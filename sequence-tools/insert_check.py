#!/usr/bin/env python3
"""
FASTQ Sequence Matcher with Command-line Sequences - Improved Version

This script searches for sequences in FASTQ files based on command-line input.
Sequences can be either:
1. Flanking sequences with extraction: "ACGT-TGCA" (extracts middle sequence)
2. Direct sequences: "ACGTACGT" (counts occurrences)

Outputs:
- Combination statistics for all matches
- Combined results file with one line per read showing all match data
- UpSet plot visualization of match combinations

Dependencies:
- biopython
- pandas
- python-Levenshtein
- matplotlib
- numpy
- upsetplot (optional, falls back to bar plot if not installed)
"""

import argparse
import sys
import logging
import time
from collections import defaultdict, Counter
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
from Bio import SeqIO, Align
from Bio.Seq import Seq
import Levenshtein
import matplotlib.pyplot as plt
import numpy as np
from upsetplot import UpSet, from_contents


# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def sanitize_name(name: str) -> str:
    """Sanitize custom sequence name for use in filenames."""
    # Replace spaces with underscores
    name = name.replace(' ', '_')
    # Remove or replace problematic characters
    name = ''.join(c for c in name if c.isalnum() or c in ('_', '-'))
    return name


def validate_sequence(seq: str) -> bool:
    """Validate that sequence contains only valid nucleotides."""
    valid_chars = set('ACGTUNRYSWKMBDHV')
    return all(c in valid_chars for c in seq.upper())


def parse_sequence_spec(seq_spec: str) -> Tuple[bool, str, Optional[str]]:
    """
    Parse sequence specification to determine if it's flanking or direct.
    
    Returns:
        tuple: (is_flanking, left_seq, right_seq) or (is_flanking, full_seq, None)
    """
    if '-' in seq_spec and seq_spec.count('-') == 1:
        # Flanking sequence format
        left, right = seq_spec.split('-')
        if left and right:  # Both parts must be non-empty
            left_clean = left.strip().upper()
            right_clean = right.strip().upper()
            
            # Validate sequences
            if not validate_sequence(left_clean) or not validate_sequence(right_clean):
                raise ValueError(f"Invalid nucleotide characters in sequence: {seq_spec}")
            
            return True, left_clean, right_clean
    
    # Direct sequence format
    seq_clean = seq_spec.strip().upper()
    if not validate_sequence(seq_clean):
        raise ValueError(f"Invalid nucleotide characters in sequence: {seq_spec}")
    
    return False, seq_clean, None


def find_fuzzy_match_optimized(sequence: str, pattern: str, max_distance: int) -> List[Tuple[int, int, int]]:
    """
    Find fuzzy matches of pattern in sequence allowing up to max_distance edits.
    Optimized version with early termination.
    
    Returns list of (start_position, end_position, distance) tuples.
    """
    matches = []
    pattern_len = len(pattern)
    sequence_upper = sequence.upper()
    
    # Early termination optimization for exact matches
    exact_pos = sequence_upper.find(pattern)
    if exact_pos != -1:
        matches.append((exact_pos, exact_pos + pattern_len, 0))
        # Continue searching for more matches
        next_pos = exact_pos + 1
        while next_pos < len(sequence) - pattern_len + 1:
            exact_pos = sequence_upper.find(pattern, next_pos)
            if exact_pos != -1:
                matches.append((exact_pos, exact_pos + pattern_len, 0))
                next_pos = exact_pos + 1
            else:
                break
    
    # If we need to look for inexact matches
    if max_distance > 0:
        for i in range(len(sequence) - pattern_len + 1):
            # Skip positions where we already found exact matches
            if any(start == i and dist == 0 for start, _, dist in matches):
                continue
                
            substring = sequence_upper[i:i + pattern_len]
            
            # Quick check: if too many character differences, skip
            char_diff = sum(1 for a, b in zip(substring, pattern) if a != b)
            if char_diff > max_distance:
                continue
                
            distance = Levenshtein.distance(substring, pattern)
            
            if distance <= max_distance and distance > 0:
                matches.append((i, i + pattern_len, distance))
    
    return sorted(matches, key=lambda x: (x[2], x[0]))  # Sort by distance, then position


def extract_middle_sequence(sequence: str, left_flank: str, right_flank: str, 
                          max_lev_distance: int, min_insert: int, max_insert: int) -> Optional[str]:
    """
    Extract sequence between two flanking sequences.
    
    Returns:
        str or None: Extracted middle sequence if found
    """
    # Find matches for both flanks
    left_matches = find_fuzzy_match_optimized(sequence, left_flank, max_lev_distance)
    right_matches = find_fuzzy_match_optimized(sequence, right_flank, max_lev_distance)
    
    if not left_matches or not right_matches:
        return None
    
    # Find best match pair within insert size constraints
    best_match = None
    best_total_distance = float('inf')
    
    for left_start, left_end, left_dist in left_matches:
        for right_start, right_end, right_dist in right_matches:
            if right_start > left_end:
                insert_size = right_start - left_end
                
                if min_insert <= insert_size <= max_insert:
                    total_dist = left_dist + right_dist
                    if total_dist < best_total_distance:
                        best_total_distance = total_dist
                        best_match = sequence[left_end:right_start]
    
    return best_match


def count_fastq_reads(fastq_file: str) -> int:
    """Count total number of reads in FASTQ file for progress estimation."""
    count = 0
    with open(fastq_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                count += 1
    return count


def process_fastq_file(fastq_file: str, sequence_specs: Dict, output_prefix: str, 
                      no_progress: bool = False, skip_count: bool = False) -> Dict:
    """
    Process FASTQ file and search for all specified sequences.
    
    Returns:
        dict: Results including combination counts and extracted sequences
    """
    # Validate file exists
    if not Path(fastq_file).exists():
        raise FileNotFoundError(f"FASTQ file not found: {fastq_file}")
    
    # Get file size for estimation
    file_size = Path(fastq_file).stat().st_size
    
    # Quick count of reads if file is small enough
    total_reads_estimate = None
    if not skip_count and file_size < 1_000_000_000:  # Less than 1GB
        logger.info("Counting reads for accurate progress reporting...")
        try:
            total_reads_estimate = count_fastq_reads(fastq_file)
            logger.info(f"Found {total_reads_estimate:,} reads in file")
        except:
            logger.warning("Could not count reads, will estimate progress based on file size")
    
    # Track matches for each read
    read_matches = []
    extracted_sequences = defaultdict(list)
    all_read_data = {}  # Store all data for each read
    read_ids_ordered = []  # Maintain order of reads
    total_reads = 0
    processed_reads = 0
    
    # Progress tracking
    start_time = time.time()
    last_progress_time = start_time
    progress_interval = 5  # seconds
    
    # Estimate bytes per read from first 1000 reads
    bytes_per_read_estimate = None
    
    try:
        with open(fastq_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fastq"):
                total_reads += 1
                
                # Estimate bytes per read from first 1000 reads
                if total_reads == 1000 and not total_reads_estimate:
                    # Rough estimate: 4 lines per read, average line length
                    avg_read_bytes = file_size / total_reads * 1000
                    bytes_per_read_estimate = avg_read_bytes
                
                # Progress reporting
                current_time = time.time()
                if not no_progress and (current_time - last_progress_time >= progress_interval or total_reads % 10000 == 0):
                    elapsed_time = current_time - start_time
                    reads_per_second = total_reads / elapsed_time if elapsed_time > 0 else 0
                    
                    # Calculate progress
                    if total_reads_estimate:
                        percent_complete = (total_reads / total_reads_estimate * 100)
                    elif bytes_per_read_estimate:
                        estimated_total_reads = file_size / bytes_per_read_estimate
                        percent_complete = (total_reads / estimated_total_reads * 100)
                    else:
                        percent_complete = 0
                    
                    # Estimate time remaining
                    if percent_complete > 0 and elapsed_time > 0:
                        total_time_estimate = elapsed_time / (percent_complete / 100)
                        time_remaining = total_time_estimate - elapsed_time
                        eta_str = f", ETA: {int(time_remaining // 60)}m {int(time_remaining % 60)}s"
                    else:
                        eta_str = ""
                    
                    logger.info(f"Progress: {total_reads:,} reads processed ({percent_complete:.1f}% complete, "
                              f"{reads_per_second:.0f} reads/sec{eta_str})")
                    last_progress_time = current_time
                
                read_id = record.id
                read_ids_ordered.append(read_id)
                sequence = str(record.seq)
                
                # Initialize read data
                read_data = {
                    'read_id': read_id,
                    'matches': {},  # seq_name -> True/False
                    'extracted': {},  # seq_name -> extracted sequence or NA
                    'match_ids': {},  # seq_name -> match ID or NO_MATCH or NA
                    'metrics': {}  # seq_name -> metrics dict or None
                }
                
                # Track which sequences matched for this read
                matches_found = set()
                
                # Check each sequence specification
                for seq_name, spec in sequence_specs.items():
                    is_flanking = spec['is_flanking']
                    max_dist = spec['max_distance']
                    
                    if is_flanking:
                        # Extract middle sequence
                        left_seq = spec['left_seq']
                        right_seq = spec['right_seq']
                        min_insert = spec['min_insert']
                        max_insert = spec['max_insert']
                        
                        # Try forward strand
                        middle = extract_middle_sequence(sequence, left_seq, right_seq, 
                                                       max_dist, min_insert, max_insert)
                        
                        if middle:
                            matches_found.add(seq_name)
                            read_data['matches'][seq_name] = True
                            
                            # Reverse complement if requested
                            if spec.get('reverse_complement', False):
                                middle = str(Seq(middle).reverse_complement())
                            
                            read_data['extracted'][seq_name] = middle
                            
                            # Match to reference if available
                            if spec['reference']:
                                match_id, ref_seq, metrics = match_sequence_to_references(
                                    middle, spec['reference'], spec['match_method'], spec['match_dist'])
                                read_data['match_ids'][seq_name] = match_id if match_id else 'NO_MATCH'
                                read_data['metrics'][seq_name] = metrics
                                extracted_sequences[seq_name].append((read_id, middle, match_id, metrics))
                            else:
                                read_data['match_ids'][seq_name] = 'NA'
                                read_data['metrics'][seq_name] = None
                                extracted_sequences[seq_name].append((read_id, middle, None, None))
                        else:
                            # Try reverse complement only if forward fails
                            rev_comp = str(record.seq.reverse_complement())
                            middle = extract_middle_sequence(rev_comp, left_seq, right_seq, 
                                                           max_dist, min_insert, max_insert)
                            if middle:
                                matches_found.add(seq_name)
                                read_data['matches'][seq_name] = True
                                
                                # Reverse complement if requested
                                if spec.get('reverse_complement', False):
                                    middle = str(Seq(middle).reverse_complement())
                                
                                read_data['extracted'][seq_name] = middle
                                
                                if spec['reference']:
                                    match_id, ref_seq, metrics = match_sequence_to_references(
                                        middle, spec['reference'], spec['match_method'], spec['match_dist'])
                                    read_data['match_ids'][seq_name] = match_id if match_id else 'NO_MATCH'
                                    read_data['metrics'][seq_name] = metrics
                                    extracted_sequences[seq_name].append((read_id, middle, match_id, metrics))
                                else:
                                    read_data['match_ids'][seq_name] = 'NA'
                                    read_data['metrics'][seq_name] = None
                                    extracted_sequences[seq_name].append((read_id, middle, None, None))
                            else:
                                # No match found
                                read_data['matches'][seq_name] = False
                                read_data['extracted'][seq_name] = 'NA'
                                read_data['match_ids'][seq_name] = 'NA'
                                read_data['metrics'][seq_name] = None
                    else:
                        # Direct sequence search
                        full_seq = spec['full_seq']
                        
                        # Try forward strand
                        if search_direct_sequence(sequence, full_seq, max_dist):
                            matches_found.add(seq_name)
                            read_data['matches'][seq_name] = True
                        else:
                            # Try reverse complement only if forward fails
                            rev_comp = str(record.seq.reverse_complement())
                            if search_direct_sequence(rev_comp, full_seq, max_dist):
                                matches_found.add(seq_name)
                                read_data['matches'][seq_name] = True
                            else:
                                read_data['matches'][seq_name] = False
                        
                        # Direct sequences don't have extraction
                        read_data['extracted'][seq_name] = 'NA'
                        read_data['match_ids'][seq_name] = 'NA'
                        read_data['metrics'][seq_name] = None
                
                # Store read data
                all_read_data[read_id] = read_data
                
                # Store match combination for this read
                if matches_found:
                    read_matches.append(tuple(sorted(matches_found)))
                else:
                    read_matches.append(())
                    
                processed_reads += 1
                
    except Exception as e:
        logger.error(f"Error processing FASTQ file: {e}")
        raise
    
    # Final progress report
    total_time = time.time() - start_time
    logger.info(f"Finished processing {total_reads:,} reads in {total_time:.1f} seconds "
                f"({total_reads/total_time:.0f} reads/sec)")
    
    # Count combinations
    combination_counts = Counter(read_matches)
    
    # Generate all possible combinations
    all_seq_names = sorted(sequence_specs.keys())
    all_combinations = []
    
    for r in range(len(all_seq_names) + 1):
        for combo in combinations(all_seq_names, r):
            all_combinations.append(combo)
    
    # Create results
    results = {
        'total_reads': total_reads,
        'processed_reads': processed_reads,
        'combination_counts': combination_counts,
        'all_combinations': all_combinations,
        'extracted_sequences': dict(extracted_sequences),
        'sequence_specs': sequence_specs,
        'all_read_data': all_read_data,
        'read_ids_ordered': read_ids_ordered
    }
    
    return results


def search_direct_sequence(sequence: str, pattern: str, max_lev_distance: int) -> bool:
    """
    Search for direct sequence matches.
    
    Returns:
        bool: True if match found
    """
    matches = find_fuzzy_match_optimized(sequence, pattern, max_lev_distance)
    return len(matches) > 0


def load_reference_sequences(reference_file: str) -> Dict[str, str]:
    """
    Load reference sequences from either FASTA or TSV file.
    
    For TSV files:
    - Column 1: sequence to match
    - Column 2: ID to report
    
    Returns:
        dict: {sequence: identifier}
    """
    references = {}
    
    if not Path(reference_file).exists():
        raise FileNotFoundError(f"Reference file not found: {reference_file}")
    
    # Try to detect format
    with open(reference_file, 'r') as f:
        first_line = f.readline().strip()
    
    if first_line.startswith('>'):
        # FASTA format
        with open(reference_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_upper = str(record.seq).upper()
                if not validate_sequence(seq_upper):
                    logger.warning(f"Skipping reference sequence {record.id} due to invalid characters")
                    continue
                references[seq_upper] = record.id
    else:
        # TSV format
        try:
            df = pd.read_csv(reference_file, sep='\t', header=None, 
                           names=['sequence', 'id'], usecols=[0, 1])
            
            for _, row in df.iterrows():
                sequence = str(row['sequence']).strip().upper()
                identifier = str(row['id']).strip()
                
                if sequence and identifier:
                    if not validate_sequence(sequence):
                        logger.warning(f"Skipping reference sequence {identifier} due to invalid characters")
                        continue
                    references[sequence] = identifier
        except Exception as e:
            raise ValueError(f"Error reading TSV file: {e}")
    
    return references


def calculate_hamming_distance(seq1: str, seq2: str) -> float:
    """
    Calculate Hamming distance between two sequences of equal length.
    Returns normalized distance (0-1, where 0 is identical).
    """
    if len(seq1) != len(seq2):
        return 1.0  # Maximum distance if lengths differ
    
    if len(seq1) == 0:
        return 0.0
    
    mismatches = sum(a != b for a, b in zip(seq1, seq2))
    return mismatches / len(seq1)


def calculate_all_metrics(query_seq: str, ref_seq: str) -> Dict[str, float]:
    """
    Calculate all three distance/similarity metrics between sequences.
    
    Returns:
        dict: {'hamming': float, 'levenshtein': float, 'homology': float}
    """
    query_upper = query_seq.upper()
    ref_upper = ref_seq.upper()
    
    # Hamming distance (normalized)
    hamming = calculate_hamming_distance(query_upper, ref_upper)
    
    # Levenshtein distance (normalized)
    lev_dist = Levenshtein.distance(query_upper, ref_upper)
    max_len = max(len(query_upper), len(ref_upper))
    levenshtein = lev_dist / max_len if max_len > 0 else 0
    
    # Homology (alignment-based similarity)
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    
    alignments = aligner.align(query_upper, ref_upper)
    homology = 0
    
    if alignments:
        alignment = alignments[0]
        aligned_seq1 = str(alignment[0])
        aligned_seq2 = str(alignment[1])
        
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) 
                     if a == b and a != '-')
        total_positions = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) 
                            if a != '-' or b != '-')
        
        if total_positions > 0:
            homology = matches / total_positions
    
    return {
        'hamming': 1 - hamming,  # Convert to similarity (1 = perfect match)
        'levenshtein': 1 - levenshtein,  # Convert to similarity
        'homology': homology
    }


def match_sequence_to_references(query_seq: str, references: Dict[str, str], 
                               match_method: str = 'homology', 
                               match_dist: Optional[float] = None) -> Tuple[Optional[str], Optional[str], Optional[Dict[str, float]]]:
    """
    Match a query sequence to reference sequences using specified method.
    
    Args:
        query_seq: Query sequence
        references: Dict of {sequence: identifier}
        match_method: 'hamming', 'levenshtein', or 'homology'
        match_dist: Maximum distance (for hamming/levenshtein) or minimum score (for homology)
    
    Returns:
        tuple: (best_match_id, best_ref_seq, all_metrics) or (None, None, None)
    """
    if not query_seq:
        return None, None, None
    
    # Check for perfect match first
    query_upper = query_seq.upper()
    if query_upper in references:
        # Perfect match found - calculate all metrics
        metrics = calculate_all_metrics(query_upper, query_upper)
        return references[query_upper], query_upper, metrics
    
    # Find best match using specified method
    best_score = float('inf') if match_method in ['hamming', 'levenshtein'] else 0
    best_match_id = None
    best_ref_seq = None
    
    for ref_seq, ref_id in references.items():
        if match_method == 'hamming':
            # For Hamming, sequences must be same length
            if len(query_upper) == len(ref_seq):
                # Calculate raw hamming distance (not normalized)
                mismatches = sum(a != b for a, b in zip(query_upper, ref_seq))
                if mismatches < best_score:
                    best_score = mismatches
                    best_match_id = ref_id
                    best_ref_seq = ref_seq
        
        elif match_method == 'levenshtein':
            # Calculate raw Levenshtein distance (not normalized)
            lev_dist = Levenshtein.distance(query_upper, ref_seq)
            if lev_dist < best_score:
                best_score = lev_dist
                best_match_id = ref_id
                best_ref_seq = ref_seq
        
        elif match_method == 'homology':
            # Use alignment-based scoring (0-1, higher is better)
            aligner = Align.PairwiseAligner()
            aligner.mode = 'local'
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -1
            
            alignments = aligner.align(query_upper, ref_seq)
            
            if alignments:
                alignment = alignments[0]
                aligned_seq1 = str(alignment[0])
                aligned_seq2 = str(alignment[1])
                
                matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) 
                             if a == b and a != '-')
                total_positions = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) 
                                    if a != '-' or b != '-')
                
                if total_positions > 0:
                    score = matches / total_positions
                    if score > best_score:
                        best_score = score
                        best_match_id = ref_id
                        best_ref_seq = ref_seq
    
    # Check if best match meets criteria
    if best_match_id:
        if match_method in ['hamming', 'levenshtein'] and best_score <= match_dist:
            # For distance metrics, score must be <= threshold
            metrics = calculate_all_metrics(query_seq, best_ref_seq)
            return best_match_id, best_ref_seq, metrics
        elif match_method == 'homology' and best_score >= match_dist:
            # For similarity metrics, score must be >= threshold
            metrics = calculate_all_metrics(query_seq, best_ref_seq)
            return best_match_id, best_ref_seq, metrics
    
    return None, None, None


def generate_reports(results: Dict, output_prefix: str) -> None:
    """
    Generate output reports including combined file and upset plot.
    """
    total_reads = results['total_reads']
    combination_counts = results['combination_counts']
    all_combinations = results['all_combinations']
    sequence_specs = results['sequence_specs']
    all_read_data = results['all_read_data']
    read_ids_ordered = results['read_ids_ordered']
    
    # Report 1: Combination statistics
    output_file = f'{output_prefix}_combination_stats.txt'
    try:
        with open(output_file, 'w') as f:
            f.write("Combination\tCount\tProportion\n")
            
            for combo in all_combinations:
                count = combination_counts.get(combo, 0)
                proportion = count / total_reads if total_reads > 0 else 0
                
                if combo:
                    combo_str = '+'.join(combo)
                else:
                    combo_str = 'none'
                
                f.write(f"{combo_str}\t{count}\t{proportion:.4f}\n")
            
            # Summary statistics
            f.write(f"\n# Summary\n")
            f.write(f"Total reads: {total_reads:,}\n")
            
            reads_with_any_match = sum(count for combo, count in combination_counts.items() if combo)
            f.write(f"Reads with at least one match: {reads_with_any_match:,} ({reads_with_any_match/total_reads*100:.2f}%)\n")
            
            # Individual sequence statistics
            f.write(f"\n# Individual sequence matches\n")
            for seq_name in sorted(sequence_specs.keys()):
                count = sum(c for combo, c in combination_counts.items() if seq_name in combo)
                f.write(f"{seq_name}: {count:,} ({count/total_reads*100:.2f}%)\n")
        
        logger.info(f"Written combination statistics to {output_file}")
        
    except Exception as e:
        logger.error(f"Error writing combination statistics: {e}")
        raise
    
    # Report 2: Combined output file with all reads
    output_file = f'{output_prefix}_combined_results.txt'
    try:
        with open(output_file, 'w') as f:
            # Build header
            header_parts = ['Read_ID']
            
            # Add columns for each sequence
            for seq_name in sorted(sequence_specs.keys()):
                header_parts.append(f'{seq_name}_match')
                
                if sequence_specs[seq_name]['is_flanking']:
                    header_parts.append(f'{seq_name}_extracted')
                    
                    if sequence_specs[seq_name]['reference']:
                        header_parts.append(f'{seq_name}_match_id')
                        header_parts.append(f'{seq_name}_hamming')
                        header_parts.append(f'{seq_name}_levenshtein')
                        header_parts.append(f'{seq_name}_homology')
            
            f.write('\t'.join(header_parts) + '\n')
            
            # Write data for each read
            for read_id in read_ids_ordered:
                if read_id not in all_read_data:
                    continue
                    
                read_data = all_read_data[read_id]
                row_parts = [read_id]
                
                for seq_name in sorted(sequence_specs.keys()):
                    # Match status (Yes/No)
                    match_status = 'Yes' if read_data['matches'].get(seq_name, False) else 'No'
                    row_parts.append(match_status)
                    
                    if sequence_specs[seq_name]['is_flanking']:
                        # Extracted sequence or NA
                        extracted = read_data['extracted'].get(seq_name, 'NA')
                        row_parts.append(extracted)
                        
                        if sequence_specs[seq_name]['reference']:
                            # Match ID
                            match_id = read_data['match_ids'].get(seq_name, 'NA')
                            row_parts.append(match_id)
                            
                            # Metrics
                            metrics = read_data['metrics'].get(seq_name)
                            if metrics:
                                row_parts.append(f"{metrics['hamming']:.3f}")
                                row_parts.append(f"{metrics['levenshtein']:.3f}")
                                row_parts.append(f"{metrics['homology']:.3f}")
                            else:
                                row_parts.append('NA')
                                row_parts.append('NA')
                                row_parts.append('NA')
                
                f.write('\t'.join(row_parts) + '\n')
        
        logger.info(f"Written combined results to {output_file}")
        
    except Exception as e:
        logger.error(f"Error writing combined results: {e}")
        raise
    
    # Report 3: Generate upset plot with percentages
    try:
        # Prepare data for upset plot
        # Create a dictionary where keys are sequence names and values are sets of read IDs
        sequence_sets = {}
        for seq_name in sorted(sequence_specs.keys()):
            matching_reads = set()
            for read_id, read_data in all_read_data.items():
                if read_data['matches'].get(seq_name, False):
                    matching_reads.add(read_id)
            sequence_sets[seq_name] = matching_reads
        
        # Add a "No_Match" category for reads with no matches
        no_match_reads = set()
        for read_id, read_data in all_read_data.items():
            if not any(read_data['matches'].get(seq_name, False) for seq_name in sequence_specs.keys()):
                no_match_reads.add(read_id)
        
        if no_match_reads:
            sequence_sets['No_Match'] = no_match_reads
        
        # Only create upset plot if there are any reads
        if sequence_sets:
            # Create upset plot data
            upset_data = from_contents(sequence_sets)
            
            # Create the plot with percentages on both axes
            fig = plt.figure(figsize=(14, 10))
            upset = UpSet(upset_data, 
                         subset_size='count',
                         intersection_plot_elements=6,
                         totals_plot_elements=3,
                         show_counts='%d',
                         sort_by='cardinality',
                         show_percentages=True)
            upset.plot(fig=fig)
            
            # Get the axes to customize labels
            axes = fig.get_axes()
            
            # Find the intersection size axis (usually the main bar plot)
            for ax in axes:
                # Check if this is the intersection size axis
                if ax.get_ylabel() and 'Intersection size' in ax.get_ylabel():
                    # Remove % from y-axis tick labels but keep the values
                    y_ticks = ax.get_yticks()
                    y_labels = [f'{int(tick)}' for tick in y_ticks if tick >= 0]
                    ax.set_yticks([tick for tick in y_ticks if tick >= 0])
                    ax.set_yticklabels(y_labels, fontsize=10)
                    ax.set_ylabel('Intersection Size', fontsize=11)
                    ax.grid(False)  # Remove grid lines
                    
                    # Customize the percentage labels above bars
                    for text in ax.texts:
                        # Remove parentheses and resize
                        label = text.get_text()
                        if '(' in label and ')' in label:
                            # Extract percentage value and remove parentheses
                            clean_label = label.replace('(', '').replace(')', '')
                            text.set_text(clean_label)
                            text.set_fontsize(6)
                
                # Check if this is the set size axis
                elif ax.get_xlabel() and 'Set size' in ax.get_xlabel():
                    ax.grid(False)  # Remove grid lines
                    # Add percentages to x-axis labels
                    x_ticks = ax.get_xticks()
                    x_labels = []
                    for tick in x_ticks:
                        if tick >= 0 and total_reads > 0:
                            pct = (tick / total_reads) * 100
                            x_labels.append(f'{int(tick)}\n{pct:.1f}')
                        else:
                            x_labels.append(f'{int(tick)}')
                    ax.set_xticklabels(x_labels, fontsize=9)
                    ax.set_xlabel('Set Size\nCount and %', fontsize=11)
            
            plt.suptitle(f'{output_prefix} Match Combinations', fontsize=16, y=0.995)
            
            output_file = f'{output_prefix}_upset_plot.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Written upset plot to {output_file}")
        else:
            logger.warning("No data found - skipping upset plot generation")
            
    except ImportError:
        logger.warning("upsetplot package not installed - generating standard bar plot with percentages instead")
        
        # Fall back to bar plot with percentages
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.grid(False)  # Remove grid lines
        
        # Include all combinations, including empty (no matches)
        plot_data = [(combo, count) for combo, count in combination_counts.items()]
        plot_data.sort(key=lambda x: x[1], reverse=True)
        
        if plot_data and total_reads > 0:
            combos = ['+'.join(c) if c else 'No_Match' for c, _ in plot_data]
            percentages = [(count / total_reads * 100) for _, count in plot_data]
            
            bars = ax.bar(range(len(combos)), percentages)
            
            # Add percentage labels on top of bars without parentheses
            for i, (bar, pct) in enumerate(zip(bars, percentages)):
                if pct > 0:  # Only label non-zero bars
                    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                            f'{pct:.1f}%', ha='center', va='bottom', 
                            fontsize=6)
            
            ax.set_xlabel('Sequence Combination')
            ax.set_ylabel('Percentage of Reads')
            ax.set_title(f'{output_prefix} Match Combinations')
            ax.set_xticks(range(len(combos)))
            ax.set_xticklabels(combos, rotation=45, ha='right')
            ax.set_ylim(0, max(percentages) * 1.15 if percentages else 100)  # Add space for labels
            plt.tight_layout()
            
            output_file = f'{output_prefix}_combination_plot.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Written combination plot to {output_file}")
        else:
            logger.warning("No data found - skipping plot generation")
            
    except Exception as e:
        logger.error(f"Error generating upset plot: {e}")
        logger.warning("Upset plot generation failed - continuing without plot")


def validate_configuration(sequence_specs: Dict[str, Dict]) -> None:
    """Validate configuration parameters for all sequences."""
    for seq_name, spec in sequence_specs.items():
        # Validate insert sizes for flanking sequences
        if spec['is_flanking']:
            if spec['min_insert'] >= spec['max_insert']:
                raise ValueError(f"{seq_name}: min_insert ({spec['min_insert']}) must be less than max_insert ({spec['max_insert']})")
            
            if spec['min_insert'] < 0:
                raise ValueError(f"{seq_name}: min_insert must be non-negative")
            
            if spec['max_insert'] > 10000:
                logger.warning(f"{seq_name}: max_insert ({spec['max_insert']}) is very large, this may impact performance")
        
        # Validate match distance thresholds
        if spec['reference'] and spec['match_dist'] is not None:
            if spec['match_method'] in ['hamming', 'levenshtein']:
                if spec['match_dist'] < 0:
                    raise ValueError(f"{seq_name}: match distance must be non-negative for {spec['match_method']}")
            elif spec['match_method'] == 'homology':
                if not 0 <= spec['match_dist'] <= 1:
                    raise ValueError(f"{seq_name}: homology score must be between 0 and 1")
        
        # Validate max distance
        if spec['max_distance'] < 0:
            raise ValueError(f"{seq_name}: max Levenshtein distance must be non-negative")
        elif spec['max_distance'] > 10:
            logger.warning(f"{seq_name}: max Levenshtein distance ({spec['max_distance']}) is large, this may impact performance and specificity")


def main() -> None:
    """Main function for FASTQ sequence matcher."""
    parser = argparse.ArgumentParser(
        description='FASTQ sequence matcher with command-line sequences',
        epilog='''
Examples:
  # Flanking sequences with extraction and custom names
  python script.py input.fastq --seq1 "ATCG-TGCA" --name1 "barcode" --seq1-dist 2 --seq1-insert 10:30 --ref-file1 barcodes.tsv --match-method1 homology --match-dist1 0.8
  
  # Direct sequence matching with custom name
  python script.py input.fastq --seq2 "ATCGATCG" --name2 "promoter" --seq2-dist 1
  
  # Multiple sequences with custom names and reverse complement
  python script.py input.fastq --seq1 "ATCG-TGCA" --name1 "bc1" --seq1-insert 10:30 --match-method1 levenshtein --match-dist1 3 --seq1-rc \\
                               --seq2 "GGCCTTAA" --name2 "marker" \\
                               --seq3 "ACTG-CAGT" --name3 "bc2" --seq3-dist 3 --seq3-insert 5:25
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('fastq_file', help='Input FASTQ file')
    parser.add_argument('--output-prefix', default='output',
                        help='Prefix for output files (default: output)')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        default='INFO', help='Logging level (default: INFO)')
    parser.add_argument('--no-progress', action='store_true',
                        help='Disable progress reporting (useful for very large files)')
    parser.add_argument('--skip-count', action='store_true',
                        help='Skip initial read counting for faster startup')
    
    # Allow multiple sequence specifications
    for i in range(1, 11):  # Support up to 10 sequences
        parser.add_argument(f'--seq{i}', help=f'Sequence {i} (format: "ATCG-TGCA" or "ATCGATCG")')
        parser.add_argument(f'--name{i}', help=f'Custom name for sequence {i} (default: seq{i})')
        parser.add_argument(f'--seq{i}-dist', type=int, default=None, 
                          help=f'Max Levenshtein distance for seq{i} (default: 2)')
        parser.add_argument(f'--seq{i}-insert', help=f'Min:max insert size for flanking seq{i} (e.g., "10:30")')
        parser.add_argument(f'--ref-file{i}', help=f'Reference file for seq{i} matches (FASTA or TSV)')
        parser.add_argument(f'--match-method{i}', choices=['hamming', 'levenshtein', 'homology'],
                          default='homology', help=f'Matching method for seq{i} (default: homology)')
        parser.add_argument(f'--match-dist{i}', type=float, default=None,
                          help=f'Match distance/score for seq{i} (default: 5 for hamming/levenshtein, 0.8 for homology)')
        parser.add_argument(f'--seq{i}-rc', action='store_true',
                          help=f'Reverse complement extracted sequence before matching for seq{i}')
    
    args = parser.parse_args()
    
    # Set logging level
    logger.setLevel(getattr(logging, args.log_level))
    
    # Validate FASTQ file exists
    if not Path(args.fastq_file).exists():
        logger.error(f"FASTQ file not found: {args.fastq_file}")
        sys.exit(1)
    
    # Parse sequence specifications
    sequence_specs: Dict[str, Dict] = {}
    
    for i in range(1, 11):
        seq_arg = getattr(args, f'seq{i}')
        if seq_arg:
            # Get custom name or use default
            custom_name = getattr(args, f'name{i}')
            if custom_name:
                seq_name = sanitize_name(custom_name)
                logger.info(f"Using custom name '{seq_name}' for sequence {i}")
            else:
                seq_name = f'seq{i}'
            
            try:
                # Get distance parameter (default to 2 if not specified)
                dist = getattr(args, f'seq{i}_dist')
                if dist is None:
                    dist = 2
                    logger.info(f"Using default Levenshtein distance of 2 for {seq_name}")
                
                # Parse sequence specification
                is_flanking, seq1, seq2 = parse_sequence_spec(seq_arg)
                
                # Get match method
                match_method = getattr(args, f'match_method{i}')
                
                # Get match distance/score with appropriate default
                match_dist = getattr(args, f'match_dist{i}')
                if match_dist is None:
                    if match_method in ['hamming', 'levenshtein']:
                        match_dist = 5
                        logger.info(f"Using default match distance of 5 for {seq_name} ({match_method})")
                    else:  # homology
                        match_dist = 0.8
                        logger.info(f"Using default match score of 0.8 for {seq_name} (homology)")
                
                spec = {
                    'is_flanking': is_flanking,
                    'max_distance': dist,
                    'match_method': match_method,
                    'match_dist': match_dist,
                    'reverse_complement': getattr(args, f'seq{i}_rc')
                }
                
                if is_flanking:
                    spec['left_seq'] = seq1
                    spec['right_seq'] = seq2
                    
                    # Get insert size
                    insert_arg = getattr(args, f'seq{i}_insert')
                    if not insert_arg:
                        logger.error(f"--seq{i}-insert is required for flanking sequences")
                        sys.exit(1)
                    
                    try:
                        min_insert, max_insert = map(int, insert_arg.split(':'))
                        spec['min_insert'] = min_insert
                        spec['max_insert'] = max_insert
                    except:
                        logger.error(f"--seq{i}-insert must be in format 'min:max' (e.g., '10:30')")
                        sys.exit(1)
                else:
                    spec['full_seq'] = seq1
                    spec['min_insert'] = None
                    spec['max_insert'] = None
                
                # Load reference file if provided
                ref_file = getattr(args, f'ref_file{i}')
                if ref_file:
                    if not is_flanking:
                        logger.warning(f"Reference file for {seq_name} ignored (only used with flanking sequences)")
                        spec['reference'] = None
                    else:
                        if not Path(ref_file).exists():
                            logger.error(f"Reference file not found: {ref_file}")
                            sys.exit(1)
                        spec['reference'] = load_reference_sequences(ref_file)
                        logger.info(f"Loaded {len(spec['reference'])} reference sequences for {seq_name}")
                        logger.info(f"Using {spec['match_method']} matching for {seq_name} with threshold {match_dist}")
                        if spec['reverse_complement']:
                            logger.info(f"Will reverse complement extracted sequences for {seq_name}")
                else:
                    spec['reference'] = None
                
                sequence_specs[seq_name] = spec
                
            except ValueError as e:
                logger.error(f"Error parsing {seq_name}: {e}")
                sys.exit(1)
    
    if not sequence_specs:
        logger.error("At least one sequence must be specified")
        sys.exit(1)
    
    # Validate configuration
    try:
        validate_configuration(sequence_specs)
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        sys.exit(1)
    
    # Print configuration
    logger.info(f"Processing FASTQ file: {args.fastq_file}")
    logger.info(f"Number of sequences to search: {len(sequence_specs)}")
    
    for seq_name, spec in sequence_specs.items():
        if spec['is_flanking']:
            logger.info(f"  {seq_name}: {spec['left_seq']}-{spec['right_seq']} (flanking, dist={spec['max_distance']}, insert={spec['min_insert']}-{spec['max_insert']})")
        else:
            logger.info(f"  {seq_name}: {spec['full_seq']} (direct, dist={spec['max_distance']})")
    
    # Process FASTQ file
    logger.info("Processing FASTQ file...")
    try:
        results = process_fastq_file(args.fastq_file, sequence_specs, args.output_prefix,
                                   no_progress=args.no_progress, skip_count=args.skip_count)
    except Exception as e:
        logger.error(f"Error processing FASTQ file: {e}")
        sys.exit(1)
    
    # Generate reports
    logger.info("Generating reports...")
    try:
        generate_reports(results, args.output_prefix)
    except Exception as e:
        logger.error(f"Error generating reports: {e}")
        sys.exit(1)
    
    # Print summary statistics
    logger.info("\n" + "="*60)
    logger.info("SUMMARY STATISTICS")
    logger.info("="*60)
    logger.info(f"Total reads processed: {results['total_reads']:,}")
    
    reads_with_matches = sum(count for combo, count in results['combination_counts'].items() if combo)
    logger.info(f"Reads with at least one match: {reads_with_matches:,} ({reads_with_matches/results['total_reads']*100:.2f}%)")
    
    logger.info("\nIndividual sequence match rates:")
    for seq_name in sorted(sequence_specs.keys()):
        count = sum(c for combo, c in results['combination_counts'].items() if seq_name in combo)
        logger.info(f"  {seq_name}: {count:,} reads ({count/results['total_reads']*100:.2f}%)")
    
    logger.info("\nTop 5 most common combinations:")
    sorted_combos = sorted(results['combination_counts'].items(), key=lambda x: x[1], reverse=True)[:5]
    for combo, count in sorted_combos:
        combo_str = '+'.join(combo) if combo else 'none'
        logger.info(f"  {combo_str}: {count:,} reads ({count/results['total_reads']*100:.2f}%)")
    
    logger.info("="*60 + "\n")
    
    logger.info("Processing complete!")
    logger.info(f"Results written to:")
    logger.info(f"  - {args.output_prefix}_combination_stats.txt (combination statistics)")
    logger.info(f"  - {args.output_prefix}_combined_results.txt (all reads with match data)")
    
    # Check if upset plot was created
    if Path(f"{args.output_prefix}_upset_plot.png").exists():
        logger.info(f"  - {args.output_prefix}_upset_plot.png (upset plot visualization)")
    elif Path(f"{args.output_prefix}_combination_plot.png").exists():
        logger.info(f"  - {args.output_prefix}_combination_plot.png (combination plot visualization)")


if __name__ == "__main__":
    main()
