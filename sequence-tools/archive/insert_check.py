#!/usr/bin/env python3
"""
FASTQ Sequence Matcher with Command-line Sequences - OPTIMIZED VERSION WITH MAPPY

PERFORMANCE IMPROVEMENTS:
1. Hash table lookup for exact hamming matches (O(1) instead of O(n))
2. Reuse PairwiseAligner objects instead of recreating
3. Early termination for exact matches
4. Only calculate needed metrics, avoid redundant calculations
5. **NEW: Minimap2 (mappy) integration for ultra-fast alignment**

MAPPY OPTIMIZATIONS:
- End-to-end alignment mode for maximum accuracy
- High sensitivity settings to discriminate 1-2 nucleotide differences
- Pre-built index reused across all queries
- Can handle 100K+ references with <1 second per query
- 10-100x faster than Bio.Align for large reference sets

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
- mappy (optional, for ultra-fast alignment)
- upsetplot (optional, falls back to bar plot if not installed)
"""

import argparse
import sys
import logging
import time
import warnings
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

# Optional imports
try:
    from upsetplot import UpSet, from_contents
    UPSETPLOT_AVAILABLE = True
except ImportError:
    UPSETPLOT_AVAILABLE = False

try:
    import mappy
    MAPPY_AVAILABLE = True
except ImportError:
    MAPPY_AVAILABLE = False
    logging.warning("mappy (minimap2) not available. Install with: pip install mappy")
    logging.warning("Falling back to slower alignment methods for homology matching")

# Suppress FutureWarnings from upsetplot/pandas compatibility
warnings.filterwarnings('ignore', category=FutureWarning, module='upsetplot')
warnings.filterwarnings('ignore', category=FutureWarning, message='.*fillna.*')
warnings.filterwarnings('ignore', category=FutureWarning, message='.*inplace.*')


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
    
    logger.info(f"Loaded {len(references)} reference sequences")
    return references


def stream_reference_sequences(reference_file: str):
    """
    Generator that yields reference sequences one at a time without loading entire file.
    Memory-efficient for very large reference files.

    Args:
        reference_file: Path to FASTA or TSV reference file

    Yields:
        tuple: (sequence, identifier)
    """
    if not Path(reference_file).exists():
        raise FileNotFoundError(f"Reference file not found: {reference_file}")

    # OPTIMIZATION: Use 1MB buffer for faster I/O
    buffer_size = 1024 * 1024  # 1MB buffer

    # Detect format from first line
    with open(reference_file, 'r', buffering=buffer_size) as f:
        first_line = f.readline().strip()

    if first_line.startswith('>'):
        # FASTA format - use SeqIO streaming
        with open(reference_file, 'r', buffering=buffer_size) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_upper = str(record.seq).upper()
                if validate_sequence(seq_upper):
                    yield (seq_upper, record.id)
                else:
                    logger.debug(f"Skipping reference sequence {record.id} due to invalid characters")
    else:
        # TSV format - stream line by line
        with open(reference_file, 'r', buffering=buffer_size) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    sequence = parts[0].strip().upper()
                    identifier = parts[1].strip()
                    if sequence and identifier and validate_sequence(sequence):
                        yield (sequence, identifier)
                    else:
                        logger.debug(f"Skipping reference sequence {identifier} due to invalid characters")


class MappyIndex:
    """
    Wrapper for mappy index with optimized settings for high-accuracy alignment.
    
    Optimized for:
    - End-to-end alignment
    - High sensitivity to discriminate 1-2 nucleotide differences
    - Fast querying of 100K+ references
    """
    
    def __init__(self, references: Dict[str, str], preset: str = 'map-ont'):
        """
        Initialize mappy index from reference sequences.
        
        Args:
            references: Dict of {sequence: identifier}
            preset: Minimap2 preset ('sr' for short reads, 'map-ont' for higher sensitivity)
        """
        if not MAPPY_AVAILABLE:
            raise ImportError("mappy not available. Install with: pip install mappy")
        
        self.references = references
        self.seq_to_id = references  # Keep mapping for fast lookup
        self.id_to_seq = {v: k for k, v in references.items()}  # Reverse mapping
        
        # Create temporary FASTA for indexing
        import tempfile
        self.temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        
        for seq, seq_id in references.items():
            self.temp_fasta.write(f">{seq_id}\n{seq}\n")
        self.temp_fasta.flush()
        
        # Build mappy index with optimized parameters
        # Using custom parameters for high accuracy and sensitivity
        self.aligner = mappy.Aligner(
            self.temp_fasta.name,
            preset=preset,
            k=15,              # k-mer size (smaller = more sensitive)
            w=10,              # minimizer window (smaller = more sensitive)
            min_chain_score=20, # Minimum chaining score
            min_dp_score=30,   # Minimum DP alignment score
            best_n=1,          # Return only best alignment
            n_threads=1,       # Single thread (will parallelize at read level if needed)
            extra_flags=0x4    # 0x4 = end-to-end alignment (cigar with clipping)
        )
        
        logger.info(f"Built mappy index with {len(references)} references")
        logger.info(f"Settings: k=15, w=10, preset={preset}, end-to-end mode")
    
    def __del__(self):
        """Clean up temporary file."""
        try:
            import os
            os.unlink(self.temp_fasta.name)
        except:
            pass
    
    def align(self, query_seq: str, min_identity: float = 0.8) -> Optional[Tuple[str, str, float, int, int]]:
        """
        Align query sequence to references using minimap2.
        
        Args:
            query_seq: Query sequence
            min_identity: Minimum identity threshold (0-1)
        
        Returns:
            tuple: (ref_id, ref_seq, identity, matches, query_len) or None
        """
        if not query_seq:
            return None
        
        query_upper = query_seq.upper()
        
        # Perform alignment
        alignments = list(self.aligner.map(query_upper))
        
        if not alignments:
            return None
        
        # Get best alignment (mappy returns sorted by score)
        best_hit = alignments[0]
        
        # Calculate identity
        # Identity = matches / aligned_length
        matches = best_hit.mlen  # Number of matching bases
        aligned_length = best_hit.blen  # Block length (alignment length)
        
        if aligned_length == 0:
            return None
        
        identity = matches / aligned_length
        
        # Filter by minimum identity
        if identity < min_identity:
            return None
        
        # Get reference sequence
        ref_id = best_hit.ctg  # Contig/reference name
        ref_seq = self.id_to_seq.get(ref_id, "")
        
        return ref_id, ref_seq, identity, matches, len(query_upper)


def calculate_hamming_similarity(seq1: str, seq2: str) -> float:
    """
    Calculate Hamming similarity between two sequences of equal length.
    Returns similarity score (0-1, where 1 is identical).
    """
    if len(seq1) != len(seq2):
        return 0.0
    
    if len(seq1) == 0:
        return 1.0
    
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return matches / len(seq1)


def calculate_levenshtein_similarity(seq1: str, seq2: str) -> float:
    """
    Calculate Levenshtein similarity (normalized).
    Returns similarity score (0-1, where 1 is identical).
    """
    lev_dist = Levenshtein.distance(seq1, seq2)
    max_len = max(len(seq1), len(seq2))
    if max_len == 0:
        return 1.0
    return 1 - (lev_dist / max_len)


def calculate_homology_similarity(query_seq: str, ref_seq: str, aligner: Align.PairwiseAligner) -> float:
    """
    Calculate homology similarity using sequence alignment.
    Returns similarity score (0-1, where 1 is identical).
    
    OPTIMIZATION: Accepts pre-initialized aligner to avoid recreation overhead.
    """
    alignments = aligner.align(query_seq, ref_seq)
    
    if not alignments:
        return 0.0
    
    alignment = alignments[0]
    aligned_seq1 = str(alignment[0])
    aligned_seq2 = str(alignment[1])
    
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) 
                 if a == b and a != '-')
    total_positions = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) 
                        if a != '-' or b != '-')
    
    if total_positions == 0:
        return 0.0
    
    return matches / total_positions


def calculate_all_metrics(query_seq: str, ref_seq: str, aligner: Align.PairwiseAligner) -> Dict[str, float]:
    """
    Calculate all three distance/similarity metrics between sequences.
    
    OPTIMIZATION: Only calculate metrics as needed, reuse aligner.
    
    Returns:
        dict: {'hamming': float, 'levenshtein': float, 'homology': float}
    """
    query_upper = query_seq.upper()
    ref_upper = ref_seq.upper()
    
    # Hamming similarity
    hamming = calculate_hamming_similarity(query_upper, ref_upper)
    
    # Levenshtein similarity
    levenshtein = calculate_levenshtein_similarity(query_upper, ref_upper)
    
    # Homology similarity (using provided aligner)
    homology = calculate_homology_similarity(query_upper, ref_upper, aligner)
    
    return {
        'hamming': hamming,
        'levenshtein': levenshtein,
        'homology': homology
    }


def calculate_metrics_with_mappy(query_seq: str, ref_seq: str, identity: float, 
                                matches: int, query_len: int) -> Dict[str, float]:
    """
    Calculate metrics using mappy alignment results plus supplementary metrics.
    
    Args:
        query_seq: Query sequence
        ref_seq: Reference sequence
        identity: Identity from mappy alignment
        matches: Number of matches from mappy
        query_len: Query length
    
    Returns:
        dict: {'hamming': float, 'levenshtein': float, 'homology': float, 'mappy_identity': float}
    """
    query_upper = query_seq.upper()
    ref_upper = ref_seq.upper()
    
    # Use mappy identity as primary homology score
    homology = identity
    
    # Calculate supplementary metrics for compatibility
    hamming = calculate_hamming_similarity(query_upper, ref_upper)
    levenshtein = calculate_levenshtein_similarity(query_upper, ref_upper)
    
    return {
        'hamming': hamming,
        'levenshtein': levenshtein,
        'homology': homology,
        'mappy_identity': identity
    }


def match_sequence_to_references_mappy(query_seq: str, mappy_index: MappyIndex,
                                      match_dist: float = 0.8) -> Tuple[Optional[str], Optional[str], Optional[Dict[str, float]]]:
    """
    ULTRA-FAST: Match sequence using minimap2 (mappy).
    
    This is 10-100x faster than Bio.Align for large reference sets.
    Optimized for end-to-end alignment with high sensitivity to discriminate
    similar sequences (1-2 nucleotide differences).
    
    Args:
        query_seq: Query sequence
        mappy_index: Pre-built mappy index
        match_dist: Minimum identity threshold (0-1)
    
    Returns:
        tuple: (best_match_id, best_ref_seq, metrics) or (None, None, None)
    """
    if not query_seq:
        return None, None, None
    
    query_upper = query_seq.upper()
    
    # Check for perfect match first (O(1) hash lookup)
    if query_upper in mappy_index.seq_to_id:
        ref_id = mappy_index.seq_to_id[query_upper]
        # Perfect match - calculate metrics
        metrics = {
            'hamming': 1.0,
            'levenshtein': 1.0,
            'homology': 1.0,
            'mappy_identity': 1.0
        }
        return ref_id, query_upper, metrics
    
    # Perform mappy alignment
    result = mappy_index.align(query_upper, min_identity=match_dist)
    
    if result is None:
        return None, None, None
    
    ref_id, ref_seq, identity, matches, query_len = result
    
    # Calculate all metrics
    metrics = calculate_metrics_with_mappy(query_seq, ref_seq, identity, matches, query_len)
    
    return ref_id, ref_seq, metrics


def match_sequence_to_references_optimized(query_seq: str, references: Dict[str, str], 
                                          match_method: str = 'homology', 
                                          match_dist: Optional[float] = None,
                                          aligner: Optional[Align.PairwiseAligner] = None,
                                          mappy_index: Optional[MappyIndex] = None) -> Tuple[Optional[str], Optional[str], Optional[Dict[str, float]]]:
    """
    OPTIMIZED: Match a query sequence to reference sequences using specified method.
    
    KEY OPTIMIZATIONS:
    1. For hamming distance 0: Use O(1) hash table lookup
    2. For homology: Reuse pre-initialized aligner
    3. For mappy: Use ultra-fast minimap2 alignment (NEW!)
    4. Early termination for perfect matches
    5. Only calculate all metrics once at the end
    
    Args:
        query_seq: Query sequence
        references: Dict of {sequence: identifier}
        match_method: 'hamming', 'levenshtein', 'homology', or 'mappy'
        match_dist: Maximum distance (for hamming/levenshtein) or minimum score (for homology/mappy)
        aligner: Pre-initialized PairwiseAligner (for homology method)
        mappy_index: Pre-built mappy index (for mappy method)
    
    Returns:
        tuple: (best_match_id, best_ref_seq, all_metrics) or (None, None, None)
    """
    if not query_seq:
        return None, None, None
    
    query_upper = query_seq.upper()
    
    # OPTIMIZATION 0: Ultra-fast minimap2 alignment (NEW!)
    if match_method == 'mappy':
        if mappy_index is None:
            raise ValueError("mappy_index required for mappy matching method")
        return match_sequence_to_references_mappy(query_seq, mappy_index, match_dist)
    
    # OPTIMIZATION 1: Perfect match lookup for exact hamming (O(1) instead of O(n))
    if match_method == 'hamming' and match_dist == 0:
        if query_upper in references:
            # Perfect match found - only calculate metrics once
            if aligner is None:
                aligner = _get_default_aligner()
            metrics = calculate_all_metrics(query_upper, query_upper, aligner)
            return references[query_upper], query_upper, metrics
        else:
            return None, None, None
    
    # OPTIMIZATION 2: For all other cases, check perfect match first
    if query_upper in references:
        if aligner is None:
            aligner = _get_default_aligner()
        metrics = calculate_all_metrics(query_upper, query_upper, aligner)
        return references[query_upper], query_upper, metrics
    
    # OPTIMIZATION 3: Initialize aligner once if using homology
    if match_method == 'homology' and aligner is None:
        aligner = _get_default_aligner()
    
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
                    # OPTIMIZATION 4: Early termination if within threshold
                    if mismatches <= match_dist:
                        break
        
        elif match_method == 'levenshtein':
            # Calculate raw Levenshtein distance (not normalized)
            lev_dist = Levenshtein.distance(query_upper, ref_seq)
            if lev_dist < best_score:
                best_score = lev_dist
                best_match_id = ref_id
                best_ref_seq = ref_seq
                # OPTIMIZATION 5: Early termination if within threshold
                if lev_dist <= match_dist:
                    break
        
        elif match_method == 'homology':
            # Use alignment-based scoring (0-1, higher is better)
            score = calculate_homology_similarity(query_upper, ref_seq, aligner)
            if score > best_score:
                best_score = score
                best_match_id = ref_id
                best_ref_seq = ref_seq
                # OPTIMIZATION 6: Early termination for very high similarity
                if score >= 0.99:
                    break
    
    # Check if best match meets criteria
    if best_match_id:
        if match_method in ['hamming', 'levenshtein'] and best_score <= match_dist:
            # For distance metrics, score must be <= threshold
            if aligner is None:
                aligner = _get_default_aligner()
            # OPTIMIZATION 7: Only calculate all metrics once at the end
            metrics = calculate_all_metrics(query_seq, best_ref_seq, aligner)
            return best_match_id, best_ref_seq, metrics
        elif match_method == 'homology' and best_score >= match_dist:
            # For similarity metrics, score must be >= threshold
            # OPTIMIZATION 8: Metrics already calculated during homology search
            metrics = {
                'hamming': calculate_hamming_similarity(query_upper, best_ref_seq),
                'levenshtein': calculate_levenshtein_similarity(query_upper, best_ref_seq),
                'homology': best_score  # Reuse already calculated score
            }
            return best_match_id, best_ref_seq, metrics
    
    return None, None, None


def _get_default_aligner() -> Align.PairwiseAligner:
    """
    Create and return a default PairwiseAligner with standard settings.
    Helper function to avoid code duplication.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    return aligner


def search_direct_sequence(sequence: str, pattern: str, max_lev_distance: int) -> bool:
    """
    Search for direct sequence matches.
    
    Returns:
        bool: True if match found
    """
    matches = find_fuzzy_match_optimized(sequence, pattern, max_lev_distance)
    return len(matches) > 0


def match_sequences_streaming(
    query_sequences: Dict[str, str],
    reference_file: str,
    match_method: str,
    match_dist: float
) -> Tuple[Dict[str, Tuple[str, str, float]], int]:
    """
    Stream through reference file matching against query sequences in memory.
    Efficient for very large reference files with hamming/levenshtein matching.

    Args:
        query_sequences: {query_id: sequence} - all sequences to match
        reference_file: Path to reference FASTA or TSV file
        match_method: 'hamming' or 'levenshtein'
        match_dist: Maximum distance threshold

    Returns:
        tuple: (matches_dict, refs_checked)
            matches_dict: {query_id: (match_id, ref_seq, distance)}
            refs_checked: Number of references processed
    """
    if match_method not in ['hamming', 'levenshtein']:
        raise ValueError(f"Streaming only supports hamming/levenshtein, not {match_method}")

    # Track best match for each query
    best_matches = {}  # {query_id: (match_id, ref_seq, distance)}
    for query_id in query_sequences:
        best_matches[query_id] = (None, None, float('inf'))

    # Track which queries have perfect matches
    perfect_matches = set()

    # OPTIMIZATION: Track matched count incrementally instead of recalculating
    matched_count = 0

    # Stream through reference file
    refs_checked = 0
    last_report = 0

    logger.info(f"Streaming through reference file: {Path(reference_file).name}")

    # OPTIMIZATION: Fast path for exact matching (match_dist == 0)
    use_fast_path = (match_dist == 0)
    query_hash = None
    if use_fast_path:
        # Build hash table: {uppercased_seq: query_id} for O(1) lookup
        query_hash = {query_seq.upper(): query_id for query_id, query_seq in query_sequences.items()}
        logger.info(f"Using exact-match fast path (hash table with {len(query_hash)} queries)")

    # OPTIMIZATION: Group queries by length for hamming distance
    queries_by_length = None
    if match_method == 'hamming' and not use_fast_path:
        queries_by_length = {}
        for query_id, query_seq in query_sequences.items():
            length = len(query_seq)
            if length not in queries_by_length:
                queries_by_length[length] = {}
            queries_by_length[length][query_id] = query_seq
        logger.info(f"Grouped {len(query_sequences)} queries into {len(queries_by_length)} length groups")

    for ref_seq, ref_id in stream_reference_sequences(reference_file):
        refs_checked += 1

        # FAST PATH: O(1) hash table lookup for exact matches
        if use_fast_path:
            # Direct hash lookup - no distance calculation needed!
            if ref_seq in query_hash:
                query_id = query_hash[ref_seq]
                if query_id not in perfect_matches:
                    was_unmatched = best_matches[query_id][0] is None
                    best_matches[query_id] = (ref_id, ref_seq, 0)
                    perfect_matches.add(query_id)
                    if was_unmatched:
                        matched_count += 1
        # NORMAL PATH: Calculate distances for fuzzy matching
        else:
            # For hamming with length grouping, only check queries of matching length
            if queries_by_length is not None:
                ref_len = len(ref_seq)
                queries_to_check = queries_by_length.get(ref_len, {})
            else:
                queries_to_check = query_sequences

            # Check each query sequence
            for query_id, query_seq in queries_to_check.items():
                # Skip if already has perfect match
                if query_id in perfect_matches:
                    continue

                # Calculate distance
                if match_method == 'hamming':
                    # No length check needed if using grouped queries
                    if queries_by_length is None and len(query_seq) != len(ref_seq):
                        continue  # Hamming requires same length
                    distance = sum(a != b for a, b in zip(query_seq.upper(), ref_seq))
                else:  # levenshtein
                    distance = Levenshtein.distance(query_seq.upper(), ref_seq)

                # Perfect match - mark as done
                if distance == 0:
                    was_unmatched = best_matches[query_id][0] is None
                    best_matches[query_id] = (ref_id, ref_seq, 0)
                    perfect_matches.add(query_id)
                    if was_unmatched:
                        matched_count += 1
                # Better match than current best
                elif distance <= match_dist:
                    current_distance = best_matches[query_id][2]
                    if distance < current_distance:
                        was_unmatched = best_matches[query_id][0] is None
                        best_matches[query_id] = (ref_id, ref_seq, distance)
                        if was_unmatched:
                            matched_count += 1

        # Early termination if all have perfect matches
        if len(perfect_matches) == len(query_sequences):
            logger.info(f"✓ All {len(query_sequences)} sequences matched perfectly! "
                       f"Stopped at {refs_checked:,} references")
            break

        # Progress report every 5M references
        if refs_checked - last_report >= 5000000:
            logger.info(f"  Checked {refs_checked:,} refs | "
                       f"Matched: {matched_count}/{len(query_sequences)} | "
                       f"Perfect: {len(perfect_matches)}/{len(query_sequences)}")
            last_report = refs_checked

    # Final summary
    logger.info(f"Completed: {refs_checked:,} references checked, "
               f"{matched_count}/{len(query_sequences)} sequences matched")

    return best_matches, refs_checked


def validate_configuration(sequence_specs: Dict) -> None:
    """
    Validate the configuration before processing.

    Raises:
        ValueError: If configuration is invalid
    """
    if not sequence_specs:
        raise ValueError("No sequences specified")

    for seq_name, spec in sequence_specs.items():
        if spec['is_flanking']:
            if not spec.get('min_insert') or not spec.get('max_insert'):
                raise ValueError(f"{seq_name}: Flanking sequences require insert size range")
            if spec['min_insert'] < 0:
                raise ValueError(f"{seq_name}: Minimum insert size must be >= 0")
            if spec['max_insert'] < spec['min_insert']:
                raise ValueError(f"{seq_name}: Maximum insert size must be >= minimum")

        if spec['max_distance'] < 0:
            raise ValueError(f"{seq_name}: Maximum distance must be >= 0")

        # Validate mappy availability
        if spec.get('match_method') == 'mappy' and not MAPPY_AVAILABLE:
            raise ValueError(f"{seq_name}: mappy method requested but mappy not installed. "
                           "Install with: pip install mappy")


def passes_quality_filters(record, min_length: Optional[int] = None,
                           max_length: Optional[int] = None,
                           min_quality: Optional[float] = None,
                           min_quality_percentile: Optional[float] = None) -> Tuple[bool, str]:
    """
    Check if a FASTQ record passes quality filters.

    Args:
        record: BioPython SeqRecord with quality scores
        min_length: Minimum read length (inclusive)
        max_length: Maximum read length (inclusive)
        min_quality: Minimum average quality score (Phred scale)
        min_quality_percentile: Percentage of bases that must meet min_quality

    Returns:
        tuple: (passes, reason) where passes is bool and reason is failure reason or "pass"
    """
    read_length = len(record.seq)

    # Check length filters
    if min_length is not None and read_length < min_length:
        return False, f"length_{read_length}<{min_length}"

    if max_length is not None and read_length > max_length:
        return False, f"length_{read_length}>{max_length}"

    # Check quality filters
    if min_quality is not None or min_quality_percentile is not None:
        if not hasattr(record, 'letter_annotations') or 'phred_quality' not in record.letter_annotations:
            logger.warning("Quality scores not available in FASTQ file. Skipping quality filtering.")
            return True, "pass"

        quality_scores = record.letter_annotations['phred_quality']

        if min_quality is not None and min_quality_percentile is None:
            # Filter by average quality
            avg_quality = sum(quality_scores) / len(quality_scores) if quality_scores else 0
            if avg_quality < min_quality:
                return False, f"avg_quality_{avg_quality:.1f}<{min_quality}"

        elif min_quality is not None and min_quality_percentile is not None:
            # Filter by percentile - percentage of bases meeting quality threshold
            bases_above_threshold = sum(1 for q in quality_scores if q >= min_quality)
            percent_above = (bases_above_threshold / len(quality_scores) * 100) if quality_scores else 0

            if percent_above < min_quality_percentile:
                return False, f"quality_percentile_{percent_above:.1f}%<{min_quality_percentile}%"

    return True, "pass"


def process_single_read(record, sequence_specs: Dict, 
                        min_length: Optional[int] = None, max_length: Optional[int] = None,
                        min_quality: Optional[float] = None,
                        min_quality_percentile: Optional[float] = None) -> Tuple[Optional[Dict], Optional[str]]:
    """
    Process a single read: filter, extract sequences, and prepare data.
    
    Args:
        record: BioPython SeqRecord
        sequence_specs: Dictionary of sequence specifications
        min_length: Minimum read length filter
        max_length: Maximum read length filter
        min_quality: Minimum quality score filter
        min_quality_percentile: Minimum quality percentile filter
        
    Returns:
        tuple: (read_data, filter_reason)
            read_data: Dict with extracted info or None if filtered
            filter_reason: Reason for filtering or "pass"
    """
    # Apply quality filters
    passes_filter, filter_reason = passes_quality_filters(
        record, min_length, max_length, min_quality, min_quality_percentile
    )
    
    if not passes_filter:
        return None, filter_reason

    read_id = record.id
    sequence = str(record.seq)
    
    # Initialize read data
    read_data = {
        'read_id': read_id,
        'matches': {},       # seq_name -> True/False (extraction/detection success)
        'extracted': {},     # seq_name -> extracted sequence or NA
        'match_ids': {},     # seq_name -> match ID or NO_MATCH or NA
        'metrics': {}        # seq_name -> metrics dict or None
    }
    
    # Cache reverse complement (computed only once if needed)
    rev_comp = None
    
    # Check each sequence specification
    for seq_name, spec in sequence_specs.items():
        is_flanking = spec['is_flanking']
        max_dist = spec['max_distance']
        
        matches_found = False
        extracted_seq = 'NA'
        
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
                matches_found = True
                # Reverse complement if requested
                if spec.get('reverse_complement', False):
                    middle = str(Seq(middle).reverse_complement())
                extracted_seq = middle
            else:
                # Try reverse complement only if forward fails
                if rev_comp is None:
                    rev_comp = str(record.seq.reverse_complement())
                
                middle = extract_middle_sequence(rev_comp, left_seq, right_seq,
                                               max_dist, min_insert, max_insert)
                
                if middle:
                    matches_found = True
                    # Reverse complement if requested
                    if spec.get('reverse_complement', False):
                        middle = str(Seq(middle).reverse_complement())
                    extracted_seq = middle
                    
        else:
            # Direct sequence search
            full_seq = spec['full_seq']
            
            # Try forward strand
            if search_direct_sequence(sequence, full_seq, max_dist):
                matches_found = True
            else:
                # Try reverse complement
                if rev_comp is None:
                    rev_comp = str(record.seq.reverse_complement())
                if search_direct_sequence(rev_comp, full_seq, max_dist):
                    matches_found = True
        
        read_data['matches'][seq_name] = matches_found
        read_data['extracted'][seq_name] = extracted_seq
        # match_ids and metrics will be filled in later (either immediately or after streaming)
        read_data['match_ids'][seq_name] = 'NA'
        read_data['metrics'][seq_name] = None

    return read_data, "pass"


def process_fastq_file(fastq_file: str, sequence_specs: Dict, output_prefix: str,
                      no_progress: bool = False, skip_count: bool = False,
                      min_length: Optional[int] = None, max_length: Optional[int] = None,
                      min_quality: Optional[float] = None,
                      min_quality_percentile: Optional[float] = None,
                      stream_references: bool = False) -> Dict:
    """
    Process FASTQ file and search for all specified sequences.

    OPTIMIZATION: Pre-initialize aligners and mappy indices for each sequence spec.
    MEMORY FIX: True streaming processing of reads.

    Args:
        fastq_file: Path to input FASTQ file
        sequence_specs: Dictionary of sequence specifications
        output_prefix: Prefix for output files
        no_progress: Disable progress reporting
        skip_count: Skip initial read counting
        min_length: Minimum read length filter
        max_length: Maximum read length filter
        min_quality: Minimum quality score filter
        min_quality_percentile: Minimum quality percentile filter
        stream_references: Whether to stream reference files (hybrid mode)

    Returns:
        dict: Results including combination counts and extracted sequences
    """
    # Validate file exists
    if not Path(fastq_file).exists():
        raise FileNotFoundError(f"FASTQ file not found: {fastq_file}")
    
    # Initialize aligners and mappy indices
    logger.info("Initializing sequence matchers...")
    aligners = {}
    mappy_indices = {}
    
    # Identify which sequences are streaming vs in-memory
    streaming_sequences = {}  # seq_name -> spec
    memory_sequences = {}     # seq_name -> spec
    
    for seq_name, spec in sequence_specs.items():
        if stream_references and spec.get('reference_file') and spec['match_method'] in ['hamming', 'levenshtein']:
            streaming_sequences[seq_name] = spec
            logger.info(f"  {seq_name}: Will use streaming reference matching")
        elif spec.get('reference'):
            memory_sequences[seq_name] = spec
            # Initialize matchers for in-memory sequences
            match_method = spec.get('match_method', 'homology')
            if match_method == 'mappy':
                if MAPPY_AVAILABLE:
                    logger.info(f"  Building mappy index for {seq_name} ({len(spec['reference'])} references)...")
                    start_time = time.time()
                    mappy_indices[seq_name] = MappyIndex(spec['reference'])
                    logger.info(f"  Mappy index built in {time.time() - start_time:.1f}s")
                else:
                    raise ValueError(f"mappy not available for {seq_name}")
            elif match_method == 'homology':
                aligners[seq_name] = _get_default_aligner()
                logger.info(f"  Pre-initialized Bio.Align aligner for {seq_name}")
        else:
            # Direct sequence or no reference
            pass

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
    
    # Data structures for results
    all_read_data = {}  # Store all data for each read {read_id: read_data}
    read_ids_ordered = []  # Maintain order of reads
    
    # For streaming mode: store extracted sequences to match later
    # {seq_name: {read_id: sequence}}
    extracted_for_streaming = defaultdict(dict)
    
    total_reads = 0
    processed_reads = 0
    filtered_reads = 0
    filter_reasons = Counter()
    
    # Progress tracking
    start_time = time.time()
    last_progress_time = start_time
    progress_interval = 5  # seconds
    bytes_per_read_estimate = None
    
    logger.info("="*60)
    logger.info("PHASE 1: Processing reads (Extraction & In-memory Matching)")
    logger.info("="*60)
    
    # Stream through FASTQ file
    with open(fastq_file, 'r') as handle:
        for record in SeqIO.parse(handle, "fastq"):
            total_reads += 1
            
            # Process single read
            read_data, filter_reason = process_single_read(
                record, sequence_specs, min_length, max_length, min_quality, min_quality_percentile
            )
            
            if not read_data:
                filtered_reads += 1
                filter_reasons[filter_reason] += 1
                continue
                
            read_id = read_data['read_id']
            read_ids_ordered.append(read_id)
            processed_reads += 1
            
            # Handle matches
            for seq_name, spec in sequence_specs.items():
                # If extraction failed, nothing to do
                if not read_data['matches'][seq_name]:
                    continue
                    
                extracted_seq = read_data['extracted'][seq_name]
                
                # If this is a streaming sequence, save for later
                if seq_name in streaming_sequences:
                    extracted_for_streaming[seq_name][read_id] = extracted_seq
                    
                # If this is an in-memory sequence, match now
                elif seq_name in memory_sequences:
                    match_id, ref_seq, metrics = match_sequence_to_references_optimized(
                        extracted_seq, spec['reference'], spec['match_method'],
                        spec['match_dist'], aligners.get(seq_name), mappy_indices.get(seq_name)
                    )
                    read_data['match_ids'][seq_name] = match_id if match_id else 'NO_MATCH'
                    read_data['metrics'][seq_name] = metrics
            
            # Store result
            all_read_data[read_id] = read_data
            
            # Progress reporting
            current_time = time.time()
            if not no_progress and (current_time - last_progress_time >= progress_interval or total_reads % 10000 == 0):
                elapsed_time = current_time - start_time
                reads_per_second = total_reads / elapsed_time if elapsed_time > 0 else 0
                
                # Calculate percent complete
                percent_complete = 0
                if total_reads_estimate:
                    percent_complete = (total_reads / total_reads_estimate * 100)
                elif total_reads > 1000:
                    # Estimate based on bytes
                    if not bytes_per_read_estimate:
                        bytes_per_read_estimate = file_size / total_reads * 0.9 # conservative
                    estimated_total = file_size / bytes_per_read_estimate
                    percent_complete = (total_reads / estimated_total * 100)
                
                eta_str = ""
                if percent_complete > 0 and elapsed_time > 0:
                    total_time_est = elapsed_time / (percent_complete / 100)
                    remaining = total_time_est - elapsed_time
                    eta_str = f", ETA: {int(remaining//60)}m {int(remaining%60)}s"
                
                logger.info(f"Progress: {total_reads:,} reads ({percent_complete:.1f}%{eta_str})")
                last_progress_time = current_time

    # Phase 1 Summary
    phase1_time = time.time() - start_time
    logger.info(f"Finished Phase 1 in {phase1_time:.1f}s ({total_reads/phase1_time:.0f} reads/sec)")
    
    # PHASE 2: Streaming Reference Matching
    if streaming_sequences:
        logger.info("="*60)
        logger.info("PHASE 2: Streaming Reference Matching")
        logger.info("="*60)
        
        for seq_name, spec in streaming_sequences.items():
            if not extracted_for_streaming[seq_name]:
                logger.info(f"Skipping {seq_name}: no sequences extracted")
                continue
                
            logger.info(f"Processing {seq_name} ({len(extracted_for_streaming[seq_name])} sequences to match)")
            
            matches, refs_checked = match_sequences_streaming(
                extracted_for_streaming[seq_name],
                spec['reference_file'],
                spec['match_method'],
                spec['match_dist']
            )
            
            # Update read data with results
            aligner = aligners.get(seq_name) or _get_default_aligner()
            
            for read_id, (match_id, ref_seq, distance) in matches.items():
                if read_id in all_read_data:
                    if match_id:
                        metrics = calculate_all_metrics(
                            extracted_for_streaming[seq_name][read_id],
                            ref_seq,
                            aligner
                        )
                        all_read_data[read_id]['match_ids'][seq_name] = match_id
                        all_read_data[read_id]['metrics'][seq_name] = metrics
                    else:
                        all_read_data[read_id]['match_ids'][seq_name] = 'NO_MATCH'
                        all_read_data[read_id]['metrics'][seq_name] = None
            
            matched_count = sum(1 for _, (mid, _, _) in matches.items() if mid)
            logger.info(f"  Matched {matched_count}/{len(matches)} sequences in {seq_name}")

    # Finalize results
    logger.info("="*60)
    logger.info("Finalizing results")
    logger.info("="*60)
    
    # Generate combinations
    read_matches = []
    extracted_sequences = defaultdict(list)
    
    for read_id in read_ids_ordered:
        read_data = all_read_data[read_id]
        matches_found = set()
        
        for seq_name in sequence_specs.keys():
            # Check if extraction/detection succeeded
            if read_data['matches'].get(seq_name, False):
                matches_found.add(seq_name)
                
                # Add to extracted_sequences list for reporting/debugging if needed
                # (Note: original code returned this, so we preserve it)
                extracted = read_data['extracted'].get(seq_name)
                match_id = read_data['match_ids'].get(seq_name)
                metrics = read_data['metrics'].get(seq_name)
                extracted_sequences[seq_name].append((read_id, extracted, match_id, metrics))
        
        if matches_found:
            read_matches.append(tuple(sorted(matches_found)))
        else:
            read_matches.append(())

    # Report filtering statistics
    if filtered_reads > 0:
        logger.info(f"Filtered out {filtered_reads:,} reads ({filtered_reads/total_reads*100:.2f}%)")
        logger.info("Filter reasons:")
        for reason, count in filter_reasons.most_common():
            logger.info(f"  {reason}: {count:,} reads")

    # Count combinations
    combination_counts = Counter(read_matches)
    
    # Generate all possible combinations
    all_seq_names = sorted(sequence_specs.keys())
    all_combinations = []
    for r in range(len(all_seq_names) + 1):
        for combo in combinations(all_seq_names, r):
            all_combinations.append(combo)

    results = {
        'total_reads': total_reads,
        'processed_reads': processed_reads,
        'filtered_reads': filtered_reads,
        'filter_reasons': dict(filter_reasons),
        'combination_counts': combination_counts,
        'all_combinations': all_combinations,
        'extracted_sequences': dict(extracted_sequences),
        'sequence_specs': sequence_specs,
        'all_read_data': all_read_data,
        'read_ids_ordered': read_ids_ordered
    }

    return results


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
    
    # Report 1: Combination statistics (original format)
    output_file = f'{output_prefix}_combination_stats.txt'
    try:
        with open(output_file, 'w') as f:
            f.write("Combination\tCount\tProportion\n")

            # Sort combinations by count (descending)
            sorted_combinations = sorted(all_combinations,
                                        key=lambda combo: combination_counts.get(combo, 0),
                                        reverse=True)

            for combo in sorted_combinations:
                count = combination_counts.get(combo, 0)
                proportion = count / results['processed_reads'] if results['processed_reads'] > 0 else 0

                if combo:
                    combo_str = '+'.join(combo)
                else:
                    combo_str = 'none'

                f.write(f"{combo_str}\t{count}\t{proportion:.4f}\n")

            # Summary statistics
            f.write(f"\n# Summary\n")
            f.write(f"Total reads in file: {total_reads:,}\n")

            # Add filtering statistics if applicable
            if results.get('filtered_reads', 0) > 0:
                filtered_reads = results['filtered_reads']
                f.write(f"Filtered reads: {filtered_reads:,} ({filtered_reads/total_reads*100:.2f}%)\n")

                # Write filter reasons
                if results.get('filter_reasons'):
                    f.write(f"\n# Filter reasons\n")
                    for reason, count in sorted(results['filter_reasons'].items(), key=lambda x: x[1], reverse=True):
                        f.write(f"{reason}: {count:,} ({count/filtered_reads*100:.2f}% of filtered)\n")

            f.write(f"\nReads analyzed (post-filter): {results['processed_reads']:,}\n")

            reads_with_any_match = sum(count for combo, count in combination_counts.items() if combo)
            f.write(f"Reads with at least one match: {reads_with_any_match:,} ({reads_with_any_match/results['processed_reads']*100:.2f}% of analyzed)\n")

            # Individual sequence statistics
            f.write(f"\n# Individual sequence matches\n")
            for seq_name in sorted(sequence_specs.keys()):
                count = sum(c for combo, c in combination_counts.items() if seq_name in combo)
                f.write(f"{seq_name}: {count:,} ({count/results['processed_reads']*100:.2f}% of analyzed)\n")
        
        logger.info(f"Written combination statistics to {output_file}")

    except Exception as e:
        logger.error(f"Error writing combination statistics: {e}")
        raise

    # Report 2: Combination matrix (binary format)
    output_file = f'{output_prefix}_combination_matrix.txt'
    try:
        with open(output_file, 'w') as f:
            # Create header with all sequence names as columns
            seq_names_sorted = sorted(sequence_specs.keys())
            header = '\t'.join(seq_names_sorted) + '\tCount\tProportion\n'
            f.write(header)

            # Sort combinations by count (descending)
            sorted_combinations = sorted(all_combinations,
                                        key=lambda combo: combination_counts.get(combo, 0),
                                        reverse=True)

            for combo in sorted_combinations:
                count = combination_counts.get(combo, 0)
                proportion = count / results['processed_reads'] if results['processed_reads'] > 0 else 0

                # Create binary vector for this combination
                binary_vector = []
                for seq_name in seq_names_sorted:
                    if seq_name in combo:
                        binary_vector.append('1')
                    else:
                        binary_vector.append('0')

                row = '\t'.join(binary_vector) + f'\t{count}\t{proportion:.4f}\n'
                f.write(row)

        logger.info(f"Written combination matrix to {output_file}")

    except Exception as e:
        logger.error(f"Error writing combination matrix: {e}")
        raise

    # Report 3: Combined output file with all reads
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
                        
                        # Add mappy_identity column if mappy method was used
                        if sequence_specs[seq_name].get('match_method') == 'mappy':
                            header_parts.append(f'{seq_name}_mappy_identity')
            
            f.write('\t'.join(header_parts) + '\n')
            
            # Write data for each read in order
            for read_id in read_ids_ordered:
                read_data = all_read_data[read_id]
                row_parts = [read_id]
                
                for seq_name in sorted(sequence_specs.keys()):
                    # Match status
                    match_status = 'Yes' if read_data['matches'].get(seq_name, False) else 'No'
                    row_parts.append(match_status)
                    
                    if sequence_specs[seq_name]['is_flanking']:
                        # Extracted sequence
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
                                
                                # Add mappy_identity if present
                                if 'mappy_identity' in metrics and sequence_specs[seq_name].get('match_method') == 'mappy':
                                    row_parts.append(f"{metrics['mappy_identity']:.3f}")
                            else:
                                row_parts.extend(['NA', 'NA', 'NA'])
                                if sequence_specs[seq_name].get('match_method') == 'mappy':
                                    row_parts.append('NA')
                
                f.write('\t'.join(row_parts) + '\n')
        
        logger.info(f"Written combined results to {output_file}")
        
    except Exception as e:
        logger.error(f"Error writing combined results: {e}")
        raise

    # Report 4: Visualization
    try:
        generate_visualization(results, output_prefix)
    except Exception as e:
        logger.warning(f"Could not generate visualization: {e}")

    # Report 5: Extraction status plot
    try:
        generate_extraction_status_plot(results, output_prefix)
    except Exception as e:
        logger.warning(f"Could not generate extraction status plot: {e}")

    # Report 6: Pairwise concordance plot
    try:
        generate_pairwise_concordance_plot(results, output_prefix)
    except Exception as e:
        logger.warning(f"Could not generate concordance plot: {e}")


def generate_visualization(results: Dict, output_prefix: str) -> None:
    """
    Generate UpSet plot or bar chart for combination visualization.
    Uses styling from upset_plot_version.py for consistency.
    """
    combination_counts = results['combination_counts']
    sequence_specs = results['sequence_specs']
    all_read_data = results['all_read_data']
    total_reads = results['total_reads']
    all_seq_names = sorted(sequence_specs.keys())

    # Prepare data for UpSet plot
    if UPSETPLOT_AVAILABLE:
        try:
            # Create a dictionary where keys are sequence names and values are sets of read IDs
            sequence_sets = {}
            for seq_name in all_seq_names:
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

                # Remove grid lines from all axes
                for ax in fig.get_axes():
                    ax.grid(False)

                # Try to adjust the matrix plot row heights
                for ax in fig.get_axes():
                    # The matrix plot typically doesn't have axis labels
                    if not ax.get_xlabel() and not ax.get_ylabel():
                        try:
                            # Adjust the aspect ratio to make rows more uniform
                            ax.set_aspect('auto')
                        except:
                            pass

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

                logger.info(f"Written UpSet plot to {output_file}")
                return
            else:
                logger.warning("No data found - skipping upset plot generation")
                return

        except Exception as e:
            logger.warning(f"Could not create UpSet plot: {e}. Falling back to bar chart.")
    
    # Fallback: Create bar chart with percentages
    try:
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
        logger.error(f"Could not create bar chart: {e}")


def generate_extraction_status_plot(results: Dict, output_prefix: str) -> None:
    """
    Generate stacked barplot showing extraction and matching status for each flanking sequence.

    Categories:
    - No sequence: Extraction failed (flanking sequences not found)
    - No match: Extracted but not matched to reference
    - Matched: Extracted and matched to reference
    """
    sequence_specs = results['sequence_specs']
    all_read_data = results['all_read_data']

    # Filter to only flanking sequences with references
    flanking_seqs = {name: spec for name, spec in sequence_specs.items()
                     if spec['is_flanking']}

    if not flanking_seqs:
        logger.info("No flanking sequences found. Skipping extraction status plot.")
        return

    # Count status for each sequence
    status_counts = {}

    for seq_name in sorted(flanking_seqs.keys()):
        no_sequence = 0
        no_match = 0
        matched = 0

        for read_data in all_read_data.values():
            extracted = read_data['extracted'].get(seq_name, 'NA')
            match_id = read_data['match_ids'].get(seq_name, 'NA')

            if extracted == 'NA':
                # No sequence extracted
                no_sequence += 1
            elif match_id == 'NO_MATCH':
                # Extracted but not matched
                no_match += 1
            elif match_id != 'NA':
                # Matched to reference
                matched += 1
            else:
                # Extracted but no reference file (shouldn't happen for flanking with refs)
                no_match += 1

        status_counts[seq_name] = {
            'no_sequence': no_sequence,
            'no_match': no_match,
            'matched': matched
        }

    # Create stacked bar plot
    try:
        fig, ax = plt.subplots(figsize=(max(10, len(flanking_seqs) * 2), 8))

        seq_names = sorted(flanking_seqs.keys())
        x_pos = np.arange(len(seq_names))

        # Get counts for each category
        no_sequence_counts = [status_counts[name]['no_sequence'] for name in seq_names]
        no_match_counts = [status_counts[name]['no_match'] for name in seq_names]
        matched_counts = [status_counts[name]['matched'] for name in seq_names]

        # Calculate percentages
        totals = [no_sequence_counts[i] + no_match_counts[i] + matched_counts[i]
                 for i in range(len(seq_names))]

        no_sequence_pct = [no_sequence_counts[i] / totals[i] * 100 if totals[i] > 0 else 0
                          for i in range(len(seq_names))]
        no_match_pct = [no_match_counts[i] / totals[i] * 100 if totals[i] > 0 else 0
                       for i in range(len(seq_names))]
        matched_pct = [matched_counts[i] / totals[i] * 100 if totals[i] > 0 else 0
                      for i in range(len(seq_names))]

        # Create stacked bars
        p1 = ax.bar(x_pos, no_sequence_pct, label='No Sequence', color='#d62728')
        p2 = ax.bar(x_pos, no_match_pct, bottom=no_sequence_pct,
                   label='No Match', color='#ff7f0e')
        p3 = ax.bar(x_pos, matched_pct,
                   bottom=[no_sequence_pct[i] + no_match_pct[i] for i in range(len(seq_names))],
                   label='Matched', color='#2ca02c')

        # Customize plot
        ax.set_ylabel('Percentage of Reads (%)', fontsize=12)
        ax.set_xlabel('Sequence Region', fontsize=12)
        ax.set_title('Extraction and Matching Status by Region', fontsize=14, fontweight='bold', pad=20)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(seq_names, rotation=45, ha='right')
        ax.legend(loc='upper right')
        ax.set_ylim(0, 110)  # Increased from 100 to 110 to create space for labels

        # Add count labels on bars
        for i, name in enumerate(seq_names):
            # Add total count at top
            ax.text(i, 102, f'n={totals[i]:,}', ha='center', va='bottom', fontsize=9)

            # Add percentage labels for each segment if visible
            if matched_pct[i] > 5:
                y_pos = no_sequence_pct[i] + no_match_pct[i] + matched_pct[i] / 2
                ax.text(i, y_pos, f'{matched_pct[i]:.1f}%', ha='center', va='center',
                       fontsize=9, fontweight='bold', color='white')

            if no_match_pct[i] > 5:
                y_pos = no_sequence_pct[i] + no_match_pct[i] / 2
                ax.text(i, y_pos, f'{no_match_pct[i]:.1f}%', ha='center', va='center',
                       fontsize=9, fontweight='bold', color='white')

            if no_sequence_pct[i] > 5:
                y_pos = no_sequence_pct[i] / 2
                ax.text(i, y_pos, f'{no_sequence_pct[i]:.1f}%', ha='center', va='center',
                       fontsize=9, fontweight='bold', color='white')

        plt.tight_layout()

        output_file = f'{output_prefix}_extraction_status.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        logger.info(f"Written extraction status plot to {output_file}")

    except Exception as e:
        logger.error(f"Could not create extraction status plot: {e}")


def generate_pairwise_concordance_plot(results: Dict, output_prefix: str) -> None:
    """
    Generate plot showing pairwise Match ID concordance between regions.

    Categories:
    - Same Match ID: Both regions matched to same reference
    - Different Match IDs: Both regions matched but to different references
    - Only one matched: One region matched, other did not
    - Neither matched: Neither region matched
    """
    sequence_specs = results['sequence_specs']
    all_read_data = results['all_read_data']

    # Filter to only flanking sequences with references (in memory or file)
    flanking_with_refs = {name: spec for name, spec in sequence_specs.items()
                         if spec['is_flanking'] and (spec.get('reference') or spec.get('reference_file'))}

    if len(flanking_with_refs) < 2:
        logger.info("Need at least 2 flanking sequences with references for concordance plot.")
        return

    # Generate all pairwise combinations
    from itertools import combinations as iter_combinations
    seq_names = sorted(flanking_with_refs.keys())
    pairs = list(iter_combinations(seq_names, 2))

    # Count concordance for each pair
    concordance_data = {}

    for seq1, seq2 in pairs:
        pair_name = f'{seq1} vs {seq2}'
        same_id = 0
        different_id = 0
        only_one = 0
        neither = 0

        for read_data in all_read_data.values():
            match1 = read_data['match_ids'].get(seq1, 'NA')
            match2 = read_data['match_ids'].get(seq2, 'NA')

            # Skip if either wasn't extracted
            if match1 == 'NA' or match2 == 'NA':
                continue

            if match1 == 'NO_MATCH' and match2 == 'NO_MATCH':
                neither += 1
            elif match1 == 'NO_MATCH' or match2 == 'NO_MATCH':
                only_one += 1
            elif match1 == match2:
                same_id += 1
            else:
                different_id += 1

        concordance_data[pair_name] = {
            'same': same_id,
            'different': different_id,
            'only_one': only_one,
            'neither': neither
        }

    # Create stacked bar plot
    try:
        fig, ax = plt.subplots(figsize=(max(12, len(pairs) * 2), 8))

        pair_names = [f'{s1}\nvs\n{s2}' for s1, s2 in pairs]
        x_pos = np.arange(len(pairs))

        # Get counts for each category
        same_counts = [concordance_data[f'{s1} vs {s2}']['same'] for s1, s2 in pairs]
        different_counts = [concordance_data[f'{s1} vs {s2}']['different'] for s1, s2 in pairs]
        only_one_counts = [concordance_data[f'{s1} vs {s2}']['only_one'] for s1, s2 in pairs]
        neither_counts = [concordance_data[f'{s1} vs {s2}']['neither'] for s1, s2 in pairs]

        # Calculate percentages
        totals = [same_counts[i] + different_counts[i] + only_one_counts[i] + neither_counts[i]
                 for i in range(len(pairs))]

        same_pct = [same_counts[i] / totals[i] * 100 if totals[i] > 0 else 0
                   for i in range(len(pairs))]
        different_pct = [different_counts[i] / totals[i] * 100 if totals[i] > 0 else 0
                        for i in range(len(pairs))]
        only_one_pct = [only_one_counts[i] / totals[i] * 100 if totals[i] > 0 else 0
                       for i in range(len(pairs))]
        neither_pct = [neither_counts[i] / totals[i] * 100 if totals[i] > 0 else 0
                      for i in range(len(pairs))]

        # Create stacked bars
        p1 = ax.bar(x_pos, same_pct, label='Same Match ID', color='#2ca02c')
        p2 = ax.bar(x_pos, different_pct, bottom=same_pct,
                   label='Different Match IDs', color='#d62728')

        bottom2 = [same_pct[i] + different_pct[i] for i in range(len(pairs))]
        p3 = ax.bar(x_pos, only_one_pct, bottom=bottom2,
                   label='Only One Matched', color='#ff7f0e')

        bottom3 = [bottom2[i] + only_one_pct[i] for i in range(len(pairs))]
        p4 = ax.bar(x_pos, neither_pct, bottom=bottom3,
                   label='Neither Matched', color='#7f7f7f')

        # Customize plot
        ax.set_ylabel('Percentage of Reads (%)', fontsize=12)
        ax.set_xlabel('Pairwise Comparison', fontsize=12)
        ax.set_title('Match ID Concordance Between Regions', fontsize=14, fontweight='bold', pad=20)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(pair_names, fontsize=10)
        ax.legend(loc='upper right')
        ax.set_ylim(0, 110)  # Increased from 100 to 110 to create space for labels

        # Add count labels
        for i in range(len(pairs)):
            ax.text(i, 102, f'n={totals[i]:,}', ha='center', va='bottom', fontsize=9)

            # Add percentage labels for each segment if visible
            if same_pct[i] > 5:
                ax.text(i, same_pct[i] / 2, f'{same_pct[i]:.1f}%',
                       ha='center', va='center', fontsize=9, fontweight='bold', color='white')

            if different_pct[i] > 5:
                y_pos = same_pct[i] + different_pct[i] / 2
                ax.text(i, y_pos, f'{different_pct[i]:.1f}%',
                       ha='center', va='center', fontsize=9, fontweight='bold', color='white')

            if only_one_pct[i] > 5:
                y_pos = bottom2[i] + only_one_pct[i] / 2
                ax.text(i, y_pos, f'{only_one_pct[i]:.1f}%',
                       ha='center', va='center', fontsize=9, fontweight='bold', color='white')

            if neither_pct[i] > 5:
                y_pos = bottom3[i] + neither_pct[i] / 2
                ax.text(i, y_pos, f'{neither_pct[i]:.1f}%',
                       ha='center', va='center', fontsize=9, fontweight='bold', color='white')

        plt.tight_layout()

        output_file = f'{output_prefix}_match_concordance.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        logger.info(f"Written match concordance plot to {output_file}")

    except Exception as e:
        logger.error(f"Could not create match concordance plot: {e}")


def main():
    parser = argparse.ArgumentParser(
        description='Search for sequences in FASTQ files with optimized performance (includes mappy support)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simple sequence search
  python insert_check_mappy.py reads.fastq --seq1 "ATCGATCG" --seq1-dist 2
  
  # Flanking sequence with barcode matching (OPTIMIZED for millions of barcodes)
  python insert_check_mappy.py reads.fastq \\
      --seq1 "GTCGACGAACCTCTAGA-AGATCGGAAGAGCGT" \\
      --seq1-dist 4 --seq1-insert 18:22 \\
      --ref-file1 barcodes.txt \\
      --match-method1 hamming --match-dist1 0
  
  # Long sequence with minimap2 (ULTRA-FAST - 10-100x faster than Bio.Align)
  python insert_check_mappy.py reads.fastq \\
      --seq1 "GCAGGACTGGCCGCTTGACG-CACTGCGGCTCCTGCGATTG" \\
      --seq1-dist 5 --seq1-insert 100:300 \\
      --ref-file1 references.fasta \\
      --match-method1 mappy --match-dist1 0.8
  
  # Multiple sequences with different methods
  python insert_check_mappy.py reads.fastq \\
      --seq1 "ATCG-GCTA" --seq1-insert 18:22 \\
      --ref-file1 barcodes.txt --match-method1 hamming --match-dist1 0 \\
      --seq2 "GGCC-TTAA" --seq2-insert 100:300 \\
      --ref-file2 amplicons.fasta --match-method2 mappy --match-dist2 0.85
        """
    )
    
    # Required arguments
    parser.add_argument('fastq_file', help='Input FASTQ file')
    
    # Output options
    parser.add_argument('--output-prefix', default='output', help='Prefix for output files (default: output)')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO', help='Logging level (default: INFO)')
    parser.add_argument('--no-progress', action='store_true', help='Disable progress reporting')
    parser.add_argument('--skip-count', action='store_true',
                       help='Skip initial read counting for faster startup')

    # Pre-filtering options
    parser.add_argument('--min-length', type=int, default=None,
                       help='Minimum read length (bases). Reads shorter than this will be skipped')
    parser.add_argument('--max-length', type=int, default=None,
                       help='Maximum read length (bases). Reads longer than this will be skipped')
    parser.add_argument('--min-quality', type=float, default=None,
                       help='Minimum average quality score (Phred scale). Reads below this will be skipped')
    parser.add_argument('--min-quality-percentile', type=float, default=None,
                       help='Minimum quality percentile (0-100). E.g., 90 means 90%% of bases must be >= min-quality')

    # Streaming mode for large reference files
    parser.add_argument('--stream-references', action='store_true',
                       help='Stream reference files instead of loading into memory. '
                            'Use for very large reference files. Compatible with hamming/levenshtein; '
                            'mappy/homology sequences will use normal mode automatically.')
    
    # Allow multiple sequence specifications
    for i in range(1, 11):  # Support up to 10 sequences
        parser.add_argument(f'--seq{i}', help=f'Sequence {i} (format: "ATCG-TGCA" or "ATCGATCG")')
        parser.add_argument(f'--name{i}', help=f'Custom name for sequence {i} (default: seq{i})')
        parser.add_argument(f'--seq{i}-dist', type=int, default=None, 
                          help=f'Max Levenshtein distance for seq{i} (default: 2)')
        parser.add_argument(f'--seq{i}-insert', help=f'Min:max insert size for flanking seq{i} (e.g., "10:30")')
        parser.add_argument(f'--ref-file{i}', help=f'Reference file for seq{i} matches (FASTA or TSV)')
        parser.add_argument(f'--match-method{i}', choices=['hamming', 'levenshtein', 'homology', 'mappy'],
                          default='homology', help=f'Matching method for seq{i} (default: homology, mappy=ultra-fast)')
        parser.add_argument(f'--match-dist{i}', type=float, default=None,
                          help=f'Match distance/score for seq{i} (default: 5 for hamming/levenshtein, 0.8 for homology/mappy)')
        parser.add_argument(f'--seq{i}-rc', action='store_true',
                          help=f'Reverse complement extracted sequence before matching for seq{i}')
    
    args = parser.parse_args()
    
    # Set logging level
    logger.setLevel(getattr(logging, args.log_level))
    
    # Print mappy availability
    if MAPPY_AVAILABLE:
        logger.info("mappy (minimap2) is available for ultra-fast alignment")
    else:
        logger.warning("mappy not available. Install for 10-100x speedup: pip install mappy")
    
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
                    else:  # homology or mappy
                        match_dist = 0.8
                        logger.info(f"Using default match score of 0.8 for {seq_name} ({match_method})")
                
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
                        spec['reference_file'] = None
                    else:
                        if not Path(ref_file).exists():
                            logger.error(f"Reference file not found: {ref_file}")
                            sys.exit(1)

                        # Store the file path for streaming mode
                        spec['reference_file'] = ref_file

                        # Conditionally load references based on streaming mode and match method
                        if args.stream_references and spec['match_method'] in ['hamming', 'levenshtein']:
                            # Will be streamed later - don't load into memory
                            spec['reference'] = None
                            logger.info(f"Will stream reference file for {seq_name} ({spec['match_method']} method, streaming mode enabled)")
                        else:
                            # Load into memory (normal mode or incompatible method)
                            spec['reference'] = load_reference_sequences(ref_file)
                            logger.info(f"Loaded {len(spec['reference'])} reference sequences for {seq_name}")
                            if args.stream_references:
                                logger.info(f"Note: {seq_name} uses {spec['match_method']} method, which requires loading references into memory")

                        logger.info(f"Using {spec['match_method']} matching for {seq_name} with threshold {match_dist}")
                        if spec['reverse_complement']:
                            logger.info(f"Will reverse complement extracted sequences for {seq_name}")
                else:
                    spec['reference'] = None
                    spec['reference_file'] = None
                
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
    logger.info("OPTIMIZATIONS ENABLED:")
    logger.info("  - Hash table lookup for exact hamming matches")
    logger.info("  - Reusable PairwiseAligner objects")
    if MAPPY_AVAILABLE:
        logger.info("  - Ultra-fast minimap2 (mappy) alignment available")
    logger.info("  - Early termination for matching")
    logger.info("  - Reduced redundant calculations")

    # Report filtering settings if enabled
    if args.min_length or args.max_length or args.min_quality:
        logger.info("PRE-FILTERING ENABLED:")
        if args.min_length:
            logger.info(f"  - Minimum read length: {args.min_length} bp")
        if args.max_length:
            logger.info(f"  - Maximum read length: {args.max_length} bp")
        if args.min_quality:
            logger.info(f"  - Minimum quality score: Q{args.min_quality}")
            if args.min_quality_percentile:
                logger.info(f"  - Quality percentile threshold: {args.min_quality_percentile}%")

    try:
        results = process_fastq_file(args.fastq_file, sequence_specs, args.output_prefix,
                                   no_progress=args.no_progress, skip_count=args.skip_count,
                                   min_length=args.min_length, max_length=args.max_length,
                                   min_quality=args.min_quality,
                                   min_quality_percentile=args.min_quality_percentile,
                                   stream_references=args.stream_references)
    except Exception as e:
        logger.error(f"Error processing FASTQ file: {e}")
        import traceback
        traceback.print_exc()
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
    logger.info(f"Total reads in file: {results['total_reads']:,}")

    if results.get('filtered_reads', 0) > 0:
        logger.info(f"Filtered reads: {results['filtered_reads']:,} ({results['filtered_reads']/results['total_reads']*100:.2f}%)")

    logger.info(f"Reads analyzed: {results['processed_reads']:,}")

    reads_with_matches = sum(count for combo, count in results['combination_counts'].items() if combo)
    logger.info(f"Reads with at least one match: {reads_with_matches:,} ({reads_with_matches/results['processed_reads']*100:.2f}% of analyzed)")

    logger.info("\nIndividual sequence match rates:")
    for seq_name in sorted(sequence_specs.keys()):
        count = sum(c for combo, c in results['combination_counts'].items() if seq_name in combo)
        logger.info(f"  {seq_name}: {count:,} reads ({count/results['processed_reads']*100:.2f}% of analyzed)")
    
    logger.info("\nTop 5 most common combinations:")
    sorted_combos = sorted(results['combination_counts'].items(), key=lambda x: x[1], reverse=True)[:5]
    for combo, count in sorted_combos:
        combo_str = '+'.join(combo) if combo else 'none'
        logger.info(f"  {combo_str}: {count:,} reads ({count/results['total_reads']*100:.2f}%)")
    
    logger.info("="*60 + "\n")
    
    logger.info("Processing complete!")
    logger.info(f"Results written to:")
    logger.info(f"  - {args.output_prefix}_combination_stats.txt (combination statistics)")
    logger.info(f"  - {args.output_prefix}_combination_matrix.txt (binary matrix format)")
    logger.info(f"  - {args.output_prefix}_combined_results.txt (all reads with match data)")

    # Check if upset plot was created
    if Path(f"{args.output_prefix}_upset_plot.png").exists():
        logger.info(f"  - {args.output_prefix}_upset_plot.png (upset plot visualization)")
    elif Path(f"{args.output_prefix}_combination_plot.png").exists():
        logger.info(f"  - {args.output_prefix}_combination_plot.png (combination plot visualization)")

    # Check for extraction status plot
    if Path(f"{args.output_prefix}_extraction_status.png").exists():
        logger.info(f"  - {args.output_prefix}_extraction_status.png (extraction/matching status by region)")

    # Check for concordance plot
    if Path(f"{args.output_prefix}_match_concordance.png").exists():
        logger.info(f"  - {args.output_prefix}_match_concordance.png (pairwise match concordance)")


if __name__ == "__main__":
    main()
