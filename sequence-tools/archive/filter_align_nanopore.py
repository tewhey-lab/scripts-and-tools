#!/usr/bin/env python3
"""
Nanopore Library Alignment Tool
Specialized for aligning reads to plasmid libraries where N regions represent variable inserts
Handles the critical issue of NOT penalizing alignments to N regions
"""

import argparse
import sys
import os
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import logging
from tqdm import tqdm
import json
import re

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

try:
    import mappy as mp
    HAS_MINIMAP2 = True
except ImportError:
    HAS_MINIMAP2 = False
    logger.warning("mappy not installed. Install with: pip install mappy")


class LibraryAligner:
    """
    Specialized aligner for library screening that properly handles N regions
    """
    
    def __init__(self, reference_fasta, min_length=0, min_quality=0, 
                 flank_size=50, min_flank_identity=0.9, circular=True):
        """
        Args:
            reference_fasta: FASTA file with reference sequences
            min_length: Minimum read length
            min_quality: Minimum mean quality score  
            flank_size: Size of flanking regions (used for insert extraction context)
            min_flank_identity: Minimum identity required (not used in new strategy but kept for API compat)
            circular: Whether plasmids are circular
        """
        self.reference_fasta = reference_fasta
        self.min_length = min_length
        self.min_quality = min_quality
        self.flank_size = flank_size
        self.min_flank_identity = min_flank_identity
        self.circular = circular
        
        self.references = {}
        self.backbones = {} # Ref name -> backbone sequence (no Ns)
        self.n_regions = {} # Ref name -> list of (start, end)
        self.backbone_map = {} # Ref name -> list mapping backbone indices to original indices
        self.aligners = {} # Ref name -> mappy.Aligner
        
        self.filtered_reads = []
        self.assignments = defaultdict(list)
        self.stats = defaultdict(int)
        
        self.load_references()
    
    def load_references(self):
        """Load references, identify N regions, and prepare aligners"""
        logger.info(f"Loading references from {self.reference_fasta}")
        
        for record in SeqIO.parse(self.reference_fasta, "fasta"):
            seq_str = str(record.seq).upper()
            self.references[record.id] = seq_str
            
            # Find N regions
            n_regions = []
            for match in re.finditer(r'N+', seq_str):
                n_regions.append((match.start(), match.end()))
            
            self.n_regions[record.id] = n_regions
            
            # Create backbone (remove Ns) and mapping
            backbone = []
            mapping = [] # backbone_index -> original_index
            
            current_orig_idx = 0
            for start, end in n_regions:
                # Add sequence before N region
                chunk = seq_str[current_orig_idx:start]
                backbone.append(chunk)
                mapping.extend(range(current_orig_idx, start))
                
                current_orig_idx = end
            
            # Add remaining sequence
            chunk = seq_str[current_orig_idx:]
            backbone.append(chunk)
            mapping.extend(range(current_orig_idx, len(seq_str)))
            
            backbone_seq = "".join(backbone)
            self.backbones[record.id] = backbone_seq
            self.backbone_map[record.id] = mapping
            
            if n_regions:
                logger.info(f"  {record.id}: {len(seq_str)} bp with {len(n_regions)} N region(s)")
                logger.info(f"    Backbone length: {len(backbone_seq)} bp")
            else:
                logger.info(f"  {record.id}: {len(seq_str)} bp (no N regions)")
            
            # Initialize aligner for this backbone
            if HAS_MINIMAP2:
                target_seq = backbone_seq
                if self.circular:
                    # Double the backbone to handle circularity
                    target_seq = backbone_seq + backbone_seq
                
                try:
                    self.aligners[record.id] = mp.Aligner(seq=target_seq, preset="map-ont", best_n=1)
                except Exception as e:
                    logger.error(f"Failed to initialize aligner for {record.id}: {e}")

    def calculate_quality(self, record):
        """Calculate mean quality score"""
        if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
            return np.mean(record.letter_annotations['phred_quality'])
        return 20  # Default
    
    def filter_read(self, record):
        """Check if read passes basic filters"""
        if len(record.seq) < self.min_length:
            return False
        if self.calculate_quality(record) < self.min_quality:
            return False
        return True
    
    def align_to_reference_smart(self, read_seq, ref_name, ref_seq):
        """
        Align read to backbone and detect inserts corresponding to N regions
        """
        if ref_name not in self.aligners:
            return None
            
        aligner = self.aligners[ref_name]
        read_str = str(read_seq)
        
        try:
            # Map read to backbone
            alignments = list(aligner.map(read_str))
            if not alignments:
                return None
            
            best = alignments[0]
            
            # Check for large insertions in the read that map to N regions
            # We need to map the backbone coordinates back to original coordinates
            # to see if the insertion happens at an N-region boundary.
            
            # Parse CIGAR to find insertions
            # best.cigar is list of (len, op) where op: 0=M, 1=I, 2=D
            
            current_ref_pos = best.r_st
            current_read_pos = best.q_st
            
            max_score = best.mlen
            detected_insert = None
            spans_n_region = False
            
            # If circular, normalize ref_pos
            backbone_len = len(self.backbones[ref_name])
            
            # Iterate through CIGAR
            for length, op in best.cigar:
                if op == 0: # Match/Mismatch
                    current_ref_pos += length
                    current_read_pos += length
                elif op == 1: # Insertion in read (potential N region insert)
                    # Check if this insertion aligns with an N region in the original sequence
                    # We check if the current_ref_pos corresponds to a break in the backbone map
                    
                    # Normalize pos if circular
                    norm_pos = current_ref_pos % backbone_len
                    
                    # Check if norm_pos is at the start of an N region
                    # In our backbone map, if we are at index i, the original index is map[i].
                    # If map[i] and map[i-1] are not contiguous, there was a gap (N region).
                    
                    is_n_junction = False
                    n_region_idx = -1
                    
                    # Check if we are at a junction
                    # We need to be careful with boundary conditions
                    if 0 < norm_pos < backbone_len:
                        orig_before = self.backbone_map[ref_name][norm_pos - 1]
                        orig_after = self.backbone_map[ref_name][norm_pos]
                        
                        if orig_after - orig_before > 1:
                            is_n_junction = True
                            # Find which N region this is
                            for i, (n_start, n_end) in enumerate(self.n_regions[ref_name]):
                                if n_start == orig_before + 1:
                                    n_region_idx = i
                                    break
                    
                    if is_n_junction:
                        insert_seq = read_str[current_read_pos : current_read_pos + length]
                        detected_insert = insert_seq
                        spans_n_region = True
                        
                        # Boost score: add back the insertion penalty (usually large for long inserts)
                        # and add a bonus for matching the N-region
                        # Minimap2 gap open/ext costs are roughly 4+2*k or similar.
                        # We'll just add length * 2 as a bonus, similar to previous logic
                        max_score += length * 2
                        
                    current_read_pos += length
                    
                elif op == 2: # Deletion from ref
                    current_ref_pos += length
            
            return {
                'ref_name': ref_name,
                'ref_start': best.r_st % backbone_len,
                'ref_end': best.r_en % backbone_len,
                'query_start': best.q_st,
                'query_end': best.q_en,
                'strand': '+' if best.strand == 1 else '-',
                'identity': best.mlen / best.blen if best.blen > 0 else 0,
                'score': max_score,
                'mapq': best.mapq,
                'has_n_region': bool(self.n_regions[ref_name]),
                'spans_n_region': spans_n_region,
                'insert_sequence': detected_insert
            }
            
        except Exception as e:
            logger.debug(f"Alignment failed: {e}")
            return None
    
    def assign_read(self, record):
        """
        Assign read to best reference using N-aware scoring
        """
        read_seq = str(record.seq)
        alignments = []
        
        # Align to each reference
        for ref_name, ref_seq in self.references.items():
            alignment = self.align_to_reference_smart(read_seq, ref_name, ref_seq)
            if alignment:
                alignments.append(alignment)
        
        if not alignments:
            return None
        
        # Sort by score (which has been adjusted for N regions)
        alignments.sort(key=lambda x: x['score'], reverse=True)
        best = alignments[0]
        
        # Check if there's a clear winner
        if len(alignments) > 1:
            second_best = alignments[1]
            score_ratio = best['score'] / (second_best['score'] + 1)
            
            if score_ratio < 1.1:  # Ambiguous assignment
                logger.debug(f"Ambiguous assignment for {record.id}: "
                           f"{best['ref_name']} ({best['score']}) vs "
                           f"{second_best['ref_name']} ({second_best['score']})")
                self.stats['ambiguous'] += 1
        
        return best
    
    def process_reads(self, input_fastq):
        """Process all reads in the FASTQ file"""
        logger.info(f"Processing {input_fastq}")
        
        for record in tqdm(SeqIO.parse(input_fastq, "fastq"), desc="Processing reads"):
            self.stats['total'] += 1
            
            # Basic filtering
            if not self.filter_read(record):
                self.stats['filtered_out'] += 1
                continue
            
            self.stats['passed_filter'] += 1
            self.filtered_reads.append(record)
            
            # Assign to reference
            assignment = self.assign_read(record)
            
            if assignment:
                self.stats['assigned'] += 1
                self.stats[f'assigned_to_{assignment["ref_name"]}'] += 1
                
                if assignment.get('spans_n_region'):
                    self.stats[f'spans_n_in_{assignment["ref_name"]}'] += 1
                
                # Store assignment
                self.assignments[assignment['ref_name']].append({
                    'read_id': record.id,
                    'record': record,
                    'alignment': assignment
                })
            else:
                self.stats['unassigned'] += 1
    
    def write_outputs(self, output_prefix):
        """Write output files"""
        # Filtered reads
        output_file = f"{output_prefix}_filtered.fastq"
        with open(output_file, 'w') as f:
            SeqIO.write(self.filtered_reads, f, "fastq")
        logger.info(f"Wrote {len(self.filtered_reads)} filtered reads to {output_file}")
        
        # Reads per reference
        for ref_name, assignments in self.assignments.items():
            safe_name = ref_name.replace(':', '_').replace('/', '_')
            ref_file = f"{output_prefix}_{safe_name}_reads.fastq"
            records = [a['record'] for a in assignments]
            with open(ref_file, 'w') as f:
                SeqIO.write(records, f, "fastq")
            logger.info(f"Wrote {len(records)} reads assigned to {ref_name}")
            
            # Write insert sequences if applicable
            inserts = [a['alignment'].get('insert_sequence') 
                      for a in assignments 
                      if a['alignment'].get('insert_sequence')]
            
            if inserts:
                insert_file = f"{output_prefix}_{safe_name}_inserts.fasta"
                with open(insert_file, 'w') as f:
                    for i, insert in enumerate(inserts):
                        f.write(f">{assignments[i]['read_id']}\n{insert}\n")
                logger.info(f"  Extracted {len(inserts)} insert sequences")
        
        # Summary statistics
        self.write_summary(output_prefix)
    
    def write_summary(self, output_prefix):
        """Write summary statistics"""
        summary_file = f"{output_prefix}_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("="*60 + "\n")
            f.write("LIBRARY ALIGNMENT SUMMARY\n")
            f.write("="*60 + "\n\n")
            
            f.write("Overall Statistics:\n")
            f.write(f"  Total reads: {self.stats['total']}\n")
            f.write(f"  Passed filters: {self.stats['passed_filter']} "
                   f"({100*self.stats['passed_filter']/max(1,self.stats['total']):.1f}%)\n")
            f.write(f"  Assigned reads: {self.stats['assigned']} "
                   f"({100*self.stats['assigned']/max(1,self.stats['passed_filter']):.1f}%)\n")
            f.write(f"  Unassigned: {self.stats['unassigned']}\n")
            f.write(f"  Ambiguous: {self.stats.get('ambiguous', 0)}\n")
            
            f.write("\nAssignments per Reference:\n")
            for ref_name in self.references.keys():
                count = self.stats.get(f'assigned_to_{ref_name}', 0)
                spans_n = self.stats.get(f'spans_n_in_{ref_name}', 0)
                pct = 100 * count / max(1, self.stats['assigned'])
                
                f.write(f"\n  {ref_name}:\n")
                f.write(f"    Total: {count} ({pct:.1f}%)\n")
                
                if self.n_regions.get(ref_name):
                    f.write(f"    Spans N region: {spans_n} "
                           f"({100*spans_n/max(1,count):.1f}% of assigned)\n")
                    
                    # Report N region details
                    for n_start, n_end in self.n_regions[ref_name]:
                        f.write(f"    N region: {n_start}-{n_end} ({n_end-n_start} bp)\n")
        
        logger.info(f"Summary written to {summary_file}")
        
        # Print to console
        with open(summary_file, 'r') as f:
            print(f.read())


def main():
    parser = argparse.ArgumentParser(
        description="Align nanopore reads to plasmid library with N-aware scoring",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("input_fastq", help="Input FASTQ file")
    parser.add_argument("-r", "--reference", required=True,
                       help="Reference FASTA with library variants (N = variable region)")
    parser.add_argument("-o", "--output-prefix", default="library_alignment",
                       help="Output file prefix")
    parser.add_argument("-l", "--min-length", type=int, default=0,
                       help="Minimum read length")
    parser.add_argument("-q", "--min-quality", type=float, default=0,
                       help="Minimum mean quality score")
    parser.add_argument("-f", "--flank-size", type=int, default=50,
                       help="Size of flanking region to use for mapping")
    parser.add_argument("--min-flank-identity", type=float, default=0.9,
                       help="Minimum identity in flanking regions")
    parser.add_argument("-c", "--circular", action="store_true", default=True,
                       help="Treat plasmids as circular")
    
    args = parser.parse_args()
    
    if not HAS_MINIMAP2:
        logger.error("This tool requires minimap2. Install with: pip install mappy")
        sys.exit(1)
    
    # Process
    aligner = LibraryAligner(
        reference_fasta=args.reference,
        min_length=args.min_length,
        min_quality=args.min_quality,
        flank_size=args.flank_size,
        min_flank_identity=args.min_flank_identity,
        circular=args.circular
    )
    
    aligner.process_reads(args.input_fastq)
    aligner.write_outputs(args.output_prefix)
    
    logger.info("Complete!")


if __name__ == "__main__":
    main()
