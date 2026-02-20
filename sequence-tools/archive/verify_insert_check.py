import unittest
import tempfile
import os
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

# Import the module to test
import insert_check

# Disable logging for tests
logging.getLogger('insert_check').setLevel(logging.CRITICAL)

class TestInsertCheckRefactor(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.fastq_file = os.path.join(self.test_dir, "test.fastq")
        self.ref_file = os.path.join(self.test_dir, "ref.fasta")
        
        # Create synthetic data
        # Read 1: Perfect match to Ref1 (flanking)
        # Read 2: No match
        # Read 3: Match to Ref2 (direct)
        
        self.left_flank = "ATCG"
        self.right_flank = "GCTA"
        self.insert1 = "AAAA" # Ref1
        self.insert2 = "TTTT" # Ref2
        
        reads = [
            SeqRecord(Seq(f"NNN{self.left_flank}{self.insert1}{self.right_flank}NNN"), id="read1", description=""),
            SeqRecord(Seq("NNNNNNNNNNNNNNNNNNNN"), id="read2", description=""),
            SeqRecord(Seq(f"NNN{self.left_flank}{self.insert2}{self.right_flank}NNN"), id="read3", description="")
        ]
        
        # Add dummy quality scores (required for FASTQ)
        for read in reads:
            read.letter_annotations["phred_quality"] = [30] * len(read)
            
        SeqIO.write(reads, self.fastq_file, "fastq")
        
        refs = [
            SeqRecord(Seq(self.insert1), id="Ref1", description=""),
            SeqRecord(Seq(self.insert2), id="Ref2", description="")
        ]
        SeqIO.write(refs, self.ref_file, "fasta")
        
        self.sequence_specs = {
            'seq1': {
                'is_flanking': True,
                'left_seq': self.left_flank,
                'right_seq': self.right_flank,
                'min_insert': 2,
                'max_insert': 10,
                'max_distance': 1,
                'match_method': 'hamming',
                'match_dist': 0,
                'reverse_complement': False,
                'reference': {self.insert1: "Ref1", self.insert2: "Ref2"},
                'reference_file': self.ref_file
            }
        }

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_normal_mode(self):
        """Test normal in-memory mode"""
        results = insert_check.process_fastq_file(
            self.fastq_file, 
            self.sequence_specs, 
            os.path.join(self.test_dir, "output_normal"),
            stream_references=False
        )
        
        self.assertEqual(results['total_reads'], 3)
        self.assertEqual(results['processed_reads'], 3)
        
        # Check read 1
        r1 = results['all_read_data']['read1']
        self.assertTrue(r1['matches']['seq1'])
        self.assertEqual(r1['extracted']['seq1'], self.insert1)
        self.assertEqual(r1['match_ids']['seq1'], "Ref1")
        
        # Check read 2
        r2 = results['all_read_data']['read2']
        self.assertFalse(r2['matches']['seq1'])
        
        # Check read 3
        r3 = results['all_read_data']['read3']
        self.assertTrue(r3['matches']['seq1'])
        self.assertEqual(r3['extracted']['seq1'], self.insert2)
        self.assertEqual(r3['match_ids']['seq1'], "Ref2")

    def test_streaming_mode(self):
        """Test streaming reference mode"""
        # In streaming mode, 'reference' should be None for streaming seqs
        specs = self.sequence_specs.copy()
        specs['seq1']['reference'] = None # Simulate what main() does
        
        results = insert_check.process_fastq_file(
            self.fastq_file, 
            specs, 
            os.path.join(self.test_dir, "output_stream"),
            stream_references=True
        )
        
        self.assertEqual(results['total_reads'], 3)
        self.assertEqual(results['processed_reads'], 3)
        
        # Check read 1
        r1 = results['all_read_data']['read1']
        self.assertTrue(r1['matches']['seq1'])
        self.assertEqual(r1['extracted']['seq1'], self.insert1)
        self.assertEqual(r1['match_ids']['seq1'], "Ref1")
        
        # Check read 2
        r2 = results['all_read_data']['read2']
        self.assertFalse(r2['matches']['seq1'])
        
        # Check read 3
        r3 = results['all_read_data']['read3']
        self.assertTrue(r3['matches']['seq1'])
        self.assertEqual(r3['extracted']['seq1'], self.insert2)
        self.assertEqual(r3['match_ids']['seq1'], "Ref2")

if __name__ == '__main__':
    unittest.main()
