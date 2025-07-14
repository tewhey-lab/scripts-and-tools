# FASTQ Sequence Matcher and Insert Checker

Tool for searching and extracting specific sequences from FASTQ files with support for fuzzy matching, flanking sequence extraction, and matching to reference files.

## Installation

### Prerequisites

```bash
#PIP INSTALL
# Required packages
pip install biopython pandas python-Levenshtein matplotlib numpy 
# Optional for UpSet plots (falls back to bar chart if not installed)
pip install upsetplot

#CONDA INSTALL
conda install -c conda-forge biopython pandas python-levenshtein matplotlib numpy
# Optional: for UpSet plots
conda install -c conda-forge upsetplot
```

## Usage

### Basic Command Structure

```bash
python insert_check.py input.fastq [options]
```

### Required Parameters

- `input.fastq`: Input FASTQ file to process
- At least one sequence for matching (`--seq1`, `--seq2`, etc.)

### Optional Parameters

#### General Options
- `--output-prefix PREFIX`: Prefix for output files (default: "output")
- `--log-level {DEBUG,INFO,WARNING,ERROR}`: Logging verbosity (default: INFO)
- `--no-progress`: Disable progress reporting
- `--skip-count`: Skip initial read counting for faster startup

#### Per-Sequence Options (replace # with 1-10)
- `--seq# "SEQUENCE"`: Sequence to search (required)
- `--seq#-dist N`: Maximum Levenshtein distance for matching (default: 2)
- `--seq#-insert MIN:MAX`: Insert size range for flanking sequences (e.g., "10:30")
- `--ref-file# FILE`: Reference file for matching (FASTA or TSV format)
- `--match-method# {hamming,levenshtein,homology}`: Matching algorithm for insert matching (default: homology)
- `--match-dist# N`: Match threshold for insert matching (default: 5 for distance metrics, 0.8 for homology)
- `--seq#-rc`: Reverse complement extracted sequences before matching

### Sequence Specification Formats

1. **Direct Sequence Search**: `"ATCGATCG"`
   - Searches for exact or fuzzy matches of the sequence
   - Counts occurrences without extraction

2. **Flanking Sequence Search**: `"ATCG-TGCA"`
   - Searches for the flanking sequences (ATCG and TGCA)
   - Extracts the variable region between them
   - Useful for barcode extraction or amplicon analysis
     
## Examples

### Example 1: Simple Direct Sequence Search
```bash
python fastq_sequence_matcher.py reads.fastq \
    --seq1 "TCCTCAGTCGCGATCGAACA" \
    --seq1-dist 3
```

### Example 2: Flanking Sequence with Barcode Extraction
```bash
python insert_check.py reads.fastq \
    --seq1 "GTCGACGAACCTCTAGA-AGATCGGAAGAGCGT" \
    --seq1-dist 4 \
    --seq1-insert 18:22 \
    --ref-file1 barcodes.ct.parsed \
    --match-method1 hamming \
    --match-dist1 0
```

### Example 3: Multiple Sequences with Different Parameters
```bash
python insert_check.py reads.fastq \
    --seq1 "GTCGACGAACCTCTAGA-AGATCGGAAGAGCGT" \
    --seq1-dist 4 \
    --seq1-insert 18:22 \
    --ref-file1 barcodes.ct.parsed \
    --match-method1 hamming \
    --match-dist1 0 \
    --seq2 "TCCTCAGTCGCGATCGAACA" \
    --seq2-dist 3 \
    --seq3 "GCAGGACTGGCCGCTTGACG-CACTGCGGCTCCTGCGATTG" \
    --seq3-dist 5 \
    --seq3-insert 100:300 \
    --ref-file3 OL54_reference.fasta \
    --match-method3 homology \
    --match-dist3 0.8 \
    --output-prefix OutputFileName
```

## Output Files

### 1. `{prefix}_combination_stats.txt`
Summary statistics showing:
- All possible sequence match combinations and their frequencies
- Total read counts and proportions
- Individual sequence match rates

Example:
```
Combination     Count   Proportion
none           245     0.0245
seq1           1823    0.1823
seq2           892     0.0892
seq1+seq2      4521    0.4521
seq1+seq2+seq3 2519    0.2519

# Summary
Total reads: 10000
Reads with at least one match: 9755 (97.55%)

# Individual sequence matches
seq1: 8863 (88.63%)
seq2: 7932 (79.32%)
seq3: 2519 (25.19%)
```

### 2. `{prefix}_combined_results.txt`
Tab-delimited file with one row per read containing:
- `Read_ID`: Sequence identifier from FASTQ
- For each sequence:
  - `seq#_match`: Yes/No indicating if sequence was found
  - `seq#_extracted`: Extracted sequence (for flanking) or NA
  - `seq#_match_id`: Best matching reference ID or NO_MATCH/NA
  - `seq#_hamming/levenshtein/homology`: Similarity scores (0-1)

Example:
```
Read_ID    seq1_match  seq1_extracted  seq1_match_id  seq1_hamming  seq1_levenshtein  seq1_homology  seq2_match
M01234:1   Yes         ATCGATCGATCGATCG BC001         1.000         0.950            0.980          No
M01234:2   No          NA              NA             NA            NA               NA             Yes
M01234:3   Yes         GCTAGCTAGCTAGCTA NO_MATCH      0.000         0.000            0.000          Yes
```

### 3. `{prefix}_upset_plot.png` or `{prefix}_combination_plot.png`
- **UpSet Plot**: Interactive-style visualization showing set intersections and sizes
- **Bar Chart** (fallback): Shows percentage of reads for each combination
- Includes "No_Match" category for reads with no sequence matches

## Reference File Formats

### FASTA Format
```
>barcode001
ATCGATCGATCGATCGATCG
>barcode002
GCTAGCTAGCTAGCTAGCTA
```

### TSV Format
Two columns: sequence and identifier (no header)
```
ATCGATCGATCGATCGATCG    barcode001
GCTAGCTAGCTAGCTAGCTA    barcode002
```

## Advanced Features

### Matching Methods Explained

1. **Hamming Distance** (`hamming`)
   - Counts position-wise mismatches
   - Requires sequences of equal length
   - Good for barcodes with known fixed length

2. **Levenshtein Distance** (`levenshtein`)
   - Allows insertions, deletions, and substitutions
   - Works with sequences of different lengths
   - More flexible but computationally intensive

3. **Homology** (`homology`)
   - Uses local sequence alignment
   - Best for biological sequence similarity
   - Returns similarity score (0-1)

### Performance Tips

- For large files (>1GB), use `--skip-count` to start processing immediately
- Use `--no-progress` when running in batch scripts
- Hamming distance is fastest for fixed-length sequences
- Set appropriate error thresholds to balance sensitivity and specificity

### Troubleshooting

1. **Memory Issues**: Process files in chunks or on a high-memory system
2. **Slow Processing**: 
   - Reduce maximum Levenshtein distance
   - Use more specific sequences
   - Consider using Hamming instead of homology for simple matches
3. **No Matches Found**: 
   - Check sequence orientation (try reverse complement)
   - Increase error tolerance
   - Verify sequence format and case

## License

GNU GENERAL PUBLIC LICENSE v3

## License

Code was partially built using Claude Opus 4. Use with caution.

## Contact

info@tewheylab.org
