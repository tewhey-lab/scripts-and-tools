# FASTQ Sequence Matcher and Insert Checker

Tool for searching and extracting specific sequences from FASTQ files with support for fuzzy matching, flanking sequence extraction, and matching to reference files.

## Key Features

- **Sequence Matching**: Direct sequence search or flanking sequence extraction with fuzzy matching
- **Pre-filtering**: Optional read length and quality score filtering to improve processing speed
- **Multiple Algorithms**: Hamming, Levenshtein, homology-based matching, or ultra-fast minimap2 (mappy)
- **Comprehensive Reports**: Matrix-format statistics, detailed per-read results, and multiple visualizations
- **Advanced Visualizations** (300 DPI, publication-ready):
  - Enhanced UpSet plots with dual-axis information (counts + percentages)
  - Extraction/matching status by region (stacked barplots)
  - Pairwise match concordance analysis
  - Professional styling: grid-free, clean labels, optimized spacing, no overlapping text
- **Performance Optimized**: Hash table lookups, early termination, and optional mappy integration

## Recent Enhancements

**Version Updates:**
- ✨ **Dual output formats**:
  - Human-readable combination statistics (`_combination_stats.txt`)
  - Binary matrix format for computational analysis (`_combination_matrix.txt`)
- 📊 **Enhanced UpSet plots**: Dual-axis labels showing both counts and percentages
- 🎨 **Professional styling**:
  - Grid-free plots with clean labels
  - Optimized title padding prevents text overlap
  - Proper spacing between all plot elements
- 📈 **New visualizations**:
  - Extraction status plots (per-region success rates)
  - Match concordance plots (pairwise barcode analysis)
- 🔍 **Pre-filtering**: Optional read length and quality filtering
- 🚀 **Performance**: Improved data representation and plotting efficiency

## Installation

### Prerequisites

```bash
#PIP INSTALL
# Required packages
pip install biopython pandas python-Levenshtein matplotlib numpy
# Recommended for enhanced UpSet plots with dual-axis labels (falls back to bar chart if not installed)
pip install upsetplot

#CONDA INSTALL
conda install -c conda-forge biopython pandas python-levenshtein matplotlib numpy
# Recommended: for enhanced UpSet plots
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

#### Pre-filtering Options
- `--min-length N`: Minimum read length in bases. Reads shorter than this will be skipped
- `--max-length N`: Maximum read length in bases. Reads longer than this will be skipped
- `--min-quality Q`: Minimum average quality score (Phred scale). Reads below this will be skipped
- `--min-quality-percentile P`: Minimum quality percentile (0-100). E.g., 90 means 90% of bases must be >= min-quality. Requires `--min-quality` to be set

#### Per-Sequence Options (replace # with 1-10)
- `--seq# "SEQUENCE"`: Sequence to search (required)
- `--name# NAME`: Custom name for sequence # (default: seq#). Names are sanitized for use in filenames
- `--seq#-dist N`: Maximum Levenshtein distance for matching (default: 2)
- `--seq#-insert MIN:MAX`: Insert size range for flanking sequences (e.g., "10:30")
- `--ref-file# FILE`: Reference file for matching (FASTA or TSV format)
- `--match-method# {hamming,levenshtein,homology,mappy}`: Matching algorithm for insert matching (default: homology, mappy=ultra-fast)
- `--match-dist# N`: Match threshold for insert matching (default: 5 for distance metrics, 0.8 for homology/mappy)
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
python insert_check.py reads.fastq \
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

### Example 4: Pre-filtering for Read Length and Quality
```bash
# Filter reads: length 50-500bp, average quality Q30+
python insert_check.py reads.fastq \
    --seq1 "GTCGACGAACCTCTAGA-AGATCGGAAGAGCGT" \
    --seq1-dist 4 \
    --seq1-insert 18:22 \
    --ref-file1 barcodes.txt \
    --min-length 50 \
    --max-length 500 \
    --min-quality 30
```

### Example 5: Strict Quality Filtering with Percentile
```bash
# Require 95% of bases to have Q30 or higher
python insert_check.py reads.fastq \
    --seq1 "ATCGATCG-GCTAGCTA" \
    --seq1-insert 18:22 \
    --min-quality 30 \
    --min-quality-percentile 95
```

## Output Files

The tool generates comprehensive reports and visualizations to analyze your sequencing data:

**Text Reports:**
1. Combination statistics (human-readable format)
2. Combination matrix (binary 1/0 format for statistical analysis)
3. Per-read detailed results

**Visualizations:**
4. UpSet plot or combination bar chart (sequence match combinations)
5. Extraction status plot (success rates per region)
6. Match concordance plot (pairwise barcode concordance)

All plots are publication-ready at 300 DPI with professional styling.

---

### 1. `{prefix}_combination_stats.txt`
Summary statistics showing:
- All possible sequence match combinations and their frequencies
- Total read counts and proportions
- Pre-filtering statistics (if filters were applied)
- Individual sequence match rates

Example:
```
Combination	Count	Proportion
AmpR+GFP+Ori+bc+oligo	1015	0.4299
AmpR+GFP+Ori+oligo	287	0.1216
GFP+Ori+bc	198	0.0839
none	89	0.0377

# Summary
Total reads in file: 12000
Filtered reads: 2000 (16.67%)

# Filter reasons
length_45<50: 1200 (60.00% of filtered)
avg_quality_25.3<30: 800 (40.00% of filtered)

Reads analyzed (post-filter): 10000
Reads with at least one match: 9755 (97.55% of analyzed)

# Individual sequence matches
AmpR: 8863 (88.63% of analyzed)
GFP: 7932 (79.32% of analyzed)
Ori: 7850 (78.50% of analyzed)
bc: 5200 (52.00% of analyzed)
oligo: 4100 (41.00% of analyzed)
```

**Use cases:**
- Quick overview of most common combinations
- Summary statistics for reporting
- Human-readable format

### 2. `{prefix}_combination_matrix.txt` **(NEW)**
Binary matrix format with 1/0 columns for each sequence - optimized for computational analysis.

Example:
```
AmpR	GFP	Ori	bc	oligo	Count	Proportion
1	1	1	1	1	1015	0.4299
1	1	1	0	1	287	0.1216
0	1	1	1	0	198	0.0839
0	0	0	0	0	89	0.0377
```

**Features:**
- Each sequence gets its own binary column (1 = present, 0 = absent)
- Easy to import into R, Python pandas, or Excel
- Perfect for statistical analysis and modeling
- Rows sorted by count (descending)
- Includes count and proportion for each combination

**Use cases:**
- Import into statistical software (R, Python, MATLAB)
- Machine learning feature matrices
- Correlation analysis between sequences
- Custom data analysis pipelines

### 3. `{prefix}_combined_results.txt`
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

### 4. `{prefix}_upset_plot.png` or `{prefix}_combination_plot.png`
**UpSet Plot** (primary visualization):
- Interactive-style visualization showing set intersections and sizes
- **Enhanced styling**: Clean design without grid lines, professional appearance
- **Dual-axis information**: Shows both counts AND percentages
  - Intersection bars: Labeled with clean percentages (e.g., "5.2%")
  - Set size bars: Dual-line labels showing count and percentage (e.g., "1234\n12.3")
- **Optimized layout**: Larger 14x10 figure size for better readability
- **Smart sorting**: Sorted by cardinality for easy interpretation
- Includes "No_Match" category for reads with no sequence matches

**Bar Chart** (fallback if upsetplot not installed):
- Shows percentage of reads for each combination
- Clean percentage labels without parentheses
- Grid-free professional styling
- Sorted by frequency (descending)

### 5. `{prefix}_extraction_status.png` **(NEW)**
Stacked barplot showing extraction and matching status for each flanking sequence region.

**Categories:**
- **No Sequence** (red): Flanking sequences not found, extraction failed
- **No Match** (orange): Sequence extracted but not matched to any reference
- **Matched** (green): Sequence extracted and successfully matched to reference

**Features:**
- Shows percentage breakdown for each region
- Includes total read counts above each bar (format: `n=1234`)
- Percentage labels displayed on visible segments (>5%)
- Professional spacing: optimized title padding prevents label overlap
- Only generated when flanking sequences are present

**Use cases:**
- Quality control: Identify regions with low extraction rates
- Optimization: Adjust flanking sequence parameters or error thresholds
- Troubleshooting: Diagnose why certain regions aren't matching

### 6. `{prefix}_match_concordance.png` **(NEW)**
Pairwise comparison plot showing Match ID concordance between different regions.

**Categories:**
- **Same Match ID** (green): Both regions matched to the same reference ID
- **Different Match IDs** (red): Both regions matched but to different references (discordant)
- **Only One Matched** (orange): One region matched, the other didn't
- **Neither Matched** (gray): Neither region matched to any reference

**Features:**
- All pairwise combinations of regions with references
- Percentage breakdown with count labels (format: `n=5678`)
- Identifies concordant vs. discordant matches
- Professional spacing: optimized title padding prevents label overlap
- Multi-line x-axis labels for clarity (format: `region1\nvs\nregion2`)
- Only generated when ≥2 flanking sequences with references exist

**Use cases:**
- Validate dual-barcoding schemes
- Detect chimeric reads (different barcodes in same read)
- Quality control for paired amplicon designs
- Identify cross-contamination or mis-priming events

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
- Use pre-filtering (`--min-length`, `--max-length`, `--min-quality`) to reduce processing time by eliminating low-quality reads early
- Pre-filtering is applied before sequence matching, improving overall throughput

### Troubleshooting

1. **Memory Issues**: Process files in chunks or on a high-memory system
2. **Slow Processing**:
   - Reduce maximum Levenshtein distance
   - Use more specific sequences
   - Consider using Hamming instead of homology for simple matches
   - Apply pre-filtering to eliminate poor quality reads early
3. **No Matches Found**:
   - Check sequence orientation (try reverse complement)
   - Increase error tolerance
   - Verify sequence format and case
   - Check if pre-filtering is too stringent (review filter statistics in output)
4. **Too Many Reads Filtered**:
   - Review filter reasons in the `_combination_stats.txt` file
   - Adjust `--min-quality`, `--min-length`, or `--max-length` thresholds
   - Consider using `--min-quality-percentile` for more flexible quality filtering
5. **Quality Filtering Not Working**:
   - Ensure your FASTQ file contains quality scores (some FASTA files may be mislabeled)
   - Check that quality scores are in standard Phred format

## License

GNU GENERAL PUBLIC LICENSE v3

## License

Code was partially built using Claude Opus 4. Use with caution.

## Contact

info@tewheylab.org
