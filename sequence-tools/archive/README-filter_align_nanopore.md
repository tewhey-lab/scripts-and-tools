# Nanopore Library Alignment Tool (`filter_align_nanopore.py`)

A specialized tool for aligning Nanopore reads to plasmid libraries where variable inserts are represented by `N` regions in the reference sequence.

## Why this tool?
Standard aligners penalize mismatches. If you align a read with a specific insert (e.g., `ATGC...`) to a reference with `NNNN...`, standard aligners will treat every position as a mismatch or try to introduce gaps, resulting in poor alignment scores or failed alignments.

This tool solves this by:
1.  **Backbone Alignment**: It temporarily removes the `N` regions to create a "backbone" reference.
2.  **Insert Detection**: It aligns reads to this backbone using `minimap2` and detects the insert as a large insertion in the read at the specific junction where the `N` region was.
3.  **Insert Extraction**: It extracts the sequence corresponding to the `N` region for downstream analysis.

## Dependencies

- python >= 3.6
- `mappy` (minimap2 python binding)
- `biopython`
- `numpy`
- `tqdm`

Install with:
```bash
pip install mappy biopython numpy tqdm
```

## Usage

```bash
python3 filter_align_nanopore.py input.fastq -r reference.fasta [options]
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `input_fastq` | Path to input FASTQ file (Nanopore reads) | Required |
| `-r`, `--reference` | Path to reference FASTA. Variable regions must be marked with `N`s. | Required |
| `-o`, `--output-prefix` | Prefix for output files | `library_alignment` |
| `-l`, `--min-length` | Minimum read length to process | 0 |
| `-q`, `--min-quality` | Minimum mean quality score (Phred) | 0 |
| `-f`, `--flank-size` | Size of flanking region context (for logging) | 50 |
| `-c`, `--circular` | Treat plasmids as circular | True |

## Input Format

**Reference FASTA**:
Must contain `N` characters where the variable insert is expected.
```fasta
>plasmid_A
ATGCATGCATGC...NNNNNNNNNNNNNNNNNNNN...ATGCATGC
```

## Outputs

1.  `{prefix}_filtered.fastq`: Reads that passed length/quality filters.
2.  `{prefix}_{ref_name}_reads.fastq`: Reads assigned to each reference.
3.  `{prefix}_{ref_name}_inserts.fasta`: Extracted insert sequences for each reference.
4.  `{prefix}_summary.txt`: Detailed statistics on assignment rates and N-region spanning.

## Algorithm Details

The tool uses `mappy` (minimap2) with the `map-ont` preset.
1.  **Preprocessing**: For each reference, `N` regions are identified and removed to create a "backbone" sequence.
2.  **Alignment**: Reads are mapped to these backbones.
3.  **Scoring**: 
    - Reads are assigned to the backbone with the highest alignment score.
    - If a read has an insertion at the coordinate corresponding to the `N` region, it is marked as "spanning the N region".
    - The alignment score is boosted by the length of this insert to ensure reads with inserts are not penalized compared to the backbone.
