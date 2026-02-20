# insert_extract.py

A command-line tool for Oxford Nanopore Technology (ONT) long-read sequencing data. It aligns reads to a single circular or linear plasmid reference, identifies contiguous stretches of ambiguous N bases in the reference, extracts the read subsequence that spans each N region (the "insert"), and matches those inserts against user-provided lookup files. Results are reported per read in a TSV file and summarised in a multi-page PDF report.

---

## Contents

- [Background](#background)
- [Requirements](#requirements)
- [Installation](#installation)
- [How It Works](#how-it-works)
- [Usage](#usage)
- [Arguments](#arguments)
- [Lookup File Formats](#lookup-file-formats)
- [Worked Example: MPRAv3 Library](#worked-example-mprav3-library)
- [Output Files](#output-files)
- [PDF Report Pages](#pdf-report-pages)
- [Notes and Caveats](#notes-and-caveats)

---

## Background

Massively-parallel reporter assays (MPRAs) and similar synthetic-biology constructs frequently place variable elements — regulatory sequences, barcodes, guide RNAs — inside a defined plasmid backbone. The variable positions are represented as runs of `N` bases in the reference sequence. This tool answers the question: *which specific sequence was cloned into each N region in each individual molecule?*

The key design decisions:

- A single reference FASTA with one or more N regions is used; every N stretch is detected automatically.
- Each N region is assigned a separate lookup file by the user.
- Inserts are extracted using minimap2 CS-string walking, correctly accounting for read orientation, boundary insertions/deletions, and the circular topology of plasmid backbones.
- Short inserts (< 50 bp) are matched by Levenshtein neighbourhood lookup rather than k-mer seeding, tolerating sequencing errors in barcodes.

---

## Requirements

| Package | Purpose |
|---|---|
| Python 3.8+ | |
| [mappy](https://github.com/lh3/minimap2/tree/master/python) | minimap2 Python bindings for read-to-reference and insert-to-FASTA alignment |
| numpy | numerical operations |
| pandas | tabular output |
| matplotlib | PDF report |

---

## Installation

```bash
pip install mappy numpy pandas matplotlib
```

Or with conda:

```bash
conda install -c bioconda mappy
conda install numpy pandas matplotlib
```

---

## How It Works

### 1. Reference scanning

On startup, the tool scans the reference FASTA for contiguous runs of `N` bases and reports each region with its 0-based coordinates and length. These are numbered 1, 2, 3, … in order of position. Each N region must be assigned a lookup file via `--n-region`.

### 2. Read alignment

Every read passing quality and length filters is aligned to the reference using minimap2 (`map-ont` preset). When `--circular` is enabled, the reference is doubled internally so that reads spanning the origin of a circular plasmid produce a single contiguous alignment rather than two split alignments. Only the highest-scoring alignment (by matching bases, `mlen`) is used.

### 3. Insert extraction

The CS string of the best alignment is walked token by token to locate the precise query (read) coordinates that correspond to the N region on the reference. Because the N bases in the reference do not constrain the aligner, the extracted insert can be longer or shorter than the N region itself — for example, a 200 bp N region may harbour a 220 bp ORF insert, with the extra 20 bp appearing as insertion tokens at the N-region boundaries in the CS string.

For minus-strand alignments, the CS string walks the reverse complement of the query. The tool converts CS-walk coordinates back to original read coordinates before slicing, so both strands are handled correctly.

### 4. Insert matching

Two strategies are used depending on the lookup file type and insert length.

**FASTA lookups (`.fa`, `.fasta`, `.fna`):**
Inserts are aligned to the indexed FASTA with minimap2. For N regions >= 50 bp the `map-ont` preset is used; for N regions < 50 bp the sequences are pre-hashed and Levenshtein matching is used instead (see below), because minimap2 k-mer seeding is unreliable for sequences shorter than its default k-mer size.
A match is called when the alignment covers at least `--min-insert-match-frac` of the insert length (default 0.70).

**TSV lookups (all other extensions):**
The file is loaded into a hash dictionary `{sequence: ID}`. Each insert is queried by Levenshtein neighbourhood lookup: all sequences within `--max-edit-dist` edits of the insert (and its reverse complement) are generated using BFS and checked against the hash in O(1) per variant. This is far more reliable than k-mer-based methods for 20 bp barcodes with sequencing errors.
The default edit distance is automatic: `max(1, round(0.05 × insert_length))`, giving distance 1 for 20 bp barcodes (allowing one substitution, insertion, or deletion). Set `--max-edit-dist 0` to require exact matching.

### 5. Cross-checking

With `--cross-check`, reads that obtained an ID match in two or more N regions are evaluated for agreement: the matched IDs from each region are compared pairwise. Because each N region independently identifies the same library member (ORF + barcode), high agreement is the expected outcome in a well-constructed library.

---

## Usage

```
python3 insert_extract.py \
    --fastq <reads.fastq[.gz]> \
    --reference <plasmid.fasta> \
    --n-region INDEX:FILE [--n-region INDEX:FILE ...] \
    [options]
```

`--n-region` must be specified at least once. `INDEX` is the 1-based position of the N stretch in the reference (first stretch = 1). `FILE` is the lookup file for that stretch.

---

## Arguments

| Argument | Default | Description |
|---|---|---|
| `--fastq` | *(required)* | ONT FASTQ file, plain or gzip-compressed |
| `--reference` | *(required)* | Single-sequence plasmid FASTA containing N-base stretches |
| `--n-region INDEX:FILE` | *(required, repeatable)* | Assign the Nth N-stretch to a lookup file. May be given multiple times for multiple N regions |
| `--tsv-seq-col` | `1` | 1-based column index of the sequence in TSV lookup files |
| `--tsv-id-col` | `2` | 1-based column index of the ID in TSV lookup files |
| `--cross-check` | off | Report pairwise ID agreement across N regions within each read |
| `--circular` | off | Treat reference as circular. Doubles the reference sequence internally so reads spanning the plasmid origin align as a single hit. Recommended for all plasmid data |
| `--min-qual` | `10.0` | Minimum average Phred quality score to accept a read |
| `--min-length` | `2000` | Minimum read length in bp |
| `--min-insert-match-frac` | `0.70` | Minimum fraction of insert bases aligned to call a FASTA lookup match |
| `--max-edit-dist` | `-1` (auto) | Maximum Levenshtein edit distance for short-insert matching (< 50 bp). `-1` = automatic (5% of insert length, minimum 1). `0` = exact match only |
| `--threads` | `4` | Threads passed to minimap2 |
| `--output` | `insert_results.tsv` | Per-read results table |
| `--plots-output` | `insert_report.pdf` | Multi-page PDF report |
| `--no-plots` | off | Skip PDF generation |

---

## Lookup File Formats

### FASTA lookup

Used for longer inserts (>= 50 bp) such as ORFs, enhancers, or sgRNA sequences. Standard FASTA format; extensions `.fa`, `.fasta`, or `.fna` trigger this mode. The sequence header (first word after `>`) is used as the reported ID.

```
>chr1:20788501-20788720
ATGCGATCGATCGATCG...
>random_2291
GCTAGCTAGCTAGCTAG...
```

Indexing is done once at startup with minimap2 (`map-ont` preset, `best_n=5`).

### TSV lookup

Used for short sequences such as barcodes or UMIs. Any file extension other than `.fa`/`.fasta`/`.fna` is treated as TSV. The default layout is:

```
ACGTACGTACGTACGTACGT    chr10:125878091-125878310
TTGCAATCGGAACCTTGACG    random_15795
...
```

Column indices are configurable with `--tsv-seq-col` and `--tsv-id-col`. The entire file is loaded into a hash dictionary at startup. For very large barcode tables (100 M+ lines, ~8 GB RAM) loading takes 2–3 minutes; plan memory accordingly.

---

## Worked Example: MPRAv3 Library

This example reflects a real MPRA library built on the pGL4.23 backbone. The plasmid contains:

- **N region 1** (positions 1943–2143, 200 bp): ORF insert site. Each molecule carries one of ~150,000 candidate regulatory sequences (~220 bp each, slightly longer than the N region).
- **N region 2** (positions 3048–3068, 20 bp): Barcode site. Each molecule carries a 20 bp barcode from a 116 M-entry barcode table.

Both inserts are present on the same molecule, so a read spanning the full plasmid should yield a matched ORF and a matched barcode that agree on library identity.

### Reference preparation

The reference FASTA (`pgl4-23mprav3gfp.fasta`) contains the backbone with Ns at the two variable positions:

```
>pGL4.23:MPRAv3:gfp
ATGCG...NNNNN(200)...GCTA...NNNNN(20)...TGCA
```

### Run command

```bash
python3 insert_extract.py \
    --fastq 2462529_OL57GFP-Final-Scaleup.fastq \
    --reference pgl4-23mprav3gfp.fasta \
    --n-region 1:OL57_reference.fasta \
    --n-region 2:OL57.merged.match.enh.mapped.barcode.ct.parsed \
    --circular \
    --cross-check \
    --min-qual 10 \
    --min-length 2000 \
    --threads 4 \
    --output insert_results.tsv \
    --plots-output insert_report.pdf
```

**Lookup file notes:**

- `OL57_reference.fasta` — FASTA with 149,670 ORF sequences (220 bp each). Matched by minimap2 with `map-ont` preset; a match requires >= 70% of insert bases aligned.
- `OL57.merged.match.enh.mapped.barcode.ct.parsed` — two-column TSV with 116 M barcodes (column 1: 20 bp sequence, column 2: genomic ID). Matched by Levenshtein neighbourhood lookup with automatic edit distance of 1 (5% of 20 bp).

**Why `--circular`:** The barcode site (N region 2) is only 607 bp from the end of the linearised reference. Without `--circular`, reads whose alignment spans the plasmid origin are split into two separate hits and only the longer one is used, causing the barcode region to be missed for approximately 40% of reads.

### Expected output summary

```
Total reads:               1,462
Passed filters:              601
Mapped to reference:         475

N-region 1 (pos 1,943–2,143, 200 bp Ns):
  Inserts extracted:  396
  Avg insert length:  209.4 bp
  Matched to lookup:  382    (96.5% of extracted)
  Unique IDs matched: 334

N-region 2 (pos 3,048–3,068, 20 bp Ns):
  Inserts extracted:  430
  Avg insert length:  19.9 bp
  Matched to lookup:  333    (77.4% of extracted)
  Unique IDs matched: 302

Cross-check (N1 vs N2):
  Reads with both matches: 272
  Same ID:                 265 / 272  (97.4%)
```

---

## Output Files

### `insert_results.tsv`

One row per mapped read. Columns:

| Column | Description |
|---|---|
| `read_id` | Read name from FASTQ header |
| `read_length` | Read length in bp |
| `avg_qual` | Average Phred quality score |
| `ref_start` | Alignment start on reference (0-based) |
| `ref_end` | Alignment end on reference |
| `strand` | Alignment strand (`+` or `-`) |
| `match_frac` | Fraction of read bases matching the reference |
| `n{i}_insert_seq` | Extracted insert sequence for N region i |
| `n{i}_insert_len` | Insert length in bp (0 if not extracted) |
| `n{i}_match_id` | ID of the best lookup match (empty if none) |
| `n{i}_match_frac` | Match fraction for the insert-to-lookup alignment |
| `n{i}_match_strand` | Strand of the insert match |

Columns are repeated for each N region that was assigned a lookup file (n1, n2, …).

---

## PDF Report Pages

### Page 1 — Processing Funnel

A stacked bar chart showing read counts at each stage of the pipeline: QC passed, mapped to reference, then per N region: insert extracted and insert matched. Each bar is split into a coloured (passing) and grey (dropout) segment against the QC-passed baseline. Counts and percentages are annotated on each bar.

### Page 2 — ID Matches and Agreement

Scoped to reads with at least one ID match. For each N region, two bars are shown:

- **Has match** — number of reads with an ID match for that region, as a fraction of all any-match reads.
- **Agrees with N{j}** — number of reads where both this region and region j were matched and returned the same ID, as a fraction of this region's own match count (e.g. `265/333 = 79.6%`).

If more than two N regions are assigned, one agreement bar is drawn per pairwise comparison.

### Page 3 — Alignment Score Distributions

One histogram per N region showing the distribution of insert-to-lookup match fractions across all extracted inserts, including those that did not match (fraction = 0). A red dashed vertical line marks the match fraction cutoff (`--min-insert-match-frac`). The title states the total number matched, the number with no alignment at all, and the total extracted.

### Page 4 — Insert Length Distributions

One histogram per N region. Matched and unmatched inserts are overlaid in different colours. A red dashed line shows the N-region size in the reference; a blue dashed line shows the mean insert length. This page is useful for diagnosing size bias or systematic extraction errors.

### Page 5 — UpSet Plot

Shows all unique combinations of four boolean categories across the 475 mapped reads: N1 insert extracted, N1 ID matched, N2 insert extracted, N2 ID matched. Bar height = read count for that combination; filled dots in the matrix below indicate which categories are active. Total category sizes are annotated on the right.

---

## Notes and Caveats

**Circular plasmids.** Always use `--circular` when sequencing plasmids. Reads that span the cloning origin will otherwise have their alignment split, and N regions near the end of the linearised sequence will be missed.

**Minus-strand reads.** The CS string produced by minimap2 for reverse-strand alignments describes the alignment of the reverse complement of the query against the forward reference. The tool accounts for this when converting CS-walk coordinates back to original read positions, ensuring both strands are extracted correctly.

**Insert length vs N-region length.** The N region in the reference is a placeholder of fixed length, but the actual insert in the molecule may differ. The extraction algorithm walks the CS string to find the precise boundaries of the insert, including sequence that minimap2 reports as insertions at the N-region flanks. The insert length distribution (page 4) will show how much variation there is.

**TSV memory usage.** Loading a 116 M-entry barcode table requires approximately 8 GB of RAM. A second reverse-complement dictionary is not built; instead, the reverse complement of each query insert is computed on the fly (the number of queries is small). Loading takes 2–3 minutes on a standard desktop.

**Levenshtein matching for barcodes.** At the default automatic edit distance of 1 for 20 bp barcodes, approximately 164 forward variants and 164 reverse-complement variants are generated per query and looked up in the hash table. This is fast (< 1 ms per read) and raises the match rate substantially compared to exact matching. Increasing `--max-edit-dist` to 2 generates roughly 10,000–20,000 variants per query; this is still fast per read but may increase the rate of false-positive matches if barcodes are closely related.

**FASTA lookup for short inserts.** If a FASTA file is assigned to an N region shorter than 50 bp, minimap2 is not used; sequences are loaded into a hash dictionary and Levenshtein matching is applied instead. This is because the default minimap2 k-mer size (21) exceeds the insert length, making k-mer seeding impossible.

**Multiple N regions.** The script supports any number of N regions in the reference, each independently assigned a lookup file. N regions without an assigned lookup are silently skipped.
