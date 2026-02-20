# plasmid_align.py

A command-line tool for assigning Oxford Nanopore Technology (ONT) reads to a library of plasmid
reference sequences. Each read is aligned to every reference using minimap2 (via the `mappy`
Python bindings) and assigned to its best-matching plasmid, or flagged as ambiguous when two
references score similarly.

---

## Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Arguments](#arguments)
- [Scoring](#scoring)
- [N-masked references](#n-masked-references)
- [Circular mode](#circular-mode)
- [Output files](#output-files)
- [Visual report](#visual-report)
- [Design notes and limitations](#design-notes-and-limitations)

---

## Overview

The script streams reads from a FASTQ file, applies quality and length filters, aligns each
passing read against all provided reference FASTA files with the `map-ont` minimap2 preset, and
assigns each read to the plasmid with the highest **match fraction** score.  Reads whose top two
scores fall within a user-defined threshold of each other are called **ambiguous** rather than
assigned.

---

## Requirements

All dependencies are available in a conda environment.  Create one if needed:

```bash
conda create -n sandbox python=3.10 mappy pandas numpy matplotlib
conda activate sandbox
```

Or run directly through an existing environment:

```bash
conda run -n sandbox python3 plasmid_align.py --help
```

| Package | Purpose |
|---|---|
| `mappy` | Python bindings for minimap2; performs alignment |
| `numpy` | Array operations, N-region cumulative sums |
| `pandas` | TSV output construction |
| `matplotlib` | PDF report generation |

---

## Quick start

The example below matches the dataset in this directory: a GFP scaleup ONT run aligned against
three plasmid references, all of which are circular.

```bash
conda run -n sandbox python3 plasmid_align.py \
    --fastq 2462529_OL57GFP-Final-Scaleup.fastq \
    --fasta pgl4-23lucxba.fasta \
            pgl4-23mprav3_deltagfp.fasta \
            pgl4-23mprav3gfp.fasta \
    --circular \
    --output results.tsv \
    --summary summary.tsv \
    --plots-output report.pdf
```

Expected runtime: ~30 seconds for 600 passing reads against three ~3 kbp references on a
4-thread laptop.

---

## Arguments

| Flag | Default | Description |
|---|---|---|
| `--fastq` | *(required)* | ONT FASTQ file.  Plain text or gzip-compressed (`.fastq.gz`). |
| `--fasta` | *(required)* | One or more plasmid FASTA files.  Each file may contain multiple sequences separated by `>` headers.  Pass multiple files as space-separated arguments. |
| `--min-qual` | `10.0` | Minimum average Phred quality score.  Reads below this threshold are excluded before alignment. |
| `--min-length` | `2000` | Minimum read length in bp.  Reads shorter than this are excluded.  For full-plasmid ONT reads a value of 2000–3000 is typical. |
| `--ambiguity-delta` | `0.05` | A read is called ambiguous when the gap between the best and second-best match fraction is less than `delta * best_score`.  Increase to flag more reads as ambiguous; decrease for stricter uniqueness. |
| `--circular` | off | Treat all references as circular.  Each sequence is internally doubled before indexing so that reads spanning the linearisation junction align correctly.  Almost always required for plasmid data. |
| `--threads` | `4` | Thread count passed to minimap2 index building. |
| `--output` | `plasmid_assignments.tsv` | Per-read TSV output file (see [Output files](#output-files)). |
| `--summary` | `plasmid_summary.tsv` | Per-plasmid summary TSV output file. |
| `--plots-output` | `plasmid_alignment_report.pdf` | Path for the multi-page PDF visual report. |
| `--no-plots` | off | Skip PDF generation entirely (faster; useful for batch runs). |
| `--fill-n` | off | Attempt to fill N-masked regions in each reference using homologous sequence from another reference in the library before building aligners.  See [N-masked references](#n-masked-references). |

---

## Scoring

Each read is scored against every reference using a **match fraction**:

```
match_frac = identity * query_coverage
           = (effective_mlen / effective_blen) * (q_span / read_len)
```

where:

- `q_span` — number of query bases in the alignment (`q_en - q_st`)
- `read_len` — total length of the read
- `effective_mlen / effective_blen` — alignment identity after excluding N-position contributions
  (see below)

Using `read_len` as the sole denominator makes the score purely read-centric.  Backbone-only reads
that align with equal quality to multiple references receive the same score for each and are
correctly flagged as ambiguous.  Earlier versions used `max(read_len, orig_ref_len)` as the
denominator, which introduced a systematic bias toward smaller references: a reference 30% shorter
than another would outscore it for any backbone-only read regardless of actual alignment quality.

### Ambiguity threshold

A read is called **ambiguous** when:

```
best_score - second_score  <  ambiguity_delta * best_score
```

The `best_plasmid` field in the output still records which reference had the highest score, but
`is_ambiguous` is `True` and the read is not counted as a unique assignment.

---

## N-masked references

Some references contain runs of `N` bases representing masked or unknown sequence.  minimap2
treats reference `N` positions as wildcards that match any query base, which has two effects:

1. **N positions inflate both `mlen` and `blen`** — because they match any base, they are counted
   as matches in `mlen` and as aligned columns in `blen`, artificially elevating computed identity.
2. **Large N runs create anchor gaps** — a 200 bp N run suppresses minimizer seeding across that
   region, potentially routing circular alignments around the masked region and reducing coverage.

### N-aware identity (always active)

The script pre-computes a cumulative N-count array for every reference.  For each alignment, the
number of N bases within the aligned reference range is subtracted from both `mlen` and `blen`
before computing identity:

```
effective_mlen = max(0, mlen - n_in_hit)
effective_blen = max(1, blen - n_in_hit)
identity       = effective_mlen / effective_blen
```

This prevents N regions from artificially inflating the identity of references that contain them.
The per-read `{plasmid}__n_in_hit` column in the output TSV records how many N bases fell within
each alignment, allowing downstream inspection.

### `--fill-n` (optional, requires verification)

When enabled, the script attempts to replace N runs in each reference with homologous sequence
drawn from the most similar other reference in the library before building aligners.

**Algorithm:**
1. For each reference containing N runs, find the donor reference with the highest total aligned
   bases (mlen) against it.
2. For each N run `[n_st, n_en)`, extract 40 bp flanking sequences from the target reference.
3. Locate each flank in the donor using a vectorised sliding-window mismatch search (tolerates up
   to 3 mismatches; N bases in the flank are treated as wildcards).
4. Extract the donor sequence between the end of the left-flank match and the start of the
   right-flank match.
5. Skip the fill if the extracted donor sequence is itself all N (i.e. the donor also has an N run
   at the corresponding position) or if the coordinates are implausible.
6. Print a verification table of all fills applied.

**Important:** `--fill-n` requires that donor and target references have N runs at *different*
positions and that the flanking sequences are conserved enough to locate.  If all references in
the library share N runs at the same coordinates, no fills will be applied (the script reports
why each fill was skipped).  Always verify fills biologically before using in production.

---

## Circular mode

With `--circular`, each reference is concatenated with itself before indexing
(`doubled_seq = seq + seq`).  Reads that span the linearisation junction produce two partial
alignments in minimap2 output; the script detects these pairs and merges them into a single
junction-spanning result.  The `is_junction_spanning` column in the output records whether the
best hit for that read was a merged pair or a single continuous alignment that crossed the
junction.

Pileup plots in the PDF report display reads in original coordinates `[0, orig_len)`, with the
junction position marked at `x = 0` and junction-spanning reads outlined in purple.

---

## Output files

### Per-read TSV (`--output`)

One row per read that passed quality and length filters (including unmapped reads).

| Column | Description |
|---|---|
| `read_id` | Read name from the FASTQ header |
| `read_length` | Read length in bp |
| `avg_qual` | Average Phred quality score |
| `best_plasmid` | Name of the highest-scoring reference, or `unmapped` |
| `is_ambiguous` | `True` if the top two scores are within `ambiguity_delta` of each other |
| `is_junction_spanning` | `True` if the best alignment spans the circular junction |
| `second_plasmid` | Name of the second-highest-scoring reference |
| `reason` | Short human-readable explanation of the assignment decision |
| `{plasmid}__match_frac` | Match fraction score for each reference |
| `{plasmid}__identity` | N-corrected alignment identity for each reference |
| `{plasmid}__query_cov` | Fraction of the read covered by the best alignment to each reference |
| `{plasmid}__mapq` | Mapping quality for the best alignment to each reference |
| `{plasmid}__n_in_hit` | Number of reference N bases within the aligned range (diagnostic) |

Plasmid name prefixes in column names have spaces, `/`, and `:` replaced with `_`.

### Per-plasmid summary TSV (`--summary`)

One row per reference, sorted by unique assigned read count.

| Column | Description |
|---|---|
| `plasmid` | Reference name |
| `assigned_reads` | Uniquely assigned reads (not ambiguous) |
| `ambiguous_reads` | Reads where this reference was the best hit but the read was ambiguous |
| `total_reads_best_hit` | Sum of the two above |
| `median_match_frac` | Median match fraction across all passing reads |
| `median_identity` | Median N-corrected identity across all passing reads |
| `median_query_cov` | Median query coverage across all passing reads |

---

## Visual report

The PDF report (`--plots-output`) contains the following pages:

**Page 1 — Read quality overview (pre-filter)**
Scatter plot of all reads (read length vs. average quality), coloured by filter status or final
plasmid assignment.  Marginal histograms on the top and right axes.

**Page 2 — Read disposition overview**
- Stacked bar chart showing raw → filtered → uniquely assigned read counts.
- Horizontal bar chart of unique and ambiguous read counts per plasmid.
- Read length distribution stacked by assignment.

**Page 3 — Alignment quality**
- Violin plot of match fraction distributions per plasmid (uniquely assigned reads only).
- Scatter plot of best vs. second-best match fraction per read (reads near the diagonal are
  ambiguous).
- Identity vs. query coverage scatter (uniquely assigned reads).
- Read length vs. match fraction scatter (uniquely assigned reads).

**Page 4 — Cross-alignment heatmap**
Median match fraction of reads assigned to one plasmid when aligned to every other plasmid.
High off-diagonal values indicate shared sequence (backbone) driving ambiguity.

**Page 5 — Per-plasmid match fraction distributions**
For each reference, overlaid density histograms comparing match fractions for reads assigned to
that reference vs. all other reads.  A well-separated distribution indicates reliable assignment.

**Pages 6+ — Per-plasmid pileup plots**
One page per reference with at least one assigned read.  Each page shows:
- A coverage depth histogram at the top, with N-region positions shaded grey.
- A packed read-track pileup below, coloured by match type:
  - Green bars — matching bases
  - Red ticks — mismatches
  - Dark grey bars — deletions
  - `I{n}` labels — insertions
  - Read border colour: blue (unique), orange (ambiguous)
  - Junction-spanning reads are labelled with purple text at their right edge: `0` if the
    read aligns to the end with no soft-clip, or `{N}bp` for the number of soft-clipped bases

---

## Design notes and limitations

**Aligner preset.** The `map-ont` minimap2 preset is used for all alignments.  It is tuned for
noisy ONT reads (high indel rate, 1–15 kbp).  It is not appropriate for Illumina short-read data.

**Reference size bias.** The match fraction denominator is `read_len` (not `max(read_len,
orig_ref_len)`).  This means that a read covering only the shared backbone of all references in
the library will receive the same score for each and will be called ambiguous.  This is the
correct behaviour — backbone-only reads genuinely cannot be assigned to a specific plasmid.

**N-region anchor gaps.** A 200 bp N run in a reference suppresses minimizer seeding across a
~400 bp window (the N run plus its surrounding k-mer shadow), which can cause circular aligners to
route around that region.  The N-aware identity correction mitigates the score effect, but it does
not restore coverage through the gap.  The `--fill-n` flag addresses this at the reference level
when biologically appropriate.

**Multi-sequence FASTA files.** Each `>header` record in a multi-sequence FASTA is treated as a
separate reference, named `{filestem}:{header}`.  Single-record FASTA files use `{filestem}` as
the reference name.

**Memory.** Read sequences for assigned reads are cached in memory for pileup plot generation.
For very large FASTQ files this can be significant.  Use `--no-plots` to skip caching and plot
generation.
