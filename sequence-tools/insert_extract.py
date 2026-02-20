#!/usr/bin/env python3
"""
insert_extract.py

Align ONT reads to a single plasmid reference, extract the subsequences
that fall within ambiguous N-base stretches, and align those inserts to
user-provided lookup references (FASTA or TSV).

Usage example
─────────────
    python3 insert_extract.py \
        --fastq 2462529_OL57GFP-Final-Scaleup.fastq \
        --reference pgl4-23mprav3gfp.fasta \
        --n-region 1:OL57_reference.fasta \
        --n-region 2:OL57.merged.match.enh.mapped.barcode.ct.parsed \
        --cross-check \
        [--min-qual 10] [--min-length 2000] [--threads 8] \
        [--output insert_results.tsv]

The --n-region flag maps the Nth stretch of Ns (1-based) in the reference
to a lookup file.  Lookup files can be:
  - FASTA (.fa/.fasta/.fna) — sequences are indexed and aligned with minimap2
  - TSV (anything else)     — column 1 = sequence, column 2 = ID

Requirements: mappy, pandas, numpy
"""

import argparse
import sys
import re
import gzip
import time
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import mappy
except ImportError:
    sys.exit("ERROR: mappy not found. Install with: pip install mappy")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract inserts from N regions of a plasmid reference "
                    "and align them to user-provided lookup files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--fastq", required=True,
                   help="ONT FASTQ file (plain or .gz)")
    p.add_argument("--reference", required=True,
                   help="Single plasmid FASTA file with N-base regions")
    p.add_argument("--n-region", action="append", required=True,
                   metavar="INDEX:FILE",
                   help="Map the Nth N-stretch (1-based) to a lookup file. "
                        "e.g. --n-region 1:OL57_reference.fasta "
                        "--n-region 2:barcodes.tsv. "
                        "Can be specified multiple times.")
    p.add_argument("--tsv-seq-col", type=int, default=1,
                   help="1-based column index for sequence in TSV lookup "
                        "(default: 1 = first column)")
    p.add_argument("--tsv-id-col", type=int, default=2,
                   help="1-based column index for ID in TSV lookup "
                        "(default: 2 = second column)")
    p.add_argument("--cross-check", action="store_true",
                   help="Evaluate whether IDs from different N-region "
                        "matches agree within the same read")
    p.add_argument("--min-qual", type=float, default=10.0,
                   help="Minimum average Phred quality score")
    p.add_argument("--min-length", type=int, default=2000,
                   help="Minimum read length (bp)")
    p.add_argument("--min-insert-match-frac", type=float, default=0.70,
                   help="Minimum fraction of insert bases matching to call "
                        "an insert aligned")
    p.add_argument("--max-edit-dist", type=int, default=-1, metavar="N",
                   help="Maximum Levenshtein edit distance for matching short "
                        "inserts (< 50 bp). -1 = auto (5%% of insert length, "
                        "min 1). 0 = exact match only (hash/minimap2). "
                        "Applies to TSV lookups and FASTA lookups on short "
                        "N-regions.")
    p.add_argument("--threads", type=int, default=4,
                   help="Threads for minimap2")
    p.add_argument("--circular", action="store_true",
                   help="Treat reference as circular (double sequence)")
    p.add_argument("--output", default="insert_results.tsv",
                   help="Output TSV file")
    p.add_argument("--plots-output", default="insert_report.pdf",
                   help="Output PDF file for plots")
    p.add_argument("--no-plots", action="store_true",
                   help="Skip plot generation")
    return p.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# FASTQ I/O
# ─────────────────────────────────────────────────────────────────────────────

def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def iter_fastq(path):
    """Yield (name, seq, qual_string) for every record."""
    with _open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip()
            fh.readline()          # +
            qual = fh.readline().rstrip()
            name = header[1:].split()[0]
            yield name, seq, qual


def avg_phred(qual_string):
    return float(np.mean([ord(c) - 33 for c in qual_string]))


# ─────────────────────────────────────────────────────────────────────────────
# FASTA loading
# ─────────────────────────────────────────────────────────────────────────────

def load_fasta(path):
    """Return list of (header, sequence) from a FASTA file."""
    records = []
    header, seqparts = None, []
    with _open(str(path)) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seqparts)))
                header = line[1:].split()[0]
                seqparts = []
            else:
                seqparts.append(line)
        if header is not None:
            records.append((header, "".join(seqparts)))
    return records


# ─────────────────────────────────────────────────────────────────────────────
# N-region detection
# ─────────────────────────────────────────────────────────────────────────────

def find_n_regions(seq):
    """Return list of (start, end) for contiguous N stretches in seq."""
    regions = []
    in_n = False
    n_start = 0
    for i, base in enumerate(seq.upper()):
        if base == 'N' and not in_n:
            in_n = True
            n_start = i
        elif base != 'N' and in_n:
            regions.append((n_start, i))
            in_n = False
    if in_n:
        regions.append((n_start, len(seq)))
    return regions


# ─────────────────────────────────────────────────────────────────────────────
# Lookup file loading
# ─────────────────────────────────────────────────────────────────────────────

def is_fasta_file(path):
    """Check if file looks like a FASTA by extension."""
    ext = Path(path).suffix.lower()
    return ext in ('.fa', '.fasta', '.fna')


_RC_TABLE = str.maketrans('ACGTacgt', 'TGCAtgca')


def revcomp(seq):
    """Return reverse complement of a DNA sequence."""
    return seq.translate(_RC_TABLE)[::-1]


_DNA_ALPHABET = ("A", "C", "G", "T")


def generate_variants(seq, max_dist):
    """
    Return a dict {variant_seq: min_edit_distance} for all sequences
    reachable from seq within Levenshtein distance max_dist.

    Uses BFS over the edit graph (each edge = one substitution, insertion,
    or deletion), so BFS depth == Levenshtein distance.
    """
    result = {seq: 0}
    if max_dist == 0:
        return result
    frontier = {seq}
    for dist in range(1, max_dist + 1):
        new_variants = {}
        for s in frontier:
            n = len(s)
            # Substitutions
            for i in range(n):
                for b in _DNA_ALPHABET:
                    if b != s[i]:
                        v = s[:i] + b + s[i + 1:]
                        if v not in result:
                            new_variants[v] = dist
            # Deletions
            for i in range(n):
                v = s[:i] + s[i + 1:]
                if v not in result:
                    new_variants[v] = dist
            # Insertions (including at the end)
            for i in range(n + 1):
                for b in _DNA_ALPHABET:
                    v = s[:i] + b + s[i:]
                    if v not in result:
                        new_variants[v] = dist
        result.update(new_variants)
        frontier = set(new_variants)
    return result


def match_insert_levenshtein(insert_seq, seq_dict, max_dist):
    """
    Match an insert against a sequence dict using Levenshtein neighbourhood
    lookup.  Generates all variants of the insert (and its RC) within
    max_dist edits and checks each against the hash dict.

    Returns (best_id, match_frac, strand) or (None, 0.0, None).
    match_frac = 1 - best_dist / insert_len (approximate).
    """
    if len(insert_seq) < 1:
        return None, 0.0, None

    insert_upper = insert_seq.upper()
    insert_rc = revcomp(insert_upper)
    insert_len = len(insert_upper)

    best_id = None
    best_dist = max_dist + 1
    best_strand = None

    for seq, strand in [(insert_upper, "+"), (insert_rc, "-")]:
        variants = generate_variants(seq, max_dist)   # {variant: dist}
        for variant, dist in variants.items():
            if dist < best_dist and variant in seq_dict:
                best_id = seq_dict[variant]
                best_dist = dist
                best_strand = strand

    if best_id is not None:
        match_frac = max(0.0, 1.0 - best_dist / max(1, insert_len))
        return best_id, round(match_frac, 4), best_strand
    return None, 0.0, None


def load_tsv_lookup(tsv_path, seq_col=1, id_col=2):
    """
    Load a TSV lookup file into a dict {sequence_upper: id}.
    Returns fwd_lookup dict only (RC is computed on-the-fly per query).
    """
    seq_idx = seq_col - 1
    id_idx = id_col - 1
    fwd_lookup = {}
    n_lines = 0
    with open(tsv_path) as fh:
        for line in fh:
            n_lines += 1
            if n_lines % 10_000_000 == 0:
                print(f"      … {n_lines:,} lines loaded …", flush=True)
            line = line.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) <= max(seq_idx, id_idx):
                continue
            seq = cols[seq_idx].upper()
            sid = cols[id_idx]
            if seq not in fwd_lookup:
                fwd_lookup[seq] = sid
    print(f"      {n_lines:,} total lines → {len(fwd_lookup):,} unique sequences",
          flush=True)
    return fwd_lookup


def build_fasta_aligner(path, n_threads, n_region_len=200):
    """
    Build a mappy aligner from a FASTA file.
    Chooses preset based on expected insert length.
    Returns aligner.
    """
    if n_region_len <= 100:
        aln = mappy.Aligner(path, preset="sr",
                            n_threads=n_threads, best_n=5)
    else:
        aln = mappy.Aligner(path, preset="map-ont",
                            n_threads=n_threads, best_n=5)
    if not aln:
        sys.exit(f"ERROR: failed to build index for {path}")
    return aln


def load_fasta_as_dict(path):
    """Load a FASTA file into a dict {sequence_upper: name} for hash lookup."""
    records = load_fasta(path)
    return {seq.upper(): name for name, seq in records}


# ─────────────────────────────────────────────────────────────────────────────
# cs-string parsing for extracting query subsequence at reference interval
# ─────────────────────────────────────────────────────────────────────────────

_CS_RE = re.compile(r':[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+')


def extract_insert_by_flanks(cs_string, r_st, q_st, interval_start,
                             interval_end):
    """
    Extract the query subsequence corresponding to the insert between two
    flanking regions on the reference.

    The reference has Ns at [interval_start, interval_end), but the actual
    insert in the read may be longer or shorter than the N stretch.
    minimap2 aligns the extra bases as insertions at the boundaries.

    Strategy: walk the cs string and track query positions.  The insert
    starts at the query position where the reference first enters the N
    region (including any insertions at that boundary) and ends where the
    reference leaves the N region.

    We find:
      - q_left:  the last query position consumed while ref < interval_start
                 (i.e. the query pos just after the last flank base)
      - q_right: the query position when ref first reaches interval_end

    Returns (query_start, query_end) or None.
    """
    ref_pos = r_st
    query_pos = q_st

    # Track the query position at the moment the ref crosses each boundary.
    # q_left_end: query pos after consuming the last token whose ref_pos
    #             is < interval_start (so the insert starts here)
    # q_right_start: query pos at the first token whose ref_pos >= interval_end
    q_left_end = None
    q_right_start = None

    tokens = _CS_RE.findall(cs_string)

    for token in tokens:
        ch = token[0]

        # At the start of each token, check if ref_pos has reached a
        # boundary.  This catches cases where the previous token ended
        # exactly at the boundary.
        if ref_pos == interval_start and q_left_end is None:
            q_left_end = query_pos
        if ref_pos == interval_end and q_right_start is None:
            q_right_start = query_pos

        if ch == ':':
            n = int(token[1:])
            block_r_st = ref_pos
            block_r_en = ref_pos + n
            block_q_st = query_pos

            # If this match block contains interval_start
            if block_r_st < interval_start < block_r_en:
                offset = interval_start - block_r_st
                q_left_end = block_q_st + offset

            # If this match block contains interval_end
            if block_r_st < interval_end < block_r_en:
                offset = interval_end - block_r_st
                q_right_start = block_q_st + offset

            ref_pos += n
            query_pos += n

        elif ch == '*':
            ref_pos += 1
            query_pos += 1

        elif ch == '+':
            ins_len = len(token) - 1
            query_pos += ins_len

        elif ch == '-':
            del_len = len(token) - 1

            # If deletion spans interval_start
            if ref_pos < interval_start < ref_pos + del_len:
                if q_left_end is None:
                    q_left_end = query_pos
            # If deletion spans interval_end
            if ref_pos < interval_end < ref_pos + del_len:
                if q_right_start is None:
                    q_right_start = query_pos

            ref_pos += del_len

    # Final boundary check after all tokens
    if ref_pos >= interval_start and q_left_end is None:
        q_left_end = query_pos
    if ref_pos >= interval_end and q_right_start is None:
        q_right_start = query_pos

    # If the alignment ended before reaching interval_end, q_right_start
    # stays None — the read doesn't span the full N region on the right.
    # Still return what we have if we at least entered the N region.
    if q_left_end is not None and q_right_start is None:
        # Alignment ended inside the N region — use query_pos at end
        q_right_start = query_pos

    if q_left_end is not None and q_right_start is not None:
        if q_left_end < q_right_start:
            return (q_left_end, q_right_start)
    return None


# ─────────────────────────────────────────────────────────────────────────────
# Insert alignment against lookups
# ─────────────────────────────────────────────────────────────────────────────

def align_insert_fasta(insert_seq, aligner, min_match_frac=0.70):
    """
    Align an insert sequence to a FASTA lookup using minimap2.
    Returns (best_target_name, match_frac, strand) or (None, 0, None).
    """
    insert_len = len(insert_seq)
    if insert_len < 10:
        return None, 0.0, None

    best_name = None
    best_mlen = 0
    best_strand = None

    hits = list(aligner.map(insert_seq, cs=True))
    for h in hits:
        if h.mlen > best_mlen:
            best_mlen = h.mlen
            best_name = h.ctg
            best_strand = "+" if h.strand == 1 else "-"

    match_frac = best_mlen / insert_len if insert_len > 0 else 0.0
    if match_frac >= min_match_frac:
        return best_name, round(match_frac, 4), best_strand
    return None, round(match_frac, 4), best_strand


def match_insert_tsv(insert_seq, fwd_lookup):
    """
    Match an insert to a TSV barcode lookup using exact or near-exact
    matching.  Checks both forward and reverse-complement of the insert.

    Strategy:
    1. Exact match of full insert (forward, then RC of insert)
    2. If insert is longer than barcode length, sliding window (fwd + RC)

    Returns (best_id, match_frac, strand) or (None, 0, None).
    """
    insert_len = len(insert_seq)
    if insert_len < 5:
        return None, 0.0, None

    insert_upper = insert_seq.upper()
    insert_rc = revcomp(insert_upper)

    # Determine expected barcode length from lookup
    sample_key = next(iter(fwd_lookup))
    bc_len = len(sample_key)

    # 1. Exact match — full insert
    if insert_upper in fwd_lookup:
        return fwd_lookup[insert_upper], 1.0, "+"
    if insert_rc in fwd_lookup:
        return fwd_lookup[insert_rc], 1.0, "-"

    # 2. Sliding window if insert is longer than barcode
    if insert_len > bc_len:
        for start in range(insert_len - bc_len + 1):
            subseq_fwd = insert_upper[start:start + bc_len]
            if subseq_fwd in fwd_lookup:
                frac = bc_len / insert_len
                return fwd_lookup[subseq_fwd], round(frac, 4), "+"
            subseq_rc = insert_rc[start:start + bc_len]
            if subseq_rc in fwd_lookup:
                frac = bc_len / insert_len
                return fwd_lookup[subseq_rc], round(frac, 4), "-"

    return None, 0.0, None


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

def generate_plots(results_df, n_regions, region_lookups, cross_check,
                   pdf_path, n_total, n_passed, n_mapped,
                   min_match_frac=0.70, region_lookup_types=None):
    """Generate a multi-page PDF report."""
    plt.rcParams.update({
        "font.size": 9, "axes.spines.top": False,
        "axes.spines.right": False, "figure.dpi": 150,
    })

    region_indices = sorted(region_lookups.keys())
    n_assigned = len(region_indices)
    cmap = plt.get_cmap("tab10")
    region_colors = {i: cmap(k % 10) for k, i in enumerate(region_indices)}

    # Pre-compute per-region insert and match counts
    region_n_insert = {}
    region_n_match = {}
    for idx in region_indices:
        col_len = f"n{idx}_insert_len"
        col_id = f"n{idx}_match_id"
        region_n_insert[idx] = int(
            (results_df[col_len] > 0).sum()
            if col_len in results_df.columns else 0)
        region_n_match[idx] = int(
            (results_df[col_id].notna() & (results_df[col_id] != "")).sum()
            if col_id in results_df.columns else 0)

    with PdfPages(pdf_path) as pdf:

        # ── PAGE 1: Processing funnel — stacked bar chart ────────────────
        # Each bar's total height = n_total.
        # Lower (colored) segment = reads passing that stage.
        # Upper (light gray) segment = dropout.
        # Stages: [All reads, QC passed, Mapped, per-region Insert, Match]
        funnel_stages = []   # (label, count, color)
        funnel_stages.append(("QC\npassed",    n_passed, "#555555"))
        funnel_stages.append(("Mapped\nto ref", n_mapped, "#333333"))
        for idx in region_indices:
            ns, ne = n_regions[idx - 1]
            n_len = ne - ns
            funnel_stages.append(
                (f"N{idx} insert\n({n_len} bp)", region_n_insert[idx],
                 region_colors[idx]))
            funnel_stages.append(
                (f"N{idx} match", region_n_match[idx],
                 region_colors[idx]))

        stage_labels  = [s[0] for s in funnel_stages]
        stage_counts  = [s[1] for s in funnel_stages]
        stage_colors  = [s[2] for s in funnel_stages]
        stage_dropout = [n_passed - c for c in stage_counts]

        fig, ax = plt.subplots(figsize=(max(8, len(funnel_stages) * 1.4 + 2), 6))
        x = np.arange(len(funnel_stages))

        ax.bar(x, stage_counts,  color=stage_colors,  edgecolor="white",
               width=0.6, alpha=0.90, label="Reads passing stage")
        ax.bar(x, stage_dropout, color="#eeeeee", edgecolor="white",
               width=0.6, bottom=stage_counts, label="Did not pass")

        # Annotate each bar: count + % of QC-passed reads
        for xi, (cnt, clr) in enumerate(zip(stage_counts, stage_colors)):
            pct = 100.0 * cnt / n_passed if n_passed > 0 else 0.0
            y_txt = cnt + n_passed * 0.01
            ax.text(xi, y_txt, f"{cnt:,}\n({pct:.1f}%)",
                    ha="center", va="bottom", fontsize=8, fontweight="bold",
                    color="#222222")

        ax.set_xticks(x)
        ax.set_xticklabels(stage_labels, fontsize=9)
        ax.set_ylabel("Number of reads")
        ax.set_ylim(0, n_passed * 1.22)
        ax.set_title("Read Processing Funnel", fontsize=13, fontweight="bold")
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ── PAGE 2: ID match + pairwise agreement ───────────────────────
        # Subset: reads with an ID match in at least one N-region.
        id_col = {idx: f"n{idx}_match_id" for idx in region_indices}
        any_match_mask = pd.Series(False, index=results_df.index)
        for col in id_col.values():
            if col in results_df.columns:
                any_match_mask |= (results_df[col].notna()
                                   & (results_df[col] != ""))
        subset = results_df[any_match_mask].copy()
        n_any = len(subset)

        if n_any > 0:
            # Compute n_match_i and pairwise_agree[i][j] for all pairs.
            n_match_per = {}
            pairwise_agree = {idx: {} for idx in region_indices}
            for idx in region_indices:
                col_i = id_col[idx]
                if col_i not in subset.columns:
                    n_match_per[idx] = 0
                    for jdx in region_indices:
                        if jdx != idx:
                            pairwise_agree[idx][jdx] = 0
                    continue
                has_i = subset[col_i] != ""
                n_match_per[idx] = int(has_i.sum())
                for jdx in region_indices:
                    if jdx == idx:
                        continue
                    col_j = id_col[jdx]
                    if col_j not in subset.columns:
                        pairwise_agree[idx][jdx] = 0
                        continue
                    has_j = subset[col_j] != ""
                    both  = has_i & has_j
                    same  = int((subset.loc[both, col_i]
                                 == subset.loc[both, col_j]).sum())
                    pairwise_agree[idx][jdx] = same

            # Layout: one group per region; each group has (1 + n_peers) bars.
            n_peers        = n_assigned - 1
            n_bars_per_grp = 1 + n_peers
            bw             = 0.25            # bar width
            grp_gap        = 0.35            # extra gap between groups
            grp_w          = n_bars_per_grp * bw + grp_gap

            all_x, all_h, all_c, all_a = [], [], [], []
            tick_bar_x, tick_bar_labels  = [], []   # per-bar x-axis labels
            grp_center_x, grp_labels     = [], []   # group separator labels
            max_h = 1

            for g, idx in enumerate(region_indices):
                ns, ne   = n_regions[idx - 1]
                n_len    = ne - ns
                grp_x0   = g * grp_w
                others   = [j for j in region_indices if j != idx]

                # Bar 0: match count for this region
                all_x.append(grp_x0)
                all_h.append(n_match_per[idx])
                all_c.append(region_colors[idx])
                all_a.append(0.90)
                tick_bar_x.append(grp_x0)
                tick_bar_labels.append("Has\nmatch")
                max_h = max(max_h, n_match_per[idx])

                # Bars 1…: pairwise agreement with each other region
                for k, jdx in enumerate(others):
                    agree = pairwise_agree[idx].get(jdx, 0)
                    all_x.append(grp_x0 + (k + 1) * bw)
                    all_h.append(agree)
                    all_c.append(region_colors[idx])
                    all_a.append(0.40)
                    tick_bar_x.append(grp_x0 + (k + 1) * bw)
                    tick_bar_labels.append(f"Agrees\nN{jdx}")
                    max_h = max(max_h, agree)

                grp_center_x.append(grp_x0 + (n_bars_per_grp - 1) * bw / 2)
                grp_labels.append(f"N{idx} ({n_len} bp)")

            fig, ax = plt.subplots(
                figsize=(max(8, len(all_x) * 1.1 + 2), 6))
            fig.suptitle(
                f"ID Matches and Region Agreement\n"
                f"(base: {n_any:,} reads with ≥1 ID match, "
                f"{100*n_any/n_passed:.1f}% of QC-passed)",
                fontsize=12, fontweight="bold")

            bars = ax.bar(all_x, all_h, width=bw * 0.9,
                          color=all_c, alpha=1.0, edgecolor="white")
            for bar, alpha in zip(bars, all_a):
                bar.set_alpha(alpha)

            # Annotate bars
            for xi, h, alpha, g_idx, bar_lbl in zip(
                    all_x, all_h, all_a,
                    [idx for idx in region_indices
                     for _ in range(n_bars_per_grp)],
                    tick_bar_labels):
                n_mi = n_match_per[g_idx]
                # "has match" bar: % of n_any; agree bar: % of n_match_i
                if "match" in bar_lbl.lower():
                    pct  = 100.0 * h / n_any if n_any > 0 else 0.0
                    frac = f"{h}/{n_any}\n({pct:.1f}%)"
                else:
                    pct  = 100.0 * h / n_mi if n_mi > 0 else 0.0
                    frac = f"{h}/{n_mi}\n({pct:.1f}%)"
                ax.text(xi + bw * 0.45, h + max_h * 0.01,
                        frac, ha="center", va="bottom",
                        fontsize=7.5, fontweight="bold" if alpha > 0.6 else "normal")

            # Per-bar x-ticks (small, rotated)
            ax.set_xticks(tick_bar_x)
            ax.set_xticklabels(tick_bar_labels, fontsize=8)

            # Group labels below the x-axis
            ax2 = ax.twiny()
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(grp_center_x)
            ax2.set_xticklabels(grp_labels, fontsize=10, fontweight="bold")
            ax2.tick_params(top=False, bottom=True, labeltop=False,
                            labelbottom=True, pad=28)
            ax2.spines["top"].set_visible(False)
            ax2.spines["bottom"].set_visible(False)

            # Faint vertical separators between groups
            for g in range(1, n_assigned):
                sep_x = g * grp_w - grp_gap / 2
                ax.axvline(sep_x, color="#cccccc", lw=0.8, ls="--")

            ax.set_ylabel("Number of reads")
            ax.set_ylim(0, max_h * 1.28)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
        else:
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.text(0.5, 0.5, "No ID matches found in any N-region.",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=12, color="#888888")
            ax.axis("off")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # ── PAGE 3: alignment score distributions ────────────────────────
        fig, axes_arr = plt.subplots(1, n_assigned,
                                     figsize=(6 * n_assigned, 5),
                                     squeeze=False)
        fig.suptitle("Insert Alignment Score Distributions", fontsize=13,
                     fontweight="bold")

        for k, idx in enumerate(region_indices):
            ax = axes_arr[0][k]
            ns, ne   = n_regions[idx - 1]
            n_len    = ne - ns
            col_len  = f"n{idx}_insert_len"
            col_frac = f"n{idx}_match_frac"
            col_id   = f"n{idx}_match_id"
            ltype    = (region_lookup_types or {}).get(idx, "unknown")

            if col_frac not in results_df.columns or col_len not in results_df.columns:
                ax.set_visible(False)
                continue

            extracted = results_df[col_len] > 0
            fracs     = results_df.loc[extracted, col_frac].fillna(0.0)
            n_extracted = int(extracted.sum())

            matched_mask = (results_df[col_id].notna()
                            & (results_df[col_id] != ""))
            n_matched    = int((extracted & matched_mask).sum())
            n_no_align   = int((extracted & (results_df[col_frac] == 0.0)).sum())

            if fracs.empty:
                ax.text(0.5, 0.5, "No inserts", transform=ax.transAxes,
                        ha="center", va="center")
                ax.set_title(f"N-region {idx} ({n_len} bp Ns)")
                continue

            ltype_label = ltype.replace("_", " ")
            title_base  = (
                f"N{idx} ({n_len} bp Ns) — {ltype_label}\n"
                f"Matched: {n_matched}  |  No alignment: {n_no_align}  |  "
                f"Extracted: {n_extracted}"
            )

            if ltype == "tsv":
                # ── Discrete edit-distance bar chart ──────────────────────
                insert_lens = results_df.loc[extracted, col_len]
                # Reconstruct edit distance from stored match_frac and length
                # edit_dist = round((1 - frac) * insert_len), 0 for exact match
                edit_dists = []
                for frac_val, ins_len in zip(fracs, insert_lens):
                    if frac_val == 0.0 and not results_df.loc[
                            fracs.index[fracs == frac_val].intersection(
                                insert_lens.index)].empty:
                        # frac==0 could mean no-match (id is NaN) or dist=0
                        pass  # handled below via matched_mask
                    ed = round((1.0 - frac_val) * max(1, ins_len))
                    edit_dists.append(ed)

                # Rebuild properly: iterate rows so we pair correctly
                edit_dists = []
                for row_idx in results_df.index[extracted]:
                    fv  = results_df.at[row_idx, col_frac]
                    il  = results_df.at[row_idx, col_len]
                    mid = results_df.at[row_idx, col_id]
                    if pd.isna(mid) or mid == "":
                        edit_dists.append(None)   # no match
                    else:
                        ed = int(round((1.0 - fv) * max(1, il)))
                        edit_dists.append(ed)

                max_ed = max((e for e in edit_dists if e is not None), default=0)
                ed_counts = {e: 0 for e in range(max_ed + 1)}
                no_match_count = 0
                for e in edit_dists:
                    if e is None:
                        no_match_count += 1
                    else:
                        ed_counts[e] = ed_counts.get(e, 0) + 1

                ed_labels  = [str(e) for e in range(max_ed + 1)] + ["no match"]
                ed_values  = [ed_counts[e] for e in range(max_ed + 1)] + [no_match_count]
                bar_colors = [region_colors[idx]] * (max_ed + 1) + ["#aaaaaa"]
                x_pos      = np.arange(len(ed_labels))

                ax.bar(x_pos, ed_values, color=bar_colors, edgecolor="white",
                       alpha=0.85)
                ax.set_xticks(x_pos)
                ax.set_xticklabels(ed_labels, fontsize=9)
                ax.set_xlabel("Edit distance (Levenshtein)")
                ax.set_ylabel("Count")
                ax.set_title(title_base, fontsize=9)

            else:
                # ── Continuous match-fraction histogram ────────────────────
                bins = np.linspace(0, 1, 42)  # 41 bins of 0.025 width
                ax.hist(fracs, bins=bins, color=region_colors[idx],
                        alpha=0.75, edgecolor="white")
                ax.axvline(min_match_frac, ls="--", color="#d62728", lw=1.5,
                           label=f"Cutoff ({min_match_frac:.2f})")
                ax.set_xlabel("Match fraction (best alignment)")
                ax.set_ylabel("Count")
                ax.set_xlim(0, 1)
                ax.legend(fontsize=7)
                ax.set_title(title_base, fontsize=9)

        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ── PAGE 4: insert length distributions ──────────────────────────
        fig, axes_arr = plt.subplots(1, n_assigned, figsize=(6 * n_assigned, 5),
                                     squeeze=False)
        fig.suptitle("Insert Length Distributions", fontsize=13,
                     fontweight="bold")

        for k, idx in enumerate(region_indices):
            ax = axes_arr[0][k]
            ns, ne = n_regions[idx - 1]
            n_len = ne - ns
            col_len = f"n{idx}_insert_len"
            col_id = f"n{idx}_match_id"

            if col_len not in results_df.columns:
                ax.set_visible(False)
                continue

            lengths_all = results_df.loc[results_df[col_len] > 0, col_len]
            matched_mask = results_df[col_id].notna() & (results_df[col_id] != "")
            lengths_matched = results_df.loc[
                (results_df[col_len] > 0) & matched_mask, col_len]
            lengths_unmatched = results_df.loc[
                (results_df[col_len] > 0) & ~matched_mask, col_len]

            if lengths_all.empty:
                ax.text(0.5, 0.5, "No inserts", transform=ax.transAxes,
                        ha="center", va="center")
                ax.set_title(f"N-region {idx} ({n_len} bp Ns)")
                continue

            lo = max(0, lengths_all.min() - 5)
            hi = lengths_all.max() + 5
            bins = np.arange(lo, hi + 1, max(1, (hi - lo) // 40))

            if not lengths_matched.empty:
                ax.hist(lengths_matched, bins=bins,
                        color=region_colors[idx], alpha=0.7,
                        label=f"Matched (n={len(lengths_matched)})",
                        edgecolor="white")
            if not lengths_unmatched.empty:
                ax.hist(lengths_unmatched, bins=bins,
                        color="#cccccc", alpha=0.6,
                        label=f"Unmatched (n={len(lengths_unmatched)})",
                        edgecolor="white")

            ax.axvline(n_len, ls="--", color="red", lw=1.2,
                       label=f"N-region size ({n_len} bp)")
            if not lengths_all.empty:
                ax.axvline(lengths_all.mean(), ls="--", color="blue", lw=1.0,
                           label=f"Mean ({lengths_all.mean():.1f} bp)")

            ax.set_xlabel("Insert length (bp)")
            ax.set_ylabel("Count")
            ax.set_title(f"N-region {idx} ({n_len} bp Ns)")
            ax.legend(fontsize=7)

        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ── PAGE 5: UpSet-style intersection plot ─────────────────────
        # Categories per read: has insert from each N-region, has match
        # from each N-region.  Show all non-empty intersections.

        # Build boolean columns
        cat_labels = []
        cat_masks = []

        for idx in region_indices:
            ns, ne = n_regions[idx - 1]
            n_len = ne - ns
            col_len = f"n{idx}_insert_len"
            col_id = f"n{idx}_match_id"

            if col_len in results_df.columns:
                ins_mask = (results_df[col_len] > 0).values
                cat_labels.append(f"N{idx} insert\n({n_len}bp)")
                cat_masks.append(ins_mask)

            if col_id in results_df.columns:
                match_mask = (results_df[col_id].notna()
                              & (results_df[col_id] != "")).values
                cat_labels.append(f"N{idx} ID match")
                cat_masks.append(match_mask)

        n_cats = len(cat_labels)
        n_reads = len(results_df)

        # Enumerate every unique combination present in the data
        combo_key = []
        for i in range(n_reads):
            combo_key.append(tuple(bool(m[i]) for m in cat_masks))

        from collections import Counter
        combo_counts = Counter(combo_key)

        # Remove the all-False combination (mapped but no insert at all)
        all_false = tuple(False for _ in range(n_cats))

        # Sort combos: by number of set bits descending, then count desc
        combos_sorted = sorted(
            combo_counts.items(),
            key=lambda x: (sum(x[0]), x[1]),
            reverse=True,
        )

        # Filter out empty intersection if present, keep it for display
        combo_labels = []
        combo_vals = []
        combo_bits = []
        for bits, count in combos_sorted:
            combo_labels.append(
                " + ".join(cat_labels[j].replace("\n", " ")
                           for j in range(n_cats) if bits[j])
                if any(bits) else "No insert"
            )
            combo_vals.append(count)
            combo_bits.append(bits)

        n_combos = len(combo_vals)

        # UpSet layout: bottom = dot matrix, top = bar chart
        fig_h = max(5, 2.5 + n_cats * 0.4 + 0.3)
        fig = plt.figure(figsize=(max(8, n_combos * 0.8 + 2), fig_h))

        # Grid: bar chart on top (3 parts), dot matrix on bottom (n_cats parts)
        gs = fig.add_gridspec(nrows=2, ncols=1,
                              height_ratios=[3, max(1.5, n_cats * 0.5)],
                              hspace=0.05)
        ax_bar = fig.add_subplot(gs[0])
        ax_dot = fig.add_subplot(gs[1], sharex=ax_bar)

        x = np.arange(n_combos)

        # Color bars by whether all categories are set
        bar_colors = []
        for bits in combo_bits:
            if all(bits):
                bar_colors.append("#2166ac")
            elif any(bits):
                bar_colors.append("#74add1")
            else:
                bar_colors.append("#cccccc")

        ax_bar.bar(x, combo_vals, color=bar_colors, edgecolor="white",
                   width=0.6)
        for xi, val in zip(x, combo_vals):
            ax_bar.text(xi, val + max(combo_vals) * 0.015, str(val),
                        ha="center", va="bottom", fontsize=8,
                        fontweight="bold")
        ax_bar.set_ylabel("Number of reads")
        ax_bar.set_title("Read Classification — UpSet Plot\n"
                         f"(all {n_reads} reads mapped to reference)",
                         fontsize=12, fontweight="bold")
        ax_bar.set_xlim(-0.6, n_combos - 0.4)
        ax_bar.tick_params(bottom=False, labelbottom=False)

        # Dot matrix
        for row_idx in range(n_cats):
            for col_idx in range(n_combos):
                is_set = combo_bits[col_idx][row_idx]
                color = "#333333" if is_set else "#dddddd"
                size = 80 if is_set else 30
                ax_dot.scatter(col_idx, row_idx, s=size, c=color,
                               zorder=3, edgecolors="none")

            # Draw connecting lines for set dots in each column
        for col_idx in range(n_combos):
            set_rows = [r for r in range(n_cats)
                        if combo_bits[col_idx][r]]
            if len(set_rows) >= 2:
                ax_dot.plot([col_idx, col_idx],
                            [min(set_rows), max(set_rows)],
                            color="#333333", lw=1.5, zorder=2)

        ax_dot.set_yticks(range(n_cats))
        ax_dot.set_yticklabels(cat_labels, fontsize=8)
        ax_dot.set_ylim(-0.5, n_cats - 0.5)
        ax_dot.invert_yaxis()
        ax_dot.tick_params(bottom=False, labelbottom=False)
        ax_dot.set_xlim(-0.6, n_combos - 0.4)

        # Light grid
        for row_idx in range(n_cats):
            ax_dot.axhline(row_idx, color="#eeeeee", lw=0.8, zorder=1)

        ax_dot.spines["left"].set_visible(False)
        ax_dot.spines["bottom"].set_visible(False)

        # Set size bars on the right side of the dot matrix
        # (total count per category)
        cat_totals = [int(m.sum()) for m in cat_masks]
        for row_idx, total in enumerate(cat_totals):
            ax_dot.text(n_combos - 0.2, row_idx,
                        f"  n={total:,}", va="center", fontsize=7,
                        color="#555555")

        fig.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.05)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"  Plots written to: {pdf_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # ── validate inputs ──────────────────────────────────────────────────────
    if not Path(args.fastq).exists():
        sys.exit(f"ERROR: FASTQ file not found: {args.fastq}")
    if not Path(args.reference).exists():
        sys.exit(f"ERROR: Reference file not found: {args.reference}")

    # ── load reference ───────────────────────────────────────────────────────
    ref_records = load_fasta(args.reference)
    if not ref_records:
        sys.exit("ERROR: no sequences in reference FASTA")
    ref_name, ref_seq = ref_records[0]
    orig_ref_len = len(ref_seq)

    # ── find N regions ───────────────────────────────────────────────────────
    n_regions = find_n_regions(ref_seq)
    if not n_regions:
        sys.exit("ERROR: no N-base stretches found in reference")

    # ── parse --n-region arguments ───────────────────────────────────────────
    region_lookups = {}  # {1-based index: file path}
    for spec in args.n_region:
        parts = spec.split(":", 1)
        if len(parts) != 2:
            sys.exit(f"ERROR: --n-region must be INDEX:FILE, got: {spec}")
        idx = int(parts[0])
        fpath = parts[1]
        if not Path(fpath).exists():
            sys.exit(f"ERROR: lookup file not found: {fpath}")
        if idx < 1 or idx > len(n_regions):
            sys.exit(f"ERROR: N-region index {idx} out of range "
                     f"(reference has {len(n_regions)} N-stretches)")
        region_lookups[idx] = fpath

    # ── print header ─────────────────────────────────────────────────────────
    print("=" * 70)
    print("  insert_extract.py")
    print("=" * 70)
    print(f"  Reference:      {args.reference}")
    print(f"  Sequence:       {ref_name} ({orig_ref_len:,} bp)")
    print(f"  FASTQ:          {args.fastq}")
    print(f"  Circular mode:  {args.circular}")
    print(f"  Min quality:    Q{args.min_qual}")
    print(f"  Min length:     {args.min_length} bp")
    print(f"  Threads:        {args.threads}")
    _edit_label = ("auto (5% of insert len)" if args.max_edit_dist == -1
                   else ("disabled" if args.max_edit_dist == 0
                         else str(args.max_edit_dist)))
    print(f"  Max edit dist:  {_edit_label} (for inserts < 50 bp)")
    print()

    print(f"  N-base stretches found: {len(n_regions)}")
    print(f"  {'#':<4} {'Start':>8} {'End':>8} {'Length':>8}  Lookup file")
    print(f"  {'─'*4} {'─'*8} {'─'*8} {'─'*8}  {'─'*40}")
    for i, (ns, ne) in enumerate(n_regions, 1):
        n_len = ne - ns
        lfile = region_lookups.get(i, "(not assigned)")
        ftype = ""
        if i in region_lookups:
            ftype = " [FASTA]" if is_fasta_file(region_lookups[i]) else " [TSV]"
        print(f"  {i:<4} {ns:>8,} {ne:>8,} {n_len:>8,}  {lfile}{ftype}")
    print()

    # ── handle circular mode ─────────────────────────────────────────────────
    if args.circular:
        align_seq = ref_seq + ref_seq
    else:
        align_seq = ref_seq

    # ── build reference aligner ──────────────────────────────────────────────
    print("Building reference aligner …")
    ref_aligner = mappy.Aligner(seq=align_seq, preset="map-ont",
                                n_threads=args.threads, best_n=5)
    if not ref_aligner:
        sys.exit("ERROR: failed to build reference aligner")

    # ── load lookup files ────────────────────────────────────────────────────
    # lookup_data[idx]:
    #   ("fasta",       aligner)    — minimap2 for long inserts (>= 50 bp)
    #   ("fasta_short", seq_dict)   — hash+Levenshtein for short FASTA inserts
    #   ("tsv",         fwd_lookup) — hash+Levenshtein or exact for TSV
    lookup_data = {}
    for idx, fpath in region_lookups.items():
        n_start, n_end = n_regions[idx - 1]
        n_len = n_end - n_start
        print(f"  Loading lookup for N-region {idx}: {fpath} …")
        if is_fasta_file(fpath):
            if n_len < 50:
                # minimap2 can't seed reads shorter than its k-mer; use hash
                seq_dict = load_fasta_as_dict(fpath)
                print(f"    FASTA (short N-region): {len(seq_dict):,} sequences "
                      f"indexed (hash + Levenshtein, preset=map-ont)")
                lookup_data[idx] = ("fasta_short", seq_dict)
            else:
                aln = build_fasta_aligner(fpath, args.threads,
                                          n_region_len=n_len)
                records = load_fasta(fpath)
                print(f"    FASTA: {len(records):,} sequences indexed (preset=map-ont)")
                lookup_data[idx] = ("fasta", aln)
        else:
            print(f"    Loading TSV barcode lookup …")
            fwd_lk = load_tsv_lookup(fpath, args.tsv_seq_col,
                                     args.tsv_id_col)
            method = ("hash + Levenshtein" if args.max_edit_dist != 0
                      else "hash exact")
            print(f"    TSV: {len(fwd_lk):,} unique barcodes ({method})")
            lookup_data[idx] = ("tsv", fwd_lk)
    print()

    # ── stream reads and process ─────────────────────────────────────────────
    print("Processing reads …")
    t0 = time.time()

    rows = []
    n_total = 0
    n_passed = 0
    n_mapped = 0
    insert_lengths = {i: [] for i in region_lookups}
    match_counts = {i: defaultdict(int) for i in region_lookups}
    n_matched = {i: 0 for i in region_lookups}

    for read_name, read_seq, read_qual in iter_fastq(args.fastq):
        n_total += 1
        if len(read_seq) < args.min_length:
            continue
        q = avg_phred(read_qual)
        if q < args.min_qual:
            continue
        n_passed += 1

        if n_passed % 500 == 0:
            elapsed = time.time() - t0
            print(f"  … {n_passed:,} reads processed ({elapsed:.0f}s)",
                  flush=True)

        # Align read to reference
        hits = list(ref_aligner.map(read_seq, cs=True))
        if not hits:
            continue

        # Take the best hit by mlen
        best_hit = max(hits, key=lambda h: h.mlen)
        if not hasattr(best_hit, 'cs') or best_hit.cs is None:
            continue

        n_mapped += 1
        read_len = len(read_seq)

        # For circular references, wrap coordinates
        r_st = best_hit.r_st
        r_en = best_hit.r_en

        row = {
            "read_id": read_name,
            "read_length": read_len,
            "avg_qual": round(q, 2),
            "ref_start": r_st % orig_ref_len if args.circular else r_st,
            "ref_end": r_en % orig_ref_len if args.circular else r_en,
            "strand": "+" if best_hit.strand == 1 else "-",
            "match_frac": round(best_hit.mlen / read_len, 4),
        }

        # Extract inserts for each N region
        for region_idx in region_lookups:
            n_start, n_end = n_regions[region_idx - 1]
            n_len = n_end - n_start

            # Also check in second copy for circular
            intervals_to_check = [(n_start, n_end)]
            if args.circular:
                intervals_to_check.append(
                    (n_start + orig_ref_len, n_end + orig_ref_len))

            insert_seq = None
            for iv_start, iv_end in intervals_to_check:
                result = extract_insert_by_flanks(
                    best_hit.cs, best_hit.r_st, best_hit.q_st,
                    iv_start, iv_end)
                if result is not None:
                    qs, qe = result
                    if best_hit.strand == -1:
                        # CS string walks the revcomp of the query, so the
                        # walk-space positions must be mirrored to get the
                        # correct original-read coordinates.
                        orig_qs = best_hit.q_en + best_hit.q_st - qe
                        orig_qe = best_hit.q_en + best_hit.q_st - qs
                        qs, qe = orig_qs, orig_qe
                    if 0 <= qs < qe <= read_len:
                        insert_seq = read_seq[qs:qe]
                        break

            prefix = f"n{region_idx}"
            if insert_seq is None or len(insert_seq) < 5:
                row[f"{prefix}_insert_seq"] = ""
                row[f"{prefix}_insert_len"] = 0
                row[f"{prefix}_match_id"] = ""
                row[f"{prefix}_match_frac"] = 0.0
                row[f"{prefix}_match_strand"] = ""
                continue

            insert_lengths[region_idx].append(len(insert_seq))
            row[f"{prefix}_insert_seq"] = insert_seq
            row[f"{prefix}_insert_len"] = len(insert_seq)

            # Align insert to lookup
            ltype, lref = lookup_data[region_idx]
            insert_len = len(insert_seq)

            # Compute effective edit distance for this insert
            if args.max_edit_dist == -1:   # auto
                if insert_len < 50:
                    eff_edit_dist = max(1, round(0.05 * insert_len))
                else:
                    eff_edit_dist = 0
            else:
                eff_edit_dist = args.max_edit_dist

            if ltype == "fasta":
                match_id, mfrac, mstrand = align_insert_fasta(
                    insert_seq, lref, args.min_insert_match_frac)
            elif ltype == "fasta_short":
                match_id, mfrac, mstrand = match_insert_levenshtein(
                    insert_seq, lref, eff_edit_dist)
            else:  # tsv
                if eff_edit_dist > 0:
                    match_id, mfrac, mstrand = match_insert_levenshtein(
                        insert_seq, lref, eff_edit_dist)
                else:
                    match_id, mfrac, mstrand = match_insert_tsv(
                        insert_seq, lref)

            row[f"{prefix}_match_id"] = match_id if match_id else ""
            row[f"{prefix}_match_frac"] = mfrac
            row[f"{prefix}_match_strand"] = mstrand if mstrand else ""

            if match_id:
                n_matched[region_idx] += 1
                match_counts[region_idx][match_id] += 1

        rows.append(row)

    elapsed = time.time() - t0

    # ── write results ────────────────────────────────────────────────────────
    results_df = pd.DataFrame(rows)
    results_df.to_csv(args.output, sep="\t", index=False)
    print(f"\nResults written to: {args.output}")

    # ── summary ──────────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"  Total reads:          {n_total:>10,}")
    print(f"  Passed filters:       {n_passed:>10,}")
    print(f"  Mapped to reference:  {n_mapped:>10,}")
    print(f"  Processing time:      {elapsed:>10.1f}s")
    print()

    for region_idx in sorted(region_lookups.keys()):
        n_start, n_end = n_regions[region_idx - 1]
        n_len = n_end - n_start
        lengths = insert_lengths[region_idx]
        counts = match_counts[region_idx]

        print(f"  N-region {region_idx} "
              f"(pos {n_start:,}–{n_end:,}, {n_len} bp Ns):")
        print(f"    Lookup file:        {region_lookups[region_idx]}")
        if lengths:
            print(f"    Inserts extracted:  {len(lengths):,}")
            print(f"    Avg insert length:  {np.mean(lengths):.1f} bp")
            print(f"    Median insert len:  {np.median(lengths):.1f} bp")
            print(f"    Min/Max insert len: {min(lengths)}/{max(lengths)} bp")
        else:
            print(f"    Inserts extracted:  0")
        print(f"    Matched to lookup:  {n_matched[region_idx]:,}")
        if counts:
            print(f"    Unique IDs matched: {len(counts):,}")
            top_n = sorted(counts.items(), key=lambda x: x[1], reverse=True)
            print(f"    Top 10 matches:")
            for sid, ct in top_n[:10]:
                print(f"      {sid:50s}  {ct:>6,}")
        print()

    # ── cross-check ──────────────────────────────────────────────────────────
    if args.cross_check and len(region_lookups) >= 2:
        print("=" * 70)
        print("  CROSS-CHECK: ID agreement across N-regions within reads")
        print("=" * 70)

        region_indices = sorted(region_lookups.keys())
        id_cols = [f"n{i}_match_id" for i in region_indices]

        # Only consider reads where all regions have a match
        has_all = results_df.copy()
        for col in id_cols:
            if col in has_all.columns:
                has_all = has_all[has_all[col] != ""]

        if len(has_all) == 0:
            print("  No reads with matches in all assigned N-regions.")
        else:
            print(f"  Reads with matches in all {len(region_indices)} "
                  f"N-regions: {len(has_all):,}")

            # Pairwise comparison
            for i in range(len(region_indices)):
                for j in range(i + 1, len(region_indices)):
                    ri = region_indices[i]
                    rj = region_indices[j]
                    col_i = f"n{ri}_match_id"
                    col_j = f"n{rj}_match_id"

                    if col_i not in has_all.columns or col_j not in has_all.columns:
                        continue

                    same = (has_all[col_i] == has_all[col_j]).sum()
                    diff = len(has_all) - same
                    pct_same = 100 * same / len(has_all) if len(has_all) > 0 else 0

                    print(f"\n  N-region {ri} vs N-region {rj}:")
                    print(f"    Same ID:       {same:>8,} ({pct_same:.1f}%)")
                    print(f"    Different ID:  {diff:>8,} ({100-pct_same:.1f}%)")

                    # Show top mismatched pairs
                    mismatched = has_all[has_all[col_i] != has_all[col_j]]
                    if len(mismatched) > 0:
                        pair_counts = (mismatched
                                       .groupby([col_i, col_j])
                                       .size()
                                       .reset_index(name="count")
                                       .sort_values("count", ascending=False))
                        print(f"    Top mismatched ID pairs:")
                        for _, prow in pair_counts.head(10).iterrows():
                            print(f"      {prow[col_i]:30s} | "
                                  f"{prow[col_j]:30s}  {prow['count']:>5,}")

        print()

    # ── plots ─────────────────────────────────────────────────────────────────
    if not args.no_plots:
        print("Generating plots …")
        region_lookup_types = {idx: ltype
                               for idx, (ltype, _) in lookup_data.items()}
        generate_plots(results_df, n_regions, region_lookups,
                       args.cross_check, args.plots_output,
                       n_total=n_total, n_passed=n_passed, n_mapped=n_mapped,
                       min_match_frac=args.min_insert_match_frac,
                       region_lookup_types=region_lookup_types)

    print("=" * 70)
    print("  Done.")
    print("=" * 70)


if __name__ == "__main__":
    main()
