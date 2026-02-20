#!/usr/bin/env python3
"""
plasmid_align.py

Align Oxford Nanopore reads to a library of plasmid FASTA sequences.
For each read, determines the most likely plasmid assignment or flags
the read as ambiguous when two references score similarly.

Plasmid references may contain stretches of N bases; minimap2 (via mappy)
handles these natively — N positions in the reference are treated as
wildcards and do not penalise the alignment score.

Usage
-----
    python3 plasmid_align.py \\
        --fastq 2462529_OL57GFP-Final-Scaleup.fastq \\
        --fasta-dir . \\
        [--min-qual 10] [--min-length 2000] \\
        [--ambiguity-delta 0.05] [--circular] \\
        [--threads 8] [--output results.tsv]

Requirements (conda run -n sandbox python3 plasmid_align.py ...)
    mappy, pandas, numpy, matplotlib
"""

import argparse
import sys
import os
import re
import gzip
import time
from pathlib import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_rgba

try:
    import mappy
except ImportError:
    sys.exit("ERROR: mappy not found. Run inside the sandbox conda env:\n"
             "  conda run -n sandbox python3 plasmid_align.py ...")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Align ONT reads to a plasmid library and assign each "
                    "read to its most likely plasmid.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--fastq", required=True,
                   help="ONT FASTQ file (plain or .gz)")
    p.add_argument("--fasta", required=True, nargs="+", metavar="FASTA",
                   help="One or more plasmid FASTA files. Each file may "
                        "contain multiple sequences.")
    p.add_argument("--min-qual", type=float, default=10.0,
                   help="Minimum average Phred quality score to keep a read")
    p.add_argument("--min-length", type=int, default=2000,
                   help="Minimum read length (bp) to keep a read")
    p.add_argument("--ambiguity-delta", type=float, default=0.05,
                   help="If the top two plasmid scores differ by less than "
                        "this fraction (of the best score), call the read "
                        "ambiguous. Range 0–1.")
    p.add_argument("--circular", action="store_true",
                   help="Treat plasmid references as circular. Each "
                        "reference sequence is internally doubled so that "
                        "reads spanning the linearisation junction can align.")
    p.add_argument("--threads", type=int, default=4,
                   help="Number of threads passed to minimap2")
    p.add_argument("--output", default="plasmid_assignments.tsv",
                   help="Output TSV file for per-read results")
    p.add_argument("--summary", default="plasmid_summary.tsv",
                   help="Output TSV file for per-plasmid summary")
    p.add_argument("--plots-output", default="plasmid_alignment_report.pdf",
                   help="Output PDF file for visual report")
    p.add_argument("--no-plots", action="store_true",
                   help="Skip plot generation")
    p.add_argument("--fill-n", action="store_true", default=False,
                   help="Fill N-masked regions in each reference using the most "
                        "homologous sequence from another reference in the library. "
                        "Improves alignment coverage through masked regions. "
                        "Verify fills are biologically correct before use in production.")
    return p.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# FASTQ I/O + filtering
# ─────────────────────────────────────────────────────────────────────────────

def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def iter_fastq(path):
    """Yield (name, seq, qual_string) for every record in a FASTQ file."""
    with _open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().rstrip()
            fh.readline()          # +
            qual = fh.readline().rstrip()
            name = header[1:].split()[0]   # strip @ and description
            yield name, seq, qual


def avg_phred(qual_string):
    """Mean Phred Q score from an ASCII quality string."""
    return float(np.frombuffer(qual_string.encode('ascii'), dtype=np.uint8).mean() - 33)




# ─────────────────────────────────────────────────────────────────────────────
# FASTA loading
# ─────────────────────────────────────────────────────────────────────────────

def load_fasta(path):
    """Return a list of (header, sequence) tuples from a FASTA file."""
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



def load_plasmid_library(fasta_paths, circular=False):
    """
    Load plasmid sequences from a list of FASTA file paths.
    Returns (plasmids, orig_lens, n_cumsums) where:
        plasmids  = {plasmid_name: sequence_string}
        orig_lens = {plasmid_name: original_length_before_doubling}
        n_cumsums = {plasmid_name: cumulative-N-count array of length orig_len+1}
                    n_cumsums[name][b] - n_cumsums[name][a] = N count in [a, b)

    If circular=True each sequence is doubled (concatenated with itself)
    so that reads spanning the linearisation junction can align.
    """
    plasmids  = {}
    orig_lens = {}
    n_cumsums = {}
    for fpath in (Path(p) for p in fasta_paths):
        records = load_fasta(fpath)
        for header, seq in records:
            name = f"{fpath.stem}:{header}" if len(records) > 1 else fpath.stem
            if not seq:
                print(f"  WARNING: {name} has empty sequence, skipping.")
                continue
            orig_len = len(seq)
            orig_lens[name] = orig_len
            # Pre-compute cumulative N count on the original (non-doubled) sequence
            cs = np.zeros(orig_len + 1, dtype=np.int32)
            cs[1:] = np.cumsum(
                np.frombuffer(seq.upper().encode(), dtype=np.uint8) == ord('N'))
            n_cumsums[name] = cs
            if circular:
                seq = seq + seq          # double for circular alignment
            plasmids[name] = seq
    return plasmids, orig_lens, n_cumsums


def _n_in_range(cumsum, a, b, orig_len):
    """
    Count N bases in [a, b) where coordinates are on the doubled reference
    (length 2*orig_len) and cumsum is built on the original (length orig_len).
    Handles all three cases: within first copy, within second copy, and spanning
    the junction.
    """
    # Normalise coordinates that fall entirely in the second copy
    if a >= orig_len:
        a -= orig_len
        b -= orig_len
    if b <= orig_len:
        return int(cumsum[b] - cumsum[a])
    # Spans junction: [a, orig_len) in copy 1 + [0, b-orig_len) in copy 2
    part1 = int(cumsum[orig_len] - cumsum[a])
    part2 = int(cumsum[min(b - orig_len, orig_len)] - cumsum[0])
    return part1 + part2


# ─────────────────────────────────────────────────────────────────────────────
# Alignment
# ─────────────────────────────────────────────────────────────────────────────

def build_aligners(plasmid_seqs, n_threads, circular=False):
    """
    Build one mappy.Aligner per plasmid.
    preset="map-ont"  – tuned for noisy ONT long reads
    n_threads         – passed through for index building
    extra_flags       – 0 (default); N bases in reference are wildcards
    best_n            – 20 in circular mode to capture split junction hits
    """
    best_n = 20 if circular else 10
    aligners = {}
    for name, seq in plasmid_seqs.items():
        aln = mappy.Aligner(seq=seq, preset="map-ont", best_n=best_n,
                            n_threads=n_threads)
        if not aln:
            print(f"  WARNING: failed to build index for {name}, skipping.")
            continue
        aligners[name] = aln
    return aligners


def _find_junction_pair(hits, orig_len):
    """
    Find the best junction-spanning pair of hits on a doubled circular
    reference.  Returns (ha, hb) sorted by reference position, or None
    if no valid junction merge beats the best single hit.
    """
    if len(hits) < 2:
        return None

    best_single_mlen = max(h.mlen for h in hits)
    L = orig_len

    best_pair = None
    best_pair_mlen = 0

    for i in range(len(hits)):
        for j in range(i + 1, len(hits)):
            h1, h2 = hits[i], hits[j]

            if h1.strand != h2.strand:
                continue

            q_ov = max(0, min(h1.q_en, h2.q_en) - max(h1.q_st, h2.q_st))
            shorter_q = min(h1.q_en - h1.q_st, h2.q_en - h2.q_st)
            if shorter_q > 0 and q_ov / shorter_q > 0.20:
                continue

            ha, hb = (h1, h2) if h1.r_st < h2.r_st else (h2, h1)
            tol = max(100, L * 0.05)
            if not ((ha.r_en > L - tol and ha.r_st < L) and
                    (hb.r_st >= L - tol and hb.r_st < L + tol)):
                continue

            m = h1.mlen + h2.mlen
            if m > best_pair_mlen:
                best_pair_mlen = m
                best_pair = (ha, hb)

    if best_pair and best_pair_mlen > best_single_mlen:
        return best_pair
    return None


def merge_circular_hits(hits, read_len, orig_len, n_cumsum=None):
    """
    Detect and merge complementary split alignments that span the circular
    junction of a doubled reference.

    In a doubled reference of length 2*L (where L = orig_len), the junction
    sits at position L.  A junction-spanning read produces two hits:
      - one ending near L  (tail of copy 1)
      - one starting near L (head of copy 2)
    Together they cover the full read.

    n_cumsum: pre-computed cumulative N array for N-aware identity scoring.
              When provided, N positions are excluded from blen before computing
              identity so that masked regions do not inflate blen and reduce
              the computed identity.

    Returns a merged result dict (with is_junction_spanning=True) if a valid
    merge beats the best single hit; otherwise returns None.
    """
    pair = _find_junction_pair(hits, orig_len)
    if pair is None:
        return None

    ha, hb = pair
    L = orig_len
    merged_mlen = ha.mlen + hb.mlen
    merged_blen = ha.blen + hb.blen
    merged_q_st = min(ha.q_st, hb.q_st)
    merged_q_en = max(ha.q_en, hb.q_en)
    q_span = merged_q_en - merged_q_st
    q_cov = q_span / read_len

    # N-aware identity: subtract N positions from both mlen and blen.
    # minimap2 treats reference N's as wildcards (matching any query base),
    # so they inflate both mlen (match count) and blen (alignment length).
    # Subtracting n_in_hit from both removes the N-induced inflation.
    n_ha = _n_in_range(n_cumsum, ha.r_st, ha.r_en, L) if n_cumsum is not None else 0
    n_hb = _n_in_range(n_cumsum, hb.r_st, hb.r_en, L) if n_cumsum is not None else 0
    n_in_merged = n_ha + n_hb
    effective_mlen = max(0, merged_mlen - n_in_merged)
    effective_blen = max(1, merged_blen - n_in_merged)
    ident = effective_mlen / effective_blen

    match_frac = ident * q_cov   # identity × (q_span / read_len)
    mapq = max(ha.mapq, hb.mapq)

    # Wrap reference coordinates to [0, L)
    seg1 = {"r_st": ha.r_st % L, "r_en": min(ha.r_en, L) % L or L}
    seg2 = {"r_st": hb.r_st % L, "r_en": hb.r_en % L or L}

    return {
        "mlen":         merged_mlen,
        "blen":         merged_blen,
        "q_st":         merged_q_st,
        "q_en":         merged_q_en,
        "query_cov":    round(q_cov, 4),
        "identity":     round(ident, 4),
        "match_frac":   round(match_frac, 4),
        "mapq":         mapq,
        "strand":       "+" if ha.strand == 1 else "-",
        "n_in_hit":     n_in_merged,
        "is_junction_spanning": True,
        "junction_segments": [seg1, seg2],
    }


def best_hit_for_plasmid(hits, read_len, circular_len=None, orig_ref_len=None,
                         n_cumsum=None):
    """
    From all minimap2 hits to one plasmid, return the metrics of the
    single best alignment (highest mlen).

    If circular_len is provided (original plasmid length before doubling),
    attempts to merge split junction-spanning hits first.

    Score used for plasmid assignment:
        match_frac = identity * query_cov
                   = (effective_mlen / effective_blen) * (q_span / read_len)

    where effective_mlen/blen exclude N-position contributions (see n_cumsum).

    Using read_len as the sole denominator makes the score purely read-centric:
    backbone-only reads that align equally well to all references score the same
    for each and are correctly flagged ambiguous, while insert-specific reads
    score high only against their own reference.  The previous denominator
    max(read_len, orig_ref_len) created a systematic size-bias: a reference that
    is 30% smaller (e.g. lucxba 2550 bp vs gfp 3655 bp) would win for any
    backbone-only read regardless of actual alignment quality.

    n_cumsum: pre-computed cumulative N array.  When provided, N positions
    within the aligned range are excluded from both mlen and blen before
    computing identity, so that masked regions do not influence the score.
    """
    if not hits:
        return None

    # Try merging split junction hits in circular mode
    if circular_len is not None:
        merged = merge_circular_hits(hits, read_len, circular_len,
                                     n_cumsum=n_cumsum)
        if merged is not None:
            return merged

    best = None
    for h in hits:
        if best is None or h.mlen > best.mlen:
            best = h
    if best is None:
        return None

    q_span = best.q_en - best.q_st
    q_cov  = q_span / read_len

    # N-aware identity: subtract N positions from both mlen and blen.
    # minimap2 treats reference N's as wildcards (matching any query base),
    # so they inflate both mlen (match count) and blen (alignment length).
    # Subtracting n_in_hit from both removes the N-induced inflation and
    # prevents computed identity from exceeding 1.0.
    n_in_hit = 0
    if n_cumsum is not None and orig_ref_len is not None:
        n_in_hit = _n_in_range(n_cumsum, best.r_st, best.r_en, orig_ref_len)
    effective_mlen = max(0, best.mlen - n_in_hit)
    effective_blen = max(1, best.blen - n_in_hit)
    ident  = effective_mlen / effective_blen if effective_blen > 0 else 0.0

    match_frac = ident * q_cov   # identity × (q_span / read_len); primary ranking score
    mapq   = best.mapq

    # A single continuous alignment can still span the junction
    spans_junction = (circular_len is not None and
                      best.r_st < circular_len and best.r_en > circular_len)

    return {
        "mlen":         best.mlen,
        "blen":         best.blen,
        "q_st":         best.q_st,
        "q_en":         best.q_en,
        "query_cov":    round(q_cov,  4),
        "identity":     round(ident,  4),
        "match_frac":   round(match_frac, 4),
        "mapq":         mapq,
        "strand":       "+" if best.strand == 1 else "-",
        "n_in_hit":     n_in_hit,
        "is_junction_spanning": spans_junction,
    }


def align_read(read_seq, aligners, circular=False, orig_lens=None,
               n_cumsums=None):
    """
    Align one read against all plasmid aligners.
    Returns dict: {plasmid_name: hit_metrics_or_None}

    orig_lens is always used: as circular_len for junction-merge detection
    (when circular=True) and as orig_ref_len for length-normalised scoring
    in both modes.

    n_cumsums: per-reference cumulative N arrays from load_plasmid_library().
    When provided, N positions are excluded from blen before computing identity
    so that masked regions do not penalise references containing N runs.
    """
    read_len = len(read_seq)
    results  = {}
    for name, aln in aligners.items():
        hits = list(aln.map(read_seq))
        circ_len     = orig_lens.get(name) if (circular and orig_lens) else None
        orig_ref_len = orig_lens.get(name) if orig_lens else None
        n_cumsum     = n_cumsums.get(name) if n_cumsums else None
        results[name] = best_hit_for_plasmid(hits, read_len,
                                             circular_len=circ_len,
                                             orig_ref_len=orig_ref_len,
                                             n_cumsum=n_cumsum)
    return results


def classify_read(plasmid_hits, ambiguity_delta):
    """
    Given per-plasmid hit dicts, return:
        best_plasmid        – name of top-scoring plasmid (or 'unmapped')
        is_ambiguous        – bool
        second_plasmid      – name of runner-up (or None)
        reason              – short explanation string
    Ranking is by match_frac (matching bases / max(read_len, orig_ref_len)).
    """
    scored = {
        name: info["match_frac"]
        for name, info in plasmid_hits.items()
        if info is not None and info["match_frac"] > 0
    }

    if not scored:
        return "unmapped", False, None, "no alignments to any plasmid"

    ranked = sorted(scored.items(), key=lambda x: x[1], reverse=True)
    best_name, best_score = ranked[0]

    if len(ranked) == 1:
        return best_name, False, None, "unique alignment"

    second_name, second_score = ranked[1]

    # Ambiguity: gap between top two is smaller than delta * best_score
    gap = best_score - second_score
    if gap < ambiguity_delta * best_score:
        reason = (f"ambiguous: {best_name} ({best_score:.3f}) vs "
                  f"{second_name} ({second_score:.3f}), "
                  f"delta={gap:.3f}")
        return best_name, True, second_name, reason

    return best_name, False, second_name, "clear best alignment"


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

def _short(name, maxlen=22):
    """Truncate long plasmid names for axis labels."""
    return name if len(name) <= maxlen else name[:maxlen - 1] + "…"


def generate_readqc_plot(pre_filter_df, plasmid_names, min_length, min_qual, pdf):
    """
    Scatter of ALL reads (pre-filter): read length (x) vs average quality (y),
    with stacked marginal histograms on the top and right axes.
    Points coloured by filter status or plasmid assignment.
    """
    n_plasmids = len(plasmid_names)
    cmap = plt.get_cmap("tab20" if n_plasmids > 10 else "tab10")
    pal  = {p: cmap(i % 20) for i, p in enumerate(plasmid_names)}
    pal["unmapped"]  = (0.55, 0.55, 0.55, 1.0)
    pal["ambiguous"] = (0.95, 0.70, 0.10, 1.0)
    slabels = {p: _short(p) for p in plasmid_names}

    # Category order — filtered reads drawn first (background), plasmids last (foreground)
    cat_order = [
        ("filtered_short",   "#cccccc", "Filtered — too short"),
        ("filtered_lowqual", "#d4a056", "Filtered — low quality"),
        ("unmapped",         pal["unmapped"], "Unmapped"),
        ("ambiguous",        pal["ambiguous"], "Ambiguous"),
    ] + [(p, pal[p], slabels[p]) for p in plasmid_names]

    plt.rcParams.update({
        "font.size": 9, "axes.spines.top": False,
        "axes.spines.right": False, "figure.dpi": 150,
    })

    # Work in log2 space for the x-axis so rare long reads don't dominate
    log2_len  = np.log2(pre_filter_df["read_length"].clip(lower=1))
    log2_min  = np.log2(max(pre_filter_df["read_length"].min(), 1))
    log2_max  = np.log2(pre_filter_df["read_length"].max())
    log2_minl = np.log2(min_length)

    fig = plt.figure(figsize=(12, 9))
    gs  = fig.add_gridspec(2, 2,
                            width_ratios=[4, 1], height_ratios=[1, 5],
                            hspace=0.05, wspace=0.05)
    ax_top   = fig.add_subplot(gs[0, 0])
    ax_main  = fig.add_subplot(gs[1, 0], sharex=ax_top)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
    fig.add_subplot(gs[0, 1]).set_visible(False)   # blank corner
    fig.suptitle("Read Quality Overview — All Reads (pre-filter)",
                 fontsize=13, fontweight="bold")

    # ── main scatter (x = log2 read length) ──────────────────────────────────
    for cat_key, col, label in cat_order:
        mask = pre_filter_df["_scatter_cat"] == cat_key
        if not mask.any():
            continue
        ax_main.scatter(log2_len[mask], pre_filter_df.loc[mask, "avg_qual"],
                        c=[col], s=6, alpha=0.5,
                        label=f"{label} (n={mask.sum():,})",
                        rasterized=True, linewidths=0)
    ax_main.axvline(log2_minl, ls="--", color="#444444", lw=1.0, alpha=0.7,
                    label=f"Min length ({min_length:,} bp)")
    ax_main.axhline(min_qual,  ls="--", color="#444444", lw=1.0, alpha=0.7,
                    label=f"Min quality (Q{min_qual})")
    ax_main.set_xlabel("Read length (bp, log₂ scale)")
    ax_main.set_ylabel("Average Phred quality (Q)")
    ax_main.legend(fontsize=7, markerscale=2, loc="upper right", framealpha=0.85)
    # Tick at every power-of-2 that falls in range; label in bp
    tick_powers = np.arange(np.floor(log2_min), np.ceil(log2_max) + 1)
    ax_main.set_xticks(tick_powers)
    ax_main.set_xticklabels(
        [f"{int(2**p/1000)}k" if 2**p >= 1000 else str(int(2**p))
         for p in tick_powers], fontsize=8)
    ax_main.set_xlim(log2_min - 0.1, log2_max + 0.1)

    # ── top marginal: read length histogram (log2-spaced bins) ───────────────
    bins_x     = np.linspace(log2_min, log2_max, 60)
    nonempty_x = [(log2_len[pre_filter_df["_scatter_cat"] == k].values, col)
                  for k, col, _ in cat_order
                  if (pre_filter_df["_scatter_cat"] == k).sum() > 0]
    ax_top.hist([d for d, _ in nonempty_x], bins=bins_x,
                color=[c for _, c in nonempty_x],
                stacked=True, alpha=0.85, edgecolor="none")
    ax_top.axvline(log2_minl, ls="--", color="#444444", lw=1.0, alpha=0.7)
    ax_top.set_ylabel("Count")
    ax_top.spines["bottom"].set_visible(False)
    plt.setp(ax_top.get_xticklabels(), visible=False)

    # ── right marginal: quality histogram (horizontal) ────────────────────────
    bins_y     = np.linspace(pre_filter_df["avg_qual"].min(),
                              pre_filter_df["avg_qual"].max(), 40)
    nonempty_y = [(pre_filter_df[pre_filter_df["_scatter_cat"] == k]["avg_qual"].values, col)
                  for k, col, _ in cat_order
                  if (pre_filter_df["_scatter_cat"] == k).sum() > 0]
    ax_right.hist([d for d, _ in nonempty_y], bins=bins_y,
                  color=[c for _, c in nonempty_y],
                  stacked=True, alpha=0.85, edgecolor="none",
                  orientation="horizontal")
    ax_right.axhline(min_qual, ls="--", color="#444444", lw=1.0, alpha=0.7)
    ax_right.set_xlabel("Count")
    ax_right.spines["left"].set_visible(False)
    plt.setp(ax_right.get_yticklabels(), visible=False)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def generate_plots(results_df, plasmid_names, counts, pdf, orig_lens=None):
    """
    Write a multi-page PDF report to plots_path.

    counts dict keys: n_total, n_short, n_lowqual, n_passed,
                      n_unmapped, n_ambiguous
    """
    n_plasmids = len(plasmid_names)
    cmap       = plt.get_cmap("tab20" if n_plasmids > 10 else "tab10")
    pal        = {p: cmap(i % 20) for i, p in enumerate(plasmid_names)}
    pal["unmapped"]  = (0.55, 0.55, 0.55, 1.0)
    pal["ambiguous"] = (0.95, 0.70, 0.10, 1.0)

    # Convenience: short label map
    slabels = {p: _short(p) for p in plasmid_names}

    # Pre-build match_frac column names
    mf_cols = {p: p.replace(" ", "_").replace("/", "_").replace(":", "_")
                  + "__match_frac"
               for p in plasmid_names}
    id_cols = {p: p.replace(" ", "_").replace("/", "_").replace(":", "_")
                  + "__identity"
               for p in plasmid_names}
    qc_cols = {p: p.replace(" ", "_").replace("/", "_").replace(":", "_")
                  + "__query_cov"
               for p in plasmid_names}

    plt.rcParams.update({
        "font.size": 9, "axes.spines.top": False,
        "axes.spines.right": False, "figure.dpi": 150,
    })

    if True:  # keep indent level consistent; pdf is passed in

        # ── PAGE 1: read disposition overview ────────────────────────────────
        fig = plt.figure(figsize=(13, 9))
        gs1 = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.35,
                               height_ratios=[1, 1])
        ax00 = fig.add_subplot(gs1[0, 0])
        ax01 = fig.add_subplot(gs1[0, 1])
        ax10 = fig.add_subplot(gs1[1, :])   # read length spans full bottom row
        fig.suptitle("Plasmid Alignment — Read Overview", fontsize=13,
                     fontweight="bold")

        # [0,0] Stacked bar: read fate at every stage
        ax = ax00
        n_unique = counts["n_passed"] - counts["n_unmapped"] - counts["n_ambiguous"]
        stages   = ["Raw reads", "Passed filters", "Uniquely assigned"]
        totals   = [counts["n_total"], counts["n_passed"], n_unique]
        colors_s = ["#4393c3", "#74add1", "#2166ac"]
        bars = ax.barh(stages, totals, color=colors_s, edgecolor="white")
        for bar, val in zip(bars, totals):
            ax.text(val + max(totals) * 0.01, bar.get_y() + bar.get_height() / 2,
                    f"{val:,}", va="center", fontsize=8)
        # Overlay filtered portions
        ax.barh(["Passed filters"],
                [counts["n_short"] + counts["n_lowqual"]],
                left=[n_unique + counts["n_unmapped"] + counts["n_ambiguous"]],
                color="#d9d9d9", edgecolor="white", label="Filtered")
        ax.set_xlabel("Read count")
        ax.set_title("Read Disposition")
        ax.set_xlim(0, max(totals) * 1.15)

        # annotation boxes
        ax2_text = (
            f"Filtered short: {counts['n_short']:,}\n"
            f"Filtered low-Q: {counts['n_lowqual']:,}\n"
            f"Unmapped:       {counts['n_unmapped']:,}\n"
            f"Ambiguous:      {counts['n_ambiguous']:,}"
        )
        ax.text(0.98, 0.05, ax2_text, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#aaaaaa"))

        # [0,1] Horizontal bar: unique reads per plasmid + single Ambiguous bar
        ax = ax01
        unique_counts = (results_df[~results_df["is_ambiguous"]]
                         ["best_plasmid"].value_counts()
                         .reindex(plasmid_names, fill_value=0))
        ambig_counts  = (results_df[results_df["is_ambiguous"]]
                         ["best_plasmid"].value_counts()
                         .reindex(plasmid_names, fill_value=0))
        total_ambiguous = counts["n_ambiguous"]
        sort_idx  = unique_counts.values.argsort()
        y_pos     = np.arange(n_plasmids)
        bar_cols  = [pal[p] for p in unique_counts.index[sort_idx]]
        # Unique reads per plasmid (solid plasmid-coloured bars, no stacking)
        ax.barh(y_pos, unique_counts.values[sort_idx],
                color=bar_cols, edgecolor="white", label="Unique")
        # Single "Ambiguous" bar at top (all ambiguous reads combined)
        ambig_y = n_plasmids
        ax.barh(ambig_y, total_ambiguous,
                color=pal["ambiguous"], edgecolor="white",
                alpha=0.85, label="Ambiguous")
        all_y   = list(y_pos) + [ambig_y]
        all_lbl = [slabels[p] for p in unique_counts.index[sort_idx]] + ["Ambiguous"]
        ax.set_yticks(all_y)
        ax.set_yticklabels(all_lbl, fontsize=8)
        # Annotate plasmid bars with unique + (unique+ambiguous-best-hit) counts
        n_unique_total   = counts["n_passed"] - counts["n_unmapped"] - counts["n_ambiguous"]
        n_assigned_total = counts["n_passed"] - counts["n_unmapped"]
        bar_xlim_max = max(
            (unique_counts.values + ambig_counts.values).max() if len(unique_counts) else 0,
            total_ambiguous,
        ) * 1.45
        for yi, p in enumerate(unique_counts.index[sort_idx]):
            u   = unique_counts[p]
            am  = ambig_counts[p]
            tot = u + am
            pct_u   = 100 * u   / n_unique_total   if n_unique_total   else 0
            pct_tot = 100 * tot / n_assigned_total if n_assigned_total else 0
            ax.text(u + bar_xlim_max * 0.01, yi,
                    f"u={u:,} ({pct_u:.1f}% uniq)  u+a={tot:,} ({pct_tot:.1f}% assigned)",
                    va="center", fontsize=7, color="#333333")
        # Annotate ambiguous bar
        if n_assigned_total:
            pct_a = 100 * total_ambiguous / n_assigned_total
            ax.text(total_ambiguous + bar_xlim_max * 0.01, ambig_y,
                    f"{total_ambiguous:,} ({pct_a:.1f}% assigned)",
                    va="center", fontsize=7, color="#333333")
        ax.set_xlim(0, bar_xlim_max)
        ax.set_xlabel("Read count")
        ax.set_title("Reads per Plasmid")
        ax.legend(fontsize=8, loc="lower right")

        # [1, :] Read length distribution — spans full bottom row
        ax = ax10
        # Build category column: plasmid name, "ambiguous", or "unmapped"
        def read_cat(row):
            if row["best_plasmid"] == "unmapped":
                return "unmapped"
            if row["is_ambiguous"]:
                return "ambiguous"
            return row["best_plasmid"]

        results_df["_cat"] = results_df.apply(read_cat, axis=1)
        max_ref_len  = max(orig_lens.values()) if orig_lens else results_df["read_length"].max()
        xmax         = 2 * max_ref_len
        n_clipped    = (results_df["read_length"] > xmax).sum()
        bins = np.linspace(results_df["read_length"].min(), xmax, 50)

        # Order: plasmids (unique reads only), ambiguous, unmapped
        cat_order  = list(plasmid_names) + ["ambiguous", "unmapped"]
        cat_colors = [pal[c] for c in cat_order]
        cat_data   = [results_df.loc[results_df["_cat"] == c, "read_length"].values
                      for c in cat_order]
        cat_labels = [f"{slabels.get(c, c)} (n={len(d):,})"
                      for c, d in zip(cat_order, cat_data)]

        # Drop empty categories so matplotlib stacks correctly
        nonempty = [(d, col, lbl)
                    for d, col, lbl in zip(cat_data, cat_colors, cat_labels)
                    if len(d) > 0]
        ax.hist([t[0] for t in nonempty], bins=bins,
                color=[t[1] for t in nonempty],
                label=[t[2] for t in nonempty],
                stacked=True, alpha=0.85, edgecolor="none")
        ax.set_xlim(bins[0], xmax)
        ax.set_xlabel("Read length (bp)")
        ax.set_ylabel("Read count")
        title = "Read Length Distribution by Assignment"
        if n_clipped:
            title += f"\n({n_clipped:,} read{'s' if n_clipped>1 else ''} >{int(xmax/1000)}k bp not shown)"
        ax.set_title(title, fontsize=9)
        ax.legend(fontsize=7, loc="upper right", ncol=2)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(
            lambda x, _: f"{int(x/1000)}k" if x >= 1000 else str(int(x))))

        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ── PAGE 2: per-plasmid alignment quality ────────────────────────────
        fig, axes = plt.subplots(2, 2, figsize=(13, 9))
        fig.suptitle("Plasmid Alignment — Alignment Quality", fontsize=13,
                     fontweight="bold")

        # [0,0] Match fraction violin — distribution across uniquely assigned reads
        ax = axes[0][0]
        vdata, vlabels, vcolors = [], [], []
        for p in plasmid_names:
            col = mf_cols[p]
            if col not in results_df.columns:
                continue
            sub = results_df.loc[results_df["_cat"] == p, col].dropna()
            if len(sub) < 2:
                continue
            vdata.append(sub.values)
            vlabels.append(slabels[p])
            vcolors.append(pal[p])

        if vdata:
            parts = ax.violinplot(vdata, positions=range(len(vdata)),
                                  showmedians=True, showextrema=False,
                                  widths=0.7)
            for body, col in zip(parts["bodies"], vcolors):
                body.set_facecolor(col); body.set_alpha(0.7)
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(1.5)
            ax.set_xticks(range(len(vlabels)))
            ax.set_xticklabels(vlabels, rotation=35, ha="right", fontsize=8)
        ax.set_ylabel("Match fraction (matching bases / read length)")
        ax.set_ylim(0, 1)
        ax.set_title("Match Fraction — Uniquely Assigned Reads")
        ax.axhline(0.5, ls="--", color="grey", lw=0.8, alpha=0.5)

        # [0,1] Top-2 score scatter: best vs second-best match_frac per read
        ax = axes[0][1]
        # For each read, find best and second-best match_frac across plasmids
        mf_matrix = results_df[[mf_cols[p] for p in plasmid_names
                                  if mf_cols[p] in results_df.columns]].values
        sorted_mf = np.sort(mf_matrix, axis=1)[:, ::-1]  # descending per row
        best_mf   = sorted_mf[:, 0]
        second_mf = sorted_mf[:, 1] if sorted_mf.shape[1] > 1 else np.zeros(len(sorted_mf))

        has_aln   = best_mf > 0
        ambig_mask = results_df["is_ambiguous"].values
        unmap_mask = (results_df["best_plasmid"] == "unmapped").values

        for mask, label, col, alpha, size in [
            (has_aln & ~ambig_mask & ~unmap_mask, "Unique",    "#2166ac", 0.5, 8),
            (has_aln & ambig_mask,                "Ambiguous", "#f4a582", 0.8, 12),
            (~has_aln | unmap_mask,               "Unmapped",  "#aaaaaa", 0.4, 6),
        ]:
            if mask.sum() > 0:
                ax.scatter(best_mf[mask], second_mf[mask],
                           c=col, alpha=alpha, s=size, label=f"{label} (n={mask.sum():,})")

        lim = max(best_mf.max(), second_mf.max()) * 1.05
        ax.plot([0, lim], [0, lim], ls="--", color="grey", lw=1,
                label="Equal score (diagonal)")
        ax.set_xlabel("Best plasmid match fraction")
        ax.set_ylabel("Second-best plasmid match fraction")
        ax.set_title("Top-2 Plasmid Scores per Read\n"
                     "(reads near diagonal → ambiguous)")
        ax.legend(fontsize=8)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)

        # [1,0] Identity vs query coverage scatter (uniquely assigned reads)
        ax = axes[1][0]
        for p in plasmid_names:
            sub = results_df[results_df["_cat"] == p]
            if sub.empty:
                continue
            ic = id_cols[p]; qc = qc_cols[p]
            if ic not in results_df.columns or qc not in results_df.columns:
                continue
            ax.scatter(sub[qc], sub[ic], c=[pal[p]], s=8, alpha=0.5,
                       label=slabels[p])
        ax.set_xlabel("Query coverage (fraction of read aligned)")
        ax.set_ylabel("Alignment identity (matching / aligned bases)")
        ax.set_title("Identity vs Query Coverage\n(uniquely assigned reads)")
        ax.legend(fontsize=7, ncol=2, markerscale=2)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)

        # [1,1] Read length vs match fraction (uniquely assigned reads)
        ax = axes[1][1]
        for p in plasmid_names:
            sub = results_df[results_df["_cat"] == p]
            if sub.empty:
                continue
            mf = mf_cols[p]
            if mf not in results_df.columns:
                continue
            ax.scatter(sub["read_length"], sub[mf], c=[pal[p]], s=8,
                       alpha=0.5, label=slabels[p])
        ax.set_xlabel("Read length (bp)")
        ax.set_ylabel("Match fraction")
        ax.set_title("Read Length vs Match Fraction\n(uniquely assigned reads)")
        ax.legend(fontsize=7, ncol=2, markerscale=2)
        ax.set_ylim(0, 1)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(
            lambda x, _: f"{int(x/1000)}k" if x >= 1000 else str(int(x))))

        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # ── PAGE 3: cross-alignment heatmap ──────────────────────────────────
        # For reads uniquely assigned to plasmid R, what is the median
        # match_frac to every other plasmid C?
        # Reveals shared sequence (backbone etc.) driving ambiguity.
        valid_plasmids = [p for p in plasmid_names
                          if mf_cols[p] in results_df.columns
                          and (results_df["_cat"] == p).sum() > 0]

        if len(valid_plasmids) >= 2:
            heatmap = np.zeros((len(valid_plasmids), len(valid_plasmids)))
            for r_idx, row_p in enumerate(valid_plasmids):
                sub = results_df[results_df["_cat"] == row_p]
                for c_idx, col_p in enumerate(valid_plasmids):
                    mc = mf_cols[col_p]
                    if mc in sub.columns:
                        heatmap[r_idx, c_idx] = sub[mc].median()

            fig, ax = plt.subplots(figsize=(max(7, len(valid_plasmids) * 0.9 + 2),
                                            max(5, len(valid_plasmids) * 0.8 + 2)))
            fig.suptitle("Cross-Alignment Heatmap\n"
                         "Median match fraction of reads assigned to row-plasmid "
                         "when aligned to column-plasmid\n"
                         "(off-diagonal values reveal shared sequence / backbone)",
                         fontsize=11, fontweight="bold")

            im = ax.imshow(heatmap, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
            plt.colorbar(im, ax=ax, label="Median match fraction", shrink=0.8)

            tick_labels = [_short(p, 28) for p in valid_plasmids]
            ax.set_xticks(range(len(valid_plasmids)))
            ax.set_yticks(range(len(valid_plasmids)))
            ax.set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=8)
            ax.set_yticklabels(tick_labels, fontsize=8)
            ax.set_xlabel("Aligned to (column plasmid)")
            ax.set_ylabel("Assigned reads from (row plasmid)")

            for r in range(len(valid_plasmids)):
                for c in range(len(valid_plasmids)):
                    val = heatmap[r, c]
                    ax.text(c, r, f"{val:.2f}", ha="center", va="center",
                            fontsize=7,
                            color="white" if val > 0.6 else "black")

            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # ── PAGE 4: per-plasmid match_frac distributions (all reads) ─────────
        # Shows how well each plasmid discriminates its own reads from others.
        fig, axes_flat = plt.subplots(
            max(1, (n_plasmids + 1) // 2), 2,
            figsize=(13, max(4, n_plasmids * 1.4)),
            squeeze=False,
        )
        fig.suptitle("Per-Plasmid Match Fraction — All Reads\n"
                     "(blue = reads assigned to this plasmid, "
                     "grey = all other reads)",
                     fontsize=12, fontweight="bold")

        for idx, p in enumerate(plasmid_names):
            ax   = axes_flat[idx // 2][idx % 2]
            mcol = mf_cols[p]
            if mcol not in results_df.columns:
                ax.set_visible(False)
                continue
            assigned_vals = results_df.loc[results_df["_cat"] == p, mcol]
            other_vals    = results_df.loc[results_df["_cat"] != p, mcol]
            bins = np.linspace(0, 1, 41)
            ax.hist(other_vals, bins=bins, color="#aaaaaa", alpha=0.55,
                    density=True, label=f"Other reads (n={len(other_vals):,})")
            if not assigned_vals.empty:
                ax.hist(assigned_vals, bins=bins, color=pal[p], alpha=0.75,
                        density=True,
                        label=f"Assigned reads (n={len(assigned_vals):,})")
            ax.set_title(slabels[p], fontsize=9)
            ax.set_xlabel("Match fraction")
            ax.set_ylabel("Density")
            ax.legend(fontsize=7)

        # Hide unused subplots
        total_axes = axes_flat.shape[0] * 2
        for idx in range(n_plasmids, total_axes):
            axes_flat[idx // 2][idx % 2].set_visible(False)

        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        pass  # metadata written by caller


# ─────────────────────────────────────────────────────────────────────────────
# Pileup plots
# ─────────────────────────────────────────────────────────────────────────────

# Regex to tokenise a minimap2 cs string
_CS_RE = re.compile(r':[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+')

N_MAX_PILEUP = 500   # reads shown per plasmid pileup panel


def parse_cs(cs_string, r_st):
    """
    Convert a minimap2 cs string into a list of typed segments in reference
    coordinates.  Each segment is a dict with key 'type':
        'match'     — r_st, r_en
        'mismatch'  — r_pos  (single reference base)
        'insertion' — r_pos, size  (does not consume reference bases)
        'deletion'  — r_st, r_en
    """
    segments = []
    ref_pos  = r_st
    for token in _CS_RE.findall(cs_string):
        ch = token[0]
        if ch == ':':                          # matches
            n = int(token[1:])
            segments.append(dict(type='match', r_st=ref_pos, r_en=ref_pos + n))
            ref_pos += n
        elif ch == '*':                        # single mismatch
            segments.append(dict(type='mismatch', r_pos=ref_pos))
            ref_pos += 1
        elif ch == '+':                        # insertion (no ref advance)
            segments.append(dict(type='insertion', r_pos=ref_pos,
                                 size=len(token) - 1))
        elif ch == '-':                        # deletion
            n = len(token) - 1
            segments.append(dict(type='deletion', r_st=ref_pos,
                                 r_en=ref_pos + n))
            ref_pos += n
    return segments


def pack_into_rows(intervals, gap=8):
    """
    Greedy interval packing — assign each read to the first row where it
    fits without overlapping.  Returns (row_assignments list, n_rows int).
    """
    order      = sorted(range(len(intervals)), key=lambda i: intervals[i][0])
    row_ends   = []
    assignment = [0] * len(intervals)
    for i in order:
        start, end = intervals[i]
        placed = False
        for row_idx, row_end in enumerate(row_ends):
            if start >= row_end + gap:
                assignment[i] = row_idx
                row_ends[row_idx] = end
                placed = True
                break
        if not placed:
            assignment[i] = len(row_ends)
            row_ends.append(end)
    return assignment, len(row_ends)


def _n_runs(seq):
    """Return list of (start, end) for runs of N/n in seq."""
    runs, in_n, n_st = [], False, 0
    for i, b in enumerate(seq.upper()):
        if b == 'N' and not in_n:
            in_n, n_st = True, i
        elif b != 'N' and in_n:
            runs.append((n_st, i)); in_n = False
    if in_n:
        runs.append((n_st, len(seq)))
    return runs


def wrap_hit_to_segments(cs_segments, r_st, r_en, orig_len):
    """
    If an alignment on a doubled reference spans the junction at position L
    (orig_len), split the parsed cs segments into two display groups with
    coordinates wrapped to [0, L).

    Returns a list of dicts:
        [{'r_st': int, 'r_en': int, 'segments': list, 'is_wrapped': bool}, ...]

    Entirely within copy 1 → one group, not wrapped.
    Entirely within copy 2 → one group, coords shifted by −L.
    Spans junction → two groups: [r_st, L) and [0, r_en − L).
    """
    L = orig_len

    if r_en <= L:
        return [{'r_st': r_st, 'r_en': r_en, 'segments': cs_segments,
                 'is_wrapped': False}]

    if r_st >= L:
        shifted = []
        for seg in cs_segments:
            s = dict(seg)
            for key in ('r_st', 'r_en', 'r_pos'):
                if key in s:
                    s[key] -= L
            shifted.append(s)
        return [{'r_st': r_st - L, 'r_en': r_en - L, 'segments': shifted,
                 'is_wrapped': False}]

    # Spans junction — split at L
    before, after = [], []
    for seg in cs_segments:
        t = seg['type']
        if t == 'match':
            if seg['r_en'] <= L:
                before.append(seg)
            elif seg['r_st'] >= L:
                after.append({'type': 'match',
                              'r_st': seg['r_st'] - L, 'r_en': seg['r_en'] - L})
            else:
                before.append({'type': 'match', 'r_st': seg['r_st'], 'r_en': L})
                after.append({'type': 'match', 'r_st': 0,
                              'r_en': seg['r_en'] - L})
        elif t == 'mismatch':
            if seg['r_pos'] < L:
                before.append(seg)
            else:
                after.append({'type': 'mismatch', 'r_pos': seg['r_pos'] - L})
        elif t == 'insertion':
            if seg['r_pos'] < L:
                before.append(seg)
            else:
                after.append({'type': 'insertion', 'r_pos': seg['r_pos'] - L,
                              'size': seg['size']})
        elif t == 'deletion':
            if seg['r_en'] <= L:
                before.append(seg)
            elif seg['r_st'] >= L:
                after.append({'type': 'deletion',
                              'r_st': seg['r_st'] - L, 'r_en': seg['r_en'] - L})
            else:
                before.append({'type': 'deletion', 'r_st': seg['r_st'],
                               'r_en': L})
                after.append({'type': 'deletion', 'r_st': 0,
                              'r_en': seg['r_en'] - L})

    result = []
    if before:
        result.append({'r_st': r_st, 'r_en': L, 'segments': before,
                       'is_wrapped': False})
    if after:
        result.append({'r_st': 0, 'r_en': r_en - L, 'segments': after,
                       'is_wrapped': True})
    return result


def generate_pileup_plots(results_df, plasmid_seqs_orig, aligners, read_seqs,
                           pdf, circular=False, orig_lens=None):
    """
    One PDF page per plasmid with assigned reads.

    Layout per page
    ───────────────
    Top strip  : coverage depth histogram with N-region shading
    Main panel : per-read alignment tracks, packed into rows
                 Green  = matching bases
                 Red    = mismatches (vertical tick)
                 I{n}   = insertion (text above read, dark-red)
                 ██████ dark grey = deletion
                 Read border: blue (unique) | orange (ambiguous) | purple (junction)

    In circular mode aligners are built on DOUBLED sequences so that
    junction-spanning reads align fully.  Display coordinates are wrapped
    back to [0, orig_len) and a dashed junction line is drawn at x = 0.
    """
    assigned = results_df[results_df["best_plasmid"] != "unmapped"]

    for plasmid_name in plasmid_seqs_orig:
        if plasmid_name not in aligners:
            continue
        aln = aligners[plasmid_name]
        ref_seq = plasmid_seqs_orig[plasmid_name]
        ref_len = len(ref_seq)             # original (non-doubled) length
        pgroup  = assigned[assigned["best_plasmid"] == plasmid_name].copy()
        if pgroup.empty:
            continue

        total_reads = len(pgroup)
        note        = ""
        if total_reads > N_MAX_PILEUP:
            pgroup = pgroup.sample(N_MAX_PILEUP, random_state=42)
            note   = f", showing {N_MAX_PILEUP} of {total_reads}"

        print(f"  Pileup: {plasmid_name}  ({len(pgroup)} reads{note}) …",
              flush=True)

        # ── re-align with cs=True ────────────────────────────────────────────
        read_data = []
        for _, row in pgroup.iterrows():
            seq = read_seqs.get(row["read_id"])
            if seq is None:
                continue
            hits = list(aln.map(seq, cs=True))
            if not hits:
                continue

            is_ambig = bool(row["is_ambiguous"])
            is_junct = False

            rlen = int(row["read_length"])
            if circular:
                L = ref_len
                pair = _find_junction_pair(hits, ref_len)
                if pair:
                    ha, hb = pair
                    display_groups = []
                    for h in (ha, hb):
                        if not hasattr(h, 'cs') or h.cs is None:
                            continue
                        raw_segs = parse_cs(h.cs, h.r_st)
                        groups = wrap_hit_to_segments(raw_segs, h.r_st,
                                                     h.r_en, L)
                        display_groups.extend(groups)

                    if display_groups:
                        max_q_en = max(ha.q_en, hb.q_en)
                        read_data.append(dict(
                            display_groups = display_groups,
                            is_ambiguous   = is_ambig,
                            is_junction    = True,
                            read_length    = rlen,
                            max_q_en       = max_q_en,
                        ))
                    continue

                # Use best single hit, wrap if needed
                h = max(hits, key=lambda x: x.mlen)
                if not hasattr(h, 'cs') or h.cs is None:
                    continue
                raw_segs = parse_cs(h.cs, h.r_st)
                display_groups = wrap_hit_to_segments(raw_segs, h.r_st,
                                                     h.r_en, L)
                # A single continuous hit can still span the junction
                single_spans_junc = (h.r_st < L and h.r_en > L)
                read_data.append(dict(
                    display_groups = display_groups,
                    is_ambiguous   = is_ambig,
                    is_junction    = single_spans_junc,
                    read_length    = rlen,
                    max_q_en       = h.q_en,
                ))
            else:
                # Non-circular: original single-hit path
                h = max(hits, key=lambda x: x.mlen)
                if not hasattr(h, 'cs') or h.cs is None:
                    continue
                read_data.append(dict(
                    display_groups = [{'r_st': h.r_st, 'r_en': h.r_en,
                                       'segments': parse_cs(h.cs, h.r_st),
                                       'is_wrapped': False}],
                    is_ambiguous   = is_ambig,
                    is_junction    = False,
                    read_length    = rlen,
                    max_q_en       = h.q_en,
                ))

        if not read_data:
            continue

        # ── packing ──────────────────────────────────────────────────────────
        # For each read compute its display interval(s).  Junction-spanning
        # reads that wrap around are packed as full-width.
        intervals = []
        for rd in read_data:
            if rd["is_junction"]:
                intervals.append((0, ref_len))    # full-width packing
            else:
                all_st = min(g['r_st'] for g in rd['display_groups'])
                all_en = max(g['r_en'] for g in rd['display_groups'])
                intervals.append((all_st, all_en))

        # Sort by start position for stable packing
        sort_order = sorted(range(len(intervals)),
                            key=lambda i: intervals[i][0])
        read_data  = [read_data[i] for i in sort_order]
        intervals  = [intervals[i] for i in sort_order]

        row_assign, n_rows = pack_into_rows(intervals)
        n_regions          = _n_runs(ref_seq)

        # ── figure sizing ────────────────────────────────────────────────────
        ROW_H_IN  = 0.16            # inches per read row
        pile_h_in = max(2.0, n_rows * ROW_H_IN)
        fig_h     = 2.5 + pile_h_in
        fig, (ax_cov, ax_pile) = plt.subplots(
            2, 1, figsize=(14, fig_h),
            gridspec_kw={"height_ratios": [2.0, pile_h_in]},
        )
        circ_tag = " (circular)" if circular else ""
        fig.suptitle(
            f"Pileup — {plasmid_name}{circ_tag}   "
            f"({len(read_data)} reads{note}, {ref_len:,} bp)",
            fontsize=11, fontweight="bold",
        )

        # ── coverage strip ───────────────────────────────────────────────────
        cov = np.zeros(ref_len, dtype=np.int32)
        for rd in read_data:
            for grp in rd['display_groups']:
                st = max(0, grp['r_st'])
                en = min(ref_len, grp['r_en'])
                if st < en:
                    cov[st:en] += 1
        ax_cov.fill_between(np.arange(ref_len), cov,
                            color="#2166ac", alpha=0.75, step="mid")
        for nst, nen in n_regions:
            ax_cov.axvspan(nst, nen, color="#d0d0d0", alpha=0.6,
                           label="N region" if nst == n_regions[0][0] else "")
        if n_regions:
            ax_cov.legend(fontsize=7, loc="upper right")
        ax_cov.set_xlim(0, ref_len)
        ax_cov.set_ylabel("Depth")
        ax_cov.set_title("Coverage depth")
        ax_cov.xaxis.set_visible(False)
        ax_cov.spines["top"].set_visible(False)
        ax_cov.spines["right"].set_visible(False)

        # Junction line in circular mode
        if circular:
            ax_cov.axvline(0, ls="--", color="#7b3294", lw=1.0, alpha=0.7)

        # ── pileup panel ─────────────────────────────────────────────────────
        H      = 0.70   # bar height in data (row) units
        INS_FS = max(2.5, min(4.5, 180 / max(n_rows, 1)))  # insertion font size

        ax_pile.set_xlim(0, ref_len)
        ax_pile.set_ylim(-0.5, n_rows + 0.5)
        ax_pile.set_yticks([])
        ax_pile.set_xlabel("Reference position (bp)")
        ax_pile.set_ylabel("Reads (packed rows)")
        ax_pile.spines["top"].set_visible(False)
        ax_pile.spines["right"].set_visible(False)

        for nst, nen in n_regions:
            ax_pile.axvspan(nst, nen, color="#ebebeb", alpha=0.7, zorder=0)

        # Junction line in circular mode
        if circular:
            ax_pile.axvline(0, ls="--", color="#7b3294", lw=1.0, alpha=0.7,
                            label="Junction (pos 0)")

        # Batch draw helpers — collect rectangles for efficiency
        match_rects, del_rects = [], []
        mm_xs, mm_ys           = [], []
        ins_labels             = []

        junc_labels = []   # (x, yc, text) for junction soft-clip annotations

        for rd, ri in zip(read_data, row_assign):
            yc = ri + 0.5
            if rd["is_ambiguous"]:
                bord = "#f4a582"
                alph = 0.55
            else:
                bord = "#2166ac"
                alph = 0.85

            if rd["is_junction"]:
                sc = rd["read_length"] - rd["max_q_en"]
                sc_label = "0" if sc <= 0 else f"{sc}bp"
                junc_labels.append((ref_len, yc, sc_label))

            # Draw background bar(s) for each display group
            for grp in rd['display_groups']:
                w = grp['r_en'] - grp['r_st']
                if w <= 0:
                    continue
                ax_pile.broken_barh(
                    [(grp['r_st'], w)],
                    (yc - H/2, H),
                    facecolors="#e4e4e4", edgecolors=bord,
                    linewidth=0.4,
                    zorder=1,
                )

                for seg in grp['segments']:
                    t = seg["type"]
                    if t == "match":
                        sw = seg["r_en"] - seg["r_st"]
                        if sw > 0:
                            match_rects.append(
                                (seg["r_st"], yc - H/2, sw, H, alph))
                    elif t == "mismatch":
                        mm_xs.append(seg["r_pos"] + 0.5)
                        mm_ys.append((yc - H/2, yc + H/2))
                    elif t == "deletion":
                        sw = seg["r_en"] - seg["r_st"]
                        if sw > 0:
                            del_rects.append(
                                (seg["r_st"], yc - H/4, sw, H/2))
                    elif t == "insertion":
                        ins_labels.append(
                            (seg["r_pos"], yc + H/2, seg["size"]))

        # Draw collected primitives
        for (x, yb, w, h_bar, alph) in match_rects:
            ax_pile.broken_barh([(x, w)], (yb, h_bar),
                                facecolors="#4dac26", edgecolors="none",
                                alpha=alph, zorder=2)

        for (x, yb, w, h_bar) in del_rects:
            ax_pile.broken_barh([(x, w)], (yb, h_bar),
                                facecolors="#404040", edgecolors="none",
                                zorder=3)

        for x, (y0, y1) in zip(mm_xs, mm_ys):
            ax_pile.plot([x, x], [y0, y1],
                         color="#d7191c", lw=0.8, zorder=4)

        for (x, y_top, size) in ins_labels:
            ax_pile.text(x, y_top + 0.02, f"I{size}",
                         fontsize=INS_FS, color="#8b0000",
                         ha="center", va="bottom", zorder=5,
                         clip_on=True)

        # Junction soft-clip annotations (purple text at right edge)
        junc_fs = max(4.5, min(6.0, INS_FS * 1.1))
        for (x, yc, label) in junc_labels:
            ax_pile.text(x - 1, yc, label,
                         fontsize=junc_fs, color="#7b3294",
                         ha="right", va="center", zorder=6,
                         clip_on=True)

        # Legend
        legend_els = [
            Patch(facecolor="#4dac26",  edgecolor="none", label="Match"),
            Line2D([0], [0], color="#d7191c", lw=1.5,    label="Mismatch"),
            Patch(facecolor="#404040",  edgecolor="none", label="Deletion"),
            Line2D([0], [0], color="none", marker="",
                   label="I{n} = insertion (dark-red text)"),
            Patch(facecolor="#e4e4e4",  edgecolor="#2166ac", linewidth=1,
                  label="Unique read"),
            Patch(facecolor="#e4e4e4",  edgecolor="#f4a582", linewidth=1,
                  label="Ambiguous read"),
        ]
        if circular:
            legend_els.append(
                Line2D([0], [0], color="#7b3294", lw=2,
                       label="Purple text = junction-spanning (0 = no soft-clip)"))
        ax_pile.legend(handles=legend_els, fontsize=7, loc="upper right",
                       ncol=3, framealpha=0.9)

        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# N-fill
# ─────────────────────────────────────────────────────────────────────────────

def _revcomp(seq):
    """Reverse-complement a DNA sequence."""
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]


def _locate_flank(flank, donor_seq, max_mm=3):
    """
    Find the best-matching position of `flank` in `donor_seq` using a
    vectorised sliding-window comparison (no mappy needed; works for any
    query length including very short flanks).

    N bases in the flank are treated as wildcards (not counted as mismatches).
    Returns (r_st, r_en) of the best match in donor_seq, or None if the
    best match has more than max_mm mismatches at non-N positions.
    """
    flen = len(flank)
    if len(donor_seq) < flen:
        return None
    farr = np.frombuffer(flank.upper().encode(),     dtype=np.uint8)
    darr = np.frombuffer(donor_seq.upper().encode(), dtype=np.uint8)
    non_n = farr != ord('N')
    if non_n.sum() < 10:      # flank almost entirely N; unusable
        return None
    n_windows = len(darr) - flen + 1
    if n_windows <= 0:
        return None
    windows = np.lib.stride_tricks.sliding_window_view(darr, flen)  # (n, flen)
    mm = np.sum(windows[:, non_n] != farr[non_n], axis=1)
    best_i = int(np.argmin(mm))
    if int(mm[best_i]) > max_mm:
        return None
    return (best_i, best_i + flen)


def fill_n_regions(plasmids_orig, n_cumsums):
    """
    Fill N-masked regions in each reference using the most homologous sequence
    from another reference in the library.

    Algorithm per target reference:
      1. Find the donor reference with the highest total mlen against this one.
      2. For each N run [n_st, n_en):
           - Extract the 40 bp flanking sequences (left and right of the N run).
           - Locate each flank in the donor via sliding-window mismatch search
             (tolerates up to 3 mismatches; treats N bases as wildcards).
           - Extract donor sequence between end-of-left-flank and
             start-of-right-flank.
           - Skip if the extracted sequence is itself all-N (donor also masked).
           - Replace the N run with the extracted fill sequence.
      3. Return updated plasmids_orig dict and print a verification table.

    Parameters
    ----------
    plasmids_orig : dict {name: sequence_string}   (original, non-doubled)
    n_cumsums     : dict {name: cumsum array}       (from load_plasmid_library)

    Returns
    -------
    updated plasmids_orig dict (same keys; filled sequences where applicable)
    """
    names = list(plasmids_orig.keys())
    if len(names) < 2:
        print("  --fill-n: fewer than 2 references; no donor available. Skipping.")
        return plasmids_orig

    print("\n  --fill-n: scanning for N regions to fill …")
    updated      = dict(plasmids_orig)
    fills_applied = []
    FLANK = 40

    for tgt_name, tgt_seq in plasmids_orig.items():
        n_runs_list = _n_runs(tgt_seq)
        if not n_runs_list:
            continue

        # ── find best donor ───────────────────────────────────────────────────
        best_donor_name = None
        best_mlen       = -1
        for src_name, src_seq in plasmids_orig.items():
            if src_name == tgt_name:
                continue
            try:
                aln = mappy.Aligner(seq=src_seq, preset="map-ont", best_n=5)
                if not aln:
                    continue
                total_mlen = sum(h.mlen for h in list(aln.map(tgt_seq))[:5])
                if total_mlen > best_mlen:
                    best_mlen       = total_mlen
                    best_donor_name = src_name
            except Exception:
                continue

        if best_donor_name is None:
            print(f"    WARNING: no donor found for {tgt_name}. Skipping fills.")
            continue

        donor_seq = plasmids_orig[best_donor_name]
        print(f"    {tgt_name}  <-  donor: {best_donor_name}")

        # ── process N runs ────────────────────────────────────────────────────
        parts    = []
        prev_end = 0
        n_filled = 0

        for n_st, n_en in sorted(n_runs_list):
            run_len     = n_en - n_st
            left_flank  = tgt_seq[max(0, n_st - FLANK) : n_st]
            right_flank = tgt_seq[n_en : min(len(tgt_seq), n_en + FLANK)]

            fill_seq    = None
            skip_reason = None

            if len(left_flank) < 10 or len(right_flank) < 10:
                skip_reason = "flanking sequence too short"
            else:
                l_loc = _locate_flank(left_flank,  donor_seq)
                r_loc = _locate_flank(right_flank, donor_seq)
                if l_loc is None:
                    skip_reason = "left flank not found in donor (>3 mismatches)"
                elif r_loc is None:
                    skip_reason = "right flank not found in donor (>3 mismatches)"
                else:
                    fill_st = l_loc[1]     # donor pos right after left flank
                    fill_en = r_loc[0]     # donor pos just before right flank
                    if fill_st >= fill_en:
                        skip_reason = (f"bad fill coords ({fill_st} >= {fill_en})"
                                       f" — flanks overlap or order inverted in donor")
                    elif fill_en - fill_st > run_len * 4:
                        skip_reason = (f"fill ({fill_en - fill_st} bp) implausibly "
                                       f"long vs N run ({run_len} bp)")
                    else:
                        candidate = donor_seq[fill_st:fill_en]
                        n_in_fill = candidate.upper().count('N')
                        if n_in_fill == len(candidate):
                            skip_reason = (f"donor also has N's at coords "
                                           f"{fill_st}-{fill_en} — no fill possible")
                        elif n_in_fill > len(candidate) * 0.5:
                            skip_reason = (f"donor fill >50% N at {fill_st}-"
                                           f"{fill_en}; skipping")
                        else:
                            fill_seq = candidate

            parts.append(tgt_seq[prev_end:n_st])
            if fill_seq is not None:
                parts.append(fill_seq)
                fills_applied.append({
                    "target":       tgt_name,
                    "donor":        best_donor_name,
                    "n_st":         n_st,
                    "n_en":         n_en,
                    "run_len":      run_len,
                    "fill_len":     len(fill_seq),
                    "donor_coords": f"{fill_st}-{fill_en}",
                })
                n_filled += 1
            else:
                parts.append(tgt_seq[n_st:n_en])   # keep N run unchanged
                print(f"    SKIP N-run [{n_st}-{n_en}) ({run_len} bp): {skip_reason}")

            prev_end = n_en

        parts.append(tgt_seq[prev_end:])
        new_seq = "".join(parts)

        if n_filled > 0:
            updated[tgt_name] = new_seq
            print(f"    Applied {n_filled} fill(s) to {tgt_name}  "
                  f"(new length: {len(new_seq):,} bp, was {len(tgt_seq):,} bp)")

    # ── print verification table ──────────────────────────────────────────────
    if fills_applied:
        print()
        hdr = (f"  {'Target':<42} {'Donor':<42} "
               f"{'N-run':<20} {'Fill':>5}  Donor coords")
        print(hdr)
        print("  " + "─" * (len(hdr) - 2))
        for f in fills_applied:
            run_str = f"[{f['n_st']}-{f['n_en']}) {f['run_len']} bp"
            print(f"  {f['target']:<42} {f['donor']:<42} "
                  f"{run_str:<20} {f['fill_len']:>5} bp  {f['donor_coords']}")
        print()
        print("  *** Verify fills are biologically correct before production use. ***")
    else:
        print("  No N-fills applied (all N runs skipped).")
    print()

    return updated


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def run_pipeline(fastq_path, fasta_paths, min_qual=10.0, min_length=2000,
                 ambiguity_delta=0.05, circular=False, threads=4,
                 output_dir=".", no_plots=False, fill_n=False,
                 progress_callback=None):
    """Run the full alignment pipeline programmatically.

    Parameters
    ----------
    fastq_path : str
        Path to ONT FASTQ file (plain or .gz).
    fasta_paths : list[str]
        One or more plasmid FASTA file paths.
    min_qual : float
        Minimum average Phred quality score.
    min_length : int
        Minimum read length (bp).
    ambiguity_delta : float
        Score-difference threshold for ambiguous calls (0–1).
    circular : bool
        Treat plasmid references as circular.
    threads : int
        Number of minimap2 alignment threads.
    output_dir : str
        Directory for output files.
    no_plots : bool
        Skip plot generation.
    fill_n : bool
        Fill N-masked regions in references using the most homologous
        donor sequence from the library.
    progress_callback : callable or None
        Called as ``progress_callback(n_passed, status_message)`` during
        alignment so callers (e.g. Streamlit) can update a progress bar.

    Returns
    -------
    dict with keys:
        results_df    – pandas DataFrame of per-read assignments
        summary_df    – pandas DataFrame of per-plasmid summary
        plots_path    – str path to PDF report (or None if no_plots)
        log           – list[str] of log messages
        output_path   – str path to per-read TSV
        summary_path  – str path to summary TSV
    """
    log = []
    def _log(msg):
        log.append(msg)

    output_path  = os.path.join(output_dir, "plasmid_assignments.tsv")
    summary_path = os.path.join(output_dir, "plasmid_summary.tsv")
    plots_path   = os.path.join(output_dir, "plasmid_alignment_report.pdf") \
                   if not no_plots else None

    # ── validate inputs ──────────────────────────────────────────────────────
    if not Path(fastq_path).exists():
        raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")
    missing = [f for f in fasta_paths if not Path(f).exists()]
    if missing:
        raise FileNotFoundError("FASTA file(s) not found:\n" +
                                "\n".join(f"  {f}" for f in missing))

    _log("=" * 60)
    _log("  plasmid_align.py")
    _log("=" * 60)
    _log(f"  FASTQ:           {fastq_path}")
    _log(f"  FASTA files:     {len(fasta_paths)}")
    for f in fasta_paths:
        _log(f"                   {f}")
    _log(f"  Min avg quality: Q{min_qual}")
    _log(f"  Min length:      {min_length} bp")
    _log(f"  Ambiguity delta: {ambiguity_delta}")
    _log(f"  Circular mode:   {circular}")
    _log(f"  Fill N regions:  {fill_n}")
    _log(f"  Threads:         {threads}")
    _log("")

    # ── load plasmid library ─────────────────────────────────────────────────
    _log("Loading plasmid sequences …")
    plasmids, orig_lens, n_cumsums = load_plasmid_library(
        fasta_paths, circular=circular)
    if not plasmids:
        raise ValueError("No valid plasmid sequences loaded.")
    for name, seq in plasmids.items():
        effective_len = len(seq) // 2 if circular else len(seq)
        n_count = seq.upper().count("N")
        pct_n = 100 * n_count / len(seq)
        _log(f"  {name:40s}  {effective_len:>8,} bp  "
             f"N bases: {n_count:,} ({pct_n:.1f}%)")

    # ── optional N-fill ──────────────────────────────────────────────────────
    if fill_n:
        plasmids_orig_seqs = {name: seq[:orig_lens[name]]
                              for name, seq in plasmids.items()}
        plasmids_orig_filled = fill_n_regions(plasmids_orig_seqs, n_cumsums)
        # Rebuild plasmids (doubled if circular), orig_lens, and n_cumsums
        # from the filled sequences so all downstream code uses filled refs.
        plasmids      = {}
        new_orig_lens = {}
        new_n_cumsums = {}
        for name, seq in plasmids_orig_filled.items():
            new_orig_lens[name] = len(seq)
            cs = np.zeros(len(seq) + 1, dtype=np.int32)
            cs[1:] = np.cumsum(
                np.frombuffer(seq.upper().encode(), dtype=np.uint8) == ord('N'))
            new_n_cumsums[name] = cs
            plasmids[name] = (seq + seq) if circular else seq
        orig_lens = new_orig_lens
        n_cumsums = new_n_cumsums
        _log("Building minimap2 indices from filled sequences …")
    else:
        _log("\nBuilding minimap2 indices (map-ont preset) …")

    t0 = time.time()
    aligners = build_aligners(plasmids, threads, circular=circular)
    _log(f"  {len(aligners)} indices built in {time.time()-t0:.1f}s")

    # ── stream + align reads ─────────────────────────────────────────────────
    _log("\nAligning reads …")
    rows             = []
    read_seqs        = {}    # cache for pileup: {read_id: sequence}
    pre_filter_rows  = []    # all reads (incl. filtered) for QC scatter plot
    n_total       = 0
    n_short       = 0
    n_lowqual     = 0
    n_passed      = 0
    n_unmapped    = 0
    n_ambiguous   = 0

    t0 = time.time()

    for name, seq, qual in iter_fastq(fastq_path):
        n_total += 1
        rlen  = len(seq)
        avg_q = avg_phred(qual)   # compute for all reads (needed for QC scatter)
        if rlen < min_length:
            n_short += 1
            pre_filter_rows.append({"read_id": name, "read_length": rlen,
                                    "avg_qual": round(avg_q, 2),
                                    "_scatter_cat": "filtered_short"})
            continue
        if avg_q < min_qual:
            n_lowqual += 1
            pre_filter_rows.append({"read_id": name, "read_length": rlen,
                                    "avg_qual": round(avg_q, 2),
                                    "_scatter_cat": "filtered_lowqual"})
            continue
        pre_filter_rows.append({"read_id": name, "read_length": rlen,
                                 "avg_qual": round(avg_q, 2),
                                 "_scatter_cat": None})  # filled after alignment
        n_passed += 1
        if n_passed % 500 == 0:
            elapsed = time.time() - t0
            msg = f"    … {n_passed} reads aligned ({elapsed:.0f}s)"
            _log(msg)
            if progress_callback:
                progress_callback(n_passed, msg)

        hit_map = align_read(seq, aligners, circular=circular,
                             orig_lens=orig_lens, n_cumsums=n_cumsums)
        best, ambig, second, reason = classify_read(hit_map, ambiguity_delta)

        if best == "unmapped":
            n_unmapped += 1
        if ambig:
            n_ambiguous += 1

        # Cache sequence for pileup generation
        if best != "unmapped":
            read_seqs[name] = seq

        # Determine if the best-plasmid hit is junction-spanning
        best_hit = hit_map.get(best) if best != "unmapped" else None
        is_jspan = (best_hit.get("is_junction_spanning", False)
                    if best_hit else False)

        row = {
            "read_id":              name,
            "read_length":          len(seq),
            "avg_qual":             round(avg_q, 2),
            "best_plasmid":         best,
            "is_ambiguous":         ambig,
            "is_junction_spanning": is_jspan,
            "second_plasmid":       second if second else "",
            "reason":               reason,
        }
        # per-plasmid columns
        for pname in plasmids:
            h = hit_map.get(pname)
            pfx = pname.replace(" ", "_").replace("/", "_").replace(":", "_")
            if h:
                row[f"{pfx}__match_frac"] = h["match_frac"]
                row[f"{pfx}__identity"]   = h["identity"]
                row[f"{pfx}__query_cov"]  = h["query_cov"]
                row[f"{pfx}__mapq"]       = h["mapq"]
                row[f"{pfx}__n_in_hit"]   = h.get("n_in_hit", 0)
            else:
                row[f"{pfx}__match_frac"] = 0.0
                row[f"{pfx}__identity"]   = 0.0
                row[f"{pfx}__query_cov"]  = 0.0
                row[f"{pfx}__mapq"]       = 0
                row[f"{pfx}__n_in_hit"]   = 0

        rows.append(row)

    # ── annotate pre-filter rows with final assignment category ─────────────
    assign_map = {row["read_id"]: (row["best_plasmid"], row["is_ambiguous"])
                  for row in rows}
    for pr in pre_filter_rows:
        if pr["_scatter_cat"] is not None:
            continue
        bp, is_ambig = assign_map.get(pr["read_id"], ("unmapped", False))
        if bp == "unmapped":
            pr["_scatter_cat"] = "unmapped"
        elif is_ambig:
            pr["_scatter_cat"] = "ambiguous"
        else:
            pr["_scatter_cat"] = bp
    pre_filter_df = pd.DataFrame(pre_filter_rows)

    # ── write per-read results ───────────────────────────────────────────────
    results_df = pd.DataFrame(rows)
    results_df.to_csv(output_path, sep="\t", index=False)
    _log(f"\nPer-read results → {output_path}")

    # ── per-plasmid summary ──────────────────────────────────────────────────
    summary_rows = []
    for pname in plasmids:
        assigned   = results_df[
            (results_df["best_plasmid"] == pname) & ~results_df["is_ambiguous"]
        ]
        ambig_sub  = results_df[
            (results_df["best_plasmid"] == pname) & results_df["is_ambiguous"]
        ]
        pfx = pname.replace(" ", "_").replace("/", "_").replace(":", "_")
        mf_col = f"{pfx}__match_frac"
        summary_rows.append({
            "plasmid":                pname,
            "assigned_reads":         len(assigned),
            "ambiguous_reads":        len(ambig_sub),
            "total_reads_best_hit":   len(assigned) + len(ambig_sub),
            "median_match_frac":      round(results_df[mf_col].median(), 4)
                                      if mf_col in results_df else None,
            "median_identity":        round(results_df[f"{pfx}__identity"].median(), 4)
                                      if f"{pfx}__identity" in results_df else None,
            "median_query_cov":       round(results_df[f"{pfx}__query_cov"].median(), 4)
                                      if f"{pfx}__query_cov" in results_df else None,
        })

    summary_df = pd.DataFrame(summary_rows).sort_values(
        "assigned_reads", ascending=False
    )
    summary_df.to_csv(summary_path, sep="\t", index=False)
    _log(f"Per-plasmid summary → {summary_path}")

    # ── log summary ──────────────────────────────────────────────────────────
    _log("")
    _log("=" * 60)
    _log("  SUMMARY")
    _log("=" * 60)
    _log(f"  Total reads in FASTQ:       {n_total:>8,}")
    _log(f"  Filtered (too short):       {n_short:>8,}")
    _log(f"  Filtered (low quality):     {n_lowqual:>8,}")
    _log(f"  Reads passed filters:       {n_passed:>8,}")
    _log(f"  Unmapped:                   {n_unmapped:>8,}")
    _log(f"  Ambiguous:                  {n_ambiguous:>8,}")
    _log(f"  Uniquely assigned:          {n_passed - n_unmapped - n_ambiguous:>8,}")
    if circular and "is_junction_spanning" in results_df.columns:
        n_jspan = results_df["is_junction_spanning"].sum()
        _log(f"  Junction-spanning:          {n_jspan:>8,}")
    _log("")
    _log("  Per-plasmid assignments (unique):")
    for _, srow in summary_df.iterrows():
        _log(f"    {srow['plasmid']:45s}  "
             f"unique={srow['assigned_reads']:>5}  "
             f"ambig={srow['ambiguous_reads']:>4}")
    _log("=" * 60)

    # ── plots ────────────────────────────────────────────────────────────────
    if not no_plots:
        _log("\nGenerating visual report …")
        if progress_callback:
            progress_callback(n_passed, "Generating plots …")
        counts = dict(
            n_total    = n_total,
            n_short    = n_short,
            n_lowqual  = n_lowqual,
            n_passed   = n_passed,
            n_unmapped = n_unmapped,
            n_ambiguous= n_ambiguous,
        )
        with PdfPages(plots_path) as pdf:
            generate_readqc_plot(pre_filter_df, list(plasmids.keys()),
                                 min_length, min_qual, pdf)
            generate_plots(results_df, list(plasmids.keys()), counts, pdf,
                           orig_lens=orig_lens)
            _log("  Generating per-plasmid pileup plots …")
            plasmid_seqs_orig = {name: seq[:len(seq)//2] if circular else seq
                                 for name, seq in plasmids.items()}
            generate_pileup_plots(
                results_df, plasmid_seqs_orig,
                aligners, read_seqs, pdf,
                circular=circular, orig_lens=orig_lens,
            )
            d = pdf.infodict()
            d["Title"]   = "Plasmid Alignment Report"
            d["Subject"] = "ONT read alignment to plasmid library"
        _log(f"Visual report  → {plots_path}")

    return {
        "results_df":   results_df,
        "summary_df":   summary_df,
        "plots_path":   plots_path,
        "log":          log,
        "output_path":  output_path,
        "summary_path": summary_path,
    }


def main():
    args = parse_args()

    # Use a temp directory for pipeline outputs, then move to user-specified paths
    import tempfile
    import shutil
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            result = run_pipeline(
                fastq_path=args.fastq,
                fasta_paths=args.fasta,
                min_qual=args.min_qual,
                min_length=args.min_length,
                ambiguity_delta=args.ambiguity_delta,
                circular=args.circular,
                threads=args.threads,
                output_dir=tmpdir,
                no_plots=args.no_plots,
                fill_n=args.fill_n,
                progress_callback=lambda n, msg: print(msg, flush=True),
            )
        except (FileNotFoundError, ValueError) as exc:
            sys.exit(f"ERROR: {exc}")

        # Print the collected log to stdout
        for line in result["log"]:
            print(line)

        # Move outputs to user-specified paths
        shutil.move(result["output_path"], args.output)
        shutil.move(result["summary_path"], args.summary)
        if not args.no_plots and result["plots_path"]:
            shutil.move(result["plots_path"], args.plots_output)


if __name__ == "__main__":
    main()
