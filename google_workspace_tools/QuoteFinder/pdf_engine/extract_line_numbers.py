#!/usr/bin/env python3
"""
extract_line_numbers.py
=======================
Read a Google-Docs-exported PDF that has "Show line numbers" turned on, and
emit  rendered_lines.json = [{"n": <line number>, "text": "<line text>"}, ...]
in document order.

This feeds the "PDF/JSON (Exact Line #)" engine of the QuoteFinder Apps Script:
it reads the *real* numbers Google printed in the left margin, so it is correct
for any layout. No Google login / API access is needed here — you feed it the PDF.

------------------------------------------------------------------------------
HOW TO PRODUCE THE PDF (important):
  In the manuscript doc:  Tools/Format menu -> turn ON "Show line numbers"
  (continuous, whole document), then  File -> Download -> PDF Document (.pdf).
  The downloaded PDF must visibly show the line numbers in the left margin.

USAGE (local):
  pip install pdfplumber
  python extract_line_numbers.py manuscript.pdf rendered_lines.json
  # then upload rendered_lines.json to your Google Drive, and point the
  # Apps Script at it:  QuoteFinder -> "Set PDF/JSON line-number file…"

USAGE (Google Colab):
  # cell 1
  !pip -q install pdfplumber
  # cell 2: use the file picker to upload manuscript.pdf
  from google.colab import files; files.upload()
  # cell 3
  !python extract_line_numbers.py manuscript.pdf rendered_lines.json
  # cell 4: download the result, then drag it into Google Drive
  from google.colab import files; files.download('rendered_lines.json')
------------------------------------------------------------------------------
"""

import sys
import re
import json

try:
    import pdfplumber
except ImportError:
    sys.exit("pdfplumber is not installed.  Run:  pip install pdfplumber")

NUM_RE = re.compile(r'^\d{1,6}$')
TOP_TOL = 3.0          # points: words within this vertical distance share a line
MARGIN_GAP = 2.0       # points: a number this far left of the text column = line #


def extract(pdf_path):
    lines = []
    warnings = []
    with pdfplumber.open(pdf_path) as pdf:
        for pageno, page in enumerate(pdf.pages, start=1):
            words = page.extract_words(use_text_flow=False, keep_blank_chars=False)
            if not words:
                continue

            # The text column's left edge = leftmost x0 among NON-numeric words.
            # Line numbers are numeric tokens sitting to the LEFT of that edge.
            non_num_x = [w["x0"] for w in words if not NUM_RE.match(w["text"])]
            if not non_num_x:
                continue
            body_left = min(non_num_x)

            # Cluster words into visual lines by their top (y) coordinate.
            words.sort(key=lambda w: (round(w["top"], 1), w["x0"]))
            groups, cur, cur_top = [], [], None
            for w in words:
                if cur_top is None or abs(w["top"] - cur_top) <= TOP_TOL:
                    cur.append(w)
                    cur_top = w["top"] if cur_top is None else cur_top
                else:
                    groups.append(cur)
                    cur, cur_top = [w], w["top"]
            if cur:
                groups.append(cur)

            for g in groups:
                g.sort(key=lambda w: w["x0"])
                num = None
                body_words = []
                for w in g:
                    if num is None and NUM_RE.match(w["text"]) and w["x0"] < body_left - MARGIN_GAP:
                        num = int(w["text"])          # the margin line number
                    else:
                        body_words.append(w)
                text = " ".join(w["text"] for w in body_words).strip()
                if num is not None:
                    lines.append({"n": num, "text": text})
                elif text:
                    # a visual line with text but no detected margin number —
                    # attach it to the previous numbered line so the quote text
                    # still matches (rare: e.g. number column mis-detected).
                    if lines:
                        lines[-1]["text"] = (lines[-1]["text"] + " " + text).strip()
                    else:
                        warnings.append("Page %d: text before any line number." % pageno)

    # ---- sanity checks (printed, not fatal) ---------------------------------
    if not lines:
        sys.exit("No line numbers were detected. Is 'Show line numbers' ON in the "
                 "PDF, and is it a single-column layout?")
    nums = [l["n"] for l in lines]
    if nums != sorted(nums):
        warnings.append("Line numbers are not strictly increasing — the doc may "
                        "have a complex layout; spot-check the output.")
    # drop accidental exact duplicates (same n AND same text)
    seen, deduped = set(), []
    for l in lines:
        key = (l["n"], l["text"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(l)

    return deduped, warnings, (min(nums), max(nums))


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python extract_line_numbers.py <input.pdf> [output.json]")
    pdf_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else "rendered_lines.json"

    lines, warnings, (lo, hi) = extract(pdf_path)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(lines, f, ensure_ascii=False)

    print("Wrote %d lines to %s  (line numbers %d…%d)" % (len(lines), out_path, lo, hi))
    print("First lines:")
    for l in lines[:3]:
        print("  %5d | %s" % (l["n"], (l["text"][:70])))
    print("Last lines:")
    for l in lines[-3:]:
        print("  %5d | %s" % (l["n"], (l["text"][:70])))
    for w in warnings:
        print("WARNING:", w)


if __name__ == "__main__":
    main()
