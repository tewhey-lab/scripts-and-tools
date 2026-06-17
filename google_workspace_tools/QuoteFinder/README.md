# QuoteFinder — refresh manuscript line references in a response document

You have two Google Docs:

- a **primary / manuscript** doc with Google Docs **line numbers** turned on (continuous), and
- a **response** doc containing markers like:

  ```
  @AUTOLINE123: "some quoted text from the manuscript"
  @AUTOLINE45-58: "a longer quote that wraps across several lines"
  ```

For each marker, QuoteFinder finds the quoted text in the manuscript, rewrites the
marker to the **current** line number(s), and copies basic formatting
(**font color, bold, italic, superscript/subscript**) from the manuscript onto
the quote in the response doc.

Single rendered line → `@AUTOLINE<start>:` Spans several lines → `@AUTOLINE<start>-<end>:`

(The in-document marker token stays `@AUTOLINE` — only the app/menu is named QuoteFinder.)

---

## The one hard part (read this)

Google **does not expose its rendered line numbers** to any script or API — they
exist only when the page is laid out. So there are two ways to get them, and
this repo gives you **both** so you can compare on your real document:

| | **gDoc (Estimated Line #)** | **PDF/JSON (Exact Line #)** |
|---|---|---|
| Setup | none | run a small Python script once per revision |
| Per‑use | one click in the doc | one click in the doc (after refreshing the JSON) |
| Accuracy | exact for ordinary flowing text; **flags** anything near tables/images instead of guessing | exact for any layout (reads Google's real numbers) |
| Login | none beyond the normal Apps Script prompt | none (you hand it the PDF) |

Everything else (finding quotes, copying formatting, rewriting markers) is
identical for both. Since the PDF/JSON numbers are ground truth, use them when you
can; **Preview matches** shows what either engine would write before you commit.

---

## Install (one time, ~3 min)

1. Open your **response** doc → **Extensions → Apps Script**.
2. Delete the default `Code.gs`, paste in [`apps_script/Code.gs`](apps_script/Code.gs), **Save**.
3. Add the sidebar UI: **＋ → HTML**, name it exactly **`ReviewSidebar`** (no `.html`),
   replace its contents with [`apps_script/ReviewSidebar.html`](apps_script/ReviewSidebar.html), **Save**.
4. Reload the response doc. A new **QuoteFinder** menu appears.
5. **QuoteFinder → Settings…** → paste the manuscript doc's URL or ID (and, for the
   exact engine, the PDF/JSON line-number file). Save.
6. First run will show a Google authorization prompt (it reads the manuscript and edits this doc). Approve it.
7. *(Optional)* Only needed for the marker-style option **“Use the document’s Normal text style”**
   (Settings → *Marker text style*). The default style mode is a fixed **set style**
   (Arial 11 #000000, regular), which needs nothing extra. To use the Normal-text option, enable
   the **Google Docs API** advanced service: in the Apps Script editor, next to **Services** click
   **＋**, pick **Google Docs API**, identifier **`Docs`**, **Add**. This **Add** action only
   appends the service — it's the safe path. The repo's
   [`apps_script/appsscript.json`](apps_script/appsscript.json) shows the exact entry
   (`enabledAdvancedServices` → `Docs`/`docs`/`v1`); add just that block if you edit the manifest
   by hand — **don't paste the whole file over an existing manifest**, since the manifest is
   replaced wholesale and you'd wipe your `timeZone`, `oauthScopes`, and any other settings. Then
   re-run once to approve the added permission.

The **QuoteFinder** menu:
- **Update references (all at once)** ▸ gDoc (Estimated Line #) / PDF/JSON (Exact Line #) — fire-and-forget; flags imperfect ones in the doc.
- **Review & update (step through imperfect)** ▸ gDoc (Estimated Line #) / PDF/JSON (Exact Line #) — the interactive panel (below).
- **Preview matches** — dry run, no changes (uses the exact PDF/JSON numbers if a file is loaded, else the estimate).

---

## Method 1 — gDoc (Estimated Line #), zero setup

**QuoteFinder → Update references (all at once) → gDoc (Estimated Line #).**

- Updated markers get their new numbers.
- **Yellow flag** = low‑confidence number (a table/image/unsupported font sits
  somewhere above it, so the absolute count may be off — check these against the
  manuscript). **Red flag** = quote not found (typo or edited wording). These are
  background shades the tool paints to mark results — not scope input.
- Re‑run any time. **Clear all QuoteFinder flags** removes the marks.

> This engine reproduces Google's wrapping using standard Arial/Times metrics.
> For plain single‑column manuscript text it is typically exact; it is *not*
> guaranteed exact once tables, inline images, or unusual fonts appear above a
> quote — that's what the yellow flag is for, and what PDF/JSON (Exact Line #) fixes.

---

## Method 2 — PDF/JSON (Exact Line #), from the PDF

Do this whenever the manuscript changes:

1. In the **manuscript** doc: turn **ON** “Show line numbers” (continuous, whole
   document), then **File → Download → PDF Document (.pdf)**. The PDF must show
   the numbers in the left margin.
2. Turn that PDF into `rendered_lines.json` — pick whichever is easier:
   - **No install (browser):** open [`pdf_engine/pdf_to_json.html`](pdf_engine/pdf_to_json.html)
     by double-clicking it, drop the PDF on the page, and click **Download
     rendered_lines.json**. (Needs internet — it loads PDF.js from a CDN.)
   - **Python (local or Colab):** see the header of
     [`pdf_engine/extract_line_numbers.py`](pdf_engine/extract_line_numbers.py):
     ```bash
     pip install pdfplumber
     python extract_line_numbers.py manuscript.pdf rendered_lines.json
     ```
   Both produce the same JSON.
3. **Upload `rendered_lines.json` to your Google Drive** (just drag it in).
4. In the response doc: **QuoteFinder → Settings…** → put the file’s name
   (`rendered_lines.json`) or its Drive ID in the *PDF/JSON line-number file* field.
   (Only needed once unless you rename it.)
5. **QuoteFinder → Update references (all at once) → PDF/JSON (Exact Line #).**

> **How matching works (robust by design):** the gDoc layout gives a *ballpark* line
> number; the PDF then **fine‑tunes** it. For each quote the tool searches the PDF
> text near that ballpark — exact first, then a **fuzzy** match within a window — so
> small PDF‑extraction quirks (e.g. `ﬁ` ligatures, superscript spacing) don't break
> it. If the PDF genuinely can't pin a quote, it **falls back to the gDoc estimate**
> and flags it **yellow / "estimated"** (so a quote found in the doc never comes back
> as "no usable line number" — it just may be flagged for a manual check). Tunable via
> `PDF_WINDOW` and `PDF_FUZZY_MIN` in `Code.gs`.

---

## Reviewing matches (the review panel)

When a quote no longer matches the manuscript word-for-word (it was edited,
shortened, or has a typo), fire-and-forget mode just flags it. To resolve those
case-by-case: **QuoteFinder → Review & update (step through imperfect) ▸ gDoc
(Estimated Line #)** (or **PDF/JSON (Exact Line #)**). A wide, **draggable** panel
opens that **floats over the doc but doesn't block it** — you can still scroll and
edit while it's open.

- Three **category tabs** at the top: **Changes** (default — anything that will change
  on apply: edited/typo'd text, not-found/malformed markers, **and** exact matches whose
  line number is moving), **Line # only** (a subset — exact-text matches whose *only*
  change is the line number, incl. low-confidence "estimated" lines), and **All** (every
  marker in scope, including the no-op perfect ones). The tabs only filter what you step
  through; they never change what **Apply** does.
- The **decision and navigation buttons sit above the quote text** (Use manuscript / Skip,
  and ◀ Prev / Next ▶), so they stay put even when a long quote grows the card.
- Each card shows a **match %** badge, the manuscript line(s) the candidate sits on,
  **your current quote** and **the manuscript text now** *side by side, with their
  real formatting* (color, bold, italic, super/subscript), and a **word-level diff**.
- For each, choose **Use manuscript ✓** / **Update ✓** or **Skip**. Defaults:
  perfect matches → refresh number + formatting; high-confidence fuzzy (≥90%) → use;
  everything else → skip. Use **◀ / ▶** to move and revisit.
- **Apply** commits everything in one pass — "*N of M will change*" — then rescans
  so anything you skipped is still there. Nothing is written until you click Apply.
- **⚡ Apply this one now ✓** commits *just the current card* immediately and shows
  "✓ Applied — line N." It re-finds the marker fresh on the server (by its order in
  the doc), so committing one **can't** shift the positions of the others. Note it
  re-reads the manuscript each time, so it costs about as much as a full Apply — use
  it to land one or two and watch them; use **Apply** for doing many at once.
- Items with **no match found** can only be skipped (there's nothing to pull in).

**Re-reading:** the panel scans when it opens. Edited the doc? Click **↻ Rescan** to
re-read (your Use/Skip picks are kept). To change which lines are **in scope**, reopen
the panel with a different drag-selection — Rescan keeps the selection you opened with
(see *Scanning a selection* below).

> **Width / docked vs floating:** the panel is a wide floating dialog by default
> (`REVIEW_WIDTH` ≈ 820 px in `Code.gs`). Set `REVIEW_PANEL: 'sidebar'` to use the
> classic docked sidebar instead — but note Google **locks sidebars to ~300 px**,
> so that can't be widened.

---

## Preview matches (dry run)

**QuoteFinder → Preview matches (dry run, no changes)** opens a table
`# | old | quote | new #` showing the line number each marker *would* get — with
**nothing written**. It uses the **PDF/JSON exact numbers when a line-number file
is loaded** (ground truth), otherwise the gDoc estimate; the heading names which
one, and quotes with no verbatim match are flagged. (Fuzzy candidates aren't shown
here — that's what **Review & update** is for.)

---

## Scanning a subset — select the lines

To process only some quotes, **drag-select the quote line(s)** and run a command — the
selection *is* the scope. **No toggle, no highlighting.** Nothing selected → the whole
document. Selection is your live cursor selection, so it's read instantly (no save‑lag).

- **Update references** / **Preview matches** read your **live selection** at run time.
  The report says e.g. *"Scope: 6 selected quote(s)."* Nothing selected → whole doc.
- **Review & update** can't read a live selection (it's a floating dialog), so it
  **captures your selection when you open the panel** and uses it for the whole
  session (scan, Rescan, Apply). To change the subset, **reopen** the panel with a
  different selection. The panel shows a blue **"selection: N quote(s)"** badge.
- It's **coarse by line**: selecting any part of a line includes that line's marker.
- **Identical quotes are matched by text.** If the same quote appears on two lines and
  you select one, *both* are processed (they'd get the same line number anyway). The
  reports show the in-scope count vs the distinct-selected count so this is visible,
  never silent. (Likewise, several bare/unfinished `@AUTOLINE` hooks can group together.)
- A selection that contains **no markers** processes **nothing** (and says so) — it
  does *not* fall back to the whole document.

## Marker format (configurable)

A marker is simply the **hook** followed, within a short distance, by a **quote**:

```
<HOOK> … <OPEN>quote<CLOSE>        default:   @AUTOLINE: "…"      (also: @AUTOLINE "…", @AUTOLINE12-15: "…")
```

QuoteFinder **only looks for the hook + a nearby quote** — it does *not* require (or
parse) a number or colon. Anything between the hook and the opening quote (an old
`12-15:`, a stray space, etc.) is ignored on input and replaced on output.

Edit these in **QuoteFinder → Settings…** (no code needed — saved per document), or
as defaults at the top of `Code.gs`:
- **`HOOK`** — the token that starts a marker (default `@AUTOLINE`). Any characters
  (regex-special ones are handled).
- **`QUOTE_OPEN` / `QUOTE_CLOSE`** — the marks that bracket the quote (default `"`/`"`;
  straight or curly double quotes both accepted). Can be anything, e.g. `«`/`»` or `[[`/`]]`.
- **`OPEN_LOOKAHEAD`** — how far **past the hook** to look for the **opening** mark
  (default 20 characters). Whitespace is always allowed; this bounds any *other*
  characters (like an old `12-15:`) between the hook and where the quote starts.
- **`MAX_LOOKAHEAD`** — the **quote body** length: how far from the **opening** mark to
  scan for the **closing** mark (default 4000). For a split quote (`[...]`, below) this
  covers the whole span, so it runs through to the **last** piece's closing mark.

So there are two windows: a small one (`OPEN_LOOKAHEAD`) to find where the quote
*starts*, and a large one (`MAX_LOOKAHEAD`) for the quote's *length*.

Notes:
- Matching ignores case and collapses whitespace; en/em dashes are normalized.
- On **auto Update** and **Preview** the quote text is left exactly as you wrote it
  (only its formatting and the line number are updated). Only **Review** can change the
  text, and only when you pick **“Use manuscript”** on a fuzzy match — then the body is
  replaced with the manuscript wording (see *Reviewing matches* and *Split quotes*).
  The marker itself is rewritten as **`<HOOK><start>: `** or **`<HOOK><start>-<end>: `**
  (the found line number(s), colon, **one space**) regardless of what was there before.
- A hook with **no quote within `OPEN_LOOKAHEAD`** is flagged as malformed (red), so
  you can spot it.

### Split quotes (elisions)

A quote that contains **`[...]`** (also `[…]` or `[. . .]`) is treated as a quote
with an omission: each side is located **independently**, and the reported range
runs from the **start of the first piece to the end of the last piece**.

```
@AUTOLINE0: "the cells were treated [...] and then imaged"
   →  finds "the cells were treated" (say line 12) and "and then imaged" (line 18)
   →  @AUTOLINE12-18: "the cells were treated [...] and then imaged"
```

On **auto Update** and **Preview**, your text (every piece *and* the `[...]`
separators) is kept exactly as written — only the line number and per-piece formatting
are updated.

In **Review**, picking **“Use manuscript”** on a fuzzy split quote rewrites only the
pieces that *differ* from the manuscript to the manuscript wording (e.g. dropping a stray
trailing period), while **exact pieces, the `[...]` separators, and the spacing around
each piece are preserved**. An all-exact split quote is left untouched. If a corrected
piece's manuscript text happens to contain its own `[...]` or a closing quote mark, that
piece is kept as-you-wrote-it instead (splicing it in would break the marker on the next
scan). One small caveat: any character formatting you applied *to the `[...]` separators
themselves* is not preserved across the rewrite (the separators keep their text, not their
styling); the pieces' own formatting is copied piece-by-piece as before.

### Malformed markers (flagged)

If a hook + number + colon is found but the quote is **missing**, **unterminated**
(no closing mark within `MAX_LOOKAHEAD`), or **empty**, the marker can't be matched.
QuoteFinder flags it instead of silently ignoring it, and **shows the hook plus the
~20 characters that follow it** so you can find the broken marker (and the passage it
was pointing at) in the document:
- **Auto** runs **highlight the hook red** (like a not-found quote) and the report
  lists each one as `MALFORMED (reason): @AUTOLINE… <next 20 chars> …`.
- **Review** lists it with a red **“malformed”** badge and the reason, and shows
  `@AUTOLINE` + the next 20 characters under *“Marker in this doc.”* It can only be
  skipped (and its hook is highlighted red on Apply). Fix the marker in the document.
- **Preview matches** shows the same hook + 20 characters in its first column.

## Settings

The document IDs, the marker-format keys (`HOOK`, `QUOTE_OPEN`, `QUOTE_CLOSE`,
`OPEN_LOOKAHEAD`, `MAX_LOOKAHEAD`), and the **marker text style** (below) are editable
from **QuoteFinder → Settings…** (saved per document; an empty field reverts to the
default). The rest are defaults at the top of `Code.gs`:

| Key | Default | Meaning |
|---|---|---|
| `HOOK` | `'@AUTOLINE'` | marker prefix that starts a reference |
| `QUOTE_OPEN` / `QUOTE_CLOSE` | `'"'` / `'"'` | marks that bracket the quote (curly accepted; any literal allowed) |
| `OPEN_LOOKAHEAD` | `20` | max chars between the colon and the **opening** mark (whitespace always allowed) |
| `MAX_LOOKAHEAD` | `4000` | quote-body length: **opening** mark → **closing** mark |
| `COPY_FORMATTING` | `true` | copy color/bold/italic/super‑sub onto the quote |
| `ALWAYS_RANGE` | `false` | `true` → always emit `a-b` even for one line |
| `CASE_INSENSITIVE` | `true` | match quotes ignoring case |
| `MARKER_STYLE_MODE` | `'set'` | how the marker scaffold (hook + line # + colon + quote marks — **not** the quote) is styled when updating: `'insert'` = keep the formatting where it sits (no restyling), `'named'` = the document-wide *Normal text* style (needs the Docs API — see Install step 7), `'set'` = the explicit values below. Pick it with the radio buttons in Settings → *Marker text style* |
| `MARKER_FONT_FAMILY` / `MARKER_FONT_SIZE` / `MARKER_FONT_COLOR` | `'Arial'` / `11` / `'#000000'` | the **set** style's font/size/color (blank fields fall back to these defaults) |
| `MARKER_BOLD` / `MARKER_ITALIC` | `false` / `false` | the **set** style's weight (default is regular) |
| `PDF_WINDOW` | `40` | PDF/JSON engine: rendered lines to search each side of the gDoc estimate |
| `PDF_FUZZY_MIN` | `0.5` | PDF/JSON engine: min word-similarity to accept a fuzzy PDF refinement |
| `REVIEW_PANEL` | `'dialog'` | review UI: `'dialog'` (wide, floating, draggable) or `'sidebar'` (docked, ~300 px) |
| `REVIEW_WIDTH` / `REVIEW_HEIGHT` | `820` / `780` | review dialog size in px (dialog mode only) |

## Known limits / things that get flagged

- A quote that appears **more than once** in the manuscript → uses the first
  occurrence and reports “multiple matches”.
- A quote that **spans paragraphs** in the manuscript → line number still works
  (PDF/JSON engine) but formatting copy is skipped and reported.
- Two‑column layouts and line‑numbering that **restarts per page** aren’t
  supported by the gDoc (Estimated Line #) engine (the PDF/JSON engine assumes single‑column).
- Font color copy treats the manuscript’s default as black (`#000000`).

> Heads up: this was written to run in your Google environment; I couldn’t
> execute it here. Try it on a **copy** of the response doc first.
