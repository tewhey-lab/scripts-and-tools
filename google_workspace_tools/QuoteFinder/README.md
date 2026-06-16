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

The **QuoteFinder** menu:
- **Update references (all at once)** ▸ gDoc (Estimated Line #) / PDF/JSON (Exact Line #) — fire-and-forget; flags imperfect ones in the doc.
- **Review & update (step through imperfect)** ▸ gDoc (Estimated Line #) / PDF/JSON (Exact Line #) — the interactive panel (below).
- **Preview matches** — dry run, no changes (uses the exact PDF/JSON numbers if a file is loaded, else the estimate).

---

## Method 1 — gDoc (Estimated Line #), zero setup

**QuoteFinder → Update references (all at once) → gDoc (Estimated Line #).**

- Updated markers get their new numbers.
- **Yellow** highlight = low‑confidence number (a table/image/unsupported font sits
  somewhere above it, so the absolute count may be off — check these against the
  manuscript). **Red** = quote not found (typo or edited wording).
- Re‑run any time. **Clear all QuoteFinder highlights** removes the marks.

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

---

## Reviewing matches (the review panel)

When a quote no longer matches the manuscript word-for-word (it was edited,
shortened, or has a typo), fire-and-forget mode just flags it. To resolve those
case-by-case: **QuoteFinder → Review & update (step through imperfect) ▸ gDoc
(Estimated Line #)** (or **PDF/JSON (Exact Line #)**). A wide, **draggable** panel
opens that **floats over the doc but doesn't block it** — you can still scroll and
edit while it's open.

- A **view toggle** at the top: **Only changes** (default — fuzzy / not-found
  markers) or **All** (also lists the perfect matches so you can eyeball/veto them).
- Each card shows a **match %** badge, the manuscript line(s) the candidate sits on,
  **your current quote** and **the manuscript text now** *side by side, with their
  real formatting* (color, bold, italic, super/subscript), and a **word-level diff**.
- For each, choose **Use manuscript ✓** / **Update ✓** or **Skip**. Defaults:
  perfect matches → refresh number + formatting; high-confidence fuzzy (≥90%) → use;
  everything else → skip. Use **◀ / ▶** to move and revisit.
- **Apply** commits everything in one pass — "*N of M will change*" — then rescans
  so anything you skipped is still there. Nothing is written until you click Apply.
- Items with **no match found** can only be skipped (there's nothing to pull in).

> **Width / docked vs floating:** the panel is a wide floating dialog by default
> (`REVIEW_WIDTH` ≈ 820 px in `Code.gs`). Set `REVIEW_PANEL: 'sidebar'` to use the
> classic docked sidebar instead — but note Google **locks sidebars to ~300 px**,
> so that can't be widened.
>
> Avoid hand-editing the markers while the panel is open — click **Rescan** if you
> do, so its positions stay in sync.

---

## Preview matches (dry run)

**QuoteFinder → Preview matches (dry run, no changes)** opens a table
`# | old | quote | new #` showing the line number each marker *would* get — with
**nothing written**. It uses the **PDF/JSON exact numbers when a line-number file
is loaded** (ground truth), otherwise the gDoc estimate; the heading names which
one, and quotes with no verbatim match are flagged. (Fuzzy candidates aren't shown
here — that's what **Review & update** is for.)

---

## Scanning only highlighted markers

To work on just a subset, use the **scope toggle** in the QuoteFinder menu. It
always names the current mode and updates immediately (no reload):
- **○ Scope: WHOLE document — click to scan highlighted text only**
- **◉ Scope: HIGHLIGHTED text only — click to scan the whole document**

With highlighted scope on, **every** command — Update, Review, Preview — ignores any
`@AUTOLINE` marker unless its marker/quote text has a **highlight (background) color**.

- Highlight the marker or its quote with the highlighter tool (any color). Scattered
  markers are fine — it's not a single drag-selection, and it persists across runs.
- The tool's own yellow/red flag shades don't count as "highlighted," so prior runs
  won't accidentally widen the scope.
- In this mode the tool **won't add or clear background colors**, so your highlights
  are left intact (results are reported in the summary / sidebar instead).
- Turn it off to go back to the whole document.

> "Highlighted" here means a **highlight color**, not the blue drag-selection. Ask
> if you'd prefer it to use the current selection instead.

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
- The quote text is left exactly as you wrote it (only its formatting is updated).
  The marker is rewritten as **`<HOOK><start>: `** or **`<HOOK><start>-<end>: `** (the
  found line number(s), colon, **one space**) regardless of what was there before.
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

Your elided text is kept as-is (it is **never** replaced by the manuscript text);
formatting is copied piece-by-piece. Works in auto, preview, and review.

### Malformed markers (flagged)

If a hook + number + colon is found but the quote is **missing**, **unterminated**
(no closing mark within `MAX_LOOKAHEAD`), or **empty**, the marker can't be matched.
QuoteFinder flags it instead of silently ignoring it:
- **Auto** runs **highlight the hook red** (like a not-found quote) and report a count.
- **Review** lists it with a red **“malformed”** badge and the reason; it can only be
  skipped (and its hook is highlighted red on Apply). Fix the marker in the document.

## Settings

The document IDs and the marker-format keys (`HOOK`, `QUOTE_OPEN`, `QUOTE_CLOSE`,
`OPEN_LOOKAHEAD`, `MAX_LOOKAHEAD`) are editable from **QuoteFinder → Settings…** (saved
per document; an empty field reverts to the default). The rest are defaults at the top of `Code.gs`:

| Key | Default | Meaning |
|---|---|---|
| `HOOK` | `'@AUTOLINE'` | marker prefix that starts a reference |
| `QUOTE_OPEN` / `QUOTE_CLOSE` | `'"'` / `'"'` | marks that bracket the quote (curly accepted; any literal allowed) |
| `OPEN_LOOKAHEAD` | `20` | max chars between the colon and the **opening** mark (whitespace always allowed) |
| `MAX_LOOKAHEAD` | `4000` | quote-body length: **opening** mark → **closing** mark |
| `COPY_FORMATTING` | `true` | copy color/bold/italic/super‑sub onto the quote |
| `ALWAYS_RANGE` | `false` | `true` → always emit `a-b` even for one line |
| `CASE_INSENSITIVE` | `true` | match quotes ignoring case |
| `HIGHLIGHTED_ONLY` | `false` | scan only highlighted markers (use the menu toggle, not this) |
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
