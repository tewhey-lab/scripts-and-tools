/**
 * QuoteFinder — refresh @AUTOLINE line references in a response document
 * =====================================================================
 *
 * WHAT IT DOES
 *   In the *response* doc you have markers like:
 *        @AUTOLINE123: "some quoted text from the manuscript"
 *        @AUTOLINE45-58: "a longer quote that spans several lines"
 *   For each marker this script:
 *     1. reads the quoted text (text inside the quotation marks),
 *     2. finds that exact text in the *primary* (manuscript) doc,
 *     3. works out the current line number(s) of that text,
 *     4. rewrites the marker  ->  @AUTOLINE<start>:  or  @AUTOLINE<start>-<end>:
 *     5. copies basic formatting (font color, bold, italic, super/subscript)
 *        from the manuscript text onto the quote in the response doc.
 *   (The in-document marker token stays @AUTOLINE; only the app is QuoteFinder.)
 *
 * TWO WAYS TO GET THE LINE NUMBER (you can try both and compare):
 *   • gDoc (Estimated Line #)  – computed in here by reproducing Google Docs'
 *                   line wrapping. Zero setup; exact for ordinary flowing text,
 *                   and it FLAGS anything near tables/images instead of guessing.
 *   • PDF/JSON (Exact Line #)  – reads Google's real rendered line numbers from a
 *                   JSON file produced by pdf_engine/extract_line_numbers.py.
 *                   Exact for any layout. See README for the 4-step recipe.
 *
 * SETUP (one time)
 *   1. Open your RESPONSE doc → Extensions → Apps Script.
 *   2. Paste this whole file in (replacing the default Code.gs). Save.
 *   3. Reload the response doc. A new menu "QuoteFinder" appears.
 *   4. QuoteFinder → "Settings…", paste the manuscript doc URL/ID (and, for the
 *      PDF/JSON method, the line-number file name or Drive ID). Save.
 *
 * NOTE: the first run will ask you to authorize the script (it needs to read
 *       the manuscript doc and edit this one). That's a standard Google prompt.
 */

// ----------------------------------------------------------------------------
// Configuration (defaults; the menu setters override these via doc properties)
// ----------------------------------------------------------------------------
var CONFIG = {
  PRIMARY_DOC_ID: '',                          // set via menu, or hardcode here
  EXACT_MAP_FILE_NAME: 'rendered_lines.json',  // JSON produced by the PDF engine
  EXACT_MAP_FILE_ID: '',                       // optional: pin an exact Drive file id
  // ---- marker format (define your convention here) ----
  HOOK: '@AUTOLINE',                           // the hook/prefix that starts a marker
  QUOTE_OPEN: '"',                             // mark that opens the quote (double quote; curly is accepted)
  QUOTE_CLOSE: '"',                            // mark that closes the quote
  OPEN_LOOKAHEAD: 20,                          // max chars between the hook and the OPENING quote mark
  MAX_LOOKAHEAD: 4000,                         // max chars of quote body (OPENING mark -> CLOSING mark)
  // -----------------------------------------------------
  COPY_FORMATTING: true,                       // copy color/bold/italic/super-sub
  ALWAYS_RANGE: false,                         // true => always write <hook>a-b
  CASE_INSENSITIVE: true,                      // match quotes ignoring case
  HIGHLIGHTED_ONLY: false,                     // scan only markers in highlighted text
  REVIEW_PANEL: 'dialog',                      // 'dialog' (wide, floating, non-blocking) or 'sidebar' (narrow, docked)
  REVIEW_WIDTH: 820,                           // review dialog width in px  (ignored for sidebar)
  REVIEW_HEIGHT: 780,                          // review dialog height in px (ignored for sidebar)
  LOWCONF_COLOR: '#fff2cc',                    // highlight: low-confidence number
  NOTFOUND_COLOR: '#f4cccc'                    // highlight: quote not found
};

// ----------------------------------------------------------------------------
// Menu
// ----------------------------------------------------------------------------
function onOpen() { buildMenu_(); }

// Builds (and rebuilds) the QuoteFinder menu. Called on open AND after toggling
// the scan scope, so the scope label reflects the current mode without a reload.
function buildMenu_() {
  var ui = DocumentApp.getUi();
  var hl = false;
  try { hl = (PropertiesService.getDocumentProperties().getProperty('HIGHLIGHTED_ONLY') === 'true'); } catch (e) {}
  var scopeLabel = hl
    ? '◉ Scope: HIGHLIGHTED text only — click to scan the whole document'
    : '○ Scope: WHOLE document — click to scan highlighted text only';
  var update = ui.createMenu('Update references (all at once)')
    .addItem('gDoc (Estimated Line #)', 'runSimulation')
    .addItem('PDF/JSON (Exact Line #)', 'runExact');
  var review = ui.createMenu('Review & update (step through imperfect)')
    .addItem('gDoc (Estimated Line #)', 'openReviewSim')
    .addItem('PDF/JSON (Exact Line #)', 'openReviewExact');
  ui.createMenu('QuoteFinder')
    .addSubMenu(update)
    .addSubMenu(review)
    .addSeparator()
    .addItem('Preview matches (dry run, no changes)', 'previewMatches')
    .addSeparator()
    .addItem(scopeLabel, 'toggleHighlightScope')
    .addItem('Settings…', 'showSettings')
    .addItem('Clear all QuoteFinder highlights', 'clearHighlights')
    .addToUi();
}

// Toggle the "scan highlighted text only" scope (persisted per document).
function toggleHighlightScope() {
  var p = PropertiesService.getDocumentProperties();
  var now = (p.getProperty('HIGHLIGHTED_ONLY') === 'true');
  p.setProperty('HIGHLIGHTED_ONLY', now ? 'false' : 'true');
  buildMenu_();   // refresh the menu label immediately — no document reload needed
  DocumentApp.getUi().alert('Scan scope is now: ' + (now ? 'THE WHOLE DOCUMENT' : 'HIGHLIGHTED TEXT ONLY') + '.');
}

function openReviewSim()   { openReview_('sim'); }
function openReviewExact() { openReview_('exact'); }
function openReview_(mode) {
  var c = cfg();
  if (!c.PRIMARY_DOC_ID) { DocumentApp.getUi().alert('Set the primary doc first (QuoteFinder → Settings…).'); return; }
  var title = 'QuoteFinder review — ' + (mode === 'sim' ? 'gDoc (Estimated Line #)' : 'PDF/JSON (Exact Line #)');
  var t = HtmlService.createTemplateFromFile('ReviewSidebar');
  t.mode = mode;
  var html = t.evaluate().setTitle(title);
  var ui = DocumentApp.getUi();
  // A native sidebar is locked to ~300px; a modeless dialog can be wide AND
  // still leaves the document editable (non-blocking) and draggable.
  if (c.REVIEW_PANEL === 'sidebar') ui.showSidebar(html);
  else ui.showModelessDialog(html.setWidth(c.REVIEW_WIDTH).setHeight(c.REVIEW_HEIGHT), title);
}

function cfg() {
  var p = PropertiesService.getDocumentProperties();
  var c = JSON.parse(JSON.stringify(CONFIG));
  c.PRIMARY_DOC_ID      = p.getProperty('PRIMARY_DOC_ID')      || c.PRIMARY_DOC_ID;
  c.EXACT_MAP_FILE_NAME = p.getProperty('EXACT_MAP_FILE_NAME') || c.EXACT_MAP_FILE_NAME;
  c.EXACT_MAP_FILE_ID   = p.getProperty('EXACT_MAP_FILE_ID')   || c.EXACT_MAP_FILE_ID;
  c.HIGHLIGHTED_ONLY    = (p.getProperty('HIGHLIGHTED_ONLY') === 'true');
  // marker format — editable via the "Settings…" dialog (falls back to CONFIG)
  var s;
  s = p.getProperty('HOOK');           if (s) c.HOOK = s;
  s = p.getProperty('QUOTE_OPEN');     if (s) c.QUOTE_OPEN = s;
  s = p.getProperty('QUOTE_CLOSE');    if (s) c.QUOTE_CLOSE = s;
  s = p.getProperty('OPEN_LOOKAHEAD'); if (s != null && s !== '' && !isNaN(parseInt(s, 10))) c.OPEN_LOOKAHEAD = parseInt(s, 10);
  s = p.getProperty('MAX_LOOKAHEAD');  if (s != null && s !== '' && !isNaN(parseInt(s, 10))) c.MAX_LOOKAHEAD = parseInt(s, 10);
  return c;
}

// Interactive settings (no code editing needed). Stored per document.
function showSettings() {
  var c = cfg();
  function a(s) { return String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;'); }
  var emf = c.EXACT_MAP_FILE_ID || c.EXACT_MAP_FILE_NAME || '';
  var h =
    '<style>body{font:13px Arial;margin:10px;color:#222}h3{font-size:13px;margin:14px 0 4px;color:#1a73e8;border-bottom:1px solid #e0e0e0;padding-bottom:2px}' +
    'label{display:block;margin:8px 0 2px;font-weight:bold}input{width:100%;padding:5px;box-sizing:border-box}' +
    '.hint{color:#777;font-weight:normal;font-size:12px}.row2{display:flex;gap:10px}.row2>div{flex:1}' +
    '.btns{margin-top:14px;display:flex;gap:8px}button{font:13px Arial;padding:7px 12px;border:1px solid #1a73e8;background:#1a73e8;color:#fff;border-radius:5px;cursor:pointer}' +
    'button.sec{background:#fff;color:#1a73e8}#msg{margin-top:10px;font-size:12px;min-height:14px}</style>' +
    '<h3>Documents</h3>' +
    '<label>Primary (manuscript) doc <span class="hint">— paste its URL or ID</span></label><input id="pid" value="' + a(c.PRIMARY_DOC_ID) + '">' +
    '<label>PDF/JSON line-number file <span class="hint">— Drive file name or ID (for the Exact engine)</span></label><input id="emf" value="' + a(emf) + '">' +
    '<h3>Marker format</h3>' +
    '<label>Hook (marker prefix) <span class="hint">— e.g. @AUTOLINE</span></label><input id="hook" value="' + a(c.HOOK) + '">' +
    '<div class="row2"><div><label>Opening quote mark</label><input id="qo" value="' + a(c.QUOTE_OPEN) + '"></div>' +
    '<div><label>Closing quote mark</label><input id="qc" value="' + a(c.QUOTE_CLOSE) + '"></div></div>' +
    '<div class="row2"><div><label>Look-ahead to opening quote <span class="hint">(chars after the hook)</span></label>' +
    '<input id="ola" type="number" min="0" value="' + a(c.OPEN_LOOKAHEAD) + '"></div>' +
    '<div><label>Max quote length <span class="hint">(opening &rarr; closing)</span></label>' +
    '<input id="mla" type="number" min="1" value="' + a(c.MAX_LOOKAHEAD) + '"></div></div>' +
    '<div id="msg"></div>' +
    '<div class="btns"><button onclick="save()">Save</button>' +
    '<button class="sec" onclick="resetFmt()">Reset marker format</button>' +
    '<button class="sec" onclick="google.script.host.close()">Cancel</button></div>' +
    '<script>' +
    'function v(id){return document.getElementById(id).value;}' +
    'function s(id,val){document.getElementById(id).value=val;}' +
    'function resetFmt(){s("hook","@AUTOLINE");s("qo",String.fromCharCode(34));s("qc",String.fromCharCode(34));s("ola","20");s("mla","4000");}' +
    'function save(){document.getElementById("msg").textContent="Saving…";' +
    'google.script.run.withSuccessHandler(function(){google.script.host.close();})' +
    '.withFailureHandler(function(e){document.getElementById("msg").textContent="Error: "+(e&&e.message?e.message:e);})' +
    '.saveSettings(JSON.stringify({PRIMARY_DOC_ID:v("pid"),EXACT_MAP:v("emf"),HOOK:v("hook"),QUOTE_OPEN:v("qo"),QUOTE_CLOSE:v("qc"),OPEN_LOOKAHEAD:v("ola"),MAX_LOOKAHEAD:v("mla")}));}' +
    '</script>';
  DocumentApp.getUi().showModalDialog(HtmlService.createHtmlOutput(h).setWidth(460).setHeight(560), 'QuoteFinder — settings');
}

// Persist settings; an empty field clears that override (reverts to default).
function saveSettings(json) {
  var o = JSON.parse(json || '{}');
  var p = PropertiesService.getDocumentProperties();
  function set(k, v) { v = (v == null) ? '' : String(v).trim(); if (v === '') p.deleteProperty(k); else p.setProperty(k, v); }
  // documents
  var pid = String(o.PRIMARY_DOC_ID == null ? '' : o.PRIMARY_DOC_ID).trim();
  if (!pid) p.deleteProperty('PRIMARY_DOC_ID'); else p.setProperty('PRIMARY_DOC_ID', extractDocId_(pid) || pid);
  var emf = String(o.EXACT_MAP == null ? '' : o.EXACT_MAP).trim();
  if (!emf) { p.deleteProperty('EXACT_MAP_FILE_ID'); p.deleteProperty('EXACT_MAP_FILE_NAME'); }
  else if (/^[A-Za-z0-9_\-]{20,}$/.test(emf)) { p.setProperty('EXACT_MAP_FILE_ID', emf); p.deleteProperty('EXACT_MAP_FILE_NAME'); }
  else { p.setProperty('EXACT_MAP_FILE_NAME', emf); p.deleteProperty('EXACT_MAP_FILE_ID'); }
  // marker format
  set('HOOK', o.HOOK);
  set('QUOTE_OPEN', o.QUOTE_OPEN);
  set('QUOTE_CLOSE', o.QUOTE_CLOSE);
  set('OPEN_LOOKAHEAD', o.OPEN_LOOKAHEAD);
  set('MAX_LOOKAHEAD', o.MAX_LOOKAHEAD);
  return 'ok';
}

// ----------------------------------------------------------------------------
// Scope: is a marker inside highlighted text? (background color present, and not
// one of the tool's own flag shades). Used when HIGHLIGHTED_ONLY is on.
// ----------------------------------------------------------------------------
function isHighlighted_(textEl, start, endIncl) {
  var lo = String(CONFIG.LOWCONF_COLOR).toLowerCase(), nf = String(CONFIG.NOTFOUND_COLOR).toLowerCase();
  for (var o = start; o <= endIncl; o++) {
    var bgc; try { bgc = textEl.getBackgroundColor(o); } catch (e) { bgc = null; }
    if (bgc) { var b = bgc.toLowerCase(); if (b !== lo && b !== nf) return true; }
  }
  return false;
}
function markerInScope_(c, textEl, mk) {
  if (!c.HIGHLIGHTED_ONLY) return true;
  var len = textEl.getText().length;
  var endIncl = Math.min(len - 1, mk.quoteRawEnd);
  if (endIncl < mk.matchStart) endIncl = mk.matchStart;
  return isHighlighted_(textEl, mk.matchStart, endIncl);
}

function extractDocId_(s) {
  var m = s.match(/\/d\/([A-Za-z0-9_\-]+)/);
  if (m) return m[1];
  if (/^[A-Za-z0-9_\-]{20,}$/.test(s)) return s;
  return '';
}

// ----------------------------------------------------------------------------
// Public entry points
// ----------------------------------------------------------------------------
function runSimulation()     { runUpdate_('sim'); }
function runExact()          { runUpdate_('exact'); }

// ============================================================================
// Text normalization (used for matching). Returns the normalized string and a
// map[normIndex] = rawIndex so we can translate matches back to real offsets.
// ============================================================================
function isBreakCode_(code) { return code === 10 || code === 13 || code === 11; } // \n \r \v
function isSpaceCode_(code) {  // spaces (NOT line breaks)
  return code === 32 || code === 9 || code === 160 || code === 0xFFFC ||
         code === 0x202F || code === 0x205F || code === 0x3000 ||
         (code >= 0x2000 && code <= 0x200B);
}
function isWsCode_(code) { return isBreakCode_(code) || isSpaceCode_(code); }

function normalize_(raw, caseInsensitive) {
  var norm = [], map = [];
  var prevSpace = true; // treat leading whitespace as collapsible
  for (var i = 0; i < raw.length; i++) {
    var code = raw.charCodeAt(i);
    if (isWsCode_(code)) {                       // any whitespace/break -> one space
      if (!prevSpace) { norm.push(' '); map.push(i); prevSpace = true; }
      continue;
    }
    var ch = raw.charAt(i);
    if (ch === '“' || ch === '”' || ch === '„' || ch === '«' || ch === '»') ch = '"';
    else if (ch === '‘' || ch === '’' || ch === '‚' || ch === '′') ch = "'";
    else if (ch === '–' || ch === '—' || ch === '−') ch = '-';
    norm.push(caseInsensitive ? ch.toLowerCase() : ch);
    map.push(i);
    prevSpace = false;
  }
  if (norm.length && norm[norm.length - 1] === ' ') { norm.pop(); map.pop(); } // trim trailing
  return { norm: norm.join(''), map: map };
}

// ============================================================================
// Font metrics (Adobe Core-14 widths, units/1000 em). Used by the SIMULATION
// engine to reproduce line wrapping. Google's Arial/Times metrics are very
// close to Helvetica/Times-Roman but not identical — hence "best effort".
// ============================================================================
var W_HELV = [278,278,355,556,556,889,667,191,333,333,389,584,278,333,278,278,556,556,556,556,556,556,556,556,556,556,278,278,584,584,584,556,1015,667,667,722,722,667,611,778,722,278,500,667,556,833,722,778,667,778,722,667,611,722,667,944,667,667,611,278,278,278,469,556,333,556,556,500,556,556,278,556,556,222,222,500,222,833,556,556,556,556,333,500,278,556,500,722,500,500,500,334,260,334,584];
var W_TIMES = [250,333,408,500,500,833,778,180,333,333,500,564,250,333,250,278,500,500,500,500,500,500,500,500,500,500,278,278,564,564,564,444,921,722,667,667,722,611,556,722,722,333,389,722,611,889,722,722,556,722,667,556,611,722,722,944,722,722,611,333,278,333,469,500,333,444,500,444,500,444,333,500,500,278,278,500,278,778,500,500,500,500,333,389,278,500,500,722,500,500,444,480,200,480,541];

function metricSet_(family) {
  family = (family || 'Arial').toLowerCase();
  if (/courier|mono|consol/.test(family)) return 'courier';
  if (/times|serif|georgia|cambria|garamond|minion|palatino|book antiqua/.test(family)) return 'times';
  return 'helv'; // arial, helvetica, calibri, verdana, roboto, etc. (approx)
}
function isExactFamily_(family) {
  family = (family || '').toLowerCase();
  return /arial|helvetica|times new roman|times|courier new|courier/.test(family);
}
function charWidth_(set, ch) {
  var code = ch.charCodeAt(0);
  if (set === 'courier') return 600;
  if (code >= 32 && code <= 126) return (set === 'times' ? W_TIMES : W_HELV)[code - 32];
  return set === 'times' ? 500 : 556; // non-ASCII fallback ~ average letter width
}

// ============================================================================
// Primary-document index: collect paragraph units, normalize each, and (for
// the simulation engine) lay each out to know which visual line every
// character falls on.
// ============================================================================
function collectUnits_(el, out, inTable) {
  var t = el.getType();
  if (t === DocumentApp.ElementType.PARAGRAPH || t === DocumentApp.ElementType.LIST_ITEM) {
    out.push({ el: el, inTable: !!inTable });
    return;
  }
  if (t === DocumentApp.ElementType.TABLE) {
    var tbl = el.asTable();
    for (var r = 0; r < tbl.getNumRows(); r++) {
      var row = tbl.getRow(r);
      for (var c = 0; c < row.getNumCells(); c++) {
        var cell = row.getCell(c);
        for (var k = 0; k < cell.getNumChildren(); k++) collectUnits_(cell.getChild(k), out, true);
      }
    }
    return;
  }
  if (el.getNumChildren) {
    for (var i = 0; i < el.getNumChildren(); i++) collectUnits_(el.getChild(i), out, inTable);
  }
}

function paragraphRuns_(textEl) {
  var runs = [];
  var len = textEl.getText().length;
  if (len === 0) return [{ start: 0, family: 'Arial', size: 11 }];
  var idx;
  try { idx = textEl.getTextAttributeIndices(); } catch (e) { idx = [0]; }
  if (!idx || idx.length === 0) idx = [0];
  for (var i = 0; i < idx.length; i++) {
    var o = idx[i], fam, sz;
    try { fam = textEl.getFontFamily(o); } catch (e) { fam = null; }
    try { sz = textEl.getFontSize(o); } catch (e) { sz = null; }
    runs.push({ start: o, family: fam || 'Arial', size: sz || 11 });
  }
  return runs;
}

function hasInlineImage_(para) {
  if (!para.getNumChildren) return false;
  for (var i = 0; i < para.getNumChildren(); i++) {
    if (para.getChild(i).getType() === DocumentApp.ElementType.INLINE_IMAGE) return true;
  }
  return false;
}

// Lay out one paragraph -> visual line starts (raw char offsets) + approx flag.
function layoutParagraph_(textEl, rawText, usableFirst, usableRest) {
  var lineStarts = [0];
  var approximated = false;
  if (rawText.length === 0) return { lineStarts: lineStarts, approximated: approximated };

  var runs = paragraphRuns_(textEl);
  function fontAt(i) {
    var lo = 0;
    for (var r = 0; r < runs.length; r++) { if (runs[r].start <= i) lo = r; else break; }
    return runs[lo];
  }
  function widthOf(i) {
    var f = fontAt(i);
    var set = metricSet_(f.family);
    if (!isExactFamily_(f.family)) approximated = true;
    return charWidth_(set, rawText.charAt(i)) * f.size / 1000.0;
  }

  var cur = 0;            // width used on current line (points)
  var lineIdx = 0;
  function usable() { return lineIdx === 0 ? usableFirst : usableRest; }

  var i = 0, n = rawText.length;
  while (i < n) {
    var code = rawText.charCodeAt(i);
    if (isBreakCode_(code)) {                 // forced line break
      lineIdx++; lineStarts.push(i + 1); cur = 0; i++; continue;
    }
    if (isSpaceCode_(code)) {                 // run of spaces: rides the current line
      var gw = 0, j = i;
      while (j < n && isSpaceCode_(rawText.charCodeAt(j))) { gw += widthOf(j); j++; }
      if (cur > 0) cur += gw;
      i = j; continue;
    }
    // a word = run of non-space, non-break chars
    var ws = i, ww = 0;
    while (i < n) {
      var cc = rawText.charCodeAt(i);
      if (isBreakCode_(cc) || isSpaceCode_(cc)) break;
      ww += widthOf(i); i++;
    }
    if (cur === 0) {
      cur = ww;                               // first word always placed
      if (ww > usable() + 0.01) {             // word longer than the line (approx break)
        approximated = true;
        var extra = Math.floor(ww / usable());
        for (var e = 0; e < extra; e++) { lineIdx++; lineStarts.push(ws); }
        cur = ww - extra * usable();
      }
    } else if (cur + ww <= usable() + 0.01) {
      cur += ww;
    } else {
      lineIdx++; lineStarts.push(ws); cur = ww;
      if (ww > usable() + 0.01) {
        approximated = true;
        var ex2 = Math.floor(ww / usable());
        for (var e2 = 0; e2 < ex2; e2++) { lineIdx++; lineStarts.push(ws); }
        cur = ww - ex2 * usable();
      }
    }
  }
  return { lineStarts: lineStarts, approximated: approximated };
}

function numOr0_(thunk) { try { var v = thunk(); return v ? v : 0; } catch (e) { return 0; } }
function lineForRaw_(lineStarts, rawOff) {
  var L = 0;
  for (var i = 0; i < lineStarts.length; i++) { if (lineStarts[i] <= rawOff) L = i; else break; }
  return L;
}

function buildPrimaryIndex_(c, withSim) {
  var doc = DocumentApp.openById(c.PRIMARY_DOC_ID);
  var body = doc.getBody();
  var units = [];
  for (var i = 0; i < body.getNumChildren(); i++) collectUnits_(body.getChild(i), units, false);

  var pageW = body.getPageWidth(), mL = body.getMarginLeft(), mR = body.getMarginRight();

  var paras = [];        // {textEl, raw, norm, map, startLine, lineStarts, lowConf}
  var simConcat = [];
  var simLineNo = [];
  var lineCounter = 0;
  var unmodeledSeen = false;

  for (var u = 0; u < units.length; u++) {
    var para = units[u].el;
    var textEl = para.editAsText();
    var raw = para.getText();
    var nm = normalize_(raw, c.CASE_INSENSITIVE);
    var rec = { textEl: textEl, raw: raw, norm: nm.norm, map: nm.map };

    if (withSim) {
      var indentS = numOr0_(function () { return para.getIndentStart(); });
      var indentE = numOr0_(function () { return para.getIndentEnd(); });
      var indentF = numOr0_(function () { return para.getIndentFirstLine(); });
      var uFirst = Math.max(40, pageW - mL - mR - indentF - indentE);
      var uRest  = Math.max(40, pageW - mL - mR - indentS - indentE);
      var lay = layoutParagraph_(textEl, raw, uFirst, uRest);
      var lowConf = units[u].inTable || hasInlineImage_(para) || lay.approximated || unmodeledSeen;
      if (units[u].inTable || hasInlineImage_(para)) unmodeledSeen = true;
      rec.startLine = lineCounter + 1;
      rec.lineStarts = lay.lineStarts;
      rec.lowConf = lowConf;
      for (var k = 0; k < rec.norm.length; k++) {
        var rawOff = rec.map[k];
        var L = lineForRaw_(rec.lineStarts, rawOff);
        simConcat.push(rec.norm.charAt(k));
        simLineNo.push(rec.startLine + L);
      }
      simConcat.push('\n'); // paragraph separator, attributed to this para's last line
      simLineNo.push(rec.startLine + (rec.lineStarts.length - 1));
      lineCounter += rec.lineStarts.length;
    }
    rec.idx = paras.length;
    paras.push(rec);
  }

  return {
    paras: paras,
    sim: withSim ? { concat: simConcat.join(''), lineNo: simLineNo, totalLines: lineCounter, unmodeledSeen: unmodeledSeen } : null
  };
}

// ============================================================================
// Exact engine: load rendered_lines.json from Drive -> concat + per-char line#
// ============================================================================
function loadExactMap_(c) {
  var file = null;
  if (c.EXACT_MAP_FILE_ID) {
    file = DriveApp.getFileById(c.EXACT_MAP_FILE_ID);
  } else {
    var it = DriveApp.getFilesByName(c.EXACT_MAP_FILE_NAME);
    var newest = null;
    while (it.hasNext()) { var f = it.next(); if (!newest || f.getLastUpdated() > newest.getLastUpdated()) newest = f; }
    file = newest;
  }
  if (!file) throw new Error('Could not find the PDF/JSON line-number file. Run extract_line_numbers.py, ' +
    'upload its JSON to your Drive, then set it in QuoteFinder → Settings….');
  var data = JSON.parse(file.getBlob().getDataAsString());
  if (!data || !data.length) throw new Error('The exact line-number file is empty or malformed.');
  return data; // [{n: <lineNumber>, text: "<line text>"}, ...] in document order
}

function buildExactIndex_(c, renderedLines) {
  // Join rendered lines with a SPACE (not a break): text wraps across rendered
  // lines at spaces, so a quote that spans several lines must read as one flow.
  // Each glue space is attributed to the line it follows.
  var concat = [], lineNo = [];
  for (var i = 0; i < renderedLines.length; i++) {
    var n = renderedLines[i].n;
    var nm = normalize_(String(renderedLines[i].text || ''), c.CASE_INSENSITIVE);
    for (var k = 0; k < nm.norm.length; k++) { concat.push(nm.norm.charAt(k)); lineNo.push(n); }
    concat.push(' '); lineNo.push(n);
  }
  return { concat: concat.join(''), lineNo: lineNo };
}

// ============================================================================
// Locate a quote -> {found, startLine, endLine, count}
// ============================================================================
function locate_(index, qnorm) {
  if (!qnorm) return { found: false };
  var first = index.concat.indexOf(qnorm);
  if (first < 0) return { found: false };
  var count = 0, from = 0, p;
  while ((p = index.concat.indexOf(qnorm, from)) >= 0) { count++; from = p + 1; }
  return { found: true, startLine: index.lineNo[first], endLine: index.lineNo[first + qnorm.length - 1], count: count };
}

// A quote containing an elision such as  [...]  /  […]  /  [. . .]  is a SPLIT
// quote: each side is matched independently and the reported range runs from the
// start of the first piece to the end of the last piece.
var ELLIPSIS_RE = /\[\s*(?:\.\s*){2,}\]|\[\s*…\s*\]/g;
function hasElision_(s) { ELLIPSIS_RE.lastIndex = 0; return ELLIPSIS_RE.test(s); }
function splitQuote_(quoteRaw) {                 // -> [{raw, offset}] pieces with real content
  ELLIPSIS_RE.lastIndex = 0;
  var pieces = [], last = 0, m;
  while ((m = ELLIPSIS_RE.exec(quoteRaw)) !== null) {
    pieces.push({ raw: quoteRaw.slice(last, m.index), offset: last });
    last = m.index + m[0].length;
  }
  pieces.push({ raw: quoteRaw.slice(last), offset: last });
  return pieces.filter(function (p) { return /\S/.test(p.raw); });
}

// Locate a (possibly split) quote in a concat index -> first-piece start .. last-piece end.
function lookupLineRange_(index, quoteRaw, c) {
  var pieces = splitQuote_(quoteRaw);
  if (!pieces.length) return { found: false };
  var first = null, last = null, count = 0;
  for (var i = 0; i < pieces.length; i++) {
    var info = locate_(index, normalize_(pieces[i].raw, c.CASE_INSENSITIVE).norm);
    if (!info.found) return { found: false, failedPiece: i };
    if (i === 0) first = info;
    last = info;
    count += info.count;
  }
  return { found: true, startLine: first.startLine, endLine: last.endLine, count: count, pieces: pieces.length };
}

// Copy formatting piece-by-piece (each piece from its own manuscript match).
function copyPartsFormatting_(primary, textEl, quoteRawStart, quoteRaw, c) {
  var pieces = splitQuote_(quoteRaw), missing = 0;
  for (var i = 0; i < pieces.length; i++) {
    var pn = normalize_(pieces[i].raw, c.CASE_INSENSITIVE).norm;
    if (!pn) continue;
    var src = findFormatSource_(primary.paras, pn);
    if (src) {
      try { copyFormatting_(src.para.textEl, src.para.map, src.normPos, textEl, quoteRawStart + pieces[i].offset, pieces[i].raw, pn, c.CASE_INSENSITIVE); }
      catch (e) { missing++; }
    } else missing++;
  }
  return { missing: missing };
}

function findFormatSource_(paras, qnorm) {
  for (var i = 0; i < paras.length; i++) {
    var pos = paras[i].norm.indexOf(qnorm);
    if (pos >= 0) return { para: paras[i], normPos: pos };
  }
  return null;
}

// ============================================================================
// Response-document marker parsing.  Format (all parts configurable in CONFIG):
//   <HOOK><num>[-<num>] : <QUOTE_OPEN>quoted text<QUOTE_CLOSE>
//   default:  @AUTOLINE12: "the quote"   /   @AUTOLINE12-15: "the quote"
// Numbers may be any length; the quote scan is bounded by MAX_LOOKAHEAD chars.
// ============================================================================
function escapeRegex_(s) { return String(s).replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); }
function bracketRe_(s, isOpen) {
  if (s === '"') return isOpen ? '["“]' : '["”]';   // accept straight OR curly double quotes
  if (s === "'") return isOpen ? "['‘]" : "['’]";   // straight OR curly single quotes
  return escapeRegex_(s);                            // any other literal mark(s), e.g. « » [ ]
}
function openLookahead_(c) { var v = parseInt(c.OPEN_LOOKAHEAD, 10); return (c.OPEN_LOOKAHEAD != null && v >= 0) ? v : 20; }
// A marker is simply the HOOK followed (within OPEN_LOOKAHEAD chars) by a quote.
// Any old number/colon between the hook and the quote is ignored on input and
// rewritten on output — we don't parse it.
function markerRe_(c) {
  var hook = escapeRegex_(c.HOOK || '@AUTOLINE');
  var max = Math.max(0, parseInt(c.MAX_LOOKAHEAD, 10) || 4000);  // OPENING -> CLOSING (quote body)
  var ola = openLookahead_(c);                                   // hook -> OPENING mark
  var openRe  = bracketRe_(c.QUOTE_OPEN  || '"', true);
  var closeRe = bracketRe_(c.QUOTE_CLOSE || '"', false);
  return new RegExp(
    hook + '\\s*[\\s\\S]{0,' + ola + '}?(' + openRe + ')([\\s\\S]{0,' + max + '}?)(' + closeRe + ')',
    'g');
}

// Finds bare hook occurrences, so a hook with no nearby quote can be flagged.
function hookRe_(c) { return new RegExp(escapeRegex_(c.HOOK || '@AUTOLINE'), 'g'); }
function openTestRe_(c) { return new RegExp('^\\s*[\\s\\S]{0,' + openLookahead_(c) + '}?(?:' + bracketRe_(c.QUOTE_OPEN || '"', true) + ')'); }
function markerScanner_(c) { return { full: markerRe_(c), hook: hookRe_(c), openTest: openTestRe_(c), ola: openLookahead_(c) }; }

// scan = markerScanner_(cfg()) built once by the caller and reused per paragraph.
// Returns entries with kind 'ok' (complete marker) or 'bad' (hook present but the
// quote is missing / unterminated / empty — flagged, never auto-updated).
function parseMarkersInPara_(para, scan) {
  var raw = para.getText();
  var found = [], spans = [];

  scan.full.lastIndex = 0;                          // 1) complete markers: hook + quote within range
  var m;
  while ((m = scan.full.exec(raw)) !== null) {
    var whole = m[0], start = m.index;
    var openLen = m[1].length, quote = m[2], closeLen = m[3].length;
    var quoteRawStart = start + (whole.length - closeLen - quote.length); // quote sits at the match end
    var prefixLen = whole.length - closeLen - quote.length - openLen;     // hook (+ any old num/colon), up to open mark
    spans.push([start, start + whole.length]);
    var e = {
      matchStart: start, prefixLen: prefixLen,
      oldNum: whole.slice(0, prefixLen).replace(/\s+$/, ''),  // original prefix (hook + any old number/colon)
      quoteRaw: quote, quoteRawStart: quoteRawStart, quoteRawEnd: quoteRawStart + quote.length
    };
    if (!quote || !/\S/.test(quote)) { e.kind = 'bad'; e.reason = 'empty quote'; e.quoteRawEnd = start + prefixLen; }
    else e.kind = 'ok';
    found.push(e);
  }

  scan.hook.lastIndex = 0;                           // 2) hooks with no usable quote within range
  var p;
  while ((p = scan.hook.exec(raw)) !== null) {
    var inside = false;                              // ignore hooks that sit inside a valid quote
    for (var k = 0; k < spans.length; k++) { if (p.index >= spans[k][0] && p.index < spans[k][1]) { inside = true; break; } }
    if (inside) continue;
    var hk = p[0];
    found.push({
      kind: 'bad',
      matchStart: p.index, prefixLen: hk.length,
      oldNum: hk,
      quoteRaw: '', quoteRawStart: p.index + hk.length, quoteRawEnd: p.index + hk.length,
      reason: scan.openTest.test(raw.slice(p.index + hk.length)) ? 'opening quote but no closing quote' : ('no quotation within ' + scan.ola + ' characters')
    });
  }

  found.sort(function (a, b) { return a.matchStart - b.matchStart; }); // document order (stable ordinals)
  return found;
}

// ============================================================================
// Formatting copy (manuscript -> response quote)
// ============================================================================
function readStyle_(t, o) {
  return {
    bold: t.isBold(o), italic: t.isItalic(o),
    color: t.getForegroundColor(o) || '#000000',
    align: t.getTextAlignment(o) || DocumentApp.TextAlignment.NORMAL
  };
}
function eqStyle_(a, b) { return a && b && a.bold === b.bold && a.italic === b.italic && a.color === b.color && a.align === b.align; }
function writeStyle_(t, a, b, st) {
  if (b < a) return;
  t.setBold(a, b, st.bold);
  t.setItalic(a, b, st.italic);
  t.setForegroundColor(a, b, st.color);
  t.setTextAlignment(a, b, st.align);
}

// Copy formatting char-by-char from a manuscript range to the response quote,
// aligning by normalized index (handles whitespace/case differences).
function copyFormatting_(srcTextEl, srcMap, srcNormPos, dstTextEl, dstQuoteRawStart, dstQuoteRaw, qnorm, caseInsensitive) {
  var dn = normalize_(dstQuoteRaw, caseInsensitive);
  var qlen = Math.min(qnorm.length, dn.norm.length);
  var runStartDst = -1, runEndDst = -1, runStyle = null;
  for (var k = 0; k < qlen; k++) {
    var srcRaw = srcMap[srcNormPos + k];
    var dstRaw = dstQuoteRawStart + dn.map[k];
    var st = readStyle_(srcTextEl, srcRaw);
    if (runStyle && eqStyle_(runStyle, st) && dstRaw === runEndDst + 1) {
      runEndDst = dstRaw;
    } else {
      if (runStyle) writeStyle_(dstTextEl, runStartDst, runEndDst, runStyle);
      runStartDst = dstRaw; runEndDst = dstRaw; runStyle = st;
    }
  }
  if (runStyle) writeStyle_(dstTextEl, runStartDst, runEndDst, runStyle);
}

// Copy formatting from src[srcStart..srcEndIncl] onto dst starting at dstStart
// (identical length: used after inserting the manuscript text verbatim).
function copyRangeIdentity_(srcTextEl, srcStart, srcEndIncl, dstTextEl, dstStart) {
  var n = srcEndIncl - srcStart;
  if (n < 0) return;
  var runStart = 0, cur = readStyle_(srcTextEl, srcStart);
  for (var i = 1; i <= n; i++) {
    var st = readStyle_(srcTextEl, srcStart + i);
    if (!eqStyle_(st, cur)) { writeStyle_(dstTextEl, dstStart + runStart, dstStart + i - 1, cur); runStart = i; cur = st; }
  }
  writeStyle_(dstTextEl, dstStart + runStart, dstStart + n, cur);
}

// safe background-color set/clear (the offset overload rejects null directly)
function bg_(t, a, b, color) {
  try {
    if (color === null) {
      var attr = {};
      attr[DocumentApp.Attribute.BACKGROUND_COLOR] = null;
      t.setAttributes(a, b, attr);
    } else {
      t.setBackgroundColor(a, b, color);
    }
  } catch (e) {}
}

// ============================================================================
// Main update routine
// ============================================================================
function runUpdate_(mode) {
  var ui = DocumentApp.getUi();
  var c = cfg();
  if (!c.PRIMARY_DOC_ID) { ui.alert('Set the primary doc first (QuoteFinder → Settings…).'); return; }

  var primary, exactIndex = null, simIndex = null;
  try {
    primary = buildPrimaryIndex_(c, mode === 'sim');
    if (mode === 'sim') simIndex = primary.sim;
    if (mode === 'exact') exactIndex = buildExactIndex_(c, loadExactMap_(c));
  } catch (e) { ui.alert('Setup problem:\n' + e.message); return; }

  var responseBody = DocumentApp.getActiveDocument().getBody();
  var units = [];
  for (var i = 0; i < responseBody.getNumChildren(); i++) collectUnits_(responseBody.getChild(i), units, false);

  var stats = { updated: 0, notFound: 0, lowConf: 0, multi: 0, noFormat: 0, malformed: 0 };
  var problems = [];
  var re = markerScanner_(c);

  for (var u = 0; u < units.length; u++) {
    var para = units[u].el;
    var markers = parseMarkersInPara_(para, re);
    if (!markers.length) continue;
    markers.sort(function (a, b) { return b.matchStart - a.matchStart; }); // right-to-left
    var textEl = para.editAsText();

    for (var mi = 0; mi < markers.length; mi++) {
      var mk = markers[mi];
      if (!markerInScope_(c, textEl, mk)) continue;
      if (mk.kind === 'bad') {                 // hook present, but quote missing/unterminated/empty
        stats.malformed++;
        if (!c.HIGHLIGHTED_ONLY) bg_(textEl, mk.matchStart, mk.matchStart + mk.prefixLen - 1, c.NOTFOUND_COLOR);
        problems.push('MALFORMED (' + mk.reason + '): ' + mk.oldNum);
        continue;
      }
      // line number(s): split quotes ([...]) span first-piece-start .. last-piece-end
      var lineInfo = lookupLineRange_(mode === 'sim' ? simIndex : exactIndex, mk.quoteRaw, c);

      if (!lineInfo.found) {
        stats.notFound++;
        if (!c.HIGHLIGHTED_ONLY) bg_(textEl, mk.matchStart, mk.matchStart + mk.prefixLen - 1, c.NOTFOUND_COLOR);
        problems.push('NOT FOUND: "' + snippet_(mk.quoteRaw) + '"');
        continue;
      }
      if (lineInfo.count > 1) { stats.multi++; problems.push('MULTIPLE (' + lineInfo.count + ') matches, used first: "' + snippet_(mk.quoteRaw) + '"'); }

      if (c.COPY_FORMATTING) {
        var fmt = copyPartsFormatting_(primary, textEl, mk.quoteRawStart, mk.quoteRaw, c);
        if (fmt.missing) { stats.noFormat++; problems.push('formatting skipped (' + fmt.missing + ' piece(s) not found in one paragraph): "' + snippet_(mk.quoteRaw) + '"'); }
      }

      var newNum = (c.ALWAYS_RANGE || lineInfo.startLine !== lineInfo.endLine)
        ? (lineInfo.startLine + '-' + lineInfo.endLine)
        : ('' + lineInfo.startLine);
      var newPrefix = c.HOOK + newNum + ': ';

      textEl.deleteText(mk.matchStart, mk.matchStart + mk.prefixLen - 1);
      textEl.insertText(mk.matchStart, newPrefix);

      var fp0 = splitQuote_(mk.quoteRaw)[0];
      var low = (mode === 'sim') && fp0 && lineNumberLowConf_(primary, normalize_(fp0.raw, c.CASE_INSENSITIVE).norm);
      if (low) stats.lowConf++;
      // In highlighted-only mode, never touch backgrounds (preserves your highlights).
      if (!c.HIGHLIGHTED_ONLY) bg_(textEl, mk.matchStart, mk.matchStart + newPrefix.length - 1, low ? c.LOWCONF_COLOR : null);
      stats.updated++;
    }
  }

  var scopeNote = c.HIGHLIGHTED_ONLY ? ' [highlighted text only]' : '';
  var msg = 'QuoteFinder (' + (mode === 'sim' ? 'gDoc (Estimated Line #)' : 'PDF/JSON (Exact Line #)') + ')' + scopeNote + ' done.\n\n' +
    'Updated: ' + stats.updated + '\n' +
    (stats.lowConf ? 'Low-confidence' + (c.HIGHLIGHTED_ONLY ? '' : ' (yellow)') + ': ' + stats.lowConf + '\n' : '') +
    (stats.notFound ? 'Not found' + (c.HIGHLIGHTED_ONLY ? '' : ' (red)') + ': ' + stats.notFound + '\n' : '') +
    (stats.malformed ? 'Malformed marker' + (c.HIGHLIGHTED_ONLY ? '' : ' (red)') + ': ' + stats.malformed + '\n' : '') +
    (stats.multi ? 'Multiple matches: ' + stats.multi + '\n' : '') +
    (stats.noFormat ? 'Formatting skipped: ' + stats.noFormat + '\n' : '');
  if (c.HIGHLIGHTED_ONLY && stats.updated === 0 && stats.notFound === 0)
    msg += '\nNo highlighted @AUTOLINE markers were found. Give the marker/quote text a highlight color, or turn off QuoteFinder → "Scan highlighted text only".\n';
  if (problems.length) msg += '\nDetails:\n• ' + problems.slice(0, 25).join('\n• ') + (problems.length > 25 ? '\n…(' + (problems.length - 25) + ' more)' : '');
  ui.alert(msg);
}

function lineNumberLowConf_(primary, qn) {
  if (primary.sim && primary.sim.unmodeledSeen) return true;
  for (var i = 0; i < primary.paras.length; i++) {
    var p = primary.paras[i];
    if (p.norm.indexOf(qn) >= 0) return !!p.lowConf;
  }
  return false;
}

function snippet_(s) { s = s.replace(/\s+/g, ' ').trim(); return s.length > 60 ? s.slice(0, 60) + '…' : s; }

// ============================================================================
// Dry-run: preview the line number each marker would get (no document changes).
// Uses the PDF/JSON exact engine when a line-number file is loaded; otherwise
// falls back to the gDoc estimate. Mirrors what "Update references" would write.
// ============================================================================
function previewMatches() {
  var ui = DocumentApp.getUi();
  var c = cfg();
  if (!c.PRIMARY_DOC_ID) { ui.alert('Set the primary doc first (QuoteFinder → Settings…).'); return; }

  var index = null, label = '', note = '';
  try {                                       // prefer the exact (PDF/JSON) engine — ground truth
    index = buildExactIndex_(c, loadExactMap_(c));
    label = 'PDF/JSON (Exact Line #)';
  } catch (e) {                               // fall back to the gDoc estimate
    try { index = buildPrimaryIndex_(c, true).sim; }
    catch (e2) { ui.alert('Setup problem:\n' + e2.message); return; }
    label = 'gDoc (Estimated Line #)';
    note = 'Using estimated numbers — no usable PDF/JSON line-number file. ' +
           'Set one in QuoteFinder → Settings… for exact numbers.';
  }

  var rows = collectPreviewRows_(c, index);
  showPreviewHtml_(rows, label, note, c.HIGHLIGHTED_ONLY);
}

function collectPreviewRows_(c, index) {
  var body = DocumentApp.getActiveDocument().getBody();
  var units = [];
  for (var i = 0; i < body.getNumChildren(); i++) collectUnits_(body.getChild(i), units, false);
  var rows = [], re = markerScanner_(c);
  for (var u = 0; u < units.length; u++) {
    var textEl = units[u].el.editAsText();
    var markers = parseMarkersInPara_(units[u].el, re);
    for (var mi = 0; mi < markers.length; mi++) {
      var mk = markers[mi];
      if (!markerInScope_(c, textEl, mk)) continue;
      if (mk.kind === 'bad') { rows.push({ old: mk.oldNum, quote: '⚠ ' + mk.reason, lines: 'malformed', found: false }); continue; }
      var info = lookupLineRange_(index, mk.quoteRaw, c);   // split quotes -> first..last range
      rows.push({ old: mk.oldNum, quote: snippet_(mk.quoteRaw), lines: fmtLines_(info), found: !!info.found });
    }
  }
  return rows;
}

function fmtLines_(info) {
  if (!info) return '—';
  if (!info.found) return 'not found';
  return info.startLine === info.endLine ? ('' + info.startLine) : (info.startLine + '-' + info.endLine);
}

function showPreviewHtml_(rows, engineLabel, note, highlightedOnly) {
  var h = ['<style>body{font:13px Arial;margin:8px}table{border-collapse:collapse;width:100%}',
    'td,th{border:1px solid #ccc;padding:3px 6px;text-align:left;vertical-align:top}',
    'th{background:#f0f0f0}.no{background:#f4cccc}</style>'];
  h.push('<p><b>' + rows.length + ' marker(s)' + (highlightedOnly ? ' — highlighted only' : '') + '.</b> ' +
    'Line numbers from <b>' + esc_(engineLabel) + '</b>. Dry run — nothing was changed.</p>');
  if (note) h.push('<p style="color:#b06000">' + esc_(note) + '</p>');
  h.push('<table><tr><th>#</th><th>old</th><th>quote</th><th>new # (' + esc_(engineLabel) + ')</th></tr>');
  for (var i = 0; i < rows.length; i++) {
    var r = rows[i];
    h.push('<tr><td>' + (i + 1) + '</td><td>' + esc_(r.old) + '</td><td>' + esc_(r.quote) + '</td>' +
      '<td' + (r.found ? '' : ' class="no"') + '>' + esc_(r.lines) + '</td></tr>');
  }
  h.push('</table>');
  DocumentApp.getUi().showModalDialog(HtmlService.createHtmlOutput(h.join('')).setWidth(680).setHeight(560), 'QuoteFinder — preview (dry run)');
}
function esc_(s) { return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;'); }

// ============================================================================
// Clear highlights left by previous runs
// ============================================================================
function clearHighlights() {
  var body = DocumentApp.getActiveDocument().getBody();
  var units = [];
  for (var i = 0; i < body.getNumChildren(); i++) collectUnits_(body.getChild(i), units, false);
  var re = markerScanner_(cfg());
  for (var u = 0; u < units.length; u++) {
    var markers = parseMarkersInPara_(units[u].el, re);
    var t = units[u].el.editAsText();
    for (var mi = 0; mi < markers.length; mi++) {
      var mk = markers[mi];
      bg_(t, mk.matchStart, mk.matchStart + mk.prefixLen - 1, null);
    }
  }
  DocumentApp.getUi().alert('Highlights cleared.');
}

// ============================================================================
// ============================================================================
// REVIEW WORKFLOW  (step through imperfect matches in a sidebar)
// ============================================================================
// ============================================================================

// ---- fuzzy matching --------------------------------------------------------
function lcsLen_(a, b) {
  var m = a.length, n = b.length;
  if (!m || !n) return 0;
  var dp = []; for (var j = 0; j <= n; j++) dp.push(0);
  for (var i = 1; i <= m; i++) {
    var prev = 0;
    for (var k = 1; k <= n; k++) {
      var tmp = dp[k];
      dp[k] = (a[i - 1] === b[k - 1]) ? prev + 1 : Math.max(dp[k], dp[k - 1]);
      prev = tmp;
    }
  }
  return dp[n];
}
function words_(s) { return s.split(' ').filter(function (w) { return w.length > 0; }); }
function wordSim_(a, b) {
  var wa = words_(a), wb = words_(b);
  if (!wa.length && !wb.length) return 1;
  if (!wa.length || !wb.length) return 0;
  var l = lcsLen_(wa, wb);
  return (2 * l) / (wa.length + wb.length);
}

// Best approximate span of qnorm inside one paragraph's normalized text.
function fuzzyInParagraph_(pn, qnorm, ws) {
  var ix = pn.indexOf(qnorm);
  if (ix >= 0) return { found: true, S: ix, E: ix + qnorm.length, sim: 1 };
  var S = -1;
  for (var p = ws.length; p >= 1; p--) {                 // longest matching prefix (words)
    var pos = pn.indexOf(ws.slice(0, p).join(' '));
    if (pos >= 0) { S = pos; break; }
  }
  var E = -1;
  for (var s = ws.length; s >= 1; s--) {                 // longest matching suffix (words)
    var sub = ws.slice(ws.length - s).join(' ');
    var pos2 = (S >= 0) ? pn.indexOf(sub, S) : pn.indexOf(sub);
    if (pos2 >= 0) { E = pos2 + sub.length; break; }
  }
  if (S < 0 && E < 0) return { found: false };
  if (S < 0) S = Math.max(0, E - qnorm.length);
  if (E < 0) E = Math.min(pn.length, S + qnorm.length);
  if (E <= S) E = Math.min(pn.length, S + qnorm.length);
  if (E - S > qnorm.length * 2 + 10) E = Math.min(pn.length, S + qnorm.length); // guard runaway span
  return { found: true, S: S, E: E, sim: wordSim_(qnorm, pn.substring(S, E)) };
}

// Locate a quote in the manuscript: exact, else best fuzzy candidate.
// Returns the manuscript paragraph + raw offsets (for text/formatting) and the
// candidate's line numbers (from the chosen engine).
function locateCandidate_(primary, mode, exactIndex, qnorm) {
  if (!qnorm) return { status: 'notfound' };
  for (var i = 0; i < primary.paras.length; i++) {           // exact in a paragraph
    var pos = primary.paras[i].norm.indexOf(qnorm);
    if (pos >= 0) return finishCandidate_(primary, mode, exactIndex, primary.paras[i], pos, pos + qnorm.length, qnorm, 1, 'exact');
  }
  var ws = words_(qnorm);
  if (!ws.length) return { status: 'notfound' };
  var probe = ws.slice().sort(function (a, b) { return b.length - a.length; })[0]; // distinctive word
  var best = null;
  for (var j = 0; j < primary.paras.length; j++) {
    var pn = primary.paras[j].norm;
    if (probe.length >= 4 && pn.indexOf(probe) < 0) continue; // quick filter
    var f = fuzzyInParagraph_(pn, qnorm, ws);
    if (f.found && (!best || f.sim > best.sim)) best = { idx: j, S: f.S, E: f.E, sim: f.sim };
  }
  if (!best) {                                                // fallback: scan all
    for (var k = 0; k < primary.paras.length; k++) {
      var f2 = fuzzyInParagraph_(primary.paras[k].norm, qnorm, ws);
      if (f2.found && (!best || f2.sim > best.sim)) best = { idx: k, S: f2.S, E: f2.E, sim: f2.sim };
    }
  }
  if (!best || best.sim < 0.25) return { status: 'notfound' };
  return finishCandidate_(primary, mode, exactIndex, primary.paras[best.idx], best.S, best.E, qnorm, best.sim, 'fuzzy');
}

function finishCandidate_(primary, mode, exactIndex, para, S, E, qnorm, sim, status) {
  var candNorm = para.norm.substring(S, E);
  var lines = candidateLines_(mode, exactIndex, para, S, E, candNorm);
  return {
    status: status, similarity: sim, srcParaIndex: para.idx, srcNormPos: S,
    srcRawStart: para.map[S], srcRawEndIncl: para.map[E - 1],
    startLine: lines.start, endLine: lines.end, candNorm: candNorm
  };
}
function candidateLines_(mode, exactIndex, para, S, E, candNorm) {
  if (mode === 'sim') {
    return {
      start: para.startLine + lineForRaw_(para.lineStarts, para.map[S]),
      end:   para.startLine + lineForRaw_(para.lineStarts, para.map[E - 1])
    };
  }
  var loc = locate_(exactIndex, candNorm);             // exact engine: find the text in the PDF map
  return loc.found ? { start: loc.startLine, end: loc.endLine } : { start: null, end: null };
}

// Split-quote aware candidate: a quote with [...] is matched piece-by-piece; the
// reported range spans the first piece's start to the last piece's end.
function locateCandidate2_(primary, mode, exactIndex, quoteRaw, c) {
  if (!hasElision_(quoteRaw)) return locateCandidate_(primary, mode, exactIndex, normalize_(quoteRaw, c.CASE_INSENSITIVE).norm);
  var pieces = splitQuote_(quoteRaw);
  if (!pieces.length) return { status: 'notfound' };
  var cands = [];
  for (var i = 0; i < pieces.length; i++) {
    var cand = locateCandidate_(primary, mode, exactIndex, normalize_(pieces[i].raw, c.CASE_INSENSITIVE).norm);
    if (cand.status === 'notfound') return { status: 'notfound' };
    cands.push(cand);
  }
  var first = cands[0], last = cands[cands.length - 1], allExact = true, minSim = 1;
  for (var k = 0; k < cands.length; k++) {
    if (!(cands[k].status === 'exact' && cands[k].startLine != null)) allExact = false;
    if (cands[k].similarity < minSim) minSim = cands[k].similarity;
  }
  return {
    status: (allExact && first.startLine != null && last.endLine != null) ? 'exact' : 'fuzzy',
    similarity: minSim, split: true, pieces: cands,
    srcParaIndex: first.srcParaIndex, srcNormPos: first.srcNormPos,
    srcRawStart: first.srcRawStart, srcRawEndIncl: first.srcRawEndIncl,
    startLine: first.startLine, endLine: last.endLine
  };
}
function candNewHtml_(primary, cand) {
  if (cand.split) return cand.pieces.map(function (pc) { return rangeToHtml_(primary.paras[pc.srcParaIndex].textEl, pc.srcRawStart, pc.srcRawEndIncl); }).join(' <span class="muted">[…]</span> ');
  return rangeToHtml_(primary.paras[cand.srcParaIndex].textEl, cand.srcRawStart, cand.srcRawEndIncl);
}
function candNewText_(primary, cand) {
  if (cand.split) return cand.pieces.map(function (pc) { return primary.paras[pc.srcParaIndex].raw.substring(pc.srcRawStart, pc.srcRawEndIncl + 1); }).join(' [...] ');
  return primary.paras[cand.srcParaIndex].raw.substring(cand.srcRawStart, cand.srcRawEndIncl + 1);
}

// ---- formatted-HTML + diff rendering --------------------------------------
function styleHtmlAt_(t, o) {
  var a = t.getTextAlignment(o);
  var al = (a === DocumentApp.TextAlignment.SUPERSCRIPT) ? 'SUP'
         : (a === DocumentApp.TextAlignment.SUBSCRIPT) ? 'SUB' : 'NORM';
  return { bold: !!t.isBold(o), italic: !!t.isItalic(o), color: t.getForegroundColor(o) || '#000000', align: al };
}
function eqHtmlStyle_(a, b) { return a.bold === b.bold && a.italic === b.italic && a.color === b.color && a.align === b.align; }
function rangeToHtml_(textEl, start, endIncl) {
  if (endIncl < start) return '<em>(empty)</em>';
  var T = textEl.getText();
  var html = [];
  function emit(a, b, st) {
    var open = [], close = [];
    if (st.align === 'SUP') { open.push('<sup>'); close.unshift('</sup>'); }
    else if (st.align === 'SUB') { open.push('<sub>'); close.unshift('</sub>'); }
    if (st.bold) { open.push('<b>'); close.unshift('</b>'); }
    if (st.italic) { open.push('<i>'); close.unshift('</i>'); }
    if (st.color && st.color !== '#000000') { open.push('<span style="color:' + st.color + '">'); close.unshift('</span>'); }
    html.push(open.join('') + esc_(T.substring(a, b + 1)) + close.join(''));
  }
  var runStart = start, cur = styleHtmlAt_(textEl, start);
  for (var o = start + 1; o <= endIncl; o++) {
    var st = styleHtmlAt_(textEl, o);
    if (!eqHtmlStyle_(st, cur)) { emit(runStart, o - 1, cur); runStart = o; cur = st; }
  }
  emit(runStart, endIncl, cur);
  return html.join('');
}
function collapseWs_(s) { return String(s).replace(/\s+/g, ' ').trim(); }
function wordDiffHtml_(oldS, newS) {                    // word-level diff, manuscript vs yours
  var A = collapseWs_(oldS).split(' '), B = collapseWs_(newS).split(' ');
  var m = A.length, n = B.length;
  var dp = []; for (var i = 0; i <= m; i++) { dp.push([]); for (var j = 0; j <= n; j++) dp[i].push(0); }
  for (i = 1; i <= m; i++) for (j = 1; j <= n; j++)
    dp[i][j] = (A[i - 1] === B[j - 1]) ? dp[i - 1][j - 1] + 1 : Math.max(dp[i - 1][j], dp[i][j - 1]);
  var out = []; i = m; j = n;
  while (i > 0 && j > 0) {
    if (A[i - 1] === B[j - 1]) { out.unshift(esc_(A[i - 1])); i--; j--; }
    else if (dp[i - 1][j] >= dp[i][j - 1]) { out.unshift('<del>' + esc_(A[i - 1]) + '</del>'); i--; }
    else { out.unshift('<ins>' + esc_(B[j - 1]) + '</ins>'); j--; }
  }
  while (i > 0) { out.unshift('<del>' + esc_(A[i - 1]) + '</del>'); i--; }
  while (j > 0) { out.unshift('<ins>' + esc_(B[j - 1]) + '</ins>'); j--; }
  return out.join(' ');
}

// ---- server entry points called from the sidebar --------------------------
function buildReviewData(mode) {
  var c = cfg();
  var primary = buildPrimaryIndex_(c, mode === 'sim');
  var exactIndex = (mode === 'exact') ? buildExactIndex_(c, loadExactMap_(c)) : null;

  var body = DocumentApp.getActiveDocument().getBody();
  var units = [];
  for (var i = 0; i < body.getNumChildren(); i++) collectUnits_(body.getChild(i), units, false);

  var ordinal = 0, total = 0, exactCount = 0, malformedCount = 0, items = [];
  var re = markerScanner_(c);
  for (var u = 0; u < units.length; u++) {
    var para = units[u].el, textEl = para.editAsText();
    var markers = parseMarkersInPara_(para, re);
    for (var mi = 0; mi < markers.length; mi++) {
      var mk = markers[mi];
      if (!markerInScope_(c, textEl, mk)) { ordinal++; continue; } // out of scope: keep ordinals aligned
      total++;
      if (mk.kind === 'bad') {                       // hook present, quote missing/unterminated/empty
        malformedCount++;
        items.push({
          ordinal: ordinal, oldNum: mk.oldNum, status: 'malformed', reason: mk.reason,
          similarity: 0, startLine: null, endLine: null,
          oldHtml: rangeToHtml_(textEl, mk.matchStart, mk.matchStart + mk.prefixLen - 1), newHtml: '', diffHtml: ''
        });
        ordinal++; continue;
      }
      var cand = locateCandidate2_(primary, mode, exactIndex, mk.quoteRaw, c);  // split-quote aware
      var isExact = (cand.status === 'exact' && cand.startLine != null);
      if (isExact) exactCount++;

      var hasCand = (cand.status !== 'notfound');
      items.push({                                  // EVERY in-scope marker (client filters the view)
        ordinal: ordinal,
        oldNum: mk.oldNum,
        status: isExact ? 'exact' : cand.status,    // 'exact' | 'fuzzy' | 'notfound'
        similarity: Math.round((cand.similarity || 0) * 100),
        startLine: hasCand ? cand.startLine : null,
        endLine: hasCand ? cand.endLine : null,
        split: !!cand.split,
        oldHtml: rangeToHtml_(textEl, mk.quoteRawStart, mk.quoteRawEnd - 1),
        newHtml: hasCand ? candNewHtml_(primary, cand) : '',
        diffHtml: hasCand ? wordDiffHtml_(mk.quoteRaw, candNewText_(primary, cand)) : ''
      });
      ordinal++;
    }
  }
  return { mode: mode, total: total, exactCount: exactCount, imperfectCount: total - exactCount,
           malformedCount: malformedCount, items: items, scope: c.HIGHLIGHTED_ONLY ? 'highlighted' : 'all' };
}

// Apply one pass; edits within each paragraph go right-to-left so offsets stay
// valid. decisionsJson = JSON map { ordinal: 'use' | 'skip' } for reviewed
// markers. Default when an ordinal is absent: refresh exact matches, skip the
// rest (so "only changes" view still auto-applies the perfect ones).
function applyReviewDecisions(mode, decisionsJson) {
  var c = cfg();
  var dec = {};
  try { dec = JSON.parse(decisionsJson || '{}') || {}; } catch (e) { dec = {}; }

  var primary = buildPrimaryIndex_(c, mode === 'sim');
  var exactIndex = (mode === 'exact') ? buildExactIndex_(c, loadExactMap_(c)) : null;

  var body = DocumentApp.getActiveDocument().getBody();
  var units = [];
  for (var i = 0; i < body.getNumChildren(); i++) collectUnits_(body.getChild(i), units, false);

  var n = { refreshed: 0, replaced: 0, skipped: 0, flagged: 0 };
  var ordinal = 0, re = markerScanner_(c);
  for (var u = 0; u < units.length; u++) {
    var para = units[u].el, textEl = para.editAsText();
    var markers = parseMarkersInPara_(para, re);
    markers.forEach(function (mk) { mk.__ord = ordinal++; });          // ordinals in document order
    markers.slice().sort(function (a, b) { return b.matchStart - a.matchStart; }) // edit right-to-left
      .forEach(function (mk) {
        if (!markerInScope_(c, textEl, mk)) return;                    // leave out-of-scope untouched
        if (mk.kind === 'bad') {                                       // flag malformed markers red, never edit
          if (!c.HIGHLIGHTED_ONLY) bg_(textEl, mk.matchStart, mk.matchStart + mk.prefixLen - 1, c.NOTFOUND_COLOR);
          n.flagged++; return;
        }
        var cand = locateCandidate2_(primary, mode, exactIndex, mk.quoteRaw, c);  // split-quote aware
        var isExact = (cand.status === 'exact' && cand.startLine != null);
        var d = dec[mk.__ord];                                         // 'use' | 'skip' | undefined
        if (d === 'skip') { n.skipped++; return; }
        var useIt = (d === 'use') || (d === undefined && isExact);     // default: auto-apply exact only
        if (!useIt) { n.skipped++; return; }
        if (isExact) { applyMarker_(textEl, mk, primary, cand, c, false); n.refreshed++; }      // number + formatting
        else if (cand.status === 'fuzzy' && cand.startLine != null) { applyMarker_(textEl, mk, primary, cand, c, true); n.replaced++; } // text + formatting + number
        else { n.skipped++; }                                          // 'use' but no usable line #
      });
  }
  return 'Applied: ' + n.refreshed + ' refreshed, ' + n.replaced + ' replaced; ' + n.skipped + ' skipped' +
    (n.flagged ? ', ' + n.flagged + ' flagged (malformed)' : '') + '.';
}

// replaceText=true -> swap the quote for the manuscript text + its formatting.
// Always refresh the @AUTOLINE number. Order: quote edits first, prefix last.
function applyMarker_(textEl, mk, primary, cand, c, replaceText) {
  if (cand.split || hasElision_(mk.quoteRaw)) {        // split quote: keep the user's elided text, copy formatting per piece
    if (c.COPY_FORMATTING) copyPartsFormatting_(primary, textEl, mk.quoteRawStart, mk.quoteRaw, c);
  } else if (replaceText) {
    var src = primary.paras[cand.srcParaIndex];
    var newTextRaw = src.raw.substring(cand.srcRawStart, cand.srcRawEndIncl + 1);
    textEl.deleteText(mk.quoteRawStart, mk.quoteRawEnd - 1);
    textEl.insertText(mk.quoteRawStart, newTextRaw);
    if (c.COPY_FORMATTING) copyRangeIdentity_(src.textEl, cand.srcRawStart, cand.srcRawEndIncl, textEl, mk.quoteRawStart);
  } else if (c.COPY_FORMATTING) {
    var src2 = primary.paras[cand.srcParaIndex];
    try {
      copyFormatting_(src2.textEl, src2.map, cand.srcNormPos, textEl, mk.quoteRawStart, mk.quoteRaw,
        normalize_(mk.quoteRaw, c.CASE_INSENSITIVE).norm, c.CASE_INSENSITIVE);
    } catch (e) {}
  }
  var numStr = (c.ALWAYS_RANGE || cand.startLine !== cand.endLine) ? (cand.startLine + '-' + cand.endLine) : ('' + cand.startLine);
  var newPrefix = c.HOOK + numStr + ': ';
  textEl.deleteText(mk.matchStart, mk.matchStart + mk.prefixLen - 1);
  textEl.insertText(mk.matchStart, newPrefix);
}
