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
  PDF_WINDOW: 40,                              // exact engine: rendered lines to search each side of the gDoc estimate
  PDF_FUZZY_MIN: 0.5,                          // exact engine: min word-similarity to accept a fuzzy PDF refinement
  REVIEW_PANEL: 'dialog',                      // 'dialog' (wide, floating, non-blocking) or 'sidebar' (narrow, docked)
  REVIEW_WIDTH: 820,                           // review dialog width in px  (ignored for sidebar)
  REVIEW_HEIGHT: 780,                          // review dialog height in px (ignored for sidebar)
  LOWCONF_COLOR: '#fff2cc',                    // highlight: low-confidence number
  NOTFOUND_COLOR: '#f4cccc',                   // highlight: quote not found
  // Marker scaffold style — the hook + line number(s) + colon + quotation marks (NOT the
  // quoted text). Choose how it's styled when a marker is updated (set in Settings):
  //   'insert' = keep the formatting where it sits (no restyling — the original behavior)
  //   'named'  = the document-wide "Normal text" style (needs the Docs API advanced service)
  //   'set'    = the explicit font/size/color/weight below (the default)
  MARKER_STYLE_MODE: 'set',                    // 'insert' | 'named' | 'set'
  MARKER_FONT_FAMILY: 'Arial',                 // 'set' mode font (blank => Arial)
  MARKER_FONT_SIZE: 11,                        // 'set' mode size (0/blank => 11)
  MARKER_FONT_COLOR: '#000000',                // 'set' mode color (blank => #000000)
  MARKER_BOLD: false,                          // 'set' mode: bold
  MARKER_ITALIC: false                         // 'set' mode: italic
};

// ----------------------------------------------------------------------------
// Menu
// ----------------------------------------------------------------------------
function onOpen() { buildMenu_(); }

// Builds the QuoteFinder menu. To scope to a subset, drag-select the quote line(s)
// and run a command — selection is the scope (no menu toggle needed).
function buildMenu_() {
  var ui = DocumentApp.getUi();
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
    .addItem('Settings…', 'showSettings')
    .addItem('Clear all QuoteFinder flags', 'clearHighlights')
    .addToUi();
}

function openReviewSim()   { openReview_('sim'); }
function openReviewExact() { openReview_('exact'); }
function openReview_(mode) {
  var c = cfg();
  if (!c.PRIMARY_DOC_ID) { DocumentApp.getUi().alert('Set the primary doc first (QuoteFinder → Settings…).'); return; }
  // This is the last code that runs with the live selection still available, so snapshot
  // it here for the whole panel session (the modeless dialog can't read getSelection()).
  // No selection -> writes null, which CLEARS any stale snapshot.
  writeSelectionSnapshot_(selectionSigSet_(c));
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

// Normalize a color to '#rrggbb' (strip '#', lowercase, expand 3-digit hex); '' if not a valid hex.
function normHexColor_(s) {
  s = String(s == null ? '' : s).trim().toLowerCase().replace(/^#/, '');
  if (/^[0-9a-f]{3}$/.test(s)) s = s.charAt(0) + s.charAt(0) + s.charAt(1) + s.charAt(1) + s.charAt(2) + s.charAt(2);
  return /^[0-9a-f]{6}$/.test(s) ? ('#' + s) : '';
}

function cfg() {
  var p = PropertiesService.getDocumentProperties();
  var c = JSON.parse(JSON.stringify(CONFIG));
  c.PRIMARY_DOC_ID      = p.getProperty('PRIMARY_DOC_ID')      || c.PRIMARY_DOC_ID;
  c.EXACT_MAP_FILE_NAME = p.getProperty('EXACT_MAP_FILE_NAME') || c.EXACT_MAP_FILE_NAME;
  c.EXACT_MAP_FILE_ID   = p.getProperty('EXACT_MAP_FILE_ID')   || c.EXACT_MAP_FILE_ID;
  // marker format — editable via the "Settings…" dialog (falls back to CONFIG)
  var s;
  s = p.getProperty('HOOK');           if (s) c.HOOK = s;
  s = p.getProperty('QUOTE_OPEN');     if (s) c.QUOTE_OPEN = s;
  s = p.getProperty('QUOTE_CLOSE');    if (s) c.QUOTE_CLOSE = s;
  s = p.getProperty('OPEN_LOOKAHEAD'); if (s != null && s !== '' && !isNaN(parseInt(s, 10))) c.OPEN_LOOKAHEAD = parseInt(s, 10);
  s = p.getProperty('MAX_LOOKAHEAD');  if (s != null && s !== '' && !isNaN(parseInt(s, 10))) c.MAX_LOOKAHEAD = parseInt(s, 10);
  // marker scaffold style
  s = p.getProperty('MARKER_STYLE_MODE'); if (s === 'insert' || s === 'named' || s === 'set') c.MARKER_STYLE_MODE = s;
  s = p.getProperty('MARKER_FONT_FAMILY'); if (s != null) c.MARKER_FONT_FAMILY = s;
  s = p.getProperty('MARKER_FONT_SIZE');  if (s != null && s !== '' && !isNaN(parseFloat(s))) c.MARKER_FONT_SIZE = parseFloat(s);
  s = p.getProperty('MARKER_FONT_COLOR'); if (s != null) { var nc = normHexColor_(s); if (nc) c.MARKER_FONT_COLOR = nc; }   // ignore a legacy/non-hex value -> keep the default
  s = p.getProperty('MARKER_BOLD');   if (s != null && s !== '') c.MARKER_BOLD = (s === 'true');
  s = p.getProperty('MARKER_ITALIC'); if (s != null && s !== '') c.MARKER_ITALIC = (s === 'true');
  return c;
}

// Interactive settings (no code editing needed). Stored per document.
function showSettings() {
  var c = cfg();
  function a(s) { return String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;'); }
  var emf = c.EXACT_MAP_FILE_ID || c.EXACT_MAP_FILE_NAME || '';
  // For the "Use the document's Normal text style" option: what the Docs API reports, and whether
  // the advanced service is even enabled (so we can tell the user if that option won't work yet).
  var named = namedNormalStyle_();
  var namedDesc;
  if (!docsApiAvailable_()) {
    namedDesc = '<span style="color:#b00">Google Docs API not enabled</span> — turn it on (Install step 7) or this option won’t restyle.';
  } else if (named.family === '' && named.size === '' && named.color === '') {
    namedDesc = 'Docs API on, but the Normal text style couldn’t be read.';
  } else {
    var nd = [];
    if (named.family) nd.push('<b>' + a(named.family) + '</b>');
    if (named.size !== '' && named.size != null) nd.push('<b>' + a(named.size) + '&nbsp;pt</b>');
    if (named.color) nd.push((/^#[0-9a-f]{3}([0-9a-f]{3})?$/i.test(named.color) ? '<span style="display:inline-block;width:10px;height:10px;border:1px solid #bbb;vertical-align:middle;margin-right:3px;background:' + a(named.color) + '"></span>' : '') + '<b>' + a(named.color) + '</b>');
    namedDesc = 'currently ' + nd.join(', ');
  }
  var mode = c.MARKER_STYLE_MODE || 'set';
  function ck(m) { return mode === m ? ' checked' : ''; }
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
    '<h3>Marker text style <span class="hint">— the hook, line #, colon &amp; quote marks (not the quoted text)</span></h3>' +
    '<label style="font-weight:normal"><input type="radio" name="mode" value="insert" style="width:auto;margin-right:6px"' + ck('insert') + ' onchange="syncMode()">' +
    'Keep the formatting where it&rsquo;s inserted <span class="hint">(no restyling)</span></label>' +
    '<label style="font-weight:normal"><input type="radio" name="mode" value="named" style="width:auto;margin-right:6px"' + ck('named') + ' onchange="syncMode()">' +
    'Use the document&rsquo;s <b>Normal text</b> style <span class="hint">— ' + namedDesc + '</span></label>' +
    '<label style="font-weight:normal"><input type="radio" name="mode" value="set" style="width:auto;margin-right:6px"' + ck('set') + ' onchange="syncMode()">' +
    'Use a set style:</label>' +
    '<div id="setbox" style="padding-left:22px">' +
    '<div class="row2"><div><label>Font</label><input id="mff" value="' + a(c.MARKER_FONT_FAMILY) + '" placeholder="Arial"></div>' +
    '<div><label>Size</label><input id="mfs" type="number" min="1" value="' + (c.MARKER_FONT_SIZE ? a(c.MARKER_FONT_SIZE) : '') + '" placeholder="11"></div>' +
    '<div><label>Color</label><input id="mfc" value="' + a(c.MARKER_FONT_COLOR) + '" placeholder="#000000"></div></div>' +
    '<label style="font-weight:normal;margin-top:6px"><input type="checkbox" id="mb" style="width:auto;margin-right:6px"' + (c.MARKER_BOLD ? ' checked' : '') + '>Bold' +
    '<input type="checkbox" id="mi" style="width:auto;margin:0 6px 0 16px"' + (c.MARKER_ITALIC ? ' checked' : '') + '>Italic</label>' +
    '<div class="hint" style="margin-top:4px"><a href="#" onclick="resetStyle();return false;">Reset to Arial, 11, #000000, regular</a></div>' +
    '</div>' +
    '<div id="msg"></div>' +
    '<div class="btns"><button onclick="save()">Save</button>' +
    '<button class="sec" onclick="resetFmt()">Reset marker format</button>' +
    '<button class="sec" onclick="google.script.host.close()">Cancel</button></div>' +
    '<script>' +
    'function v(id){return document.getElementById(id).value;}' +
    'function s(id,val){document.getElementById(id).value=val;}' +
    'function mode(){var r=document.getElementsByName("mode");for(var i=0;i<r.length;i++)if(r[i].checked)return r[i].value;return "set";}' +
    'function syncMode(){document.getElementById("setbox").style.opacity=(mode()==="set")?"1":"0.45";}' +
    'function resetFmt(){s("hook","@AUTOLINE");s("qo",String.fromCharCode(34));s("qc",String.fromCharCode(34));s("ola","20");s("mla","4000");}' +
    'function resetStyle(){s("mff","Arial");s("mfs","11");s("mfc","#000000");document.getElementById("mb").checked=false;document.getElementById("mi").checked=false;}' +
    'function save(){document.getElementById("msg").textContent="Saving…";' +
    'google.script.run.withSuccessHandler(function(){google.script.host.close();})' +
    '.withFailureHandler(function(e){document.getElementById("msg").textContent="Error: "+(e&&e.message?e.message:e);})' +
    '.saveSettings(JSON.stringify({PRIMARY_DOC_ID:v("pid"),EXACT_MAP:v("emf"),HOOK:v("hook"),QUOTE_OPEN:v("qo"),QUOTE_CLOSE:v("qc"),OPEN_LOOKAHEAD:v("ola"),MAX_LOOKAHEAD:v("mla"),MARKER_STYLE_MODE:mode(),MARKER_FONT_FAMILY:v("mff"),MARKER_FONT_SIZE:v("mfs"),MARKER_FONT_COLOR:v("mfc"),MARKER_BOLD:document.getElementById("mb").checked,MARKER_ITALIC:document.getElementById("mi").checked}));}' +
    'syncMode();' +
    '</script>';
  DocumentApp.getUi().showModalDialog(HtmlService.createHtmlOutput(h).setWidth(460).setHeight(680), 'QuoteFinder — settings');
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
  // marker scaffold style
  var mode = (o.MARKER_STYLE_MODE === 'insert' || o.MARKER_STYLE_MODE === 'named' || o.MARKER_STYLE_MODE === 'set') ? o.MARKER_STYLE_MODE : 'set';
  p.setProperty('MARKER_STYLE_MODE', mode);
  set('MARKER_FONT_FAMILY', o.MARKER_FONT_FAMILY);      // blank => 'set' mode uses Arial
  var mfs = String(o.MARKER_FONT_SIZE == null ? '' : o.MARKER_FONT_SIZE).trim();
  if (!mfs || isNaN(parseFloat(mfs))) p.deleteProperty('MARKER_FONT_SIZE'); else p.setProperty('MARKER_FONT_SIZE', String(parseFloat(mfs)));  // blank => 11
  var mfc = normHexColor_(o.MARKER_FONT_COLOR);         // only a valid hex survives; anything else clears -> #000000 default
  if (!mfc) p.deleteProperty('MARKER_FONT_COLOR'); else p.setProperty('MARKER_FONT_COLOR', mfc);
  p.setProperty('MARKER_BOLD', o.MARKER_BOLD ? 'true' : 'false');
  p.setProperty('MARKER_ITALIC', o.MARKER_ITALIC ? 'true' : 'false');
  // clean up the obsolete pre-mode toggle so a stale value can't linger
  p.deleteProperty('STYLE_MARKER');
  return 'ok';
}

// ---- selection-based scope ---------------------------------------------------
// Position-independent per-marker signature (survives the offset-shifting edits
// Apply makes): ok markers -> 'Q:'+normalized quote; malformed -> 'B:'+context.
// Capture and check MUST both go through this with the same cfg() so they agree.
function markerSig_(mk, c) {
  return mk.kind === 'bad' ? ('B:' + (mk.context || '')) : ('Q:' + normalize_(mk.quoteRaw, c.CASE_INSENSITIVE).norm);
}
// Walk up from a selection element to its containing paragraph/list-item (the same
// unit types collectUnits_ enumerates), so any selected run/cell resolves correctly.
function paragraphOf_(el) {
  for (var e = el, g = 0; e && g < 60; g++) {
    var t; try { t = e.getType(); } catch (x) { return null; }
    if (t === DocumentApp.ElementType.PARAGRAPH || t === DocumentApp.ElementType.LIST_ITEM) return e;
    try { e = e.getParent(); } catch (x) { return null; }
  }
  return null;
}
// Signatures of the markers in the LIVE selection. Coarse by paragraph: any selected
// element pulls in ALL markers of its paragraph. Returns:
//   null      -> no selection at all (caller falls through to the whole document)
//   {} or set -> a selection EXISTS (even if it hit no markers => process nothing)
function selectionSigSet_(c) {
  var sel; try { sel = DocumentApp.getActiveDocument().getSelection(); } catch (e) { sel = null; }
  if (!sel) return null;
  var els; try { els = sel.getRangeElements(); } catch (e) { els = null; }
  if (!els || !els.length) return {};
  var re = markerScanner_(c), sigs = {}, seen = {};
  for (var i = 0; i < els.length; i++) {
    var para = paragraphOf_(els[i].getElement());
    if (!para) continue;
    var key; try { key = para.getText(); } catch (e2) { key = null; }
    if (key != null) { if (seen[key]) continue; seen[key] = true; }   // skip re-parsing same paragraph text
    var ms = parseMarkersInPara_(para, re);
    for (var j = 0; j < ms.length; j++) sigs[markerSig_(ms[j], c)] = true;
  }
  return sigs;
}
function sigCount_(sigs) { return sigs ? Object.keys(sigs).length : 0; }

// Snapshot the selection for the Review session (a modeless dialog can't read the
// live selection). Chunked to stay under the ~9KB-per-property limit. Writing null
// CLEARS it, so a stale snapshot can never survive into a new no-selection session.
function writeSelectionSnapshot_(sigs) {
  var p = PropertiesService.getDocumentProperties();
  var old = parseInt(p.getProperty('SELECTION_SNAPSHOT_N') || '0', 10) || 0;
  for (var i = 0; i < old; i++) p.deleteProperty('SELECTION_SNAPSHOT_' + i);
  p.deleteProperty('SELECTION_SNAPSHOT_N');
  if (sigs == null) return;
  var json = JSON.stringify(Object.keys(sigs)), CH = 8000, parts = [];  // empty set -> "[]" (one chunk, truthy on read)
  for (var s = 0; s < json.length; s += CH) parts.push(json.substring(s, s + CH));
  for (var k = 0; k < parts.length; k++) p.setProperty('SELECTION_SNAPSHOT_' + k, parts[k]);
  p.setProperty('SELECTION_SNAPSHOT_N', String(parts.length));
}
function readSelectionSnapshot_() {
  var p = PropertiesService.getDocumentProperties();
  var n = p.getProperty('SELECTION_SNAPSHOT_N');
  if (n == null) return null;
  n = parseInt(n, 10) || 0;
  var json = ''; for (var i = 0; i < n; i++) json += (p.getProperty('SELECTION_SNAPSHOT_' + i) || '');
  try { var arr = JSON.parse(json || '[]'), set = {}; for (var j = 0; j < arr.length; j++) set[arr[j]] = true; return set; }
  catch (e) { return null; }
}

// Scope: a drag-selection (if any) limits to its markers; otherwise the whole document.
function markerInScope_(c, textEl, mk) {
  if (c._selectionSigs) return c._selectionSigs[markerSig_(mk, c)] === true;
  return true;
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

var LIGATURES_ = { 'ﬀ': 'ff', 'ﬁ': 'fi', 'ﬂ': 'fl', 'ﬃ': 'ffi', 'ﬄ': 'ffl', 'ﬅ': 'ft', 'ﬆ': 'st' };
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
    var lig = LIGATURES_[ch];                    // ﬁ/ﬂ/… (common PDF artifacts) -> their letters
    if (lig) {                                   // expand to multiple chars, all mapped to this raw offset
      for (var g = 0; g < lig.length; g++) { norm.push(caseInsensitive ? lig.charAt(g).toLowerCase() : lig.charAt(g)); map.push(i); }
      prevSpace = false;
      continue;
    }
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
      var indentS = 0, indentE = 0, indentF = 0;
      try { indentS = para.getIndentStart() || 0; } catch (e) {}
      try { indentE = para.getIndentEnd() || 0; } catch (e) {}
      try { indentF = para.getIndentFirstLine() || 0; } catch (e) {}
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

// ---- PDF refinement: gDoc gives the ballpark line; the PDF fine-tunes it --------
// Char-index bounds [lo,hi) of all chars whose PDF line is within ±win of centerLine.
// lineNoArr is non-decreasing (lines appended in document order), so use lower/upper
// bounds — tolerant of repeats (each line spans many chars) and gaps (missing line #s).
function lineWindowRange_(lineNoArr, centerLine, win) {
  var n = lineNoArr.length;
  if (!n) return { lo: 0, hi: 0 };
  var loLine = centerLine - win, hiLine = centerLine + win, a, b, mid;
  a = 0; b = n; while (a < b) { mid = (a + b) >> 1; if (lineNoArr[mid] < loLine) a = mid + 1; else b = mid; }
  var lo = a;
  a = 0; b = n; while (a < b) { mid = (a + b) >> 1; if (lineNoArr[mid] <= hiLine) a = mid + 1; else b = mid; }
  var hi = a;
  if (hi < lo) hi = lo;
  return { lo: lo, hi: hi };
}
// Given a gDoc estimate (estStartLine/estEndLine), return the PDF line range for a quote
// piece: (1) any EXACT occurrence anywhere, pick the one NEAREST the estimate; else
// (2) FUZZY match within ±PDF_WINDOW lines of the estimate; else (3) the gDoc estimate
// itself (estimated:true). So a quote found in the doc never lacks a line number.
function refinePieceWithPdf_(exactIndex, pieceNorm, estStartLine, estEndLine, c) {
  if (!exactIndex || !exactIndex.concat || !pieceNorm) return { start: estStartLine, end: estEndLine, estimated: true };
  var win = (c && c.PDF_WINDOW != null) ? c.PDF_WINDOW : 40;
  var fmin = (c && c.PDF_FUZZY_MIN != null) ? c.PDF_FUZZY_MIN : 0.5;
  var est0 = (estStartLine == null) ? 0 : estStartLine;   // no gDoc estimate (cross-paragraph) -> take first occurrence
  var lineNo = exactIndex.lineNo, concat = exactIndex.concat;
  // (1) exact occurrences (global), nearest to the estimate
  var bestPos = -1, bestDist = Infinity, from = 0, p;
  while ((p = concat.indexOf(pieceNorm, from)) >= 0) {
    var d = Math.abs(lineNo[p] - est0);
    if (d < bestDist) { bestDist = d; bestPos = p; }
    from = p + 1;
  }
  if (bestPos >= 0) return { start: lineNo[bestPos], end: lineNo[bestPos + pieceNorm.length - 1], estimated: false };
  // (2) fuzzy within the window around the estimate
  var rng = lineWindowRange_(lineNo, est0, win);
  if (rng.hi > rng.lo) {
    var slice = concat.substring(rng.lo, rng.hi);
    var f = fuzzyInParagraph_(slice, pieceNorm, words_(pieceNorm));
    if (f.found && f.sim >= fmin) {
      var s = rng.lo + f.S, e = rng.lo + f.E - 1;
      if (e < s) e = s; if (e >= lineNo.length) e = lineNo.length - 1;
      return { start: lineNo[s], end: lineNo[e], estimated: false };
    }
  }
  // (3) gDoc fallback
  return { start: estStartLine, end: estEndLine, estimated: true };
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

// Locate a (possibly split) quote -> first-piece start .. last-piece end. The gDoc/sim
// index decides EXISTENCE and the ballpark line; in exact mode each piece is then
// fine-tuned against the PDF (estimated:true if the PDF couldn't pin it). Used by the
// auto Update and Preview paths. `exactIndex` may be null (sim mode / no PDF file).
function lookupLineRange_(simIndex, exactIndex, mode, quoteRaw, c) {
  var pieces = splitQuote_(quoteRaw);
  if (!pieces.length) return { found: false };
  var firstStart = null, lastEnd = null, count = 0, anyEst = false;
  for (var i = 0; i < pieces.length; i++) {
    var pn = normalize_(pieces[i].raw, c.CASE_INSENSITIVE).norm;
    var est = locate_(simIndex, pn);                 // gDoc estimate (ballpark + existence)
    var ps, pe, pest;
    if (mode === 'exact' && exactIndex) {
      // The PDF concat joins lines with spaces, so it can find a quote that crosses a
      // hard paragraph break (which the sim concat can't). So: refine via the PDF, using
      // the sim estimate as the ballpark when available; if the PDF can't pin it AND the
      // sim never located it either, the piece is genuinely absent.
      var r = refinePieceWithPdf_(exactIndex, pn, est.found ? est.startLine : null, est.found ? est.endLine : null, c);
      if (r.estimated && !est.found) return { found: false, failedPiece: i };
      ps = r.start; pe = r.end; pest = r.estimated;
      count += (est.found ? est.count : 1);
    } else {
      if (!est.found) return { found: false, failedPiece: i };
      ps = est.startLine; pe = est.endLine; pest = false;
      count += est.count;
    }
    if (i === 0) firstStart = ps;
    lastEnd = pe;
    if (pest) anyEst = true;
  }
  return { found: true, startLine: firstStart, endLine: lastEnd, count: count, pieces: pieces.length, estimated: anyEst };
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

// True when a manuscript core can be spliced into the marker body WITHOUT breaking the
// marker on a later re-parse. It must contain neither a [...] elision (which splitQuote_
// would read as a phantom separator, inflating the piece count and silently defeating the
// 1:1 alignment guard) nor the CLOSING quote mark (which markerRe_ would treat as the end
// of the quote, truncating the marker). When unsafe we keep the user's piece verbatim.
function coreInjectable_(core, c) {
  if (hasElision_(core)) return false;
  var closeRe = new RegExp(bracketRe_((c && c.QUOTE_CLOSE) || '"', false));
  return !closeRe.test(core);
}

// "Use manuscript" on a split [...] quote: rebuild the quote body so each FUZZY piece
// is corrected to its manuscript wording (the same text the diff compared against),
// while KEEPING exact pieces verbatim and preserving the user's [...] separators and the
// whitespace around each piece. pieceCands[i] (from locateCandidate2_) lines up 1:1 with
// splitQuote_(quoteRaw)[i]. Returns the new body, or null if it can't rebuild safely.
function rebuildSplitQuote_(quoteRaw, pieceCands, primary, c) {
  var pieces = splitQuote_(quoteRaw);
  if (!pieces.length || !pieceCands || pieces.length !== pieceCands.length) return null;
  var out = '', cursor = 0;
  for (var i = 0; i < pieces.length; i++) {
    out += quoteRaw.substring(cursor, pieces[i].offset);     // inter-piece text incl. the [...] separators
    var raw = pieces[i].raw, pc = pieceCands[i];
    var fuzzy = pc && pc.status === 'fuzzy' && pc.srcParaIndex != null && primary.paras[pc.srcParaIndex] &&
                typeof pc.srcRawStart === 'number' && typeof pc.srcRawEndIncl === 'number' &&
                pc.srcRawEndIncl >= pc.srcRawStart;
    var core = fuzzy ? primary.paras[pc.srcParaIndex].raw.substring(pc.srcRawStart, pc.srcRawEndIncl + 1) : null;
    if (fuzzy && coreInjectable_(core, c)) {
      var lead = (raw.match(/^\s*/) || [''])[0];             // keep the user's spacing around the piece
      var trail = (raw.match(/\s*$/) || [''])[0];
      out += lead + core + trail;                            // manuscript wording for the differing piece
    } else {
      out += raw;                                            // exact / unresolved / unsafe core: keep the user's text
    }
    cursor = pieces[i].offset + raw.length;
  }
  out += quoteRaw.substring(cursor);
  return out;
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
function markerScanner_(c) { return { full: markerRe_(c), hook: hookRe_(c), openTest: openTestRe_(c), ola: openLookahead_(c), hookLen: (c.HOOK || '@AUTOLINE').length }; }
// The ~20 characters right after the hook — shown for malformed markers so the
// user can locate the broken marker in the document.
function afterHook_(raw, hookStart, hookLen) { return raw.substring(hookStart + hookLen, hookStart + hookLen + 20).replace(/\s+/g, ' '); }

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
      matchStart: start, prefixLen: prefixLen, openLen: openLen, closeLen: closeLen,
      oldNum: whole.slice(0, prefixLen).replace(/\s+$/, ''),  // original prefix (hook + any old number/colon)
      quoteRaw: quote, quoteRawStart: quoteRawStart, quoteRawEnd: quoteRawStart + quote.length
    };
    if (!quote || !/\S/.test(quote)) { e.kind = 'bad'; e.reason = 'empty quote'; e.quoteRawEnd = start + prefixLen; e.context = afterHook_(raw, start, scan.hookLen); }
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
      oldNum: hk, context: afterHook_(raw, p.index, scan.hookLen),
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
// ---- styles read once per formatting RUN, not per character (the speed win) ----
// Uses the same getters as before (isBold/getForegroundColor/…) — just at run
// boundaries — so behavior matches the known-good per-character version exactly.
function styleRuns_(t) {
  var idx; try { idx = t.getTextAttributeIndices(); } catch (e) { idx = [0]; }
  if (!idx || !idx.length) idx = [0];
  var runs = [];
  for (var i = 0; i < idx.length; i++) {
    var o = idx[i], bold = false, ital = false, col = null, al = null;
    try { bold = t.isBold(o); } catch (e) {}
    try { ital = t.isItalic(o); } catch (e) {}
    try { col = t.getForegroundColor(o); } catch (e) {}
    try { al = t.getTextAlignment(o); } catch (e) {}
    runs.push({ start: o, bold: !!bold, italic: !!ital, color: col || '#000000', align: al || DocumentApp.TextAlignment.NORMAL });
  }
  return runs;
}
function styleAt_(runs, o) {            // run whose start <= o (runs sorted ascending)
  var lo = 0;
  for (var i = 0; i < runs.length; i++) { if (runs[i].start <= o) lo = i; else break; }
  return runs[lo];
}
function eqStyle_(a, b) { return a && b && a.bold === b.bold && a.italic === b.italic && a.color === b.color && a.align === b.align; }
function writeStyle_(t, a, b, st) {     // one write per run, using the original setters
  if (b < a) return;
  try { t.setBold(a, b, st.bold); } catch (e) {}
  try { t.setItalic(a, b, st.italic); } catch (e) {}
  try { t.setForegroundColor(a, b, st.color); } catch (e) {}
  try { t.setTextAlignment(a, b, st.align); } catch (e) {}
}

// Copy formatting char-by-char from a manuscript range to the response quote,
// aligning by normalized index. Source styles are read once per run (cached).
function copyFormatting_(srcTextEl, srcMap, srcNormPos, dstTextEl, dstQuoteRawStart, dstQuoteRaw, qnorm, caseInsensitive) {
  var dn = normalize_(dstQuoteRaw, caseInsensitive);
  var qlen = Math.min(qnorm.length, dn.norm.length);
  var srcRuns = styleRuns_(srcTextEl);
  var runStartDst = -1, runEndDst = -1, runStyle = null;
  for (var k = 0; k < qlen; k++) {
    var srcRaw = srcMap[srcNormPos + k];
    var dstRaw = dstQuoteRawStart + dn.map[k];
    var st = styleAt_(srcRuns, srcRaw);
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
  var runs = styleRuns_(srcTextEl);
  var runStart = 0, cur = styleAt_(runs, srcStart);
  for (var i = 1; i <= n; i++) {
    var st = styleAt_(runs, srcStart + i);
    if (!eqStyle_(st, cur)) { writeStyle_(dstTextEl, dstStart + runStart, dstStart + i - 1, cur); runStart = i; cur = st; }
  }
  writeStyle_(dstTextEl, dstStart + runStart, dstStart + n, cur);
}

// ---- marker scaffold style (hook + number + colon + quotation marks) ----------
// Convert a Docs API rgbColor ({red,green,blue} floats 0..1; missing component => 0; the empty
// object {} that Docs returns for pure black => #000000) to a #rrggbb string.
function rgbToHex_(rgb) {
  rgb = rgb || {};
  function h(v) { v = Math.round((v || 0) * 255); if (v < 0) v = 0; if (v > 255) v = 255; var s = v.toString(16); return s.length < 2 ? '0' + s : s; }
  return '#' + h(rgb.red) + h(rgb.green) + h(rgb.blue);
}

// Is the Docs API advanced service enabled? (Needed for the 'named' style mode.)
function docsApiAvailable_() { return (typeof Docs !== 'undefined' && Docs && Docs.Documents) ? true : false; }

// Read the document-wide "Normal text" NAMED STYLE via the Docs Advanced Service (Docs API v1)
// — the definition the user sets in Format > Paragraph styles > Normal text, which applies
// document-wide. Returns {family,size,color} with '' for anything unreadable (and all-blank when
// the service isn't enabled or the style is absent). Used by the 'named' mode + the Settings dialog.
function namedNormalStyle_() {
  var out = { family: '', size: '', color: '' };
  if (!docsApiAvailable_()) return out;
  try {
    var id = DocumentApp.getActiveDocument().getId();
    var doc = Docs.Documents.get(id, { fields: 'namedStyles' });   // partial response: just the styles
    var styles = doc && doc.namedStyles && doc.namedStyles.styles;
    if (!styles) return out;
    var ts = null;
    for (var i = 0; i < styles.length; i++) {
      if (styles[i] && styles[i].namedStyleType === 'NORMAL_TEXT') { ts = styles[i].textStyle; break; }
    }
    if (!ts) return out;
    if (ts.weightedFontFamily && ts.weightedFontFamily.fontFamily) out.family = ts.weightedFontFamily.fontFamily;
    if (ts.fontSize && ts.fontSize.magnitude != null) out.size = ts.fontSize.magnitude;
    if (ts.foregroundColor && ts.foregroundColor.color && ts.foregroundColor.color.rgbColor)
      out.color = rgbToHex_(ts.foregroundColor.color.rgbColor);   // skip theme-only colors (no rgbColor)
  } catch (e) {}
  return out;
}

// Resolve the scaffold style for the marker (hook + line number(s) + colon + quotation marks —
// NOT the quoted text). Stashed once per run on c._markerStyle. The mode is the user's explicit
// choice from Settings:
//   'insert' -> null   (no restyling — the scaffold keeps the formatting where it sits)
//   'named'  -> the document-wide Normal text style (Docs API), or null if unavailable/empty
//              (so it harmlessly behaves like 'insert' until the Docs API is enabled)
//   'set'    -> the explicit font/size/color + bold/italic (defaults: Arial / 11 / #000000 / regular)
function resolveMarkerStyle_(c) {
  var mode = c.MARKER_STYLE_MODE || 'set';
  if (mode === 'insert') return null;
  if (mode === 'named') {
    var nrm = namedNormalStyle_();
    var fam = nrm.family || null;
    var size = (nrm.size !== '' && nrm.size != null) ? nrm.size : null;
    var col = nrm.color || null;
    if (fam == null && size == null && col == null) return null;   // API off / style empty -> keep insert formatting
    return { family: fam, size: size, color: col };
  }
  return {                                                          // 'set'
    family: c.MARKER_FONT_FAMILY || 'Arial',
    size: (c.MARKER_FONT_SIZE && c.MARKER_FONT_SIZE > 0) ? c.MARKER_FONT_SIZE : 11,
    color: c.MARKER_FONT_COLOR || '#000000',
    bold: !!c.MARKER_BOLD,
    italic: !!c.MARKER_ITALIC
  };
}
function applyTextStyle_(textEl, a, b, style) {
  if (!style || a == null || b == null || b < a) return;
  try { if (style.family) textEl.setFontFamily(a, b, style.family); } catch (e) {}
  try { if (style.size)   textEl.setFontSize(a, b, style.size); } catch (e) {}
  try { if (style.color)  textEl.setForegroundColor(a, b, style.color); } catch (e) {}
  try { if (typeof style.bold === 'boolean')   textEl.setBold(a, b, style.bold); } catch (e) {}     // 'set' mode: force regular/bold
  try { if (typeof style.italic === 'boolean') textEl.setItalic(a, b, style.italic); } catch (e) {} // 'set' mode: force regular/italic
}
// Style the scaffold of a marker AT ITS CURRENT LAYOUT (after any prefix rewrite):
// prefix + opening mark = [matchStart, openEnd]; closing mark = [closeStart, closeEnd].
function styleMarkerScaffold_(textEl, matchStart, openEnd, closeStart, closeEnd, style) {
  applyTextStyle_(textEl, matchStart, openEnd, style);
  applyTextStyle_(textEl, closeStart, closeEnd, style);
}

// Set or clear a background color on a text range. The range overload accepts a
// hex string OR null (to clear); setAttributes is a fallback. (setAttributes with a
// null value can silently no-op, which is why the range overload is tried first.)
function bg_(t, a, b, color) {
  if (a == null || b == null || b < a) return;
  try { t.setBackgroundColor(a, b, color); return; } catch (e) {}
  try { var x = {}; x[DocumentApp.Attribute.BACKGROUND_COLOR] = color; t.setAttributes(a, b, x); } catch (e2) {}
}

function isFlagShade_(bg, c) {              // is this background one of OUR flag colors?
  if (!bg) return false;
  bg = String(bg).toLowerCase();
  return bg === String(c.LOWCONF_COLOR).toLowerCase() || bg === String(c.NOTFOUND_COLOR).toLowerCase();
}

// Clear the tool's own yellow/red flags from every marker prefix (leaves any
// highlights YOU applied untouched). Returns how many were cleared.
function clearFlags_(c) {
  var body = DocumentApp.getActiveDocument().getBody();
  var units = [];
  for (var i = 0; i < body.getNumChildren(); i++) collectUnits_(body.getChild(i), units, false);
  var re = markerScanner_(c), n = 0;
  for (var u = 0; u < units.length; u++) {
    var t = units[u].el.editAsText();
    var ms = parseMarkersInPara_(units[u].el, re);
    for (var k = 0; k < ms.length; k++) {
      var mk = ms[k], bgc;
      try { bgc = t.getBackgroundColor(mk.matchStart); } catch (e) { bgc = null; }
      if (isFlagShade_(bgc, c)) { bg_(t, mk.matchStart, mk.matchStart + mk.prefixLen - 1, null); n++; }
    }
  }
  return n;
}

// ============================================================================
// Main update routine
// ============================================================================
function runUpdate_(mode) {
  var ui = DocumentApp.getUi();
  var c = cfg();
  if (!c.PRIMARY_DOC_ID) { ui.alert('Set the primary doc first (QuoteFinder → Settings…).'); return; }
  c._selectionSigs = selectionSigSet_(c);          // LIVE selection scope (menu command)
  c._markerStyle = resolveMarkerStyle_(c);         // scaffold style (hook/number/colon/quote marks)

  var primary, exactIndex = null, simIndex = null;
  try {
    primary = buildPrimaryIndex_(c, true);       // sim layout always (gDoc estimate; exact mode refines it)
    simIndex = primary.sim;
    if (mode === 'exact') exactIndex = buildExactIndex_(c, loadExactMap_(c));
  } catch (e) { ui.alert('Setup problem:\n' + e.message); return; }

  if (!c._selectionSigs) clearFlags_(c);   // clean slate (whole-doc only; selection leaves other markers alone)

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
        bg_(textEl, mk.matchStart, mk.matchStart + mk.prefixLen - 1, c.NOTFOUND_COLOR);
        problems.push('MALFORMED (' + mk.reason + '): ' + c.HOOK + mk.context + ' …');
        continue;
      }
      // line number(s): split quotes ([...]) span first-piece-start .. last-piece-end
      var lineInfo = lookupLineRange_(simIndex, exactIndex, mode, mk.quoteRaw, c);

      if (!lineInfo.found) {
        stats.notFound++;
        bg_(textEl, mk.matchStart, mk.matchStart + mk.prefixLen - 1, c.NOTFOUND_COLOR);
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
      var low = (mode === 'sim')
        ? (fp0 && lineNumberLowConf_(primary, normalize_(fp0.raw, c.CASE_INSENSITIVE).norm))   // sim: table/image drift
        : !!lineInfo.estimated;                                                                // exact: PDF couldn't pin it
      if (low) stats.lowConf++;
      bg_(textEl, mk.matchStart, mk.matchStart + newPrefix.length - 1, low ? c.LOWCONF_COLOR : null);
      if (c._markerStyle) {                          // style scaffold; body unchanged here, so use the original length
        var openLen = (mk.openLen != null) ? mk.openLen : (mk.quoteRawStart - mk.matchStart - mk.prefixLen);
        var closeLen = (mk.closeLen != null) ? mk.closeLen : 1;
        var base = mk.matchStart + newPrefix.length, bodyLen = mk.quoteRawEnd - mk.quoteRawStart;
        var cStart = base + openLen + bodyLen;
        styleMarkerScaffold_(textEl, mk.matchStart, base + openLen - 1, cStart, cStart + closeLen - 1, c._markerStyle);
      }
      stats.updated++;
    }
  }

  var scopeNote = c._selectionSigs ? (' [Scope: ' + sigCount_(c._selectionSigs) + ' selected quote(s)]') : '';
  var msg = 'QuoteFinder (' + (mode === 'sim' ? 'gDoc (Estimated Line #)' : 'PDF/JSON (Exact Line #)') + ')' + scopeNote + ' done.\n\n' +
    'Updated: ' + stats.updated + '\n' +
    (stats.lowConf ? 'Low-confidence (yellow): ' + stats.lowConf + '\n' : '') +
    (stats.notFound ? 'Not found (red): ' + stats.notFound + '\n' : '') +
    (stats.malformed ? 'Malformed marker (red): ' + stats.malformed + '\n' : '') +
    (stats.multi ? 'Multiple matches: ' + stats.multi + '\n' : '') +
    (stats.noFormat ? 'Formatting skipped: ' + stats.noFormat + '\n' : '');
  if (c._selectionSigs && sigCount_(c._selectionSigs) === 0)
    msg += '\nYour selection contained no @AUTOLINE markers, so nothing was changed. Select the quote line(s) you want, then run again.\n';
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
  c._selectionSigs = selectionSigSet_(c);          // LIVE selection scope (menu command)

  var simIndex, label, note = '', mode;
  try { simIndex = buildPrimaryIndex_(c, true).sim; }              // gDoc estimate — always
  catch (e2) { ui.alert('Setup problem:\n' + e2.message); return; }
  var exactIndex = null;
  try {                                                           // prefer the exact (PDF/JSON) engine — ground truth
    exactIndex = buildExactIndex_(c, loadExactMap_(c)); mode = 'exact'; label = 'PDF/JSON (Exact Line #)';
  } catch (e) {                                                   // no usable PDF/JSON file -> gDoc estimate
    mode = 'sim'; label = 'gDoc (Estimated Line #)';
    note = 'Using estimated numbers — no usable PDF/JSON line-number file. ' +
           'Set one in QuoteFinder → Settings… for exact numbers.';
  }

  var rows = collectPreviewRows_(c, simIndex, exactIndex, mode);
  showPreviewHtml_(rows, label, note, sigCount_(c._selectionSigs), c._selectionSigs != null);
}

function collectPreviewRows_(c, simIndex, exactIndex, mode) {
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
      if (mk.kind === 'bad') { rows.push({ old: c.HOOK + mk.context + ' …', quote: '⚠ ' + mk.reason, lines: 'malformed', found: false }); continue; }
      var info = lookupLineRange_(simIndex, exactIndex, mode, mk.quoteRaw, c);   // split quotes -> first..last range
      rows.push({ old: mk.oldNum, quote: snippet_(mk.quoteRaw), lines: fmtLines_(info), found: !!info.found, estimated: !!info.estimated });
    }
  }
  return rows;
}

function fmtLines_(info) {
  if (!info) return '—';
  if (!info.found) return 'not found';
  return info.startLine === info.endLine ? ('' + info.startLine) : (info.startLine + '-' + info.endLine);
}

function showPreviewHtml_(rows, engineLabel, note, selCount, selActive) {
  var h = ['<style>body{font:13px Arial;margin:8px}table{border-collapse:collapse;width:100%}',
    'td,th{border:1px solid #ccc;padding:3px 6px;text-align:left;vertical-align:top}',
    'th{background:#f0f0f0}.no{background:#f4cccc}.est{background:#fff2cc}</style>'];
  var scopeTxt = selActive ? (' — Scope: ' + (selCount || 0) + ' selected quote(s)') : '';
  h.push('<p><b>' + rows.length + ' marker(s)' + scopeTxt + '.</b> ' +
    'Line numbers from <b>' + esc_(engineLabel) + '</b>. Dry run — nothing was changed. ' +
    '<span style="background:#fff2cc">Yellow</span> = PDF couldn\'t pin it, gDoc estimate used.</p>');
  if (selActive && !rows.length) h.push('<p style="color:#b06000">Your selection contained no @AUTOLINE markers. Select the quote line(s) you want, then run Preview again.</p>');
  if (note) h.push('<p style="color:#b06000">' + esc_(note) + '</p>');
  h.push('<table><tr><th>#</th><th>old</th><th>quote</th><th>new # (' + esc_(engineLabel) + ')</th></tr>');
  for (var i = 0; i < rows.length; i++) {
    var r = rows[i];
    var cls = !r.found ? ' class="no"' : (r.estimated ? ' class="est"' : '');
    h.push('<tr><td>' + (i + 1) + '</td><td>' + esc_(r.old) + '</td><td>' + esc_(r.quote) + '</td>' +
      '<td' + cls + '>' + esc_(r.lines) + (r.estimated ? ' (est)' : '') + '</td></tr>');
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
function locateCandidate_(primary, mode, exactIndex, qnorm, c) {
  if (!qnorm) return { status: 'notfound' };
  for (var i = 0; i < primary.paras.length; i++) {           // exact in a paragraph
    var pos = primary.paras[i].norm.indexOf(qnorm);
    if (pos >= 0) return finishCandidate_(primary, mode, exactIndex, primary.paras[i], pos, pos + qnorm.length, qnorm, 1, 'exact', c);
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
  return finishCandidate_(primary, mode, exactIndex, primary.paras[best.idx], best.S, best.E, qnorm, best.sim, 'fuzzy', c);
}

function finishCandidate_(primary, mode, exactIndex, para, S, E, qnorm, sim, status, c) {
  var candNorm = para.norm.substring(S, E);
  var lines = candidateLines_(mode, exactIndex, para, S, E, candNorm, c);
  return {
    status: status, similarity: sim, srcParaIndex: para.idx, srcNormPos: S,
    srcRawStart: para.map[S], srcRawEndIncl: para.map[E - 1],
    startLine: lines.start, endLine: lines.end, candNorm: candNorm,
    lineEstimated: !!lines.estimated
  };
}
// gDoc/sim estimate from the doc layout; in exact mode, fine-tune via the PDF.
function candidateLines_(mode, exactIndex, para, S, E, candNorm, c) {
  var estStart = para.startLine + lineForRaw_(para.lineStarts, para.map[S]);
  var estEnd   = para.startLine + lineForRaw_(para.lineStarts, para.map[E - 1]);
  if (mode === 'sim') return { start: estStart, end: estEnd, estimated: false };
  var r = refinePieceWithPdf_(exactIndex, candNorm, estStart, estEnd, c);
  return { start: r.start, end: r.end, estimated: r.estimated };
}

// Split-quote aware candidate: a quote with [...] is matched piece-by-piece; the
// reported range spans the first piece's start to the last piece's end.
function locateCandidate2_(primary, mode, exactIndex, quoteRaw, c) {
  if (!hasElision_(quoteRaw)) return locateCandidate_(primary, mode, exactIndex, normalize_(quoteRaw, c.CASE_INSENSITIVE).norm, c);
  var pieces = splitQuote_(quoteRaw);
  if (!pieces.length) return { status: 'notfound' };
  var cands = [];
  for (var i = 0; i < pieces.length; i++) {
    var cand = locateCandidate_(primary, mode, exactIndex, normalize_(pieces[i].raw, c.CASE_INSENSITIVE).norm, c);
    if (cand.status === 'notfound') return { status: 'notfound' };
    cands.push(cand);
  }
  var first = cands[0], last = cands[cands.length - 1], allExact = true, minSim = 1, anyEst = false;
  for (var k = 0; k < cands.length; k++) {
    if (!(cands[k].status === 'exact' && cands[k].startLine != null)) allExact = false;
    if (cands[k].similarity < minSim) minSim = cands[k].similarity;
    if (cands[k].lineEstimated) anyEst = true;
  }
  return {
    status: (allExact && first.startLine != null && last.endLine != null) ? 'exact' : 'fuzzy',
    similarity: minSim, split: true, pieces: cands, lineEstimated: anyEst,
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

// ---- formatted-HTML + diff rendering (styles read once per run) ------------
function rangeToHtml_(textEl, start, endIncl) {
  if (endIncl < start) return '<em>(empty)</em>';
  var T = textEl.getText();
  var runs = styleRuns_(textEl);
  var html = [];
  function emit(a, b, st) {
    var open = [], close = [];
    if (st.align === DocumentApp.TextAlignment.SUPERSCRIPT) { open.push('<sup>'); close.unshift('</sup>'); }
    else if (st.align === DocumentApp.TextAlignment.SUBSCRIPT) { open.push('<sub>'); close.unshift('</sub>'); }
    if (st.bold) { open.push('<b>'); close.unshift('</b>'); }
    if (st.italic) { open.push('<i>'); close.unshift('</i>'); }
    if (st.color && st.color !== '#000000') { open.push('<span style="color:' + st.color + '">'); close.unshift('</span>'); }
    html.push(open.join('') + esc_(T.substring(a, b + 1)) + close.join(''));
  }
  var runStart = start, cur = styleAt_(runs, start);
  for (var o = start + 1; o <= endIncl; o++) {
    var st = styleAt_(runs, o);
    if (!eqStyle_(st, cur)) { emit(runStart, o - 1, cur); runStart = o; cur = st; }
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

// Parse the old line number(s) out of a marker's oldNum prefix (HOOK + optional n|a-b + colon/ws;
// see where oldNum is set). Returns { hadNumber, start, end }. No digits (bare hook, "@AUTOLINE:")
// => hadNumber:false (a first-time numbering counts as a change); 0 is a real number.
function parseOldNum_(oldNum, c) {
  var hook = (c && c.HOOK) || '@AUTOLINE';
  var s = String(oldNum == null ? '' : oldNum);
  if (s.indexOf(hook) === 0) s = s.slice(hook.length);     // strip the hook prefix
  var m = s.match(/(\d+)\s*-\s*(\d+)|(\d+)/);               // range first, else a single number
  if (!m) return { hadNumber: false, start: null, end: null };
  if (m[1] != null) return { hadNumber: true, start: parseInt(m[1], 10), end: parseInt(m[2], 10) };
  var n = parseInt(m[3], 10);
  return { hadNumber: true, start: n, end: n };
}

// ---- server entry points called from the sidebar --------------------------
function buildReviewData(mode) {
  var c = cfg();
  c._selectionSigs = readSelectionSnapshot_();      // selection captured when the panel opened (session-fixed)
  var primary = buildPrimaryIndex_(c, true);   // always build the sim layout (gDoc estimate, even in exact mode)
  var exactIndex = (mode === 'exact') ? buildExactIndex_(c, loadExactMap_(c)) : null;

  var body = DocumentApp.getActiveDocument().getBody();
  var units = [];
  for (var i = 0; i < body.getNumChildren(); i++) collectUnits_(body.getChild(i), units, false);

  var ordinal = 0, total = 0, exactCount = 0, malformedCount = 0, items = [];
  var allMarkers = 0;                             // diagnostic: total markers seen in the doc
  var re = markerScanner_(c);
  for (var u = 0; u < units.length; u++) {
    var para = units[u].el, textEl = para.editAsText();
    var markers = parseMarkersInPara_(para, re);
    for (var mi = 0; mi < markers.length; mi++) {
      var mk = markers[mi];
      allMarkers++;
      var inScope = c._selectionSigs ? markerInScope_(c, textEl, mk) : true;   // selection limits; else whole doc
      if (!inScope) { ordinal++; continue; }       // out of scope: keep ordinals aligned (stable IDs)
      total++;
      if (mk.kind === 'bad') {                       // hook present, quote missing/unterminated/empty
        malformedCount++;
        items.push({
          ordinal: ordinal, oldNum: mk.oldNum, status: 'malformed', reason: mk.reason,
          context: c.HOOK + mk.context, similarity: 0, startLine: null, endLine: null,
          numberChanged: false,
          oldHtml: esc_(c.HOOK + mk.context) + ' <span class="muted">…</span>', newHtml: '', diffHtml: ''
        });
        ordinal++; continue;
      }
      var cand = locateCandidate2_(primary, mode, exactIndex, mk.quoteRaw, c);  // split-quote aware
      var isExact = (cand.status === 'exact' && cand.startLine != null);
      var isEstimated = !!cand.lineEstimated;       // PDF couldn't pin it -> gDoc estimate used
      if (isExact && !isEstimated) exactCount++;     // an estimated line isn't "perfect" — surface it

      var hasCand = (cand.status !== 'notfound');
      var oldN = parseOldNum_(mk.oldNum, c);          // is the line NUMBER actually moving? (for the "Line # only" tab)
      var numberChanged = hasCand && cand.startLine != null &&
        (!oldN.hadNumber || oldN.start !== cand.startLine || oldN.end !== cand.endLine);
      items.push({                                  // EVERY in-scope marker (client filters the view)
        ordinal: ordinal,
        oldNum: mk.oldNum,
        status: isExact ? 'exact' : cand.status,    // 'exact' | 'fuzzy' | 'notfound'
        similarity: Math.round((cand.similarity || 0) * 100),
        startLine: hasCand ? cand.startLine : null,
        endLine: hasCand ? cand.endLine : null,
        split: !!cand.split,
        numberChanged: numberChanged,
        lineEstimated: hasCand ? isEstimated : false,
        oldHtml: rangeToHtml_(textEl, mk.quoteRawStart, mk.quoteRawEnd - 1),
        newHtml: hasCand ? candNewHtml_(primary, cand) : '',
        diffHtml: hasCand ? wordDiffHtml_(mk.quoteRaw, candNewText_(primary, cand)) : ''
      });
      ordinal++;
    }
  }
  return { mode: mode, total: total, exactCount: exactCount, imperfectCount: total - exactCount,
           malformedCount: malformedCount, items: items,
           scope: c._selectionSigs ? 'selection' : 'all',
           selectionCount: sigCount_(c._selectionSigs), allMarkers: allMarkers };
}

// Apply one pass; edits within each paragraph go right-to-left so offsets stay
// valid. decisionsJson = JSON map { ordinal: 'use' | 'skip' } for reviewed
// markers. Default when an ordinal is absent: refresh exact matches, skip the
// rest (so "only changes" view still auto-applies the perfect ones).
function applyReviewDecisions(mode, decisionsJson) {
  var c = cfg();
  c._selectionSigs = readSelectionSnapshot_();      // same snapshot the panel scanned with
  c._markerStyle = resolveMarkerStyle_(c);          // scaffold style (hook/number/colon/quote marks)
  var dec = {};
  try { dec = JSON.parse(decisionsJson || '{}') || {}; } catch (e) { dec = {}; }

  var primary = buildPrimaryIndex_(c, true);   // always build the sim layout (gDoc estimate, even in exact mode)
  var exactIndex = (mode === 'exact') ? buildExactIndex_(c, loadExactMap_(c)) : null;

  if (!c._selectionSigs) clearFlags_(c);   // whole-doc only; selection leaves other markers alone

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
          bg_(textEl, mk.matchStart, mk.matchStart + mk.prefixLen - 1, c.NOTFOUND_COLOR);
          n.flagged++; return;
        }
        // Selection mode skips the global clearFlags_, so clear an in-scope OK marker's
        // stale red/yellow shade here (the applied path also self-cleans; this also covers
        // the skip / use-but-unresolved paths, matching whole-doc behavior).
        if (c._selectionSigs) {
          var pbg; try { pbg = textEl.getBackgroundColor(mk.matchStart); } catch (e) { pbg = null; }
          if (isFlagShade_(pbg, c)) bg_(textEl, mk.matchStart, mk.matchStart + mk.prefixLen - 1, null);
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
  var newBodyLen = mk.quoteRawEnd - mk.quoteRawStart;  // unchanged unless we replace the text (fuzzy)
  if (cand.split || hasElision_(mk.quoteRaw)) {        // split quote (with [...])
    var did = false;
    if (replaceText && cand.split && cand.pieces) {    // "Use manuscript": correct each fuzzy piece, keep the [...]
      var rebuilt = rebuildSplitQuote_(mk.quoteRaw, cand.pieces, primary, c);
      if (rebuilt != null && rebuilt !== mk.quoteRaw) {
        newBodyLen = rebuilt.length;
        textEl.deleteText(mk.quoteRawStart, mk.quoteRawEnd - 1);
        textEl.insertText(mk.quoteRawStart, rebuilt);
        if (c.COPY_FORMATTING) copyPartsFormatting_(primary, textEl, mk.quoteRawStart, rebuilt, c);
        did = true;
      }
    }
    if (!did && c.COPY_FORMATTING) copyPartsFormatting_(primary, textEl, mk.quoteRawStart, mk.quoteRaw, c);
  } else if (replaceText) {
    var src = primary.paras[cand.srcParaIndex];
    var newTextRaw = src.raw.substring(cand.srcRawStart, cand.srcRawEndIncl + 1);
    newBodyLen = newTextRaw.length;
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
  // An estimated line (PDF couldn't pin it) is flagged yellow so you know to check it;
  // otherwise drop any stale flag shade inherited by the rewritten prefix (never a
  // highlight YOU applied), so a cleanly-resolved marker self-cleans.
  if (cand.lineEstimated) {
    bg_(textEl, mk.matchStart, mk.matchStart + newPrefix.length - 1, c.LOWCONF_COLOR);
  } else {
    var pbg; try { pbg = textEl.getBackgroundColor(mk.matchStart); } catch (e) { pbg = null; }
    if (isFlagShade_(pbg, c)) bg_(textEl, mk.matchStart, mk.matchStart + newPrefix.length - 1, null);
  }
  // style the scaffold (hook + number + colon + quote marks) — current layout after the rewrite
  if (c._markerStyle) {
    var openLen = mk.openLen != null ? mk.openLen : (mk.quoteRawStart - mk.matchStart - mk.prefixLen);
    var closeLen = mk.closeLen != null ? mk.closeLen : 1;
    var base = mk.matchStart + newPrefix.length;       // first char after the rewritten prefix = opening mark
    var closeStart = base + openLen + newBodyLen;
    styleMarkerScaffold_(textEl, mk.matchStart, base + openLen - 1, closeStart, closeStart + closeLen - 1, c._markerStyle);
  }
}

// Apply ONE marker, found fresh by its ordinal (Nth marker in document order) so a
// prior edit can't have moved it — same engine as the bulk Apply, scoped to one.
// Returns {ok, msg, startLine, endLine}. Note: like Apply, this re-reads the
// manuscript, so it costs about as much per click as a full Apply.
function applyOne(mode, ordinal) {
  var c = cfg();
  c._selectionSigs = readSelectionSnapshot_();      // same snapshot the panel scanned with
  c._markerStyle = resolveMarkerStyle_(c);          // scaffold style (hook/number/colon/quote marks)
  var primary = buildPrimaryIndex_(c, true);   // always build the sim layout (gDoc estimate, even in exact mode)
  var exactIndex = (mode === 'exact') ? buildExactIndex_(c, loadExactMap_(c)) : null;

  var body = DocumentApp.getActiveDocument().getBody();
  var units = [];
  for (var i = 0; i < body.getNumChildren(); i++) collectUnits_(body.getChild(i), units, false);

  var ord = 0, re = markerScanner_(c);
  for (var u = 0; u < units.length; u++) {
    var para = units[u].el, textEl = para.editAsText();
    var markers = parseMarkersInPara_(para, re);
    for (var mi = 0; mi < markers.length; mi++) {
      var mk = markers[mi];
      if (ord++ !== ordinal) continue;                 // not this one (ordinals match buildReviewData)
      if (!markerInScope_(c, textEl, mk)) return { ok: false, msg: 'That marker is out of the current scope.' };
      if (mk.kind === 'bad') return { ok: false, msg: 'Malformed marker — fix it in the document.' };
      var cand = locateCandidate2_(primary, mode, exactIndex, mk.quoteRaw, c);
      if (cand.status === 'notfound' || cand.startLine == null)
        return { ok: false, msg: 'No manuscript match / line number for that quote.' };
      var replace = (cand.status === 'fuzzy');
      applyMarker_(textEl, mk, primary, cand, c, replace);   // edits this one marker only, then returns
      var lineStr = (cand.startLine === cand.endLine) ? ('' + cand.startLine) : (cand.startLine + '–' + cand.endLine);
      return { ok: true, startLine: cand.startLine, endLine: cand.endLine, estimated: !!cand.lineEstimated,
               msg: (replace ? 'Replaced text + number + formatting' : 'Refreshed number + formatting') + ' — line ' + lineStr +
                    (cand.lineEstimated ? ' (estimated — PDF couldn\'t pin it)' : '') + '.' };
    }
  }
  return { ok: false, msg: 'Could not find that marker — click ↻ Rescan (did the markers change?).' };
}
