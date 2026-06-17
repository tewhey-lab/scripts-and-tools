# Interactive qPCR Data Visualization

`qPCR_viewer.html` is a self-contained, browser-based app for the interactive
visualization of qPCR data from QuantStudio Excel exports. Everything runs
locally in your web browser — there is nothing to install, and no data ever
leaves your machine.

## Features

- **Configurable dashboard**: The plots live in a grid of panels you can
  **move** (drag the panel header), **resize** (drag any edge or corner), and
  **remove** (the × button). Use **+ Add view** to add another panel and
  **Reset** to restore the default layout. Each panel has its own dropdown to
  switch which view it shows.
- **Views offered match your file**: The app inspects the workbook and only
  offers the views it can actually plot. Up to six are available:
    1.  Melt Curve: Fluorescence vs. Temperature
    2.  Melt Curve: Derivative (−dF/dT) vs. Temperature
    3.  Amplification: Rn vs. Cycle
    4.  Amplification: ΔRn vs. Cycle
    5.  Amplification: Absolute Fluorescence vs. Cycle *
    6.  Amplification: ΔFluorescence vs. Cycle *

    \* The two absolute-fluorescence views require raw per-cycle dye signal,
    which lives in QuantStudio's **Multicomponent Data** sheet. When that sheet
    is present, the app uses it (and labels the panel with the reporter dye it
    detected, e.g. *SYBR*); when it is absent, those two views are simply not
    offered, since the standard Amplification Data sheet only contains the
    normalized `Rn` / `Delta Rn`.
- **Filtering**: Toggle which traces are displayed by **Target Name**,
  **Sample Name**, and **Well ID**, with "All" / "None" buttons. **Smart
  Groups** auto-detects common patterns in your sample names so you can show
  related samples in one click.
- **Display options**: Adjustable line width, an optional CT threshold line,
  and a log-scale toggle for the amplification y-axis. CT values are marked on
  the amplification curves.
- **Plot interaction**: Hover for detailed tooltips (well, sample, target, CT,
  Tm), click-and-drag to zoom, and double-click a panel to reset its view.
- **Dark mode**: A light/dark theme toggle.
- **Export**: Save every panel as a PNG or SVG image.

## How to Use

1.  **Open `qPCR_viewer.html`** in any modern web browser (Chrome, Firefox,
    Safari, or Edge). Double-click the file, or drag it into a browser window.

2.  **Load your data** by dragging your QuantStudio Excel export (`.xls` or
    `.xlsx`) onto the upload area, or by clicking it to browse for the file.

3.  The dashboard renders automatically with one panel per available view.
    Rearrange, resize, add, or remove panels to taste; use the sidebar to
    filter traces and the toolbar to export images or toggle dark mode.

> **Note:** The app loads its charting, spreadsheet, and dashboard libraries
> (Plotly, SheetJS, Gridstack) from public CDNs, so an internet connection is
> required the first time you open it.
