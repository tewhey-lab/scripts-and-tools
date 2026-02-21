"""Streamlit web UI for plasmid_align.py."""

import os
import tempfile
import requests
import streamlit as st

from plasmid_align import run_pipeline

# ── GitHub plasmid library ───────────────────────────────────────────────────
GITHUB_API_URL = (
    "https://api.github.com/repos/tewhey-lab/scripts-and-tools"
    "/contents/sequence-tools/plasmids"
)
GITHUB_RAW_BASE = (
    "https://raw.githubusercontent.com/tewhey-lab/scripts-and-tools"
    "/main/sequence-tools/plasmids"
)


@st.cache_data(ttl=600)
def _fetch_plasmid_catalog():
    """Return a dict {display_name: filename} of FASTA files in the repo."""
    try:
        resp = requests.get(GITHUB_API_URL, timeout=10)
        resp.raise_for_status()
        entries = resp.json()
        catalog = {}
        for entry in entries:
            name = entry["name"]
            if name.lower().endswith((".fasta", ".fa", ".fna")):
                display = os.path.splitext(name)[0]
                catalog[display] = name
        return catalog
    except Exception:
        return {}


def _download_plasmid(filename, dest_dir):
    """Download a single FASTA file from the repo and return its local path."""
    url = f"{GITHUB_RAW_BASE}/{filename}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    local_path = os.path.join(dest_dir, filename)
    with open(local_path, "wb") as f:
        f.write(resp.content)
    return local_path


def _pdf_to_images(pdf_path):
    """Convert each page of a PDF to a PNG byte buffer using PyMuPDF."""
    images = []
    try:
        import fitz  # PyMuPDF
        doc = fitz.open(pdf_path)
        for page in doc:
            pix = page.get_pixmap(dpi=150)
            images.append(pix.tobytes("png"))
        doc.close()
    except ImportError:
        pass  # PyMuPDF not installed; fall back to PDF download only
    return images


# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(page_title="Plasmid Aligner", layout="wide")
st.title("Plasmid Aligner")
st.markdown(
    "Align Oxford Nanopore reads to a library of plasmid FASTA sequences. "
    "Upload your files, adjust parameters, and click **Run**."
)

# ── Sidebar: file uploads + parameters ───────────────────────────────────────
with st.sidebar:
    st.header("Input files")
    fastq_file = st.file_uploader(
        "FASTQ file", type=["fastq", "fq", "gz"],
        help="Oxford Nanopore FASTQ file (plain or gzip-compressed).",
    )

    # ── Plasmid reference selection ──────────────────────────────────────────
    st.subheader("Plasmid references")
    plasmid_catalog = _fetch_plasmid_catalog()

    selected_plasmids = []
    if plasmid_catalog:
        selected_plasmids = st.multiselect(
            "Select from library",
            options=sorted(plasmid_catalog.keys()),
            help="Choose one or more plasmid references from the lab library.",
        )
    else:
        st.caption("⚠️ Could not load plasmid library from GitHub.")

    fasta_files = st.file_uploader(
        "Or upload custom FASTA(s)",
        type=["fasta", "fa", "fna"],
        accept_multiple_files=True,
        help="Upload additional plasmid FASTA files. "
             "These are combined with any selections above.",
    )

    has_references = len(selected_plasmids) > 0 or len(fasta_files) > 0

    st.header("Parameters")
    min_qual = st.number_input(
        "Min avg quality (Phred)", value=10.0, min_value=0.0, step=1.0,
        help="Reads below this average quality are filtered out.",
    )
    min_length = st.number_input(
        "Min read length (bp)", value=2000, min_value=0, step=100,
        help="Reads shorter than this are filtered out.",
    )
    ambiguity_delta = st.number_input(
        "Ambiguity delta", value=0.05, min_value=0.0, max_value=1.0, step=0.01,
        help="If top two plasmid scores differ by less than this fraction, "
             "the read is called ambiguous.",
    )
    circular = st.checkbox(
        "Circular mode", value=True,
        help="Double each reference to capture reads spanning the "
             "linearisation junction.",
    )
    fill_n = st.checkbox(
        "Fill N regions", value=True,
        help="Fill N-masked regions in references using the most homologous "
             "sequence from another reference in the library. Improves "
             "alignment coverage through masked regions.",
    )
    threads = st.slider("Threads", min_value=1, max_value=4, value=2)

# ── Main area: run button + results ──────────────────────────────────────────
can_run = fastq_file is not None and has_references

if st.button("Run Alignment", disabled=not can_run, type="primary"):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write uploaded FASTQ to disk
        fastq_path = os.path.join(tmpdir, fastq_file.name)
        with open(fastq_path, "wb") as f:
            f.write(fastq_file.getbuffer())

        fasta_paths = []

        # Download selected library plasmids
        if selected_plasmids:
            with st.spinner("Downloading plasmid references …"):
                for display_name in selected_plasmids:
                    filename = plasmid_catalog[display_name]
                    p = _download_plasmid(filename, tmpdir)
                    fasta_paths.append(p)

        # Write uploaded FASTA files to disk
        for uf in fasta_files:
            p = os.path.join(tmpdir, uf.name)
            with open(p, "wb") as f:
                f.write(uf.getbuffer())
            fasta_paths.append(p)

        output_dir = os.path.join(tmpdir, "results")
        os.makedirs(output_dir)

        progress_bar = st.progress(0, text="Starting alignment …")
        status_text = st.empty()

        def _progress(n_passed, msg):
            status_text.text(msg)

        try:
            result = run_pipeline(
                fastq_path=fastq_path,
                fasta_paths=fasta_paths,
                min_qual=min_qual,
                min_length=min_length,
                ambiguity_delta=ambiguity_delta,
                circular=circular,
                threads=threads,
                output_dir=output_dir,
                no_plots=False,
                fill_n=fill_n,
                progress_callback=_progress,
            )
        except Exception as exc:
            st.error(f"Pipeline failed: {exc}")
            st.stop()

        progress_bar.progress(100, text="Done!")

        # ── Display results ──────────────────────────────────────────────
        st.success("Alignment complete!")

        # ── Mapping statistics ───────────────────────────────────────────
        counts = result["counts"]
        summary_df = result["summary_df"]
        n_unique = counts["n_passed"] - counts["n_unmapped"] - counts["n_ambiguous"]

        st.subheader("Mapping statistics")
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total reads", f"{counts['n_total']:,}")
        m2.metric("Passed filters", f"{counts['n_passed']:,}")
        pct_unique = 100 * n_unique / counts["n_passed"] if counts["n_passed"] else 0
        m3.metric("Uniquely assigned", f"{n_unique:,} ({pct_unique:.1f}%)")
        pct_ambig = 100 * counts["n_ambiguous"] / counts["n_passed"] if counts["n_passed"] else 0
        m4.metric("Ambiguous", f"{counts['n_ambiguous']:,} ({pct_ambig:.1f}%)")

        # Per-plasmid breakdown with percentages
        st.subheader("Per-plasmid summary")
        display_df = summary_df.copy()
        if n_unique > 0:
            display_df["% of unique"] = (
                display_df["assigned_reads"] / n_unique * 100
            ).round(1)
        n_assigned_total = counts["n_passed"] - counts["n_unmapped"]
        if n_assigned_total > 0:
            display_df["% of all assigned"] = (
                (display_df["assigned_reads"] + display_df["ambiguous_reads"])
                / n_assigned_total * 100
            ).round(1)
        st.dataframe(display_df, use_container_width=True)

        # ── Download buttons ─────────────────────────────────────────────
        col1, col2, col3 = st.columns(3)

        with open(result["output_path"], "rb") as f:
            col1.download_button(
                "Download assignments TSV",
                data=f.read(),
                file_name="plasmid_assignments.tsv",
                mime="text/tab-separated-values",
            )

        with open(result["summary_path"], "rb") as f:
            col2.download_button(
                "Download summary TSV",
                data=f.read(),
                file_name="plasmid_summary.tsv",
                mime="text/tab-separated-values",
            )

        if result["plots_path"] and os.path.exists(result["plots_path"]):
            with open(result["plots_path"], "rb") as f:
                pdf_bytes = f.read()
            col3.download_button(
                "Download PDF report",
                data=pdf_bytes,
                file_name="plasmid_alignment_report.pdf",
                mime="application/pdf",
            )

        # ── Inline plot images ───────────────────────────────────────────
        if result["plots_path"] and os.path.exists(result["plots_path"]):
            st.subheader("Visual report")
            page_images = _pdf_to_images(result["plots_path"])
            if page_images:
                for img_bytes in page_images:
                    st.image(img_bytes, use_container_width=True)
            else:
                st.caption(
                    "Install PyMuPDF (`pymupdf`) to display plots inline. "
                    "The PDF is still available for download above."
                )

        # Log
        with st.expander("Pipeline log", expanded=False):
            st.text("\n".join(result["log"]))

elif not can_run:
    st.info(
        "Upload a FASTQ file and select at least one plasmid reference "
        "(from the library or by uploading) to get started."
    )
