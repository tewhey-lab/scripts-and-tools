"""Streamlit web UI for plasmid_align.py."""

import os
import tempfile
import streamlit as st

from plasmid_align import run_pipeline

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
    fasta_files = st.file_uploader(
        "FASTA reference(s)", type=["fasta", "fa", "fna"],
        accept_multiple_files=True,
        help="One or more plasmid FASTA files. Each may contain multiple sequences.",
    )

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
        "Circular mode",
        help="Double each reference to capture reads spanning the "
             "linearisation junction.",
    )
    fill_n = st.checkbox(
        "Fill N regions",
        help="Fill N-masked regions in references using the most homologous "
             "sequence from another reference in the library. Improves "
             "alignment coverage through masked regions.",
    )
    threads = st.slider("Threads", min_value=1, max_value=4, value=2)

# ── Main area: run button + results ──────────────────────────────────────────
can_run = fastq_file is not None and len(fasta_files) > 0

if st.button("Run Alignment", disabled=not can_run, type="primary"):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write uploaded files to disk so mappy can access them
        fastq_path = os.path.join(tmpdir, fastq_file.name)
        with open(fastq_path, "wb") as f:
            f.write(fastq_file.getbuffer())

        fasta_paths = []
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

        # Log
        with st.expander("Pipeline log", expanded=False):
            st.text("\n".join(result["log"]))

        # Summary table
        st.subheader("Per-plasmid summary")
        st.dataframe(result["summary_df"], use_container_width=True)

        # Download buttons
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
                col3.download_button(
                    "Download PDF report",
                    data=f.read(),
                    file_name="plasmid_alignment_report.pdf",
                    mime="application/pdf",
                )

elif not can_run:
    st.info("Upload a FASTQ file and at least one FASTA reference to get started.")
