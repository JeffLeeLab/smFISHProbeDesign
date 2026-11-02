"""smFISH Probe Design page — refactored from original app.py."""

import os
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

# Ensure the probedesign package is importable from the repo root
_repo_root = Path(__file__).resolve().parent.parent.parent
if str(_repo_root / "src") not in sys.path:
    sys.path.insert(0, str(_repo_root / "src"))

# Add streamlit_app to path for utils import
_app_dir = Path(__file__).resolve().parent.parent
if str(_app_dir) not in sys.path:
    sys.path.insert(0, str(_app_dir))

from utils import (
    validate_fasta_text,
    save_fasta_to_temp,
    save_uploaded_file_to_temp,
    cleanup_temp_files,
    check_prerequisites,
    run_design,
    run_batch,
    discover_fasta_files,
    format_batch_summary,
    write_batch_summary,
    package_batch_results_zip,
    BatchResult,
    clean_name_for_prefix,
    extract_first_fasta_header,
    mask_warnings,
    idt_plate_for_run,
    IDT_BUTTON_HELP,
    XLSX_MIME,
)
from probedesign.output import (
    format_oligos_content,
    format_seq_content,
    format_hits_content,
)
from probedesign import __version__
from probedesign.provenance import PASTED_INPUT_NAME


# ── Session state initialization ─────────────────────────────────────────────

for key, default in {
    "smfish_single_result": None,
    "smfish_batch_results": None,
    "smfish_batch_output_dir": None,
    "smfish_temp_files": [],
}.items():
    if key not in st.session_state:
        st.session_state[key] = default


# ── Sidebar — parameters panel ───────────────────────────────────────────────

with st.sidebar:
    st.header("smFISH Probe Design")
    st.caption(f"ProbeDesign v{__version__}")

    mode = st.radio(
        "Mode",
        ["Single sequence", "Batch"],
        key="smfish_mode_radio",
        horizontal=True,
    )

    st.write("")

    # Basic parameters
    species = st.selectbox(
        "Species",
        ["human", "mouse", "elegans", "drosophila", "rat"],
        key="smfish_species",
    )

    n_probes = st.number_input(
        "Number of probes",
        min_value=1, value=48, step=1,
        key="smfish_n_probes",
    )
    oligo_mode = st.radio(
        "Oligo length mode",
        ["Fixed", "Mixed range"],
        key="smfish_oligo_mode",
        horizontal=True,
        help="Fixed: all probes same length. Mixed range: probes vary between min-max length.",
    )

    mixed_lengths = None
    if oligo_mode == "Fixed":
        oligo_length = st.number_input(
            "Oligo length (bp)",
            min_value=10, max_value=60, value=20, step=1,
            key="smfish_oligo_length",
        )
    else:
        ocol1, ocol2 = st.columns(2)
        with ocol1:
            mixed_min = st.number_input(
                "Min length (bp)",
                min_value=10, max_value=60, value=18, step=1,
                key="smfish_mixed_min",
            )
        with ocol2:
            mixed_max = st.number_input(
                "Max length (bp)",
                min_value=10, max_value=60, value=22, step=1,
                key="smfish_mixed_max",
            )
        if mixed_min > mixed_max:
            st.warning("Min length must not be greater than max length.")
        oligo_length = mixed_min
        mixed_lengths = (mixed_min, mixed_max)
    spacer_length = st.number_input(
        "Spacer length (bp)",
        min_value=0, value=3, step=1,
        key="smfish_spacer_length",
    )
    target_gibbs = st.number_input(
        "Target Gibbs FE (kcal/mol)",
        min_value=-100.0, max_value=0.0, value=-23.0, step=0.5,
        format="%.1f",
        key="smfish_target_gibbs",
    )

    col1, col2 = st.columns(2)
    with col1:
        gibbs_min = st.number_input(
            "Gibbs min",
            min_value=-200.0, max_value=0.0, value=-26.0, step=0.5,
            format="%.1f",
            key="smfish_gibbs_min",
        )
    with col2:
        gibbs_max = st.number_input(
            "Gibbs max",
            min_value=-200.0, max_value=0.0, value=-20.0, step=0.5,
            format="%.1f",
            key="smfish_gibbs_max",
        )

    st.divider()

    # Low-complexity filter
    st.markdown("Low-complexity filter")
    col_hp, col_di = st.columns(2)
    with col_hp:
        hp_threshold = st.number_input(
            "Homopolymer repeats",
            min_value=3, max_value=10, value=5, step=1,
            key="smfish_hp_threshold",
            help="Mask probes containing single-nucleotide runs of this length or longer (e.g. AAAAA at threshold 5)",
        )
    with col_di:
        di_threshold = st.number_input(
            "Dinucleotide repeats",
            min_value=3, max_value=10, value=3, step=1,
            key="smfish_di_threshold",
            help="Mask probes containing dinucleotide motifs repeated this many times or more (e.g. ATATAT at threshold 3)",
        )

    st.divider()

    # Bowtie masking
    st.markdown("Off-target filter")
    pseudogene_mask_on = st.checkbox(
        "Pseudogene mask", value=True, key="smfish_pseudogene_mask",
        help="Masking against pseudogenes"
    )
    genome_mask_on = st.checkbox(
        "Genome mask", value=True, key="smfish_genome_mask",
        help="Masking against the whole genome (more comprehensive but slower than pseudogene masking"
    )

    index_dir = None
    if pseudogene_mask_on or genome_mask_on:
        index_dir = st.text_input(
            "Bowtie index directory",
            value=str(_repo_root / "bowtie_indexes"),
            key="smfish_index_dir",
            help="Directory containing bowtie indexes for the selected species",
        )

    save_bowtie_raw = st.checkbox(
        "Save raw bowtie output",
        value=False,
        key="smfish_save_bowtie_raw",
        help="Include full bowtie alignment files in downloads",
    )

    # RepeatMasker
    with st.expander("Advanced RepeatMasker options", expanded=False):
        repeatmask_mode = st.radio(
            "Repeat masking",
            ["None", "Auto (RepeatMasker)", "Provide file"],
            key="smfish_repeatmask_mode",
        )

        repeatmask_upload = None
        if repeatmask_mode == "Provide file":
            repeatmask_upload = st.file_uploader(
                "Repeat mask FASTA",
                type=["fa", "fasta", "txt"],
                key="smfish_repeatmask_file_upload",
            )


# ── Helper functions ─────────────────────────────────────────────────────────

_RM_MODE_MAP = {
    "None": "none",
    "Auto (RepeatMasker)": "auto",
    "Provide file": "file",
}
rm_mode_key = _RM_MODE_MAP[repeatmask_mode]


def _run_prereq_checks() -> bool:
    prereq = check_prerequisites(
        pseudogene_mask=pseudogene_mask_on,
        genome_mask=genome_mask_on,
        repeatmask_mode=rm_mode_key,
        species=species,
        index_dir=index_dir or "",
    )
    ok = True
    if (pseudogene_mask_on or genome_mask_on) and not prereq.bowtie_ok:
        st.error("Bowtie is not installed or not found. Masking requires bowtie.")
        ok = False
    if rm_mode_key == "auto" and not prereq.repeatmasker_ok:
        st.error("RepeatMasker is not installed. Auto repeat masking requires RepeatMasker.")
        ok = False
    for msg in prereq.missing_indexes:
        st.warning(msg)
    return ok


def _resolve_repeatmask_file():
    if rm_mode_key == "file":
        if repeatmask_upload is None:
            st.warning("Please upload a repeat mask file.")
            return "__STOP__"
        tmp = save_uploaded_file_to_temp(repeatmask_upload)
        st.session_state["smfish_temp_files"].append(tmp)
        return str(tmp)
    return None


def _repeatmask_upload_name():
    """The uploaded repeat-mask file's own name, for the run header."""
    if rm_mode_key == "file" and repeatmask_upload is not None:
        return repeatmask_upload.name
    return None


def _show_single_results():
    run_result = st.session_state.get("smfish_single_result")
    if run_result is None:
        return

    if run_result.error:
        st.error(run_result.error)

    for warning in mask_warnings(run_result.result):
        st.warning(warning)

    if run_result.stdout_log.strip():
        with st.expander("Design log"):
            st.code(run_result.stdout_log, language=None)

    if run_result.result is not None and run_result.result.probes:
        result = run_result.result
        prefix = result.template_name

        mcol1, mcol2 = st.columns(2)
        mcol1.metric("Probes found", len(result.probes))
        mcol2.metric("Score", f"{result.score:.4f}")

        st.subheader("Designed Probes")
        probe_data = []
        probe_lengths = set(len(p.sequence) for p in result.probes)
        show_length = len(probe_lengths) > 1
        for p in result.probes:
            row = {"Index": p.index}
            if show_length:
                row["Length"] = len(p.sequence)
            row.update({
                "GC%": p.gc_percent,
                "Tm": p.tm,
                "Gibbs FE": p.gibbs_fe,
                "Sequence": p.sequence,
                "Name": p.name,
            })
            probe_data.append(row)
        df = pd.DataFrame(probe_data)
        st.dataframe(df, width='stretch', hide_index=True)

        seq_content = format_seq_content(result, mask_seqs=result.mask_strings)
        with st.expander("Sequence alignment"):
            st.code(seq_content, language=None)

        st.subheader("Download Files")
        oligos_content = format_oligos_content(result)

        dcol1, dcol2, dcol3 = st.columns(3)
        with dcol1:
            st.download_button(
                label="Download oligos.txt",
                data=oligos_content,
                file_name=f"{prefix}_oligos.txt",
                mime="text/plain",
                key="smfish_dl_oligos",
            )
        with dcol2:
            st.download_button(
                label="Download seq.txt",
                data=seq_content,
                file_name=f"{prefix}_seq.txt",
                mime="text/plain",
                key="smfish_dl_seq",
            )
        with dcol3:
            try:
                idt_name, idt_xlsx, idt_warnings = idt_plate_for_run(oligos_content)
            except ImportError as e:
                st.caption(str(e))
            else:
                st.download_button(
                    label="Download IDT plate (.xlsx)",
                    data=idt_xlsx,
                    file_name=idt_name,
                    mime=XLSX_MIME,
                    key="smfish_dl_idt",
                    help=IDT_BUTTON_HELP,
                )
                for w in idt_warnings:
                    st.caption(w)

        hits_files = format_hits_content(result)
        for suffix, content in hits_files.items():
            st.download_button(
                label=f"Download {suffix}",
                data=content,
                file_name=f"{prefix}_{suffix}",
                mime="text/plain",
                key=f"smfish_dl_{suffix}",
            )

        if result.bowtie_pseudogene_raw:
            st.download_button(
                label="Download bowtie_pseudogene_raw.txt",
                data=result.bowtie_pseudogene_raw,
                file_name=f"{prefix}_bowtie_pseudogene_raw.txt",
                mime="text/plain",
                key="smfish_dl_pseudo_raw",
            )
        if result.bowtie_genome_raw:
            st.download_button(
                label="Download bowtie_genome_raw.txt",
                data=result.bowtie_genome_raw,
                file_name=f"{prefix}_bowtie_genome_raw.txt",
                mime="text/plain",
                key="smfish_dl_genome_raw",
            )


# ══════════════════════════════════════════════════════════════════════════════
#  SINGLE MODE
# ══════════════════════════════════════════════════════════════════════════════

if mode == "Single sequence":
    st.title("ProbeDesign")
    st.caption("Design oligonucleotide probes for single molecule RNA FISH")

    output_name = st.text_input(
        "Output name prefix (optional, auto-derived from filename)",
        value="",
        key="smfish_single_output_name",
    )

    if gibbs_min >= gibbs_max:
        st.warning("Gibbs min should be less than Gibbs max.")

    upload_tab, paste_tab = st.tabs(["Upload FASTA", "Paste FASTA"])

    with upload_tab:
        uploaded_file = st.file_uploader(
            "Choose a FASTA file",
            type=["fa", "fasta", "txt"],
            key="smfish_single_upload",
        )
        run_from_upload = st.button("Design Probes", key="smfish_single_run_upload", type="primary")

    with paste_tab:
        pasted_text = st.text_area(
            "Paste FASTA-formatted sequence",
            height=200,
            placeholder=">gene_name\nATCGATCG...",
            key="smfish_single_paste",
        )
        run_from_paste = st.button("Design Probes", key="smfish_single_run_paste", type="primary")

    if run_from_upload or run_from_paste:
        cleanup_temp_files(st.session_state["smfish_temp_files"])
        st.session_state["smfish_temp_files"] = []
        st.session_state["smfish_single_result"] = None

        input_path = None
        auto_name = None
        input_name = None

        if run_from_upload:
            if uploaded_file is None:
                st.warning("Please upload a FASTA file.")
                st.stop()
            tmp_path = save_uploaded_file_to_temp(uploaded_file)
            st.session_state["smfish_temp_files"].append(tmp_path)
            input_path = str(tmp_path)
            auto_name = clean_name_for_prefix(Path(uploaded_file.name).stem)
            input_name = uploaded_file.name
        else:
            if not pasted_text.strip():
                st.warning("Please paste a FASTA sequence.")
                st.stop()
            valid, err = validate_fasta_text(pasted_text)
            if not valid:
                st.error(f"Invalid FASTA: {err}")
                st.stop()
            tmp_path = save_fasta_to_temp(pasted_text)
            st.session_state["smfish_temp_files"].append(tmp_path)
            input_path = str(tmp_path)
            header = extract_first_fasta_header(pasted_text)
            auto_name = clean_name_for_prefix(header) if header else "probe"
            input_name = PASTED_INPUT_NAME

        rm_file = _resolve_repeatmask_file()
        if rm_file == "__STOP__":
            st.stop()

        if not _run_prereq_checks():
            st.stop()

        resolved_name = output_name.strip() if output_name.strip() else auto_name

        with st.spinner("Running probe design..."):
            run_result = run_design(
                input_path=input_path,
                n_probes=n_probes,
                oligo_length=oligo_length,
                spacer_length=spacer_length,
                target_gibbs=target_gibbs,
                allowable_gibbs=(gibbs_min, gibbs_max),
                output_name=resolved_name,
                species=species,
                pseudogene_mask=pseudogene_mask_on,
                genome_mask=genome_mask_on,
                index_dir=index_dir,
                repeatmask_mode=rm_mode_key,
                repeatmask_file=rm_file,
                input_name=input_name,
                repeatmask_name=_repeatmask_upload_name(),
                save_bowtie_raw=save_bowtie_raw,
                mixed_lengths=mixed_lengths,
                hp_threshold=hp_threshold,
                di_threshold=di_threshold,
            )

        st.session_state["smfish_single_result"] = run_result

    _show_single_results()


# ══════════════════════════════════════════════════════════════════════════════
#  BATCH MODE
# ══════════════════════════════════════════════════════════════════════════════

else:
    st.title("ProbeDesign — Batch Mode")
    st.caption("Design probes for multiple sequences at once")

    input_method = st.radio(
        "Input method",
        ["Upload files", "Directory path"],
        key="smfish_batch_input_method",
        horizontal=True,
    )

    input_dir_path = None
    uploaded_files = None

    if input_method == "Directory path":
        input_dir_path = st.text_input(
            "Directory containing .fa/.fasta files",
            key="smfish_batch_input_dir",
        )
        st.caption("Input filesystem path. Use 'Upload files' if connecting remotely.")
    else:
        uploaded_files = st.file_uploader(
            "Upload FASTA files",
            type=["fa", "fasta", "txt"],
            accept_multiple_files=True,
            key="smfish_batch_upload",
        )

    with st.expander("Auto-save outputs to Local/Server directory (optional)", expanded=False):
        st.caption("Leave blank to download results as a ZIP file (recommended for remote access).")
        batch_output_dir = st.text_input(
            "Output directory path on server",
            value="",
            key="smfish_batch_output_dir_input",
        )

    st.write("")

    if gibbs_min >= gibbs_max:
        st.warning("Gibbs min should be less than Gibbs max.")

    if st.button("Run Batch", key="smfish_batch_run_btn", type="primary"):
        cleanup_temp_files(st.session_state["smfish_temp_files"])
        st.session_state["smfish_temp_files"] = []
        st.session_state["smfish_batch_results"] = None
        st.session_state["smfish_batch_output_dir"] = None

        fasta_paths = []
        batch_name_overrides = {}
        batch_input_names = {}
        if input_method == "Directory path":
            if not input_dir_path or not Path(input_dir_path).is_dir():
                st.error(f"Directory not found: {input_dir_path}")
                st.stop()
            fasta_paths = discover_fasta_files(input_dir_path)
            if not fasta_paths:
                st.warning("No .fa or .fasta files found in directory.")
                st.stop()
        else:
            if not uploaded_files:
                st.warning("Please upload at least one FASTA file.")
                st.stop()
            for uf in uploaded_files:
                tmp = save_uploaded_file_to_temp(uf)
                st.session_state["smfish_temp_files"].append(tmp)
                fasta_paths.append(tmp)
                batch_name_overrides[tmp] = clean_name_for_prefix(Path(uf.name).stem)
                batch_input_names[tmp] = uf.name

        rm_file = _resolve_repeatmask_file()
        if rm_file == "__STOP__":
            st.stop()

        if not _run_prereq_checks():
            st.stop()

        params = {
            "n_probes": n_probes,
            "oligo_length": oligo_length,
            "spacer_length": spacer_length,
            "target_gibbs": target_gibbs,
            "allowable_gibbs": (gibbs_min, gibbs_max),
            "species": species,
            "pseudogene_mask": pseudogene_mask_on,
            "genome_mask": genome_mask_on,
            "index_dir": index_dir,
            "repeatmask_mode": rm_mode_key,
            "repeatmask_file": rm_file,
            "repeatmask_name": _repeatmask_upload_name(),
            "save_bowtie_raw": save_bowtie_raw,
            "mixed_lengths": mixed_lengths,
            "hp_threshold": hp_threshold,
            "di_threshold": di_threshold,
        }

        progress_bar = st.progress(0, text="Starting batch...")

        def _progress_cb(current, total, filename):
            if total > 0:
                progress_bar.progress(
                    current / total,
                    text=f"Processing {filename} ({current + 1}/{total})..."
                    if current < total
                    else "Batch complete!",
                )

        batch_results = run_batch(
            fasta_paths=fasta_paths,
            output_dir=batch_output_dir.strip() or None,
            params=params,
            progress_callback=_progress_cb,
            name_overrides=batch_name_overrides if batch_name_overrides else None,
            input_names=batch_input_names or None,
        )

        if batch_output_dir.strip():
            write_batch_summary(batch_results, batch_output_dir)

        progress_bar.progress(1.0, text="Batch complete!")
        st.session_state["smfish_batch_results"] = batch_results
        st.session_state["smfish_batch_output_dir"] = batch_output_dir.strip() or None

    batch_results = st.session_state.get("smfish_batch_results")
    saved_output_dir = st.session_state.get("smfish_batch_output_dir")

    if batch_results is not None:
        st.subheader("Batch Summary")

        summary_data = []
        for br in batch_results:
            summary_data.append({
                "File": br.filename,
                "Probes": br.n_probes_found,
                "Score": f"{br.score:.4f}" if br.score != float('inf') else "N/A",
                "Status": br.status,
                "Error": br.error_msg or "",
                "Warnings": "; ".join(mask_warnings(br.result)),
            })
        df = pd.DataFrame(summary_data)
        masking_failed = [br.filename for br in batch_results if mask_warnings(br.result)]
        if masking_failed:
            st.warning(
                "Masking failed and was not applied for: "
                + ", ".join(masking_failed) + ". See the Warnings column."
            )
        st.dataframe(df, width='stretch', hide_index=True)

        zip_bytes = package_batch_results_zip(batch_results)
        dcol1, dcol2 = st.columns(2)
        with dcol1:
            st.download_button(
                "Download all results (.zip)",
                data=zip_bytes,
                file_name="probedesign_batch_results.zip",
                mime="application/zip",
                key="smfish_dl_batch_zip",
            )
        with dcol2:
            summary_text = format_batch_summary(batch_results)
            st.download_button(
                "Download batch_summary.tsv",
                data=summary_text,
                file_name="batch_summary.tsv",
                mime="text/tab-separated-values",
                key="smfish_dl_batch_summary",
            )

        if saved_output_dir:
            st.info(f"Output files also written to server: {saved_output_dir}")

        for idx, br in enumerate(batch_results):
            if br.result is not None and br.result.probes:
                with st.expander(f"{br.filename} — {br.n_probes_found} probes (score: {br.score:.4f})"):
                    br_lengths = set(len(p.sequence) for p in br.result.probes)
                    br_show_length = len(br_lengths) > 1
                    probe_data = []
                    for p in br.result.probes:
                        row = {"Index": p.index}
                        if br_show_length:
                            row["Length"] = len(p.sequence)
                        row.update({
                            "GC%": p.gc_percent,
                            "Tm": p.tm,
                            "Gibbs FE": p.gibbs_fe,
                            "Sequence": p.sequence,
                        })
                        probe_data.append(row)
                    st.dataframe(
                        pd.DataFrame(probe_data),
                        width='stretch',
                        hide_index=True,
                    )
            elif br.error_msg:
                with st.expander(f"{br.filename} — FAILED"):
                    st.error(br.error_msg)
