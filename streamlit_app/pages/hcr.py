"""HCR Split-Initiator Probe Design page."""

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
    run_hcr_design,
    run_hcr_batch,
    discover_fasta_files,
    format_hcr_batch_summary,
    package_hcr_batch_results_zip,
    HCRBatchResult,
    clean_name_for_prefix,
    extract_first_fasta_header,
)
from probedesign.hcr_output import (
    format_hcr_oligos,
    format_hcr_seq,
    format_hcr_hits,
)
from probedesign.hcr import VALID_AMPLIFIERS


# ── Session state initialization ─────────────────────────────────────────────

for key, default in {
    "hcr_single_result": None,
    "hcr_batch_results": None,
    "hcr_batch_output_dir": None,
    "hcr_temp_files": [],
}.items():
    if key not in st.session_state:
        st.session_state[key] = default


# ── Sidebar — parameters panel ───────────────────────────────────────────────

with st.sidebar:
    st.header("HCR Probe Design")

    mode = st.radio(
        "Mode",
        ["Single sequence", "Batch"],
        key="hcr_mode_radio",
        horizontal=True,
    )

    st.write("")

    # Basic parameters
    species = st.selectbox(
        "Species",
        ["human", "mouse", "elegans", "drosophila", "rat"],
        key="hcr_species",
    )

    amplifier = st.selectbox(
        "Amplifier",
        VALID_AMPLIFIERS,
        key="hcr_amplifier",
        help="HCR amplifier system. B1-B5 are standard; B7+ use WW (degenerate) spacers.",
    )

    n_pairs = st.number_input(
        "Number of pairs",
        min_value=1, value=30, step=1,
        key="hcr_n_pairs",
    )

    pair_spacing = st.number_input(
        "Pair spacing (nt)",
        min_value=0, value=2, step=1,
        key="hcr_pair_spacing",
        help="Minimum gap between consecutive 52-nt pair blocks",
    )

    st.divider()

    # Thermodynamic parameters
    st.markdown("Thermodynamic parameters")
    target_gibbs = st.number_input(
        "Target ΔG per half (kcal/mol)",
        min_value=-100.0, max_value=0.0, value=-31.0, step=0.5,
        format="%.1f",
        key="hcr_target_gibbs",
    )

    col1, col2 = st.columns(2)
    with col1:
        gibbs_min = st.number_input(
            "Strict ΔG min",
            min_value=-200.0, max_value=0.0, value=-35.0, step=0.5,
            format="%.1f",
            key="hcr_gibbs_min",
        )
    with col2:
        gibbs_max = st.number_input(
            "Strict ΔG max",
            min_value=-200.0, max_value=0.0, value=-27.0, step=0.5,
            format="%.1f",
            key="hcr_gibbs_max",
        )

    asymmetric_gibbs = st.checkbox(
        "Asymmetric ΔG leniency",
        value=False,
        key="hcr_asymmetric_gibbs",
        help="Allow one half to have relaxed ΔG range while the other stays strict",
    )

    lenient_gibbs_min = -42.0
    if asymmetric_gibbs:
        lenient_gibbs_min = st.number_input(
            "Lenient ΔG floor",
            min_value=-200.0, max_value=0.0, value=-42.0, step=0.5,
            format="%.1f",
            key="hcr_lenient_gibbs_min",
        )

    st.divider()

    # Low-complexity filter
    st.markdown("Low-complexity filter")
    col_hp, col_di = st.columns(2)
    with col_hp:
        hp_threshold = st.number_input(
            "Homopolymer repeats",
            min_value=3, max_value=10, value=5, step=1,
            key="hcr_hp_threshold",
            help="Mask probes containing single-nucleotide runs of this length or longer",
        )
    with col_di:
        di_threshold = st.number_input(
            "Dinucleotide repeats",
            min_value=3, max_value=10, value=3, step=1,
            key="hcr_di_threshold",
            help="Mask probes containing dinucleotide motifs repeated this many times or more",
        )

    st.divider()

    # Off-target filter
    st.markdown("Off-target filter")
    pseudogene_mask_on = st.checkbox(
        "Pseudogene mask", value=True, key="hcr_pseudogene_mask",
        help="Masking against pseudogenes"
    )
    genome_mask_on = st.checkbox(
        "Genome mask", value=True, key="hcr_genome_mask",
        help="Masking against the whole genome"
    )

    asymmetric_bowtie = False
    if pseudogene_mask_on or genome_mask_on:
        asymmetric_bowtie = st.checkbox(
            "Asymmetric bowtie",
            value=False,
            key="hcr_asymmetric_bowtie",
            help="Only require one half to be clean of off-target hits",
        )

    index_dir = None
    if pseudogene_mask_on or genome_mask_on:
        index_dir = st.text_input(
            "Bowtie index directory",
            value=str(_repo_root / "bowtie_indexes"),
            key="hcr_index_dir",
        )

    save_bowtie_raw = st.checkbox(
        "Save raw bowtie output",
        value=False,
        key="hcr_save_bowtie_raw",
    )

    # RepeatMasker
    with st.expander("Advanced RepeatMasker options", expanded=False):
        repeatmask_mode = st.radio(
            "Repeat masking",
            ["None", "Auto (RepeatMasker)", "Provide file"],
            key="hcr_repeatmask_mode",
        )
        repeatmask_upload = None
        if repeatmask_mode == "Provide file":
            repeatmask_upload = st.file_uploader(
                "Repeat mask FASTA",
                type=["fa", "fasta", "txt"],
                key="hcr_repeatmask_file_upload",
            )

    # Spacer resolution (advanced)
    with st.expander("Spacer options", expanded=False):
        resolve_spacer_input = st.text_input(
            "Resolve WW spacer (leave blank for IUPAC default)",
            value="",
            key="hcr_resolve_spacer",
            help="For amplifiers B7+, replace WW with explicit bases (e.g. AA, TT)",
        )
        resolve_spacer = resolve_spacer_input.strip().upper() if resolve_spacer_input.strip() else None


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
        st.error("RepeatMasker is not installed.")
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
        st.session_state["hcr_temp_files"].append(tmp)
        return str(tmp)
    return None


def _show_hcr_results():
    run_result = st.session_state.get("hcr_single_result")
    if run_result is None:
        return

    if run_result.error:
        st.error(run_result.error)

    if run_result.stdout_log.strip():
        with st.expander("Design log"):
            st.code(run_result.stdout_log, language=None)

    if run_result.result is not None and run_result.result.pairs:
        result = run_result.result
        prefix = result.template_name

        mcol1, mcol2, mcol3 = st.columns(3)
        mcol1.metric("Pairs found", len(result.pairs))
        mcol2.metric("Score", f"{result.score:.4f}")
        mcol3.metric("Mode", result.asymmetric_mode)

        # Probe pair table
        st.subheader("Designed Probe Pairs")
        pair_data = []
        for pair in result.pairs:
            p1_id = f"{2 * pair.pair_index - 1:02d}"
            p2_id = f"{2 * pair.pair_index:02d}"
            pair_data.append({
                "Pair": pair.pair_index,
                "Half": "P1",
                "Start": pair.right.position + 1,
                "GC%": pair.right.gc_percent,
                "Tm": pair.right.tm,
                "Gibbs": pair.right.gibbs_fe,
                "Strict": "yes" if pair.right.is_strict else "no",
                "Oligo": pair.right.oligo_seq,
                "Name": f"{result.template_name}_HCR{result.amplifier}_{p1_id}",
            })
            pair_data.append({
                "Pair": pair.pair_index,
                "Half": "P2",
                "Start": pair.left.position + 1,
                "GC%": pair.left.gc_percent,
                "Tm": pair.left.tm,
                "Gibbs": pair.left.gibbs_fe,
                "Strict": "yes" if pair.left.is_strict else "no",
                "Oligo": pair.left.oligo_seq,
                "Name": f"{result.template_name}_HCR{result.amplifier}_{p2_id}",
            })
        df = pd.DataFrame(pair_data)
        st.dataframe(df, width='stretch', hide_index=True)

        # Sequence alignment
        seq_content = format_hcr_seq(result)
        with st.expander("Sequence alignment"):
            st.code(seq_content, language=None)

        # Download buttons
        st.subheader("Download Files")
        oligos_content = format_hcr_oligos(result)
        hits_content = format_hcr_hits(result)

        dcol1, dcol2, dcol3 = st.columns(3)
        with dcol1:
            st.download_button(
                label="Download HCR_oligos.txt",
                data=oligos_content,
                file_name=f"{prefix}_HCR_oligos.txt",
                mime="text/plain",
                key="hcr_dl_oligos",
            )
        with dcol2:
            st.download_button(
                label="Download HCR_seq.txt",
                data=seq_content,
                file_name=f"{prefix}_HCR_seq.txt",
                mime="text/plain",
                key="hcr_dl_seq",
            )
        with dcol3:
            st.download_button(
                label="Download HCR_hits.txt",
                data=hits_content,
                file_name=f"{prefix}_HCR_hits.txt",
                mime="text/plain",
                key="hcr_dl_hits",
            )


# ══════════════════════════════════════════════════════════════════════════════
#  SINGLE MODE
# ══════════════════════════════════════════════════════════════════════════════

if mode == "Single sequence":
    st.title("HCR ProbeDesign")
    st.caption("Design HCR v3 split-initiator probe pairs for RNA FISH")

    output_name = st.text_input(
        "Output name prefix (optional, auto-derived from filename)",
        value="",
        key="hcr_single_output_name",
    )

    if gibbs_min >= gibbs_max:
        st.warning("Strict ΔG min should be less than max.")

    upload_tab, paste_tab = st.tabs(["Upload FASTA", "Paste FASTA"])

    with upload_tab:
        uploaded_file = st.file_uploader(
            "Choose a FASTA file",
            type=["fa", "fasta", "txt"],
            key="hcr_single_upload",
        )
        run_from_upload = st.button("Design HCR Probes", key="hcr_single_run_upload", type="primary")

    with paste_tab:
        pasted_text = st.text_area(
            "Paste FASTA-formatted sequence",
            height=200,
            placeholder=">gene_name\nATCGATCG...",
            key="hcr_single_paste",
        )
        run_from_paste = st.button("Design HCR Probes", key="hcr_single_run_paste", type="primary")

    if run_from_upload or run_from_paste:
        cleanup_temp_files(st.session_state["hcr_temp_files"])
        st.session_state["hcr_temp_files"] = []
        st.session_state["hcr_single_result"] = None

        input_path = None
        auto_name = None

        if run_from_upload:
            if uploaded_file is None:
                st.warning("Please upload a FASTA file.")
                st.stop()
            tmp_path = save_uploaded_file_to_temp(uploaded_file)
            st.session_state["hcr_temp_files"].append(tmp_path)
            input_path = str(tmp_path)
            auto_name = clean_name_for_prefix(Path(uploaded_file.name).stem)
        else:
            if not pasted_text.strip():
                st.warning("Please paste a FASTA sequence.")
                st.stop()
            valid, err = validate_fasta_text(pasted_text)
            if not valid:
                st.error(f"Invalid FASTA: {err}")
                st.stop()
            tmp_path = save_fasta_to_temp(pasted_text)
            st.session_state["hcr_temp_files"].append(tmp_path)
            input_path = str(tmp_path)
            header = extract_first_fasta_header(pasted_text)
            auto_name = clean_name_for_prefix(header) if header else "probe"

        rm_file = _resolve_repeatmask_file()
        if rm_file == "__STOP__":
            st.stop()

        if not _run_prereq_checks():
            st.stop()

        resolved_name = output_name.strip() if output_name.strip() else auto_name

        with st.spinner("Designing HCR probe pairs..."):
            run_result = run_hcr_design(
                input_path=input_path,
                n_pairs=n_pairs,
                amplifier=amplifier,
                target_gibbs=target_gibbs,
                strict_range=(gibbs_min, gibbs_max),
                asymmetric_gibbs=asymmetric_gibbs,
                lenient_gibbs_min=lenient_gibbs_min,
                asymmetric_bowtie=asymmetric_bowtie,
                pair_spacing=pair_spacing,
                species=species,
                pseudogene_mask=pseudogene_mask_on,
                genome_mask=genome_mask_on,
                index_dir=index_dir,
                repeatmask_file=rm_file,
                repeatmask_mode=rm_mode_key,
                hp_threshold=hp_threshold,
                di_threshold=di_threshold,
                resolve_spacer=resolve_spacer,
                output_name=resolved_name,
                save_bowtie_raw=save_bowtie_raw,
            )

        st.session_state["hcr_single_result"] = run_result

    _show_hcr_results()


# ══════════════════════════════════════════════════════════════════════════════
#  BATCH MODE
# ══════════════════════════════════════════════════════════════════════════════

else:
    st.title("HCR ProbeDesign — Batch Mode")
    st.caption("Design HCR probe pairs for multiple sequences at once")

    input_method = st.radio(
        "Input method",
        ["Upload files", "Directory path"],
        key="hcr_batch_input_method",
        horizontal=True,
    )

    input_dir_path = None
    uploaded_files = None

    if input_method == "Directory path":
        input_dir_path = st.text_input(
            "Directory containing .fa/.fasta files",
            key="hcr_batch_input_dir",
        )
    else:
        uploaded_files = st.file_uploader(
            "Upload FASTA files",
            type=["fa", "fasta", "txt"],
            accept_multiple_files=True,
            key="hcr_batch_upload",
        )

    with st.expander("Auto-save outputs to Local/Server directory (optional)", expanded=False):
        st.caption("Leave blank to download results as a ZIP file.")
        batch_output_dir = st.text_input(
            "Output directory path on server",
            value="",
            key="hcr_batch_output_dir_input",
        )

    st.write("")

    if gibbs_min >= gibbs_max:
        st.warning("Strict ΔG min should be less than max.")

    if st.button("Run HCR Batch", key="hcr_batch_run_btn", type="primary"):
        cleanup_temp_files(st.session_state["hcr_temp_files"])
        st.session_state["hcr_temp_files"] = []
        st.session_state["hcr_batch_results"] = None
        st.session_state["hcr_batch_output_dir"] = None

        fasta_paths = []
        batch_name_overrides = {}
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
                st.session_state["hcr_temp_files"].append(tmp)
                fasta_paths.append(tmp)
                batch_name_overrides[tmp] = clean_name_for_prefix(Path(uf.name).stem)

        rm_file = _resolve_repeatmask_file()
        if rm_file == "__STOP__":
            st.stop()

        if not _run_prereq_checks():
            st.stop()

        params = {
            "n_pairs": n_pairs,
            "amplifier": amplifier,
            "target_gibbs": target_gibbs,
            "strict_range": (gibbs_min, gibbs_max),
            "asymmetric_gibbs": asymmetric_gibbs,
            "lenient_gibbs_min": lenient_gibbs_min,
            "asymmetric_bowtie": asymmetric_bowtie,
            "pair_spacing": pair_spacing,
            "species": species,
            "pseudogene_mask": pseudogene_mask_on,
            "genome_mask": genome_mask_on,
            "index_dir": index_dir,
            "repeatmask_mode": rm_mode_key,
            "repeatmask_file": rm_file,
            "hp_threshold": hp_threshold,
            "di_threshold": di_threshold,
            "resolve_spacer": resolve_spacer,
            "save_bowtie_raw": save_bowtie_raw,
        }

        progress_bar = st.progress(0, text="Starting HCR batch...")

        def _progress_cb(current, total, filename):
            if total > 0:
                progress_bar.progress(
                    current / total,
                    text=f"Processing {filename} ({current + 1}/{total})..."
                    if current < total
                    else "Batch complete!",
                )

        batch_results = run_hcr_batch(
            fasta_paths=fasta_paths,
            output_dir=batch_output_dir.strip() or None,
            params=params,
            progress_callback=_progress_cb,
            name_overrides=batch_name_overrides if batch_name_overrides else None,
        )

        progress_bar.progress(1.0, text="Batch complete!")
        st.session_state["hcr_batch_results"] = batch_results
        st.session_state["hcr_batch_output_dir"] = batch_output_dir.strip() or None

    batch_results = st.session_state.get("hcr_batch_results")
    saved_output_dir = st.session_state.get("hcr_batch_output_dir")

    if batch_results is not None:
        st.subheader("HCR Batch Summary")

        summary_data = []
        for br in batch_results:
            summary_data.append({
                "File": br.filename,
                "Pairs": br.n_pairs_found,
                "Score": f"{br.score:.4f}" if br.score != float('inf') else "N/A",
                "Status": br.status,
                "Error": br.error_msg or "",
            })
        df = pd.DataFrame(summary_data)
        st.dataframe(df, width='stretch', hide_index=True)

        zip_bytes = package_hcr_batch_results_zip(batch_results)
        dcol1, dcol2 = st.columns(2)
        with dcol1:
            st.download_button(
                "Download all results (.zip)",
                data=zip_bytes,
                file_name="hcr_batch_results.zip",
                mime="application/zip",
                key="hcr_dl_batch_zip",
            )
        with dcol2:
            summary_text = format_hcr_batch_summary(batch_results)
            st.download_button(
                "Download hcr_batch_summary.tsv",
                data=summary_text,
                file_name="hcr_batch_summary.tsv",
                mime="text/tab-separated-values",
                key="hcr_dl_batch_summary",
            )

        if saved_output_dir:
            st.info(f"Output files also written to server: {saved_output_dir}")

        for idx, br in enumerate(batch_results):
            if br.result is not None and br.result.pairs:
                with st.expander(f"{br.filename} — {br.n_pairs_found} pairs (score: {br.score:.4f})"):
                    pair_data = []
                    for pair in br.result.pairs:
                        pair_data.append({
                            "Pair": pair.pair_index,
                            "Half": "P1",
                            "GC%": pair.right.gc_percent,
                            "Tm": pair.right.tm,
                            "Gibbs": pair.right.gibbs_fe,
                            "Strict": "yes" if pair.right.is_strict else "no",
                            "Oligo": pair.right.oligo_seq,
                        })
                        pair_data.append({
                            "Pair": pair.pair_index,
                            "Half": "P2",
                            "GC%": pair.left.gc_percent,
                            "Tm": pair.left.tm,
                            "Gibbs": pair.left.gibbs_fe,
                            "Strict": "yes" if pair.left.is_strict else "no",
                            "Oligo": pair.left.oligo_seq,
                        })
                    st.dataframe(
                        pd.DataFrame(pair_data),
                        width='stretch',
                        hide_index=True,
                    )
            elif br.error_msg:
                with st.expander(f"{br.filename} — FAILED"):
                    st.error(br.error_msg)
