"""Backend helpers for the Streamlit ProbeDesign app."""

import io
import os
import re
import tempfile
import zipfile
from contextlib import redirect_stdout
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from probedesign.core import design_probes, ProbeDesignResult
from probedesign.output import (
    write_output_files,
    format_oligos_content,
    format_seq_content,
    format_hits_content,
)
from probedesign.hcr import design_hcr_probes, HCRDesignResult
from probedesign.hcr_output import (
    write_hcr_output_files,
    format_hcr_oligos,
    format_hcr_seq,
    format_hcr_hits,
)
from probedesign.provenance import failed_masks
from probedesign.idt import (
    build_plates,
    parse_oligos_text,
    plate_filename,
    write_plates_xlsx,
)
from probedesign.masking import (
    GENOME_INDEXES,
    PSEUDOGENE_INDEXES,
    find_bowtie,
    find_repeatmasker,
    run_repeatmasker,
)

XLSX_MIME = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
IDT_BUTTON_HELP = (
    "Excel file for IDT's 96-well plate upload, this probe set only. "
    "To put several probe sets on one order, use the IDT Plate tab."
)


# --- IDT plate export ---


def idt_plate_for_run(oligos_text: str) -> Tuple[str, bytes, List[str]]:
    """IDT plate order for one run: (file name, .xlsx bytes, warnings).

    Takes the run's _oligos.txt / _HCR_oligos.txt content, so the plate
    carries exactly the names and sequences of the downloaded oligos file.
    Raises ImportError if openpyxl is missing.
    """
    plates, warnings = build_plates([parse_oligos_text(oligos_text)])
    return plate_filename(plates), write_plates_xlsx(plates), warnings


# --- Name helpers ---


def mask_warnings(result) -> List[str]:
    """One message per requested bowtie mask that failed and was not applied."""
    if result is None:
        return []
    names = {"P": "Pseudogene", "B": "Genome"}
    return [
        f"{names.get(m.letter, m.letter)} masking failed and was not applied: {m.error}"
        for m in failed_masks(result)
    ]


def _label_run_info(result, input_name, repeatmask_mode, repeatmask_name) -> None:
    """Swap temp-file names in result.run_info for the names the user knows."""
    if not result.run_info:
        return
    if input_name:
        result.run_info["input_name"] = input_name
    if repeatmask_mode == "auto":
        # design saw RepeatMasker's temp output as a plain repeat-mask file
        result.run_info["repeatmask_mode"] = "auto"
        result.run_info["repeatmask_name"] = None
    elif repeatmask_name:
        result.run_info["repeatmask_name"] = repeatmask_name


def clean_name_for_prefix(name: str) -> str:
    """Convert a filename stem or FASTA header into a clean output name prefix.

    - Strips leading '>' (FASTA header marker) and surrounding whitespace
    - Replaces whitespace and non-alphanumeric/dot/dash characters with '_'
    - Squeezes consecutive underscores
    - Strips leading/trailing underscores
    - Falls back to 'probe' if the result would be empty

    Args:
        name: Raw string (filename stem or FASTA header text)

    Returns:
        Cleaned string suitable for use as an output name prefix.
    """
    name = name.lstrip(">").strip()
    name = re.sub(r'[^\w.\-]', '_', name)
    name = re.sub(r'_+', '_', name)
    name = name.strip('_')
    return name if name else "probe"


def extract_first_fasta_header(text: str) -> str:
    """Return the text of the first FASTA header line (everything after '>').

    Args:
        text: Raw FASTA text

    Returns:
        Header text (stripped), or empty string if no header found.
    """
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            return line[1:].strip()
    return ""


# --- FASTA validation and temp file helpers ---


def validate_fasta_text(text: str) -> Tuple[bool, str]:
    """Validate pasted FASTA text.

    Checks:
    - At least one header line starting with '>'
    - At least one non-empty sequence line after each header
    - Sequence lines contain only valid bases (ACGTNacgtn and whitespace)

    Returns:
        (is_valid, error_message). error_message is "" on success.
    """
    text = text.strip()
    if not text:
        return False, "Empty input."

    lines = text.split("\n")

    # Must start with a header
    if not lines[0].strip().startswith(">"):
        return False, "FASTA must start with a header line beginning with '>'."

    found_header = False
    found_seq_after_header = False

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if found_header and not found_seq_after_header:
                return False, "Header line has no sequence following it."
            found_header = True
            found_seq_after_header = False
        else:
            if not found_header:
                return False, "Sequence data found before any header line."
            # Check for valid characters (ACGTN only, case-insensitive)
            if not re.match(r'^[ACGTNacgtn]+$', line):
                invalid = re.findall(r'[^ACGTNacgtn]', line)
                return False, f"Invalid characters in sequence: {set(invalid)}"
            found_seq_after_header = True

    if not found_seq_after_header:
        return False, "No sequence data found after the last header."

    return True, ""


def save_fasta_to_temp(text: str) -> Path:
    """Write validated FASTA text to a named temp file.

    The file is NOT auto-deleted — the caller is responsible for cleanup.

    Returns:
        Path to the temp file.
    """
    f = tempfile.NamedTemporaryFile(
        mode='w', suffix='.fa', delete=False, prefix='probedesign_'
    )
    f.write(text)
    f.close()
    return Path(f.name)


def save_uploaded_file_to_temp(uploaded_file) -> Path:
    """Write a Streamlit UploadedFile to a named temp file.

    Returns:
        Path to the temp file.
    """
    content = uploaded_file.getvalue().decode('utf-8')
    f = tempfile.NamedTemporaryFile(
        mode='w', suffix='.fa', delete=False, prefix='probedesign_'
    )
    f.write(content)
    f.close()
    return Path(f.name)


def cleanup_temp_files(paths: List[Path]) -> None:
    """Delete temporary files. Silently ignores missing files."""
    for p in paths:
        try:
            if p.exists():
                p.unlink()
        except OSError:
            pass


# --- Prerequisite checks ---


@dataclass
class PrereqStatus:
    bowtie_ok: bool = False
    bowtie_path: Optional[str] = None
    repeatmasker_ok: bool = False
    repeatmasker_path: Optional[str] = None
    index_dir_exists: bool = False
    missing_indexes: List[str] = field(default_factory=list)


def check_prerequisites(
    pseudogene_mask: bool,
    genome_mask: bool,
    repeatmask_mode: str,
    species: str,
    index_dir: str,
) -> PrereqStatus:
    """Check that all required external tools and indexes are available.

    Args:
        pseudogene_mask: Whether pseudogene masking is enabled
        genome_mask: Whether genome masking is enabled
        repeatmask_mode: One of "none", "auto", "file"
        species: Species for masking databases
        index_dir: Path to bowtie indexes directory

    Returns:
        PrereqStatus with detailed availability info
    """
    status = PrereqStatus()

    # Check bowtie
    if pseudogene_mask or genome_mask:
        try:
            status.bowtie_path = find_bowtie()
            status.bowtie_ok = True
        except (FileNotFoundError, Exception):
            status.bowtie_ok = False

    # Check RepeatMasker
    if repeatmask_mode == "auto":
        try:
            status.repeatmasker_path = find_repeatmasker()
            status.repeatmasker_ok = True
        except (FileNotFoundError, Exception):
            status.repeatmasker_ok = False

    # Check index directory and files
    # Bowtie 1.3+ can read both .ebwt (bowtie1) and .bt2 (bowtie2) indexes
    if pseudogene_mask or genome_mask:
        if index_dir and Path(index_dir).is_dir():
            status.index_dir_exists = True

            if pseudogene_mask:
                db_name = PSEUDOGENE_INDEXES.get(species.lower(), "")
                if db_name:
                    d = Path(index_dir)
                    if not (d / f"{db_name}.1.ebwt").exists() and not (d / f"{db_name}.1.bt2").exists():
                        status.missing_indexes.append(
                            f"Pseudogene index not found: {db_name} in {index_dir}"
                        )

            if genome_mask:
                db_name = GENOME_INDEXES.get(species.lower(), "")
                if db_name:
                    d = Path(index_dir)
                    if not (d / f"{db_name}.1.ebwt").exists() and not (d / f"{db_name}.1.bt2").exists():
                        status.missing_indexes.append(
                            f"Genome index not found: {db_name} in {index_dir}"
                        )
        else:
            status.index_dir_exists = False
            if index_dir:
                status.missing_indexes.append(
                    f"Index directory does not exist: {index_dir}"
                )
            else:
                status.missing_indexes.append(
                    "No index directory specified. Required for pseudogene/genome masking."
                )

    return status


# --- Design runner ---


@dataclass
class DesignRunResult:
    result: Optional[ProbeDesignResult] = None
    stdout_log: str = ""
    error: Optional[str] = None


def run_design(
    input_path: str,
    n_probes: int = 48,
    oligo_length: int = 20,
    spacer_length: int = 3,
    target_gibbs: float = -23.0,
    allowable_gibbs: Tuple[float, float] = (-26.0, -20.0),
    output_name: Optional[str] = None,
    species: str = "human",
    pseudogene_mask: bool = False,
    genome_mask: bool = False,
    index_dir: Optional[str] = None,
    repeatmask_mode: str = "none",
    repeatmask_file: Optional[str] = None,
    save_bowtie_raw: bool = False,
    mixed_lengths: Optional[Tuple[int, int]] = None,
    hp_threshold: int = 5,
    di_threshold: int = 3,
    input_name: Optional[str] = None,
    repeatmask_name: Optional[str] = None,
) -> DesignRunResult:
    """Run probe design, replicating the CLI orchestration logic.

    1. If repeatmask_mode == "auto", calls run_repeatmasker() first.
    2. Calls design_probes() with the resolved repeatmask_file.
    3. Captures stdout (design_probes prints progress messages).
    4. Catches all exceptions and returns them as error string.

    input_name / repeatmask_name are the names the user knows (upload names);
    they replace the temp-file names in the result's run_info header.

    Returns:
        DesignRunResult with either result or error populated.
    """
    stdout_capture = io.StringIO()
    try:
        actual_repeatmask_file = repeatmask_file

        # Step 1: RepeatMasker (mirrors cli.py lines 144-162)
        if repeatmask_mode == "auto":
            with redirect_stdout(stdout_capture):
                print(f"Running RepeatMasker on {input_path}...")
            masked_file = run_repeatmasker(Path(input_path), species=species)
            actual_repeatmask_file = str(masked_file)
            with redirect_stdout(stdout_capture):
                print(f"RepeatMasker output: {masked_file}")

        # Step 2: design with stdout capture
        with redirect_stdout(stdout_capture):
            result = design_probes(
                input_file=input_path,
                n_probes=n_probes,
                oligo_length=oligo_length,
                spacer_length=spacer_length,
                target_gibbs=target_gibbs,
                allowable_gibbs=allowable_gibbs,
                output_name=output_name,
                species=species,
                pseudogene_mask=pseudogene_mask,
                genome_mask=genome_mask,
                index_dir=index_dir,
                repeatmask_file=actual_repeatmask_file,
                save_bowtie_raw=save_bowtie_raw,
                mixed_lengths=mixed_lengths,
                hp_threshold=hp_threshold,
                di_threshold=di_threshold,
            )
        _label_run_info(result, input_name, repeatmask_mode, repeatmask_name)

        if not result.probes:
            return DesignRunResult(
                result=result,
                stdout_log=stdout_capture.getvalue(),
                error="No valid probes found for this sequence. "
                      "Try adjusting Gibbs FE range or oligo length.",
            )

        return DesignRunResult(
            result=result,
            stdout_log=stdout_capture.getvalue(),
            error=None,
        )

    except Exception as e:
        return DesignRunResult(
            result=None,
            stdout_log=stdout_capture.getvalue(),
            error=f"{type(e).__name__}: {e}",
        )


# --- Batch runner ---


@dataclass
class BatchResult:
    filename: str
    n_probes_found: int = 0
    score: float = float('inf')
    status: str = "pending"  # "success", "error", "no_probes"
    error_msg: Optional[str] = None
    result: Optional[ProbeDesignResult] = None


def discover_fasta_files(directory: str) -> List[Path]:
    """Find all .fa/.fasta files in a directory (non-recursive).

    Returns:
        Sorted list of FASTA file paths.
    """
    d = Path(directory)
    files = list(d.glob("*.fa")) + list(d.glob("*.fasta"))
    return sorted(set(files), key=lambda p: p.name)


def run_batch(
    fasta_paths: List[Path],
    output_dir: Optional[str] = None,
    params: dict = None,
    progress_callback: Optional[Callable] = None,
    name_overrides: Optional[Dict[Path, str]] = None,
    input_names: Optional[Dict[Path, str]] = None,
) -> List[BatchResult]:
    """Run probe design on a list of FASTA files.

    Optionally writes output files to output_dir for each successful run.
    If output_dir is None, results are only held in memory (retrieve via
    package_batch_results_zip).

    Args:
        fasta_paths: List of FASTA file paths to process
        output_dir: Optional server-side directory for output files. When None,
            no files are written to disk.
        params: Dict of design parameters (passed to run_design)
        progress_callback: Optional callable(current_index, total, filename)
        name_overrides: Optional dict mapping Path -> output name
        input_names: Optional dict mapping Path -> upload name (for run_info)

    Returns:
        List of BatchResult for summary table.
    """
    if params is None:
        params = {}

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    batch_results = []

    for i, fasta_path in enumerate(fasta_paths):
        resolved_output_name = (name_overrides or {}).get(fasta_path, fasta_path.stem)
        fname = resolved_output_name

        if progress_callback:
            progress_callback(i, len(fasta_paths), fname)
        run_result = run_design(
            input_path=str(fasta_path),
            n_probes=params.get("n_probes", 48),
            oligo_length=params.get("oligo_length", 20),
            spacer_length=params.get("spacer_length", 3),
            target_gibbs=params.get("target_gibbs", -23.0),
            allowable_gibbs=params.get("allowable_gibbs", (-26.0, -20.0)),
            output_name=resolved_output_name,
            species=params.get("species", "human"),
            pseudogene_mask=params.get("pseudogene_mask", False),
            genome_mask=params.get("genome_mask", False),
            index_dir=params.get("index_dir"),
            repeatmask_mode=params.get("repeatmask_mode", "none"),
            repeatmask_file=params.get("repeatmask_file"),
            save_bowtie_raw=params.get("save_bowtie_raw", False),
            mixed_lengths=params.get("mixed_lengths"),
            hp_threshold=params.get("hp_threshold", 5),
            di_threshold=params.get("di_threshold", 3),
            input_name=(input_names or {}).get(fasta_path),
            repeatmask_name=params.get("repeatmask_name"),
        )

        if run_result.result is not None and run_result.result.probes:
            if output_dir:
                write_output_files(
                    run_result.result,
                    resolved_output_name,
                    mask_seqs=run_result.result.mask_strings,
                    output_dir=output_dir,
                )
            batch_results.append(BatchResult(
                filename=fname,
                n_probes_found=len(run_result.result.probes),
                score=run_result.result.score,
                status="success",
                error_msg=None,
                result=run_result.result,
            ))
        else:
            batch_results.append(BatchResult(
                filename=fname,
                n_probes_found=0,
                score=float('inf'),
                status="error" if run_result.error else "no_probes",
                error_msg=run_result.error,
                result=run_result.result,
            ))

    if progress_callback:
        progress_callback(len(fasta_paths), len(fasta_paths), "Done")

    return batch_results


def format_batch_summary(results: List[BatchResult]) -> str:
    """Format batch results as a TSV summary string.

    Columns: Filename, Probes_Found, Score, Status, Error, Warnings

    Args:
        results: List of BatchResult from run_batch

    Returns:
        TSV-formatted summary string
    """
    lines = ["Filename\tProbes_Found\tScore\tStatus\tError\tWarnings"]
    for br in results:
        score_str = f"{br.score:.4f}" if br.score != float('inf') else "N/A"
        error_str = br.error_msg or ""
        warnings = "; ".join(mask_warnings(br.result))
        lines.append(f"{br.filename}\t{br.n_probes_found}\t{score_str}\t{br.status}\t{error_str}\t{warnings}")
    return "\n".join(lines) + "\n"


def write_batch_summary(results: List[BatchResult], output_dir: str) -> str:
    """Write batch summary TSV to the output directory.

    Args:
        results: List of BatchResult from run_batch
        output_dir: Directory to write the summary file

    Returns:
        Path to the written summary file
    """
    summary = format_batch_summary(results)
    filepath = os.path.join(output_dir, "batch_summary.tsv")
    with open(filepath, 'w') as f:
        f.write(summary)
    return filepath


def package_batch_results_zip(results: List[BatchResult]) -> bytes:
    """Package all batch results into a ZIP archive (in memory).

    The ZIP contains per-gene output files generated from in-memory data,
    plus a batch_summary.tsv. Only successful runs (with probes) are included.

    ZIP structure:
        {prefix}_oligos.txt
        {prefix}_seq.txt
        {prefix}_{bowtie_hit_suffix}   (for each hit file, if masking was used)
        {prefix}_bowtie_pseudogene_raw.txt  (if raw bowtie output saved)
        {prefix}_bowtie_genome_raw.txt      (if raw bowtie output saved)
        batch_summary.tsv

    Args:
        results: List of BatchResult from run_batch

    Returns:
        ZIP file contents as bytes.
    """
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zf:
        for br in results:
            if br.result is None or not br.result.probes:
                continue
            prefix = br.filename

            zf.writestr(f"{prefix}_oligos.txt", format_oligos_content(br.result))
            zf.writestr(
                f"{prefix}_seq.txt",
                format_seq_content(br.result, mask_seqs=br.result.mask_strings),
            )

            for suffix, content in format_hits_content(br.result).items():
                zf.writestr(f"{prefix}_{suffix}", content)

            if br.result.bowtie_pseudogene_raw:
                zf.writestr(
                    f"{prefix}_bowtie_pseudogene_raw.txt",
                    br.result.bowtie_pseudogene_raw,
                )
            if br.result.bowtie_genome_raw:
                zf.writestr(
                    f"{prefix}_bowtie_genome_raw.txt",
                    br.result.bowtie_genome_raw,
                )

        zf.writestr("batch_summary.tsv", format_batch_summary(results))

    return buf.getvalue()


# --- HCR Design runner ---


@dataclass
class HCRDesignRunResult:
    result: Optional[HCRDesignResult] = None
    stdout_log: str = ""
    error: Optional[str] = None


def run_hcr_design(
    input_path: str,
    n_pairs: int = 30,
    amplifier: str = "B1",
    target_gibbs: float = -31.0,
    strict_range: Tuple[float, float] = (-35.0, -27.0),
    asymmetric_gibbs: bool = False,
    lenient_gibbs_min: float = -42.0,
    asymmetric_bowtie: bool = False,
    pair_spacing: int = 3,
    species: str = "human",
    pseudogene_mask: bool = False,
    genome_mask: bool = False,
    index_dir: Optional[str] = None,
    repeatmask_file: Optional[str] = None,
    repeatmask_mode: str = "none",
    hp_threshold: int = 5,
    di_threshold: int = 3,
    resolve_spacer: Optional[str] = None,
    output_name: Optional[str] = None,
    save_bowtie_raw: bool = False,
    input_name: Optional[str] = None,
    repeatmask_name: Optional[str] = None,
) -> HCRDesignRunResult:
    """Run HCR probe design with stdout capture and error handling.

    Analogous to run_design() for smFISH, including input_name / repeatmask_name.

    Returns:
        HCRDesignRunResult with either result or error populated.
    """
    stdout_capture = io.StringIO()
    try:
        actual_repeatmask_file = repeatmask_file

        if repeatmask_mode == "auto":
            with redirect_stdout(stdout_capture):
                print(f"Running RepeatMasker on {input_path}...")
            masked_file = run_repeatmasker(Path(input_path), species=species)
            actual_repeatmask_file = str(masked_file)
            with redirect_stdout(stdout_capture):
                print(f"RepeatMasker output: {masked_file}")

        with redirect_stdout(stdout_capture):
            result = design_hcr_probes(
                input_file=input_path,
                n_pairs=n_pairs,
                amplifier=amplifier,
                target_gibbs=target_gibbs,
                strict_range=strict_range,
                asymmetric_gibbs=asymmetric_gibbs,
                lenient_gibbs_min=lenient_gibbs_min,
                asymmetric_bowtie=asymmetric_bowtie,
                pair_spacing=pair_spacing,
                species=species,
                pseudogene_mask=pseudogene_mask,
                genome_mask=genome_mask,
                index_dir=index_dir,
                repeatmask_file=actual_repeatmask_file,
                hp_threshold=hp_threshold,
                di_threshold=di_threshold,
                resolve_spacer=resolve_spacer,
                output_name=output_name,
                save_bowtie_raw=save_bowtie_raw,
            )
        _label_run_info(result, input_name, repeatmask_mode, repeatmask_name)

        if not result.pairs:
            return HCRDesignRunResult(
                result=result,
                stdout_log=stdout_capture.getvalue(),
                error="No valid HCR probe pairs found. "
                      "Try adjusting ΔG range or enabling asymmetric mode.",
            )

        return HCRDesignRunResult(
            result=result,
            stdout_log=stdout_capture.getvalue(),
            error=None,
        )

    except Exception as e:
        return HCRDesignRunResult(
            result=None,
            stdout_log=stdout_capture.getvalue(),
            error=f"{type(e).__name__}: {e}",
        )


@dataclass
class HCRBatchResult:
    filename: str
    n_pairs_found: int = 0
    score: float = float('inf')
    status: str = "pending"
    error_msg: Optional[str] = None
    result: Optional[HCRDesignResult] = None


def run_hcr_batch(
    fasta_paths: List[Path],
    output_dir: Optional[str] = None,
    params: dict = None,
    progress_callback: Optional[Callable] = None,
    name_overrides: Optional[Dict[Path, str]] = None,
    input_names: Optional[Dict[Path, str]] = None,
) -> List[HCRBatchResult]:
    """Run HCR probe design on a list of FASTA files.

    Analogous to run_batch() for smFISH.

    Args:
        fasta_paths: List of FASTA file paths to process
        output_dir: Optional directory for output files
        params: Dict of HCR design parameters
        progress_callback: Optional callable(current_index, total, filename)
        name_overrides: Optional dict mapping Path -> output name
        input_names: Optional dict mapping Path -> upload name (for run_info)

    Returns:
        List of HCRBatchResult
    """
    if params is None:
        params = {}

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    batch_results = []

    for i, fasta_path in enumerate(fasta_paths):
        resolved_output_name = (name_overrides or {}).get(fasta_path, fasta_path.stem)
        fname = resolved_output_name

        if progress_callback:
            progress_callback(i, len(fasta_paths), fname)

        run_result = run_hcr_design(
            input_path=str(fasta_path),
            n_pairs=params.get("n_pairs", 30),
            amplifier=params.get("amplifier", "B1"),
            target_gibbs=params.get("target_gibbs", -31.0),
            strict_range=params.get("strict_range", (-35.0, -27.0)),
            asymmetric_gibbs=params.get("asymmetric_gibbs", False),
            lenient_gibbs_min=params.get("lenient_gibbs_min", -42.0),
            asymmetric_bowtie=params.get("asymmetric_bowtie", False),
            pair_spacing=params.get("pair_spacing", 3),
            species=params.get("species", "human"),
            pseudogene_mask=params.get("pseudogene_mask", False),
            genome_mask=params.get("genome_mask", False),
            index_dir=params.get("index_dir"),
            repeatmask_mode=params.get("repeatmask_mode", "none"),
            repeatmask_file=params.get("repeatmask_file"),
            hp_threshold=params.get("hp_threshold", 5),
            di_threshold=params.get("di_threshold", 3),
            resolve_spacer=params.get("resolve_spacer"),
            output_name=resolved_output_name,
            save_bowtie_raw=params.get("save_bowtie_raw", False),
            input_name=(input_names or {}).get(fasta_path),
            repeatmask_name=params.get("repeatmask_name"),
        )

        if run_result.result is not None and run_result.result.pairs:
            if output_dir:
                write_hcr_output_files(
                    run_result.result,
                    resolved_output_name,
                    output_dir=output_dir,
                )
            batch_results.append(HCRBatchResult(
                filename=fname,
                n_pairs_found=len(run_result.result.pairs),
                score=run_result.result.score,
                status="success",
                error_msg=None,
                result=run_result.result,
            ))
        else:
            batch_results.append(HCRBatchResult(
                filename=fname,
                n_pairs_found=0,
                score=float('inf'),
                status="error" if run_result.error else "no_pairs",
                error_msg=run_result.error,
                result=run_result.result,
            ))

    if progress_callback:
        progress_callback(len(fasta_paths), len(fasta_paths), "Done")

    return batch_results


def format_hcr_batch_summary(results: List[HCRBatchResult]) -> str:
    """Format HCR batch results as a TSV summary string."""
    lines = ["Filename\tPairs_Found\tScore\tStatus\tError\tWarnings"]
    for br in results:
        score_str = f"{br.score:.4f}" if br.score != float('inf') else "N/A"
        error_str = br.error_msg or ""
        warnings = "; ".join(mask_warnings(br.result))
        lines.append(f"{br.filename}\t{br.n_pairs_found}\t{score_str}\t{br.status}\t{error_str}\t{warnings}")
    return "\n".join(lines) + "\n"


def package_hcr_batch_results_zip(results: List[HCRBatchResult]) -> bytes:
    """Package all HCR batch results into a ZIP archive (in memory)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zf:
        for br in results:
            if br.result is None or not br.result.pairs:
                continue
            prefix = br.filename
            zf.writestr(f"{prefix}_HCR_oligos.txt", format_hcr_oligos(br.result))
            zf.writestr(f"{prefix}_HCR_seq.txt", format_hcr_seq(br.result))
            zf.writestr(f"{prefix}_HCR_hits.txt", format_hcr_hits(br.result))

        zf.writestr("hcr_batch_summary.tsv", format_hcr_batch_summary(results))

    return buf.getvalue()
