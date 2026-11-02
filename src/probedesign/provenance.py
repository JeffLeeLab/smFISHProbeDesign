"""Run provenance: what a design used, written as a '#' header in _seq.txt files.

design_probes() and design_hcr_probes() fill result.run_info with the effective
parameters, input details, software versions and one MaskRecord per bowtie
mask. format_header() turns that into '#' lines that sit above the sequence
overlay, and build_command() gives the equivalent CLI command.

The header names files by basename (or upload name) only, so it never shows
temp paths or the server's directory layout.
"""

import functools
import math
import shlex
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from . import __version__
from .masking import (
    DEFAULT_INDEX_DIR,
    GENOME_ASSEMBLIES,
    GENOME_INDEXES,
    PSEUDOGENE_INDEXES,
    bowtie_version,
    index_fingerprint,
)

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Input name the app uses for a pasted sequence (also used in the command line)
PASTED_INPUT_NAME = "pasted_sequence.fa"

@dataclass
class MaskRecord:
    """Outcome of one bowtie mask (pseudogene "P" or genome "B")."""
    letter: str
    status: str  # "applied" or "failed"
    index: str
    assembly: Optional[str] = None
    fingerprint: Optional[str] = None
    error: Optional[str] = None


# ── Collecting ───────────────────────────────────────────────────────────────

def _git_revision() -> Optional[str]:
    """Short commit hash of the checkout, with '-dirty' for uncommitted changes."""
    if not (_REPO_ROOT / ".git").exists():
        return None
    try:
        def git(*args):
            return subprocess.run(
                ["git", "-C", str(_REPO_ROOT), *args],
                capture_output=True, text=True, check=True, timeout=2,
            ).stdout.strip()
        rev = git("rev-parse", "--short", "HEAD")
        dirty = git("status", "--porcelain", "--untracked-files=no")
    except (OSError, subprocess.SubprocessError):
        return None
    return rev + ("-dirty" if dirty else "")


@functools.lru_cache(maxsize=None)
def software_info() -> dict:
    """Package version and git revision (cached: the loaded code can't change)."""
    return {"version": __version__, "git": _git_revision()}


def new_run_info(
    tool: str,
    input_file: str,
    n_entries: int,
    seq: str,
    params: dict,
    repeatmask_file: Optional[str],
    index_dir: Optional[str],
) -> dict:
    """Start the run_info dict for a design; masks are appended as they run."""
    if repeatmask_file:
        repeatmask_mode = "file"
    elif "n" in seq:
        repeatmask_mode = "input_n"
    else:
        repeatmask_mode = "none"
    uses_bowtie = params.get("pseudogene_mask") or params.get("genome_mask")
    return {
        "tool": tool,
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        **software_info(),
        "input_name": Path(input_file).name,
        "n_entries": n_entries,
        "length_nt": len(seq) - seq.count(">"),
        "params": params,
        "repeatmask_mode": repeatmask_mode,
        "repeatmask_name": Path(repeatmask_file).name if repeatmask_file else None,
        "index_dir": str(index_dir) if index_dir else None,
        "bowtie": bowtie_version() if uses_bowtie else None,
        "masks": [],
    }


def bowtie_mask_record(
    letter: str,
    species: str,
    index_dir: Optional[Path],
    error: Optional[BaseException] = None,
) -> MaskRecord:
    """Describe a pseudogene ("P") or genome ("B") mask that ran or failed."""
    species = species.lower()
    table = PSEUDOGENE_INDEXES if letter == "P" else GENOME_INDEXES
    index = table.get(species, "?")
    message = None
    if error is not None:
        lines = str(error).strip().splitlines()
        message = (lines[0] if lines else type(error).__name__)[:200]
    return MaskRecord(
        letter=letter,
        status="failed" if error is not None else "applied",
        index=index,
        assembly=GENOME_ASSEMBLIES.get(species) if letter == "B" else None,
        fingerprint=index_fingerprint(index_dir, index),
        error=message,
    )


def failed_masks(result) -> List[MaskRecord]:
    """Bowtie masks that were requested but not applied."""
    return [m for m in result.run_info.get("masks", []) if m.status == "failed"]


def display_index_dir(index_dir: Optional[str]) -> str:
    """'default' for <repo>/bowtie_indexes (however it was given), else the path."""
    if not index_dir:
        return "default"
    try:
        if Path(index_dir).resolve() == DEFAULT_INDEX_DIR.resolve():
            return "default"
    except OSError:
        pass
    return str(index_dir)


# ── Command line ─────────────────────────────────────────────────────────────

def _fmt(value) -> str:
    if isinstance(value, float) and value.is_integer():
        return f"{value:.1f}"
    return str(value)


def _option(flag: str, value) -> str:
    text = _fmt(value)
    # '=' keeps negative numbers from looking like options
    if text.startswith("-"):
        return f"{flag}={shlex.quote(text)}"
    return f"{flag} {shlex.quote(text)}"


def _range(pair) -> str:
    return f"{_fmt(pair[0])},{_fmt(pair[1])}"


def build_command(result) -> str:
    """The `probedesign design` / `hcrdesign` command that reproduces this run.

    Uses long flags only and lists only options that differ from the CLI
    defaults. File arguments are the names recorded in run_info, so the
    command runs from the directory that holds those files.
    """
    info = result.run_info
    p = info["params"]
    hcr = info["tool"] == "HCR"
    parts = ["probedesign", "hcrdesign" if hcr else "design", shlex.quote(info["input_name"])]

    def add(flag, value, default):
        if value != default:
            parts.append(_option(flag, value))

    def add_flag(flag, on):
        if on:
            parts.append(flag)

    if hcr:
        add("--pairs", p["n_pairs"], 30)
        add("--amplifier", p["amplifier"], "B1")
        add("--target-gibbs", p["target_gibbs"], -31.0)
        if tuple(p["strict_range"]) != (-35.0, -27.0):
            parts.append(_option("--allowable-gibbs", _range(p["strict_range"])))
        add_flag("--asymmetric-gibbs", p["asymmetric_gibbs"])
        add("--lenient-gibbs-min", p["lenient_gibbs_min"], -42.0)
        add_flag("--asymmetric-bowtie", p["asymmetric_bowtie"])
        add("--pair-spacing", p["pair_spacing"], 2)
    else:
        add("--probes", p["n_probes"], 48)
        if p["mixed_lengths"]:
            parts.append(_option("--oligo-length", f"{p['mixed_lengths'][0]}-{p['mixed_lengths'][1]}"))
        else:
            add("--oligo-length", p["oligo_length"], 20)
        add("--spacer-length", p["spacer_length"], 2)
        add("--target-gibbs", p["target_gibbs"], -23.0)
        if tuple(p["allowable_gibbs"]) != (-26.0, -20.0):
            parts.append(_option("--allowable-gibbs", _range(p["allowable_gibbs"])))

    add("--species", p["species"], "human")
    add_flag("--pseudogene-mask", p["pseudogene_mask"])
    add_flag("--genome-mask", p["genome_mask"])
    if (p["pseudogene_mask"] or p["genome_mask"]) and display_index_dir(info["index_dir"]) != "default":
        parts.append(_option("--index-dir", info["index_dir"]))

    note = ""
    if info["repeatmask_mode"] == "file":
        parts.append(_option("--repeatmask-file", info["repeatmask_name"]))
    elif info["repeatmask_mode"] == "auto":
        if hcr:
            # hcrdesign has no --repeatmask: pass RepeatMasker's output instead
            parts.append(_option("--repeatmask-file", f"{info['input_name']}.masked"))
            note = "   (run RepeatMasker first: hcrdesign has no --repeatmask)"
        else:
            parts.append("--repeatmask")

    add("--hp-threshold", p["hp_threshold"], 5)
    add("--di-threshold", p["di_threshold"], 3)
    if hcr:
        add("--resolve-spacer", p["resolve_spacer"], None)
    add_flag("--save-bowtie-raw", p["save_bowtie_raw"])
    if result.template_name != Path(info["input_name"]).stem:
        parts.append(_option("--output", result.template_name))

    return " ".join(parts) + note


# ── Header ───────────────────────────────────────────────────────────────────

def _line(label: str = "", value="") -> str:
    if not label:
        return "#"
    return f"# {label + ':':<17} {value}"


def _mask_line(result, letter: str, requested: bool) -> str:
    if not requested:
        return "off"
    record = next((m for m in result.run_info["masks"] if m.letter == letter), None)
    if record is None:
        return "requested, but no result recorded"
    if record.status == "failed":
        return f"FAILED, not applied: {record.error}"
    fields = ["applied", record.index]
    if record.assembly:
        fields.append(record.assembly)
    fields.append(f"sha {record.fingerprint}" if record.fingerprint else "sha unknown")
    return "  ".join(fields)


def _repeat_line(info: dict) -> str:
    mode = info["repeatmask_mode"]
    if mode == "file":
        return f"file {info['repeatmask_name']}"
    if mode == "auto":
        return "RepeatMasker (auto)"
    if mode == "input_n":
        return "N's in input"
    return "none"


def _coverage(result) -> str:
    total = result.run_info["length_nt"] or 1
    parts = []
    for kind, mask in zip(result.mask_kinds, result.mask_strings):
        parts.append(f"{kind} {100 * mask.count(kind) / total:.1f}%")
    return "  ".join(parts) or "none"


def format_header(result) -> str:
    """'#' header for _seq.txt / _HCR_seq.txt; '' if the result has no run_info."""
    info = result.run_info
    if not info:
        return ""
    p = info["params"]
    hcr = info["tool"] == "HCR"
    git = info["git"] or "unknown"
    entries = info["n_entries"]

    lines = [
        f"# ProbeDesign {info['version']} (git {git}) · {info['tool']} · {info['timestamp']}",
        f"# command: {build_command(result)}",
        _line(),
        "# == Input ==",
        _line("file", "pasted sequence" if info["input_name"] == PASTED_INPUT_NAME else info["input_name"]),
        _line("name", result.template_name),
        _line("length", f"{info['length_nt']:,} nt ({entries} {'entry' if entries == 1 else 'entries'})"),
        _line(),
        "# == Design parameters ==",
    ]
    low_complexity = f"homopolymer >= {p['hp_threshold']}, dinucleotide repeat >= {p['di_threshold']}"
    if hcr:
        lenient = f"on, lenient half down to {_fmt(p['lenient_gibbs_min'])} kcal/mol"
        lines += [
            _line("pairs", f"{p['n_pairs']} requested"),
            _line("amplifier", p["amplifier"]),
            _line("target dG", f"{_fmt(p['target_gibbs'])} kcal/mol per 25-nt half"),
            _line("strict dG range", f"{_fmt(p['strict_range'][0])} to {_fmt(p['strict_range'][1])} kcal/mol"),
            _line("asymmetric dG", lenient if p["asymmetric_gibbs"] else "off"),
            _line("asymmetric mask", "on (one half may overlap P/B masks)" if p["asymmetric_bowtie"] else "off"),
            _line("pair spacing", p["pair_spacing"]),
            _line("species", p["species"]),
            _line("low complexity", low_complexity),
            _line("resolve spacer", p["resolve_spacer"] or "none (WW)"),
        ]
    else:
        if p["mixed_lengths"]:
            length = f"{p['mixed_lengths'][0]}-{p['mixed_lengths'][1]} (mixed)"
        else:
            length = p["oligo_length"]
        lines += [
            _line("probes", f"{p['n_probes']} requested"),
            _line("oligo length", length),
            _line("spacer", p["spacer_length"]),
            _line("target dG", f"{_fmt(p['target_gibbs'])} kcal/mol"),
            _line("allowable dG", f"{_fmt(p['allowable_gibbs'][0])} to {_fmt(p['allowable_gibbs'][1])} kcal/mol"),
            _line("species", p["species"]),
            _line("low complexity", low_complexity),
        ]

    lines += [
        _line(),
        "# == Masking ==",
        _line("repeat mask", _repeat_line(info)),
        _line("pseudogene mask", _mask_line(result, "P", p["pseudogene_mask"])),
        _line("genome mask", _mask_line(result, "B", p["genome_mask"])),
    ]
    if p["pseudogene_mask"] or p["genome_mask"]:
        lines += [
            _line("index dir", display_index_dir(info["index_dir"])),
            _line("bowtie", info["bowtie"] or "not found"),
        ]

    found = len(result.pairs) if hcr else len(result.probes)
    wanted = p["n_pairs"] if hcr else p["n_probes"]
    score = f"{result.score:.4f}" if math.isfinite(result.score) else "n/a"
    masks = ", ".join(f"{kind} mask" for kind in result.mask_kinds)
    block = (f"position, sequence, {masks}, pairs, labels" if hcr
             else f"sequence, {masks}, probes, labels").replace(", , ", ", ")
    lines += [
        _line(),
        "# == Result ==",
        _line("pairs found" if hcr else "probes found", f"{found} / {wanted}"),
        _line("score", score),
        _line("masked", _coverage(result)),
        _line(),
        "# == Sequence overlay ==",
        f"# each block: {block}",
    ]
    return "\n".join(lines) + "\n"
