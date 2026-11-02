"""Tests for the run header that _seq.txt / _HCR_seq.txt start with.

The header is a block of '#' lines (see probedesign/provenance.py) recording
the parameters, input, masks and software behind a design. It sits above the
sequence overlay, which must stay exactly as it was, except that no overlay
line may start with '#' (a label wrapped at its '#' loses it). These tests cover the
header text, the `command:` line (build_command), the per-mask records
(applied / failed, index fingerprints) and display_index_dir.

Designs that need no bowtie (CDKN1A with a repeat-mask file, EIF1 for HCR)
run in well under a second. Tests that need bowtie skip cleanly when bowtie
or the index files are missing, the same way test_smfish_regression.py does.
"""

import contextlib
import copy
import dataclasses
import hashlib
import inspect
import io
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pytest
from click.testing import CliRunner

# Ensure probedesign package is importable
_repo_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_repo_root / "src"))

from probedesign import __version__
from probedesign.cli import main
from probedesign.core import design_probes
from probedesign.fasta import read_fasta, write_fasta
from probedesign.hcr import design_hcr_probes
from probedesign.hcr_output import format_hcr_oligos, format_hcr_seq, write_hcr_output_files
from probedesign.masking import (
    DEFAULT_INDEX_DIR,
    GENOME_INDEXES,
    PSEUDOGENE_INDEXES,
    bowtie_version,
    find_bowtie,
    index_fingerprint,
)
from probedesign.output import format_seq_content, write_seq_file
from probedesign.provenance import (
    PASTED_INPUT_NAME,
    bowtie_mask_record,
    build_command,
    display_index_dir,
    failed_masks,
    format_header,
    software_info,
)


REGRESSION = _repo_root / "tests" / "regression"
INDEX_DIR = _repo_root / "bowtie_indexes"

CDKN1A_DIR = REGRESSION / "CDKN1A_32"
CDKN1A_FA = CDKN1A_DIR / "CDKN1A.fa"
CDKN1A_REPEATMASKED = CDKN1A_DIR / "CDKN1A_repeatmasked.fa"
EIF1_DIR = REGRESSION / "EIF1_CDS_HCR"
EIF1_FA = EIF1_DIR / "EIF1_Exons.fasta"
KRT19_FA = REGRESSION / "KRT19_withUTRs" / "KRT19_withUTRs.fa"
MIXED_FA = REGRESSION / "mixed_length_test" / "mixed_gc.fa"

HUMAN_GENOME_INDEX = GENOME_INDEXES["human"]

# Switches the L (low-complexity) mask off, as in test_smfish_regression.py
NO_L_MASK = {"hp_threshold": 10, "di_threshold": 10}

# Header lines whose value is not aligned to column 20 (long free-form text)
UNALIGNED_LABELS = ("command", "each block")


# ── Helpers ──────────────────────────────────────────────────────────────────

def _quiet(func, *args, **kwargs):
    """Call func with its progress prints swallowed."""
    with contextlib.redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


def _smfish(input_file=CDKN1A_FA, **kwargs):
    """design_probes() on CDKN1A (or another input) with the repeat-mask file by default."""
    if input_file == CDKN1A_FA:
        kwargs.setdefault("repeatmask_file", str(CDKN1A_REPEATMASKED))
    return _quiet(design_probes, str(input_file), **kwargs)


def _hcr(input_file=EIF1_FA, **kwargs):
    """design_hcr_probes() on EIF1 (or another input)."""
    return _quiet(design_hcr_probes, str(input_file), **kwargs)


def _render(result):
    """The _seq.txt / _HCR_seq.txt text for either kind of result."""
    if hasattr(result, "pairs"):
        return format_hcr_seq(result)
    return format_seq_content(result, mask_seqs=result.mask_strings)


def _header_lines(text):
    """The leading '#' lines of an output file; no '#' line may follow them."""
    lines = text.splitlines()
    n = 0
    while n < len(lines) and lines[n].startswith("#"):
        n += 1
    assert not any(line.startswith("#") for line in lines[n:]), "'#' line after the header"
    return lines[:n]


def _without_header(text):
    """The text with every '#' line removed (what the overlay alone looks like)."""
    return "".join(line for line in text.splitlines(keepends=True) if not line.startswith("#"))


def _fields(header):
    """{label: value} for the aligned '# label:   value' lines, checking the column."""
    fields = {}
    for line in header:
        m = re.match(r"# ([A-Za-z][A-Za-z ]*?):( +)(\S.*)$", line)
        if not m or m.group(1) in UNALIGNED_LABELS:
            continue
        assert m.start(3) == 20, f"value does not start in column 20: {line!r}"
        fields[m.group(1)] = m.group(3)
    return fields


def _header_command(header):
    """The text after '# command: '."""
    (line,) = [l for l in header if l.startswith("# command: ")]
    return line[len("# command: "):]


def _command_tokens(command):
    """Shell tokens of a header command."""
    return shlex.split(command)


def _command_options(tokens):
    """{flag: value or True} for the options after `probedesign <tool> <input>`."""
    options = {}
    i = 0
    while i < len(tokens):
        token = tokens[i]
        assert token.startswith("--"), f"expected a long flag, got {token!r}"
        if "=" in token:
            flag, value = token.split("=", 1)
        elif i + 1 < len(tokens) and not tokens[i + 1].startswith("--"):
            flag, value = token, tokens[i + 1]
            i += 1
        else:
            flag, value = token, True
        options[flag] = value
        i += 1
    return options


def _drop_output(tokens):
    """Tokens without `--output NAME`."""
    kept = []
    skip = False
    for token in tokens:
        if skip:
            skip = False
        elif token == "--output":
            skip = True
        else:
            kept.append(token)
    return kept


def _overlay_blocks(text):
    """The overlay (header removed) as a list of blocks, each a list of lines.

    A block is the lines between blank lines: for smFISH one sequence line,
    one line per mask, the probe alignment and the label line.
    """
    blocks, current = [], []
    for line in _without_header(text).splitlines():
        if line == "":
            blocks.append(current)
            current = []
        else:
            current.append(line)
    assert not current, "overlay does not end with a blank line"
    return blocks


def _oligo_sequences(oligos_file):
    """Probe sequences in file order (second-to-last column of *_oligos.txt)."""
    with open(oligos_file) as f:
        return [line.rstrip("\n").split("\t")[-2] for line in f if line.strip()]


def _hcr_oligo_rows(text):
    """(pair, half, full_oligo) per data row of an _HCR_oligos.txt text."""
    rows = [line.split("\t") for line in text.splitlines()[1:] if line.strip()]
    return [(row[0], row[1], row[-2]) for row in rows]


def _write_three_entry_fasta(path):
    """Write mixed_gc.fa's sequence as three entries; returns the sequence length."""
    _, seqs = read_fasta(str(MIXED_FA))
    full = "".join(seqs)
    cuts = [len(full) // 3, 2 * len(full) // 3]
    write_fasta(str(path), ["a", "b", "c"], [full[:cuts[0]], full[cuts[0]:cuts[1]], full[cuts[1]:]])
    return len(full)


def _require_indices(*indices):
    """Skip the test unless every (index_name, extension) pair is on disk."""
    suffixes = ("1", "2", "3", "4", "rev.1", "rev.2")
    missing = [
        f"{name}.{suffix}.{ext}"
        for name, ext in indices
        for suffix in suffixes
        if not (INDEX_DIR / f"{name}.{suffix}.{ext}").exists()
    ]
    if missing:
        pytest.skip(f"bowtie index files missing in {INDEX_DIR}: {', '.join(missing)}")


# ── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def bowtie_env(monkeypatch):
    """Put the env's bin dir first on PATH, then skip if bowtie is still missing."""
    env_bin = str(Path(sys.executable).parent)
    monkeypatch.setenv("PATH", env_bin + os.pathsep + os.environ.get("PATH", ""))
    if shutil.which("bowtie") is None:
        pytest.skip("bowtie not found on PATH")


@pytest.fixture
def human_indices(bowtie_env):
    """Skip unless the human pseudogene (.ebwt) and genome (.bt2) indices exist."""
    _require_indices(("humanPseudo", "ebwt"), (HUMAN_GENOME_INDEX, "bt2"))


@pytest.fixture
def no_bowtie(monkeypatch, tmp_path):
    """Hide bowtie: PATH is an empty directory and BOWTIEHOME is unset.

    find_bowtie() also looks at /usr/bin/bowtie and ~/bowtie/bowtie, so skip if
    either exists. Returns the empty PATH directory.
    """
    for fixed in ("/usr/bin/bowtie", os.path.expanduser("~/bowtie/bowtie")):
        if os.path.exists(fixed):
            pytest.skip(f"bowtie found at {fixed}")
    software_info()  # cache the git lookup while git is still on PATH
    bin_dir = tmp_path / "empty_bin"
    bin_dir.mkdir()
    monkeypatch.setenv("PATH", str(bin_dir))
    monkeypatch.delenv("BOWTIEHOME", raising=False)
    return bin_dir


# ── CDKN1A through the CLI ───────────────────────────────────────────────────

def test_cdkn1a_cli_seq_file_starts_with_header(tmp_path):
    """`design` writes a '#' header with aligned values, basenames only and no absolute paths."""
    prefix = tmp_path / "cdkn1a"
    result = CliRunner().invoke(main, [
        "design", str(CDKN1A_FA),
        "--probes", "32",
        "--repeatmask-file", str(CDKN1A_REPEATMASKED),
        "--hp-threshold", "10",
        "--di-threshold", "10",
        "-o", str(prefix),
        "--quiet",
    ])
    assert result.exit_code == 0, f"{result.output}\n{result.exception!r}"

    text = (tmp_path / "cdkn1a_seq.txt").read_text()
    lines = text.splitlines()
    assert lines[0].startswith("# ProbeDesign ")
    header = _header_lines(text)
    assert len(header) > 10
    assert "# repeat mask:      file CDKN1A_repeatmasked.fa" in lines
    assert "# probes found:     32 / 32" in lines
    assert "# file:             CDKN1A.fa" in lines
    assert "# name:             cdkn1a" in lines

    # Long flags, only the non-default options, basenames, --output because the name differs
    tokens = _command_tokens(_header_command(header))
    assert tokens[:3] == ["probedesign", "design", "CDKN1A.fa"]
    assert _command_options(tokens[3:]) == {
        "--probes": "32",
        "--repeatmask-file": "CDKN1A_repeatmasked.fa",
        "--hp-threshold": "10",
        "--di-threshold": "10",
        "--output": "cdkn1a",
    }

    for forbidden in {str(tmp_path), str(tmp_path.resolve()), str(_repo_root), str(_repo_root.resolve())}:
        assert forbidden not in text


# ── The header is only an addition ───────────────────────────────────────────

def test_smfish_header_is_only_an_addition(tmp_path):
    """Removing '#' lines leaves exactly the overlay of a result without run_info."""
    result = _smfish(n_probes=32, **NO_L_MASK)
    assert result.run_info and len(result.mask_strings) >= 2
    bare = dataclasses.replace(result, run_info={})

    with_header = format_seq_content(result, mask_seqs=result.mask_strings)
    overlay = format_seq_content(bare, mask_seqs=result.mask_strings)

    assert format_header(bare) == ""
    assert not any(line.startswith("#") for line in overlay.splitlines())
    assert _without_header(with_header) == overlay
    assert with_header == format_header(result) + overlay
    assert _header_lines(with_header)  # non-empty, contiguous, at the top

    # write_seq_file() writes the same text
    path = tmp_path / "out_seq.txt"
    write_seq_file(result, str(path), mask_seqs=result.mask_strings)
    assert path.read_text() == with_header


def test_hcr_header_is_only_an_addition(tmp_path):
    """Same for _HCR_seq.txt, with an R mask line in the overlay."""
    result = _hcr(CDKN1A_FA, n_pairs=8, repeatmask_file=str(CDKN1A_REPEATMASKED))
    assert result.run_info and result.pairs and result.mask_kinds == ["R"]
    bare = dataclasses.replace(result, run_info={})

    with_header = format_hcr_seq(result)
    overlay = format_hcr_seq(bare)

    assert format_header(bare) == ""
    assert not any(line.startswith("#") for line in overlay.splitlines())
    assert _without_header(with_header) == overlay
    assert with_header == format_header(result) + overlay
    assert _header_lines(with_header)

    # write_hcr_output_files() writes the same text
    write_hcr_output_files(result, str(tmp_path / "out"))
    assert (tmp_path / "out_HCR_seq.txt").read_text() == with_header


def test_label_wrapped_at_its_hash_does_not_start_a_row_with_hash():
    """A label that wraps just before its '#' ("Prb|# 2,...") leaves no '#' line after the header."""
    result = _smfish(n_probes=32, **NO_L_MASK)
    probe = result.probes[1]
    width = probe.position + len("Prb")  # the second block starts on this label's '#'

    text = format_seq_content(result, mask_seqs=result.mask_strings, line_width=width)
    assert _header_lines(text)  # fails if any later line starts with '#'
    blocks = _overlay_blocks(text)
    assert blocks[0][-1].endswith("Prb")
    assert blocks[1][-1].startswith(f"  {probe.index},Pos ")


def test_overlay_description_matches_block_layout():
    """The `# each block:` line names as many rows as a block of the overlay has."""
    cases = [
        _smfish(n_probes=10),
        _hcr(EIF1_FA, n_pairs=5),                                                   # no masks
        _hcr(CDKN1A_FA, n_pairs=5, repeatmask_file=str(CDKN1A_REPEATMASKED)),       # one mask
    ]
    for result in cases:
        text = _render(result)
        (line,) = [l for l in _header_lines(text) if l.startswith("# each block: ")]
        described = line[len("# each block: "):].split(", ")
        blocks = _overlay_blocks(text)
        assert blocks
        assert {len(block) for block in blocks} == {len(described)}, line
        assert len(result.mask_kinds) == len(result.mask_strings)


# ── Header content ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("tool", ["smfish", "hcr"])
def test_title_line(tool):
    """First line: version, git revision, tool and an ISO timestamp with a UTC offset."""
    result = _smfish(n_probes=5) if tool == "smfish" else _hcr(n_pairs=5)
    title = format_header(result).splitlines()[0]
    m = re.fullmatch(
        r"# ProbeDesign (\S+) \(git (\S+)\) \| (smFISH|HCR) \| (\S+)", title
    )
    assert m, title
    assert m.group(1) == __version__
    assert m.group(3) == ("smFISH" if tool == "smfish" else "HCR")
    stamp = datetime.fromisoformat(m.group(4))
    assert stamp.tzinfo is not None
    assert abs((datetime.now().astimezone() - stamp).total_seconds()) < 86400


def test_smfish_header_fields_match_design():
    """Every aligned field reflects the design that was run."""
    result = _smfish(
        n_probes=20, mixed_lengths=(18, 22), spacer_length=3,
        target_gibbs=-24.0, allowable_gibbs=(-27.0, -19.0),
        species="mouse", hp_threshold=9, di_threshold=8,
    )
    fields = _fields(_header_lines(format_header(result)))
    seq = result.input_sequence
    length = len(seq) - seq.count(">")

    assert fields["file"] == "CDKN1A.fa"
    assert fields["name"] == "CDKN1A"
    assert fields["length"] == f"{length:,} nt (1 entry)"
    assert fields["probes"] == "20 requested"
    assert fields["oligo length"] == "18-22 (mixed)"
    assert fields["spacer"] == "3"
    assert fields["target dG"] == "-24.0 kcal/mol"
    assert fields["allowable dG"] == "-27.0 to -19.0 kcal/mol"
    assert fields["species"] == "mouse"
    assert fields["low complexity"] == "homopolymer >= 9, dinucleotide repeat >= 8"
    assert fields["repeat mask"] == "file CDKN1A_repeatmasked.fa"
    assert fields["pseudogene mask"] == "off"
    assert fields["genome mask"] == "off"
    assert "index dir" not in fields and "bowtie" not in fields
    assert fields["probes found"] == f"{len(result.probes)} / 20"
    assert fields["score"] == f"{result.score:.4f}"

    # Coverage: one "<letter> <percent>%" per mask line, in mask order
    covered = re.findall(r"(\w) (\d+\.\d)%", fields["masked"])
    assert [kind for kind, _ in covered] == result.mask_kinds
    for (kind, percent), mask in zip(covered, result.mask_strings):
        assert float(percent) == pytest.approx(100 * mask.count(kind) / length, abs=0.051)


def test_hcr_header_fields_match_design():
    """Every aligned HCR field reflects the design, including the non-default options."""
    result = _hcr(
        n_pairs=8, amplifier="B3", target_gibbs=-30.0, strict_range=(-34.0, -28.0),
        asymmetric_gibbs=True, lenient_gibbs_min=-40.0, asymmetric_bowtie=True,
        pair_spacing=4, resolve_spacer="AA", hp_threshold=6, di_threshold=4,
    )
    fields = _fields(_header_lines(format_header(result)))
    seq = result.input_sequence

    assert fields["file"] == "EIF1_Exons.fasta"
    assert fields["name"] == "EIF1_Exons"
    assert fields["length"] == f"{len(seq) - seq.count('>'):,} nt (1 entry)"
    assert fields["pairs"] == "8 requested"
    assert fields["amplifier"] == "B3"
    assert fields["target dG"] == "-30.0 kcal/mol per 25-nt half"
    assert fields["strict dG range"] == "-34.0 to -28.0 kcal/mol"
    assert fields["asymmetric dG"] == "on, lenient half down to -40.0 kcal/mol"
    assert fields["asymmetric mask"] == "on (one half may overlap P/B masks)"
    assert fields["pair spacing"] == "4"
    assert fields["species"] == "human"
    assert fields["low complexity"] == "homopolymer >= 6, dinucleotide repeat >= 4"
    assert fields["resolve spacer"] == "AA"
    assert fields["repeat mask"] == "none"
    assert fields["pseudogene mask"] == "off" and fields["genome mask"] == "off"
    assert fields["pairs found"] == f"{len(result.pairs)} / 8"
    assert fields["masked"] == "none"

    # Defaults read as off / none
    default_fields = _fields(_header_lines(format_header(_hcr(n_pairs=5))))
    assert default_fields["asymmetric dG"] == "off"
    assert default_fields["asymmetric mask"] == "off"
    assert default_fields["resolve spacer"] == "none (WW)"


def test_repeat_mask_modes():
    """The repeat mask line says whether the mask came from a file, N's in the input, or nothing."""
    assert _fields(_header_lines(format_header(_smfish(n_probes=5))))["repeat mask"] == \
        "file CDKN1A_repeatmasked.fa"
    plain = _smfish(n_probes=5, repeatmask_file=None)
    assert _fields(_header_lines(format_header(plain)))["repeat mask"] == "none"


def test_input_with_n_reports_n_mask(tmp_path):
    """A FASTA that already contains N's is reported as "N's in input"."""
    fasta = tmp_path / "with_n.fa"
    _, seqs = read_fasta(str(CDKN1A_REPEATMASKED))
    write_fasta(str(fasta), ["masked"], seqs)
    result = _quiet(design_probes, str(fasta), n_probes=5)
    assert result.mask_kinds[0] == "R"
    assert _fields(_header_lines(format_header(result)))["repeat mask"] == "N's in input"


def test_pasted_sequence_is_named_as_such(tmp_path):
    """The app's pasted-sequence input shows as "pasted sequence" but keeps its name in the command."""
    fasta = tmp_path / PASTED_INPUT_NAME
    shutil.copy(CDKN1A_FA, fasta)
    result = _quiet(design_probes, str(fasta), n_probes=5)
    header = _header_lines(format_header(result))
    assert _fields(header)["file"] == "pasted sequence"
    assert _command_tokens(_header_command(header))[2] == PASTED_INPUT_NAME


def test_multi_entry_input_is_counted(tmp_path):
    """Three FASTA entries give '3 entries' and a length without the '>' junction markers."""
    fasta = tmp_path / "three.fa"
    n_nt = _write_three_entry_fasta(fasta)

    result = _quiet(design_probes, str(fasta), n_probes=10)
    assert result.input_sequence.count(">") == 3
    fields = _fields(_header_lines(format_header(result)))
    assert fields["length"] == f"{n_nt:,} nt (3 entries)"


# ── Mixed-length overlay ─────────────────────────────────────────────────────

@pytest.mark.parametrize("entries", [1, 3])
def test_mixed_length_overlay_matches_probe_lengths(entries, tmp_path):
    """In the overlay each probe's alignment run is as long as that probe (18-22 nt mix)."""
    fasta = MIXED_FA
    if entries > 1:
        # Same sequence in several entries, so the overlay contains '>' junction markers
        fasta = tmp_path / "split.fa"
        _write_three_entry_fasta(fasta)

    result = _quiet(design_probes, str(fasta), n_probes=20, spacer_length=2, mixed_lengths=(18, 22))
    probes = result.probes
    lengths = {len(p.sequence) for p in probes}
    assert len(lengths) > 1, f"expected a genuine mix of lengths, got {lengths}"
    assert all(18 <= n <= 22 for n in lengths)

    text = format_seq_content(result, mask_seqs=result.mask_strings)
    n_masks = len(result.mask_strings)
    blocks = _overlay_blocks(text)
    assert all(len(block) == n_masks + 3 for block in blocks)
    assert all(len(block[0]) == 110 for block in blocks[:-1])
    assert 0 < len(blocks[-1][0]) <= 110

    # Re-join the blocks into whole-sequence rows
    seq = "".join(block[0] for block in blocks)
    masks = ["".join(block[1 + k] for block in blocks) for k in range(n_masks)]
    align = "".join(block[1 + n_masks] for block in blocks)
    assert seq == result.input_sequence
    assert masks == result.mask_strings
    assert len(align) == len(seq)

    for probe in probes:
        run = []
        i = probe.position
        while i < len(seq):
            if seq[i] == ">":
                assert align[i] == " "  # junction markers carry no alignment
            elif align[i] != " ":
                run.append(align[i])
            else:
                break
            i += 1
        assert len(run) == len(probe.sequence), (
            f"probe {probe.index} (len {len(probe.sequence)}) has an alignment run of {len(run)}"
        )
        # The alignment shows the complement of the template, i.e. the probe read backwards
        assert "".join(run)[::-1].lower() == probe.sequence.lower()

    # Nothing is drawn outside the probes
    assert sum(c != " " for c in align) == sum(len(p.sequence) for p in probes)


# ── The command: round trip ──────────────────────────────────────────────────

def _run_command(command, workdir, out_dir, monkeypatch):
    """Run a header `command:` through the CLI from workdir; returns the output prefix.

    The file arguments in the command are basenames, so it has to run from the
    directory that holds them. --output is pointed into out_dir (or added).
    """
    tokens = shlex.split(command)
    assert tokens[0] == "probedesign"
    args = tokens[1:]
    if "--output" in args:
        i = args.index("--output")
        args[i + 1] = str(out_dir / args[i + 1])
        prefix = Path(args[i + 1])
    else:
        prefix = out_dir / "rt"
        args += ["--output", str(prefix)]
    monkeypatch.chdir(workdir)
    result = CliRunner().invoke(main, args + ["--quiet"])
    assert result.exit_code == 0, f"{command}\n{result.output}\n{result.exception!r}"
    return prefix


SMFISH_ROUND_TRIPS = {
    "defaults": dict(n_probes=30, **NO_L_MASK),
    "mixed-spacer-named": dict(n_probes=30, mixed_lengths=(17, 23), spacer_length=2,
                               output_name="myprobes", **NO_L_MASK),
    "fixed-gibbs": dict(n_probes=20, oligo_length=22, target_gibbs=-24.5,
                        allowable_gibbs=(-27.5, -19.5), hp_threshold=7, di_threshold=4),
    "mixed-spacer": dict(n_probes=30, mixed_lengths=(18, 22), spacer_length=3, **NO_L_MASK),
}


@pytest.mark.parametrize("case", SMFISH_ROUND_TRIPS)
def test_smfish_command_round_trip(case, tmp_path, monkeypatch):
    """Running the header's command again gives the same probes and the same overlay."""
    result = _smfish(**SMFISH_ROUND_TRIPS[case])
    assert len(result.probes) >= 10

    command = build_command(result)
    prefix = _run_command(command, CDKN1A_DIR, tmp_path, monkeypatch)

    assert _oligo_sequences(f"{prefix}_oligos.txt") == [p.sequence for p in result.probes]

    rt_text = Path(f"{prefix}_seq.txt").read_text()
    assert _without_header(rt_text) == _without_header(_render(result))
    # The rerun describes itself with the same command (apart from where it wrote)
    assert _drop_output(_command_tokens(_header_command(_header_lines(rt_text)))) == \
        _drop_output(_command_tokens(command))


HCR_ROUND_TRIPS = {
    "basic": dict(n_pairs=8, amplifier="B3"),
    "named": dict(n_pairs=8, amplifier="B3", output_name="myhcr"),
    "many-options": dict(
        n_pairs=6, amplifier="B9", target_gibbs=-30.0, strict_range=(-34.0, -28.0),
        asymmetric_gibbs=True, lenient_gibbs_min=-40.0, asymmetric_bowtie=True,
        pair_spacing=5, resolve_spacer="AA", hp_threshold=6, di_threshold=4,
    ),
}


@pytest.mark.parametrize("case", HCR_ROUND_TRIPS)
def test_hcr_command_round_trip(case, tmp_path, monkeypatch):
    """Running the header's hcrdesign command again gives the same pairs and the same overlay."""
    result = _hcr(**HCR_ROUND_TRIPS[case])
    assert len(result.pairs) >= 3

    command = build_command(result)
    assert command.startswith("probedesign hcrdesign EIF1_Exons.fasta")
    prefix = _run_command(command, EIF1_DIR, tmp_path, monkeypatch)

    rt_oligos = Path(f"{prefix}_HCR_oligos.txt").read_text()
    rows = _hcr_oligo_rows(rt_oligos)
    assert len({pair for pair, _, _ in rows}) == len(result.pairs)
    assert rows == _hcr_oligo_rows(format_hcr_oligos(result))

    rt_text = Path(f"{prefix}_HCR_seq.txt").read_text()
    assert _without_header(rt_text) == _without_header(format_hcr_seq(result))
    assert _drop_output(_command_tokens(_header_command(_header_lines(rt_text)))) == \
        _drop_output(_command_tokens(command))


# ── The command: flag drift ──────────────────────────────────────────────────

# Every parameter of a design set to a non-default value ...
SMFISH_NON_DEFAULT = dict(
    n_probes=7, oligo_length=22, mixed_lengths=(16, 24), spacer_length=5,
    target_gibbs=-24.5, allowable_gibbs=(-27.0, -19.0), species="mouse",
    pseudogene_mask=True, genome_mask=True, hp_threshold=7, di_threshold=5,
    save_bowtie_raw=True,
)
HCR_NON_DEFAULT = dict(
    n_pairs=11, amplifier="B9", target_gibbs=-30.0, strict_range=(-34.0, -28.0),
    asymmetric_gibbs=True, lenient_gibbs_min=-40.0, asymmetric_bowtie=True,
    pair_spacing=4, species="mouse", pseudogene_mask=True, genome_mask=True,
    hp_threshold=7, di_threshold=5, resolve_spacer="AA", save_bowtie_raw=True,
)

# ... and what the CLI should parse back from the command (click parameter names)
SMFISH_PARSED = dict(
    n_probes=7, spacer_length=5, target_gibbs=-24.5, allowable_gibbs="-27.0,-19.0",
    species="mouse", pseudogene_mask=True, genome_mask=True, hp_threshold=7,
    di_threshold=5, save_bowtie_raw=True, output="renamed",
)
HCR_PARSED = dict(
    n_pairs=11, amplifier="B9", target_gibbs=-30.0, allowable_gibbs="-34.0,-28.0",
    asymmetric_gibbs=True, lenient_gibbs_min=-40.0, asymmetric_bowtie=True,
    pair_spacing=4, species="mouse", pseudogene_mask=True, genome_mask=True,
    hp_threshold=7, di_threshold=5, resolve_spacer="AA", save_bowtie_raw=True,
    output="renamed",
)

# Options of the CLI commands that are not design parameters
NOT_REPRODUCED = {"--quiet", "--no-pseudogene-mask", "--no-genome-mask", "--no-repeatmask"}


def _cli_options(command_name):
    """Every option string (long and short, on and off forms) of a CLI command."""
    return {
        opt
        for param in main.commands[command_name].params
        for opt in list(param.opts) + list(param.secondary_opts)
    }


def _maximal_result(tool, tmp_path, *, variant):
    """A result whose run_info has every parameter non-default (via a real design call).

    variant: "file" (repeat mask from a file) or "auto" (RepeatMasker), and for
    smFISH also "fixed" instead of "mixed" lengths: "file-fixed", "auto-mixed".
    """
    if tool == "smfish":
        result = _smfish(n_probes=3)
        non_default = dict(SMFISH_NON_DEFAULT)
        if variant.endswith("fixed"):
            non_default["mixed_lengths"] = None
    else:
        result = _hcr(n_pairs=3)
        non_default = dict(HCR_NON_DEFAULT)
    info = copy.deepcopy(result.run_info)
    assert set(non_default) == set(info["params"]), (
        "run_info['params'] changed: update the non-default values above and teach "
        "build_command() (and this test) about the new parameter"
    )
    info["params"].update(non_default)
    info["index_dir"] = str(tmp_path)
    if variant.startswith("auto"):
        info["repeatmask_mode"], info["repeatmask_name"] = "auto", None
    else:
        info["repeatmask_mode"], info["repeatmask_name"] = "file", "my_repeatmasked.fa"
    return dataclasses.replace(result, run_info=info, template_name="renamed")


FLAG_CASES = [
    ("smfish", "file-fixed"),
    ("smfish", "auto-mixed"),
    ("hcr", "file"),
    ("hcr", "auto"),
]


@pytest.mark.parametrize("tool,variant", FLAG_CASES)
def test_command_flags_exist_and_parse_back(tool, variant, tmp_path, monkeypatch):
    """Every flag in the command is a real CLI option, and the CLI reads back the same values."""
    result = _maximal_result(tool, tmp_path, variant=variant)
    command_name = "design" if tool == "smfish" else "hcrdesign"
    info = result.run_info
    tokens = _command_tokens(build_command(result))

    assert tokens[:2] == ["probedesign", command_name]
    assert tokens[2] == info["input_name"]
    flags = {t.split("=", 1)[0] for t in tokens[3:] if t.startswith("--")}
    assert flags <= _cli_options(command_name), flags - _cli_options(command_name)
    assert not any(re.fullmatch(r"-[A-Za-z]", t) for t in tokens), "short flag in command"

    # Which flags: everything that was set, and nothing else
    if tool == "smfish":
        expected_flags = {
            "--probes", "--oligo-length", "--spacer-length", "--target-gibbs",
            "--allowable-gibbs", "--species", "--pseudogene-mask", "--genome-mask",
            "--index-dir", "--hp-threshold", "--di-threshold", "--save-bowtie-raw", "--output",
        }
        expected_flags.add("--repeatmask" if variant.startswith("auto") else "--repeatmask-file")
    else:
        expected_flags = {
            "--pairs", "--amplifier", "--target-gibbs", "--allowable-gibbs",
            "--asymmetric-gibbs", "--lenient-gibbs-min", "--asymmetric-bowtie",
            "--pair-spacing", "--species", "--pseudogene-mask", "--genome-mask",
            "--index-dir", "--repeatmask-file", "--hp-threshold", "--di-threshold",
            "--resolve-spacer", "--save-bowtie-raw", "--output",
        }
    assert flags == expected_flags

    # The CLI parses the command to the values it was built from. The files it
    # names must exist (click.Path(exists=True)), so run in a directory with them.
    for name in (info["input_name"], info["repeatmask_name"], f"{info['input_name']}.masked"):
        if name:
            (tmp_path / name).touch()
    monkeypatch.chdir(tmp_path)
    ctx = main.commands[command_name].make_context(command_name, tokens[2:])

    expected = dict(SMFISH_PARSED if tool == "smfish" else HCR_PARSED, index_dir=str(tmp_path))
    if tool == "smfish":
        expected["oligo_length"] = "22" if variant.endswith("fixed") else "16-24"
        if variant.startswith("auto"):
            expected.update(repeatmask=True, repeatmask_file=None)
        else:
            expected.update(repeatmask=False, repeatmask_file="my_repeatmasked.fa")
    else:
        # hcrdesign has no --repeatmask: RepeatMasker's output is passed as a file
        expected["repeatmask_file"] = (
            f"{info['input_name']}.masked" if variant == "auto" else "my_repeatmasked.fa"
        )
    assert {key: ctx.params[key] for key in expected} == expected


@pytest.mark.parametrize("tool", ["smfish", "hcr"])
def test_every_cli_option_is_reproduced_by_the_command(tool, tmp_path):
    """A CLI option that build_command() never emits would make commands silently incomplete."""
    command_name = "design" if tool == "smfish" else "hcrdesign"
    variants = ["file-fixed", "auto-mixed"] if tool == "smfish" else ["file", "auto"]
    emitted = set()
    for variant in variants:
        tokens = _command_tokens(build_command(_maximal_result(tool, tmp_path, variant=variant)))
        emitted |= {t.split("=", 1)[0] for t in tokens if t.startswith("--")}
    long_options = {opt for opt in _cli_options(command_name) if opt.startswith("--")}
    missing = long_options - emitted - NOT_REPRODUCED
    assert not missing, f"CLI options not reproduced in the header command: {sorted(missing)}"


def test_smfish_defaults_agree_between_api_cli_and_command():
    """API defaults give a bare command, and the CLI defaults are the API defaults."""
    result = _quiet(design_probes, str(CDKN1A_FA))
    assert build_command(result) == "probedesign design CDKN1A.fa"

    cli = main.commands["design"].make_context("design", [str(CDKN1A_FA)]).params
    api = inspect.signature(design_probes).parameters
    assert int(cli["oligo_length"]) == api["oligo_length"].default
    assert tuple(map(float, cli["allowable_gibbs"].split(","))) == api["allowable_gibbs"].default
    for name in ("n_probes", "spacer_length", "target_gibbs", "species", "pseudogene_mask",
                 "genome_mask", "index_dir", "repeatmask_file", "save_bowtie_raw",
                 "hp_threshold", "di_threshold"):
        assert cli[name] == api[name].default, name
    assert api["mixed_lengths"].default is None


def test_hcr_defaults_agree_between_api_cli_and_command():
    """Same for hcrdesign."""
    result = _hcr()
    assert build_command(result) == "probedesign hcrdesign EIF1_Exons.fasta"

    cli = main.commands["hcrdesign"].make_context("hcrdesign", [str(EIF1_FA)]).params
    api = inspect.signature(design_hcr_probes).parameters
    assert tuple(map(float, cli["allowable_gibbs"].split(","))) == api["strict_range"].default
    for name in ("n_pairs", "amplifier", "target_gibbs", "asymmetric_gibbs", "lenient_gibbs_min",
                 "asymmetric_bowtie", "pair_spacing", "species", "pseudogene_mask", "genome_mask",
                 "index_dir", "repeatmask_file", "hp_threshold", "di_threshold",
                 "resolve_spacer", "save_bowtie_raw"):
        assert cli[name] == api[name].default, name


# ── Mask records (bowtie masks) ──────────────────────────────────────────────

def _check_applied_krt19_masks(result):
    """Both human masks recorded as applied, with fingerprints that match the index files."""
    records = {m.letter: m for m in result.run_info["masks"]}
    assert set(records) == {"P", "B"}
    assert failed_masks(result) == []
    assert result.mask_kinds[:2] == ["P", "B"]

    names = {"P": PSEUDOGENE_INDEXES["human"], "B": GENOME_INDEXES["human"]}
    extensions = {"P": "ebwt", "B": "bt2"}
    for letter, record in records.items():
        assert record.status == "applied" and record.error is None
        assert record.index == names[letter]
        fingerprint = index_fingerprint(DEFAULT_INDEX_DIR, names[letter])
        assert re.fullmatch(r"[0-9a-f]{12}", fingerprint)
        assert record.fingerprint == fingerprint
        # ... computed independently: SHA-256 of the index's .3 file
        third_file = DEFAULT_INDEX_DIR / f"{names[letter]}.3.{extensions[letter]}"
        assert fingerprint == hashlib.sha256(third_file.read_bytes()).hexdigest()[:12]
    assert "GRCh38" in records["B"].assembly and records["P"].assembly is None

    text = _render(result)
    fields = _fields(_header_lines(text))
    assert fields["pseudogene mask"] == (
        f"applied  {names['P']}  sha {records['P'].fingerprint}"
    )
    assert fields["genome mask"] == (
        f"applied  {names['B']}  {records['B'].assembly}  sha {records['B'].fingerprint}"
    )
    assert fields["index dir"] == "default"
    assert result.run_info["bowtie"] and fields["bowtie"] == bowtie_version()
    for absolute in (str(_repo_root), str(INDEX_DIR)):
        assert absolute not in text


def test_krt19_smfish_masks_are_recorded(human_indices):
    """KRT19 with pseudogene + genome masks: both recorded as applied, index shown as default."""
    result = _smfish(KRT19_FA, n_probes=32, pseudogene_mask=True, genome_mask=True, **NO_L_MASK)
    _check_applied_krt19_masks(result)
    fields = _fields(_header_lines(format_header(result)))
    assert fields["probes found"] == f"{len(result.probes)} / 32"


def test_krt19_hcr_masks_are_recorded(human_indices):
    """Same for hcrdesign (KRT19 leaves no valid pair under both masks: the early-return path)."""
    result = _hcr(KRT19_FA, n_pairs=10, pseudogene_mask=True, genome_mask=True)
    _check_applied_krt19_masks(result)
    fields = _fields(_header_lines(format_header(result)))
    assert fields["pairs found"] == f"{len(result.pairs)} / 10"


@pytest.mark.parametrize("tool", ["smfish", "hcr"])
def test_broken_index_marks_genome_mask_failed(tool, bowtie_env, tmp_path):
    """bowtie exits non-zero on an empty index dir: the genome mask is recorded as failed."""
    index_dir = str(tmp_path)
    if tool == "smfish":
        result = _smfish(genome_mask=True, index_dir=index_dir)
        assert result.probes
    else:
        result = _hcr(n_pairs=5, genome_mask=True, index_dir=index_dir)
        assert result.pairs

    (record,) = failed_masks(result)
    assert (record.letter, record.status) == ("B", "failed")
    assert record.index == GENOME_INDEXES["human"]
    assert record.error and record.fingerprint is None
    assert "B" not in result.mask_kinds

    text = _render(result)
    header = _header_lines(text)
    fields = _fields(header)
    assert fields["genome mask"].startswith("FAILED, not applied: ")
    assert fields["genome mask"] == f"FAILED, not applied: {record.error}"
    assert fields["pseudogene mask"] == "off"
    # A custom index dir is printed as given: on the index dir and command lines only
    assert fields["index dir"] == index_dir
    assert f"--index-dir {shlex.quote(index_dir)}" in _header_command(header)
    for line in text.splitlines():
        if index_dir in line:
            assert line.startswith(("# index dir:", "# command:")), line
    assert index_dir not in record.error


@pytest.mark.parametrize("tool", ["smfish", "hcr"])
def test_missing_bowtie_marks_pseudogene_mask_failed(tool, no_bowtie):
    """With no bowtie on PATH the pseudogene mask is recorded as failed and the design carries on."""
    if tool == "smfish":
        result = _smfish(pseudogene_mask=True)
        assert result.probes
    else:
        result = _hcr(n_pairs=5, pseudogene_mask=True)
        assert result.pairs

    (record,) = failed_masks(result)
    assert (record.letter, record.status) == ("P", "failed")
    assert record.error
    assert "P" not in result.mask_kinds
    assert result.run_info["bowtie"] is None

    fields = _fields(_header_lines(format_header(result)))
    assert fields["pseudogene mask"] == f"FAILED, not applied: {record.error}"
    assert fields["genome mask"] == "off"
    assert fields["bowtie"] == "not found"
    assert record.error.startswith("Could not find bowtie")


def test_non_executable_bowtie_is_not_used(no_bowtie, tmp_path, monkeypatch):
    """A BOWTIEHOME/bowtie that can't be executed doesn't count as found."""
    home = tmp_path / "bowtie_home"
    home.mkdir()
    (home / "bowtie").write_text("not a program\n")
    monkeypatch.setenv("BOWTIEHOME", str(home))
    with pytest.raises(FileNotFoundError):
        find_bowtie()


def test_failure_reason_never_contains_a_path():
    """Errors other than bowtie's own are reduced to their type: their text may name server paths."""
    secret = "/srv/private/home/user/bowtie"
    record = bowtie_mask_record("P", "human", None, error=PermissionError(13, "Permission denied", secret))
    assert record.status == "failed"
    assert secret not in record.error
    assert record.error.startswith("PermissionError")

    unknown = bowtie_mask_record("B", "martian", None, error=ValueError("Unknown species: martian"))
    assert (unknown.index, unknown.error) == ("?", "Unknown species: martian")


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0, reason="root can read any file")
def test_unreadable_index_has_no_fingerprint(tmp_path):
    """index_fingerprint() returns None rather than raising (it runs inside except blocks)."""
    small = tmp_path / "humanPseudo.3.ebwt"
    small.write_bytes(b"layout")
    assert index_fingerprint(tmp_path, "humanPseudo") == hashlib.sha256(b"layout").hexdigest()[:12]
    small.chmod(0)
    try:
        assert index_fingerprint(tmp_path, "humanPseudo") is None
    finally:
        small.chmod(0o644)


def test_hcr_repeatmasker_command_is_pasteable():
    """HCR + RepeatMasker (auto, app only): the command parses as shell and the note is its own line."""
    result = _hcr(n_pairs=5)
    result.run_info["repeatmask_mode"], result.run_info["repeatmask_name"] = "auto", None
    header = _header_lines(format_header(result))
    command = _header_command(header)
    assert f"--repeatmask-file {EIF1_FA.name}.masked" in command
    assert "(" not in command
    tokens = shlex.split(command)
    assert all(not t.startswith("--") or t in _cli_options("hcrdesign") or "=" in t for t in tokens)
    assert _fields(header)["command note"].startswith(f"run RepeatMasker on {EIF1_FA.name} first")
    bash = shutil.which("bash")
    if bash:
        assert subprocess.run([bash, "-n", "-c", command]).returncode == 0


def test_block_description_lists_only_drawn_masks():
    """format_seq_content() without mask_seqs draws no mask lines, and the header says so."""
    result = _smfish(n_probes=5)
    assert result.mask_kinds  # F and L at least
    header = _header_lines(format_seq_content(result))
    assert header[-1] == "# each block: sequence, probes, labels"
    with_masks = _header_lines(format_seq_content(result, mask_seqs=result.mask_strings))
    assert with_masks[-1] == "# each block: sequence, " + ", ".join(
        f"{k} mask" for k in result.mask_kinds) + ", probes, labels"


# ── display_index_dir ────────────────────────────────────────────────────────

def test_display_index_dir(tmp_path, monkeypatch):
    """The default index dir, however written, shows as 'default'; anything else as given."""
    default = str(DEFAULT_INDEX_DIR)
    assert display_index_dir(None) == "default"
    assert display_index_dir("") == "default"
    assert display_index_dir(default) == "default"
    assert display_index_dir(default + "/") == "default"
    assert display_index_dir(str(DEFAULT_INDEX_DIR / ".." / DEFAULT_INDEX_DIR.name)) == "default"

    monkeypatch.chdir(tmp_path)
    assert display_index_dir(os.path.relpath(DEFAULT_INDEX_DIR, tmp_path)) == "default"

    assert display_index_dir(str(tmp_path)) == str(tmp_path)
    assert display_index_dir(default + "_other") == default + "_other"
    assert display_index_dir("custom_indexes") == "custom_indexes"
