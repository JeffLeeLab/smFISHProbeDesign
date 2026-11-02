"""Regression tests for the smFISH probe-design pipeline (design_probes + CLI).

The reference outputs in tests/regression/ predate the low-complexity "L"
mask, so every comparison against them passes hp_threshold=10 and
di_threshold=10 to switch that filter off. Probe sets are compared as sorted
lists of sequences (the second-to-last column of a *_oligos.txt file), the
same way setup_all.sh does.

Tests that need bowtie skip cleanly when bowtie or the required index files
are missing. They prepend the running interpreter's bin directory to PATH so
the conda env's bowtie is found even if the env is not activated.
"""

import os
import shutil
import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

# Ensure probedesign package is importable
_repo_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_repo_root / "src"))

from probedesign.cli import main
from probedesign.core import design_probes
from probedesign.masking import has_homopolymer, has_dinucleotide_repeat


REGRESSION = _repo_root / "tests" / "regression"
INDEX_DIR = _repo_root / "bowtie_indexes"

CDKN1A_DIR = REGRESSION / "CDKN1A_32"
CDKN1A_FA = CDKN1A_DIR / "CDKN1A.fa"
CDKN1A_REPEATMASKED = CDKN1A_DIR / "CDKN1A_repeatmasked.fa"
CDKN1A_REF = CDKN1A_DIR / "CDKN1A_32_genomemaskoff_oligos.txt"

KRT19_DIR = REGRESSION / "KRT19_withUTRs"
MOUSE_MITF_DIR = REGRESSION / "mouseMITF"
MIXED_FA = REGRESSION / "mixed_length_test" / "mixed_gc.fa"

HUMAN_GENOME_INDEX = "GCA_000001405.15_GRCh38_no_alt_analysis_set"

# Switches the L (low-complexity) mask off so results match the old references
NO_L_MASK = {"hp_threshold": 10, "di_threshold": 10}


# ── Helpers ──────────────────────────────────────────────────────────────────

def _oligo_sequences(oligos_file):
    """Sorted probe sequences from an *_oligos.txt (second-to-last column)."""
    with open(oligos_file) as f:
        return sorted(line.rstrip("\n").split("\t")[-2] for line in f if line.strip())


def _designed_sequences(result):
    """Sorted probe sequences from a ProbeDesignResult."""
    return sorted(p.sequence for p in result.probes)


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


# ── Bowtie Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def bowtie_env(monkeypatch):
    """Put the env's bin dir first on PATH, then skip if bowtie is still missing.

    masking.find_bowtie() runs `which bowtie` at call time, so PATH must be
    set before design_probes() runs. monkeypatch restores PATH afterwards.
    """
    env_bin = str(Path(sys.executable).parent)
    monkeypatch.setenv("PATH", env_bin + os.pathsep + os.environ.get("PATH", ""))
    if shutil.which("bowtie") is None:
        pytest.skip("bowtie not found on PATH")


@pytest.fixture
def human_indices(bowtie_env):
    """Directory holding the human pseudogene (.ebwt) and genome (.bt2) indices."""
    _require_indices(("humanPseudo", "ebwt"), (HUMAN_GENOME_INDEX, "bt2"))
    return str(INDEX_DIR)


@pytest.fixture
def mouse_indices(bowtie_env):
    """Directory holding the mouse pseudogene (.ebwt) and genome (.bt2) indices."""
    _require_indices(("mousePseudo", "ebwt"), ("mm10", "bt2"))
    return str(INDEX_DIR)


# ── CDKN1A: fixed length, repeat mask from file (no bowtie) ──────────────────

def test_cdkn1a_matches_reference():
    """32 probes with the repeat-mask file reproduce the reference exactly."""
    result = design_probes(
        str(CDKN1A_FA),
        n_probes=32,
        repeatmask_file=str(CDKN1A_REPEATMASKED),
        **NO_L_MASK,
    )
    assert len(result.probes) == 32
    assert _designed_sequences(result) == _oligo_sequences(CDKN1A_REF)


def test_cdkn1a_cli_writes_outputs_matching_reference(tmp_path):
    """The `design` CLI command writes _oligos.txt and _seq.txt with the same probes."""
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

    oligos_file = tmp_path / "cdkn1a_oligos.txt"
    seq_file = tmp_path / "cdkn1a_seq.txt"
    assert oligos_file.is_file()
    assert seq_file.is_file() and seq_file.stat().st_size > 0

    # index, start, GC%, Tm, dG, sequence, name
    rows = [line.split("\t") for line in oligos_file.read_text().splitlines()]
    assert len(rows) == 32
    assert all(len(row) == 7 for row in rows)
    # -o includes tmp_path, but probe names use only its last part
    assert [row[-1] for row in rows] == [f"cdkn1a_{i}" for i in range(1, 33)]
    assert _oligo_sequences(oligos_file) == _oligo_sequences(CDKN1A_REF)


def test_cdkn1a_default_thresholds_exclude_low_complexity_probes():
    """With the L mask on (defaults), no probe has a homopolymer >=5 or dinucleotide repeat >=3."""
    kwargs = {"n_probes": 32, "repeatmask_file": str(CDKN1A_REPEATMASKED)}

    result = design_probes(str(CDKN1A_FA), **kwargs)
    assert len(result.probes) == 32
    for probe in result.probes:
        assert not has_homopolymer(probe.sequence, 5), probe.sequence
        assert not has_dinucleotide_repeat(probe.sequence, 3), probe.sequence

    # Guard against a vacuous pass: with the L mask off the design does pick such probes
    unfiltered = design_probes(str(CDKN1A_FA), **kwargs, **NO_L_MASK)
    assert any(
        has_homopolymer(p.sequence, 5) or has_dinucleotide_repeat(p.sequence, 3)
        for p in unfiltered.probes
    )


# ── Mixed-length design ──────────────────────────────────────────────────────

def test_mixed_length_probes_satisfy_constraints():
    """18-22 nt probes are sorted, spaced, and inside the allowable dG range."""
    spacer = 2
    result = design_probes(
        str(MIXED_FA),
        n_probes=20,
        spacer_length=spacer,
        allowable_gibbs=(-26.0, -20.0),
        mixed_lengths=(18, 22),
    )
    probes = result.probes
    assert 0 < len(probes) <= 20

    assert all(18 <= p.length <= 22 for p in probes)
    assert all(len(p.sequence) == p.length for p in probes)
    assert len({p.length for p in probes}) > 1, "expected a genuine mix of lengths"
    assert all(-26.0 <= p.gibbs_fe <= -20.0 for p in probes)

    # Sorted, non-overlapping, and at least `spacer` nt between consecutive probes
    assert [p.position for p in probes] == sorted(p.position for p in probes)
    for left, right in zip(probes, probes[1:]):
        gap = right.position - (left.position + left.length)
        assert gap >= spacer, f"probes {left.index}/{right.index} are only {gap} nt apart"


def test_mixed_length_with_equal_bounds_equals_fixed_length():
    """mixed_lengths=(20, 20) is identical to the fixed-length path with oligo_length=20."""
    fixed = design_probes(str(MIXED_FA), n_probes=20, oligo_length=20)
    mixed = design_probes(str(MIXED_FA), n_probes=20, mixed_lengths=(20, 20))
    assert fixed.probes, "expected at least one probe"
    assert mixed.probes == fixed.probes
    assert mixed.score == pytest.approx(fixed.score)


# ── Bowtie masking (skipped when bowtie / indices are unavailable) ───────────

def test_krt19_pseudogene_and_genome_masks_match_reference(human_indices):
    """KRT19 with pseudogene + genome masking gives the 6-probe reference set."""
    result = design_probes(
        str(KRT19_DIR / "KRT19_withUTRs.fa"),
        n_probes=32,
        pseudogene_mask=True,
        genome_mask=True,
        index_dir=human_indices,
        **NO_L_MASK,
    )
    # Masking failures only warn (missing bowtie) or give all-zero hits (bad
    # index), so check that bowtie really found hits
    assert any(result.bowtie_pseudogene_hits or [])
    assert any(any(hits) for hits in (result.bowtie_genome_hits or {}).values())
    assert len(result.probes) == 6
    assert _designed_sequences(result) == _oligo_sequences(
        KRT19_DIR / "KRT19_withUTRs_oligos.txt"
    )


def test_mouse_mitf_pseudogene_and_genome_masks_match_reference(mouse_indices):
    """mouseMITF (24 probes, mouse, repeat-mask file) matches the reference.

    Both masks are needed: pseudogene masking alone gives 22/24 matching probes.
    """
    result = design_probes(
        str(MOUSE_MITF_DIR / "mouseMITF.fa"),
        n_probes=24,
        species="mouse",
        pseudogene_mask=True,
        genome_mask=True,
        repeatmask_file=str(MOUSE_MITF_DIR / "mouseMITF_repeatmasked.fa"),
        index_dir=mouse_indices,
        **NO_L_MASK,
    )
    assert any(result.bowtie_pseudogene_hits or [])
    assert any(any(hits) for hits in (result.bowtie_genome_hits or {}).values())
    assert len(result.probes) == 24
    assert _designed_sequences(result) == _oligo_sequences(
        MOUSE_MITF_DIR / "mouseMITF_oligos.txt"
    )
