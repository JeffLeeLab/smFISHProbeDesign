"""Tests for HCR split-initiator probe design module."""

import math
import os
import sys
import tempfile
from pathlib import Path

import pytest

# Ensure probedesign package is importable
_repo_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_repo_root / "src"))

from probedesign.hcr import (
    AMPLIFIER_TABLE,
    VALID_AMPLIFIERS,
    HALF_LENGTH,
    GAP_LENGTH,
    PAIR_BLOCK_LENGTH,
    HCRProbeHalf,
    HCRProbePair,
    HCRDesignResult,
    calculate_pair_badness,
    find_best_pairs,
    attach_initiators,
    design_hcr_probes,
)
from probedesign.thermodynamics import gibbs_rna_dna
from probedesign.masking import has_homopolymer, has_dinucleotide_repeat
from probedesign.sequence import reverse_complement
from probedesign.cli import main
from click.testing import CliRunner


# ── Test Fixtures ────────────────────────────────────────────────────────────

# Fixture 1: BOTH_STRICT — Both halves balanced GC, within strict range
BOTH_STRICT_SEQ = "atcgatcgatcgatcgatcgatcgaatatcgatcgatcgatcgatcgatcga"
#                  |------left 25nt------|gap|------right 25nt------|

# Fixture 2: ASYM_RESCUE — left balanced, right GC-rich (strict fail, lenient pass)
ASYM_RESCUE_SEQ = "atcgatcgatcgatcgatcgatcgaaagcgcgatcgcgatcgcgatcgcgat"

# Fixture 3: BOTH_LENIENT_FAIL — both halves GC-rich
BOTH_LENIENT_FAIL_SEQ = "gcgcgatcgcgatcgcgatcgcgatcgcgcgatcgcgatcgcgatcgcgatc"

# Fixture 4: HOMOPOLYMER_LEFT — left half contains AAAAA
HOMOPOLYMER_LEFT_SEQ = "atcgaaaaatcgatcgatcgatcgaatatcgatcgatcgatcgatcgatcga"

# Fixture 6: N_IN_GAP — N's at gap positions 25-26
N_IN_GAP_SEQ = "atcgatcgatcgatcgatcgatcgannatcgatcgatcgatcgatcgatcga"

# Fixture 7: TOO_SHORT — 51 nt, below 52-nt minimum
TOO_SHORT_SEQ = "atcgatcgatcgatcgatcgatcgaatatcgatcgatcgatcgatcgatcg"

# Fixture 8: MULTI_PAIR — 160 nt, enough for 2+ pairs
MULTI_PAIR_SEQ = ("atcgatcgatcgatcgatcgatcga" * 6) + "atcgatcgat"

# Fixture 9: MASK_ONE_HALF — same as BOTH_STRICT, used with synthetic mask
MASK_ONE_HALF_SEQ = BOTH_STRICT_SEQ

# Fixture 10: DINUCLEOTIDE_RIGHT — right half has ATATAT
DINUCLEOTIDE_RIGHT_SEQ = "atcgatcgatcgatcgatcgatcgaaatatatatcgatcgatcgatcgatcg"


def _write_temp_fasta(seq: str, name: str = "TEST") -> str:
    """Write a sequence to a temp FASTA file, return path."""
    f = tempfile.NamedTemporaryFile(
        mode='w', suffix='.fa', delete=False, prefix='hcr_test_'
    )
    f.write(f">{name}\n{seq}\n")
    f.close()
    return f.name


def _write_multi_fasta(entries: list) -> str:
    """Write multiple FASTA entries to a temp file, return path."""
    f = tempfile.NamedTemporaryFile(
        mode='w', suffix='.fa', delete=False, prefix='hcr_test_'
    )
    for name, seq in entries:
        f.write(f">{name}\n{seq}\n")
    f.close()
    return f.name


# ── Phase 1 Unit Tests ───────────────────────────────────────────────────────

class TestBadnessCalculation:
    """Tests for calculate_pair_badness()."""

    def test_badness_strict_both_pass(self):
        """Both halves of BOTH_STRICT have finite strict badness."""
        seq = BOTH_STRICT_SEQ
        left = seq[0:25]
        right = seq[27:52]
        gibbs_left = gibbs_rna_dna(left)
        gibbs_right = gibbs_rna_dna(right)
        # Both should be within strict range (-35, -27)
        assert -35 <= gibbs_left <= -27, f"Left ΔG {gibbs_left} outside strict range"
        assert -35 <= gibbs_right <= -27, f"Right ΔG {gibbs_right} outside strict range"

        pbad, orientations = calculate_pair_badness(seq, target_gibbs=-31.0, strict_range=(-35.0, -27.0))
        assert math.isfinite(pbad[0]), "Pair at position 0 should have finite badness"

    def test_badness_strict_gc_rich_fail(self):
        """Right half of ASYM_RESCUE fails strict but passes lenient."""
        seq = ASYM_RESCUE_SEQ
        right = seq[27:52]
        gibbs_right = gibbs_rna_dna(right)
        # GC-rich right half should have ΔG more negative than -35
        assert gibbs_right < -35, f"Right ΔG {gibbs_right} should be below -35"

        # Strict should fail
        pbad_s, _ = calculate_pair_badness(seq, target_gibbs=-31.0, strict_range=(-35.0, -27.0))
        assert not math.isfinite(pbad_s[0]), "Symmetric mode should reject this pair"

        # Lenient should pass
        pbad_a, _ = calculate_pair_badness(
            seq, target_gibbs=-31.0, strict_range=(-35.0, -27.0),
            asymmetric_gibbs=True, lenient_gibbs_min=-42.0
        )
        assert math.isfinite(pbad_a[0]), "Asymmetric mode should rescue this pair"

    def test_pair_badness_symmetric_pass(self):
        """BOTH_STRICT passes in symmetric mode."""
        pbad, orientations = calculate_pair_badness(BOTH_STRICT_SEQ)
        assert math.isfinite(pbad[0])
        assert orientations[0] == "both_strict"

    def test_pair_badness_symmetric_fail(self):
        """ASYM_RESCUE fails in symmetric mode."""
        pbad, _ = calculate_pair_badness(ASYM_RESCUE_SEQ)
        assert not math.isfinite(pbad[0])

    def test_pair_badness_asymmetric_rescue(self):
        """ASYM_RESCUE passes in asymmetric mode."""
        pbad, orientations = calculate_pair_badness(
            ASYM_RESCUE_SEQ, asymmetric_gibbs=True, lenient_gibbs_min=-42.0
        )
        assert math.isfinite(pbad[0])
        assert orientations[0] == "left_strict"

    def test_pair_badness_both_lenient_rejected(self):
        """BOTH_LENIENT_FAIL rejected even in asymmetric mode."""
        pbad, _ = calculate_pair_badness(
            BOTH_LENIENT_FAIL_SEQ,
            asymmetric_gibbs=True,
            lenient_gibbs_min=-42.0,
        )
        # Check both halves actually fail strict
        left = BOTH_LENIENT_FAIL_SEQ[0:25]
        right = BOTH_LENIENT_FAIL_SEQ[27:52]
        gibbs_left = gibbs_rna_dna(left)
        gibbs_right = gibbs_rna_dna(right)
        assert gibbs_left < -35 or gibbs_left > -27, f"Left should fail strict: {gibbs_left}"
        assert gibbs_right < -35 or gibbs_right > -27, f"Right should fail strict: {gibbs_right}"
        # Both fail strict → no valid orientation → inf
        assert not math.isfinite(pbad[0])

    def test_pair_badness_best_orientation(self):
        """Asymmetric mode picks the orientation with lower combined badness."""
        pbad, orientations = calculate_pair_badness(
            ASYM_RESCUE_SEQ, asymmetric_gibbs=True, lenient_gibbs_min=-42.0
        )
        assert math.isfinite(pbad[0])
        # The left half is balanced (strict pass), right is GC-rich (strict fail)
        # So orientation should be left_strict
        assert orientations[0] == "left_strict"


class TestJunctionAndGap:
    """Tests for junction marker and N-in-gap handling."""

    def test_junction_in_gap(self):
        """Junction marker '>' at gap position invalidates pair."""
        # SEQ1 is exactly 25nt, SEQ2 starts after junction
        seq1 = "atcgatcgatcgatcgatcgatcga"  # 25 nt
        seq2 = "aatcgatcgatcgatcgatcgatcgatcg"  # 28 nt
        # After concatenation with junction: seq1 + '>' + seq2
        # The '>' falls at position 25, which is in the gap
        concatenated = seq1 + ">" + seq2
        pbad, _ = calculate_pair_badness(concatenated)
        assert not math.isfinite(pbad[0]), "Pair spanning junction should be invalid"

    def test_n_in_gap_valid(self):
        """N's at gap positions don't invalidate the pair."""
        seq = N_IN_GAP_SEQ
        # N's are at positions 25-26 (gap), but halves are clean
        pbad, _ = calculate_pair_badness(seq)
        assert math.isfinite(pbad[0]), "N's in gap should not invalidate pair"


class TestSequenceLength:
    """Tests for sequence length validation."""

    def test_sequence_too_short(self):
        """Sequence shorter than 52nt returns empty result."""
        tmp = _write_temp_fasta(TOO_SHORT_SEQ, "TOO_SHORT")
        try:
            result = design_hcr_probes(tmp)
            assert len(result.pairs) == 0, "Should find 0 pairs for sequence < 52nt"
        finally:
            os.unlink(tmp)


class TestLowComplexity:
    """Tests for L-mask (homopolymer and dinucleotide repeat filtering)."""

    def test_homopolymer_masks_pair(self):
        """Homopolymer in left half masks the pair."""
        left = HOMOPOLYMER_LEFT_SEQ[0:25]
        assert has_homopolymer(left, 5), "Left half should have homopolymer"
        pbad, _ = calculate_pair_badness(HOMOPOLYMER_LEFT_SEQ, hp_threshold=5)
        assert not math.isfinite(pbad[0]), "Pair with homopolymer should be masked"

    def test_dinucleotide_masks_pair(self):
        """Dinucleotide repeat in right half masks the pair."""
        right = DINUCLEOTIDE_RIGHT_SEQ[27:52]
        assert has_dinucleotide_repeat(right, 3), "Right half should have dinucleotide repeat"
        pbad, _ = calculate_pair_badness(DINUCLEOTIDE_RIGHT_SEQ, di_threshold=3)
        assert not math.isfinite(pbad[0]), "Pair with dinucleotide repeat should be masked"

    def test_l_mask_always_symmetric(self):
        """L-mask rejects pair even in asymmetric ΔG mode."""
        pbad, _ = calculate_pair_badness(
            HOMOPOLYMER_LEFT_SEQ,
            asymmetric_gibbs=True,
            lenient_gibbs_min=-42.0,
            hp_threshold=5,
        )
        assert not math.isfinite(pbad[0]), "L-mask should be symmetric regardless of mode"


class TestBowtieFiltering:
    """Tests for bowtie mask application (using synthetic masks)."""

    def test_mask_symmetric_rejects(self):
        """Symmetric bowtie mode rejects when one half is masked."""
        seq = MASK_ONE_HALF_SEQ
        pb_mask = [1] * 25 + [0] * 27  # left half masked
        pbad, _ = calculate_pair_badness(seq, nuc_mask_pb=pb_mask, asymmetric_bowtie=False)
        assert not math.isfinite(pbad[0])

    def test_mask_asymmetric_rescues(self):
        """Asymmetric bowtie mode rescues when one half is clean."""
        seq = MASK_ONE_HALF_SEQ
        pb_mask = [1] * 25 + [0] * 27  # left half masked, right clean
        pbad, _ = calculate_pair_badness(seq, nuc_mask_pb=pb_mask, asymmetric_bowtie=True)
        assert math.isfinite(pbad[0])

    def test_mask_both_masked_rejected(self):
        """Both halves masked → rejected even in asymmetric mode."""
        seq = MASK_ONE_HALF_SEQ
        pb_mask = [1] * 52  # both halves masked
        pbad, _ = calculate_pair_badness(seq, nuc_mask_pb=pb_mask, asymmetric_bowtie=True)
        assert not math.isfinite(pbad[0])

    def test_bowtie_mask_spanning_gap(self):
        """Mask spanning both halves through gap rejects in symmetric mode."""
        seq = BOTH_STRICT_SEQ
        # Mask covering positions 20-31 (bleeding from left through gap into right)
        pb_mask = [0] * 20 + [1] * 12 + [0] * 20
        pbad, _ = calculate_pair_badness(seq, nuc_mask_pb=pb_mask, asymmetric_bowtie=False)
        assert not math.isfinite(pbad[0])


class TestRMask:
    """Tests for R-mask (repeat) handling."""

    def test_r_mask_n_in_gap_ignored(self):
        """R-mask at gap positions is ignored — pair remains valid."""
        seq = N_IN_GAP_SEQ
        # R-mask only at gap positions 25-26
        r_mask = [0] * 25 + [1, 1] + [0] * 25
        pbad, _ = calculate_pair_badness(seq, nuc_mask_r=r_mask)
        assert math.isfinite(pbad[0]), "R-mask at gap positions should be ignored"

    def test_r_mask_always_symmetric(self):
        """R-mask rejects even with asymmetric bowtie enabled."""
        seq = MASK_ONE_HALF_SEQ
        # R-mask on left half only
        r_mask = [1] * 25 + [0] * 27
        pbad, _ = calculate_pair_badness(
            seq, nuc_mask_r=r_mask, asymmetric_bowtie=True
        )
        assert not math.isfinite(pbad[0]), "R-mask should be always symmetric"


class TestDP:
    """Tests for find_best_pairs() dynamic programming."""

    def test_dp_no_overlap(self):
        """Two pairs have sufficient spacing."""
        seq = MULTI_PAIR_SEQ
        pbad, _ = calculate_pair_badness(seq)
        score, positions = find_best_pairs(pbad, pair_spacing=2, n_pairs=5)
        assert len(positions) >= 2, f"Should find at least 2 pairs, got {len(positions)}"
        for i in range(1, len(positions)):
            assert positions[i] - positions[i - 1] >= 54, \
                f"Pairs too close: {positions[i]} - {positions[i-1]} = {positions[i] - positions[i-1]}"

    def test_dp_maximise_pairs(self):
        """DP prefers more pairs over lower badness."""
        seq = MULTI_PAIR_SEQ
        pbad, _ = calculate_pair_badness(seq)
        score, positions = find_best_pairs(pbad, pair_spacing=2, n_pairs=10)
        # With 160nt and spacing 2: pairs at 0, 54, 108 (ends at 160) = 3 pairs
        # max = (len - PAIR_BLOCK_LENGTH) // (PAIR_BLOCK_LENGTH + spacing) + 1
        max_possible = (len(seq) - PAIR_BLOCK_LENGTH) // (PAIR_BLOCK_LENGTH + 2) + 1
        assert len(positions) == max_possible, \
            f"Should find {max_possible} pairs, got {len(positions)}"


class TestInitiatorAttachment:
    """Tests for attach_initiators()."""

    def test_initiator_b1(self):
        """B1 initiator attachment produces correct oligo structure."""
        binding_left = "atcgatcgatcgatcgatcgatcga"
        binding_right = "gctagctagctagctagctagctag"  # different from left
        p1, p2 = attach_initiators(binding_left, binding_right, "B1")

        rc_left = reverse_complement(binding_left)
        rc_right = reverse_complement(binding_right)
        amp = AMPLIFIER_TABLE["B1"]

        # P1: RC(right) + spacer_b + initiator_b — binds RIGHT half, init_b points inward
        expected_tail = (amp["spacer_b"] + amp["initiator_b"]).upper()
        assert p1 == rc_right.lower() + expected_tail
        assert p1[:25].islower(), "Binding region should be lowercase"
        assert p1[25:].isupper(), "Spacer+initiator should be uppercase"

        # P2: initiator_a + spacer_a + RC(left) — binds LEFT half, init_a points inward
        expected_head = (amp["initiator_a"] + amp["spacer_a"]).upper()
        assert p2 == expected_head + rc_left.lower()

    def test_initiator_all_amplifiers(self):
        """All amplifiers produce oligos with correct structure."""
        binding = "atcgatcgatcgatcgatcgatcga"
        for amp_id in VALID_AMPLIFIERS:
            p1, p2 = attach_initiators(binding, binding, amp_id)
            assert len(p1) > 25, f"{amp_id} P1 should be longer than 25"
            assert len(p2) > 25, f"{amp_id} P2 should be longer than 25"

    def test_ww_spacer_default_iupac(self):
        """WW spacer outputs as IUPAC by default."""
        binding = "atcgatcgatcgatcgatcgatcga"
        p1, p2 = attach_initiators(binding, binding, "B9")
        # B9 has WW spacers — should appear uppercase in output
        assert "WW" in p1, "P1 should contain WW spacer"
        assert "WW" in p2, "P2 should contain WW spacer"

    def test_ww_spacer_resolved(self):
        """--resolve-spacer replaces WW with explicit bases."""
        binding = "atcgatcgatcgatcgatcgatcga"
        p1, p2 = attach_initiators(binding, binding, "B9", resolve_spacer="AA")
        assert "WW" not in p1, "WW should be resolved"
        assert "WW" not in p2, "WW should be resolved"
        # The resolved spacer should be in the output
        assert "AA" in p1.upper()
        assert "AA" in p2.upper()


class TestOutputConventions:
    """Tests for output format conventions."""

    def test_output_case_convention(self):
        """Binding region lowercase, spacer+initiator uppercase."""
        binding = "atcgatcgatcgatcgatcgatcga"
        p1, p2 = attach_initiators(binding, binding, "B1")
        # First 25 chars of P1 = binding_rc (lowercase)
        assert p1[:25].islower()
        # Rest = spacer + initiator (uppercase)
        assert p1[25:].isupper()
        # P2: first part = initiator + spacer (uppercase)
        amp = AMPLIFIER_TABLE["B1"]
        head_len = len(amp["initiator_a"]) + len(amp["spacer_a"])
        assert p2[:head_len].isupper()
        # Last 25 chars = binding_rc (lowercase)
        assert p2[-25:].islower()

    def test_output_naming(self):
        """Probe names follow convention: GENE_HCR{amp}_{NN}."""
        tmp = _write_temp_fasta(MULTI_PAIR_SEQ, "TESTGENE")
        try:
            result = design_hcr_probes(tmp, n_pairs=2, amplifier="B1", output_name="TESTGENE")
            if result.pairs:
                pair = result.pairs[0]
                # Odd=P1, even=P2 naming is applied in output module, not here
                # But pair_index should be 1-based
                assert pair.pair_index == 1
        finally:
            os.unlink(tmp)


class TestEndToEnd:
    """End-to-end integration tests."""

    def test_design_both_strict_single_pair(self):
        """Design a single pair from BOTH_STRICT sequence."""
        tmp = _write_temp_fasta(BOTH_STRICT_SEQ, "BOTH_STRICT")
        try:
            result = design_hcr_probes(tmp, n_pairs=1)
            assert len(result.pairs) == 1
            pair = result.pairs[0]
            assert pair.left.length == 25
            assert pair.right.length == 25
            assert pair.amplifier == "B1"
            assert len(pair.left.oligo_seq) > 25
            assert len(pair.right.oligo_seq) > 25
        finally:
            os.unlink(tmp)

    def test_design_multi_pair(self):
        """Design multiple pairs from longer sequence."""
        tmp = _write_temp_fasta(MULTI_PAIR_SEQ, "MULTI")
        try:
            result = design_hcr_probes(tmp, n_pairs=5, pair_spacing=2)
            assert len(result.pairs) >= 2
            # Verify no overlap
            for i in range(1, len(result.pairs)):
                gap = result.pairs[i].pair_position - result.pairs[i-1].pair_position
                assert gap >= 54
        finally:
            os.unlink(tmp)

    def test_design_asymmetric_gibbs(self):
        """Asymmetric ΔG mode rescues pairs that symmetric rejects."""
        tmp = _write_temp_fasta(ASYM_RESCUE_SEQ, "ASYM")
        try:
            res_sym = design_hcr_probes(tmp, n_pairs=1, asymmetric_gibbs=False)
            res_asym = design_hcr_probes(tmp, n_pairs=1, asymmetric_gibbs=True, lenient_gibbs_min=-42.0)
            assert len(res_sym.pairs) == 0, "Symmetric should find 0 pairs"
            assert len(res_asym.pairs) == 1, "Asymmetric should find 1 pair"
            assert res_asym.asymmetric_mode == "asymmetric_gibbs"
        finally:
            os.unlink(tmp)

    def test_design_too_short(self):
        """Too-short sequence returns empty result."""
        tmp = _write_temp_fasta(TOO_SHORT_SEQ, "SHORT")
        try:
            result = design_hcr_probes(tmp)
            assert len(result.pairs) == 0
        finally:
            os.unlink(tmp)


# ── Phase 2 Output Tests ─────────────────────────────────────────────────────

from probedesign.hcr_output import (
    format_hcr_oligos,
    format_hcr_seq,
    format_hcr_hits,
    write_hcr_output_files,
)


class TestOutputFormatting:
    """Tests for HCR output module."""

    def _make_result(self):
        """Helper to create an HCRDesignResult for output tests."""
        tmp = _write_temp_fasta(MULTI_PAIR_SEQ, "GENE")
        try:
            result = design_hcr_probes(tmp, n_pairs=2, amplifier="B1", output_name="GENE")
        finally:
            os.unlink(tmp)
        return result

    def test_oligos_tsv_header(self):
        """Oligos file has correct header."""
        result = self._make_result()
        content = format_hcr_oligos(result)
        header = content.split("\n")[0]
        assert "pair" in header
        assert "half" in header
        assert "full_oligo" in header
        assert "name" in header

    def test_oligos_naming_convention(self):
        """Oligo names follow GENE_HCR{amp}_{NN} with odd=P1, even=P2."""
        result = self._make_result()
        content = format_hcr_oligos(result)
        lines = [l for l in content.strip().split("\n") if l and not l.startswith("pair")]
        if len(lines) >= 2:
            # First pair: P1=01, P2=02
            assert "GENE_HCRB1_01" in lines[0]
            assert "GENE_HCRB1_02" in lines[1]
        if len(lines) >= 4:
            # Second pair: P1=03, P2=04
            assert "GENE_HCRB1_03" in lines[2]
            assert "GENE_HCRB1_04" in lines[3]

    def test_oligos_strict_column(self):
        """Strict column shows yes/no."""
        tmp = _write_temp_fasta(ASYM_RESCUE_SEQ, "ASYM")
        try:
            result = design_hcr_probes(
                tmp, n_pairs=1, amplifier="B1",
                asymmetric_gibbs=True, lenient_gibbs_min=-42.0,
                output_name="ASYM",
            )
        finally:
            os.unlink(tmp)
        if result.pairs:
            content = format_hcr_oligos(result)
            lines = content.strip().split("\n")[1:]  # skip header
            # One half should be "no" (lenient)
            has_no = any("\tno\t" in line for line in lines)
            assert has_no, "One half should be lenient (strict=no)"

    def test_seq_file_contains_sequence(self):
        """Seq file contains the input sequence."""
        result = self._make_result()
        content = format_hcr_seq(result)
        assert "SEQUENCE:" in content
        assert "Position:" in content
        assert "Pairs:" in content

    def test_hits_file_structure(self):
        """Hits file has correct structure even without bowtie data."""
        result = self._make_result()
        content = format_hcr_hits(result)
        assert "Pair 1" in content
        assert "P1" in content
        assert "P2" in content

    def test_write_output_files(self):
        """write_hcr_output_files creates all three files."""
        result = self._make_result()
        with tempfile.TemporaryDirectory() as tmpdir:
            write_hcr_output_files(result, "GENE", output_dir=tmpdir)
            assert os.path.exists(os.path.join(tmpdir, "GENE_HCR_oligos.txt"))
            assert os.path.exists(os.path.join(tmpdir, "GENE_HCR_seq.txt"))
            assert os.path.exists(os.path.join(tmpdir, "GENE_HCR_hits.txt"))


# ── CLI Tests ────────────────────────────────────────────────────────────────

class TestCLI:
    """hcrdesign command-line interface."""

    def test_output_dir_not_in_probe_names(self, tmp_path):
        """-o with a directory writes files there; probe names use only its last part."""
        input_fa = _repo_root / "tests" / "regression" / "EIF1_CDS_HCR" / "EIF1_Exons.fasta"
        result = CliRunner().invoke(main, [
            "hcrdesign", str(input_fa), "-p", "5", "-o", str(tmp_path / "EIF1"), "--quiet",
        ])
        assert result.exit_code == 0, f"{result.output}\n{result.exception!r}"

        rows = (tmp_path / "EIF1_HCR_oligos.txt").read_text().splitlines()[1:]
        names = [row.split("\t")[-1] for row in rows]
        assert names, "expected at least one probe pair"
        assert all(name.startswith("EIF1_HCRB1_") for name in names), names


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
