"""HCR split-initiator probe design module.

Designs HCR v3 split-initiator probe pairs for RNA FISH. Each pair consists
of two 25-nt probes separated by a 2-nt gap on the target, with half-initiator
sequences attached for HCR signal amplification.

Key design principles:
- Each 25-nt half hybridises independently (not a single 52-nt duplex)
- Asymmetric leniency: one half can have relaxed ΔG/bowtie filtering
- DP operates on 52-nt pair blocks as atomic units
"""

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .thermodynamics import gibbs_rna_dna, tm_rna_dna
from .sequence import reverse_complement, percent_gc, has_invalid_chars, clean_sequence
from .fasta import read_fasta, sequences_to_single_string, get_basename
from .masking import has_homopolymer, has_dinucleotide_repeat
from .provenance import new_run_info, bowtie_mask_record


# ── HCR amplifier initiator table ───────────────────────────────────────────
# Source: Choi et al. HCR v3.0 amplifiers (from HCR_design.R)
# Case convention preserved from R script: lowercase = RNA-binding context
# Internally normalised to uppercase for sequence operations.

AMPLIFIER_TABLE: Dict[str, Dict[str, str]] = {
    "B1":  {"initiator_a": "gAggAgggCAgCAAACgg",  "initiator_b": "gAAgAgTCTTCCTTTACg",  "spacer_a": "AA", "spacer_b": "TA"},
    "B2":  {"initiator_a": "CCTCgTAAATCCTCATCA",  "initiator_b": "ATCATCCAgTAAACCgCC",  "spacer_a": "AA", "spacer_b": "AA"},
    "B3":  {"initiator_a": "gTCCCTgCCTCTATATCT",  "initiator_b": "CCACTCAACTTTAACCCg",  "spacer_a": "TT", "spacer_b": "TT"},
    "B4":  {"initiator_a": "CCTCAACCTACCTCCAAC",  "initiator_b": "TCTCACCATATTCgCTTC",  "spacer_a": "AA", "spacer_b": "AT"},
    "B5":  {"initiator_a": "CTCACTCCCAATCTCTAT",  "initiator_b": "CTACCCTACAAATCCAAT",  "spacer_a": "AA", "spacer_b": "AA"},
    "B7":  {"initiator_a": "CTTCAACCTCCACCTACC",  "initiator_b": "TCCAATCCCTACCCTCAC",  "spacer_a": "WW", "spacer_b": "WW"},
    "B9":  {"initiator_a": "CACGTATCTACTCCACTC",  "initiator_b": "TCAGCACACTCCCAACCC",  "spacer_a": "WW", "spacer_b": "WW"},
    "B10": {"initiator_a": "CCTCAAGATACTCCTCTA",  "initiator_b": "CCTACTCGACTACCCTAG",  "spacer_a": "WW", "spacer_b": "WW"},
    "B11": {"initiator_a": "CGCTTAGATATCACTCCT",  "initiator_b": "ACGTCGACCACACTCATC",  "spacer_a": "WW", "spacer_b": "WW"},
    "B13": {"initiator_a": "AGGTAACGCCTTCCTGCT",  "initiator_b": "TTATGCTCAACATACAAC",  "spacer_a": "WW", "spacer_b": "WW"},
    "B14": {"initiator_a": "AATGTCAATAGCGAGCGA",  "initiator_b": "CCCTATATTTCTGCACAG",  "spacer_a": "WW", "spacer_b": "WW"},
    "B15": {"initiator_a": "CAGATTAACACACCACAA",  "initiator_b": "GGTATCTCGAACACTCTC",  "spacer_a": "WW", "spacer_b": "WW"},
    "B17": {"initiator_a": "CGATTGTTTGTTGTGGAC",  "initiator_b": "GCATGCTAATCGGATGAG",  "spacer_a": "WW", "spacer_b": "WW"},
}

VALID_AMPLIFIERS = list(AMPLIFIER_TABLE.keys())

# Constants
HALF_LENGTH = 25
GAP_LENGTH = 2
PAIR_BLOCK_LENGTH = HALF_LENGTH * 2 + GAP_LENGTH  # 52


# ── Data structures ─────────────────────────────────────────────────────────

@dataclass
class HCRProbeHalf:
    """One half of an HCR split-initiator probe pair."""
    position: int          # Start position on target (0-indexed)
    length: int            # Always 25
    binding_seq: str       # 25-nt target-binding sequence (sense strand)
    binding_rc: str        # Reverse complement (what is synthesised)
    gc_percent: float
    tm: float
    gibbs_fe: float
    is_strict: bool        # True if this half passed strict filtering
    oligo_seq: str         # Full synthesised oligo (binding_rc + spacer + initiator)


@dataclass
class HCRProbePair:
    """A complete HCR split-initiator probe pair."""
    pair_index: int         # 1-based pair number
    left: HCRProbeHalf      # Binds 5' portion of target window (P2 in output — carries initiator_a)
    right: HCRProbeHalf     # Binds 3' portion of target window (P1 in output — carries initiator_b)
    pair_position: int      # Start of the 52-nt window on target
    combined_badness: float
    amplifier: str          # e.g. "B1"


@dataclass
class HCRDesignResult:
    """Result of HCR probe design."""
    pairs: List[HCRProbePair]
    score: float
    input_sequence: str
    template_name: str
    amplifier: str
    mask_strings: List[str] = field(default_factory=list)
    asymmetric_mode: str = "symmetric"  # "symmetric", "asymmetric_gibbs", "asymmetric_bowtie", "asymmetric_both"
    bowtie_pseudogene_hits: Optional[List[int]] = None
    bowtie_genome_hits: Optional[dict] = None
    bowtie_pseudogene_raw: Optional[str] = None
    bowtie_genome_raw: Optional[str] = None
    # Parallel to mask_strings: the mask letter of each line (R, P, B)
    mask_kinds: List[str] = field(default_factory=list)
    # Parameters, input and mask provenance; see provenance.py
    run_info: dict = field(default_factory=dict)


# ── Core functions ───────────────────────────────────────────────────────────

def calculate_pair_badness(
    seq: str,
    target_gibbs: float = -31.0,
    strict_range: Tuple[float, float] = (-35.0, -27.0),
    asymmetric_gibbs: bool = False,
    lenient_gibbs_min: float = -42.0,
    hp_threshold: int = 5,
    di_threshold: int = 3,
    nuc_mask_r: Optional[List[int]] = None,
    nuc_mask_pb: Optional[List[int]] = None,
    asymmetric_bowtie: bool = False,
) -> Tuple[List[float], List[Optional[str]]]:
    """Calculate pair-level badness for all valid positions.

    Returns badness and orientation for each pair start position p, where
    the pair occupies seq[p:p+52] (left=p:p+25, gap=p+25:p+27, right=p+27:p+52).

    Args:
        seq: Full concatenated sequence (lowercase, may contain '>' junctions)
        target_gibbs: Target ΔG per 25-nt half (kcal/mol)
        strict_range: (min, max) strict ΔG range
        asymmetric_gibbs: Enable asymmetric ΔG leniency
        lenient_gibbs_min: Floor ΔG for lenient half
        hp_threshold: Homopolymer run threshold for L-mask
        di_threshold: Dinucleotide repeat threshold for L-mask
        nuc_mask_r: Nucleotide-level repeat mask (always symmetric)
        nuc_mask_pb: Nucleotide-level pseudogene+genome mask (potentially asymmetric)
        asymmetric_bowtie: Enable asymmetric bowtie filtering

    Returns:
        Tuple of (pair_badness, orientations) where:
        - pair_badness[p] = combined badness for pair at position p (inf if invalid)
        - orientations[p] = "left_strict" or "right_strict" or None
    """
    from .core import calculate_badness

    n = len(seq)
    n_pairs = n - PAIR_BLOCK_LENGTH + 1
    if n_pairs <= 0:
        return [], []

    # Step 0: compute 25-nt badness arrays
    badness_strict = calculate_badness(seq, HALF_LENGTH, target_gibbs, strict_range)
    lenient_range = (lenient_gibbs_min, strict_range[1])
    if asymmetric_gibbs:
        badness_lenient = calculate_badness(seq, HALF_LENGTH, target_gibbs, lenient_range)
    else:
        badness_lenient = None

    # Apply L-mask to badness arrays (always symmetric)
    for i in range(len(badness_strict)):
        if badness_strict[i] == float('inf'):
            continue
        sub = seq[i:i + HALF_LENGTH]
        if has_homopolymer(sub, hp_threshold) or has_dinucleotide_repeat(sub, di_threshold):
            badness_strict[i] = float('inf')
            if badness_lenient is not None and i < len(badness_lenient):
                badness_lenient[i] = float('inf')

    if badness_lenient is not None:
        for i in range(len(badness_lenient)):
            if badness_lenient[i] == float('inf'):
                continue
            if badness_strict[i] != float('inf'):
                continue  # already handled above
            sub = seq[i:i + HALF_LENGTH]
            if has_homopolymer(sub, hp_threshold) or has_dinucleotide_repeat(sub, di_threshold):
                badness_lenient[i] = float('inf')

    pair_badness = [float('inf')] * n_pairs
    orientations: List[Optional[str]] = [None] * n_pairs

    for p in range(n_pairs):
        left_pos = p
        right_pos = p + HALF_LENGTH + GAP_LENGTH  # p + 27

        # R-mask: always symmetric — if either half is in a repeat, pair is dead
        if nuc_mask_r is not None:
            if any(nuc_mask_r[p:p + HALF_LENGTH]) or any(nuc_mask_r[right_pos:right_pos + HALF_LENGTH]):
                continue  # pair_badness[p] stays inf

        # Junction marker in gap
        if '>' in seq[p + HALF_LENGTH:p + HALF_LENGTH + GAP_LENGTH]:
            continue  # pair_badness[p] stays inf

        # Bounds check for badness arrays
        if left_pos >= len(badness_strict) or right_pos >= len(badness_strict):
            continue

        # ΔG filter (independent)
        if asymmetric_gibbs and badness_lenient is not None:
            gibbs_a = badness_strict[left_pos] + badness_lenient[right_pos]   # left strict
            gibbs_b = badness_lenient[left_pos] + badness_strict[right_pos]   # right strict
            if gibbs_a <= gibbs_b:
                gibbs_score = gibbs_a
                orientation = "left_strict"
            else:
                gibbs_score = gibbs_b
                orientation = "right_strict"
        else:
            gibbs_score = badness_strict[left_pos] + badness_strict[right_pos]
            orientation = "both_strict"

        if not math.isfinite(gibbs_score):
            continue  # pair_badness[p] stays inf

        # Bowtie filter (independent of ΔG orientation)
        if nuc_mask_pb is not None:
            left_offmasked = any(nuc_mask_pb[p:p + HALF_LENGTH])
            right_offmasked = any(nuc_mask_pb[right_pos:right_pos + HALF_LENGTH])
            if asymmetric_bowtie:
                bowtie_ok = (not left_offmasked) or (not right_offmasked)
            else:
                bowtie_ok = (not left_offmasked) and (not right_offmasked)
            if not bowtie_ok:
                continue  # pair_badness[p] stays inf

        pair_badness[p] = gibbs_score
        orientations[p] = orientation

    return pair_badness, orientations


def find_best_pairs(
    pair_badness: List[float],
    pair_spacing: int = 2,
    n_pairs: int = 30,
) -> Tuple[float, List[int]]:
    """Find optimal pair positions using dynamic programming.

    Maximises pair count first, then minimises average badness.

    Args:
        pair_badness: Per-position pair badness (inf = invalid)
        pair_spacing: Minimum gap between consecutive 52-nt blocks
        n_pairs: Maximum number of pairs to find

    Returns:
        Tuple of (score, positions) where positions is a list of pair start positions.
        Returns (inf, []) if no valid pairs exist.
    """
    n = len(pair_badness)
    if n == 0:
        return float('inf'), []

    block_step = PAIR_BLOCK_LENGTH + pair_spacing  # 52 + spacing

    # dp[p][k] = best average badness for k+1 pairs with last pair at or before p
    dp = [[float('inf')] * n_pairs for _ in range(n)]
    tracker = [[None] * n_pairs for _ in range(n)]

    # Initialize
    if math.isfinite(pair_badness[0]):
        dp[0][0] = pair_badness[0]
        tracker[0][0] = 0

    for p in range(1, n):
        # Propagate from previous position
        for k in range(n_pairs):
            if dp[p - 1][k] < dp[p][k]:
                dp[p][k] = dp[p - 1][k]
                tracker[p][k] = tracker[p - 1][k]

        # Try placing pair at position p
        if not math.isfinite(pair_badness[p]):
            continue

        for k in range(n_pairs):
            if k == 0:
                score = pair_badness[p]
            else:
                prev_p = p - block_step
                if prev_p < 0 or tracker[prev_p][k - 1] is None:
                    continue
                # Running average
                score = (dp[prev_p][k - 1] * k + pair_badness[p]) / (k + 1)

            if score < dp[p][k]:
                dp[p][k] = score
                tracker[p][k] = p

    # Find best solution: maximise k first, then minimise score
    best_k = -1
    best_score = float('inf')
    for k in range(n_pairs - 1, -1, -1):
        if tracker[n - 1][k] is not None and math.isfinite(dp[n - 1][k]):
            best_k = k
            best_score = dp[n - 1][k]
            break

    if best_k < 0:
        return float('inf'), []

    # Backtrack
    positions = []
    p = n - 1
    k = best_k
    while k >= 0:
        pos = tracker[p][k]
        positions.append(pos)
        p = pos - block_step
        k -= 1

    positions.reverse()
    return best_score, positions


def attach_initiators(
    binding_seq_left: str,
    binding_seq_right: str,
    amplifier: str = "B1",
    resolve_spacer: Optional[str] = None,
) -> Tuple[str, str]:
    """Attach initiator sequences to binding regions.

    P1 (odd probes, binds RIGHT half):  5'—RC(right)—spacer_b—initiator_b—3'
    P2 (even probes, binds LEFT half):  5'—initiator_a—spacer_a—RC(left)—3'

    Both initiator tails point INWARD toward the 2-nt gap for cooperative
    HCR triggering. Case convention: binding region lowercase, spacer+initiator UPPERCASE.

    Args:
        binding_seq_left: 25-nt sense strand of left half
        binding_seq_right: 25-nt sense strand of right half
        amplifier: Amplifier ID (e.g. "B1")
        resolve_spacer: If provided, replace WW spacer with this (e.g. "AA")

    Returns:
        Tuple of (p1_oligo, p2_oligo) with case convention applied.
    """
    amp = AMPLIFIER_TABLE[amplifier]
    spacer_a = amp["spacer_a"]
    spacer_b = amp["spacer_b"]

    if resolve_spacer is not None:
        if spacer_a == "WW":
            spacer_a = resolve_spacer
        if spacer_b == "WW":
            spacer_b = resolve_spacer

    rc_left = reverse_complement(binding_seq_left)
    rc_right = reverse_complement(binding_seq_right)

    # P1 (odd probes): 5'—RC(right)—spacer_b—initiator_b—3'
    # Binds RIGHT half of target; initiator_b at 3' end points inward toward gap
    p1_oligo = rc_right.lower() + (spacer_b + amp["initiator_b"]).upper()

    # P2 (even probes): 5'—initiator_a—spacer_a—RC(left)—3'
    # Binds LEFT half of target; initiator_a at 5' end points inward toward gap
    p2_oligo = (amp["initiator_a"] + spacer_a).upper() + rc_left.lower()

    return p1_oligo, p2_oligo


def design_hcr_probes(
    input_file: str,
    n_pairs: int = 30,
    amplifier: str = "B1",
    target_gibbs: float = -31.0,
    strict_range: Tuple[float, float] = (-35.0, -27.0),
    asymmetric_gibbs: bool = False,
    lenient_gibbs_min: float = -42.0,
    asymmetric_bowtie: bool = False,
    pair_spacing: int = 2,
    species: str = "human",
    pseudogene_mask: bool = False,
    genome_mask: bool = False,
    index_dir: Optional[str] = None,
    repeatmask_file: Optional[str] = None,
    hp_threshold: int = 5,
    di_threshold: int = 3,
    resolve_spacer: Optional[str] = None,
    output_name: Optional[str] = None,
    save_bowtie_raw: bool = False,
) -> HCRDesignResult:
    """Design HCR split-initiator probe pairs for a target sequence.

    Main entry point for HCR probe design.

    Args:
        input_file: Path to input FASTA file
        n_pairs: Number of probe pairs to design
        amplifier: HCR amplifier (B1-B5, B7, B9-B11, B13-B15, B17)
        target_gibbs: Target ΔG per 25-nt half (kcal/mol)
        strict_range: Strict ΔG range (min, max)
        asymmetric_gibbs: Enable asymmetric ΔG leniency
        lenient_gibbs_min: Floor ΔG for lenient half
        asymmetric_bowtie: Enable asymmetric bowtie filtering
        pair_spacing: Minimum gap between consecutive 52-nt blocks
        species: Species for masking databases
        pseudogene_mask: Enable pseudogene masking
        genome_mask: Enable genome masking
        index_dir: Directory containing bowtie indexes
        repeatmask_file: Path to RepeatMasker output file
        hp_threshold: Homopolymer run threshold
        di_threshold: Dinucleotide repeat threshold
        resolve_spacer: Resolve WW spacer to explicit bases
        output_name: Output file prefix
        save_bowtie_raw: Save raw bowtie output

    Returns:
        HCRDesignResult containing designed probe pairs
    """
    if amplifier not in AMPLIFIER_TABLE:
        raise ValueError(f"Unknown amplifier: {amplifier}. Valid: {VALID_AMPLIFIERS}")

    # Read input sequence
    headers, seqs = read_fasta(input_file)
    full_seq = sequences_to_single_string(seqs, mark_junctions=True)
    seq = clean_sequence(full_seq)

    if output_name is None:
        output_name = get_basename(input_file)

    run_info = new_run_info(
        "HCR", input_file, len(seqs), seq,
        params={
            "n_pairs": n_pairs,
            "amplifier": amplifier,
            "target_gibbs": target_gibbs,
            "strict_range": tuple(strict_range),
            "asymmetric_gibbs": asymmetric_gibbs,
            "lenient_gibbs_min": lenient_gibbs_min,
            "asymmetric_bowtie": asymmetric_bowtie,
            "pair_spacing": pair_spacing,
            "species": species,
            "pseudogene_mask": pseudogene_mask,
            "genome_mask": genome_mask,
            "hp_threshold": hp_threshold,
            "di_threshold": di_threshold,
            "resolve_spacer": resolve_spacer,
            "save_bowtie_raw": save_bowtie_raw,
        },
        repeatmask_file=repeatmask_file,
        index_dir=index_dir,
    )

    # Check minimum sequence length
    if len(seq) < PAIR_BLOCK_LENGTH:
        print(f"Sequence too short for HCR probe design (minimum {PAIR_BLOCK_LENGTH} nt, got {len(seq)} nt)")
        return HCRDesignResult(
            pairs=[], score=float('inf'), input_sequence=seq,
            template_name=output_name, amplifier=amplifier,
            run_info=run_info,
        )

    # Determine asymmetric mode string
    if asymmetric_gibbs and asymmetric_bowtie:
        asym_mode = "asymmetric_both"
    elif asymmetric_gibbs:
        asym_mode = "asymmetric_gibbs"
    elif asymmetric_bowtie:
        asym_mode = "asymmetric_bowtie"
    else:
        asym_mode = "symmetric"

    # Build masks
    mask_strings = []
    mask_kinds = []
    nuc_mask_r = None  # repeat mask
    nuc_mask_pb = None  # pseudogene + genome mask
    has_n_masking = 'n' in seq

    pseudo_hits = None
    genome_hits_dict = None
    pseudo_raw_text = None
    genome_raw_text = None

    # R-mask from repeatmask file or N's in input
    if repeatmask_file:
        _, rm_seqs = read_fasta(repeatmask_file)
        rm_full_seq = sequences_to_single_string(rm_seqs, mark_junctions=True)
        rm_seq = clean_sequence(rm_full_seq)
        nuc_mask_r = [1 if rm_seq[i].lower() == 'n' else 0 for i in range(len(rm_seq))]
        rstr = "".join("R" if nuc_mask_r[i] else seq[i] for i in range(len(seq)))
        mask_strings.append(rstr)
        mask_kinds.append("R")
        print(f"Repeat masking (from file): {sum(nuc_mask_r)} positions masked")
    elif has_n_masking:
        nuc_mask_r = [1 if seq[i] == 'n' else 0 for i in range(len(seq))]
        rstr = "".join("R" if nuc_mask_r[i] else seq[i] for i in range(len(seq)))
        mask_strings.append(rstr)
        mask_kinds.append("R")
        print(f"Repeat masking (from N's): {sum(nuc_mask_r)} positions masked")

    # P/B masks from bowtie
    if pseudogene_mask or genome_mask:
        from .masking import (
            pseudogene_mask as get_pseudogene_mask,
            genome_mask as get_genome_mask,
        )
        idx_dir = Path(index_dir) if index_dir else None
        nuc_mask_pb = [0] * len(seq)

        if pseudogene_mask:
            try:
                pmask, pseudo_hits, pseudo_raw = get_pseudogene_mask(seq, species, idx_dir)
                if not save_bowtie_raw:
                    pseudo_raw = None
                pseudo_raw_text = pseudo_raw
                for i, v in enumerate(pmask):
                    nuc_mask_pb[i] |= v
                pstr = "".join("P" if pmask[i] else seq[i] for i in range(len(seq)))
                mask_strings.append(pstr)
                mask_kinds.append("P")
                run_info["masks"].append(bowtie_mask_record("P", species, idx_dir))
                print(f"Pseudogene masking: {sum(pmask)} positions masked")
            except Exception as e:
                run_info["masks"].append(bowtie_mask_record("P", species, idx_dir, error=e))
                print(f"Warning: Pseudogene masking failed: {e}")

        if genome_mask:
            try:
                gmask, genome_hits_dict, genome_raw = get_genome_mask(seq, species, idx_dir)
                if not save_bowtie_raw:
                    genome_raw = None
                genome_raw_text = genome_raw
                for i, v in enumerate(gmask):
                    nuc_mask_pb[i] |= v
                gstr = "".join("B" if gmask[i] else seq[i] for i in range(len(seq)))
                mask_strings.append(gstr)
                mask_kinds.append("B")
                run_info["masks"].append(bowtie_mask_record("B", species, idx_dir))
                print(f"Genome masking: {sum(gmask)} positions masked")
            except Exception as e:
                run_info["masks"].append(bowtie_mask_record("B", species, idx_dir, error=e))
                print(f"Warning: Genome masking failed: {e}")

    # Calculate pair badness
    pbad, orientations = calculate_pair_badness(
        seq=seq,
        target_gibbs=target_gibbs,
        strict_range=strict_range,
        asymmetric_gibbs=asymmetric_gibbs,
        lenient_gibbs_min=lenient_gibbs_min,
        hp_threshold=hp_threshold,
        di_threshold=di_threshold,
        nuc_mask_r=nuc_mask_r,
        nuc_mask_pb=nuc_mask_pb,
        asymmetric_bowtie=asymmetric_bowtie,
    )

    # Report filter stats
    n_total = len(pbad)
    n_valid = sum(1 for b in pbad if math.isfinite(b))
    print(f"Valid pair positions: {n_valid}/{n_total}")

    if n_valid == 0:
        print("No valid pair positions found. Try adjusting ΔG range or enabling asymmetric mode.")
        return HCRDesignResult(
            pairs=[], score=float('inf'), input_sequence=seq,
            template_name=output_name, amplifier=amplifier,
            mask_strings=mask_strings, asymmetric_mode=asym_mode,
            mask_kinds=mask_kinds, run_info=run_info,
            bowtie_pseudogene_hits=pseudo_hits,
            bowtie_genome_hits=genome_hits_dict,
            bowtie_pseudogene_raw=pseudo_raw_text,
            bowtie_genome_raw=genome_raw_text,
        )

    # DP to find optimal pair positions
    score, positions = find_best_pairs(pbad, pair_spacing, n_pairs)

    if not positions:
        print("DP found no valid placement.")
        return HCRDesignResult(
            pairs=[], score=float('inf'), input_sequence=seq,
            template_name=output_name, amplifier=amplifier,
            mask_strings=mask_strings, asymmetric_mode=asym_mode,
            mask_kinds=mask_kinds, run_info=run_info,
        )

    if len(positions) < n_pairs:
        print(f"Warning: Only {len(positions)} pairs found (requested {n_pairs})")

    # Build probe pairs
    pairs = []
    for i, p in enumerate(positions):
        left_seq = seq[p:p + HALF_LENGTH]
        right_pos = p + HALF_LENGTH + GAP_LENGTH
        right_seq = seq[right_pos:right_pos + HALF_LENGTH]

        orientation = orientations[p]
        left_strict = orientation in ("left_strict", "both_strict")
        right_strict = orientation in ("right_strict", "both_strict")

        p1_oligo, p2_oligo = attach_initiators(
            left_seq, right_seq, amplifier, resolve_spacer
        )

        left_half = HCRProbeHalf(
            position=p,
            length=HALF_LENGTH,
            binding_seq=left_seq,
            binding_rc=reverse_complement(left_seq),
            gc_percent=round(percent_gc(left_seq) * 100),
            tm=round(tm_rna_dna(left_seq), 1),
            gibbs_fe=round(gibbs_rna_dna(left_seq), 1),
            is_strict=left_strict,
            oligo_seq=p2_oligo,   # P2 binds left half (carries initiator_a)
        )

        right_half = HCRProbeHalf(
            position=right_pos,
            length=HALF_LENGTH,
            binding_seq=right_seq,
            binding_rc=reverse_complement(right_seq),
            gc_percent=round(percent_gc(right_seq) * 100),
            tm=round(tm_rna_dna(right_seq), 1),
            gibbs_fe=round(gibbs_rna_dna(right_seq), 1),
            is_strict=right_strict,
            oligo_seq=p1_oligo,   # P1 binds right half (carries initiator_b)
        )

        pair = HCRProbePair(
            pair_index=i + 1,
            left=left_half,
            right=right_half,
            pair_position=p,
            combined_badness=pbad[p],
            amplifier=amplifier,
        )
        pairs.append(pair)

    print(f"Designed {len(pairs)} HCR probe pairs (amplifier {amplifier}, score {score:.4f})")

    return HCRDesignResult(
        pairs=pairs,
        score=score,
        input_sequence=seq,
        template_name=output_name,
        amplifier=amplifier,
        mask_strings=mask_strings,
        mask_kinds=mask_kinds,
        run_info=run_info,
        asymmetric_mode=asym_mode,
        bowtie_pseudogene_hits=pseudo_hits,
        bowtie_genome_hits=genome_hits_dict,
        bowtie_pseudogene_raw=pseudo_raw_text,
        bowtie_genome_raw=genome_raw_text,
    )
