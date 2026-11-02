"""Core probe design algorithms.

This module implements the main probe design logic:
1. Calculate "badness" scores for each position based on thermodynamic properties
2. Use dynamic programming to find optimal probe placements
3. Generate probe designs from the results
"""

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

from .thermodynamics import gibbs_rna_dna, tm_rna_dna
from .sequence import (
    reverse_complement,
    percent_gc,
    has_invalid_chars,
    clean_sequence,
)
from .fasta import read_fasta, sequences_to_single_string, get_basename
from .provenance import new_run_info, bowtie_mask_record


@dataclass
class Probe:
    """Represents a designed oligonucleotide probe."""
    index: int           # 1-based probe number
    position: int        # 0-based start position in template
    length: int          # Probe length in bp
    sequence: str        # Probe sequence (reverse complement of template)
    gc_percent: float    # GC content as percentage (0-100)
    tm: float            # Melting temperature (°C)
    gibbs_fe: float      # Gibbs free energy (kcal/mol)
    name: str            # Probe identifier


@dataclass
class ProbeDesignResult:
    """Result of probe design."""
    probes: List[Probe]
    score: float
    input_sequence: str
    template_name: str
    mask: Optional[List[int]] = None
    mask_strings: List[str] = field(default_factory=list)
    bowtie_pseudogene_hits: Optional[List[int]] = None
    bowtie_genome_hits: Optional[dict] = None
    bowtie_pseudogene_raw: Optional[str] = None
    bowtie_genome_raw: Optional[str] = None
    # Parallel to mask_strings: the mask letter of each line (R, P, B, F, L)
    mask_kinds: List[str] = field(default_factory=list)
    # Parameters, input and mask provenance; see provenance.py
    run_info: dict = field(default_factory=dict)


def calculate_badness(
    seq: str,
    oligo_len: int,
    target_gibbs: float,
    allowable_range: Tuple[float, float],
) -> List[float]:
    """Calculate badness score for each position in the sequence.

    The badness is the squared difference from the target Gibbs free energy.
    Positions with invalid characters or Gibbs FE outside the allowable range
    get infinite badness.

    Args:
        seq: Input sequence (lowercase)
        oligo_len: Length of oligonucleotides
        target_gibbs: Target Gibbs free energy (kcal/mol)
        allowable_range: Tuple of (min_gibbs, max_gibbs) for valid probes

    Returns:
        List of badness scores, one per position. Length = len(seq) - oligo_len + 1
    """
    min_gibbs, max_gibbs = sorted(allowable_range)
    goodlen = len(seq) - oligo_len + 1
    badness = []

    for i in range(goodlen):
        oligo = seq[i:i + oligo_len]

        # Check for invalid/masked characters
        if has_invalid_chars(oligo):
            badness.append(float('inf'))
            continue

        # Calculate Gibbs free energy
        try:
            gibbs = gibbs_rna_dna(oligo)
        except KeyError:
            # Invalid dinucleotide (shouldn't happen with clean sequence)
            badness.append(float('inf'))
            continue

        # Check if within allowable range
        if gibbs < min_gibbs or gibbs > max_gibbs:
            badness.append(float('inf'))
            continue

        # Badness is squared distance from target
        badness.append((gibbs - target_gibbs) ** 2)

    return badness


def calculate_badness_mixed(
    seq: str,
    min_len: int,
    max_len: int,
    target_gibbs: float,
    allowable_range: Tuple[float, float],
) -> List[List[float]]:
    """Calculate badness scores for mixed-length probes.

    For each position, calculates badness for every probe length in
    [min_len, max_len]. This enables the mixed-length DP to consider
    all valid lengths at each position.

    Args:
        seq: Input sequence (lowercase)
        min_len: Minimum probe length
        max_len: Maximum probe length
        target_gibbs: Target Gibbs free energy (kcal/mol)
        allowable_range: Tuple of (min_gibbs, max_gibbs) for valid probes

    Returns:
        2D list: badness_2d[pos][len_idx] where len_idx = L - min_len.
        Length of outer list = len(seq) - min_len + 1.
        Length of inner list = max_len - min_len + 1.
    """
    min_gibbs, max_gibbs = sorted(allowable_range)
    max_goodlen = len(seq) - min_len + 1
    n_lengths = max_len - min_len + 1
    badness_2d = []

    for i in range(max_goodlen):
        row = []
        for L_idx in range(n_lengths):
            L = min_len + L_idx

            # Check if probe would extend past sequence end
            if i + L > len(seq):
                row.append(float('inf'))
                continue

            oligo = seq[i:i + L]

            # Check for invalid/masked characters
            if has_invalid_chars(oligo):
                row.append(float('inf'))
                continue

            # Calculate Gibbs free energy
            try:
                gibbs = gibbs_rna_dna(oligo)
            except KeyError:
                row.append(float('inf'))
                continue

            # Check if within allowable range
            if gibbs < min_gibbs or gibbs > max_gibbs:
                row.append(float('inf'))
                continue

            # Badness is squared distance from target
            row.append((gibbs - target_gibbs) ** 2)

        badness_2d.append(row)

    return badness_2d


def _calc_score(old_score: float, k: int, new_badness: float) -> float:
    """Calculate running average score.

    This implements the scoring function from the MATLAB code:
    score = (old_score * k + new_badness) / (k + 1)

    Args:
        old_score: Previous average score for k probes
        k: Number of probes so far (0-indexed, so k=0 means first probe)
        new_badness: Badness of the new probe

    Returns:
        New average score
    """
    return (old_score * k + new_badness) / (k + 1)


def find_best_probes(
    badness: List[float],
    seq_len: int,
    oligo_len: int,
    spacer_len: int,
    n_probes: int,
) -> List[Tuple[float, List[int]]]:
    """Find optimal probe positions using dynamic programming.

    This implements the MATLAB find_best_matches algorithm.

    Args:
        badness: List of badness scores per position
        seq_len: Length of input sequence
        oligo_len: Length of oligonucleotides
        spacer_len: Minimum spacing between probes
        n_probes: Maximum number of probes to find

    Returns:
        List of (score, positions) tuples for 1..n_probes solutions.
        Each positions list contains 0-based start positions.
    """
    goodlen = len(badness)
    probe_spacer_len = oligo_len + spacer_len

    # Initialize DP tables
    # bmsf_pos[x][k] = position of probe k+1 when ending at x
    # bmsf_sco[x][k] = best score for k+1 probes ending at or before x
    bmsf_pos = [[None] * n_probes for _ in range(goodlen)]
    bmsf_sco = [[float('inf')] * n_probes for _ in range(goodlen)]

    # Initialize first position
    bmsf_pos[0][0] = 0
    bmsf_sco[0][0] = badness[0]

    # Fill DP tables
    for x in range(1, goodlen):
        # Copy previous best solutions
        for k in range(n_probes):
            bmsf_pos[x][k] = bmsf_pos[x - 1][k]
            bmsf_sco[x][k] = bmsf_sco[x - 1][k]

        # Try placing probe at position x
        for k in range(n_probes):
            potential_score = float('inf')

            if k == 0:
                # First probe: just use this position's badness
                potential_score = badness[x]
            else:
                # Later probes: need to check if we can place after previous
                prev_x = x - probe_spacer_len
                if prev_x >= 0 and bmsf_pos[prev_x][k - 1] is not None:
                    potential_score = _calc_score(
                        bmsf_sco[prev_x][k - 1], k, badness[x]
                    )

            if potential_score < bmsf_sco[x][k]:
                bmsf_pos[x][k] = x
                bmsf_sco[x][k] = potential_score

    # Backtrack to find probe positions
    results = []
    for k in range(n_probes):
        x = goodlen - 1
        curr_k = k

        if bmsf_pos[x][curr_k] is None:
            continue

        score = bmsf_sco[x][curr_k]

        # Only include solutions with reasonable scores
        if score >= 1_000_000:
            continue

        positions = []
        for _ in range(k + 1):
            pos = bmsf_pos[x][curr_k]
            positions.append(pos)
            x = pos - probe_spacer_len
            curr_k -= 1

        positions.reverse()
        results.append((score, positions))

    return results


def find_best_probes_mixed(
    badness_2d: List[List[float]],
    seq_len: int,
    min_len: int,
    max_len: int,
    spacer_len: int,
    n_probes: int,
) -> List[Tuple[float, List[Tuple[int, int]]]]:
    """Find optimal probe positions with mixed-length probes using DP.

    Uses an end-position-indexed DP to handle variable-length probes.
    The DP considers all valid lengths at every position and finds the
    globally optimal combination of probe placements and lengths.

    State: dp[e][k] = best average badness for k+1 probes with last probe
    ending at or before position e.

    Args:
        badness_2d: 2D badness array [pos][len_idx] from calculate_badness_mixed()
        seq_len: Length of input sequence
        min_len: Minimum probe length
        max_len: Maximum probe length
        spacer_len: Minimum gap between probes
        n_probes: Maximum number of probes to find

    Returns:
        List of (score, [(start, length), ...]) tuples for 1..n_probes solutions.
    """
    n_lengths = max_len - min_len + 1
    max_goodlen = len(badness_2d)  # number of valid start positions

    # DP indexed by end position (0..seq_len-1)
    # dp[e][k] = best average badness for k+1 probes ending at or before e
    dp = [[float('inf')] * n_probes for _ in range(seq_len)]
    # tracker[e][k] = (start_pos, length) of the last probe in best solution
    tracker = [[None] * n_probes for _ in range(seq_len)]

    for e in range(seq_len):
        # Step 1: Propagate from previous end position
        if e > 0:
            for k in range(n_probes):
                if dp[e - 1][k] < dp[e][k]:
                    dp[e][k] = dp[e - 1][k]
                    tracker[e][k] = tracker[e - 1][k]

        # Step 2: Try placing a probe ending at exactly position e
        for L_idx in range(n_lengths):
            L = min_len + L_idx
            x = e - L + 1  # start position
            if x < 0 or x >= max_goodlen:
                continue
            b = badness_2d[x][L_idx]
            if b == float('inf'):
                continue

            for k in range(n_probes):
                if k == 0:
                    score = b
                else:
                    prev_end = x - spacer_len - 1
                    if prev_end < 0 or tracker[prev_end][k - 1] is None:
                        continue
                    score = _calc_score(dp[prev_end][k - 1], k, b)

                if score < dp[e][k]:
                    dp[e][k] = score
                    tracker[e][k] = (x, L)

    # Backtrack to find probe positions for each number of probes
    results = []
    for k in range(n_probes):
        e = seq_len - 1

        if tracker[e][k] is None:
            continue

        score = dp[e][k]

        # Only include solutions with reasonable scores
        if score >= 1_000_000:
            continue

        positions_and_lengths = []
        curr_k = k
        for _ in range(k + 1):
            if tracker[e][curr_k] is None:
                break
            x, L = tracker[e][curr_k]
            positions_and_lengths.append((x, L))
            e = x - spacer_len - 1
            curr_k -= 1

        positions_and_lengths.reverse()
        results.append((score, positions_and_lengths))

    return results


def design_probes(
    input_file: str,
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
    repeatmask_file: Optional[str] = None,
    save_bowtie_raw: bool = False,
    mixed_lengths: Optional[Tuple[int, int]] = None,
    hp_threshold: int = 5,
    di_threshold: int = 3,
) -> ProbeDesignResult:
    """Design oligonucleotide probes for a target sequence.

    Main entry point for probe design. Reads a FASTA file and designs
    optimal probes based on thermodynamic properties.

    Args:
        input_file: Path to input FASTA file
        n_probes: Number of probes to design (default 48)
        oligo_length: Length of each oligonucleotide (default 20, used when
            mixed_lengths is None)
        spacer_length: Minimum gap between probes (default 3)
        target_gibbs: Target Gibbs free energy in kcal/mol (default -23)
        allowable_gibbs: (min, max) Gibbs FE range (default (-26, -20))
        output_name: Base name for output files (default: derived from input)
        species: Species for masking databases (default: "human")
        pseudogene_mask: Whether to mask pseudogene alignments (default: False)
        genome_mask: Whether to mask repetitive genomic regions (default: False)
        index_dir: Directory containing bowtie indexes (default: auto-detect)
        repeatmask_file: Path to repeatmask FASTA file (default: None)
        save_bowtie_raw: Whether to save raw bowtie output (default: False)
        mixed_lengths: Optional (min_len, max_len) tuple for mixed-length
            probe design. When provided, probes of varying lengths within
            this range are considered. When None, uses fixed oligo_length.

    Returns:
        ProbeDesignResult containing the designed probes
    """
    use_mixed = mixed_lengths is not None
    if use_mixed:
        min_len, max_len = mixed_lengths
        if min_len > max_len:
            raise ValueError(f"mixed_lengths min ({min_len}) is greater than max ({max_len})")
    else:
        min_len = max_len = oligo_length

    # Read input sequence
    headers, seqs = read_fasta(input_file)

    # Concatenate multi-entry FASTA into single sequence (with '>' marking junctions)
    full_seq = sequences_to_single_string(seqs, mark_junctions=True)

    # Check for N's in the sequence (indicates pre-repeatmasked input)
    has_n_masking = 'n' in full_seq.lower()

    # Clean sequence (lowercase, keep only actgn>)
    seq = clean_sequence(full_seq)

    # Determine output name
    if output_name is None:
        output_name = get_basename(input_file)

    run_info = new_run_info(
        "smFISH", input_file, len(seqs), seq,
        params={
            "n_probes": n_probes,
            "oligo_length": oligo_length,
            "mixed_lengths": tuple(mixed_lengths) if use_mixed else None,
            "spacer_length": spacer_length,
            "target_gibbs": target_gibbs,
            "allowable_gibbs": tuple(allowable_gibbs),
            "species": species,
            "pseudogene_mask": pseudogene_mask,
            "genome_mask": genome_mask,
            "hp_threshold": hp_threshold,
            "di_threshold": di_threshold,
            "save_bowtie_raw": save_bowtie_raw,
        },
        repeatmask_file=repeatmask_file,
        index_dir=index_dir,
    )

    # Calculate badness scores (thermodynamic filtering)
    if use_mixed:
        badness_2d = calculate_badness_mixed(
            seq, min_len, max_len, target_gibbs, allowable_gibbs
        )
        # For F mask: position is valid if ANY length gives finite badness
        max_goodlen = len(badness_2d)
        n_lengths = max_len - min_len + 1
    else:
        badness = calculate_badness(
            seq, oligo_length, target_gibbs, allowable_gibbs
        )
        goodlen = len(badness)

    # Apply masking if requested
    full_mask = [0] * len(seq)
    mask_strings = []
    mask_kinds = []

    # Handle repeat masking from separate file or from N's in input
    if repeatmask_file:
        # Read repeatmasked sequence from separate file
        _, rm_seqs = read_fasta(repeatmask_file)
        rm_full_seq = sequences_to_single_string(rm_seqs, mark_junctions=True)
        rm_seq = clean_sequence(rm_full_seq)

        # Create repeat mask from N positions in the repeatmask file
        rmask = [1 if rm_seq[i].lower() == 'n' else 0 for i in range(len(rm_seq))]
        rstr = "".join("R" if rmask[i] else seq[i] for i in range(len(seq)))
        mask_strings.append(rstr)
        mask_kinds.append("R")
        for i, v in enumerate(rmask):
            full_mask[i] += v
        print(f"Repeat masking (from file): {sum(rmask)} positions masked")
    elif has_n_masking:
        # Create repeat mask from N positions in input file
        rmask = [1 if seq[i].lower() == 'n' else 0 for i in range(len(seq))]
        rstr = "".join("R" if rmask[i] else seq[i] for i in range(len(seq)))
        mask_strings.append(rstr)
        mask_kinds.append("R")
        for i, v in enumerate(rmask):
            full_mask[i] += v
        print(f"Repeat masking (from N's): {sum(rmask)} positions masked")

    pseudo_raw_text = None
    genome_raw_text = None
    pseudo_hits = None
    genome_hits_dict = None

    if pseudogene_mask or genome_mask:
        from .masking import (
            pseudogene_mask as get_pseudogene_mask,
            genome_mask as get_genome_mask,
        )

        idx_dir = Path(index_dir) if index_dir else None

        if pseudogene_mask:
            try:
                pmask, pseudo_hits, pseudo_raw = get_pseudogene_mask(seq, species, idx_dir)
                if not save_bowtie_raw:
                    pseudo_raw = None
                pseudo_raw_text = pseudo_raw
                for i, v in enumerate(pmask):
                    full_mask[i] += v
                # Create visualization string
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
                    full_mask[i] += v
                # Create visualization string
                gstr = "".join("B" if gmask[i] else seq[i] for i in range(len(seq)))
                mask_strings.append(gstr)
                mask_kinds.append("B")
                run_info["masks"].append(bowtie_mask_record("B", species, idx_dir))
                print(f"Genome masking: {sum(gmask)} positions masked")
            except Exception as e:
                run_info["masks"].append(bowtie_mask_record("B", species, idx_dir, error=e))
                print(f"Warning: Genome masking failed: {e}")

    # Create F mask and apply masking to badness
    if use_mixed:
        # F mask for mixed-length: show 'F' only if NO length gives finite badness
        fstr_parts = []
        for i in range(len(seq)):
            valid_for_any = False
            if i < max_goodlen:
                for L_idx in range(n_lengths):
                    if badness_2d[i][L_idx] != float('inf'):
                        valid_for_any = True
                        break
            if valid_for_any:
                fstr_parts.append(seq[i])
            else:
                fstr_parts.append('F')
        fstr = "".join(fstr_parts)
        mask_strings.append(fstr)
        mask_kinds.append("F")

        # --- L-mask: low-complexity filtering (always active, per-length) ---
        from .masking import has_homopolymer, has_dinucleotide_repeat
        l_masked_count = 0
        lstr_parts = []
        for i in range(len(seq)):
            if i < max_goodlen:
                any_length_masked = False
                for L_idx in range(n_lengths):
                    L = min_len + L_idx
                    if i + L > len(seq):
                        continue
                    if badness_2d[i][L_idx] == float('inf'):
                        continue  # already excluded by F-mask
                    probe_sub = seq[i:i + L]
                    if (has_homopolymer(probe_sub, hp_threshold)
                            or has_dinucleotide_repeat(probe_sub, di_threshold)):
                        badness_2d[i][L_idx] = float('inf')
                        l_masked_count += 1
                        any_length_masked = True
                # Visualization: 'L' only if ALL lengths are now dead
                # and at least one was killed by L-mask (not just F-mask)
                any_viable = any(
                    badness_2d[i][L_idx] != float('inf')
                    for L_idx in range(n_lengths)
                    if i + min_len + L_idx <= len(seq)
                )
                lstr_parts.append('L' if not any_viable and any_length_masked else seq[i])
            else:
                lstr_parts.append(seq[i])
        lstr = "".join(lstr_parts)
        mask_strings.append(lstr)
        mask_kinds.append("L")
        if l_masked_count > 0:
            print(f"Low-complexity masking: {l_masked_count} (position,length) pairs masked "
                  f"(homopolymer>={hp_threshold}, dinucleotide>={di_threshold})")

        # Apply mask to 2D badness (after F mask is created)
        if any(full_mask):
            from .masking import mask_to_badness_mixed
            mask_bad_2d = mask_to_badness_mixed(full_mask, min_len, max_len)
            for i in range(len(badness_2d)):
                for L_idx in range(n_lengths):
                    if i < len(mask_bad_2d) and mask_bad_2d[i][L_idx] == float('inf'):
                        badness_2d[i][L_idx] = float('inf')

        # Find optimal probe positions using mixed-length DP
        solutions = find_best_probes_mixed(
            badness_2d, len(seq), min_len, max_len, spacer_length, n_probes
        )
    else:
        # F mask for fixed-length (original behavior)
        fstr_parts = []
        for i in range(len(seq)):
            if i < goodlen and badness[i] != float('inf'):
                fstr_parts.append(seq[i])
            else:
                fstr_parts.append('F')
        fstr = "".join(fstr_parts)
        mask_strings.append(fstr)
        mask_kinds.append("F")

        # --- L-mask: low-complexity filtering (always active, probe-level) ---
        from .masking import has_homopolymer, has_dinucleotide_repeat
        l_masked_count = 0
        lstr_parts = []
        for i in range(len(seq)):
            if i < goodlen:
                probe_sub = seq[i:i + oligo_length]
                if badness[i] != float('inf') and (
                        has_homopolymer(probe_sub, hp_threshold)
                        or has_dinucleotide_repeat(probe_sub, di_threshold)):
                    badness[i] = float('inf')
                    l_masked_count += 1
                    lstr_parts.append('L')
                else:
                    lstr_parts.append(seq[i])
            else:
                lstr_parts.append(seq[i])
        lstr = "".join(lstr_parts)
        mask_strings.append(lstr)
        mask_kinds.append("L")
        if l_masked_count > 0:
            print(f"Low-complexity masking: {l_masked_count} probe positions masked "
                  f"(homopolymer>={hp_threshold}, dinucleotide>={di_threshold})")

        # Apply mask to badness (after F mask is created)
        if any(full_mask):
            from .masking import mask_to_badness
            mask_badness = mask_to_badness(full_mask, oligo_length)
            for i in range(len(badness)):
                if mask_badness[i] == float('inf'):
                    badness[i] = float('inf')

        # Find optimal probe positions using fixed-length DP
        solutions = find_best_probes(
            badness, len(seq), oligo_length, spacer_length, n_probes
        )

    if not solutions:
        return ProbeDesignResult(
            probes=[],
            score=float('inf'),
            input_sequence=seq,
            template_name=output_name,
            mask=full_mask,
            mask_strings=mask_strings,
            mask_kinds=mask_kinds,
            run_info=run_info,
            bowtie_pseudogene_hits=pseudo_hits,
            bowtie_genome_hits=genome_hits_dict,
            bowtie_pseudogene_raw=pseudo_raw_text,
            bowtie_genome_raw=genome_raw_text,
        )

    # Use the solution with the most probes (last in list)
    best_solution = solutions[-1]

    # Generate probe objects
    probes = []
    if use_mixed:
        best_score, positions_and_lengths = best_solution
        for i, (pos, probe_len) in enumerate(positions_and_lengths):
            # Extract template region, skipping '>' characters
            template_region = ""
            j = pos
            while len(template_region) < probe_len and j < len(seq):
                if seq[j] != '>':
                    template_region += seq[j]
                j += 1

            probe_seq = reverse_complement(template_region)

            probe = Probe(
                index=i + 1,
                position=pos,
                length=len(probe_seq),
                sequence=probe_seq,
                gc_percent=round(percent_gc(probe_seq) * 100),
                tm=round(tm_rna_dna(template_region), 1),
                gibbs_fe=round(gibbs_rna_dna(template_region), 1),
                name=f"{output_name}_{i + 1}",
            )
            probes.append(probe)
    else:
        best_score, positions = best_solution
        for i, pos in enumerate(positions):
            # Extract template region, skipping '>' characters
            template_region = ""
            j = pos
            while len(template_region) < oligo_length and j < len(seq):
                if seq[j] != '>':
                    template_region += seq[j]
                j += 1

            probe_seq = reverse_complement(template_region)

            probe = Probe(
                index=i + 1,
                position=pos,
                length=len(probe_seq),
                sequence=probe_seq,
                gc_percent=round(percent_gc(probe_seq) * 100),
                tm=round(tm_rna_dna(template_region), 1),
                gibbs_fe=round(gibbs_rna_dna(template_region), 1),
                name=f"{output_name}_{i + 1}",
            )
            probes.append(probe)

    return ProbeDesignResult(
        probes=probes,
        score=best_score,
        input_sequence=seq,
        template_name=output_name,
        mask=full_mask,
        mask_strings=mask_strings,
        mask_kinds=mask_kinds,
        run_info=run_info,
        bowtie_pseudogene_hits=pseudo_hits,
        bowtie_genome_hits=genome_hits_dict,
        bowtie_pseudogene_raw=pseudo_raw_text,
        bowtie_genome_raw=genome_raw_text,
    )
