"""HCR-specific output formatting for probe pair designs.

Generates _HCR_oligos.txt, _HCR_seq.txt, and _HCR_hits.txt files.
"""

import os
from typing import Dict, List, Optional
from .hcr import HCRDesignResult, HCRProbePair, HALF_LENGTH, GAP_LENGTH, PAIR_BLOCK_LENGTH
from .provenance import format_header
from .sequence import complement


def _nucleotide_position(seq: str, pos: int) -> int:
    """Convert 0-based string position to 1-based nucleotide position.

    Adjusts for '>' junction markers in the sequence string.
    """
    return pos - seq[:pos].count('>') + 1


def format_hcr_oligos(result: HCRDesignResult) -> str:
    """Format HCR probe pairs as TSV content for _HCR_oligos.txt.

    Columns: pair, half, start, GC%, Tm, Gibbs, strict, full_oligo, name

    Args:
        result: HCRDesignResult from design_hcr_probes()

    Returns:
        TSV-formatted string
    """
    lines = ["pair\thalf\tstart\tGC%\tTm\tGibbs\tstrict\tfull_oligo\tname"]
    for pair in result.pairs:
        # Odd number = P1, even = P2
        p1_id = f"{2 * pair.pair_index - 1:02d}"
        p2_id = f"{2 * pair.pair_index:02d}"
        p1_name = f"{result.template_name}_HCR{result.amplifier}_{p1_id}"
        p2_name = f"{result.template_name}_HCR{result.amplifier}_{p2_id}"

        left_start = _nucleotide_position(result.input_sequence, pair.left.position)
        right_start = _nucleotide_position(result.input_sequence, pair.right.position)

        # P1 (odd) = right half of target (carries initiator_b)
        lines.append(
            f"{pair.pair_index}\tP1\t{right_start}\t"
            f"{pair.right.gc_percent}\t{pair.right.tm}\t{pair.right.gibbs_fe}\t"
            f"{'yes' if pair.right.is_strict else 'no'}\t"
            f"{pair.right.oligo_seq}\t{p1_name}"
        )
        # P2 (even) = left half of target (carries initiator_a)
        lines.append(
            f"{pair.pair_index}\tP2\t{left_start}\t"
            f"{pair.left.gc_percent}\t{pair.left.tm}\t{pair.left.gibbs_fe}\t"
            f"{'yes' if pair.left.is_strict else 'no'}\t"
            f"{pair.left.oligo_seq}\t{p2_name}"
        )

    return "\n".join(lines) + "\n"


def format_hcr_seq(
    result: HCRDesignResult,
    line_width: int = 110,
) -> str:
    """Format sequence alignment visualization for _HCR_seq.txt.

    Starts with the '#' run header (provenance.format_header), then shows the
    target sequence with mask lines and pair annotations.

    Args:
        result: HCRDesignResult from design_hcr_probes()
        line_width: Characters per line for wrapping

    Returns:
        Formatted visualization string
    """
    seq = result.input_sequence
    n = len(seq)

    # Build pair annotation strings
    pair_annot = [' '] * n
    pair_labels = [' '] * n

    for pair in result.pairs:
        p = pair.pair_position
        # Left half: [===LEFT===]
        for i in range(HALF_LENGTH):
            if p + i < n:
                pair_annot[p + i] = '='
        # Gap: --
        for i in range(GAP_LENGTH):
            if p + HALF_LENGTH + i < n:
                pair_annot[p + HALF_LENGTH + i] = '-'
        # Right half: [===RIGHT===]
        right_start = p + HALF_LENGTH + GAP_LENGTH
        for i in range(HALF_LENGTH):
            if right_start + i < n:
                pair_annot[right_start + i] = '='

        # Label
        left_pos = _nucleotide_position(seq, pair.left.position)
        label = f"Pair {pair.pair_index} (pos {left_pos})"
        for i, c in enumerate(label):
            if p + i < n:
                pair_labels[p + i] = c

    annot_str = ''.join(pair_annot)
    label_str = ''.join(pair_labels)

    # Mask letters come from design_hcr_probes(); guess them for results built by hand
    mask_chars = result.mask_kinds or [
        next((c for c in m if c in "RPBFL"), None) for m in result.mask_strings
    ]

    # Build output
    output = []
    for start in range(0, n, line_width):
        end = min(start + line_width, n)

        # Position header
        start_pos = _nucleotide_position(seq, start)
        end_pos = _nucleotide_position(seq, end - 1)
        output.append(f"Position: {start_pos}-{end_pos}")

        # All lines use a fixed-width prefix for alignment
        # Format: "LABEL:  DATA" where LABEL is right-justified to 8 chars
        output.append(f"{'SEQUENCE':>8s}:  {seq[start:end]}")

        # Mask lines
        for mask_char, mask_str in zip(mask_chars, result.mask_strings):
            if mask_str:
                label = f"{mask_char} mask" if mask_char else "mask"
                output.append(f"{label:>8s}:  {mask_str[start:end]}")

        # Pair annotations
        output.append(f"{'Pairs':>8s}:  {annot_str[start:end]}")
        output.append(f"{'':>8s}   {label_str[start:end]}")
        output.append("")

    return format_header(result) + "\n".join(output) + "\n"


def format_hcr_hits(result: HCRDesignResult) -> str:
    """Format bowtie hit summary for _HCR_hits.txt.

    Shows per-pair, per-half bowtie hit counts with strict/lenient annotations.

    Args:
        result: HCRDesignResult from design_hcr_probes()

    Returns:
        Formatted hit summary string
    """
    if not result.pairs:
        return "No pairs designed — no hit data.\n"

    lines = []
    for pair in result.pairs:
        for half_label, half in [("P1 (right", pair.right), ("P2 (left", pair.left)]:
            strict_label = "STRICT" if half.is_strict else "LENIENT"
            start = _nucleotide_position(result.input_sequence, half.position)
            end = start + HALF_LENGTH - 1
            lines.append(f"=== Pair {pair.pair_index}, {half_label}, pos {start}-{end}, {strict_label}) ===")

            if result.bowtie_pseudogene_hits is not None:
                # Count hits in this half's region
                hits = sum(result.bowtie_pseudogene_hits[half.position:half.position + HALF_LENGTH])
                lines.append(f"Pseudogene hits: {hits}")
            else:
                lines.append("Pseudogene hits: not checked")

            if result.bowtie_genome_hits is not None:
                # Sum hits across all tiers
                total = 0
                for tier, arr in result.bowtie_genome_hits.items():
                    total += sum(arr[half.position:half.position + HALF_LENGTH])
                lines.append(f"Genome hits: {total}")
            else:
                lines.append("Genome hits: not checked")

            if not half.is_strict:
                lines.append("  (LENIENT — hits recorded but not used for filtering)")

            lines.append("")

    return "\n".join(lines) + "\n"


def write_hcr_output_files(
    result: HCRDesignResult,
    output_prefix: str,
    output_dir: Optional[str] = None,
) -> None:
    """Write all HCR output files (_HCR_oligos.txt, _HCR_seq.txt, _HCR_hits.txt).

    Args:
        result: HCRDesignResult from design_hcr_probes()
        output_prefix: Prefix for output files
        output_dir: Directory for output files (default: CWD)
    """
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        prefix = os.path.join(output_dir, output_prefix)
    else:
        prefix = output_prefix

    oligos_content = format_hcr_oligos(result)
    with open(f"{prefix}_HCR_oligos.txt", 'w') as f:
        f.write(oligos_content)

    seq_content = format_hcr_seq(result)
    with open(f"{prefix}_HCR_seq.txt", 'w') as f:
        f.write(seq_content)

    hits_content = format_hcr_hits(result)
    with open(f"{prefix}_HCR_hits.txt", 'w') as f:
        f.write(hits_content)

    # Raw bowtie output
    if result.bowtie_pseudogene_raw is not None:
        with open(f"{prefix}_HCR_bowtie_pseudogene_raw.txt", 'w') as f:
            f.write(result.bowtie_pseudogene_raw)
    if result.bowtie_genome_raw is not None:
        with open(f"{prefix}_HCR_bowtie_genome_raw.txt", 'w') as f:
            f.write(result.bowtie_genome_raw)
