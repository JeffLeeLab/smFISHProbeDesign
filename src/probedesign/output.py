"""Output file generation for probe designs."""

import os
from typing import Dict, List, Optional
from .core import Probe, ProbeDesignResult
from .provenance import format_header
from .sequence import complement


def _nucleotide_position(seq: str, pos: int) -> int:
    """Convert 0-based string position to 1-based nucleotide position.

    Adjusts for '>' junction markers in the sequence string.

    Args:
        seq: Input sequence string (may contain '>' junction markers)
        pos: 0-based position in the string

    Returns:
        1-based nucleotide position (excluding '>' markers)
    """
    return pos - seq[:pos].count('>') + 1


def write_oligos_file(result: ProbeDesignResult, filepath: str) -> None:
    """Write probe information to a TSV file.

    Output format (tab-separated):
    index  start  GC%  Tm  GibbsFE  sequence  name

    Args:
        result: ProbeDesignResult from design_probes()
        filepath: Output file path
    """
    with open(filepath, 'w') as f:
        for probe in result.probes:
            start_pos = _nucleotide_position(result.input_sequence, probe.position)
            f.write(
                f"{probe.index}\t"
                f"{start_pos}\t"
                f"{probe.gc_percent}\t"
                f"{probe.tm}\t"
                f"{probe.gibbs_fe}\t"
                f"{probe.sequence}\t"
                f"{probe.name}\n"
            )


def write_seq_file(
    result: ProbeDesignResult,
    filepath: str,
    mask_seqs: Optional[List[str]] = None,
    line_width: int = 110,
) -> None:
    """Write sequence alignment visualization file.

    Creates a multi-line visualization showing:
    - A '#' header with the run's parameters and masks (provenance.format_header)
    - Original sequence (with '>' at exon junctions where present)
    - Masked regions (if any) - each mask shows sequence with mask chars (R, P, B, F, L)
    - Probe alignments with complementary sequences
    - Probe labels

    Format matches MATLAB output:
    - No '>' prefix added to lines - the '>' only appears where it exists in the sequence
    - Mask lines show the original sequence with masked positions replaced by mask characters

    Args:
        result: ProbeDesignResult from design_probes()
        filepath: Output file path
        mask_seqs: List of mask strings (same length as sequence, with mask chars at masked positions)
        line_width: Characters per line for wrapping
    """
    with open(filepath, 'w') as f:
        f.write(format_seq_content(result, mask_seqs, line_width))


def format_probes_table(result: ProbeDesignResult) -> str:
    """Format probes as a printable table.

    Args:
        result: ProbeDesignResult from design_probes()

    Returns:
        Formatted string table
    """
    if not result.probes:
        return "No probes found."

    lines = [
        "Index\tStart\tGC%\tTm\tGibbs\tSequence\tName",
        "-" * 80,
    ]

    for probe in result.probes:
        start_pos = _nucleotide_position(result.input_sequence, probe.position)
        lines.append(
            f"{probe.index}\t"
            f"{start_pos}\t"
            f"{probe.gc_percent}\t"
            f"{probe.tm}\t"
            f"{probe.gibbs_fe}\t"
            f"{probe.sequence}\t"
            f"{probe.name}"
        )

    return "\n".join(lines)


def _format_hits_table(
    seq: str,
    hits_data: Dict[str, List[int]],
) -> str:
    """Format bowtie hit counts as a TSV table.

    Args:
        seq: Input sequence (may contain '>' junction markers)
        hits_data: Dict mapping column names to hit count arrays,
                   e.g. {"16mer": [...]} or {"12mer":[], "14mer":[], "16mer":[]}

    Returns:
        TSV-formatted string with header and one row per nucleotide position
    """
    col_names = list(hits_data.keys())
    header = "position\t" + "\t".join(f"hits_{name}" for name in col_names)
    lines = [header]

    # Find the maximum valid range across all hit arrays
    max_len = max(len(arr) for arr in hits_data.values())

    # Iterate through array positions, skipping '>' junction markers
    nuc_pos = 0
    for i in range(max_len):
        if i < len(seq) and seq[i] == '>':
            continue  # Skip junction markers
        nuc_pos += 1  # 1-based nucleotide position

        values = []
        for name in col_names:
            arr = hits_data[name]
            values.append(str(arr[i]) if i < len(arr) else "0")

        lines.append(f"{nuc_pos}\t" + "\t".join(values))

    return "\n".join(lines) + "\n"


def write_output_files(
    result: ProbeDesignResult,
    output_prefix: str,
    mask_seqs: Optional[List[str]] = None,
    output_dir: Optional[str] = None,
) -> None:
    """Write both oligos and seq output files, plus bowtie QC files if available.

    Args:
        result: ProbeDesignResult from design_probes()
        output_prefix: Prefix for output files (will add _oligos.txt and _seq.txt)
        mask_seqs: Optional mask sequences for seq file
        output_dir: Optional directory for output files. When provided, files are
                    written to this directory. When None, files are written to CWD.
    """
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        prefix = os.path.join(output_dir, output_prefix)
    else:
        prefix = output_prefix

    write_oligos_file(result, f"{prefix}_oligos.txt")
    write_seq_file(result, f"{prefix}_seq.txt", mask_seqs)

    # Write parsed bowtie hits tables (compact default output)
    if result.bowtie_pseudogene_hits is not None:
        table = _format_hits_table(
            result.input_sequence,
            {"16mer": result.bowtie_pseudogene_hits},
        )
        with open(f"{prefix}_bowtie_pseudogene.txt", 'w') as f:
            f.write(table)

    if result.bowtie_genome_hits is not None:
        table = _format_hits_table(
            result.input_sequence,
            result.bowtie_genome_hits,
        )
        with open(f"{prefix}_bowtie_genome.txt", 'w') as f:
            f.write(table)

    # Write raw bowtie alignment files (only when --save-bowtie-raw is used)
    if result.bowtie_pseudogene_raw is not None:
        with open(f"{prefix}_bowtie_pseudogene_raw.txt", 'w') as f:
            f.write(result.bowtie_pseudogene_raw)
    if result.bowtie_genome_raw is not None:
        with open(f"{prefix}_bowtie_genome_raw.txt", 'w') as f:
            f.write(result.bowtie_genome_raw)


def format_oligos_content(result: ProbeDesignResult) -> str:
    """Return oligos file content as a string (for download buttons).

    Produces identical content to write_oligos_file() without file I/O.

    Args:
        result: ProbeDesignResult from design_probes()

    Returns:
        TSV-formatted string with probe information
    """
    lines = []
    for probe in result.probes:
        start_pos = _nucleotide_position(result.input_sequence, probe.position)
        lines.append(
            f"{probe.index}\t"
            f"{start_pos}\t"
            f"{probe.gc_percent}\t"
            f"{probe.tm}\t"
            f"{probe.gibbs_fe}\t"
            f"{probe.sequence}\t"
            f"{probe.name}"
        )
    return "\n".join(lines) + "\n" if lines else ""


def format_seq_content(
    result: ProbeDesignResult,
    mask_seqs: Optional[List[str]] = None,
    line_width: int = 110,
) -> str:
    """Return sequence alignment visualization as a string (for download buttons).

    This is the content write_seq_file() writes: the '#' run header, then the
    sequence overlay.

    Args:
        result: ProbeDesignResult from design_probes()
        mask_seqs: List of mask strings (same length as sequence)
        line_width: Characters per line for wrapping

    Returns:
        Formatted visualization string
    """
    seq = result.input_sequence

    # Build probe alignment string
    probe_align = [' '] * len(seq)
    probe_labels = [' '] * len(seq)

    for probe in result.probes:
        pos = probe.position
        oligo_len = len(probe.sequence)  # per probe: mixed-length designs vary
        actual_seq = ''
        seq_pos = pos
        chars_collected = 0
        while chars_collected < oligo_len and seq_pos < len(seq):
            if seq[seq_pos] != '>':
                actual_seq += seq[seq_pos]
                chars_collected += 1
            seq_pos += 1

        comp_seq = complement(actual_seq)

        comp_idx = 0
        for i in range(oligo_len + 10):
            if pos + i >= len(probe_align):
                break
            if seq[pos + i] == '>':
                continue
            if comp_idx < len(comp_seq):
                probe_align[pos + i] = comp_seq[comp_idx]
                comp_idx += 1

        start_pos = _nucleotide_position(seq, probe.position)
        label = f"Prb# {probe.index},Pos {start_pos},FE {probe.gibbs_fe},GC {probe.gc_percent}"
        for i, c in enumerate(label):
            if pos + i < len(probe_labels):
                probe_labels[pos + i] = c

    probe_align_str = ''.join(probe_align)
    probe_labels_str = ''.join(probe_labels)

    # Build output lines
    output_lines = []
    for start in range(0, len(seq), line_width):
        end = min(start + line_width, len(seq))

        output_lines.append(seq[start:end])

        if mask_seqs:
            for mask in mask_seqs:
                if mask:
                    output_lines.append(mask[start:end])

        output_lines.append(probe_align_str[start:end])
        output_lines.append(probe_labels_str[start:end])
        output_lines.append("")

    body = "\n".join(output_lines) + "\n" if output_lines else ""
    drawn = [k for k, m in zip(result.mask_kinds, result.mask_strings) if m] if mask_seqs else []
    return format_header(result, drawn_kinds=drawn) + body


def format_hits_content(result: ProbeDesignResult) -> Dict[str, str]:
    """Return bowtie hit file contents as a dict of {suffix: content_string}.

    Args:
        result: ProbeDesignResult from design_probes()

    Returns:
        Dict mapping filename suffix to TSV content string.
        Possible keys: "bowtie_pseudogene.txt", "bowtie_genome.txt"
    """
    files = {}
    if result.bowtie_pseudogene_hits is not None:
        files["bowtie_pseudogene.txt"] = _format_hits_table(
            result.input_sequence,
            {"16mer": result.bowtie_pseudogene_hits},
        )
    if result.bowtie_genome_hits is not None:
        files["bowtie_genome.txt"] = _format_hits_table(
            result.input_sequence,
            result.bowtie_genome_hits,
        )
    return files
