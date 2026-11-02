"""Command-line interface for ProbeDesign."""

import click
from pathlib import Path

from .core import design_probes
from .output import write_output_files, format_probes_table


@click.group()
@click.version_option()
def main():
    """ProbeDesign - Design oligonucleotide probes for FISH experiments."""
    pass


@main.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option(
    '-n', '--probes', 'n_probes',
    default=48,
    type=int,
    help='Number of probes to design (default: 48)'
)
@click.option(
    '-l', '--oligo-length',
    default='20',
    type=str,
    help='Oligo length: single int (e.g. 20) or range for mixed lengths (e.g. 18-22)'
)
@click.option(
    '-s', '--spacer-length',
    default=2,
    type=int,
    help='Minimum gap between probes (default: 2)'
)
@click.option(
    '-g', '--target-gibbs',
    default=-23.0,
    type=float,
    help='Target Gibbs free energy in kcal/mol (default: -23)'
)
@click.option(
    '--allowable-gibbs',
    default='-26,-20',
    type=str,
    help='Allowable Gibbs FE range as min,max (default: -26,-20)'
)
@click.option(
    '-o', '--output',
    default=None,
    type=str,
    help='Output file prefix (default: derived from input filename)'
)
@click.option(
    '--quiet', '-q',
    is_flag=True,
    help='Suppress output to stdout'
)
@click.option(
    '--species',
    default='human',
    type=click.Choice(['human', 'mouse', 'elegans', 'drosophila', 'rat']),
    help='Species for masking databases (default: human)'
)
@click.option(
    '--pseudogene-mask/--no-pseudogene-mask',
    default=False,
    help='Mask regions that align to pseudogenes (default: off)'
)
@click.option(
    '--genome-mask/--no-genome-mask',
    default=False,
    help='Mask repetitive genomic regions (default: off)'
)
@click.option(
    '--index-dir',
    default=None,
    type=click.Path(exists=True, file_okay=False),
    help='Directory containing bowtie indexes (default: auto-detect)'
)
@click.option(
    '--repeatmask-file',
    default=None,
    type=click.Path(exists=True),
    help='FASTA file with N\'s marking repeat regions (for manual repeat masking)'
)
@click.option(
    '--repeatmask/--no-repeatmask',
    default=False,
    help='Run RepeatMasker to automatically mask repeat regions (default: off)'
)
@click.option(
    '--save-bowtie-raw',
    is_flag=True,
    default=False,
    help='Save raw bowtie alignment output (large files)'
)
@click.option(
    '--hp-threshold',
    default=5,
    type=click.IntRange(3, 10),
    help='Homopolymer repeat threshold for L-mask (default: 5, range: 3-10). '
         'Probes containing single-nucleotide runs >= this length are excluded.'
)
@click.option(
    '--di-threshold',
    default=3,
    type=click.IntRange(3, 10),
    help='Dinucleotide repeat threshold for L-mask (default: 3, range: 3-10). '
         'Probes containing dinucleotide motifs repeated >= this many times are excluded.'
)
def design(
    input_file: str,
    n_probes: int,
    oligo_length: str,
    spacer_length: int,
    target_gibbs: float,
    allowable_gibbs: str,
    output: str,
    quiet: bool,
    species: str,
    pseudogene_mask: bool,
    genome_mask: bool,
    index_dir: str,
    repeatmask_file: str,
    repeatmask: bool,
    save_bowtie_raw: bool,
    hp_threshold: int,
    di_threshold: int,
):
    """Design probes for a target sequence.

    INPUT_FILE is a FASTA file containing the target sequence.

    Examples:

        probedesign design input.fa

        probedesign design input.fa -n 6 -l 20 -o MyProbes

        probedesign design input.fa --target-gibbs -23 --allowable-gibbs -26,-20

        probedesign design input.fa --pseudogene-mask --species human
    """
    # Parse allowable Gibbs range
    try:
        min_gibbs, max_gibbs = map(float, allowable_gibbs.split(','))
    except ValueError:
        raise click.BadParameter(
            f"Invalid format for --allowable-gibbs: {allowable_gibbs}. "
            "Expected format: min,max (e.g., -26,-20)"
        )

    # Parse oligo length: single int or range (e.g. "20" or "18-22")
    mixed_lengths = None
    parsed_oligo_length = 20  # default
    oligo_length_str = oligo_length.strip()
    if '-' in oligo_length_str:
        parts = oligo_length_str.split('-')
        if len(parts) == 2:
            try:
                lo, hi = int(parts[0]), int(parts[1])
                if lo > hi:
                    raise click.BadParameter(
                        f"Invalid range for --oligo-length: {oligo_length_str}. "
                        "Min must be <= max (e.g., 18-22)"
                    )
                if lo < 10 or hi > 60:
                    raise click.BadParameter(
                        f"Invalid range for --oligo-length: {oligo_length_str}. "
                        "Values must be between 10 and 60"
                    )
                mixed_lengths = (lo, hi)
                parsed_oligo_length = lo  # fallback; not used in mixed mode
            except ValueError:
                raise click.BadParameter(
                    f"Invalid format for --oligo-length: {oligo_length_str}. "
                    "Expected integer or range (e.g., 20 or 18-22)"
                )
        else:
            raise click.BadParameter(
                f"Invalid format for --oligo-length: {oligo_length_str}. "
                "Expected integer or range (e.g., 20 or 18-22)"
            )
    else:
        try:
            parsed_oligo_length = int(oligo_length_str)
        except ValueError:
            raise click.BadParameter(
                f"Invalid format for --oligo-length: {oligo_length_str}. "
                "Expected integer or range (e.g., 20 or 18-22)"
            )

    # Determine output prefix
    if output is None:
        output = Path(input_file).stem

    # Handle --repeatmask flag: run RepeatMasker and use output as repeatmask_file
    actual_repeatmask_file = repeatmask_file
    if repeatmask:
        if repeatmask_file:
            raise click.BadParameter(
                "Cannot use both --repeatmask and --repeatmask-file. "
                "Choose one or the other."
            )
        if not quiet:
            click.echo(f"Running RepeatMasker on {input_file}...")
        try:
            from .masking import run_repeatmasker
            masked_file = run_repeatmasker(Path(input_file), species=species)
            actual_repeatmask_file = str(masked_file)
            if not quiet:
                click.echo(f"RepeatMasker output: {masked_file}")
        except FileNotFoundError as e:
            raise click.ClickException(str(e))
        except RuntimeError as e:
            raise click.ClickException(str(e))

    # Run probe design
    if not quiet:
        if mixed_lengths:
            click.echo(f"Designing mixed-length probes ({mixed_lengths[0]}-{mixed_lengths[1]}bp) for {input_file}...")
        else:
            click.echo(f"Designing probes for {input_file}...")

    result = design_probes(
        input_file=input_file,
        n_probes=n_probes,
        oligo_length=parsed_oligo_length,
        spacer_length=spacer_length,
        target_gibbs=target_gibbs,
        allowable_gibbs=(min_gibbs, max_gibbs),
        output_name=Path(output).name,  # -o may include a directory; names don't
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

    if not result.probes:
        click.echo("ERROR: Could not find any valid probes.", err=True)
        raise SystemExit(1)

    # Write output files
    write_output_files(result, output, mask_seqs=result.mask_strings)

    # Print summary
    if not quiet:
        click.echo(f"\nFound {len(result.probes)} probes (score: {result.score:.4f})")
        click.echo(f"\nOutput files:")
        click.echo(f"  {output}_oligos.txt")
        click.echo(f"  {output}_seq.txt")
        if result.bowtie_pseudogene_hits is not None:
            click.echo(f"  {output}_bowtie_pseudogene.txt")
        if result.bowtie_genome_hits is not None:
            click.echo(f"  {output}_bowtie_genome.txt")
        if result.bowtie_pseudogene_raw is not None:
            click.echo(f"  {output}_bowtie_pseudogene_raw.txt")
        if result.bowtie_genome_raw is not None:
            click.echo(f"  {output}_bowtie_genome_raw.txt")
        click.echo(f"\nProbes:")
        click.echo(format_probes_table(result))


@main.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option(
    '-p', '--pairs', 'n_pairs',
    default=30,
    type=int,
    help='Number of probe pairs to design (default: 30)'
)
@click.option(
    '-a', '--amplifier',
    default='B1',
    type=click.Choice(['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B9', 'B10', 'B11', 'B13', 'B14', 'B15', 'B17']),
    help='HCR amplifier (default: B1)'
)
@click.option(
    '-g', '--target-gibbs',
    default=-31.0,
    type=float,
    help='Target ΔG per 25-nt half in kcal/mol (default: -31.0)'
)
@click.option(
    '--allowable-gibbs',
    default='-35,-27',
    type=str,
    help='Strict ΔG range as min,max (default: -35,-27)'
)
@click.option(
    '--asymmetric-gibbs',
    is_flag=True,
    default=False,
    help='Enable asymmetric ΔG leniency (one half strict, one lenient)'
)
@click.option(
    '--lenient-gibbs-min',
    default=-42.0,
    type=float,
    help='Floor ΔG for lenient half (default: -42.0, only with --asymmetric-gibbs)'
)
@click.option(
    '--asymmetric-bowtie',
    is_flag=True,
    default=False,
    help='Enable asymmetric bowtie filtering (at least one half clean)'
)
@click.option(
    '--pair-spacing',
    default=2,
    type=int,
    help='Minimum gap between consecutive 52-nt pair blocks (default: 2)'
)
@click.option(
    '-s', '--species',
    default='human',
    type=click.Choice(['human', 'mouse', 'elegans', 'drosophila', 'rat']),
    help='Species for masking databases (default: human)'
)
@click.option(
    '--pseudogene-mask/--no-pseudogene-mask',
    default=False,
    help='Mask pseudogene alignments (default: off)'
)
@click.option(
    '--genome-mask/--no-genome-mask',
    default=False,
    help='Mask genome alignments (default: off)'
)
@click.option(
    '--index-dir',
    default=None,
    type=click.Path(exists=True, file_okay=False),
    help='Directory containing bowtie indexes'
)
@click.option(
    '--repeatmask-file',
    default=None,
    type=click.Path(exists=True),
    help='FASTA file with N\'s marking repeat regions'
)
@click.option(
    '--hp-threshold',
    default=5,
    type=click.IntRange(3, 10),
    help='Homopolymer run threshold (default: 5)'
)
@click.option(
    '--di-threshold',
    default=3,
    type=click.IntRange(3, 10),
    help='Dinucleotide repeat threshold (default: 3)'
)
@click.option(
    '--resolve-spacer',
    default=None,
    type=str,
    help='Resolve WW spacer to explicit bases (e.g. AA, TT). Default: output WW as IUPAC'
)
@click.option(
    '-o', '--output',
    default=None,
    type=str,
    help='Output file prefix (default: derived from input filename)'
)
@click.option(
    '--save-bowtie-raw',
    is_flag=True,
    default=False,
    help='Save raw bowtie alignment output'
)
@click.option(
    '--quiet', '-q',
    is_flag=True,
    help='Suppress output to stdout'
)
def hcrdesign(
    input_file: str,
    n_pairs: int,
    amplifier: str,
    target_gibbs: float,
    allowable_gibbs: str,
    asymmetric_gibbs: bool,
    lenient_gibbs_min: float,
    asymmetric_bowtie: bool,
    pair_spacing: int,
    species: str,
    pseudogene_mask: bool,
    genome_mask: bool,
    index_dir: str,
    repeatmask_file: str,
    hp_threshold: int,
    di_threshold: int,
    resolve_spacer: str,
    output: str,
    save_bowtie_raw: bool,
    quiet: bool,
):
    """Design HCR split-initiator probe pairs for a target sequence.

    INPUT_FILE is a FASTA file containing the target sequence.

    Each probe pair consists of two 25-nt probes (P1+P2) separated by a
    2-nt gap on the target, with half-initiator sequences for HCR signal
    amplification.

    Examples:

        probedesign hcrdesign input.fa --pairs 20 --amplifier B1

        probedesign hcrdesign input.fa -p 25 -a B2 --asymmetric-gibbs

        probedesign hcrdesign input.fa -p 30 -a B9 --resolve-spacer AA
    """
    from .hcr import design_hcr_probes
    from .hcr_output import write_hcr_output_files, format_hcr_oligos

    # Parse allowable Gibbs range
    try:
        min_gibbs, max_gibbs = map(float, allowable_gibbs.split(','))
    except ValueError:
        raise click.BadParameter(
            f"Invalid format for --allowable-gibbs: {allowable_gibbs}. "
            "Expected format: min,max (e.g., -35,-27)"
        )

    # Determine output prefix
    if output is None:
        output = Path(input_file).stem

    if not quiet:
        click.echo(f"Designing HCR probe pairs for {input_file}...")
        click.echo(f"  Amplifier: {amplifier}, Pairs: {n_pairs}, Spacing: {pair_spacing}")
        if asymmetric_gibbs:
            click.echo(f"  Asymmetric ΔG: lenient floor = {lenient_gibbs_min}")
        if asymmetric_bowtie:
            click.echo(f"  Asymmetric bowtie: enabled")

    result = design_hcr_probes(
        input_file=input_file,
        n_pairs=n_pairs,
        amplifier=amplifier,
        target_gibbs=target_gibbs,
        strict_range=(min_gibbs, max_gibbs),
        asymmetric_gibbs=asymmetric_gibbs,
        lenient_gibbs_min=lenient_gibbs_min,
        asymmetric_bowtie=asymmetric_bowtie,
        pair_spacing=pair_spacing,
        species=species,
        pseudogene_mask=pseudogene_mask,
        genome_mask=genome_mask,
        index_dir=index_dir,
        repeatmask_file=repeatmask_file,
        hp_threshold=hp_threshold,
        di_threshold=di_threshold,
        resolve_spacer=resolve_spacer,
        output_name=Path(output).name,  # -o may include a directory; names don't
        save_bowtie_raw=save_bowtie_raw,
    )

    if not result.pairs:
        click.echo("ERROR: Could not find any valid HCR probe pairs.", err=True)
        raise SystemExit(1)

    # Write output files
    write_hcr_output_files(result, output)

    # Print summary
    if not quiet:
        click.echo(f"\nFound {len(result.pairs)} HCR probe pairs (score: {result.score:.4f})")
        click.echo(f"\nOutput files:")
        click.echo(f"  {output}_HCR_oligos.txt")
        click.echo(f"  {output}_HCR_seq.txt")
        click.echo(f"  {output}_HCR_hits.txt")
        click.echo(f"\nProbe pairs:")
        click.echo(format_hcr_oligos(result))


@main.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option(
    '-l', '--oligo-length',
    default=20,
    type=int,
    help='Length of oligonucleotide to analyze (default: 20)'
)
def analyze(input_file: str, oligo_length: int):
    """Analyze a sequence for probe design feasibility.

    Shows the distribution of Gibbs free energies across the sequence.
    """
    from .fasta import read_fasta, sequences_to_single_string
    from .sequence import clean_sequence, has_invalid_chars
    from .thermodynamics import gibbs_rna_dna

    headers, seqs = read_fasta(input_file)
    seq = clean_sequence(sequences_to_single_string(seqs))

    click.echo(f"Sequence length: {len(seq)}")
    click.echo(f"Analyzing with oligo length: {oligo_length}")

    gibbs_values = []
    invalid_count = 0

    for i in range(len(seq) - oligo_length + 1):
        oligo = seq[i:i + oligo_length]
        if has_invalid_chars(oligo):
            invalid_count += 1
        else:
            try:
                gibbs_values.append(gibbs_rna_dna(oligo))
            except KeyError:
                invalid_count += 1

    if gibbs_values:
        click.echo(f"\nGibbs FE distribution:")
        click.echo(f"  Min: {min(gibbs_values):.2f} kcal/mol")
        click.echo(f"  Max: {max(gibbs_values):.2f} kcal/mol")
        click.echo(f"  Mean: {sum(gibbs_values)/len(gibbs_values):.2f} kcal/mol")

    click.echo(f"\nValid positions: {len(gibbs_values)}")
    click.echo(f"Invalid positions: {invalid_count}")


if __name__ == '__main__':
    main()
