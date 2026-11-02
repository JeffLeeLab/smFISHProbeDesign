# Bowtie Installation and Setup

This document describes how Bowtie and its indexes are set up for probe design masking (pseudogene mask, genome mask). **`setup_all.sh` does all of this automatically.** Read on only for manual setup or troubleshooting.

## Overview

Bowtie is a short-read aligner. ProbeDesign uses it to find sequences that align to pseudogenes or to many genomic locations, and masks those regions during probe design.

**Important Notes**:
- This project uses **Bowtie 1** (not Bowtie 2). Bowtie 1 is optimised for short reads (≤50 bp), which matches oligo lengths.
- Install bowtie from **bioconda**, NOT homebrew. Homebrew's `bowtie` is an unrelated package.
- Check you have the right one with `bowtie --version`. It should show "bowtie-align-s version 1.x.x".
- Bowtie 1 v1.3+ reads Bowtie 2 `.bt2` index files directly, which is why the genome indexes can be downloaded pre-built rather than built locally.

## Installation

`environment.yml` installs `bowtie=1.3.1` into the `probedesign` environment. Never install it into `base`. To add it manually to an existing environment:

```bash
mamba install -n probedesign -c conda-forge -c bioconda --channel-priority strict bowtie=1.3.1
# or: conda install ...
mamba activate probedesign
bowtie --version
```

The `defaults` channel is deliberately not used; see `environment.yml`.

## Index Setup

All indexes live in `<repo>/bowtie_indexes/`. That is the default for `--index-dir` and for the app's "Bowtie index directory" field.

### Pseudogene Indexes (`--pseudogene-mask`)

These are built locally from the FASTA files in `probedesign/pseudogeneDBs/`. `setup_all.sh` builds human, mouse and drosophila. Manual equivalent:

```bash
mkdir -p bowtie_indexes && cd bowtie_indexes
bowtie-build ../probedesign/pseudogeneDBs/human.fasta      humanPseudo
bowtie-build ../probedesign/pseudogeneDBs/mouse.fasta      mousePseudo
bowtie-build ../probedesign/pseudogeneDBs/drosophila.fasta drosophilaPseudo
bowtie-build ../probedesign/pseudogeneDBs/elegans.fasta    celegansPseudo   # optional, not built by setup
```

### Genome Indexes (`--genome-mask`)

`setup_all.sh` downloads pre-built Bowtie 2 indexes from `https://genome-idx.s3.amazonaws.com/bt/`, then moves and renames the files to the prefixes the code expects:

| Species | Source zip | Index prefix | Download |
|---------|------------|--------------|----------|
| human | `GRCh38_noalt_as.zip` | `GCA_000001405.15_GRCh38_no_alt_analysis_set` | ~3.5 GB |
| mouse | `mm10.zip` | `mm10` | ~3.1 GB |
| drosophila | `BDGP6.zip` | `drosophila` | ~176 MB |

Use `./setup_all.sh --skip-genome` to skip this download.

### Index Names

The code looks up these prefixes (`masking.py`):

| Species | Pseudogene | Genome | Built by setup |
|---------|------------|--------|----------------|
| human | `humanPseudo` | `GCA_000001405.15_GRCh38_no_alt_analysis_set` | yes |
| mouse | `mousePseudo` | `mm10` | yes |
| drosophila | `drosophilaPseudo` | `drosophila` | yes |
| elegans | `celegansPseudo` | `celegans` | no |
| rat | `ratPseudo` | `rat` | no |

### Expected Layout

```
bowtie_indexes/
├── humanPseudo.{1,2,3,4,rev.1,rev.2}.ebwt                                   # built (Bowtie 1)
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.{1,2,3,4,rev.1,rev.2}.bt2    # downloaded (Bowtie 2 format)
├── mousePseudo.*.ebwt, mm10.*.bt2
└── drosophilaPseudo.*.ebwt, drosophila.*.bt2
```

## Finding the Bowtie Executable

`find_bowtie()` in `masking.py` checks, in order:
1. `bowtie` on `PATH` (the activated `probedesign` env)
2. `/usr/bin/bowtie`
3. `$BOWTIEHOME/bowtie`
4. `~/bowtie/bowtie`

If none is found, masking prints `Warning: ... masking failed: Could not find bowtie executable` and the design **continues without that mask**. Always activate the environment first. If you call the env's Python directly, put `<env>/bin` on `PATH` or set `BOWTIEHOME=<env>/bin`.

## Testing the Installation

```bash
mamba activate probedesign
echo -e ">test\nACGTACGTACGTACGT" | bowtie -f -v 0 -k 1 bowtie_indexes/humanPseudo -
```

## Troubleshooting

### "command not found: bowtie"
- Activate the environment: `mamba activate probedesign`
- Check it's installed: `mamba list -n probedesign bowtie`

### "Could not find bowtie executable"
- See [Finding the Bowtie Executable](#finding-the-bowtie-executable).

### Index not found
- Check the files are in `bowtie_indexes/` (or the directory given to `--index-dir`).
- The index name is the prefix before `.1.ebwt` / `.1.bt2`, and must match the table above.

## References

- Bowtie homepage: http://bowtie-bio.sourceforge.net/
- Bowtie manual: http://bowtie-bio.sourceforge.net/manual.shtml
- Bioconda recipe: https://bioconda.github.io/recipes/bowtie/README.html
