# ProbeDesign

> Design oligonucleotide probes for single molecule RNA FISH experiments.

This repo is a fork of https://github.com/arjunrajlaboratory/ProbeDesign from the Raj Lab.

ProbeDesign creates oligonucleotide probes for single molecule RNA FISH experiments. It selects optimal probe sequences based on:

- **Thermodynamic properties** — targets specific Gibbs free energy for RNA-DNA hybridization (Sugimoto 1995)
- **Sequence masking** — avoids cross-hybridization to pseudogenes, repeats, and repetitive genomic regions
- **Optimal spacing** — uses dynamic programming to find the best probe placements
- **Mixed-length probes** — variable-length oligos (e.g. 18-22 bp) to maximise probe count on GC-variable sequences
- **Low-complexity filtering** — excludes probes with homopolymer or dinucleotide repeats (L-mask)

## Development Status

**Completed** 
- ✅ Create streamlit GUI application
- ✅ Batch mode for CLI and GUI
- ✅ Implement homopolymer/dinucleotide filter
- ✅ Implement automated mixed-length probe design compatible with dynamic programming
- ✅ Allow bowtie result output
- ✅ Implement HCRv3 probe design algorithms from the RShiny scripts
- ✅ Host test streamlit application on Glasgow VM
- ✅ Implement Isoform Explorer — browse constitutive sequences and download FASTA for probe design (Ensembl v113, Hsap/Mmus/Dmel)


**In Development / Planned**
- [ ] Test local RepeatMasker installation
- [ ] Implement hairpin structure & inter-probe duplex formation filters

---

## Table of Contents

1. [Quick Start](#1-quick-start)
2. [Installation](#2-installation)
3. [Web App (Streamlit)](#3-web-app-streamlit)
4. [Command-Line Interface](#4-command-line-interface)
5. [Masking Options](#5-masking-options)
6. [Output Files](#6-output-files)
7. [Troubleshooting](#7-troubleshooting)
8. [Project Structure](#8-project-structure)
9. [Further Reading](#9-further-reading)

---

## 1. Quick Start

```bash
# 1. Install a package manager (mamba or conda)
#    See Section 2 for details

# 2. Clone and set up
git clone https://github.com/JYLeeSci/smFISHProbeDesign.git
cd smFISHProbeDesign
chmod +x setup_all.sh
./setup_all.sh          # ~18 min first run (downloads genome indices)

# 3. Launch the web app
mamba activate probedesign   # or: conda activate probedesign
streamlit run streamlit_app/app.py

# 4. Or use the CLI
probedesign design input.fa --probes 48 --pseudogene-mask --genome-mask
```

---

## 2. Installation

### Platform Support

| OS | Supported | Notes |
|----|-----------|-------|
| **macOS** (Apple Silicon / Intel) | ✅ | |
| **Linux** (x86_64 / aarch64) | ✅ | |
| **Windows** | ⚠️ WSL2 only | `wsl --install` in PowerShell, then use Ubuntu terminal |

### Step 1 — Install a Package Manager

You need **one** of `mamba` or `conda`. The setup script auto-detects whichever is available.

**Option A — Miniforge** (recommended):

1. Download from [github.com/conda-forge/miniforge/releases/latest](https://github.com/conda-forge/miniforge/releases/latest)
2. Run: `bash ~/Downloads/Miniforge3-*.sh`
3. Restart terminal, verify: `conda --version`

### Step 2 — Clone and Run Setup

```bash
git clone https://github.com/JYLeeSci/smFISHProbeDesign.git
cd smFISHProbeDesign
chmod +x setup_all.sh
./setup_all.sh
```

The script will:
1. Create a `probedesign` conda environment (Python 3.11 + Bowtie 1.3.1)
2. Install the `probedesign` package + Streamlit
3. Build pseudogene alignment indices (~1 min)
4. Download genome alignment indices (~12 min, ~6.7 GB)
5. Run validation tests

To skip the 6.7 GB genome download: `./setup_all.sh --skip-genome`

### Manual Setup (Advanced)

```bash
# Create environment from environment.yml
mamba env create -f environment.yml   # or: conda env create -f environment.yml
mamba activate probedesign
pip install -e .

# Update an existing environment
mamba env update -n probedesign -f environment.yml --prune

# Remove and recreate from scratch
mamba env remove -n probedesign
mamba env create -f environment.yml
```

Without `environment.yml`:

```bash
mamba create -n probedesign \
    -c conda-forge -c bioconda --channel-priority strict \
    python=3.11 "click>=8.0" bowtie=1.3.1 -y
mamba activate probedesign
pip install -e .
pip install "streamlit>=1.32" "pandas>=1.5"
```

---

## 3. Web App (Streamlit)

The primary interface for probe design is a Streamlit web application.

### Launch

```bash
mamba activate probedesign   # or conda
streamlit run streamlit_app/app.py
# Opens http://localhost:8501
```

### Single Mode

1. Upload or paste a FASTA file
2. Adjust parameters in the sidebar (probes, oligo length, masking)
3. Click **Design Probes**
4. Download `_oligos.txt` and `_seq.txt`

### Batch Mode

1. Switch to **Batch** mode
2. Provide a folder of FASTA files (or upload multiple files)
3. Set an output directory
4. Click **Run Batch** — results are written per-file with a summary TSV

### Sidebar Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Number of probes | 48 | Target probe count (no upper limit) |
| Oligo length mode | Fixed (20 bp) | Fixed length or mixed range (e.g. 18-22 bp) |
| Spacer length | 2 bp | Minimum gap between probes (no upper limit) |
| Target Gibbs FE | −23.0 kcal/mol | Optimal binding free energy |
| Allowable Gibbs range | −26 to −20 | Probes outside this are excluded |
| Species | human | Genome/pseudogene database |
| Pseudogene mask | on | Filter probes matching pseudogenes |
| Genome mask | on | Filter repetitive genome regions |
| RepeatMasker | off | Auto-detect repeats (requires extra setup) |

### Mixed-Length Probes

Select **Mixed range** under "Oligo length mode" to use variable-length probes (e.g. 18-22 bp). This maximises probe count on sequences with variable GC content — GC-rich regions get shorter probes, AT-rich regions get longer probes, all within the same Gibbs FE target range.

### Isoform Explorer

The **Isoform Explorer** tab (top navigation bar) lets you inspect the isoform architecture of any gene and identify whether a constitutive sequence exists — a region shared by all isoforms — for use as probe-design input.

**Workflow**:
1. Select species (Human, Mouse, or Drosophila) and search for a gene by name or Ensembl ID.
2. The isoform structure plot shows all transcripts coloured by Transcript Support Level (TSL) for mammals. Isoforms excluded by the adaptive TSL filter are de-emphasised.
3. If a constitutive sequence exists, copy it from the sequence box or download the FASTA.
4. Paste the FASTA into the **smFISH** or **HCR** designer — the N markers at exon–exon junctions prevent probes from spanning junctions.

**Data**: Ensembl v113 · Hsap 55,184 genes · Mmus 54,705 genes · Dmel 16,797 genes.
Mammalian data uses an adaptive TSL filter (TSL ≤ 2 where available, otherwise the best available level). The FASTA header records the per-gene TSL cutoff applied.

### TransQuant

The **TransQuant** tab (top navigation bar) embeds the [TransQuant probe-weight tool](https://github.com/JeffLeeLab/TransQuant-probe-weight) (external `transquant_w` package). Given the full genomic locus of a gene and a probe set, it computes the **probe weight factor W** (corrects transcription-site brightness for where probes bind relative to Pol II), the **gene length L**, and the **probe localisation profile N** — following the TransQuant methodology (Bahar Halpern & Itzkovitz, *Methods* 2016). It is fully self-contained — no inputs or results are shared with the smFISH / HCR / Isoform Explorer tabs.

A standalone, browser-only version (no install needed) is hosted at [jeffleelab.github.io/TransQuant-probe-weight](https://jeffleelab.github.io/TransQuant-probe-weight/).

The embedded version is pinned to a release tag in `environment.yml` (`git+https://github.com/JeffLeeLab/TransQuant-probe-weight.git@v0.2.0`). To update it: bump the `@<tag>` in `environment.yml`, then re-run `./setup_all.sh` or `pip install --upgrade "git+https://github.com/JeffLeeLab/TransQuant-probe-weight.git@<new-tag>"` inside the `probedesign` environment.

---

## 4. Command-Line Interface

### `probedesign design`

```
Usage: probedesign design [OPTIONS] INPUT_FILE

  Design probes for a target sequence.
  INPUT_FILE is a FASTA file containing the target sequence.

Options:
  -n, --probes INTEGER            Number of probes to design  [default: 48]
  -l, --oligo-length TEXT         Oligo length: single int (e.g. 20) or range
                                  for mixed lengths (e.g. 18-22)  [default: 20]
  -s, --spacer-length INTEGER     Minimum gap between probes  [default: 2]
  -g, --target-gibbs FLOAT        Target Gibbs free energy in kcal/mol  [default: -23]
  --allowable-gibbs TEXT          Allowable Gibbs FE range as min,max  [default: -26,-20]
  -o, --output TEXT               Output file prefix (default: derived from input filename)
  -q, --quiet                     Suppress output to stdout
  --species [human|mouse|elegans|drosophila|rat]
                                  Species for masking databases  [default: human]
  --pseudogene-mask / --no-pseudogene-mask
                                  Mask regions that align to pseudogenes  [default: off]
  --genome-mask / --no-genome-mask
                                  Mask repetitive genomic regions  [default: off]
  --index-dir DIRECTORY           Directory containing bowtie indexes (default: auto-detect)
  --repeatmask-file PATH          FASTA file with N's marking repeat regions (manual masking)
  --repeatmask / --no-repeatmask  Run RepeatMasker automatically to mask repeats  [default: off]
  --hp-threshold INTEGER          Homopolymer repeat threshold for L-mask (probes with
                                  single-nucleotide runs >= this length are excluded)
                                  [default: 5; range: 3-10]
  --di-threshold INTEGER          Dinucleotide repeat threshold for L-mask (probes with
                                  dinucleotide motifs repeated >= this many times are
                                  excluded)  [default: 3; range: 3-10]
  --save-bowtie-raw               Save raw bowtie alignment output (produces large files)
  --help                          Show this message and exit.
```

### `probedesign analyze`

Check the Gibbs FE distribution of a sequence before designing probes — useful for choosing
appropriate `--target-gibbs` and `--allowable-gibbs` values.

```
Usage: probedesign analyze [OPTIONS] INPUT_FILE

Options:
  -l, --oligo-length INTEGER  Length of oligonucleotide to analyze  [default: 20]
  --help                      Show this message and exit.
```

### Examples

```bash
# Check Gibbs FE distribution before designing
probedesign analyze input.fa

# Basic: 32 probes, no masking
probedesign design input.fa --probes 32

# Recommended: pseudogene + genome masking
probedesign design input.fa --probes 48 --pseudogene-mask --genome-mask

# Mixed-length probes (18-22 bp) — maximises count on GC-variable sequences
probedesign design input.fa --probes 48 -l 18-22

# Repeat masking from a pre-masked file (N's mark repeats)
probedesign design input.fa --probes 32 --repeatmask-file input_masked.fa

# HCR probes (52 bp, adjusted thermodynamic targets)
probedesign design input.fa --probes 20 -l 52 \
  --target-gibbs -60 --allowable-gibbs -80,-40

# Mouse gene, all masking, custom output prefix
probedesign design my_gene.fa --probes 48 --species mouse \
  --pseudogene-mask --genome-mask -o MyGene_probes
```

---

## 5. Masking Options

Masking removes probe candidates that would hybridize non-specifically.

### Repeat Masking (R Mask)

| Method | Flag | Use case |
|--------|------|----------|
| Auto-detect | *(none)* | FASTA already has `N`s at repeats |
| Manual file | `--repeatmask-file PATH` | Pre-masked sequence file |
| Automatic | `--repeatmask` | Runs RepeatMasker locally (see [docs/REPEATMASKER.md](docs/REPEATMASKER.md)) |

### Pseudogene Masking (P Mask) — `--pseudogene-mask`

16-mer exact alignment against curated pseudogene databases. Indices built by `setup_all.sh`.

### Genome Masking (B Mask) — `--genome-mask`

Three bowtie searches (12/14/16-mers) against the whole genome with tiered hit-count thresholds.

### Low-Complexity Filtering (L Mask) — always active

Probes containing homopolymer runs (e.g. `AAAAA`) or dinucleotide repeats (e.g. `ACACAC`) are excluded. Thresholds are adjustable via `--hp-threshold` (default 5) and `--di-threshold` (default 3).

### Supported Species

| Species | `--species` | Pseudogene | Genome |
|---------|-------------|------------|--------|
| Human | `human` | ✅ | ✅ |
| Mouse | `mouse` | ✅ | ✅ |
| Drosophila | `drosophila` | ✅ | ✅ |
| C. elegans | `elegans` | Manual | Manual |
| Rat | `rat` | Manual | Manual |

---

## 6. Output Files

### `<name>_oligos.txt`

Tab-separated, one probe per line. Order sequences from the last two columns when placing synthesis orders.

| Column | Content |
|--------|---------|
| 1 — Index | Probe number (1-based) |
| 2 — Start | Start position in target sequence (1-based, nucleotides) |
| 3 — GC% | GC content |
| 4 — Tm | Melting temperature (°C) |
| 5 — Gibbs | ΔG (kcal/mol) |
| 6 — Sequence | Probe sequence (5'→3') |
| 7 — Name | Probe label |

### `<name>_seq.txt`

Visual sequence map (110 chars per line) showing probe placement and masking. Useful for checking masking behaviour.

Each block shows:
- **Sequence line**: Original target sequence
- **R line** (when repeat masking used): `R` = masked position, base = unmasked
- **P line** (when pseudogene mask used): `P` = masked
- **B line** (when genome mask used): `B` = masked
- **F line** (always): `F` = thermodynamically excluded, base = valid probe start
- **Probe complement line**: reverse-complement sequence placed below the target
- **Probe label line**: `Prb# N, Pos X, FE Y, GC Z` for each probe

### Additional files (when masking is enabled)

| File | Content |
|------|---------|
| `<name>_bowtie_pseudogene.txt` | Per-position 16-mer hit counts from pseudogene search |
| `<name>_bowtie_genome.txt` | Per-position hit counts from genome search (12/14/16-mers) |
| `<name>_bowtie_*_raw.txt` | Raw bowtie alignment output (only with `--save-bowtie-raw`) |

---

## 7. Troubleshooting

| Problem | Solution |
|---------|----------|
| `mamba / conda: command not found` | Restart terminal after installing; re-run installer if needed |
| `bowtie: command not found` | Activate environment: `mamba activate probedesign` |
| Wrong bowtie (file manager, not aligner) | Only install via bioconda, not homebrew |
| Setup fails at pseudogene indices | Check `probedesign/pseudogeneDBs/*.fasta` exist |
| Genome download interrupted | Re-run `./setup_all.sh` (resumes automatically) |
| `streamlit: command not found` | Activate environment first |
| Port 8501 in use | `streamlit run streamlit_app/app.py --server.port 8502` |
| Windows | Use WSL2 Ubuntu terminal, not PowerShell |

---

## 8. Project Structure

```
smFISHProbeDesign/
├── src/probedesign/          # Python package
│   ├── cli.py                #   Click-based CLI
│   ├── core.py               #   Probe design algorithm (badness, DP, mixed-length)
│   ├── thermodynamics.py     #   Gibbs FE and Tm (Sugimoto 1995)
│   ├── masking.py            #   Bowtie masking (pseudogene, genome, RepeatMasker)
│   ├── sequence.py           #   Sequence utilities
│   ├── fasta.py              #   FASTA I/O
│   └── output.py             #   Output file generation
├── streamlit_app/            # Streamlit web app
│   ├── app.py                #   Main UI
│   └── utils.py              #   Backend helpers
├── probedesign/              # Legacy MATLAB code
│   ├── findprobesLocal.m     #   Original MATLAB implementation
│   └── pseudogeneDBs/        #   Pseudogene FASTA files
├── docs/                     # Documentation
│   ├── BOWTIE.md             #   Bowtie setup guide
│   ├── REPEATMASKER.md       #   RepeatMasker setup guide
│   ├── MIXED_LENGTH_PLAN.md  #   Mixed-length design plan
│   ├── probe_design_principles.md  # Algorithm details
│   └── summary.md            #   Development notes
├── test_cases/               # Validation test cases
├── bowtie_indexes/           # Genome/pseudogene indices (gitignored)
├── setup_all.sh              # Automated setup script
├── environment.yml           # Conda environment specification
├── pyproject.toml            # Python package config
├── CLAUDE.md                 # AI coding context
└── README.md                 # This file
```

---

## 9. Further Reading

| Document | Description |
|----------|-------------|
| [docs/probe_design_principles.md](docs/probe_design_principles.md) | Thermodynamic model, masking algorithm, DP optimization |
| [docs/summary.md](docs/summary.md) | Development notes, Bowtie 1/2 index compatibility |
| [docs/BOWTIE.md](docs/BOWTIE.md) | Bowtie installation and index setup |
| [docs/REPEATMASKER.md](docs/REPEATMASKER.md) | RepeatMasker installation (~56 GB database) |
| [docs/MIXED_LENGTH_PLAN.md](docs/MIXED_LENGTH_PLAN.md) | Mixed-length probe design plan and analysis |
