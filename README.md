# ProbeDesign

> Design oligonucleotide probes for single molecule RNA FISH (smFISH) and HCR v3 experiments.

This repo is a fork of https://github.com/arjunrajlaboratory/ProbeDesign from the Raj Lab.

ProbeDesign creates oligonucleotide probes for smFISH and HCR v3 split-initiator FISH. It selects optimal probe sequences based on:

- **Thermodynamic properties** — targets specific Gibbs free energy for RNA-DNA hybridization (Sugimoto 1995)
- **Sequence masking** — avoids cross-hybridization to pseudogenes, repeats, and repetitive genomic regions
- **Optimal spacing** — uses dynamic programming to find the best probe placements
- **Mixed-length probes** — variable-length oligos (e.g. 18-22 bp) to maximise probe count on GC-variable sequences
- **Low-complexity filtering** — excludes probes with homopolymer or dinucleotide repeats (L-mask)
- **HCR v3 split-initiator pairs** — two 25-nt halves per pair, 13 amplifiers (B1–B17)

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
- ✅ IDT 96-well plate order export — per-run button, IDT Plate tab to combine probe sets, `probedesign idtplate`


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
git clone https://github.com/JeffLeeLab/smFISHProbeDesign.git
cd smFISHProbeDesign
chmod +x setup_all.sh
./setup_all.sh          # ~18 min first run (downloads genome indices)

# 3. Launch the web app
mamba activate probedesign   # or: conda activate probedesign
streamlit run streamlit_app/app.py

# 4. Or use the CLI
probedesign design input.fa --probes 48 --pseudogene-mask --genome-mask      # smFISH
probedesign hcrdesign input.fa --pairs 30 -a B1 --pseudogene-mask --genome-mask  # HCR
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
git clone https://github.com/JeffLeeLab/smFISHProbeDesign.git
cd smFISHProbeDesign
chmod +x setup_all.sh
./setup_all.sh
```

The script will:
1. Create a `probedesign` conda environment (Python 3.11 + Bowtie 1.3.1)
2. Install the `probedesign` package and the web-app dependencies (Streamlit, pandas, pyarrow, matplotlib, TransQuant)
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
    python=3.11 "click>=8.0" bowtie=1.3.1 "pytest>=7.0" -y
mamba activate probedesign
pip install -e .
pip install "streamlit>=1.48" "pandas>=1.5" "pyarrow>=14.0" "matplotlib>=3.8" "openpyxl>=3.1"
pip install "git+https://github.com/JeffLeeLab/TransQuant-probe-weight.git@v0.2.0"
```

### Running Tests

```bash
python -m pytest tests/      # HCR unit, smFISH regression, run-header and IDT plate tests
```

`tests/test_smfish_regression.py` reproduces the reference designs in `tests/regression/` (CDKN1A, KRT19, mouseMITF). The KRT19 and mouseMITF tests need bowtie and the human/mouse indices, and are skipped without them. `setup_all.sh` also checks CDKN1A and KRT19 against the same references. `tests/test_provenance.py` covers the `#` run header in `_seq.txt` / `_HCR_seq.txt`, including that its command reproduces the design.

---

## 3. Web App (Streamlit)

The primary interface for probe design is a Streamlit web application with four tabs in the top navigation bar: **smFISH**, **HCR**, **Isoform Explorer** and **TransQuant**.

### Launch

```bash
mamba activate probedesign   # or conda
streamlit run streamlit_app/app.py
# Opens http://localhost:8501
```

### smFISH — Single Mode

1. Upload or paste a FASTA file
2. Adjust parameters in the sidebar (probes, oligo length, masking)
3. Click **Design Probes**
4. Download `_oligos.txt`, `_seq.txt`, the bowtie hit tables (when masking is on), and raw bowtie files (when **Save raw bowtie output** is ticked)
5. **Download IDT plate (.xlsx)** gives an IDT 96-well plate order for this probe set alone. To put several sets on one order, use the [IDT Plate](#idt-plate-builder) tab

### smFISH — Batch Mode

1. Switch **Mode** to **Batch**
2. Provide a folder of FASTA files (or upload multiple files)
3. Click **Run Batch**, then download all results as a ZIP plus `batch_summary.tsv`. Its `Warnings` column lists any mask that failed for a gene.
4. Optional: set an output directory under **Auto-save outputs** to also write results to disk on the server

### smFISH Sidebar Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Mode | Single sequence | Single sequence or Batch |
| Species | human | Genome/pseudogene database |
| Number of probes | 48 | Target probe count (no upper limit) |
| Oligo length mode | Fixed (20 bp) | Fixed length or mixed range (e.g. 18-22 bp); 10–60 bp |
| Spacer length | 2 bp | Minimum gap between probes (no upper limit) |
| Target Gibbs FE | −23.0 kcal/mol | Optimal binding free energy |
| Allowable Gibbs range | −26 to −20 | Probes outside this are excluded |
| Homopolymer / Dinucleotide repeats | 5 / 3 | Low-complexity (L-mask) thresholds, 3–10 |
| Pseudogene mask | on | Filter probes matching pseudogenes |
| Genome mask | on | Filter repetitive genome regions |
| Bowtie index directory | `<repo>/bowtie_indexes` | Shown when a bowtie mask is on |
| Save raw bowtie output | off | Add raw alignment files to the downloads |
| RepeatMasker | None | None, Auto (requires extra setup), or provide a masked FASTA |

### Mixed-Length Probes

Select **Mixed range** under "Oligo length mode" to use variable-length probes (e.g. 18-22 bp). This maximises probe count on sequences with variable GC content — GC-rich regions get shorter probes, AT-rich regions get longer probes, all within the same Gibbs FE target range.

### HCR

The **HCR** tab designs HCR v3 split-initiator probe pairs. Each pair is two 25-nt probes separated by a 2-nt gap on the target, with half-initiator sequences attached. The two halves hybridise independently, so ΔG targets apply **per 25-nt half**. Single and batch modes work as for smFISH. Downloads are `_HCR_oligos.txt`, `_HCR_seq.txt`, `_HCR_hits.txt` and an IDT plate order (`.xlsx`, full oligos with initiators); batch mode adds a ZIP and `hcr_batch_summary.tsv`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| Amplifier | B1 | B1–B5, B7, B9–B11, B13–B15, B17 (B7+ use a `WW` degenerate spacer) |
| Number of pairs | 30 | Target pair count |
| Pair spacing | 2 nt | Minimum gap between consecutive 52-nt pair blocks |
| Target ΔG per half | −31.0 kcal/mol | Optimal binding free energy per half |
| Strict ΔG range | −35 to −27 | Range each half must satisfy |
| Asymmetric ΔG leniency | off | Let one half go down to the **Lenient ΔG floor** (−42) |
| Homopolymer / Dinucleotide repeats | 5 / 3 | Low-complexity (L-mask) thresholds |
| Pseudogene / Genome mask | on / on | Bowtie off-target filters |
| Asymmetric bowtie | off | Only require one half to be free of off-target hits |
| Resolve WW spacer | blank | Replace `WW` with explicit bases (e.g. `AA`) for B7+ |

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

The embedded version is pinned to a release tag (currently `v0.2.0`). To update it, bump the `@<tag>` in all three places it appears — `environment.yml`, `setup_all.sh` and the install hint in `streamlit_app/pages/transquant.py` — then run `pip install --upgrade "git+https://github.com/JeffLeeLab/TransQuant-probe-weight.git@<new-tag>"` inside the `probedesign` environment (or re-run `./setup_all.sh`).

### IDT Plate Builder

The **IDT Plate** tab combines probe sets from several runs into one IDT 96-well plate order (the Excel file IDT's plate upload accepts).

1. Upload `_oligos.txt` / `_HCR_oligos.txt` files, or a batch results ZIP. smFISH and HCR sets can share a plate.
2. In the table, set the **Order** sets go onto plates, untick sets to leave out, and check the **Gene** column. Gene is the probe-set name without its region (`Eaat1-cmm` → `Eaat1`, `TDP-43-CDS-DmelCO` → `TDP-43`); edit it if the region wasn't stripped correctly.
3. Set **Minimum oligos per plate** in the sidebar: 24 for standard 96-well plates (default), 48 for HOTPlates.
4. Check the plate maps and download `IDT_plate-ordering_<YYYYMMDD>_<gene>_<gene>….xlsx`.

Layout rules:
- Wells fill down each column: A1–H1, then A2 …
- Each probe set starts on a new column.
- A set is never split across plates; one that doesn't fit moves to the next plate. Only a set larger than 96 oligos is split, over whole plates.
- Each plate is a worksheet named after its genes (Excel's 31-character limit; extra genes show as `+N`).
- Plates are filled in order with the fewest plates that each meet the minimum, earlier plates first. If no layout meets it (e.g. a single 16-oligo set), the app warns.

The same rules apply to the single-run **Download IDT plate** buttons and to `probedesign idtplate`.

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
                                  for mixed lengths (e.g. 18-22; range bounds
                                  10-60)  [default: 20]
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
  --index-dir DIRECTORY           Directory containing bowtie indexes
                                  [default: <repo>/bowtie_indexes]
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

### `probedesign hcrdesign`

Designs HCR v3 split-initiator probe pairs. ΔG values are **per 25-nt half**. Don't use `design -l 52` for HCR.

```
Usage: probedesign hcrdesign [OPTIONS] INPUT_FILE

Options:
  -p, --pairs INTEGER             Number of probe pairs to design  [default: 30]
  -a, --amplifier [B1|B2|B3|B4|B5|B7|B9|B10|B11|B13|B14|B15|B17]
                                  HCR amplifier  [default: B1]
  -g, --target-gibbs FLOAT        Target ΔG per 25-nt half  [default: -31.0]
  --allowable-gibbs TEXT          Strict ΔG range as min,max  [default: -35,-27]
  --asymmetric-gibbs              One half strict, one lenient
  --lenient-gibbs-min FLOAT       Floor ΔG for the lenient half  [default: -42.0]
  --asymmetric-bowtie             Require only one half to be bowtie-clean
  --pair-spacing INTEGER          Minimum gap between 52-nt pair blocks  [default: 2]
  -s, --species [human|mouse|elegans|drosophila|rat]
                                  Species for masking databases  [default: human]
  --pseudogene-mask / --no-pseudogene-mask   [default: off]
  --genome-mask / --no-genome-mask           [default: off]
  --index-dir DIRECTORY           [default: <repo>/bowtie_indexes]
  --repeatmask-file PATH          FASTA file with N's marking repeat regions
  --hp-threshold INTEGER          Homopolymer threshold  [default: 5; range: 3-10]
  --di-threshold INTEGER          Dinucleotide threshold  [default: 3; range: 3-10]
  --resolve-spacer TEXT           Replace WW spacer with explicit bases (e.g. AA)
  -o, --output TEXT               Output file prefix
  --save-bowtie-raw               Save raw bowtie alignment output
  -q, --quiet                     Suppress output to stdout
```

Note that `-s` means `--species` here but `--spacer-length` in `design`.

### `probedesign analyze`

Check the Gibbs FE distribution of a sequence before designing probes — useful for choosing
appropriate `--target-gibbs` and `--allowable-gibbs` values.

```
Usage: probedesign analyze [OPTIONS] INPUT_FILE

Options:
  -l, --oligo-length INTEGER  Length of oligonucleotide to analyze  [default: 20]
  --help                      Show this message and exit.
```

### `probedesign idtplate`

Build an IDT 96-well plate order (`.xlsx`) from oligo files of one or more runs, with the same layout rules as the [IDT Plate tab](#idt-plate-builder). Files go onto plates in the order given.

```
Usage: probedesign idtplate [OPTIONS] OLIGO_FILES...

Options:
  -o, --output-dir DIRECTORY     Directory for the plate file (default: current directory)
  --min-oligos INTEGER RANGE     IDT minimum oligos per plate: 24 for standard 96-well
                                 plates, 48 for HOTPlates (default: 24)
```

`OLIGO_FILES` are `_oligos.txt` / `_HCR_oligos.txt` files or batch-result ZIPs. Needs `openpyxl` (in `environment.yml`; or `pip install -e ".[idt]"`).

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

# Mouse gene, all masking, custom output prefix
probedesign design my_gene.fa --probes 48 --species mouse \
  --pseudogene-mask --genome-mask -o MyGene_probes

# HCR: 30 pairs for amplifier B1 with masking
probedesign hcrdesign input.fa --pairs 30 -a B1 --pseudogene-mask --genome-mask

# HCR: relax ΔG on one half when a target is hard to cover
probedesign hcrdesign input.fa -p 25 -a B2 --asymmetric-gibbs

# HCR: B7+ amplifier with WW spacer resolved to AA
probedesign hcrdesign input.fa -p 30 -a B9 --resolve-spacer AA

# IDT plate order for two probe sets (HOTPlates: 48-oligo minimum)
probedesign idtplate Eaat1-cmm_oligos.txt Eaat2-cmm_oligos.txt --min-oligos 48
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

`setup_all.sh` builds or downloads indices for human, mouse and drosophila. For elegans or rat, build the indices yourself into `bowtie_indexes/`, named `celegansPseudo` / `celegans` and `ratPseudo` / `rat`.

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

The file starts with a run header: lines beginning with `#`. It records:
- the equivalent `probedesign` command;
- the input and every design parameter;
- the version and git commit;
- for each mask, whether it was applied, which bowtie index was used, and a short fingerprint of that index (so runs on different machines can be compared);
- the probes found and the percentage covered by each mask.

A mask that was requested but failed (bowtie missing, a broken index) reads `FAILED, not applied`. `grep -v '^#'` strips the header and leaves the sequence map unchanged.

Each block shows:
- **Sequence line**: Original target sequence
- **R line** (when repeat masking used): `R` = masked position, base = unmasked
- **P line** (when pseudogene mask used): `P` = masked
- **B line** (when genome mask used): `B` = masked
- **F line** (always): `F` = thermodynamically excluded, base = valid probe start
- **L line** (always): `L` = probe start excluded by the low-complexity filter
- **Probe complement line**: reverse-complement sequence placed below the target
- **Probe label line**: `Prb# N,Pos X,FE Y,GC Z` for each probe

### Additional files (when masking is enabled)

| File | Content |
|------|---------|
| `<name>_bowtie_pseudogene.txt` | Per-position 16-mer hit counts from pseudogene search |
| `<name>_bowtie_genome.txt` | Per-position hit counts from genome search (12/14/16-mers) |
| `<name>_bowtie_*_raw.txt` | Raw bowtie alignment output (only with `--save-bowtie-raw`) |

### HCR outputs (`hcrdesign` / HCR tab)

| File | Content |
|------|---------|
| `<name>_HCR_oligos.txt` | TSV with header: `pair, half (P1/P2), start, GC%, Tm, Gibbs, strict (yes/no), full_oligo, name`. `full_oligo` includes the half-initiator and spacer. Names are `<name>_HCR<amp>_NN` (odd = P1, even = P2) |
| `<name>_HCR_seq.txt` | `#` run header (as for `_seq.txt`), then the sequence map with mask lines and pair annotations |
| `<name>_HCR_hits.txt` | Per-pair, per-half bowtie hit counts with strict/lenient annotation |
| `<name>_HCR_bowtie_*_raw.txt` | Raw bowtie output (only with `--save-bowtie-raw`; CLI only) |

### IDT plate order (`IDT_plate-ordering_<YYYYMMDD>_<genes>.xlsx`)

Excel workbook for IDT's 96-well plate upload. One worksheet per plate, named after its genes, with columns `Well Position`, `Name`, `Sequence`; empty wells are left out. Written by the **Download IDT plate** buttons, the [IDT Plate tab](#idt-plate-builder) and `probedesign idtplate`.

---

## 7. Troubleshooting

| Problem | Solution |
|---------|----------|
| `mamba / conda: command not found` | Restart terminal after installing; re-run installer if needed |
| `bowtie: command not found` | Activate environment: `mamba activate probedesign` |
| `Warning: Pseudogene/Genome masking failed: Could not find bowtie executable` | Bowtie isn't on `PATH`, so the design ran **without** that mask. Activate the environment and re-run |
| Wrong bowtie (file manager, not aligner) | Only install via bioconda, not homebrew |
| Setup fails at pseudogene indices | Check `data/pseudogeneDBs/*.fasta` exist |
| Genome download interrupted | Re-run `./setup_all.sh`. Species already downloaded are skipped; the interrupted one restarts |
| `streamlit: command not found` | Activate environment first |
| Port 8501 in use | `streamlit run streamlit_app/app.py --server.port 8502` |
| Windows | Use WSL2 Ubuntu terminal, not PowerShell |

---

## 8. Project Structure

```
smFISHProbeDesign/
├── src/probedesign/          # Python package
│   ├── cli.py                #   Click CLI: design, hcrdesign, analyze, idtplate
│   ├── core.py               #   smFISH algorithm (badness, DP, mixed-length)
│   ├── hcr.py                #   HCR split-initiator pair design
│   ├── thermodynamics.py     #   Gibbs FE and Tm (Sugimoto 1995)
│   ├── masking.py            #   Bowtie masking (pseudogene, genome), RepeatMasker, low-complexity
│   ├── sequence.py           #   Sequence utilities
│   ├── fasta.py              #   FASTA I/O
│   ├── output.py             #   smFISH output files
│   ├── hcr_output.py         #   HCR output files
│   └── idt.py                #   IDT 96-well plate order (.xlsx)
├── streamlit_app/            # Streamlit web app
│   ├── app.py                #   Entry point + top navigation
│   ├── pages/                #   smfish.py, hcr.py, isoform_explorer.py, transquant.py, idt_plate.py
│   ├── utils.py              #   Backend helpers (run/batch wrappers, validation)
│   └── isoform_data.py       #   Isoform Explorer data access + plot
├── data/
│   ├── isoform_explorer/     #   Isoform Explorer data (Ensembl v113, Parquet)
│   └── pseudogeneDBs/        #   Pseudogene FASTA files (indices built from these)
├── docs/                     # Documentation (see Further Reading)
├── tests/                    # pytest suite (HCR unit, smFISH regression, run header, IDT plate)
│   └── regression/           #   Regression inputs + reference outputs
├── legacy/                   # Archived upstream MATLAB / Python 2 code (not used)
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
| [docs/hcr_designer_plan.md](docs/hcr_designer_plan.md) | HCR split-initiator designer plan |
| [docs/nucleotide_repeat_plan.md](docs/nucleotide_repeat_plan.md) | Low-complexity (L-mask) filter plan |
| [docs/isoform_explorer/](docs/isoform_explorer/) | Isoform Explorer plan and data-build scripts |

The `*_plan.md` files are design plans for features that have already shipped. They're kept for background, so where they disagree with the code, the code wins.
