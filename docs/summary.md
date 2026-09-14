# ProbeDesign — Setup & Development Summary

> Consolidated documentation covering installation, key findings, and the Streamlit web app.
>
> **Last updated**: 2026-02-18

---

## Table of Contents

1. [Installation Steps](#1-installation-steps)
2. [Key Findings During Installation](#2-key-findings-during-installation)
3. [Streamlit Web App](#3-streamlit-web-app)
4. [Python Source Changes for the App](#4-python-source-changes-for-the-app)
5. [Validation & Test Results](#5-validation--test-results)
6. [Quick Reference](#6-quick-reference)

---

## 1. Installation Steps

All steps below are automated in [`setup_all.sh`](setup_all.sh). This section explains what each step does.

### 1.1 Mamba Environment

We use **mamba** (a lightweight conda alternative) to create an isolated environment with pinned versions of Python, Bowtie 1, and the CLI framework.

```bash
mamba create -n probedesign \
    -c conda-forge -c bioconda \
    --channel-priority strict \
    python=3.11 \
    "click>=8.0" \
    bowtie=1.3.1 \
    "pytest>=7.0" \
    -y
```

| Package | Version | Channel | Purpose |
|---------|---------|---------|---------|
| python | 3.11.x | conda-forge | Runtime |
| click | ≥ 8.0 | conda-forge | CLI framework used by `probedesign` |
| bowtie | 1.3.1 | bioconda | Short-read aligner for masking |
| pytest | ≥ 7.0 | conda-forge | Test runner |

**Why mamba?** It resolves dependencies faster than conda and works well on Apple Silicon Macs.

**Why `--channel-priority strict`?** Prevents cross-channel conflicts between `conda-forge` and `bioconda`. Without this, the solver can pick incompatible builds.

**Why bioconda for bowtie?** On macOS, `brew install bowtie` installs a *different* tool (a file manager). The bioinformatics Bowtie 1 aligner must come from bioconda.

### 1.2 Create Environment and Install Packages

The preferred way is via the environment file (handles all dependencies in one step):

```bash
# Using mamba / conda (auto-detected by setup_all.sh):
mamba env create -f environment.yml
mamba activate probedesign
pip install -e .
```

Or to create the environment manually:

```bash
mamba create -n probedesign \
    -c conda-forge -c bioconda --channel-priority strict \
    python=3.11 "click>=8.0" bowtie=1.3.1 "pytest>=7.0" -y
mamba activate probedesign
pip install -e .
pip install "streamlit>=1.32" "pandas>=1.5"
```

The `-e` (editable) flag means code changes in `src/probedesign/` take effect immediately without reinstalling. The entry point is defined in `pyproject.toml`:

```toml
[project.scripts]
probedesign = "probedesign.cli:main"
```

Streamlit and pandas are required for the web GUI.

### 1.3 Build Pseudogene Indices

The repository ships pseudogene FASTA files in `probedesign/pseudogeneDBs/`. We build Bowtie 1 indices (`.ebwt`) from them:

```bash
mkdir -p bowtie_indexes

bowtie-build probedesign/pseudogeneDBs/human.fasta       bowtie_indexes/humanPseudo
bowtie-build probedesign/pseudogeneDBs/mouse.fasta       bowtie_indexes/mousePseudo
bowtie-build probedesign/pseudogeneDBs/drosophila.fasta  bowtie_indexes/drosophilaPseudo
```

| Species | Index prefix | Sequences | Size | Build time |
|---------|-------------|-----------|------|------------|
| Human | `humanPseudo` | 18,046 | 35 MB | ~23 s |
| Mouse | `mousePseudo` | 19,086 | 39 MB | ~22 s |
| Drosophila | `drosophilaPseudo` | 2,204 | 10 MB | ~2 s |

Each index produces 6 `.ebwt` files (84 MB total).

### 1.4 Download Genome Indices

Instead of building genome indices from FASTA (1–3 hours), we download **pre-built Bowtie 2** (`.bt2`) indices from the AWS index mirror maintained by Ben Langmead's lab. Bowtie 1 (v1.3+) reads `.bt2` files directly (see [Section 2.1](#21-bowtie-1--bowtie-2-index-compatibility)).

```bash
cd bowtie_indexes

# Human GRCh38 (~3.5 GB)
curl -L -o GRCh38_noalt_as.zip \
    https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip
# ZIP extracts into a GRCh38_noalt_as/ subdirectory — move .bt2 files up and
# rename the prefix to match the code expectation:
mv GRCh38_noalt_as/*.bt2 .
rm -rf GRCh38_noalt_as GRCh38_noalt_as.zip
# Files now named: GCA_000001405.15_GRCh38_no_alt_analysis_set.{1,2,3,4,rev.1,rev.2}.bt2

# Mouse mm10 (~3.1 GB)
curl -L -o mm10.zip https://genome-idx.s3.amazonaws.com/bt/mm10.zip
unzip mm10.zip && rm mm10.zip
# Files: mm10.{1,2,3,4,rev.1,rev.2}.bt2

# Drosophila BDGP6 (~176 MB)
curl -L -o BDGP6.zip https://genome-idx.s3.amazonaws.com/bt/BDGP6.zip
unzip BDGP6.zip
# Rename prefix from BDGP6 → drosophila to match masking.py expectations:
for f in BDGP6.*.bt2; do mv "$f" "drosophila.${f#BDGP6.}"; done
rm -f BDGP6.zip
```

| Species | Target prefix | Source | Download time | Size |
|---------|--------------|--------|---------------|------|
| Human GRCh38 | `GCA_000001405.15_GRCh38_no_alt_analysis_set` | AWS pre-built Bowtie 2 | ~5 min | 3.5 GB |
| Mouse mm10 | `mm10` | AWS pre-built Bowtie 2 | ~5 min | 3.1 GB |
| Drosophila BDGP6 | `drosophila` | AWS pre-built Bowtie 2 | ~20 s | 176 MB |

**Renaming**: The AWS zips sometimes use different prefixes (e.g., `GRCh38_noalt_as`, `BDGP6`) than what `masking.py` expects. The script renames the extracted files to match the hard-coded index names in the code (see `src/probedesign/masking.py`, the `GENOME_DB` dict).

### 1.5 (Optional) RepeatMasker

RepeatMasker is only needed for the `--repeatmask` automatic mode. If you always supply a pre-masked FASTA file via `--repeatmask-file`, you can skip this.

```bash
mamba install -n probedesign -c bioconda -c conda-forge repeatmasker -y
```

You then need the Dfam database partition for your species:

| Partition | Taxonomic group | Species covered | Size (extracted) |
|-----------|----------------|-----------------|------------------|
| 0 (root) | All | Always included | ~77 MB |
| 1 | Brachycera | *Drosophila melanogaster* | ~45 GB |
| 7 | Mammalia | Human, Mouse, Rat | ~56 GB |

Download partitions into the famdb directory:

```bash
FAMDB="$MAMBA_ROOT_PREFIX/share/RepeatMasker/Libraries/famdb"
cd "$FAMDB"
curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.7.h5.gz
gunzip dfam39_full.7.h5.gz    # produces ~56 GB file; ensure ≥65 GB free
```

---

## 2. Key Findings During Installation

### 2.1 Bowtie 1 ↔ Bowtie 2 Index Compatibility

**Discovery**: Bowtie 1 (version 1.3+) can **read Bowtie 2 `.bt2` index files directly**, with no conversion needed.

**How we found it**: The original documentation and code assumed genome indices must be in `.ebwt` (Bowtie 1) format. Building `.ebwt` from human/mouse genome FASTAs takes 1–3 hours and requires significant RAM. While investigating faster alternatives, we tested whether the widely available pre-built Bowtie 2 indices (`.bt2`) from [AWS genome indexes](https://benlangmead.github.io/aws-indexes/bowtie) would work — and they did.

**Practical impact**:
- **Before**: Download genome FASTA → `bowtie-build` (1–3 hours, ~16 GB RAM for human)
- **After**: Download pre-built `.bt2` (12 minutes total for all 3 species, no build step)

**What we tested**: All 3 genome indices (human GRCh38, mouse mm10, drosophila BDGP6) were queried with Bowtie 1 using `.bt2` files. Bowtie reported hits correctly in every test. The full validation pipeline (`run_tests.sh`) also passes with `.bt2` genome indices.

### 2.2 Homebrew `bowtie` ≠ Bioinformatics Bowtie

On macOS, `brew install bowtie` installs a **file manager GUI app**, not the bioinformatics short-read aligner. The correct installation is via bioconda:

```bash
# WRONG (installs a file manager)
brew install bowtie

# CORRECT (installs the bioinformatics aligner)
mamba install -c bioconda bowtie=1.3.1
```

Verify you have the right one:

```bash
bowtie --version
# Correct output: "bowtie-align-s version 1.3.1"
```

### 2.3 Index Prefix Naming Matters

The ProbeDesign codebase has **hard-coded index prefix names** in `src/probedesign/masking.py`. When downloading pre-built indices from external sources, files must be renamed to match these expectations:

| Species | Pseudogene prefix | Genome prefix |
|---------|------------------|---------------|
| human | `humanPseudo` | `GCA_000001405.15_GRCh38_no_alt_analysis_set` |
| mouse | `mousePseudo` | `mm10` |
| drosophila | `drosophilaPseudo` | `drosophila` |

The Drosophila genome index from AWS extracts as `BDGP6.*`, which must be renamed to `drosophila.*` to match the code.

### 2.4 Channel Priority is Critical for Conda

Without `--channel-priority strict`, conda/mamba may resolve the `bowtie` package from the wrong channel and install an incompatible or incorrect version. Always use:

```bash
mamba create ... -c conda-forge -c bioconda --channel-priority strict ...
```

---

## 3. Streamlit Web App

### 3.1 Overview

A Streamlit-based web GUI was built at `streamlit_app/` to make ProbeDesign accessible without command-line expertise. The app wraps the same `design_probes()` function used by the CLI — no subprocess calls, no code duplication.

**Launch**:
```bash
mamba activate probedesign
streamlit run streamlit_app/app.py
```

Opens at `http://localhost:8501`.

### 3.2 Design Criteria

The app was built with the following goals:

1. **Direct library import** — Call `design_probes()` from `core.py` as a function, not via subprocess. This avoids shell escaping issues and gives access to the full `ProbeDesignResult` dataclass.

2. **Feature parity with CLI** — All CLI options (`--probes`, `--oligo-length`, `--target-gibbs`, `--pseudogene-mask`, `--genome-mask`, `--repeatmask`, `--repeatmask-file`, `--species`, `--save-bowtie-raw`) are exposed in the sidebar.

3. **Two modes** — "Single sequence" for individual files and "Batch" for processing a directory of FASTAs with progress tracking.

4. **Pre-flight checks** — Before running, the app verifies bowtie is installed, RepeatMasker is available (if needed), and index files exist for the selected species. Missing prerequisites show as warnings, not crashes.

5. **In-browser results** — Probe table (sortable dataframe), sequence alignment viewer (monospaced text), and download buttons for all output files (`_oligos.txt`, `_seq.txt`, bowtie hit tables, raw bowtie output).

6. **No heavyweight dependencies on the core package** — `streamlit` and `pandas` are only needed for the app, not for the CLI or library. They are not listed in `pyproject.toml`.

### 3.3 Architecture

```
streamlit_app/
├── app.py                      # Streamlit UI: sidebar params, single/batch mode panels
├── utils.py                    # Backend helpers: FASTA validation, design runner, batch runner
├── setup_all.sh                # Automated setup script
├── summary.md                  # This file — setup & development documentation
├── probe_design_principles.md  # Detailed probe generation & filtering parameters
└── README.md                   # Comprehensive end-user guide
```

**`app.py`** — Defines the sidebar (all parameters), input methods (file upload or paste), and result display (metrics, dataframe, download buttons). Two top-level modes are toggled with `st.radio`.

**`utils.py`** — Provides:
- `validate_fasta_text()` — Validates pasted FASTA (headers, valid base characters)
- `save_fasta_to_temp()` / `save_uploaded_file_to_temp()` — Write input to temp files for `design_probes()`
- `check_prerequisites()` — Checks bowtie, RepeatMasker, and index file availability
- `run_design()` — Wraps `design_probes()` with stdout capture and exception handling
- `run_batch()` — Iterates FASTA files, tracks progress, writes outputs
- `format_batch_summary()` / `write_batch_summary()` — TSV summary for batch runs

### 3.4 User Interface Layout

| Section | Location | Controls |
|---------|----------|----------|
| Mode toggle | Sidebar top | Radio: "Single sequence" / "Batch" |
| Basic params | Sidebar | Number inputs for probes, oligo length, spacer, Gibbs FE |
| Species | Sidebar | Selectbox: human, mouse, elegans, drosophila, rat |
| Masking toggles | Sidebar | Checkboxes: Pseudogene mask (on by default), Genome mask (on by default) |
| Index directory | Sidebar | Text input (auto-filled to `bowtie_indexes/`) |
| RepeatMasker | Sidebar expander | Radio: None / Auto / Provide file; file uploader |
| Input | Main panel | Tabs: "Upload FASTA" / "Paste FASTA" |
| Results | Main panel | Metrics, sortable table, sequence viewer, download buttons |

---

## 4. Python Source Changes for the App

The Streamlit app required several additions to the Python source package (`src/probedesign/`) to support in-process usage (vs CLI-only). These changes are backward-compatible — the CLI behaves identically.

### 4.1 `output.py` — Added `format_*` Functions

The original `output.py` only had `write_oligos_file()` and `write_seq_file()`, which write directly to disk. The Streamlit app needs file contents as strings (for `st.download_button`). Three new functions were added:

| Function | Purpose |
|----------|---------|
| `format_oligos_content(result)` | Returns `_oligos.txt` content as a string |
| `format_seq_content(result, mask_seqs)` | Returns `_seq.txt` content as a string |
| `format_hits_content(result)` | Returns bowtie hit table contents as a `dict[str, str]` |

These produce **identical output** to the `write_*` functions but return strings instead of writing files.

### 4.2 `output.py` — `write_output_files()` Got `output_dir` Parameter

Previously, output files were always written to the current working directory. A new `output_dir` parameter was added so the app and batch mode can specify where files go without `os.chdir()`:

```python
def write_output_files(result, output_prefix, mask_seqs=None, output_dir=None):
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        prefix = os.path.join(output_dir, output_prefix)
    else:
        prefix = output_prefix
    ...
```

### 4.3 `core.py` — `ProbeDesignResult` Extended

The `ProbeDesignResult` dataclass was extended with fields so bowtie alignment data can be passed through to the app's download buttons:

| New field | Type | Purpose |
|-----------|------|---------|
| `mask_strings` | `List[str]` | Visualization strings (R/P/B/F masks) for `_seq.txt` |
| `bowtie_pseudogene_hits` | `List[int]` | Parsed hit counts per position |
| `bowtie_genome_hits` | `dict` | Parsed hit counts per mer length |
| `bowtie_pseudogene_raw` | `str` | Raw bowtie alignment output |
| `bowtie_genome_raw` | `str` | Raw bowtie alignment output |

### 4.4 `core.py` — `save_bowtie_raw` Parameter

A `save_bowtie_raw` flag was added to `design_probes()` to optionally preserve the full bowtie alignment text. By default it is `False` (no extra memory usage), but the app and CLI `--save-bowtie-raw` can enable it.

### 4.5 Summary of Source Changes

All changes are **additive** — no existing function signatures were broken, and no existing behavior was changed. The CLI path is unaffected. The new `format_*` functions and extended dataclass simply provide the app with programmatic access to what was previously only available via file I/O.

---

## 5. Validation & Test Results

### 5.1 Test Cases

| Test | Command | Expected | Result |
|------|---------|----------|--------|
| CDKN1A_32 | `probedesign design CDKN1A.fa --probes 32 --repeatmask-file CDKN1A_repeatmasked.fa` | 32/32 probes match MATLAB (100%) | **PASS** |
| KRT19 | `probedesign design KRT19_withUTRs.fa --probes 32 --pseudogene-mask --genome-mask` | 6/6 probes match MATLAB (100%) | **PASS** |
| EIF1 HCR | `probedesign design EIF1_Exons.fasta --probes 20 -l 52 --target-gibbs -60 --allowable-gibbs -80,-40 --pseudogene-mask --genome-mask` | ≥ 15/20 probes match (≥ 75%) | **PASS (15/20, 75%)** |

### 5.2 Running Tests

```bash
# Quick validation (automated, all tests):
./streamlit_app/setup_all.sh

# Or, if already set up:
./run_tests.sh
```

### 5.3 Resource Summary

| Component | Disk usage | Setup time |
|-----------|-----------|------------|
| Mamba environment | ~500 MB | ~5 min |
| Pseudogene indices (3 species) | 84 MB | ~47 s |
| Genome indices (3 species) | 6.7 GB | ~12 min |
| RepeatMasker + Dfam partition 7 | ~56 GB | ~30 min |
| **Total (without RepeatMasker)** | **~7.3 GB** | **~18 min** |
| **Total (with RepeatMasker)** | **~63 GB** | **~48 min** |

---

## 6. Quick Reference

### Create / Activate Environment
```bash
# Create from environment.yml (first time):
mamba env create -f environment.yml
pip install -e .

# Activate (every session):
mamba activate probedesign   # or: conda activate probedesign
```

### Design Probes (CLI)
```bash
probedesign design input.fa --probes 32
probedesign design input.fa --probes 32 --pseudogene-mask --genome-mask
probedesign design input.fa --probes 32 --repeatmask-file masked.fa
```

### Launch Streamlit App
```bash
streamlit run streamlit_app/app.py
```

### Verify Installation
```bash
probedesign --version          # 0.1.0
bowtie --version               # bowtie-align-s version 1.3.1
```

### File Locations

| What | Path |
|------|------|
| Python package | `src/probedesign/` |
| CLI entry point | `src/probedesign/cli.py` |
| Streamlit app | `streamlit_app/app.py` |
| Conda environment file | `environment.yml` |
| Bowtie indices | `bowtie_indexes/` |
| Pseudogene FASTAs | `probedesign/pseudogeneDBs/` |
| Test cases | `test_cases/` |
| Setup script | `setup_all.sh` |
