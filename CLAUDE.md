## Subagent usage instruction

Use haiku subagents for general codebase exploration and simple documentation updates. Use sonnet subagents for web research or ingesting PDF articles or papers. 

## Project Overview

ProbeDesign designs oligonucleotide probes for single-molecule RNA FISH (smFISH) and HCR v3 split-initiator FISH. The product is the **Python package** (`src/probedesign/`) and the **Streamlit web app** (`streamlit_app/`). The MATLAB code (in `legacy/`) is legacy — don't modify it or chase MATLAB parity.

## Environment

**Always use the `probedesign` mamba/conda env.** Never install Python packages into the base environment. In non-interactive shells, call the env's binaries directly, e.g. `~/miniforge3/envs/probedesign/bin/python`.

```bash
./setup_all.sh                          # env from environment.yml, pip install -e ., bowtie indices
streamlit run streamlit_app/app.py      # web app
probedesign design input.fa --probes 48             # smFISH
probedesign hcrdesign input.fa --pairs 30 -a B1     # HCR
python -m pytest tests/                 # HCR unit, smFISH regression, run-header tests
```

The app pages put `src/` on `sys.path` themselves; the `probedesign` CLI needs the package installed (`pip install -e .`).

## Layout

```
src/probedesign/            Python package
streamlit_app/              Web app (app.py = top nav only; pages/ = one file per tab)
tests/                      pytest suite; regression/ = reference inputs + outputs (command.txt = original MATLAB call)
data/isoform_explorer/      Isoform Explorer data (Ensembl v113 Parquet; Hsap, Mmus, Dmel)
data/pseudogeneDBs/         Pseudogene FASTAs that setup_all.sh builds indices from
docs/                       Reference docs + implemented design plans
bowtie_indexes/             Bowtie indices (gitignored; built/downloaded by setup_all.sh)
legacy/                     Archived upstream MATLAB / Python 2 code — ignore
```

## Python package (`src/probedesign/`)

| File | Purpose |
|-|-|
| `cli.py` | Click CLI: `design` (smFISH), `hcrdesign` (HCR pairs), `analyze` (ΔG profile) |
| `core.py` | smFISH pipeline `design_probes()`: badness, masks, fixed- and mixed-length DP |
| `hcr.py` | HCR pipeline `design_hcr_probes()`: pair badness, pair DP, initiator attachment |
| `masking.py` | Bowtie pseudogene/genome masks, RepeatMasker, low-complexity checks |
| `thermodynamics.py` | ΔG / Tm (Sugimoto 1995 RNA-DNA; 0.33 M salt, 50 µM primer) |
| `output.py`, `hcr_output.py` | smFISH / HCR output files |
| `provenance.py` | `run_info` → `#` run header and equivalent CLI command for `_seq.txt` / `_HCR_seq.txt` |
| `sequence.py`, `fasta.py` | Sequence utilities; FASTA I/O (multi-entry joined with `>` junction markers) |

## Streamlit app (`streamlit_app/`)

`app.py` only builds the top nav (smFISH, HCR, Isoform Explorer, TransQuant). Each page owns its own sidebar and widgets.

| File | Purpose |
|-|-|
| `pages/smfish.py` | smFISH design, single + batch |
| `pages/hcr.py` | HCR design, single + batch |
| `pages/isoform_explorer.py`, `isoform_data.py` | Isoform browser; exports N-collapsed constitutive FASTA. No probe-design code |
| `pages/transquant.py` | Embeds external `transquant_w.ui.render()`. Independent of the other tabs |
| `utils.py` | Backend: FASTA validation, temp files, `check_prerequisites()`, run/batch wrappers, batch summaries + ZIPs |

```
pages/smfish.py → run_design() / run_batch()         (utils.py) → design_probes()     (core.py)
pages/hcr.py    → run_hcr_design() / run_hcr_batch() (utils.py) → design_hcr_probes() (hcr.py)
```

TransQuant is pinned to a git tag in **three** places — bump all together: `environment.yml`, `setup_all.sh`, `pages/transquant.py` (install hint).

## smFISH algorithm

1. Badness per position = `(ΔG − target)²`; `inf` for invalid chars or ΔG outside the allowable range. Defaults: target −23, range −26 to −20 kcal/mol.
2. Masks, each a line in `_seq.txt`:
   - **R** repeats — N's in input, `--repeatmask-file`, or RepeatMasker
   - **P** pseudogene / **B** genome — bowtie hits (genome uses 12/14/16-mer tiers)
   - **F** thermodynamic — badness is `inf`
   - **L** low-complexity, probe-level — homopolymer ≥ `--hp-threshold` (5) or dinucleotide repeat ≥ `--di-threshold` (3)
3. DP picks placements minimising mean badness. Fixed length: start-position DP.
4. Mixed length (`-l 18-22`): 2D badness `[pos][len]` (`calculate_badness_mixed()`), end-position DP over all lengths (`find_best_probes_mixed()`). F shows only where no length is finite. A single int `-l 20` runs the fixed-length path unchanged.

## HCR

Each pair is two 25-nt halves with a 2-nt gap; the 52-nt block is the DP unit. The halves hybridise independently, so ΔG is per half: target −31, strict range −35 to −27. `--asymmetric-gibbs` lets one half go down to `--lenient-gibbs-min` (−42); `--asymmetric-bowtie` requires only one half to be bowtie-clean. There are 13 amplifiers (`AMPLIFIER_TABLE` in `hcr.py`). Always use `hcrdesign`, never `design -l 52`.

## Masking setup

- **Bowtie 1** from bioconda (not homebrew, not bowtie2). v1.3+ reads `.bt2` indices directly.
- `--index-dir` defaults to `<repo>/bowtie_indexes`.
- `setup_all.sh` builds indices for human, mouse and drosophila only:

| Species | Pseudogene | Genome |
|-|-|-|
| human | `humanPseudo` | `GCA_000001405.15_GRCh38_no_alt_analysis_set` |
| mouse | `mousePseudo` | `mm10` |
| drosophila | `drosophilaPseudo` | `drosophila` |

  The CLI and app also accept `elegans` (`celegansPseudo` / `celegans`) and `rat` (`ratPseudo` / `rat`). Masking for those needs indices you build yourself.
- RepeatMasker: bioconda; mammals need Dfam partition 7 (~56 GB). See `docs/REPEATMASKER.md`.

## Outputs

- `_seq.txt` and `_HCR_seq.txt` start with `#` lines from `provenance.format_header()`:
  - the parameters, input, version and git hash;
  - each bowtie mask's status, index and fingerprint;
  - the equivalent CLI command.

  Files are named by basename or upload name, never by temp or server path. A requested mask that failed shows `FAILED`, and the app also shows a warning for it.
- **smFISH** (`output.py`):
  - `_oligos.txt`: TSV of index, start, GC%, Tm, ΔG, sequence, name.
  - `_seq.txt`: sequence blocks with R/P/B/F/L lines and probe labels.
  - `_bowtie_pseudogene.txt` / `_bowtie_genome.txt`: written when masking is on; `--save-bowtie-raw` adds `_raw` versions.
- **HCR** (`hcr_output.py`): `_HCR_oligos.txt`, `_HCR_seq.txt` and `_HCR_hits.txt`. `--save-bowtie-raw` adds `_HCR_bowtie_*_raw.txt`.

## Testing

```bash
python -m pytest tests/
```

- `test_hcr.py`: HCR unit tests.
- `test_smfish_regression.py`: smFISH against the references in `tests/regression/`:
  - CDKN1A (repeat-mask file), also through the CLI;
  - the L mask;
  - mixed length;
  - KRT19 and mouseMITF with pseudogene and genome masks. These skip if bowtie or the indices are missing.
- `test_provenance.py`: the `#` run header. It covers:
  - header fields;
  - the overlay being unchanged;
  - the mixed-length overlay;
  - CLI command round trips and flag drift;
  - mask status, fingerprints and failures;
  - the header never showing a path.
- The references predate the L mask. Tests pass `--hp-threshold 10 --di-threshold 10` (the maximum) to effectively switch it off; do the same when comparing by hand.
- `mixed_length_test/` is an input only, with no reference. `EIF1_CDS_HCR/` uses the old single-52-mer approach and is not a reference for `hcrdesign`.
- If bowtie isn't on `PATH`, masking prints a warning, records `FAILED` in the header and carries on unmasked. The tests put the interpreter's `bin/` on `PATH`. When running the CLI with the env's python directly, add `~/miniforge3/envs/probedesign/bin` yourself.

## Docs

- `docs/probe_design_principles.md`: algorithm reference.
- `docs/BOWTIE.md`, `docs/REPEATMASKER.md`: installation.
- `docs/summary.md`: dev notes.
- The other `docs/` files are design plans for features that have already shipped. They are history, not spec.

## Adding features

- **smFISH option**:
  1. Add a Click option to `design` in `cli.py`.
  2. Thread it through `design_probes()` in `core.py`.
  3. Add a widget in `pages/smfish.py`.
  4. Pass it through `run_design()` and the `run_batch()` params dict in `utils.py`.
- **HCR option**: same chain, using `hcrdesign` → `design_hcr_probes()` (`hcr.py`) → `pages/hcr.py` → `run_hcr_design()` / `run_hcr_batch()`.
- For either kind, also:
  - add the option to the `run_info` params dict in the design function;
  - add it to `build_command()` and `format_header()` in `provenance.py`.

  `test_provenance.py` checks that the command's flags exist.
- Then update `README.md` and this file.
