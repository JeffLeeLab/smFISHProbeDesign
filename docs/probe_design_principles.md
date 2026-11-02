# ProbeDesign — Probe Generation & Filtering Parameters

> **Purpose**: Document how the ProbeDesign Python pipeline generates, scores, filters, and selects oligonucleotide probes for smFISH (Sections 1–9) and HCR v3 split-initiator pairs (Section 10).
> **Code base**: `src/probedesign/` in `JeffLeeLab/smFISHProbeDesign`. Code references give the file (and function where useful). They leave out line numbers, which go stale.
> **Last updated**: 2026-11-02

---

## Table of Contents

1. [Pipeline Overview](#1-pipeline-overview)
2. [Candidate Generation (Sliding Window)](#2-candidate-generation-sliding-window)
3. [Thermodynamic Model](#3-thermodynamic-model)
4. [Thermodynamic Filtering (F Mask)](#4-thermodynamic-filtering-f-mask)
4B. [Low-Complexity Filtering (L Mask)](#4b-low-complexity-filtering-l-mask)
5. [What Is NOT Checked](#5-what-is-not-checked)
6. [Specificity Filtering (Masking)](#6-specificity-filtering-masking)
   - [6A. Repeat Masking (R Mask)](#6a-repeat-masking-r-mask)
   - [6B. Pseudogene Masking (P Mask)](#6b-pseudogene-masking-p-mask)
   - [6C. Genome Masking (B Mask)](#6c-genome-masking-b-mask)
   - [6D. Mask-to-Badness Conversion](#6d-mask-to-badness-conversion)
7. [Dynamic Programming Optimization](#7-dynamic-programming-optimization)
8. [Output: Probe Extraction & Reporting](#8-output-probe-extraction--reporting)
9. [Parameter Reference Table](#9-parameter-reference-table)
10. [HCR Split-Initiator Probes](#10-hcr-split-initiator-probes)

---

## 1. Pipeline Overview

```
Input FASTA
    |
    v
[1] Parse & Concatenate multi-entry FASTA (exon junctions marked with '>')
    |                                                       Code: fasta.py
    v
[2] Clean Sequence (lowercase; keep only a,c,t,g,n,>)
    |                                                       Code: sequence.py
    v
[3] Calculate Badness (thermodynamic scoring per position)
    |   badness[i] = (dG[i] - target)^2  or  inf           Code: core.py
    v
[4] Create F Mask (positions where badness == inf, BEFORE other masks)
    |                                                       Code: core.py
    v
[4b] Create L Mask (low-complexity: homopolymers & dinucleotide repeats)
    |                                                       Code: core.py
    v
[5] Apply Optional Masks (in order):
    |   (a) R mask — repeat regions (N's in input or separate file or RepeatMasker)
    |   (b) P mask — pseudogene alignments (bowtie, 16-mer)
    |   (c) B mask — genome alignments (bowtie, 12/14/16-mer triple-stringency)
    |                                                       Code: masking.py
    v
[6] Merge masks into badness (any masked position => inf)
    |                                                       Code: core.py
    v
[7] Dynamic Programming — find optimal probe positions
    |   minimizing average badness                          Code: core.py
    v
[8] Extract Probes — reverse complement template to get probe sequences
    |                                                       Code: core.py
    v
[9] Write Output — _oligos.txt and _seq.txt
                                                            Code: output.py
```

---

## 2. Candidate Generation (Sliding Window)

Candidates are **implicit** — every position `i` in `[0, len(seq) - oligo_length]` is a potential probe start. There is no explicit list of candidates; instead, the algorithm computes a "badness" score per position and the dynamic programming step selects which positions to use.

### How It Works

For a sequence of length `L` and oligo length `k` (default 20 bp):
- There are `L - k + 1` candidate start positions (indices `0` through `L - k`).
- At each position `i`, the candidate oligo is `seq[i : i + k]`.
- The probe sequence (what you order) is the **reverse complement** of this oligo.

### End Effects

| Position | Behavior |
|----------|----------|
| `i < 0` | Not considered (array bounds) |
| `0 <= i <= L - k` | Evaluated normally; badness is calculated |
| `i > L - k` | No candidate possible (oligo would extend past sequence end) |

There is no explicit padding or special handling at sequence ends — the sliding window simply stops when there aren't enough nucleotides remaining.

### Junction Handling (Multi-Exon FASTA)

When the input FASTA contains multiple entries (e.g., separate exons), they are concatenated with `>` characters marking junctions:

```
exon1_seq > exon2_seq > exon3_seq
```

The `>` character is treated as an invalid character. Any oligo window that contains `>` receives `badness = inf`, preventing probes from spanning exon junctions.

When extracting the actual probe sequence from a selected position, `>` characters are skipped:
```python
# core.py
while len(template_region) < oligo_length and j < len(seq):
    if seq[j] != '>':
        template_region += seq[j]
    j += 1
```

**Code reference:** `src/probedesign/fasta.py` — `sequences_to_single_string()` performs concatenation.

### Overlap Handling

There is **no hard prohibition on overlapping probes**. Instead, a minimum spacing constraint is enforced by the DP algorithm:

```
minimum distance between probe starts = oligo_length + spacer_length
```

With defaults (oligo_length=20, spacer_length=3), probe starts must be at least 23 bp apart. Since probes are 20 bp long, this creates a minimum 3 bp gap between consecutive probes. Mixed-length mode enforces the same gap between one probe's end and the next probe's start (Section 7).

**Code reference:** `src/probedesign/core.py` — `probe_spacer_len = oligo_len + spacer_len`

---

## 3. Thermodynamic Model

### Reference

Sugimoto et al. (1995), "Thermodynamic Parameters To Predict Stability of RNA/DNA Hybrid Duplexes", *Biochemistry*.

### Nearest-Neighbor Parameters

The model uses 16 dinucleotide stacking parameters for RNA-DNA hybrids:

**Enthalpy ΔH (kcal/mol)** — `thermodynamics.py`:

| | A | C | G | T |
|---|---|---|---|---|
| **A** | -7.8 | -5.9 | -9.1 | -8.3 |
| **C** | -9.0 | -9.3 | -16.3 | -7.0 |
| **G** | -5.5 | -8.0 | -12.8 | -7.8 |
| **T** | -7.8 | -8.6 | -10.4 | -11.5 |

**Entropy ΔS (cal/(mol·K))** — `thermodynamics.py`:

| | A | C | G | T |
|---|---|---|---|---|
| **A** | -21.9 | -12.3 | -23.5 | -23.9 |
| **C** | -26.1 | -23.2 | -47.1 | -19.7 |
| **G** | -13.5 | -17.1 | -31.9 | -21.6 |
| **T** | -23.2 | -22.9 | -28.4 | -36.4 |

**Initialization terms** — `thermodynamics.py`:
- ΔH_init = +1.9 kcal/mol
- ΔS_init = -3.9 cal/(mol·K)

### Gibbs Free Energy Calculation

```
For an oligo of length k:
  ΔH = Σ(ΔH[dinuc_i]) + ΔH_init    (sum over k-1 dinucleotides)
  ΔS = Σ(ΔS[dinuc_i]) + ΔS_init

  ΔG = (ΔH × 1000 - T × ΔS) / 1000    [kcal/mol]

where T = 310.15 K (37°C)
```

**Code reference:** `src/probedesign/thermodynamics.py` — `thermo_rna_dna(seq)`

### Melting Temperature Calculation

```
Tm = (ΔH × 1000) / (ΔS + R × ln(Ct/4)) - 273.15 + 16.6 × log₁₀(salt)

where:
  R    = 1.9872 cal/(mol·K)     (gas constant)
  Ct   = 50 × 10⁻⁶ M           (primer concentration, 50 µM)
  salt = 0.33 M                 (sodium concentration, 2× SSC buffer)
```

The Tm equation uses the SantaLucia (1998) PNAS form with the OligoCalc salt correction term.

**Code reference:** `src/probedesign/thermodynamics.py`

---

## 4. Thermodynamic Filtering (F Mask)

The F mask identifies positions where a probe **cannot start** due to thermodynamic reasons. This is purely based on the badness calculation, computed **before** any alignment-based masking.

### Badness Formula

For each position `i`:

```python
oligo = seq[i : i + oligo_length]

if has_invalid_chars(oligo):      # Contains N, X, M, H, P, B, or >
    badness[i] = inf
elif gibbs_rna_dna(oligo) raises KeyError:
    badness[i] = inf
elif ΔG < min_allowable or ΔG > max_allowable:
    badness[i] = inf
else:
    badness[i] = (ΔG - target_ΔG)²
```

**Key points:**
- `badness = 0` is the ideal case (ΔG exactly equals target)
- The badness function is a **squared loss** — probes far from the target ΔG are penalized quadratically
- `badness = inf` means the position is **completely excluded** (no probe can start here)
- The F mask in the output file shows `F` at positions with `inf` badness, and the sequence character at positions with finite badness
- In mixed-length mode each `(start, length)` cell is checked separately. The F line shows `F` only where **no** length in the range is acceptable, so an unmarked start may still have out-of-range lengths; the DP only ever picks a length whose own ΔG is in range

**Code reference:** `src/probedesign/core.py` — `calculate_badness()`

### Invalid Character Detection

```python
# sequence.py — has_invalid_chars()
invalid_pattern = re.compile(r'[xXmMhHpPbBnN>]')
```

Characters that trigger `inf` badness: `x`, `X`, `m`, `M`, `h`, `H`, `p`, `P`, `b`, `B`, `n`, `N`, `>`

**Code reference:** `src/probedesign/sequence.py`

---

### 4B. Low-Complexity Filtering (L Mask)

The L mask identifies candidate probes containing low-complexity subsequences — specifically homopolymer runs and dinucleotide repeats. Unlike the nucleotide-level R/P/B masks, the L mask operates at the **probe level**: each candidate probe substring is evaluated independently.

The L mask is **always active** (no CLI flag to enable/disable it). It modifies the badness array directly (setting matching entries to `inf`), using the same architecture as the F mask. It does **not** write to the `full_mask[]` array used by R/P/B masks.

#### Checks Performed

| Check | What it detects | Default threshold | CLI flag | Range |
|-|-|-|-|-|
| Homopolymer run | Consecutive identical bases (e.g. `AAAAA`) | 5 | `--hp-threshold` | 3-10 |
| Dinucleotide repeat | Repeating two-base motif (e.g. `ACACAC`) | 3 | `--di-threshold` | 3-10 |

A candidate probe is rejected (`badness = inf`) if it contains a homopolymer run of length >= `hp_threshold` **or** a dinucleotide repeat of >= `di_threshold` units.

#### Fixed-Length Mode

For a single oligo length `k`, the algorithm checks each probe `seq[i : i + k]`. If the probe fails either low-complexity check, `badness[i] = inf`.

#### Mixed-Length Mode

In mixed-length mode, each `(position, length)` pair in the 2D badness table is checked independently. A shorter probe that does not reach a low-complexity region is preserved even if a longer probe at the same start position is rejected. This mirrors the per-cell independence of the F mask in mixed mode.

#### Examples

With default thresholds (`hp_threshold=5`, `di_threshold=3`):

```
ATCGAAAAATCGATCG    → rejected (homopolymer: AAAAA, run of 5)
ATCGACACACGATCG     → rejected (dinucleotide: ACACAC, 3 repeats of AC)
ATCGAAAATCGATCGA    → accepted (AAAA is only 4, below threshold of 5)
```

---

## 5. What Is NOT Checked

The following criteria are **not used** in probe selection or filtering:

| Criterion | Status | Notes |
|-----------|--------|-------|
| **GC content (%)** | Reported, not filtered | GC% is computed and printed in output but does NOT affect the badness score or filtering |
| **Melting temperature (Tm)** | Reported, not filtered | Tm is computed and printed but does NOT affect scoring; ΔG filtering acts as an indirect Tm filter |
| **Self-dimers** | Not checked | No analysis of probe self-complementarity |
| **Hairpins / secondary structure** | Not checked | No folding prediction (e.g., no mfold/UNAfold integration) |
| **Cross-hybridization between probes** | Not checked | Probes are designed independently; no pairwise dimer check |
| **Off-target binding (transcriptome BLAST)** | Not checked directly | Only genome-level repeat filtering (bowtie) and pseudogene filtering are performed; there is no BLAST against the full transcriptome |
| **Probe GC clamp** | Not checked | No requirement for G/C at 3' end |

> **Implication**: The pipeline optimizes purely for thermodynamic uniformity (ΔG close to target) + specificity (avoid pseudogenes, repeats, and repetitive genome regions). Users who need secondary structure or cross-reactivity analysis should run a separate tool (e.g., OligoAnalyzer, NUPACK) on the output probes.

---

## 6. Specificity Filtering (Masking)

All masking is optional and controlled via CLI flags. Masked positions receive `badness = inf`, completely excluding them from probe selection.

### 6A. Repeat Masking (R Mask)

Identifies repetitive elements (SINEs, LINEs, etc.) in the target sequence.

**Three modes:**

| Mode | CLI Flag | Input | How It Works |
|------|----------|-------|-------------|
| Auto-detect | (none — automatic) | N's already in input FASTA | Positions with `n`/`N` are masked |
| Manual file | `--repeatmask-file PATH` | Separate FASTA with N's at repeat positions | Positions where masked file has `N` get masked |
| Automatic | `--repeatmask` | RepeatMasker runs on input | RepeatMasker produces `.masked` file; N positions become masked |

The `--repeatmask` and `--repeatmask-file` flags are **mutually exclusive**.

**Visualization:** In the `_seq.txt` output, the R mask line shows `R` at masked positions and the original sequence character at unmasked positions.

**Code reference:** `src/probedesign/core.py` (R mask creation), `src/probedesign/masking.py` (RepeatMasker integration)

### 6B. Pseudogene Masking (P Mask)

Identifies regions that align to known pseudogene sequences, which could cause cross-hybridization.

**CLI flag:** `--pseudogene-mask` (requires `--index-dir` or default `bowtie_indexes/`)

**Algorithm:**

1. Generate all 16-mers from the input sequence
2. Align each 16-mer against the pseudogene bowtie index using:
   ```
   bowtie -f -v 0 -k 1 --quiet <pseudogene_index> -
   ```
   - `-v 0`: zero mismatches allowed
   - `-k 1`: report up to 1 alignment (enough to detect any hit)
3. Count hits per position (0 or 1+)
4. Apply threshold: **any hit** (threshold = 0) triggers masking
5. Extend mask: each hit position extends the mask across the full 16-mer window
6. Remove short runs: keep only masked runs of **≥ 18 bp** (min_length=20, tolerance=2)

**Pseudogene index names:**

| Species | Index Prefix | FASTA Source |
|---------|-------------|--------------|
| human | `humanPseudo` | `data/pseudogeneDBs/human.fasta` |
| mouse | `mousePseudo` | `data/pseudogeneDBs/mouse.fasta` |
| drosophila | `drosophilaPseudo` | `data/pseudogeneDBs/drosophila.fasta` |
| elegans | `celegansPseudo` | `data/pseudogeneDBs/elegans.fasta` |
| rat | `ratPseudo` | No FASTA shipped (partial support) |

**Code reference:** `src/probedesign/masking.py` — `pseudogene_mask()`

### 6C. Genome Masking (B Mask)

Identifies regions that align to many locations in the genome (repetitive sequences that would cause non-specific hybridization).

**CLI flag:** `--genome-mask` (requires `--index-dir` or default `bowtie_indexes/`)

**Algorithm (triple-stringency):**

The pipeline runs **three separate bowtie searches** at different k-mer lengths and reports back positions exceeding a hit-count threshold:

| Pass | k-mer Length | Max Hits (`-k`) | Threshold | Rationale |
|------|-------------|-----------------|-----------|-----------|
| 1 | 12 bp | 5,000 | 4,000 | Catches highly repetitive short motifs |
| 2 | 14 bp | 1,000 | 500 | Catches moderately repetitive regions |
| 3 | 16 bp | 100 | 20 | Catches regions with moderate copy number |

Each pass:
1. Generate all N-mers from the input sequence
2. Align against the genome bowtie index:
   ```
   bowtie -f -v 0 -k <max_hits> --quiet <genome_index> -
   ```
3. Count total hits per position: `num_other_alignments + 1`
4. Positions with hits > threshold get masked (extended across the mer window)

The three masks are combined with **OR logic**: if any pass masks a position, it is masked.

**Genome index names:**

| Species | Index Prefix | Source |
|---------|-------------|--------|
| human | `GCA_000001405.15_GRCh38_no_alt_analysis_set` | Pre-built from JHU FTP |
| mouse | `mm10` | Pre-built from JHU FTP |
| drosophila | `drosophila` | Build from dm6 or rename pre-built |
| elegans | `celegans` | Build or download |
| rat | `rat` | Build or download |
| cow | `cow` | Build or download |

**Code reference:** `src/probedesign/masking.py` — `genome_mask()`

### 6D. Mask-to-Badness Conversion

After all masks are computed and combined (OR-union), the combined mask is converted to badness scores:

```python
# masking.py — mask_to_badness()
for each oligo start position i:
    if ANY nucleotide in seq[i : i + oligo_length] is masked:
        badness[i] = inf
```

This means a single masked nucleotide blocks all oligos that would overlap it.

The mask-derived badness is **added** to the thermodynamic badness:
```python
# core.py
for i in range(len(badness)):
    if mask_badness[i] == inf:
        badness[i] = inf
```

**Important ordering:** The F mask (in the output visualization) is computed from the **pre-masking** badness. This matches the MATLAB behavior where `mask_string(inseq, badness==inf, 'F')` is called before alignment masking is applied to the badness array.

---

## 7. Dynamic Programming Optimization

### Objective

Find the probe configuration (set of start positions) that **minimizes the average badness** across all selected probes.

### DP Formulation

```
State: bmsf[x][k]   (best solution ending at or before position x, using k+1 probes)
  - bmsf_sco[x][k] = best average score
  - bmsf_pos[x][k] = position of the (k+1)-th probe

Spacing constraint:
  probe_spacer_len = oligo_length + spacer_length
  (consecutive probe starts must be >= probe_spacer_len apart)

Initialization:
  bmsf_sco[0][0] = badness[0]
  bmsf_pos[0][0] = 0

Recurrence (for x = 1..goodlen-1, k = 0..n_probes-1):
  1. Copy forward: bmsf[x][k] = bmsf[x-1][k]

  2. Try placing probe k at position x:
     if k == 0:
       potential_score = badness[x]         # First probe: just its badness
     else:
       prev_x = x - probe_spacer_len
       if prev_x >= 0 and bmsf[prev_x][k-1] exists:
         potential_score = (bmsf_sco[prev_x][k-1] * k + badness[x]) / (k + 1)

  3. if potential_score < bmsf_sco[x][k]:
       update bmsf[x][k] = (position=x, score=potential_score)
```

### Scoring Function

The score is a **running average**:

```
score(k+1 probes) = (score(k probes) × k + badness[new_probe]) / (k + 1)
```

This means the optimization seeks to minimize the **mean** badness, not the sum. All probes contribute equally to the final score.

**Code reference:** `src/probedesign/core.py` — `_calc_score()`

### Backtracking

After filling the DP table, the algorithm backtracks from position `goodlen - 1` for each probe count `k`:
1. Record the last probe position
2. Jump back by `probe_spacer_len`
3. Repeat for `k-1`, `k-2`, ..., `0`
4. Reverse the collected positions

**Code reference:** `src/probedesign/core.py`

### Solution Selection

The DP returns solutions for 1, 2, ..., N probes. The pipeline selects the **solution with the most probes** (last in the list), provided its score is below the rejection threshold (`score < 1,000,000`).

If fewer than N probes can be placed (e.g., too many masked positions), the solution with the most feasible probes is used.

### Mixed-Length DP (`-l 18-22`)

With a length range, badness is a 2D table `badness[pos][len]` (`calculate_badness_mixed()`). Each cell is filtered independently by the F and L masks, and the P/B/R masks are applied per cell by `mask_to_badness_mixed()`. Because probes can differ in length, the DP is indexed by **end** position (`find_best_probes_mixed()`):

```
State: dp[e][k] = best average badness for k+1 probes, last probe ending at or before e

For each end position e:
  1. Copy forward: dp[e][k] = dp[e-1][k]
  2. For each length L in [min_len, max_len], start x = e - L + 1:
       k == 0: score = badness[x][L]
       k  > 0: prev_end = x - spacer_length - 1
               score = (dp[prev_end][k-1] * k + badness[x][L]) / (k + 1)
       keep the minimum in dp[e][k], recording (x, L)
```

Every length is considered at every position, so the optimum covers placements and lengths jointly. GC-rich regions tend to get shorter probes and AT-rich regions longer ones, all against the same `target_gibbs` and allowable range. A single integer `-l` uses the fixed-length DP above, unchanged.

**Code reference:** `src/probedesign/core.py` — `solutions[-1]`

---

## 8. Output: Probe Extraction & Reporting

### Probe Sequence

For each selected position `i`:
1. Extract the template region `seq[i : i + oligo_length]`, skipping `>` junction markers
2. Compute the **reverse complement** — this is the probe sequence (what you order from the oligo vendor)

**Code reference:** `src/probedesign/core.py`

### Output Files

**`_oligos.txt`** — Tab-separated, one probe per line:

| Column | Content | Example |
|--------|---------|---------|
| 1 | Index (1-based) | `1` |
| 2 | Start position (1-based, nucleotides) | `11` |
| 3 | GC% (integer) | `55` |
| 4 | Tm (°C, 1 decimal) | `66.2` |
| 5 | ΔG (kcal/mol, 1 decimal) | `-24.7` |
| 6 | Probe sequence | `gtctcagaagctgcgattcg` |
| 7 | Probe name | `KRT19_withUTRs_1` |

The reference files in `tests/regression/` lack the Start column.

**`_seq.txt`** — Visual alignment wrapped at 110 characters:

| Line | Content |
|------|---------|
| 1 | Original sequence (with `>` at junctions) |
| 2 | R mask (if repeat masking used): `R` = masked, sequence char = unmasked |
| 3 | P mask (if pseudogene masking used): `P` = masked |
| 4 | B mask (if genome masking used): `B` = masked |
| Next | F mask (always present): `F` = inf badness, sequence char = valid |
| Next | L mask (always present): `L` = low-complexity rejected, sequence char = valid |
| Below | Probe complementary sequences and labels: `Prb# N,Pos X,FE Y,GC Z` |

With masking on, `_bowtie_pseudogene.txt` / `_bowtie_genome.txt` hold per-position hit counts. `--save-bowtie-raw` adds the raw alignments.

**Code reference:** `src/probedesign/output.py`

---

## 9. Parameter Reference Table

### User-Configurable Parameters (CLI)

| Parameter | Default | CLI Flag(s) | Description | Code Location |
|-----------|---------|-------------|-------------|---------------|
| Number of probes | 48 | `-n`, `--probes` | Target number of probes to design | `cli.py`, `core.py` |
| Oligo length | 20 bp | `-l`, `--oligo-length` | Single length (e.g. `20`) or range for mixed-length mode (e.g. `18-22`, bounds 10–60) | `cli.py`, `core.py` |
| Spacer length | 3 bp | `-s`, `--spacer-length` | Minimum gap between probe starts (added to oligo length) | `cli.py`, `core.py` |
| Target ΔG | -23.0 kcal/mol | `-g`, `--target-gibbs` | Optimal Gibbs free energy for probe binding | `cli.py`, `core.py` |
| Allowable ΔG range | -26, -20 kcal/mol | `--allowable-gibbs "min,max"` | Probes outside this range get `inf` badness | `cli.py`, `core.py` |
| Output prefix | (input filename) | `-o`, `--output` | Base name for output files | `cli.py` |
| Quiet mode | off | `-q`, `--quiet` | Suppress progress messages | `cli.py` |
| Species | human | `--species` | Species for masking databases (human/mouse/elegans/drosophila/rat) | `cli.py` |
| Pseudogene mask | off | `--pseudogene-mask` | Enable pseudogene alignment masking | `cli.py` |
| Genome mask | off | `--genome-mask` | Enable repetitive genome region masking | `cli.py` |
| Index directory | `<repo>/bowtie_indexes/` | `--index-dir PATH` | Directory containing bowtie indexes | `cli.py`, `masking.py` (`DEFAULT_INDEX_DIR`) |
| Raw bowtie output | off | `--save-bowtie-raw` | Also write raw alignment files | `cli.py`, `output.py` |
| Repeat mask file | (none) | `--repeatmask-file PATH` | FASTA with N's marking repeat regions | `cli.py` |
| Auto repeat mask | off | `--repeatmask` | Run RepeatMasker automatically | `cli.py` |
| Homopolymer threshold | 5 | `--hp-threshold` | Reject probes with homopolymer runs >= this length (3-10) | `cli.py`, `core.py` |
| Dinucleotide repeat threshold | 3 | `--di-threshold` | Reject probes with dinucleotide repeats >= this many units (3-10) | `cli.py`, `core.py` |

### Hard-Coded Constants (Not Configurable via CLI)

| Constant | Value | Purpose | Code Location |
|----------|-------|---------|---------------|
| Salt concentration | 0.33 M (= 2× SSC) | Tm and ΔG calculation | `thermodynamics.py` |
| Primer concentration | 50 µM (50 × 10⁻⁶ M) | Tm calculation | `thermodynamics.py` |
| Temperature | 37°C (310.15 K) | ΔG calculation | `thermodynamics.py` |
| Gas constant R | 1.9872 cal/(mol·K) | Tm calculation | `thermodynamics.py` |
| ΔH_init | +1.9 kcal/mol | Nearest-neighbor initialization | `thermodynamics.py` |
| ΔS_init | -3.9 cal/(mol·K) | Nearest-neighbor initialization | `thermodynamics.py` |
| Pseudogene mer length | 16 bp | Bowtie alignment seed | `masking.py` |
| Pseudogene mismatches | 0 | Bowtie `-v` parameter | `masking.py` (default in `run_bowtie`) |
| Pseudogene hit threshold | 0 (any hit) | Minimum hits to trigger mask | `masking.py` |
| Pseudogene min run length | 20 bp | Minimum masked run to keep | `masking.py` |
| Pseudogene run tolerance | 2 bp | Gap allowance in run detection | `masking.py` |
| Genome 12-mer max hits | 5,000 | Bowtie `-k` parameter | `masking.py` |
| Genome 12-mer threshold | 4,000 | Hits above this → masked | `masking.py` |
| Genome 14-mer max hits | 1,000 | Bowtie `-k` parameter | `masking.py` |
| Genome 14-mer threshold | 500 | Hits above this → masked | `masking.py` |
| Genome 16-mer max hits | 100 | Bowtie `-k` parameter | `masking.py` |
| Genome 16-mer threshold | 20 | Hits above this → masked | `masking.py` |
| Score rejection threshold | 1,000,000 | DP solutions above this are dropped | `core.py` |
| Output line width | 110 chars | Sequence alignment wrapping | `output.py` |
| Invalid character regex | `[xXmMhHpPbBnN>]` | Characters that cause `inf` badness | `sequence.py` |

---

## 10. HCR Split-Initiator Probes

HCR v3 probes have their own pipeline: `probedesign hcrdesign` → `design_hcr_probes()` in `hcr.py`. Don't use `design -l 52` for HCR. That treats the pair as one 52-nt duplex, but the two halves hybridise independently.

### Pair Geometry

```
target:  5' ─[ left 25 nt ]─ gap 2 nt ─[ right 25 nt ]─ 3'     (52-nt block, start p)
           left = p..p+24,  gap = p+25..p+26,  right = p+27..p+51

P1 (odd, binds right half):  5'─ RC(right) ─ spacer_b ─ initiator_b ─3'
P2 (even, binds left half):  5'─ initiator_a ─ spacer_a ─ RC(left) ─3'
```

Both initiator tails point inward toward the gap. The initiator and spacer sequences for the 13 amplifiers (B1–B5, B7, B9–B11, B13–B15, B17) are in `AMPLIFIER_TABLE`. Amplifiers B7 and up use a degenerate `WW` spacer, which `--resolve-spacer AA` (for example) replaces with explicit bases. In the output, binding regions are lowercase and spacer + initiator are uppercase (`attach_initiators()`).

### Pair Badness (`calculate_pair_badness()`)

1. Compute 25-nt badness with the smFISH formula (Section 4) against the **strict** range (default target −31, range −35 to −27 kcal/mol per half). Apply the L mask to each half.
2. A pair at `p` is invalid (`inf`) if either half overlaps the R mask, or the 2-nt gap contains a `>` junction.
3. ΔG score = `badness(left) + badness(right)`. With `--asymmetric-gibbs`, one half may instead use the **lenient** range (`--lenient-gibbs-min`, default −42, up to the strict max). The better of left-strict and right-strict is kept, and the orientation is recorded.
4. P/B mask: by default both halves must be free of pseudogene/genome hits. `--asymmetric-bowtie` requires only one clean half.

### Pair DP (`find_best_pairs()`)

This is the same running-average DP as Section 7, applied to 52-nt blocks. Consecutive pair starts must be at least `52 + pair_spacing` apart (default spacing 3). Among the solutions it **maximises pair count first**, then minimises average badness.

### HCR Parameters

| Parameter | Default | CLI Flag |
|-----------|---------|----------|
| Pairs | 30 | `-p`, `--pairs` |
| Amplifier | B1 | `-a`, `--amplifier` |
| Target ΔG per half | −31.0 | `-g`, `--target-gibbs` |
| Strict ΔG range | −35, −27 | `--allowable-gibbs` |
| Asymmetric ΔG / lenient floor | off / −42.0 | `--asymmetric-gibbs`, `--lenient-gibbs-min` |
| Asymmetric bowtie | off | `--asymmetric-bowtie` |
| Pair spacing | 3 | `--pair-spacing` |
| Masking, L-mask thresholds | as smFISH | `--species`, `--pseudogene-mask`, `--genome-mask`, `--repeatmask-file`, `--hp-threshold`, `--di-threshold` |
| WW spacer resolution | none (IUPAC) | `--resolve-spacer` |

Thermodynamic constants (salt, primer concentration, temperature) match smFISH.

### Outputs (`hcr_output.py`)

- `_HCR_oligos.txt`: `pair, half, start, GC%, Tm, Gibbs, strict, full_oligo, name`
- `_HCR_seq.txt`: sequence map with mask lines and pair annotations
- `_HCR_hits.txt`: per-half bowtie hit counts
- `_HCR_bowtie_*_raw.txt`: only with `--save-bowtie-raw`

**Tests:** `tests/test_hcr.py`. **Design background:** [hcr_designer_plan.md](hcr_designer_plan.md).

`tests/regression/EIF1_CDS_HCR/` predates `hcrdesign`. It holds an old single-52-mer design and is not a reference for this pipeline.
