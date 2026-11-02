# L-Mask Implementation Plan: Low-Complexity Sequence Filtering

> **Status**: Implemented (rev 2, `ec32e72`). The shipped homopolymer default is **5**, not the 4 used throughout this plan; the dinucleotide default is 3 as planned. Where this plan disagrees with the code, the code wins.
> **Date**: 2026-02-24

---

## 1. Objective

Add an **L-mask** (low-complexity filter) that marks candidate probe starting positions as `inf` badness when the **probe's own substring** contains:

1. **Homopolymer runs** — stretches of a single nucleotide repeated N or more times (e.g. `AAAA`, `CCCC`). Default threshold: 4 (mask runs of 4+).
2. **Dinucleotide repeats** — alternating two-nucleotide motifs repeated N or more times (e.g. `ATATAT`, `GCGCGC`). All 12 ordered pairs considered. Default threshold: 3 (mask 3+ repeating units, i.e. 6+ nt).

The L-mask is **always active** (not optional) and sits in the pipeline immediately after the F-mask and before R/P/B masks.

---

## 2. Pipeline Position

Current pipeline order in `core.py:design_probes()`:

```
[3] Calculate badness (thermodynamic scoring)
[4] Create F mask visualization (badness == inf)
[5] Apply R mask → P mask → B mask  (optional, into full_mask[])
[6] Merge full_mask into badness via mask_to_badness()
[7] DP optimization
```

New pipeline order:

```
[3] Calculate badness (thermodynamic scoring)
[4] Create F mask visualization (badness == inf)
[4b] L mask: check each candidate probe substring for low complexity,
     set badness = inf directly (per position, per length)          ← NEW
[5] Create L mask visualization string                              ← NEW
[6] Apply R mask → P mask → B mask  (optional, into full_mask[])
[7] Merge full_mask into badness via mask_to_badness()
[8] DP optimization
```

**Critical design decision**: The L-mask operates on **probe substrings** (like the F-mask), NOT on nucleotide positions (unlike R/P/B masks). It modifies the badness array directly and does **not** go through `full_mask[]`. See Section 4 for the rationale.

---

## 3. Two Masking Architectures in the Codebase

The existing codebase uses two distinct masking architectures. Understanding why each exists is essential for choosing the right one for L-mask.

### Architecture 1: Nucleotide-Level Masks (R / P / B)

R/P/B masks flag individual **nucleotide positions** as problematic:
- R mask: "this nucleotide is in a repeat element"
- P mask: "this nucleotide is in a region that aligns to a pseudogene"
- B mask: "this nucleotide is in a highly repetitive genomic region"

These are stored in `full_mask[]` (a per-nucleotide binary array), then converted to probe-level badness via `mask_to_badness()` (fixed-length) or `mask_to_badness_mixed()` (mixed-length).

**How `mask_to_badness_mixed()` works** (`masking.py:377-412`):

```python
for i in range(max_goodlen):
    for L_idx in range(n_lengths):
        L = min_len + L_idx
        if any(mask[i:i + L]):     # ← checks the specific probe window
            badness[i][L_idx] = inf
```

This already provides **per-(position, length) sensitivity**: a probe is only killed if its window `seq[i:i+L]` overlaps a masked nucleotide. Shorter probes that don't reach the masked position survive. This is correct for R/P/B because the semantics are: "any probe touching this nucleotide inherits its problem." One problematic nucleotide anywhere in the probe window is sufficient reason to exclude it.

### Architecture 2: Probe-Level Checks (F mask)

The F mask evaluates each **candidate probe substring** independently. In `calculate_badness_mixed()` (`core.py:104-168`), each `(position, length)` pair gets its own Gibbs FE calculation on `seq[i:i+L]`. The badness is set directly — no intermediate nucleotide-level mask.

This is necessary because thermodynamic properties are properties of the **entire probe substring**, not of individual nucleotides. A nucleotide doesn't have an intrinsic "bad Gibbs FE" — only a probe of a specific length starting at a specific position does.

### Why L-mask Must Use Architecture 2 (Probe-Level)

The L-mask question is: "does this probe's own substring contain a homopolymer run ≥ threshold?" This is a **substring property**, not a nucleotide property.

Consider `tttttt` (6 t's) at positions 10-15, hp_threshold=4:

**If we used nucleotide-level masking**: flag positions 10-15. Then `mask_to_badness_mixed()` checks `any(mask[i:i+L])`:

| Probe | Window | Overlaps flagged pos? | Actual t-run in probe | Should mask? |
|-|-|-|-|-|
| pos=13, len=18 | 13-30 | Yes (pos 13,14,15) | `ttt` (3 t's) | **No** — below threshold |
| pos=8, len=18 | 8-25 | Yes (pos 10-15) | `tttttt` (6 t's) | Yes |
| pos=14, len=18 | 14-31 | Yes (pos 14,15) | `tt` (2 t's) | **No** — below threshold |

The nucleotide-level mask kills all three. Only the second deserves it.

**Root cause**: For R/P/B, one masked nucleotide in the probe window = bad probe. For L-mask, a single flagged nucleotide means nothing — you need a *complete run of threshold length* within the probe's own substring. Partial overlaps with a longer genomic run are not a problem.

**If we use probe-level checking**: evaluate each probe substring directly:

| Probe | Substring | `has_homopolymer(substr, 4)`? | Masked? |
|-|-|-|-|
| pos=13, len=18 | `tttXYZ...` | False (only 3 t's) | No ✓ |
| pos=8, len=18 | `XXttttttXX...` | True (6 t's) | Yes ✓ |
| pos=14, len=18 | `ttXYZ...` | False (only 2 t's) | No ✓ |

This maximises viable probe discovery — the stated goal.

### MATLAB Ancestor Validation

The MATLAB code (`probedesign/findprobesLocal.m`) already implements both architectures and explicitly documents the distinction. The comment at line 225 states:

```matlab
% NOTE: This is not really a mask. It gives a score per oligo.
% In order to use this, it really should be added to badness
% directly without extending the masking to the length of the oligo.
```

The badness integration (line 303-304) separates the two:

```matlab
badness = badness + mask_to_badness(fullmask, oligolen);  % R/P/B: nucleotide → probe expansion
badness = badness + GCmask + GCrunmask + ...;             % Quality: already per-oligo, added directly
```

The MATLAB `mask_oligos_with_runs()` function is the direct ancestor of our L-mask approach. It iterates per probe start position, extracts the probe substring `inseq(i:i+oligolen-1)`, and checks for runs within it — probe-level, not nucleotide-level. The MATLAB version only checks C and G runs (≥7 bp, 2 mismatches allowed). Our L-mask extends this to all 4 nucleotides and adds dinucleotide repeat detection.

`GC_badness()` follows the same probe-level pattern: per-position loop, extract `inseq(i:(i+oligolen-1))`, evaluate the substring. Same architecture as F-mask.

### Summary: Mask Architecture Selection

| Mask | Architecture | MATLAB precedent | Why |
|-|-|-|-|
| R / P / B | Nucleotide-level → `mask_to_badness()` | `fullmask` → `mask_to_badness(fullmask, oligolen)` | Nucleotide intrinsically problematic. One bad nt in probe = bad probe. |
| F | Probe-level (direct badness) | `GC_badness()` per-oligo loop | Thermodynamic properties are substring-dependent. |
| **L** | **Probe-level (direct badness)** | **`mask_oligos_with_runs()` per-oligo loop** | **Low-complexity is substring-dependent. Partial overlap with a genomic run ≠ probe problem.** |

---

## 4. Mixed-Length Mode: Per-(Position, Length) L-Masking

### The Challenge

In mixed-length mode (e.g. 18-22 bp), a probe starting at position `i` can have different lengths. Whether the probe contains a low-complexity pattern depends on **both** its start position and its length.

Example with `hp_threshold=4`, sequence region `...NNNNNNNNNNNNNNattttcg...`:

| Probe | Substring ending | Contains `tttt`? | L-masked? |
|-|-|-|-|
| pos=0, len=17 | `...NNNNNNNNNNNNNNatt` | No (2 t's) | No |
| pos=0, len=18 | `...NNNNNNNNNNNNNNattt` | No (3 t's) | No |
| pos=0, len=19 | `...NNNNNNNNNNNNNNatttt` | Yes (4 t's) | **Yes** |
| pos=0, len=20 | `...NNNNNNNNNNNNNNattttc` | Yes (4 t's) | **Yes** |

With a flat nucleotide-level mask, ALL lengths at position 0 would be masked because the probe window overlaps the flagged `tttt` nucleotide positions. This **loses the 17-nt and 18-nt probes** that are perfectly viable.

### The Solution: Direct Badness Modification Per (Position, Length)

For mixed-length mode, the L-mask check iterates over every `(position, length)` pair in the 2D badness table and checks the specific probe substring:

```python
for i in range(max_goodlen):
    for L_idx in range(n_lengths):
        L = min_len + L_idx
        if i + L > len(seq):
            continue
        probe_substr = seq[i:i+L]
        if has_homopolymer(probe_substr, hp_threshold) or has_dinuc_repeat(probe_substr, di_threshold):
            badness_2d[i][L_idx] = float('inf')
```

For fixed-length mode, the same logic applies but with a single length:

```python
for i in range(goodlen):
    probe_substr = seq[i:i+oligo_length]
    if has_homopolymer(probe_substr, hp_threshold) or has_dinuc_repeat(probe_substr, di_threshold):
        badness[i] = float('inf')
```

### Why This Is Architecturally Clean

The L-mask modifies `badness[]` / `badness_2d[]` **directly**, just like the F-mask does (F-mask: `badness = inf` when Gibbs FE is out of range). It does NOT go through `full_mask[]` → `mask_to_badness()`. This means:

- `full_mask[]` remains purely for nucleotide-level biological specificity masks (R/P/B)
- L-mask effects are baked into badness before `full_mask[]` is applied
- R/P/B masks can still add `inf` on top (harmless — inf stays inf)
- No changes to `mask_to_badness()`, `mask_to_badness_mixed()`, or the DP algorithm
- The DP sees badness values and picks optimal positions regardless of why a position is inf

### Computational Cost

For mixed-length mode: `O(max_goodlen × n_lengths × max_len)` character comparisons.

Typical case (seq_len=5000, n_lengths=5, max_len=22): ~550K operations. This is the same order as `calculate_badness_mixed()` itself. **Negligible overhead.**

For fixed-length mode: `O(goodlen × oligo_length)` — same order as `calculate_badness()`. Negligible.

No impact on DP complexity (which is `O(seq_len × n_probes)` for fixed, `O(seq_len × n_lengths × n_probes)` for mixed).

---

## 5. Low-Complexity Detection Functions

### 5A. Homopolymer Check (Per-Substring)

Check if a given string contains a homopolymer run ≥ threshold. This operates on the probe substring, not the full sequence.

```python
def has_homopolymer(subseq: str, threshold: int = 4) -> bool:
    """Return True if subseq contains a single-nucleotide run >= threshold.

    Only counts runs of a/c/g/t. Characters like '>', 'n' break runs.
    """
    if threshold < 1 or len(subseq) < threshold:
        return False
    count = 1
    for i in range(1, len(subseq)):
        if subseq[i] in 'acgt' and subseq[i] == subseq[i-1]:
            count += 1
            if count >= threshold:
                return True
        else:
            count = 1
    return False
```

Examples (threshold=4):
- `has_homopolymer("actgatttttcg")` → True  (5 t's)
- `has_homopolymer("actgatttcgat")` → False (only 3 t's)
- `has_homopolymer("actt>ttcgatc")` → False (`>` breaks the run)

### 5B. Dinucleotide Repeat Check (Per-Substring)

Check if a given string contains a dinucleotide motif repeated ≥ threshold times consecutively.

A "dinucleotide repeat of 3" means the 2-nt motif appears 3+ times = 6+ nucleotides (e.g. `ATATAT`).

```python
def has_dinucleotide_repeat(subseq: str, threshold: int = 3) -> bool:
    """Return True if subseq contains any dinucleotide motif repeated >= threshold times.

    Checks all 12 ordered dinucleotide pairs (AT, TA, AC, CA, AG, GA,
    TC, CT, TG, GT, GC, CG). Only considers a/c/g/t characters.
    """
    min_span = threshold * 2  # e.g. threshold=3 → need 6 consecutive nt
    if len(subseq) < min_span:
        return False

    for i in range(len(subseq) - 1):
        a, b = subseq[i], subseq[i+1]
        if a not in 'acgt' or b not in 'acgt' or a == b:
            continue  # skip non-nucleotides, junctions, and same-char (that's homopolymer)
        # Count consecutive repetitions of the motif (a,b) starting at position i
        count = 1
        j = i + 2
        while j + 1 < len(subseq) and subseq[j] == a and subseq[j+1] == b:
            count += 1
            j += 2
        if count >= threshold:
            return True
    return False
```

Examples (threshold=3):
- `has_dinucleotide_repeat("acgatatatcg")` → True (3× AT = `atatat`)
- `has_dinucleotide_repeat("acgatatcgat")` → False (only 2× AT = `atat`)
- `has_dinucleotide_repeat("gcgcgcgcgcg")` → True (5× GC)

### 5C. No "Combined Position Mask" Function Needed

Unlike the original plan, there is no `low_complexity_mask()` function that returns a nucleotide-level position array. The check is applied directly per probe substring in `core.py`. The helper functions `has_homopolymer()` and `has_dinucleotide_repeat()` are simple predicate functions.

---

## 6. File-by-File Implementation Plan

### 6A. `src/probedesign/masking.py` — New Functions

Add two predicate functions:

```python
def has_homopolymer(subseq: str, threshold: int = 4) -> bool:
    """True if subseq contains a single-nucleotide run >= threshold."""

def has_dinucleotide_repeat(subseq: str, threshold: int = 3) -> bool:
    """True if subseq contains a dinucleotide motif repeated >= threshold times."""
```

These are pure string-checking functions with no side effects. Place them near the top of `masking.py` (they don't depend on bowtie or any external tools).

### 6B. `src/probedesign/core.py` — Integration into `design_probes()`

**New parameters** for `design_probes()`:

```python
def design_probes(
    ...,
    hp_threshold: int = 4,        # NEW
    di_threshold: int = 3,        # NEW
) -> ProbeDesignResult:
```

**Insert L-mask block** after both the F-mask visualization string is created AND after the badness calculation, but BEFORE the R/P/B masking section. The insertion point differs slightly between the fixed-length and mixed-length branches.

#### Fixed-length branch (after F mask visualization, before `if any(full_mask):`):

```python
# --- L-mask: low-complexity filtering (always active, probe-level) ---
from .masking import has_homopolymer, has_dinucleotide_repeat
l_masked_count = 0
lstr_parts = []
for i in range(len(seq)):
    if i < goodlen:
        probe_sub = seq[i:i + oligo_length]
        is_l_masked = (
            has_homopolymer(probe_sub, hp_threshold)
            or has_dinucleotide_repeat(probe_sub, di_threshold)
        )
        if is_l_masked:
            badness[i] = float('inf')
            l_masked_count += 1
            lstr_parts.append('L')
        else:
            lstr_parts.append(seq[i])
    else:
        lstr_parts.append(seq[i])
lstr = "".join(lstr_parts)
mask_strings.append(lstr)
if l_masked_count > 0:
    print(f"Low-complexity masking: {l_masked_count} probe positions masked "
          f"(homopolymer>={hp_threshold}, dinucleotide>={di_threshold})")
```

#### Mixed-length branch (after F mask visualization, before `if any(full_mask):`):

```python
# --- L-mask: low-complexity filtering (always active, per-length) ---
from .masking import has_homopolymer, has_dinucleotide_repeat
l_masked_count = 0
l_all_masked_count = 0  # positions where ALL lengths are L-masked
lstr_parts = []
for i in range(len(seq)):
    if i < max_goodlen:
        all_lengths_masked = True
        any_length_masked = False
        for L_idx in range(n_lengths):
            L = min_len + L_idx
            if i + L > len(seq):
                continue  # already inf from calculate_badness_mixed
            if badness_2d[i][L_idx] == float('inf'):
                continue  # already excluded by F-mask, don't count
            probe_sub = seq[i:i + L]
            is_l = (
                has_homopolymer(probe_sub, hp_threshold)
                or has_dinucleotide_repeat(probe_sub, di_threshold)
            )
            if is_l:
                badness_2d[i][L_idx] = float('inf')
                l_masked_count += 1
                any_length_masked = True
            else:
                all_lengths_masked = False
        # Visualization: 'L' only if ALL lengths at this position are now dead
        # (either F-masked or L-masked — no viable length remains)
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
if l_masked_count > 0:
    print(f"Low-complexity masking: {l_masked_count} (position,length) pairs masked "
          f"(homopolymer>={hp_threshold}, dinucleotide>={di_threshold})")
```

**Key points about mixed-length visualization**:
- A position shows `L` only if the L-mask caused it to lose ALL remaining viable lengths (i.e., lengths that weren't already F-masked). This is consistent with how F-mask shows `F` only if NO length gives finite badness.
- If some lengths at a position are L-masked but others survive, the position shows the sequence character — the user can see which specific probes were selected in the oligos output.

**Mask string ordering**: After this change, `mask_strings` will be built in order:
1. R-mask string (if repeat masking active)
2. P-mask string (if pseudogene masking active)
3. B-mask string (if genome masking active)
4. F-mask string (always)
5. L-mask string (always)

The `_seq.txt` output shows: sequence → R → P → B → F → L → probes.
The two always-active quality masks (F and L) are grouped together at the bottom, closest to probe alignments. Optional biological masks (R, P, B) are grouped above them.

### 6C. `src/probedesign/cli.py` — New CLI Options

Add two new Click options to the `design` command:

```python
@click.option(
    '--hp-threshold',
    default=4,
    type=click.IntRange(3, 10),
    help='Homopolymer repeat threshold for L-mask (default: 4, range: 3-10). '
         'Probes containing single-nucleotide runs >= this length are excluded.'
)
@click.option(
    '--di-threshold',
    default=3,
    type=click.IntRange(3, 10),
    help='Dinucleotide repeat threshold for L-mask (default: 3, range: 3-10). '
         'Probes containing dinucleotide motifs repeated >= this many times are excluded.'
)
```

Pass both through to `design_probes()`:

```python
result = design_probes(
    ...,
    hp_threshold=hp_threshold,
    di_threshold=di_threshold,
)
```

Update the `design()` function signature to include `hp_threshold: int` and `di_threshold: int`.

### 6D. `streamlit_app/app.py` — New UI Section

Add a "Low-complexity filter" section in the sidebar **above** the Pseudogene mask checkbox (between the `st.divider()` after species selection and the bowtie masking section, around line 145):

```python
st.divider()

# Low-complexity filter — always active
st.markdown("**Low-complexity filter**")
col_hp, col_di = st.columns(2)
with col_hp:
    hp_threshold = st.number_input(
        "Homopolymer repeats",
        min_value=3, max_value=10, value=4, step=1,
        key="hp_threshold",
        help="Mask probes containing single-nucleotide runs of this length or longer (e.g. AAAA at threshold 4)"
    )
with col_di:
    di_threshold = st.number_input(
        "Dinucleotide repeats",
        min_value=3, max_value=10, value=3, step=1,
        key="di_threshold",
        help="Mask probes containing dinucleotide motifs repeated this many times or more (e.g. ATATAT at threshold 3)"
    )

st.divider()

# Bowtie masking (existing section follows)
```

**Defaults**: hp_threshold=4, di_threshold=3 in both CLI and GUI for consistency.

### 6E. `streamlit_app/utils.py` — Pass-Through Parameters

**`run_design()`** — Add parameters:

```python
def run_design(
    ...,
    hp_threshold: int = 4,
    di_threshold: int = 3,
) -> DesignRunResult:
```

Pass them to `design_probes()`:

```python
result = design_probes(
    ...,
    hp_threshold=hp_threshold,
    di_threshold=di_threshold,
)
```

**`run_batch()`** — Extract from params dict (same pattern as `mixed_lengths`):

```python
hp_threshold = params.get("hp_threshold", 4)
di_threshold = params.get("di_threshold", 3)
```

Pass them to `run_design()`.

**Single-mode call site** in `app.py` (line ~433):

```python
run_result = run_design(
    ...,
    hp_threshold=hp_threshold,
    di_threshold=di_threshold,
)
```

**Batch-mode params dict** in `app.py` (line ~545):

```python
params = {
    ...,
    "hp_threshold": hp_threshold,
    "di_threshold": di_threshold,
}
```

---

## 7. Output Format Changes

### 7A. `_seq.txt` Visualization

The L-mask line appears after the F-mask line. It shows `L` at **probe-starting positions** where the L-mask caused exclusion, and the original sequence character elsewhere.

Example: sequence `actgatcgtttttactgatcgatatatcg` (28 nt), oligo_length=18, hp_threshold=4, di_threshold=3.

Consider probe at position 0 (length 18): `actgatcgtttttactga` — contains `ttttt` (5 t's ≥ 4) → **L-masked**.
Consider probe at position 5 (length 18): `tcgtttttactgatcgat` — contains `ttttt` (5 t's ≥ 4) → **L-masked**.
Consider probe at position 10 (length 18): `tactgatcgatatatcg.` — wait, only 18 chars from pos 10 = `tttactgatcgatatat` ... this contains `atatat` (3× AT dinuc) → **L-masked**.
Consider probe at position 1 (length 18): `ctgatcgtttttactgat` — contains `ttttt` → **L-masked**.

```
actgatcgtttttactgatcgatatatcg     (original sequence)
actgatcgtttttactgatcgatatatcg     (R mask — none in this example)
actgatcgFFFttactgatcgatatatcg     (F mask)
LLLLLLLLLLLactgatcgatatatLLLL     (L mask — starting positions where probe is masked)
           tgactagctatatag         (probe complements)
           Prb# 1,Pos 12,...       (probe labels)
```

Note the difference from the original (incorrect) plan: L-mask shows `L` at each **starting position** where a probe of that length would contain low-complexity sequence, not at the nucleotide positions of the repeat itself. Position 11 (`a` in `actgatcg...`) is fine because an 18-nt probe starting there (`actgatcgatatatcg...`) doesn't contain enough repeats.

### 7B. Mixed-Length Visualization

In mixed-length mode, a position shows `L` only if ALL lengths at that position are dead (F-masked or L-masked). Positions where some lengths survive show the sequence character.

Example: position 0, lengths 17-20, sequence `...NNNNNNNNNNNNNNattttcg...`:

| Length | Probe substring | Has `tttt`? | L-masked? |
|-|-|-|-|
| 17 | `...NNNNNNNNNNNNNNatt` | No | No |
| 18 | `...NNNNNNNNNNNNNNattt` | No | No |
| 19 | `...NNNNNNNNNNNNNNatttt` | Yes | Yes |
| 20 | `...NNNNNNNNNNNNNNattttc` | Yes | Yes |

Position 0 visualization: shows sequence char (not `L`) because lengths 17 and 18 are still viable. The DP can select a 17- or 18-nt probe here.

### 7C. No Changes to `_oligos.txt`

The oligos file format is unchanged — it only lists selected probes. L-masked (position, length) pairs get `inf` badness and are excluded by the DP.

---

## 8. Compatibility Analysis

### 8A. Fixed-Length Code Path

- L-mask modifies `badness[i]` directly (sets to inf for L-masked probes)
- `full_mask[]` is NOT involved — L-mask does not contribute to `full_mask[]`
- The subsequent `mask_to_badness()` call for R/P/B masks applies on top (inf + inf = inf, harmless)
- No changes to DP or mask-to-badness conversion

### 8B. Mixed-Length Code Path

- L-mask modifies `badness_2d[i][L_idx]` per (position, length) pair
- Each length is checked independently — shorter probes that don't reach a repeat are preserved
- `mask_to_badness_mixed()` for R/P/B still applies on top via `full_mask[]`
- No changes to DP: the mixed-length DP (`find_best_probes_mixed`) already iterates over all (position, length) pairs and skips inf entries

### 8C. F-Mask Interaction

- F-mask visualization is computed BEFORE L-mask runs → F-mask output is unaffected
- The L-mask skips entries already at inf (from F-mask) in its counting — prevents double-reporting in the log message
- In mixed mode, the L visualization checks what's viable after both F and L, but only shows `L` if L-mask was the reason (not F-mask) → the `any_length_masked` flag ensures this

### 8D. Backward Compatibility

- Default thresholds (hp=4, di=3) mean L-mask is always active. This **will change output** for sequences containing probes with homopolymer runs ≥4 or dinucleotide repeats ≥3
- Existing test cases may need updated expected outputs
- To effectively disable: `--hp-threshold 10 --di-threshold 10`

### 8E. Test Case Impact

Check each test case for low-complexity probes:
- `CDKN1A_32/` — scan for probes containing homopolymers ≥4 / dinuc repeats ≥3
- `KRT19_withUTRs/` — same
- `EIF1_CDS_HCR/` — same (52-nt probes are more likely to contain repeats)
- If any test case is affected, expected output needs updating or the test should be re-run with documented new expected values

---

## 9. Strategy Comparison: Why Per-(Position, Length) Is the Right Choice

Three strategies were considered for mixed-length mode:

### Strategy A: Flat Nucleotide-Level Mask (REJECTED)

Flag nucleotide positions that are part of repeats in the full sequence, then use `mask_to_badness_mixed()`.

**Fatal flaw**: Over-masks. A probe that only overlaps 3 nt of a 4-nt homopolymer run gets masked even though its own substring only has 3 consecutive same-nucleotides. Loses viable probes.

Example: sequence has `tttt` at positions 15-18. A probe starting at position 16 of length 18 covers `tttcg...` — only 3 t's, below threshold. But nucleotide-level mask flags position 16, killing this probe.

### Strategy B: Separate Full L-Mask Per Length (REJECTED)

Pre-compute a full position-level mask for each possible length, using probe-substring checks. Store as `lmask[L_idx][pos]`. Apply each length's mask independently.

**Works correctly** but adds unnecessary complexity and memory. For 5 lengths and 5000 positions, that's 25,000 entries — not a memory concern, but the code is more complex than needed for no benefit.

### Strategy C: Direct Per-(Position, Length) Badness Modification (CHOSEN)

Check each `(position, length)` pair in the inner loop and set `badness_2d[i][L_idx] = inf` directly.

**Advantages**:
- Simplest code — single pass over the badness table
- No intermediate data structures
- Naturally handles mixed-length mode without special cases
- Fixed-length mode is a degenerate case (single length)
- Same computational complexity as Strategy B
- Modification of badness directly is architecturally consistent with how F-mask works

**No DP changes needed**: The DP already handles per-(position, length) inf values. It simply skips them. Masking a subset of lengths at a position doesn't break the DP — it just reduces the candidates at that position, potentially allowing shorter/longer probes to be chosen instead.

---

## 10. Documentation Updates

### 10A. `docs/probe_design_principles.md`

1. **Pipeline Overview** (Section 1): Insert `[4b] L mask` between `[4] F mask` and `[5] Apply Optional Masks`
2. **New Section 4B**: "Low-Complexity Filtering (L Mask)" — describe probe-level homopolymer and dinucleotide repeat checks, thresholds, examples, and the distinction from nucleotide-level R/P/B masks
3. **Section 5 "What Is NOT Checked"**: Remove "Homopolymer runs — Not checked" row
4. **Section 6 Specificity Filtering**: Add note that L-mask is applied before R/P/B masks, is always active, and operates at the probe level (not nucleotide level)
5. **Section 8 Output**: Update `_seq.txt` format table to include L mask line
6. **Section 9 Parameter Reference**: Add `hp_threshold` and `di_threshold` to CLI parameters table

### 10B. `README.md`

Add `--hp-threshold` and `--di-threshold` to CLI usage examples and parameter table. Mention low-complexity filter in feature list.

### 10C. `CLAUDE.md`

1. **Algorithm section**: Update step list to include L-mask
2. **Architecture table**: Add `has_homopolymer()` and `has_dinucleotide_repeat()` to masking.py description
3. **Streamlit App section**: Mention Low-complexity filter UI
4. **Masking Details**: Add L-mask description noting its probe-level (not nucleotide-level) nature

---

## 11. Implementation Order

1. **`masking.py`** — Add `has_homopolymer()` and `has_dinucleotide_repeat()` predicate functions
2. **`core.py`** — Add `hp_threshold` and `di_threshold` parameters to `design_probes()`. Insert L-mask block in both fixed-length and mixed-length branches (after F-mask visualization, before R/P/B masking)
3. **`cli.py`** — Add `--hp-threshold` and `--di-threshold` Click options. Pass to `design_probes()`
4. **Test CLI** — Run existing test cases, verify L-mask appears in output, check for regressions
5. **`streamlit_app/utils.py`** — Add parameters to `run_design()` signature and pass-through
6. **`streamlit_app/app.py`** — Add Low-complexity filter UI section in sidebar. Pass values to `run_design()` and batch `params` dict
7. **Test Streamlit** — Verify UI renders correctly, parameters flow through, results include L-mask
8. **Documentation** — Update `probe_design_principles.md`, `README.md`, `CLAUDE.md`

---

## 12. Edge Cases & Considerations

| Case | Behavior |
|-|-|
| Probe substring contains `>` (junction marker) | `>` is not in `acgt`, so it breaks homopolymer runs and doesn't form dinucleotide motifs. These probes already have `inf` badness from F-mask, so L-mask setting inf on top is harmless. |
| Probe substring contains `n` | `n` is not in `acgt`, so it breaks runs. The L-mask won't false-positive on `nnn` sequences. |
| Threshold at minimum (3) | Homopolymer: masks probes with 3+ same-nt runs. Dinucleotide: masks probes with 3+ motif repeats (6+ nt). |
| Threshold at maximum (10) | Effectively disables the filter for most biological sequences. |
| Overlapping dinucleotide patterns | `atatatatat` in a probe matches both AT-repeat (starting at even positions) and TA-repeat (starting at odd positions). Either match triggers the mask — result is the same. |
| Very short probes | If probe length < threshold (e.g., 3-nt probe with threshold 4), no homopolymer can exist → never masked. |
| Mixed-length: repeat at boundary | The critical case. E.g., `...attt` (18nt: 3 t's → OK) vs `...atttt` (19nt: 4 t's → masked). Per-(position, length) check handles this correctly, preserving the shorter probe. |
| F-mask already killed a (pos, len) pair | L-mask check still runs (cheap) and sets inf (harmless). The log message counts only L-mask-specific additions by skipping entries already at inf before L-mask. |

---

## 13. Summary of Changes by File

| File | Changes |
|-|-|
| `src/probedesign/masking.py` | Add 2 predicate functions: `has_homopolymer(subseq, threshold)`, `has_dinucleotide_repeat(subseq, threshold)` |
| `src/probedesign/core.py` | Add `hp_threshold`, `di_threshold` params to `design_probes()`. Insert L-mask block in both fixed-length and mixed-length branches: check each probe substring, set badness=inf, build visualization string, append to `mask_strings`. Does NOT touch `full_mask[]`. |
| `src/probedesign/cli.py` | Add `--hp-threshold` (IntRange 3-10, default 4) and `--di-threshold` (IntRange 3-10, default 3) options. Pass to `design_probes()`. |
| `streamlit_app/app.py` | Add "Low-complexity filter" section in sidebar with two number inputs. Pass values to `run_design()` and batch params dict. |
| `streamlit_app/utils.py` | Add `hp_threshold`, `di_threshold` to `run_design()` signature. Extract from params in `run_batch()`. Pass through to `design_probes()`. |
| `docs/probe_design_principles.md` | Add Section 4B, update pipeline diagram, update parameter tables, remove "Homopolymer — Not checked" from Section 5, update `_seq.txt` format. |
| `README.md` | Add CLI options, mention low-complexity filter in features. |
| `CLAUDE.md` | Update algorithm steps, add masking.py functions, mention Streamlit UI. |
