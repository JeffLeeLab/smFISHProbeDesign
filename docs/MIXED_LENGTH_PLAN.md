# Mixed-Length Probe Design: Implementation Plan

> **Status**: Implemented (`383cae7`, 2026-02-19). Historical design plan: where it disagrees with the code, the code wins. Current behaviour is in [probe_design_principles.md](probe_design_principles.md#7-dynamic-programming-optimization).

## 1. Critical Analysis: Does the Logic Make Sense?

**The logic is sound.** The nearest-neighbor thermodynamic model (Sugimoto 1995) sums per-dinucleotide contributions, so Gibbs FE scales roughly linearly with probe length. For any given sequence composition:
- GC-rich regions produce more negative Gibbs FE per bp → shorter probes (18bp) can meet the target range [-26, -20]
- AT-rich regions produce less negative Gibbs FE per bp → longer probes (22bp) are needed to reach the target range

With a fixed 20bp length, positions in extreme GC/AT regions get `badness=inf` (Gibbs outside allowable range). Mixed lengths 18-22 can recover many of these positions, increasing valid probe locations and maximizing probe count.

### Disadvantages

| Concern | Severity | Mitigation |
|---------|----------|------------|
| **Non-uniform Tm across probe set** | Low | smFISH hybridization at 37°C in formamide is well below Tm for all 18-22mers; the Gibbs FE constraint already ensures similar binding thermodynamics |
| **Shorter probes (18bp) have less sequence specificity** | Low | Bowtie masking uses 12-16bp seeds (shorter than any probe), so off-target filtering is already conservative |
| **No MATLAB reference for validation** | Medium | Test that when min_len==max_len, output is identical to existing fixed-length code |
| **Increased algorithmic complexity** | Low | DP complexity increases by constant factor (5x for 5 lengths), trivial for typical sequence lengths |
| **Pooled oligo synthesis incompatibility** | Low-Medium | Some synthesis platforms require uniform length; document that mixed-length probes should be ordered individually |
| **Signal uniformity across probe set** | Low | Binding energy (Gibbs FE) is constrained to the same range for all lengths, so per-molecule signal should be similar |

## 2. Implementation Compatibility with Existing Code

### Components that need modification

| Component | Current behavior | Required change |
|-----------|-----------------|-----------------|
| `calculate_badness()` in `core.py:51-100` | Single badness value per position (1D array) | New `calculate_badness_mixed()` returns 2D array `badness[pos][len_idx]` |
| `find_best_probes()` in `core.py:120-206` | DP with fixed spacing `probe_spacer_len` | New `find_best_probes_mixed()` using end-position DP with variable spacing |
| `mask_to_badness()` in `masking.py:351-374` | Checks `[i, i+oligo_length)` for each position | New `mask_to_badness_mixed()` checks per (position, length) pair |
| F mask generation in `core.py:349-356` | `goodlen = len(seq) - oligo_len + 1` | Position valid if ANY length gives finite badness |
| Probe extraction in `core.py:390-410` | Uses global `oligo_length` | Uses per-probe length from DP result |
| `output.py:73` | Gets `oligo_len` from first probe | Uses `len(probe.sequence)` per probe |
| `output.py:85-89` | Collects `oligo_len` chars | Already uses `len(probe.sequence)` indirectly via complement, but needs explicit per-probe handling |
| `cli.py:25-30` | `--oligo-length` accepts single int | Accept `"18-22"` range string OR single int |

### Components that DON'T need modification

| Component | Why unchanged |
|-----------|---------------|
| `thermodynamics.py` | `gibbs_rna_dna()` works on any sequence length; no change needed |
| `sequence.py` | Generic utilities (`reverse_complement`, `percent_gc`, `has_invalid_chars`) are length-agnostic |
| `fasta.py` | Sequence I/O is independent of probe length |
| Position-level mask generation (`masking.py:86-169`) | Bowtie alignment uses fixed n-mer lengths (12-16bp), independent of probe length |
| `hits_to_mask()` in `masking.py` | Produces position-level binary mask, not probe-level |

## 3. Bowtie Filtering Implications

**No changes needed to bowtie parameters or thresholds.**

The bowtie masking operates at the sub-probe level:
- Pseudogene masking: 16-mers with 0 mismatches → position-level binary mask
- Genome masking: 12/14/16-mers with tiered thresholds → position-level binary mask

These position-level masks are **independent of probe length**. The conversion from position-level mask to probe-level badness (`mask_to_badness`) is the only step that depends on probe length, and it needs a per-length version.

**Specificity analysis for 18-22bp probes:**
- An 18bp probe with a 16-mer pseudogene hit: 16/18 = 89% identity → stringent (correct)
- A 22bp probe with a 16-mer pseudogene hit: 16/22 = 73% identity → less stringent proportionally

This means shorter probes are *more* stringently filtered as a fraction of their length, which is appropriate since shorter probes inherently have less sequence specificity. The existing thresholds work well for the 18-22bp range.

## 4. Comprehensive Implementation Plan

### Phase 1: CLI Changes (`cli.py`)

**File**: `src/probedesign/cli.py`

Change the `--oligo-length` option from `type=int` to `type=str` with a custom parser:

```python
@click.option('-l', '--oligo-length', default='20',
    help='Oligo length: single int (e.g. 20) or range (e.g. 18-22)')
```

Parse logic:
- `"20"` → `oligo_length=20, mixed_lengths=None` (backward compatible)
- `"18-22"` → `oligo_length=None, mixed_lengths=(18, 22)`

Validation: min_len >= 15, max_len <= 30, min_len <= max_len, max_len - min_len <= 10.

Pass `mixed_lengths` tuple to `design_probes()` as a new optional parameter.

### Phase 2: Core Algorithm (`core.py`)

#### 2a. New `calculate_badness_mixed()` function

```python
def calculate_badness_mixed(
    seq: str,
    min_len: int, max_len: int,
    target_gibbs: float,
    allowable_range: Tuple[float, float],
) -> List[List[float]]:
    """Returns badness_2d[pos][len_idx] for lengths min_len..max_len."""
```

For each position `i` in `range(len(seq) - min_len + 1)`:
  For each length `L` in `range(min_len, max_len + 1)`:
    - If `i + L > len(seq)`: `inf`
    - Extract `oligo = seq[i:i+L]`
    - If `has_invalid_chars(oligo)`: `inf`
    - Compute Gibbs FE; if outside allowable range: `inf`
    - Otherwise: `(gibbs - target)^2`

#### 2b. New `find_best_probes_mixed()` function (End-Position DP)

**State**: `dp[e][k]` = best average badness for k+1 probes with last probe ending at or before position e.

**Tracker**: `tracker[e][k] = (start_pos, length)` for backtracking.

**Algorithm**:
```
for e in range(seq_len):
    # Propagate: inherit best from previous end position
    if e > 0:
        for k in range(n_probes):
            if dp[e-1][k] < dp[e][k]:
                dp[e][k] = dp[e-1][k]
                tracker[e][k] = tracker[e-1][k]

    # Try placing probe ending at exactly position e
    for L in range(min_len, max_len + 1):
        x = e - L + 1  # start position
        if x < 0: continue
        L_idx = L - min_len
        if badness_2d[x][L_idx] == inf: continue

        for k in range(n_probes):
            if k == 0:
                score = badness_2d[x][L_idx]
            else:
                prev_end = x - spacer_len - 1
                if prev_end < 0 or tracker[prev_end][k-1] is None:
                    continue
                score = _calc_score(dp[prev_end][k-1], k, badness_2d[x][L_idx])

            if score < dp[e][k]:
                dp[e][k] = score
                tracker[e][k] = (x, L)
```

**Backtracking** (for each target probe count):
```
e = seq_len - 1
positions_and_lengths = []
for curr_k in range(target_k, -1, -1):
    if tracker[e][curr_k] is None: break
    x, L = tracker[e][curr_k]
    positions_and_lengths.append((x, L))
    e = x - spacer_len - 1
positions_and_lengths.reverse()
```

**Returns**: `List[Tuple[float, List[Tuple[int, int]]]]` — list of (score, [(start, length), ...]) for 1..n_probes solutions.

**Complexity**: O(seq_len × n_lengths × n_probes) ≈ O(N × 5 × K). For typical N=2000, K=48: ~480,000 iterations.

#### 2c. Update `design_probes()` signature and logic

Add optional parameter `mixed_lengths: Optional[Tuple[int, int]] = None`.

When `mixed_lengths` is provided:
1. Call `calculate_badness_mixed()` instead of `calculate_badness()`
2. F mask: show 'F' only if ALL lengths give `inf` at that position
3. `mask_to_badness_mixed()` for mask integration
4. Call `find_best_probes_mixed()` instead of `find_best_probes()`
5. Probe extraction uses per-probe length from DP results

When `mixed_lengths` is None: existing behavior unchanged (backward compatible).

#### 2d. Update `Probe` dataclass

Add `length: int` field. Set to `len(sequence)` for both fixed and mixed modes:

```python
@dataclass
class Probe:
    index: int
    position: int
    length: int          # NEW
    sequence: str
    gc_percent: float
    tm: float
    gibbs_fe: float
    name: str
```

### Phase 3: Masking Changes (`masking.py`)

Add `mask_to_badness_mixed()`:

```python
def mask_to_badness_mixed(
    mask: List[int], min_len: int, max_len: int
) -> List[List[float]]:
    """Convert position-level mask to 2D probe-level badness."""
    max_goodlen = len(mask) - min_len + 1
    n_lengths = max_len - min_len + 1
    badness = [[0.0] * n_lengths for _ in range(max_goodlen)]
    for i in range(max_goodlen):
        for L_idx, L in enumerate(range(min_len, max_len + 1)):
            if i + L > len(mask):
                badness[i][L_idx] = float('inf')
            elif any(mask[i:i + L]):
                badness[i][L_idx] = float('inf')
    return badness
```

### Phase 4: Output Changes (`output.py`)

#### 4a. `write_oligos_file()` / `format_oligos_content()`

Add `Length` column to header. Each probe's length comes from `len(probe.sequence)` or `probe.length`.

New format: `index  start  Length  GC%  Tm  GibbsFE  sequence  name`

#### 4b. `write_seq_file()` / `format_seq_content()`

Change `oligo_len = len(result.probes[0].sequence)` to per-probe: use `len(probe.sequence)` when extracting template and placing complement in the alignment visualization.

Key line changes:
- Line 73: Remove global `oligo_len` derivation
- Lines 85-89: Use `len(probe.sequence)` per probe in the collection loop
- Line 95: Use `len(probe.sequence) + 10` per probe for placement loop

### Phase 5: Testing

#### 5a. Backward Compatibility Test

When `mixed_lengths=(20, 20)` (min==max), verify output is **identical** to `oligo_length=20`. This validates that the new DP produces the same results as the old one when only one length is available.

Run with existing test cases:
```bash
# Fixed length (existing)
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa -l 20

# Mixed length with single value (should be identical)
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa -l 20-20
```

#### 5b. Mixed-Length Integration Test

Run with existing CDKN1A test case:
```bash
# Mixed length range
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --probes 48 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa -l 18-22
```

Verify:
- All probe sequences are 18-22bp
- All probe Gibbs FE values are within [-26, -20] range
- Probe count >= what fixed-length 20 produces
- No overlap between probes (spacing constraint met)
- Probe labels include length info

#### 5c. Minimal Synthetic Test Case

Create `test_cases/mixed_length_test/mixed_gc.fa`:
```
>mixed_gc_test
[A ~200bp sequence with alternating GC-rich and AT-rich regions]
```

Run with both fixed and mixed lengths to demonstrate the benefit:
```bash
# Fixed 20bp — some regions will have no valid probes
probedesign design test_cases/mixed_length_test/mixed_gc.fa --probes 20 -l 20

# Mixed 18-22bp — more regions accessible
probedesign design test_cases/mixed_length_test/mixed_gc.fa --probes 20 -l 18-22
```

Expected: mixed-length produces more probes than fixed-length.

#### 5d. Unit Tests

Add tests in a test file:
1. `test_calculate_badness_mixed()` — verify 2D badness table for known sequences
2. `test_find_best_probes_mixed_single_length()` — verify identical to `find_best_probes()` when min==max
3. `test_find_best_probes_mixed_variable()` — verify correct probe placement with different lengths
4. `test_mask_to_badness_mixed()` — verify masking with variable probe lengths
5. `test_cli_parse_oligo_length()` — verify parsing of "20", "18-22" formats

### Implementation Order

1. `core.py`: Add `Probe.length` field, update all existing Probe instantiations
2. `core.py`: Add `calculate_badness_mixed()`
3. `masking.py`: Add `mask_to_badness_mixed()`
4. `core.py`: Add `find_best_probes_mixed()`
5. `core.py`: Update `design_probes()` to support `mixed_lengths` parameter
6. `output.py`: Update output functions for per-probe lengths
7. `cli.py`: Update `--oligo-length` parsing
8. Run backward compatibility tests (existing test cases must still pass)
9. Create synthetic test case and run mixed-length tests
10. Update CLAUDE.md with new feature documentation

### Files to Modify

| File | Changes |
|------|---------|
| `src/probedesign/core.py` | Add `Probe.length`, `calculate_badness_mixed()`, `find_best_probes_mixed()`, update `design_probes()` |
| `src/probedesign/masking.py` | Add `mask_to_badness_mixed()` |
| `src/probedesign/output.py` | Per-probe length handling in seq/oligos output |
| `src/probedesign/cli.py` | Parse `--oligo-length` as int or range string |
| `test_cases/mixed_length_test/mixed_gc.fa` | New synthetic test FASTA |
| `CLAUDE.md` | Document new feature |

### Design Decision: Same Gibbs Target for All Lengths

The same `target_gibbs` and `allowable_gibbs` values are used for ALL lengths in the range. This is deliberate:
- The Gibbs FE constraint ensures similar binding thermodynamics regardless of probe length
- A valid 18-mer and a valid 22-mer in the same probe set will have similar hybridization energy
- This means experimental performance should be consistent across the probe set

A per-length-normalized target (e.g., target_gibbs_per_bp × L) is a possible future enhancement but not needed for v1.

### Design Decision: Global DP (All Lengths Considered at Every Position)

The DP considers ALL valid lengths at every position and finds the globally optimal combination of probe placements and lengths. This is preferred over pre-selecting the best length per position because:

1. **Maximizes probe count**: The global DP can choose a shorter probe at one position to make room for an additional probe nearby, even if a longer probe had slightly better badness at that position.
2. **Negligible extra cost**: The complexity is O(N × 5 × K) vs O(N × K) — a 5x constant factor that's trivial for sequences < 10kbp.
3. **Subsumes the local approach**: When the globally optimal length happens to equal the locally best length at every position (which is common), the results are identical.
4. **No need for separate modes**: A single code path is simpler to maintain and test.

### Design Decision: CLI Syntax

Use `-l 18-22` (dash-separated range) for mixed lengths. The existing `-l 20` (single integer) continues to work unchanged. Parsing logic:
- If the string contains `-` and both parts are valid integers with the first < second: treat as range
- Otherwise, try to parse as single integer for backward compatibility

### Verification Plan

1. **Backward compatibility**: Run all existing test cases (`run_tests.sh`) — must produce identical output
2. **Single-length equivalence**: `-l 20-20` must produce identical output to `-l 20`
3. **Mixed-length validation on existing data**: Run CDKN1A with `-l 18-22` and verify:
   - All probes are 18-22bp
   - All Gibbs FE values within allowable range
   - No probe overlaps (spacing constraint satisfied)
   - Probe count >= fixed-length result
4. **Synthetic test case**: Create a sequence with alternating GC-rich/AT-rich regions to demonstrate mixed-length benefit
5. **Unit tests**: Test individual functions (badness_mixed, DP, mask_to_badness_mixed, CLI parsing)
