# HCR Split-Initiator Probe Designer — Implementation Plan

> **Status**: Implemented (`8917ef4`, 2026-02-28; initiator orientation fix `edd2446`). Historical design plan: where it disagrees with the code, the code wins (e.g. there is no `apply_hcr_masks()`). Current behaviour is in [probe_design_principles.md](probe_design_principles.md#10-hcr-split-initiator-probes).

## 1. Critical Assessment of Design Principles

### 1.1 Independent hybridisation: correct

The user's central insight is sound. The existing codebase treats HCR probes as a single 52-nt oligo with `target_gibbs=-60`, which is **thermodynamically incorrect**. The two 25-nt halves are not covalently linked on the target strand; they hybridise as independent events separated by a 2-nt gap. Each half must be evaluated at its own ΔG (target ≈ −31 kcal/mol under 5×SSC / 30% formamide). The existing approach conflates two separate duplex stabilities into one number, which could select probes where one half is too weak and the other compensates — a pair that would fail in practice.

### 1.2 Asymmetric leniency: sound with caveats

The mechanistic argument is valid: a single half-initiator cannot nucleate HCR hairpin polymerisation regardless of how stably it binds off-target. Relaxing the strict half's ΔG ceiling or bowtie filters is therefore safe for specificity.

**Caveats to address in implementation:**

1. **Cross-pair off-target coincidence.** In theory, a lenient half from pair X could bind off-target adjacent to a lenient half from pair Y at some non-target locus, reconstituting a full initiator spuriously. In practice this requires two independent 25-mers to land within 2 nt of each other on a random transcript — astronomically unlikely given transcriptome size. No special handling required, but worth noting in documentation.

2. **Which half is lenient must be decided per-pair, not globally.** The user specifies "OR-mode" (either half can be lenient, not both). The implementation should evaluate both orientations (A-strict/B-lenient and A-lenient/B-strict) for each candidate pair and pick the one with lower combined badness. This doubles the candidate search space but does not double runtime since both evaluations use pre-computed per-position data.

3. **Intramolecular structure risk at very negative ΔG.** The user proposes allowing the lenient half up to −42 kcal/mol. Very GC-rich 25-mers can form stable hairpins or G-quadruplexes. The homopolymer/dinucleotide filters partially mitigate this, but a secondary-structure check (self-complementarity scan) could be added as a future enhancement. For v1, the homopolymer filter plus a GC% ceiling (e.g. ≤76%, i.e. ≤19/25 GC) is a reasonable guard.

### 1.3 DP as 52-nt block: correct approach

Treating each probe pair as an atomic 52-nt block for placement is the right abstraction. The DP should optimise over pair positions (where a "position" is the start of the 52-nt window), with a user-configurable minimum spacing between consecutive pair blocks.

### 1.4 Bowtie relaxation for lenient half: sound

Skipping or relaxing transcriptome-wide alignment filtering for the lenient half is mechanistically justified. The implementation should still run bowtie for both halves (to report hits transparently) but only apply mask penalties to the strict half. The user should be able to choose between symmetric and asymmetric bowtie filtering.

### 1.5 Potential pitfalls

| Pitfall | Mitigation |
|-|-|
| WW spacers for amplifiers B7+ are IUPAC ambiguity codes, not resolved sequences | Output `WW` directly in oligo sequences — IDT and other vendors accept IUPAC nomenclature. Provide CLI flag `--resolve-spacer AA` for users who want explicit bases |
| Mixed-case initiator sequences in R script (lowercase = RNA-binding in output convention) | Store initiator sequences with original case in AMPLIFIER_TABLE for reference; apply UPPERCASE case convention to spacer+initiator in output (binding region lowercase) |
| Homopolymer filter must consider the full synthesised oligo (binding + spacer + initiator), not just the 25-nt binding region | Apply synthesis-quality filters to the complete oligo after initiator attachment |
| Very short transcripts may yield zero valid pairs even with asymmetric leniency | Report clearly; suggest relaxing ΔG bounds or reducing requested pair count |
| Genome/pseudogene mask levels: existing bowtie produces nucleotide-level masks on the full sequence | Reuse the existing single bowtie run on the full concatenated sequence (same as smFISH). Apply the resulting nucleotide-level mask **per-half**: for a pair at position `p`, check `mask[p:p+25]` (left) and `mask[p+27:p+52]` (right). Gap positions `p+25:p+27` are skipped. No separate bowtie runs per half needed |
| The 2-nt gap nucleotides are part of the target but not hybridised — they must not be masked (they are "free") | Masking logic and `mask_to_badness` conversion explicitly skip gap positions; only the two 25-nt binding windows contribute to pair mask evaluation |
| **Junction marker `>` in gap not caught by `calculate_badness`** | `has_invalid_chars()` rejects 25-mers containing `>`, but if `>` falls at position `p+25` or `p+26` (the gap), neither 25-nt half sees it — the pair appears valid but actually spans two concatenated sequences. **Must add explicit check**: `if '>' in seq[p+25:p+27]: pair_badness = inf` |
| N's in gap positions | N's at `p+25`/`p+26` do NOT hybridise and should NOT invalidate the pair. The R-mask marks these positions, but the per-half masking logic skips the gap, so pairs with N-only-in-gap are correctly retained |

---

## 2. Architecture

### 2.1 New files

| File | Purpose |
|-|-|
| `src/probedesign/hcr.py` | HCR-specific logic: pair badness, pair DP, initiator attachment, HCR constants |
| `src/probedesign/hcr_output.py` | HCR-specific output formatting (oligos.txt, seq.txt, hits summary) |
| `streamlit_app/pages/smfish.py` | smFISH probe design page (refactored from current `app.py`) |
| `streamlit_app/pages/hcr.py` | HCR probe design page |
| `tests/test_hcr.py` | HCR unit tests with simulated FASTA sequences (written BEFORE implementation) |

### 2.2 Modified files

| File | Change |
|-|-|
| `cli.py` | New `hcrdesign` subcommand |
| `app.py` | Refactored into multipage entry point with `st.navigation`; existing smFISH logic moves to `pages/smfish.py` |
| `utils.py` | `run_hcr_design()` wrapper function added |

### 2.3 Unchanged files

`core.py`, `thermodynamics.py`, `sequence.py`, `fasta.py`, `masking.py` — all reused as-is. The thermodynamic functions (`gibbs_rna_dna`, `tm_rna_dna`) already work on arbitrary-length sequences. `calculate_badness` from `core.py` is reused at 25-nt length. Masking functions (`pseudogene_mask`, `genome_mask`, `mask_to_badness`) produce nucleotide-level masks on the full sequence; HCR applies these per-half in `hcr.py` rather than adding new functions to `masking.py`.

---

## 3. Data Structures

```python
@dataclass
class HCRProbeHalf:
    """One half of an HCR split-initiator probe pair."""
    position: int          # Start position on target (0-indexed)
    length: int            # Always 25
    binding_seq: str       # 25-nt target-binding sequence (sense strand)
    binding_rc: str        # Reverse complement (what is synthesised)
    gc_percent: float
    tm: float
    gibbs_fe: float
    is_strict: bool        # True if this half passed strict filtering
    oligo_seq: str          # Full synthesised oligo (binding_rc + spacer + initiator)

@dataclass
class HCRProbePair:
    """A complete HCR split-initiator probe pair."""
    pair_index: int         # 1-based pair number
    left: HCRProbeHalf      # Binds 5' portion of target window (P2 in output — carries initiator_a)
    right: HCRProbeHalf     # Binds 3' portion of target window (P1 in output — carries initiator_b)
    pair_position: int      # Start of the 52-nt window on target
    combined_badness: float
    amplifier: str          # e.g. "B1"

@dataclass
class HCRDesignResult:
    """Result of HCR probe design."""
    pairs: List[HCRProbePair]
    score: float
    input_sequence: str
    template_name: str
    amplifier: str
    mask_left: Optional[List[int]]   # Per-position mask for left halves
    mask_right: Optional[List[int]]  # Per-position mask for right halves
    mask_strings: List[str]
    bowtie_hits_left: Optional[dict]
    bowtie_hits_right: Optional[dict]
    asymmetric_mode: str              # "symmetric", "asymmetric_gibbs", "asymmetric_bowtie", "asymmetric_both"
```

---

## 4. Algorithm Detail

### 4.1 Badness calculation for HCR pairs

**Step 0 — Pre-compute 25-nt badness arrays.** Reuse `calculate_badness()` from `core.py` with `oligo_length=25`. This function returns one float per start position, where each entry evaluates the 25-mer at `seq[i:i+25]`. Any 25-mer containing invalid characters (including `>` junction markers or `n`) gets `inf` badness automatically via `has_invalid_chars()`.

Compute TWO arrays:
- `badness_strict[pos]` using `target_gibbs=-31.0`, `allowable_range=(-35.0, -27.0)`
- `badness_lenient[pos]` using `target_gibbs=-31.0`, `allowable_range=(-42.0, -27.0)` (lenient floor is user-configurable via `--lenient-gibbs-min`)

In symmetric mode, only `badness_strict` is needed. In any asymmetric mode, both arrays are computed.

**Step 1 — Pair-level badness.** For each candidate pair starting at position `p` (valid range: `0` to `len(seq) - 52`):
- Left half: positions `[p, p+25)` → `badness_strict[p]` or `badness_lenient[p]`
- Gap: positions `p+25` and `p+26` — not evaluated for ΔG
- Right half: positions `[p+27, p+52)` → `badness_strict[p+27]` or `badness_lenient[p+27]`

**Step 2 — Junction marker check.** Before computing pair badness, check for junction markers in the gap:
```
if '>' in seq[p+25:p+27]:
    pair_badness[p] = inf   # pair spans two concatenated sequences
    continue
```
This catches a case that `calculate_badness()` misses: `>` at gap positions `p+25` or `p+26` is invisible to both 25-nt halves but means the pair straddles two separate sequences.

**Step 3 — Combine half badness into pair badness.**

**Symmetric mode:**
```
pair_badness[p] = badness_strict[p] + badness_strict[p+27]
```
Both halves must have finite strict badness.

**Asymmetric mode (OR-mode):**
```
option_A = badness_strict[p] + badness_lenient[p+27]    # left strict, right lenient
option_B = badness_lenient[p] + badness_strict[p+27]    # right strict, left lenient
pair_badness[p] = min(option_A, option_B)
```
At least one half must pass strict filtering. The other half must pass lenient filtering. Both halves lenient (`badness_strict` is `inf` for both) → `inf`. Record which orientation was chosen for each position.

### 4.2 Masking

**Key concept: nucleotide-level masks → pair-level application.** All mask sources (R, P, B) produce nucleotide-level binary arrays of length `len(seq)`, exactly as in smFISH. The HCR-specific step is converting these to pair-level decisions while skipping gap positions. No new masking functions are needed in `masking.py`.

**Mask-to-pair conversion (in `hcr.py`):**
For a pair at position `p`:
```python
left_masked  = any(nuc_mask[p : p+25])
right_masked = any(nuc_mask[p+27 : p+52])
# nuc_mask[p+25] and nuc_mask[p+26] are IGNORED (gap)
```

**Which masks are always symmetric vs potentially asymmetric:**

| Mask | Symmetric always? | Rationale |
|-|-|-|
| R (repeats) | Yes | Repeat regions indicate unreliable sequence — neither half should bind there |
| L (low-complexity) | Yes | Homopolymer/dinucleotide affects synthesis quality regardless of specificity role |
| F (thermodynamic) | No | Handled by the strict/lenient badness arrays in §4.1, not a separate mask step |
| P (pseudogene) | Potentially asymmetric | Off-target specificity concern — safe to relax for the lenient half |
| B (genome) | Potentially asymmetric | Same rationale as P mask |

**R mask (repeats):** Reuse existing nucleotide-level repeat mask (from N's in input or RepeatMasker file). A pair is masked if `left_masked OR right_masked`. Gap positions are not checked. This is unconditional — both halves must be repeat-free.

**P mask (pseudogene) and B mask (genome):** Run the existing `pseudogene_mask()` and `genome_mask()` functions on the full concatenated sequence — exactly one bowtie invocation per mask type, same as smFISH. These return nucleotide-level masks. Then apply per-half:

- **Symmetric bowtie mode:** `if left_masked or right_masked: pair_badness = inf`. Both halves must be off-target-free.
- **Asymmetric bowtie mode:** A pair survives if at least one half is bowtie-clean (`not left_masked or not right_masked`). The other half's bowtie hits are **recorded for reporting** but do not mask. This check is independent of which half is ΔG-strict.

**Combined evaluation — ΔG and bowtie are independent filters:**

When both `--asymmetric-gibbs` and `--asymmetric-bowtie` are enabled, the two leniency dimensions are evaluated **independently**. The ΔG filter picks its own best orientation; the bowtie filter picks its own. They do not need to agree on which half is strict. Rationale: a single half cannot nucleate HCR regardless of how stably it binds off-target, so specificity (bowtie) and binding strength (ΔG) are independent safety concerns — each needs at least one good half, but not necessarily the same one.

```python
# For each pair position p:
# Pre-computed: badness_strict[p], badness_lenient[p], badness_strict[p+27], badness_lenient[p+27]
# Pre-computed: left_offmasked = any(pb_mask[p:p+25]), right_offmasked = any(pb_mask[p+27:p+52])

# R-mask is always symmetric — if either half is in a repeat, pair is dead
if any(r_mask[p:p+25]) or any(r_mask[p+27:p+52]):
    pair_badness = inf; continue

# Junction in gap
if '>' in seq[p+25:p+27]:
    pair_badness = inf; continue

# ── ΔG filter (independent) ──
if asymmetric_gibbs:
    gibbs_A = badness_strict[p] + badness_lenient[p+27]      # left ΔG-strict
    gibbs_B = badness_lenient[p] + badness_strict[p+27]      # right ΔG-strict
    gibbs_score = min(gibbs_A, gibbs_B)                       # best orientation
else:
    gibbs_score = badness_strict[p] + badness_strict[p+27]    # both strict

if not math.isfinite(gibbs_score):
    pair_badness = inf; continue

# ── Bowtie filter (independent) ──
if asymmetric_bowtie:
    bowtie_ok = (not left_offmasked) or (not right_offmasked) # at least one half clean
else:
    bowtie_ok = (not left_offmasked) and (not right_offmasked) # both must be clean

if not bowtie_ok:
    pair_badness = inf; continue

pair_badness = gibbs_score
```

This logic handles all four mode combinations:
1. **Symmetric** (no flags): both halves strict ΔG, both halves strict bowtie
2. **Asymmetric ΔG only**: best ΔG orientation (one strict, one lenient); both halves must be bowtie-clean
3. **Asymmetric bowtie only**: both halves strict ΔG; at least one half bowtie-clean
4. **Asymmetric both**: best ΔG orientation (independent); at least one half bowtie-clean (independent) — the ΔG-strict half and the bowtie-clean half need not be the same

**L mask (low-complexity):** Applied to each 25-nt half independently, **always symmetric** (both halves checked regardless of asymmetric mode). Uses `has_homopolymer()` and `has_dinucleotide_repeat()` from `masking.py` on the binding sequence only. If a half fails, set its badness to `inf` in both strict and lenient arrays before pair evaluation.

**F mask (thermodynamic):** Implicit in the badness calculation — positions with `inf` badness are already filtered. For visualization in `_HCR_seq.txt`, a position shows 'F' only if both `badness_strict[p]` and `badness_lenient[p]` are `inf`.

### 4.3 Dynamic programming

**State:** `dp[p][k]` = minimum total pair badness achievable by placing exactly `k` non-overlapping pairs, with the last pair starting at or before position `p`.

**Transition:** For pair `k` starting at position `p`, the previous pair must end before `p - pair_spacing`. Since each pair occupies 52 nt, the previous pair's start must be ≤ `p - 52 - pair_spacing`.

```
dp[p][k] = min(
    dp[p-1][k],                                           # skip position p
    dp[p - 52 - pair_spacing][k-1] + pair_badness[p]      # place pair at p
)
```

**Objective:** Maximise pair count first (find largest `k` such that `dp[N][k] < inf`), then minimise total badness among all solutions with that `k`. This matches the existing smFISH approach.

**Backtracking:** Standard DP backtrack to recover pair positions.

**Computational cost:** O(N × K) where N = sequence length, K = max pairs. For a typical transcript (≤5000 nt) and ≤50 pairs, this is trivial (<1ms). No performance concerns.

### 4.4 Initiator attachment and orientation

After DP selects pair positions, attach initiator sequences using the amplifier lookup table.

**Probe structure (initiators point INWARD toward the 2-nt gap):**

```
P1 (odd probes, binds RIGHT half):  5'—[RC(right)]—[spacer_b]—[initiator_b]—3'
P2 (even probes, binds LEFT half):  5'—[initiator_a]—[spacer_a]—[RC(left)]—3'
```

Where `RC(left)` and `RC(right)` are the reverse complements of the left and right 25-nt
target-binding regions respectively (since probes are antisense DNA).

**Why initiators must point inward:**

In HCR v3.0, the two half-initiators on adjacent probes must be spatially proximal to
cooperatively nucleate hairpin polymerisation. When both probes hybridise to the target,
their initiator tails must converge at the 2-nt gap between binding sites. If initiators
pointed outward, they would be ~50 nt apart and unable to co-trigger HCR.


**Validated against reference R script output (Dmel cip4-exon, amplifier B1):**

Reference pair id_618 (target position 618-669) and our pair 2 (target position 620-671)
bind nearly the same target region (2 nt offset). Both produce structurally identical probes:

```
Reference (R script):
  P1 (01): ggacagatcctcgaagggtatgtccTAgAAgAgTCTTCCTTTACg
           |-----RC(right) 25nt----||--spacer_b+init_b--|
           binding (lowercase)      tail (mixed case)

  P2 (02): gAggAgggCAgCAAACggAAtggtggcgtgaaaccagattgatat
           |--init_a+spacer_a--||------RC(left) 25nt----|
           head (mixed case)    binding (lowercase)

Our output (Python):
  P1 (01): ttggacagatcctcgaagggtatgtTAGAAGAGTCTTCCTTTACG
           |-----RC(right) 25nt----||--spacer_b+init_b--|
           binding (lowercase)      tail (UPPERCASE)

  P2 (02): GAGGAGGGCAGCAAACGGAActtggtggcgtgaaaccagattgat
           |--init_a+spacer_a--||------RC(left) 25nt----|
           head (UPPERCASE)     binding (lowercase)

Structure comparison:
  P1 structure:  5'--binding_rc--spacer_b--init_b--3'     MATCH
  P2 structure:  5'--init_a--spacer_a--binding_rc--3'     MATCH
  init_b seq:    GAAGAGTCTTCCTTTACG (B1)                  MATCH
  init_a seq:    GAGGAGGGCAGCAAACGG (B1)                  MATCH
  spacer_b:      TA (B1)                                  MATCH
  spacer_a:      AA (B1)                                  MATCH
  Orientation:   Both initiators point INWARD              MATCH
  Case:          Our output uses UPPERCASE for initiators  (design choice)

 5'-IIIIIIIIIIIIIII  IIIIIIIIIIIIIII-3'       
                  S  S
                  S  S
       3'-NNNNNNNNN  NNNNNNNNNN-5'
5'-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-3' target mRNA
```

**WW spacer handling:** For amplifiers B7-B17, the R script uses `WW` (IUPAC: A or T). Implementation:
- **Default: output `WW` directly** in the oligo sequence. Oligo vendors (IDT, Sigma, etc.) routinely accept IUPAC ambiguity codes and synthesise a degenerate pool.
- CLI flag `--resolve-spacer` to override with explicit bases (e.g. `--resolve-spacer AA` or `--resolve-spacer TT`)
- When WW is output, the full oligo length comment in the output file notes "contains degenerate bases"

**Folded probe architecture — actual sequences from Dmel cip4-exon pair 2 (amplifier B1):**

When both P1 and P2 hybridise to adjacent target sites, the single-stranded spacer+initiator
tails fold away from the target and converge at the 2-nt gap. This spatial proximity enables
cooperative nucleation of HCR hairpin polymerisation.

```
================================================================================
         FOLDED HCR PROBE PAIR ON TARGET
         Dmel cip4-exon pair 2, amplifier B1
================================================================================

              HALF-INITIATORS CONVERGE AT GAP
              (single-stranded, not hybridised)
                          |
                          V
       P2 init_a (18nt)                   P1 init_b (18nt)
     5'-GAGGAGGGCAGCAAACGG                 GAAGAGTCTTCCTTTACG-3'
                         \                 /
                         A|               |T   <-- spacer_a (AA)
                         A|     2nt       |A   <-- spacer_b (TA)
                          |     gap       |
                          |      |        |
  3'-TAGTTAGACCAAAGTGCGGTGGTTC   |  TGTATGGGAAGCTCCTAGACAGGTT-5'  <-- PROBES
      |||||||||||||||||||||||||  |   |||||||||||||||||||||||||      (antisense)
  5'-ATCAATCTGGTTTCACGCCACCAAG--GG--ACATACCCTTCGAGGATCTGTCCAA-3'   <-- TARGET
      <------ LEFT_HALF ------>  gap  <------ RIGHT_HALF ------>
             (25 nt)                          (25 nt)
          pos 620-644           645-6       pos 647-671

--------------------------------------------------------------------------------
COMPONENT LEGEND
--------------------------------------------------------------------------------

TARGET mRNA (5'-->3'):
  LEFT_HALF  = ATCAATCTGGTTTCACGCCACCAAG   (25 nt, pos 620-644)
  gap        = GG                          (2 nt,  pos 645-646)
  RIGHT_HALF = ACATACCCTTCGAGGATCTGTCCAA   (25 nt, pos 647-671)

P2 OLIGO -- binds LEFT_HALF (read 5'-->3' as synthesised):
  +--------------------+----+---------------------------+
  | init_a (18 nt)     | AA | RC(LEFT_HALF) (25 nt)     |
  | GAGGAGGGCAGCAAACGG | AA | cttggtggcgtgaaaccagattgat |
  +--------------------+----+---------------------------+
  <-- 5' end                                3' end -->

P1 OLIGO -- binds RIGHT_HALF (read 5'-->3' as synthesised):
  +---------------------------+----+--------------------+
  | RC(RIGHT_HALF) (25 nt)    | TA | init_b (18 nt)     |
  | ttggacagatcctcgaagggtatgt | TA | GAAGAGTCTTCCTTTACG |
  +---------------------------+----+--------------------+
  <-- 5' end                                3' end -->

--------------------------------------------------------------------------------
HYBRIDISATION GEOMETRY
--------------------------------------------------------------------------------

  * P2's 3' end (binding region) aligns to LEFT_HALF's 5' end (anti-parallel)
  * P2's 5' end (init_a) protrudes INWARD toward the gap
  * P1's 5' end (binding region) aligns to RIGHT_HALF's 3' end (anti-parallel)
  * P1's 3' end (init_b) protrudes INWARD toward the gap
  * Both spacer+initiator tails are single-stranded and fold UP from the target
  * The two half-initiators converge above the 2-nt gap --> cooperative HCR

================================================================================
```

Full synthesised oligos (lowercase = target-binding, UPPERCASE = spacer + initiator):
```
P1 (HCRB1_01): 5'-ttggacagatcctcgaagggtatgtTAGAAGAGTCTTCCTTTACG-3'  (45 nt)
P2 (HCRB1_02): 5'-GAGGAGGGCAGCAAACGGAActtggtggcgtgaaaccagattgat-3'  (45 nt)
```

### 4.5 Post-DP validation

After probe pairs are selected:
1. Verify no pair overlaps (start positions differ by ≥ 52 + pair_spacing)
2. Apply homopolymer/dinucleotide filter to the full synthesised oligo (binding_rc + spacer + initiator) as a final check. If a full oligo fails, consider it a warning rather than a hard reject (since the initiator is fixed and cannot be changed)
3. Report per-pair summary: position, ΔG of each half, which half is strict/lenient, GC%, Tm

---

## 5. CLI Interface

### 5.1 New subcommand: `hcrdesign`

```
probedesign hcrdesign INPUT_FILE [OPTIONS]
```

### 5.2 Options

| Option | Type | Default | Description |
|-|-|-|-|
| `--pairs` / `-p` | int | 30 | Number of probe pairs to design |
| `--amplifier` / `-a` | choice | B1 | HCR amplifier (B1–B5, B7, B9–B11, B13–B15, B17) |
| `--target-gibbs` | float | -31.0 | Target ΔG for each 25-nt half (kcal/mol) |
| `--allowable-gibbs` | str | "-35,-27" | Strict ΔG range (min,max) |
| `--asymmetric-gibbs` | flag | off | Enable asymmetric ΔG leniency |
| `--lenient-gibbs-min` | float | -42.0 | Floor ΔG for lenient half — most negative allowed (only with `--asymmetric-gibbs`) |
| `--asymmetric-bowtie` | flag | off | Enable asymmetric bowtie filtering |
| `--pair-spacing` | int | 2 | Minimum gap (nt) between consecutive 52-nt pair blocks |
| `--species` / `-s` | choice | human | Species for masking |
| `--pseudogene-mask` | flag | off | Enable pseudogene bowtie masking |
| `--genome-mask` | flag | off | Enable genome bowtie masking |
| `--index-dir` | str | None | Path to bowtie index directory |
| `--repeatmask-file` | str | None | Pre-computed RepeatMasker output |
| `--hp-threshold` | int | 5 | Homopolymer run threshold |
| `--di-threshold` | int | 3 | Dinucleotide repeat threshold |
| `--resolve-spacer` | str | None | Resolve WW spacer to explicit bases (e.g. "AA", "TT"). Default: output WW as IUPAC |
| `--output-name` / `-o` | str | None | Output file prefix |

### 5.3 Example usage

```bash
# Symmetric design, 20 pairs, amplifier B1
probedesign hcrdesign input.fa --pairs 20 --amplifier B1 \
  --pseudogene-mask --genome-mask --index-dir bowtie_indexes

# Asymmetric leniency on both Gibbs and bowtie
probedesign hcrdesign input.fa --pairs 25 --amplifier B2 \
  --asymmetric-gibbs --lenient-gibbs-min -42 \
  --asymmetric-bowtie \
  --pseudogene-mask --genome-mask --index-dir bowtie_indexes

# Custom spacing, species mouse, resolve WW spacer
probedesign hcrdesign input.fa -p 30 -a B9 -s mouse \
  --pair-spacing 5 --target-gibbs -30 --resolve-spacer AA
```

---

## 6. Output Format

### 6.1 `_HCR_oligos.txt`

Tab-separated, one row per oligo (two rows per pair). Case convention: **lowercase** = RNA-binding sequence, **UPPERCASE** = spacer + initiator.

```
pair	half	start	GC%	Tm	Gibbs	strict	full_oligo	name
1	P1	69	52.0	62.1	-30.8	yes	acgtacgtacgtacgtacgtacgtaTAGAAGAGTCTTCCTTTACG	GENE_HCRB1_01
1	P2	42	48.0	59.3	-29.5	yes	GAGGAGGGCAGCAAACGGAAacgtacgtacgtacgtacgtacgta	GENE_HCRB1_02
2	P1	183	56.0	64.2	-32.1	yes	acgtacgtacgtacgtacgtacgtaTAGAAGAGTCTTCCTTTACG	GENE_HCRB1_03
2	P2	156	44.0	57.8	-28.3	no	GAGGAGGGCAGCAAACGGAAacgtacgtacgtacgtacgtacgta	GENE_HCRB1_04
```

Column definitions:
- `pair`: Pair number (1-based)
- `half`: P1 (binds right half, carries initiator_b at 3') or P2 (binds left half, carries initiator_a at 5')
- `start`: Start position of the 25-nt binding region on the target
- `GC%`, `Tm`, `Gibbs`: Thermodynamic properties of the 25-nt binding region only
- `strict`: Whether this half passed strict filtering (`yes`/`no`; `no` = lenient half in asymmetric mode)
- `full_oligo`: Complete synthesised sequence with case convention
- `name`: `{gene}_HCR{amplifier}_{NN}` where odd NN = P1, even NN = P2

### 6.2 `_HCR_seq.txt`

Visual alignment similar to smFISH `_seq.txt`, adapted for pairs:

```
Position: 1-110
SEQUENCE:  ATGCGATCGA...
R mask:    ..........R...
P mask:    ..........P...
B mask:    ...............
L mask:    ...............
Pair 1:    ..........[===LEFT===]--[===RIGHT===]..........
Pair 2:    ..............................[===LEFT===]--[===RIGHT===]...
```

Each pair is shown as a single line with left half, 2-nt gap (`--`), and right half. Position labels align with the sequence. The mask lines show the union of masks for both halves.

### 6.3 `_HCR_hits.txt`

Bowtie hit summary, split by half:

```
=== Pair 1, P1 (right, pos 69-93, STRICT) ===
Pseudogene hits: 0
Genome hits: chr1:12345 (1 mismatch), chr5:67890 (2 mismatches)

=== Pair 1, P2 (left, pos 42-66, STRICT) ===
Pseudogene hits: 0
Genome hits: none

=== Pair 2, P1 (right, pos 183-207, STRICT) ===
...
=== Pair 2, P2 (left, pos 156-180, LENIENT — hits not used for filtering) ===
Pseudogene hits: 1
Genome hits: chr3:11111 (0 mismatches)
```

This makes it transparent which half was lenient and that its hits were recorded but not used for masking.

### 6.4 Naming convention

Odd numbers = P1 (binds right half of target), even numbers = P2 (binds left half of target). This follows the R script convention: the full reverse complement of the 52-nt target block is split into first half (P1, odd) and second half (P2, even). P1 carries initiator_b and P2 carries initiator_a.

---

## 7. Streamlit Integration — Multipage Architecture

### 7.1 Navigation structure

Convert the app from a single-page layout to a multipage architecture using Streamlit's `st.Page` + `st.navigation` (available since Streamlit 1.36; installed version is 1.54).

**Top navigation:** `st.navigation(position=...)` only supports `"sidebar"` and `"hidden"` natively — there is no `position="top"`. To achieve top-bar navigation:

1. Use `st.navigation(position="hidden")` to handle page routing without the default sidebar nav.
2. Render a custom top navigation bar using `st.columns` + `st.button` (or `st.pills` if available) at the top of each page.
3. Use `st.switch_page()` to navigate between pages programmatically.

**File structure:**
```
streamlit_app/
├── app.py              # Entry point: st.navigation + top nav bar
├── pages/
│   ├── smfish.py       # smFISH probe design (refactored from current app.py)
│   └── hcr.py          # HCR probe design
├── utils.py            # Shared helpers (existing + run_hcr_design)
└── README.md
```

**`app.py` (entry point):**
```python
import streamlit as st

st.set_page_config(page_title="ProbeDesign", layout="wide")

smfish_page = st.Page("pages/smfish.py", title="smFISH Probe Design")
hcr_page = st.Page("pages/hcr.py", title="HCR Probe Design")

pg = st.navigation([smfish_page, hcr_page], position="hidden")

# Top navigation bar
nav_cols = st.columns([1, 1, 6])
with nav_cols[0]:
    if st.button("smFISH", use_container_width=True,
                 type="primary" if pg.title == "smFISH Probe Design" else "secondary"):
        st.switch_page(smfish_page)
with nav_cols[1]:
    if st.button("HCR", use_container_width=True,
                 type="primary" if pg.title == "HCR Probe Design" else "secondary"):
        st.switch_page(hcr_page)
st.divider()

pg.run()
```

### 7.2 smFISH page (`pages/smfish.py`)

Move the entire current `app.py` logic (sidebar parameters, single/batch modes, result display) into this file. No functional changes — purely a refactor.

### 7.3 HCR page (`pages/hcr.py`)

**Sidebar layout** (mirroring smFISH structure, parameters in logical order):

```
┌─ Sidebar ────────────────────────┐
│ HCR Probe Design                 │
│                                  │
│ Mode: [Single] [Batch]           │
│                                  │
│ Species: [human ▼]               │
│ Amplifier: [B1 ▼]               │
│ Number of pairs: [30]            │
│ Pair spacing (nt): [2]           │
│                                  │
│ ── Thermodynamic parameters ──   │
│ Target ΔG per half: [-31.0]      │
│ Strict ΔG min: [-35.0]          │
│ Strict ΔG max: [-27.0]          │
│                                  │
│ [ ] Asymmetric ΔG leniency       │
│   └─ Lenient ΔG floor: [-42.0]  │  ← only visible when toggled on
│                                  │
│ ── Low-complexity filter ──      │
│ Homopolymer: [5]  Dinucs: [3]    │
│                                  │
│ ── Off-target filter ──          │
│ [✓] Pseudogene mask              │
│ [✓] Genome mask                  │
│ [ ] Asymmetric bowtie            │  ← only visible when masks are on
│ Index dir: [bowtie_indexes]      │
│ [ ] Save raw bowtie              │
│                                  │
│ ── RepeatMasker (advanced) ──    │
│ [expander]                       │
└──────────────────────────────────┘
```

**Key UI decisions:**
- `Lenient ΔG floor` only appears when `Asymmetric ΔG leniency` is toggled on
- `Asymmetric bowtie` only appears when pseudogene or genome mask is enabled
- Amplifier dropdown includes B1–B5 (standard) and B7+ (labelled "experimental" — WW spacer)
- Single/batch modes mirror smFISH (upload, paste, directory path, batch zip download)

### 7.4 Backend changes in `utils.py`

New functions:
- `run_hcr_design()` — wraps the HCR pipeline, analogous to `run_design()`. Accepts HCR-specific parameters, calls `design_hcr_probes()`, captures stdout.
- `run_hcr_batch()` — iterates FASTA files calling `run_hcr_design()` per file, writes outputs, tracks progress.

### 7.5 HCR result display

- Probe pair table with columns: pair, half (P1/P2), start, GC%, Tm, Gibbs, strict/lenient, full oligo
- Download buttons for `_HCR_oligos.txt`, `_HCR_seq.txt`, `_HCR_hits.txt`
- Sequence alignment viewer showing pair positions with left/gap/right blocks
- Batch mode: summary table + per-file expandable details + ZIP download

---

## 8. Implementation Steps — Test-Gated Development

**Gate rule:** Each phase's tests must pass before proceeding to the next phase. If tests fail, fix the implementation before moving on. Tests are written FIRST for each phase, then the implementation is built to satisfy them.

### Phase 0: Test scaffolding (`tests/test_hcr.py`)

0. Create simulated FASTA test fixtures (see §9 for exact sequences)
1. Write all unit tests for Phase 1 functions as stubs — they should all FAIL initially
2. Verify test runner discovers and runs the test file

**Gate:** Test file runs, all tests discovered (expected: all fail).

### Phase 1: Core HCR module (`hcr.py`)

1. Define data structures (`HCRProbeHalf`, `HCRProbePair`, `HCRDesignResult`)
2. Define initiator lookup table (all amplifiers B1–B17 with sequences and spacers)
3. Implement `calculate_pair_badness()`:
   - Call existing `calculate_badness()` with `oligo_length=25` for strict range → `badness_strict`
   - Call again with lenient range → `badness_lenient` (only if asymmetric mode)
   - Apply L-mask to both arrays (homopolymer/dinucleotide on each 25-nt half)
   - Check for junction markers in gap positions
   - Combine into per-pair-position badness array using symmetric/asymmetric logic from §4.2
4. Implement `apply_hcr_masks()`:
   - Accept nucleotide-level masks from existing R/P/B mask functions (no new bowtie runs)
   - Convert to pair-level: check each 25-nt half window, skip gap positions `[p+25:p+27]`
   - Apply R-mask symmetrically; apply P/B masks symmetrically or asymmetrically per mode
5. Implement `find_best_pairs()` (DP):
   - State: pair starting positions, pair count
   - Transition: enforce 52-nt block + pair_spacing gap
   - Backtrack to recover positions
6. Implement `attach_initiators()`:
   - Reverse complement binding sequences
   - Attach spacer + initiator per the P1/P2 convention
   - Output WW spacer as IUPAC by default; resolve only if `--resolve-spacer` provided
7. Implement main entry point `design_hcr_probes()`:
   - Orchestrates: FASTA read → badness → masking → DP → initiator attachment → result

**Gate:** All Phase 1 unit tests pass (see §9.1 test cases).

### Phase 2: Output module (`hcr_output.py`)

8. Implement `format_hcr_oligos()` → `_HCR_oligos.txt` content
9. Implement `format_hcr_seq()` → `_HCR_seq.txt` content
10. Implement `format_hcr_hits()` → `_HCR_hits.txt` content
11. Implement `write_hcr_output_files()` that writes all three files

**Gate:** Output formatting tests pass.

### Phase 3: CLI integration

12. Add `hcrdesign` subcommand to `cli.py` with all options from §5.2
13. Wire subcommand to `design_hcr_probes()` and `write_hcr_output_files()`

**Gate:** CLI integration tests pass (end-to-end with simulated FASTA).

### Phase 4: Streamlit multipage refactor

14. Refactor `app.py` into multipage entry point with `st.navigation`
15. Move existing smFISH logic to `pages/smfish.py` — verify smFISH still works
16. Create `pages/hcr.py` with HCR-specific sidebar and main panel

**Gate:** smFISH regression test (existing functionality unbroken). HCR page renders with all controls.

### Phase 5: Streamlit HCR wiring

17. Add `run_hcr_design()` and `run_hcr_batch()` to `utils.py`
18. Wire HCR page controls to backend, add result display and downloads

**Gate:** End-to-end Streamlit HCR test with simulated input produces correct output.

### Phase 6: Integration testing

19. Run full integration tests (see §9.2)
20. Validate against existing smFISH pipeline (regression tests)

**Gate:** All integration tests pass. smFISH produces identical results to before.

---

## 9. Test Cases — Simulated Sequences and Expected Outcomes

All tests use simulated FASTA sequences with pre-computed thermodynamic properties. Tests call `gibbs_rna_dna()` at setup time to compute exact ΔG values for each half, then assert pipeline behavior based on those values. This avoids hardcoding ΔG numbers that could drift if thermodynamic parameters change.

### 9.0 Test fixtures — simulated FASTA sequences

Each fixture is a minimal sequence designed to trigger a specific condition. Approximate ΔG values shown in comments are for orientation only; tests compute exact values at runtime.

```python
# ── Fixture 1: BOTH_STRICT ──────────────────────────────────────────────────
# Both 25-nt halves have balanced GC → ΔG ≈ -29 to -31, within strict range (-35, -27)
# Expected: valid pair in symmetric mode
BOTH_STRICT = ">BOTH_STRICT\natcgatcgatcgatcgatcgatcgaatatcgatcgatcgatcgatcgatcga"
# Left:  atcgatcgatcgatcgatcgatcga  (positions 0–24)
# Gap:   at                          (positions 25–26)
# Right: atcgatcgatcgatcgatcgatcga  (positions 27–51)

# ── Fixture 2: ASYM_RESCUE ──────────────────────────────────────────────────
# Left half: balanced GC → ΔG ≈ -29 (strict pass)
# Right half: GC-rich → ΔG ≈ -38 (strict fail, lenient pass at -42 floor)
# Expected: invalid in symmetric mode, valid in asymmetric mode (left=strict, right=lenient)
ASYM_RESCUE = ">ASYM_RESCUE\natcgatcgatcgatcgatcgatcgaaagcgcgatcgcgatcgcgatcgcgat"

# ── Fixture 3: BOTH_LENIENT_FAIL ─────────────────────────────────────────────
# Both halves GC-rich → ΔG ≈ -38 (both strict fail)
# Expected: invalid in ALL modes (both-lenient is forbidden)
BOTH_LENIENT_FAIL = ">BOTH_LENIENT_FAIL\ngcgcgatcgcgatcgcgatcgcgatcgcgcgatcgcgatcgcgatcgcgatc"

# ── Fixture 4: HOMOPOLYMER_LEFT ──────────────────────────────────────────────
# Left half contains AAAAA (homopolymer run at threshold 5)
# Right half is clean balanced GC
# Expected: pair invalid due to L-mask on left half (regardless of asymmetric mode)
HOMOPOLYMER_LEFT = ">HOMOPOLYMER_LEFT\natcgaaaaatcgatcgatcgatcgaatatcgatcgatcgatcgatcgatcga"

# ── Fixture 5: JUNCTION_IN_GAP ───────────────────────────────────────────────
# Two FASTA entries where the junction marker '>' falls at gap position 25
# After concatenation: first_entry(25nt) + '>' + second_entry(26nt+)
# Expected: pair at position 0 is INVALID (spans two sequences)
JUNCTION_GAP_SEQ1 = ">SEQ1\natcgatcgatcgatcgatcgatcga"   # exactly 25 nt
JUNCTION_GAP_SEQ2 = ">SEQ2\naatcgatcgatcgatcgatcgatcgatcg"  # 29 nt
# Concatenated: atcgatcgatcgatcgatcgatcga>aatcgatcgatcgatcgatcgatcgatcg
# Position 25 = '>' (gap) → pair at p=0 must be caught by junction check

# ── Fixture 6: N_IN_GAP ─────────────────────────────────────────────────────
# N's at positions 25-26 (gap). Both halves are clean balanced GC.
# Expected: pair is VALID (N's in gap do not affect binding)
N_IN_GAP = ">N_IN_GAP\natcgatcgatcgatcgatcgatcgannatcgatcgatcgatcgatcgatcga"

# ── Fixture 7: TOO_SHORT ────────────────────────────────────────────────────
# 51 nucleotides — below the 52-nt minimum for a single pair
# Expected: error message "Sequence too short for HCR probe design (minimum 52 nt)"
TOO_SHORT = ">TOO_SHORT\natcgatcgatcgatcgatcgatcgaatatcgatcgatcgatcgatcgatcg"

# ── Fixture 8: MULTI_PAIR ───────────────────────────────────────────────────
# 160-nt sequence with balanced GC throughout, enough for 2 pairs with spacing=2
# Expected: DP finds 2 pairs, starts differ by ≥ 54 (52+2)
MULTI_PAIR = ">MULTI_PAIR\n" + "atcgatcgatcgatcgatcgatcga" * 6 + "atcgatcgat"

# ── Fixture 9: MASK_ONE_HALF ────────────────────────────────────────────────
# Used with a synthetic nucleotide-level mask that covers only the left half.
# Tests that symmetric mode rejects but asymmetric bowtie mode accepts.
# Sequence itself is balanced GC (both halves pass strict ΔG).
MASK_ONE_HALF = ">MASK_ONE_HALF\natcgatcgatcgatcgatcgatcgaatatcgatcgatcgatcgatcgatcga"
# Synthetic mask: positions 0-24 marked as 1 (masked), positions 25-51 marked as 0
# Symmetric: pair invalid (left masked)
# Asymmetric bowtie, orientation B (right=strict): pair valid (right is clean)

# ── Fixture 10: DINUCLEOTIDE_RIGHT ───────────────────────────────────────────
# Right half contains ATATAT (dinucleotide repeat at threshold 3)
# Left half is clean
# Expected: pair invalid due to L-mask on right half
DINUCLEOTIDE_RIGHT = ">DINUCLEOTIDE_RIGHT\natcgatcgatcgatcgatcgatcgaaatatatatcgatcgatcgatcgatcg"
```

### 9.1 Unit tests (Phase 1 gate)

Each test below references a specific fixture and states the exact expected outcome. Tests that depend on ΔG values compute them with `gibbs_rna_dna()` in setup.

| Test | Fixture | Procedure | Expected outcome |
|-|-|-|-|
| `test_badness_strict_both_pass` | BOTH_STRICT | Compute `badness_strict` for left (pos 0) and right (pos 27). | Both finite. Sum is finite. |
| `test_badness_strict_gc_rich_fail` | ASYM_RESCUE | Compute `badness_strict` for right half. | `inf` (ΔG outside strict range). `badness_lenient` for same half is finite. |
| `test_pair_badness_symmetric_pass` | BOTH_STRICT | Call `calculate_pair_badness(symmetric=True)` at position 0. | Finite pair badness. |
| `test_pair_badness_symmetric_fail` | ASYM_RESCUE | Call `calculate_pair_badness(symmetric=True)` at position 0. | `inf` (right half fails strict). |
| `test_pair_badness_asymmetric_rescue` | ASYM_RESCUE | Call `calculate_pair_badness(asymmetric_gibbs=True)` at position 0. | Finite pair badness. Orientation recorded as "left=strict". |
| `test_pair_badness_both_lenient_rejected` | BOTH_LENIENT_FAIL | Call `calculate_pair_badness(asymmetric_gibbs=True)` at position 0. | `inf` — both halves fail strict, so no valid orientation exists. |
| `test_pair_badness_best_orientation` | ASYM_RESCUE | Call with asymmetric mode and verify that the chosen orientation has the lower combined badness. | `min(option_A, option_B)` is selected; the "strict" label matches the half with finite strict badness. |
| `test_junction_in_gap` | JUNCTION_IN_GAP | Concatenate SEQ1+SEQ2 with `>`. Evaluate pair at position 0. | `inf` — junction marker at gap position 25 detected. |
| `test_n_in_gap_valid` | N_IN_GAP | Evaluate pair at position 0. | Finite pair badness — N's at gap positions 25-26 are ignored. |
| `test_sequence_too_short` | TOO_SHORT | Call `design_hcr_probes()`. | Raises error or returns 0 pairs with message containing "too short". |
| `test_homopolymer_masks_pair` | HOMOPOLYMER_LEFT | L-mask check on left half `seq[0:25]`. | `has_homopolymer(seq[0:25], 5)` returns True. Pair badness = `inf`. |
| `test_dinucleotide_masks_pair` | DINUCLEOTIDE_RIGHT | L-mask check on right half `seq[27:52]`. | `has_dinucleotide_repeat(seq[27:52], 3)` returns True. Pair badness = `inf`. |
| `test_l_mask_always_symmetric` | HOMOPOLYMER_LEFT | Call with `asymmetric_gibbs=True`. | Pair still `inf` — L-mask ignores asymmetric mode. |
| `test_mask_symmetric_rejects` | MASK_ONE_HALF | Apply synthetic mask `[1]*25 + [0]*27`. Symmetric bowtie mode. | Pair `inf` (left half masked). |
| `test_mask_asymmetric_rescues` | MASK_ONE_HALF | Apply synthetic mask `[1]*25 + [0]*27`. Asymmetric bowtie mode. | Pair finite (orientation B: right=strict, right is unmasked). |
| `test_mask_both_masked_rejected` | MASK_ONE_HALF | Apply synthetic mask `[1]*52`. Asymmetric bowtie mode. | Pair `inf` — no clean half exists for strict role. |
| `test_bowtie_mask_spanning_gap` | BOTH_STRICT | Apply synthetic P/B mask `[0]*20 + [1]*12 + [0]*20` (simulates a 12-mer hit at position 20, covering positions 20-31 — bleeding from left half through gap into right half). Symmetric bowtie mode. | Pair `inf` — both halves have masked nucleotides (left: positions 20-24, right: positions 27-31). Gap positions 25-26 are also marked but irrelevant. This is the realistic scenario: no bowtie mask can mark only the gap without also marking a half. |
| `test_r_mask_n_in_gap_ignored` | N_IN_GAP | Convert N's at positions 25-26 to R-mask `[0]*25 + [1]*2 + [0]*25`. Evaluate pair. | Pair finite — R-mask at gap positions is ignored. Unlike bowtie masks, R-masks CAN mark only the gap (from N's in input sequence), making gap-skipping necessary for this mask type. |
| `test_r_mask_always_symmetric` | MASK_ONE_HALF | Apply R-mask (repeat type) to left half only. Asymmetric bowtie mode. | Pair `inf` — R-mask is always symmetric, even with asymmetric bowtie. |
| `test_dp_no_overlap` | MULTI_PAIR | Run DP requesting 2 pairs with spacing=2. | Two pairs found; `start_2 - start_1 >= 54`. |
| `test_dp_maximise_pairs` | MULTI_PAIR | DP with spacing=2, requesting 5 pairs. | Returns maximum possible pairs (≤ `len(seq) // 54`); prefers more pairs over lower badness. |
| `test_initiator_B1` | (any) | Call `attach_initiators(amplifier="B1")` with different left/right binding sequences. | P1 = `RC(right) + "TA" + initiator_b_B1` (binds right half, init_b at 3'). P2 = `initiator_a_B1 + "AA" + RC(left)` (binds left half, init_a at 5'). Both initiators point inward. |
| `test_initiator_all_amplifiers` | (any) | Call `attach_initiators()` for all B1–B17. | All produce oligos with correct structure. WW spacers output as "WW" by default. |
| `test_ww_spacer_default_iupac` | (any) | Call `attach_initiators(amplifier="B9")` without `--resolve-spacer`. | Spacer in output oligo is literally "WW" (IUPAC). |
| `test_ww_spacer_resolved` | (any) | Call `attach_initiators(amplifier="B9", resolve_spacer="AA")`. | Spacer in output oligo is "AA". |
| `test_output_case_convention` | (any) | Format a probe pair into oligos.txt row. | Binding region is lowercase, spacer+initiator is uppercase. |
| `test_output_naming` | (any) | Format probe pairs with gene name "GENE" and amplifier "B1". | Oligo 1 = `GENE_HCRB1_01` (P1), oligo 2 = `GENE_HCRB1_02` (P2). Odd=P1, even=P2. |

### 9.2 Integration tests (Phase 6 gate)

| Test | Input | Expected outcome |
|-|-|-|
| `test_cdkn1a_hcr_symmetric` | `CDKN1A.fa`, 15 pairs, B1, symmetric, repeatmask file | ≥10 valid pairs; all pair ΔG within (-35, -27) for both halves |
| `test_cdkn1a_hcr_asymmetric` | Same input, asymmetric ΔG enabled | ≥ pairs found than symmetric (demonstrates leniency benefit); strict halves within (-35, -27), lenient halves within (-42, -27) |
| `test_short_transcript` | Synthetic 200-nt balanced sequence | Symmetric: 0-2 pairs; asymmetric mode: potentially more pairs (demonstrates rescue) |
| `test_gc_extreme_sequence` | Synthetic high-GC 500-nt sequence | Asymmetric mode finds more valid pairs than symmetric mode |
| `test_smfish_regression` | All existing test cases from `test_cases/` | Existing smFISH `probedesign design` produces identical results (no regression from multipage refactor) |
| `test_streamlit_smfish_unchanged` | Load smFISH page, verify all sidebar widgets present | All existing smFISH controls render; no missing widgets or crashes |

### 9.3 Edge cases

| Case | Expected behaviour | How to test |
|-|-|-|
| Sequence shorter than 52 nt | Error: "Sequence too short for HCR probe design (minimum 52 nt)" | Fixture TOO_SHORT |
| Zero valid pairs (all pairs fail all filters) | Warning with diagnostic: report count of positions failing each filter type (F/R/P/B/L/junction) | High-GC sequence with all homopolymers |
| Sequence with N's only in gap region | Pair is valid (gap is not hybridised) | Fixture N_IN_GAP |
| Sequence with N's in binding region | Half is masked (existing `has_invalid_chars` logic) | Place N at position 10 in a 52-nt sequence |
| Junction marker `>` at position p+25 | Pair is invalid (caught by gap junction check) | Fixture JUNCTION_IN_GAP |
| Junction marker `>` at position p+26 | Pair is invalid (caught by gap junction check) | Shift SEQ1 to 26nt |
| Junction marker `>` within a 25-nt half | Half badness is `inf` (caught by `has_invalid_chars` in `calculate_badness`) | Place `>` at position 12 in left half |
| All pairs fail strict on both halves | Report "0 pairs found" with filter breakdown by mask type | All-GC sequence |
| `--pairs` exceeds maximum possible | Design maximum possible, warn user explicitly | Request 100 pairs from 160-nt sequence |
| Amplifier B9 with `--resolve-spacer TT` | WW replaced with TT in output | Verify output oligo contains "TT" not "WW" |
| Amplifier B1 with `--resolve-spacer AA` | No change (B1 already uses AA) | Verify no warning, output unchanged |
| Multi-FASTA input with >2 entries | Pairs cannot span junction markers; pairs on each segment are valid | 3x 100-nt entries concatenated |

### 9.4 Validation against R script

| Test | Method |
|-|-|
| `test_r_script_parity` | Run R script (`HCR_design.R`) and Python `hcrdesign` on the same input with equivalent parameters. Compare: (a) selected pair positions match, (b) oligo sequences match, (c) initiator attachment matches. Note: exact parity may not be achievable due to different DP tiebreaking, but thermodynamic calculations and filtering should agree. |

---

## 10. Computational Performance Notes

| Step | Cost | Notes |
|-|-|-|
| Badness calculation (25-nt, ×2 arrays) | O(N) | Two passes for strict + lenient; trivial |
| Pair badness combination | O(N) | Simple array operations per pair position |
| Bowtie masking | ~2-5s per mask type | **Single bowtie run on full sequence** — identical to smFISH. No per-half bowtie runs. Nucleotide-level mask is reused; per-half application is O(N). |
| DP | O(N × K) | Negligible for typical inputs |
| Initiator attachment | O(K) | Trivial |
| Output generation | O(N + K) | Trivial |

**Total expected runtime:** Comparable to smFISH design. Bowtie runs are identical (same sequence, same k-mers, same indices). The only additional cost is a second `calculate_badness()` call for the lenient array (O(N), negligible) and the per-pair mask evaluation (O(N), negligible). Total runtime remains under 30 seconds for typical transcripts.

**No significant new computational bottlenecks are introduced.**

---

## 11. Code Quality Requirements

- Production-style, maintainable code with type hints throughout
- Docstrings for all public functions
- No changes to existing smFISH code paths — HCR is a parallel module
- Reuse existing utilities (`thermodynamics.py`, `sequence.py`, `fasta.py`, `masking.py`) without modification
- Per-half mask application logic lives in `hcr.py`, not in `masking.py`
- All initiator data stored as a frozen dict constant in `hcr.py`, not in external files
- Error messages are specific and actionable (e.g. "Pair 7 failed: left half ΔG = -22.1 below minimum -27.0")
- CLI help text includes examples

---

## 12. Out of Scope (Future Enhancements)

These are explicitly **not** part of the initial implementation:

1. **Secondary structure prediction** for probe oligos (hairpin / self-dimer checks)
2. **Cross-pair off-target analysis** (checking if two lenient halves from different pairs could co-localise)
3. **Automatic amplifier selection** (choosing best amplifier based on initiator-target interactions)
4. **Multi-amplifier design** (designing probe sets for multiple channels simultaneously)
5. **Mixed-length HCR probes** (varying the 25-nt half length) — could be a future extension
6. **BLAST-based filtering** (the R script uses BLAST; this implementation uses bowtie for consistency with smFISH pipeline)

---

## 13. Summary of Key Design Decisions

| Decision | Rationale |
|-|-|
| Separate `hcrdesign` subcommand (not a flag on `design`) | HCR has fundamentally different parameters (pairs vs probes, amplifier, split-initiator); a separate subcommand is cleaner than conditional logic |
| Independent 25-nt ΔG calculation (not 52-nt) | Thermodynamically correct: two non-covalent hybridisation events |
| OR-mode asymmetric leniency | Maximises candidate pool while maintaining at least one specific half per pair |
| Single bowtie run on full sequence (not per-half) | Reuses existing nucleotide-level masking infrastructure unchanged; per-half application is O(N) post-processing in `hcr.py` |
| Nucleotide masks → pair-level with gap skipping | Gap positions `[p+25:p+27]` are unhybridised; masking them would incorrectly reject valid pairs, especially when N's fall in the gap |
| Explicit junction check in gap | `calculate_badness()` catches `>` within 25-nt halves but not at gap positions — dedicated check prevents cross-sequence pairs |
| DP on pair positions (not individual halves) | Pairs are atomic units; optimising halves independently could produce overlapping pairs |
| WW spacer output as IUPAC (not resolved) | IDT and other vendors accept IUPAC ambiguity codes; degenerate pools are biologically appropriate; `--resolve-spacer` flag available for explicit bases |
| Initiator data as code constants | Avoids external file dependencies; sequences are stable and published |
| Lowercase/uppercase output convention | Clear visual separation of binding vs initiator regions in oligo orders |
| Odd/even naming for P1/P2 | P1 (odd) binds right half with init_b, P2 (even) binds left half with init_a — initiators point inward toward gap; matches R script convention |
| Multipage Streamlit with `st.navigation` | Clean separation between smFISH and HCR; independent sidebars; no UI clutter from conditional visibility |
| Test-gated development phases | Each phase's tests written first and must pass before proceeding; prevents accumulation of hidden bugs |
| R/L masks always symmetric, P/B masks potentially asymmetric | Repeat regions and synthesis quality affect both halves regardless; off-target specificity is the property being relaxed for the lenient half |
