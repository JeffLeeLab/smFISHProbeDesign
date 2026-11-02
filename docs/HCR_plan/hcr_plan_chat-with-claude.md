# HCR Probe Designer — Plan Review Chat Summary

> **Status**: Historical review notes for [hcr_designer_plan.md](../hcr_designer_plan.md), implemented in `8917ef4`.

This document summarises the key decisions, corrections, and clarifications made during
the review of `docs/hcr_designer_plan.md` with Claude.

---

## 1. Asymmetric ΔG and bowtie modes are independent

**Question:** When both `--asymmetric-gibbs` and `--asymmetric-bowtie` are enabled,
should the strict/lenient roles be coupled (same half strict for both) or independent?

**Decision: Independent.** The two filters evaluate different biological concerns:

- ΔG asymmetry: at least one half must bind well (strict ΔG range)
- Bowtie asymmetry: at least one half must be off-target clean

These do not need to be the same half. A single half cannot nucleate HCR amplification
regardless of how stably it binds off-target, so the two safety guarantees are orthogonal.

**Algorithm (from §4.2):**

```python
# ΔG filter (independent)
if asymmetric_gibbs:
    gibbs_A = badness_strict[p] + badness_lenient[p+27]   # left ΔG-strict
    gibbs_B = badness_lenient[p] + badness_strict[p+27]   # right ΔG-strict
    gibbs_score = min(gibbs_A, gibbs_B)
else:
    gibbs_score = badness_strict[p] + badness_strict[p+27]

# Bowtie filter (independent)
if asymmetric_bowtie:
    bowtie_ok = (not left_offmasked) or (not right_offmasked)  # at least one clean
else:
    bowtie_ok = (not left_offmasked) and (not right_offmasked) # both clean
```

The original plan incorrectly coupled the roles ("orientation A = left strict for both ΔG
and bowtie"). This was corrected.

---

## 2. Bowtie masking mechanism

**How it works:**

1. `sequence_to_nmers(seq, mer_length)` generates all sliding-window substrings of a
   fixed length (12, 14, or 16 nt) from the full target sequence — including substrings
   that span the 2-nt gap between probe halves.
2. All substrings are piped to bowtie as FASTA reads. Bowtie reports per-position hit counts.
3. `hits_to_mask()` marks a **block of `mer_length` nucleotides** starting at each
   hit position (conservative — every nucleotide in the matching k-mer is flagged).
4. Genome masking uses three tiers (12-mer/threshold 4000, 14-mer/threshold 500,
   16-mer/threshold 20), OR-combined. Pseudogene masking uses 16-mers/threshold 0.
5. The result is a **nucleotide-level binary mask** (one value per position).

**Key point:** Bowtie is entirely gap-unaware. The gap positions are just regular
nucleotides at this stage. Gap-awareness only applies downstream when evaluating pairs.

---

## 3. Gap-skipping: when it actually matters

**Question:** Can a bowtie mask cover only the 2-nt gap without touching either half?

**Answer: No, for bowtie masks.** The minimum n-mer is 12 nt. Any 12-mer that marks a
gap position also marks at least 10 positions in one or both halves. There is no
arrangement of a 12+ nt block that marks only positions 25-26 without spilling into
the left half (positions 0-24) or right half (positions 27-51).

**Gap-skipping matters for R-masks (repeat/N masks), not bowtie masks.** If the input
FASTA has N's at positions 25-26 (the gap region), the R-mask marks exactly those 2
positions without touching either half. Without gap-skipping, this would incorrectly
reject a valid pair.

**Tests corrected accordingly:**

| Test | What it covers |
|-|-|
| `test_bowtie_mask_spanning_gap` | Realistic 12-nt block (`[0]*20 + [1]*12 + [0]*20`) bleeding from left half through gap into right half — pair correctly rejected (both halves contaminated) |
| `test_r_mask_n_in_gap_ignored` | N's at gap positions produce R-mask `[0]*25 + [1]*2 + [0]*25` — pair correctly passes (gap positions skipped) |

The previously planned test `[0]*25 + [1]*2 + [0]*25` as a bowtie mask was unrealistic
and has been replaced.

---

## 4. Pair scoring and DP optimality

**Pair badness** = sum of squared ΔG deviations from target across both halves:
`(ΔG_left - target)² + (ΔG_right - target)²`

Halves outside their allowable range get `inf` badness (hard reject).

**DP priority:**
1. Maximise pair count (primary — more probes = stronger signal)
2. Minimise total badness among solutions with the same count (secondary — prefer
   uniform binding)

This lexicographic ordering is appropriate: in FISH, signal scales with probe count, so
trading one probe for marginally better ΔG on the remaining probes is almost never
worthwhile.

---

## 5. Function name corrections

The plan referenced non-existent function names. Corrected to match `thermodynamics.py`:

| Old (wrong) | Correct |
|-|-|
| `calculate_gibbs_fe` | `gibbs_rna_dna` |
| `calculate_tm` | `tm_rna_dna` |

---

## 6. Fixture sequence lengths (§9.0)

All 8 simulated FASTA test fixtures had incorrect lengths (53-55 nt instead of the
required 52 nt, or 52 nt instead of 51 for TOO_SHORT). All were corrected and verified
by running `len()` checks in Python.

Key fix: the BOTH_STRICT fixture gap is `at` (positions 25-26), both halves are the same
25-mer `atcgatcgatcgatcgatcgatcga`.

---

## 7. Streamlit navigation

`st.navigation(position="top")` does not exist in Streamlit 1.54. The workaround:

```python
pg = st.navigation([smfish_page, hcr_page], position="hidden")
# Custom top bar with st.columns + st.button + st.switch_page()
```

---

## 8. WW spacer (amplifiers B7-B17)

Default behaviour: output `WW` directly as IUPAC ambiguity code. IDT and other vendors
accept degenerate pools. `--resolve-spacer AA` / `TT` etc. available for explicit bases.

---

## 9. Key plan sections changed from original

| Section | Change |
|-|-|
| §1.5 Pitfalls | Added junction-in-gap bug, N-in-gap handling, WW IUPAC default |
| §2.3 Unchanged files | Added masking.py; corrected thermodynamic function names |
| §4.1 Badness | Two-array approach (strict + lenient), explicit gap junction check |
| §4.2 Masking | Full rewrite: single bowtie run, nucleotide→pair conversion, independent ΔG/bowtie evaluation, mask symmetry table |
| §4.4 WW spacer | Default outputs WW; `--resolve-spacer` flag for explicit bases |
| §5.2 CLI | Renamed `--lenient-gibbs-max` → `--lenient-gibbs-min`, `--initiator-spacer` → `--resolve-spacer` |
| §7 Streamlit | Multipage with `position="hidden"` + custom top nav workaround |
| §8 Phases | Test-gated development; Phase 0 (scaffold tests) before any implementation |
| §9 Tests | 10 simulated fixtures, 29 unit tests, realistic masking scenarios |
| §10 Performance | Single bowtie run (not two) — identical runtime to smFISH |
| §13 Summary | Added independent-filters decision, gap-skipping rationale, WW IUPAC |
