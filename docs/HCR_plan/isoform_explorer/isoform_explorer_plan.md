# Isoform Explorer — Implementation Plan

## Objective

- Add a new page to the existing Streamlit app (same design language) for exploring the **constitutive sequence** of an RNA target — the portion present in *all* isoforms of a gene.
- Rationale: genes with complex isoform structures should be probed against a constitutive region so the probe set hybridises every transcript. This page lets the user inspect isoform architecture, see whether a usable constitutive sequence exists, and export it.
- The exported N-collapsed constitutive sequence is the intended **input to the smFISH/HCR probe designer** (the `N` markers stop probes spanning exon–exon junctions).
- This page is **self-contained**: it does not touch the existing smFISH/HCR designer code paths and is not exposed via the CLI.

---

## Decisions (resolved 2026-06-16)

| Topic | Decision | Notes |
|-|-|-|
| Storage format | **Parquet (zstd)** | Chosen over SQLite after measurement — see below. Columnar, Python-native, compresses better than `.qs`. |
| Distribution | **Commit directly** | Files are well under GitHub's 100 MB/file limit. `.qs` and the source CSVs are removed from the repo once converted. |
| Plot rendering | **On-the-fly matplotlib** | The plot is a handful of horizontal segments per transcript — milliseconds per gene. No stored images. |
| Species at launch | **All three (Hsap, Mmus, Dmel)** | Mmus data is now present in `data/isoform_explorer/Mmus/`. |
| Mammalian filtering | **TSL only** | Use the `*_tsl_*` datasets for Hsap/Mmus (stringent transcript-support filtering). Dmel has no TSL variant. |
| Environment | **`probedesign` mamba env** | All Python/data work runs inside `mamba activate probedesign` (`/Users/jefflee/miniforge3/envs/probedesign`). Never the base CLI env. |

### TSL cutoff used to generate the constitutive sequence (mammals)

From `isoform_explorer_script.qmd` (the "Add TSL support" chunk, lines ~1132–1155), the constitutive N-collapsed sequence for **Hsap/Mmus** was computed from an *adaptive, per-gene* subset of isoforms:

1. If the gene's best (minimum) TSL ≤ 2 → keep all transcripts with **TSL ≤ 2** (TSL1 + TSL2).
2. Else if the best TSL is 3–5 → keep only transcripts at that **single best level** (`transcript_support_level == min`).
3. Else (no usable TSL annotation) → keep **all** transcripts.

So the cutoff is *not* a single fixed level. **Dmel** has no TSL and uses **all isoforms**. This provenance must travel with every export (FASTA header + a per-gene column — see below), since the effective cutoff differs gene to gene.

### Why Parquet, not SQLite (measured on Hsap, 55,184 genes)

| Format | File size |
|-|-|
| Source CSV (seq inline) | 54.2 MB |
| **Parquet, zstd** | **16.4 MB** |
| SQLite (indexed) | 66.1 MB |
| Metadata-only parquet (no seq) | 0.54 MB |

SQLite stores sequences uncompressed and ends up *larger* than the CSV, conflicting with the "compress like `.qs`, commit directly" goal. Parquet+zstd is 16 MB (smaller than the original 62 MB `.qs`) and, being columnar, lets the search dropdown read only the light `gene_id`/`gene_name`/`gene_biotype`/`constitutive_length` columns while skipping the 16 MB sequence column entirely.

---

## Data architecture

### Current source files (per species, under `data/isoform_explorer/{Hsap,Mmus,Dmel}/`)

| File | Contents | Python-readable? |
|-|-|-|
| `*_constitutive-sequence_summary.csv` | gene metadata + full N-collapsed sequence | yes (but heavy) |
| `*_exon_list.qs` | per-gene `exon_intervals` (transcript-level exon coords) | **no — R-native** |
| `*_constitutive-region_list.qs` | per-gene `constitutive_intervals` | **no — R-native** |
| `*_constitutive-sequence_list.qs` | above + extracted sequences | **no — R-native** |

`.qs` cannot be read from Python, so a **one-time R re-export** is unavoidable regardless of target format. Mammalian species use the `*_tsl_*` variants; Dmel uses the plain variants (`_ensembl_v113_proteincoding_noncoding_*`).

### Target files (per species, committed to repo)

Two Parquet files per species replace the CSV **and** all three `.qs` files:

**1. `{Sp}_summary.parquet`** — replaces the CSV. One row per gene.

| column | type | source |
|-|-|-|
| gene_id | str | summary |
| gene_name | str | summary |
| gene_biotype | str | summary |
| strand | str (`+`/`-`) | summary |
| isoform_count | int | summary (TSL-filtered count for mammals) |
| constitutive_exist | bool | summary |
| constitutive_length | int | summary |
| constitutive_seq_Ncollapsed | str | summary (the export sequence) |
| tsl_cutoff | str | **provenance of the adaptive TSL rule applied to this gene**, e.g. `TSL<=2`, `TSL==4`, or `all (no TSL)` for Dmel. Derived from `min_transcript_support_level` per the rule above. Used in the FASTA header. |

**2. `{Sp}_intervals.parquet`** — replaces the `.qs` lists; long-format, drives the plot.

| column | type | meaning |
|-|-|-|
| gene_id | str | join key |
| transcript_id | str | which isoform (`CONSTITUTIVE` for the constitutive track) |
| feature | str | `exon` or `constitutive` |
| start | int | genomic start (0-based, as in the qs/bedtools output) |
| end | int | genomic end |
| strand | str | `+`/`-` |
| transcript_support_level | int/null | per-transcript TSL (1–5, null = NA/`CONSTITUTIVE` track). **Required for TSL-coloured plotting.** Null for all Dmel transcripts. |
| tsl_included | bool | whether this transcript was in the TSL-filtered set used for the constitutive calc (lets the plot distinguish included vs excluded isoforms). Always true for Dmel. |

> The plot shows **all** isoforms (full `exon_intervals`), coloured by TSL — so the user can see which isoforms were included vs excluded by the cutoff. The TSL columns come from `exon_intervals$transcript_support_level` in the qs (present for mammals).

Sort both files by `gene_id` so pyarrow predicate-pushdown reads only the relevant row groups on per-gene lookup.

### Estimated repo footprint

~16 MB summary + ~10–15 MB intervals per species × 3 ≈ **80–90 MB total**, comparable to the single current Hsap `.qs`. Within the "acceptable bloat, commit directly" call. Hsap is the largest; Mmus/Dmel are smaller.

### R export script (to be written)

A standalone R script (e.g. `docs/HCR_plan/isoform_explorer/export_to_parquet.R`) that, per species:
0. **Verify TSL column is present (mammals).** Before flattening, confirm `exon_intervals` in the mammalian qs lists actually carries a per-row `transcript_support_level` column — the plot colouring depends on it. The qmd reads it via `min(exon_intervals$transcript_support_level)`, so it should be there, but check explicitly: e.g. `pcnc_constitutive_seq_list[[1]]$exon_intervals %>% names()` and `... %>% pull(transcript_support_level) %>% unique()` for Hsap and Mmus. If the column is missing or all-NA, the TSL-coloured plot cannot be built and the export step must first re-join TSL from the GTF/exon list. Dmel has no TSL and is expected to lack the column (uniform-colour fallback).
1. `qread()` the appropriate `*_constitutive-sequence_list.qs` (TSL variant for mammals).
2. Build the summary tibble (already exists as the CSV logic in `isoform_explorer_script.qmd`), derive the per-gene `tsl_cutoff` string from `min_transcript_support_level` using the adaptive rule above, and `arrow::write_parquet(compression = "zstd")`.
3. Flatten `exon_intervals` (carrying `transcript_support_level` and a `tsl_included` flag derived from `tsl_filtered_transcript_ids`) + `constitutive_intervals` across all genes into the long intervals tibble and write the second parquet.
4. After verifying the Python app reads them, delete the `.qs` and `.csv` source files from the repo.

Requires the R `arrow` package. This script is the only remaining R dependency; it is run once per data refresh, not by end users.

---

## Page layout

New page `streamlit_app/pages/isoform_explorer.py`, registered in `app.py` alongside the smFISH/HCR pages.

### Navigation (`app.py`)

Add a third `st.Page` and a third top-nav button. Current layout is `st.columns([1, 1, 6])` with two buttons → change to three buttons (e.g. `st.columns([1, 1, 1, 5])`), pattern-matching the existing primary/secondary highlight logic.

### Sidebar (left)

- **Species selector** — `st.selectbox(["Hsap", "Mmus", "Dmel"])` (display friendly labels, e.g. "Human (Hsap)").
- **Gene search** — searchable dropdown over genes for the selected species:
  - Populate from the light columns of `{Sp}_summary.parquet` only (read `gene_id`, `gene_name`, `gene_biotype`, `constitutive_length` — skips the sequence column).
  - Label format: `gene_name — gene_id` so both are searchable by type-ahead.
  - **Performance note to validate:** `st.selectbox` with ~55k options ships the full option list to the browser (~1 MB). If sluggish, fall back to a `st.text_input` filter that narrows the list before the selectbox. Cache the per-species metadata frame with `st.cache_data` keyed on species.

### Main panel (right)

On gene selection, show a detail report:

1. **Isoform structure plot** (matplotlib, on the fly — see below).
2. **Status messages**, driven by `constitutive_exist` / `constitutive_length`:
   - No constitutive sequence → warning: *"This gene has no constitutive sequence shared by all isoforms. Inspect the isoform structure above and consider using a MANE/reference isoform instead."*
   - Constitutive sequence < 500 nt → warning: *"This gene has a short constitutive sequence (<500 nt), which may not fit enough probes. Inspect the isoform structure and consider a MANE/reference isoform."*
   - Otherwise → success/info showing the length.
3. **Sequence text box** — `st.code` (monospace, copy button) of the N-collapsed constitutive sequence, so it can be copy-pasted directly. Read on demand for the selected gene (filtered parquet read).
4. **Download button** — FASTA of the N-collapsed constitutive sequence.

### FASTA download format

> Correction to brain-dump: the original note `>{gene_id}_{gene_name}_{constitutive-seq-N-collapsed}` would put the whole sequence in the header. Intended format:

```
>{gene_id}_{gene_name} constitutive_seq_Ncollapsed | Ensembl v113 | {tsl_cutoff}
ACGT...N...ACGT
```

`{tsl_cutoff}` is the per-gene provenance string (e.g. `TSL<=2`, `TSL==4`, or `all isoforms (no TSL)` for Dmel) so every export records exactly which isoform subset produced the constitutive sequence. Header on line 1, sequence on the following line(s). Filename e.g. `{gene_name}_{gene_id}_constitutive.fa`.

---

## Plotting (matplotlib, on the fly)

Port `plot_isoforms()` from `isoform_explorer_script.qmd` (the ggplot version, lines ~404–436) to matplotlib **and extend it with TSL colouring** (the qmd version uses a uniform `gray30` for all exons — TSL colouring is a new addition requested here):

- One horizontal lane per `transcript_id`; thin baseline segment spanning min–max coords (gray), thick segments for exons.
- **Exon colour encodes TSL level (mammals):** colour each isoform's exons by its `transcript_support_level` using a fixed discrete palette, e.g. TSL1 → darkest, TSL2, TSL3, TSL4, TSL5 → lightest, NA → gray. Include a legend mapping colour → TSL. Isoforms *excluded* by the cutoff (`tsl_included == false`) should be visually de-emphasised (e.g. reduced alpha or a hatch) so it's clear which isoforms drove the constitutive region.
- **Dmel (no TSL):** all `transcript_support_level` are null → fall back to a single uniform exon colour (the original `gray30`), no TSL legend.
- A separate `CONSTITUTIVE` lane with thick segments in the app accent colour `#bb5a38` (matches `primaryColor`).
- Title: `{gene_id} {gene_name} ({strand} strand, constitutive {N} nt)`.
- X axis: chromosomal coordinate; Y axis: transcript IDs.
- Strand handling: mirror the qmd (minus-strand genes arranged by descending start) so layout matches the R reference.
- Render via `st.pyplot(fig)`; close the figure after to avoid leaks.
- Style to match the app theme (light background `#f4f3ed`, text `#3d3a2a`); keep it simple and legible for genes with many isoforms.

Data: filtered read of `{Sp}_intervals.parquet` where `gene_id == selected` (predicate pushdown), including the `transcript_support_level` / `tsl_included` columns. Cache the plot/data per `(species, gene_id)`.

---

## Implementation steps

1. **R export script** — write and run `export_to_parquet.R`; produce the two parquet files per species; verify schemas. (Blocking for everything else.)
2. **Dependencies** — add `matplotlib` and `pyarrow` to `environment.yml` (pip section; pandas already declared). `pyarrow` is needed for `pandas.read_parquet`.
3. **Data-access helpers** — small module (e.g. `streamlit_app/isoform_data.py`):
   - `load_species_metadata(species)` → cached light DataFrame for the dropdown.
   - `get_gene_record(species, gene_id)` → metadata + sequence (filtered read).
   - `get_gene_intervals(species, gene_id)` → intervals DataFrame for the plot.
4. **Plot function** — `plot_isoforms(intervals_df, record)` → matplotlib figure.
5. **Page** — `streamlit_app/pages/isoform_explorer.py` wiring sidebar + main panel.
6. **Navigation** — register the page and add the third nav button in `app.py`.
7. **Cleanup** — remove `.qs` and `*_summary.csv` from `data/isoform_explorer/` once the app reads parquet successfully.
8. **Docs** — update `CLAUDE.md` (repository structure + new page) and root `README.md`.

---

## Dependencies to add

- `matplotlib` — plotting (not currently a dependency).
- `pyarrow` — parquet I/O (present in the dev env but not declared in `environment.yml`).

No new heavy services; no DuckDB/SQLite needed.

---

## Testing

- **Run everything inside the project env:** `mamba activate probedesign` (or `/Users/jefflee/miniforge3/envs/probedesign/bin/python`). Never install into or run from the base CLI env.
- Spin up a minimal Streamlit server scoped to this page for quick iteration during development.
- Self-testing/automation: use the **playwright-cli** skill to drive the page (select species, search a known gene, assert plot + sequence + download render).
- Spot-check genes against the qmd reference plots (e.g. Dmel `Myc`/`pros`/`Syp`; a known no-constitutive gene; a short-constitutive gene) to confirm the matplotlib port matches the ggplot original.
- Validate dropdown responsiveness with the full Hsap set (~55k genes) and decide whether the text-filter fallback is needed.

---

## Open questions / future

- **Dropdown UX at 55k genes** — confirm `st.selectbox` is acceptable or implement the text-filter fallback (decision deferred to first real test).
- **Cross-link to designer** — future nicety: a "Send to probe designer" button that hands the constitutive FASTA to the smFISH/HCR page directly (out of scope for v1; the page stays self-contained for now).
- **Data versioning** — files are pinned to Ensembl v113; note the version in the FASTA header and on the page so exports are traceable.

---

## Design language

- Reuse the existing theme (`streamlit_app/.streamlit/config.toml`): warm palette, Styrene B font, accent `#bb5a38`.
- Match the smFISH/HCR pages: `st.header` in sidebar, `st.title` + `st.caption` at top of main panel, `st.dataframe`/`st.code`/`st.download_button`/`st.expander` for content blocks.
