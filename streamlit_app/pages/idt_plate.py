"""IDT Plate Builder — combine designed probe sets into one IDT plate order.

Takes _oligos.txt / _HCR_oligos.txt files (or batch-result zips) from the
smFISH and HCR tabs and builds the Excel file IDT's 96-well plate upload
accepts. Layout rules live in probedesign.idt.build_plates().
"""

import sys
import zipfile
from pathlib import Path

import pandas as pd
import streamlit as st

# Ensure the probedesign package is importable from the repo root
_repo_root = Path(__file__).resolve().parent.parent.parent
if str(_repo_root / "src") not in sys.path:
    sys.path.insert(0, str(_repo_root / "src"))

# Add streamlit_app to path for utils import
_app_dir = Path(__file__).resolve().parent.parent
if str(_app_dir) not in sys.path:
    sys.path.insert(0, str(_app_dir))

from utils import XLSX_MIME
from probedesign import __version__
from probedesign.idt import (
    DEFAULT_MIN_OLIGOS,
    N_COLUMNS,
    ROWS,
    ProbeSet,
    build_plates,
    plate_filename,
    read_oligo_sets,
    write_plates_xlsx,
)

# Translucent fills so set colours read on light and dark themes
_SET_COLOURS = [
    "rgba(31,119,180,0.30)", "rgba(255,127,14,0.30)", "rgba(44,160,44,0.30)",
    "rgba(214,39,40,0.30)", "rgba(148,103,189,0.30)", "rgba(140,86,75,0.30)",
    "rgba(227,119,194,0.30)", "rgba(188,189,34,0.30)",
]


# ── Sidebar ──────────────────────────────────────────────────────────────────

with st.sidebar:
    st.header("IDT Plate")
    st.caption(f"ProbeDesign v{__version__}")
    min_oligos = st.number_input(
        "Minimum oligos per plate",
        min_value=0,
        max_value=96,
        value=DEFAULT_MIN_OLIGOS,
        step=1,
        help="IDT's minimum for a 96-well plate order: 24 for standard plates, "
             "48 for HOTPlates. Sets are packed so every plate meets it where possible.",
    )
    st.markdown(
        "**Layout rules**\n"
        "- Wells fill down each column: A1–H1, then A2 …\n"
        "- Each probe set starts on a new column.\n"
        "- A set is never split across plates; one that doesn't fit moves to the next plate.\n"
        "- One worksheet per plate, named after its genes."
    )


def _plate_grid(plate):
    """8 x 12 well grid of gene labels plus a matching grid of set colours."""
    labels = pd.DataFrame("", index=list(ROWS), columns=range(1, N_COLUMNS + 1))
    colours = labels.copy()
    wells = iter(plate.wells)
    for k, pset in enumerate(plate.sets):
        for _ in pset.oligos:
            well = next(wells)[0]
            row, col = well[0], int(well[1:])
            labels.loc[row, col] = pset.gene
            colours.loc[row, col] = f"background-color: {_SET_COLOURS[k % len(_SET_COLOURS)]}"
    return labels.style.apply(lambda _: colours, axis=None)


def _set_summary(plate) -> str:
    parts, start = [], 0
    for pset in plate.sets:
        first, last = plate.wells[start][0], plate.wells[start + len(pset.oligos) - 1][0]
        parts.append(f"{pset.name} ({len(pset.oligos)}, {first}–{last})")
        start += len(pset.oligos)
    return " · ".join(parts)


# ── Main ─────────────────────────────────────────────────────────────────────

st.title("IDT Plate Builder")
st.caption("Combine probe sets from smFISH and HCR runs into one IDT 96-well plate order")

uploads = st.file_uploader(
    "Oligo files",
    type=["txt", "zip"],
    accept_multiple_files=True,
    help="_oligos.txt or _HCR_oligos.txt files from the smFISH / HCR tabs, "
         "or a batch results .zip.",
)

if not uploads:
    st.info("Upload the oligo files of the probe sets to order together.")
    st.stop()

sets = []
for uf in uploads:
    try:
        sets.extend(read_oligo_sets(uf.name, uf.getvalue()))
    except (ValueError, UnicodeDecodeError, zipfile.BadZipFile) as e:
        st.error(f"{uf.name}: {e}")
if not sets:
    st.stop()

st.subheader("Probe sets")
st.caption("Sets go onto plates in Order. Gene names the worksheets and the file; "
           "edit it if the region wasn't stripped correctly.")
table = pd.DataFrame({
    "Order": range(1, len(sets) + 1),
    "Include": True,
    "Probe set": [s.name for s in sets],
    "Gene": [s.gene for s in sets],
    "Oligos": [len(s.oligos) for s in sets],
})
upload_sig = abs(hash(tuple((uf.name, uf.size) for uf in uploads)))
edited = st.data_editor(
    table,
    key=f"idt_sets_{upload_sig}",
    hide_index=True,
    width="stretch",
    disabled=["Probe set", "Oligos"],
    column_config={
        "Order": st.column_config.NumberColumn(min_value=1, step=1),
        "Include": st.column_config.CheckboxColumn(),
        "Gene": st.column_config.TextColumn(
            help="Probe-set name without its region (cmm, CDS, exon …)"),
    },
)

edited = edited.assign(_i=range(len(sets)), Order=edited["Order"].fillna(len(sets) + 1))
chosen = [
    ProbeSet(sets[i].name, sets[i].oligos, (gene or "").strip() or sets[i].gene)
    for i, include, gene in edited.sort_values(["Order", "_i"])[["_i", "Include", "Gene"]].itertuples(index=False)
    if include
]
if not chosen:
    st.info("Include at least one probe set.")
    st.stop()

plates, warnings = build_plates(chosen, min_oligos=int(min_oligos))
for w in warnings:
    st.warning(w)

filename = plate_filename(plates)
try:
    xlsx = write_plates_xlsx(plates)
except ImportError as e:
    st.error(str(e))
else:
    st.download_button(
        label=f"Download {filename}",
        data=xlsx,
        file_name=filename,
        mime=XLSX_MIME,
        type="primary",
        key="idt_dl_plate",
    )

st.subheader(f"{len(plates)} plate{'s' if len(plates) != 1 else ''}")
for plate in plates:
    st.markdown(f"**{plate.sheet_name}** — {len(plate.wells)} oligos")
    st.caption(_set_summary(plate))
    st.dataframe(_plate_grid(plate), width="stretch")
