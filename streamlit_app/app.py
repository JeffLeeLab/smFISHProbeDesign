"""Streamlit GUI for ProbeDesign — multipage entry point.

Builds the top navigation bar and routes between the pages in pages/
(smFISH, HCR, Isoform Explorer, TransQuant, IDT Plate) using st.navigation.
Each page owns its own sidebar.
"""

import sys
from pathlib import Path

import streamlit as st

# Ensure the probedesign package is importable from the repo root
_repo_root = Path(__file__).resolve().parent.parent
if str(_repo_root / "src") not in sys.path:
    sys.path.insert(0, str(_repo_root / "src"))

# Add streamlit_app directory to path for utils import
_app_dir = Path(__file__).resolve().parent
if str(_app_dir) not in sys.path:
    sys.path.insert(0, str(_app_dir))


# ── Page config ──────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="ProbeDesign",
    layout="wide",
)

# ── Navigation ───────────────────────────────────────────────────────────────

smfish_page    = st.Page("pages/smfish.py",          title="smFISH Probe Design")
hcr_page       = st.Page("pages/hcr.py",             title="HCR Probe Design")
isoform_page   = st.Page("pages/isoform_explorer.py", title="Isoform Explorer")
transquant_page = st.Page("pages/transquant.py",      title="TransQuant")
idt_page       = st.Page("pages/idt_plate.py",        title="IDT Plate Builder")

pg = st.navigation(
    [smfish_page, hcr_page, isoform_page, transquant_page, idt_page], position="hidden"
)

# Top navigation bar: one row of buttons, each as wide as its label
with st.container(horizontal=True):
    for label, page in [
        ("smFISH", smfish_page),
        ("HCR", hcr_page),
        ("Isoform Explorer", isoform_page),
        ("TransQuant", transquant_page),
        ("IDT Plate", idt_page),
    ]:
        if st.button(label, type="primary" if pg.title == page.title else "secondary"):
            st.switch_page(page)
st.divider()

pg.run()
