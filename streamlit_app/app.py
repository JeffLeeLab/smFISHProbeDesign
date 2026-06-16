"""Streamlit GUI for ProbeDesign — multipage entry point.

Routes between smFISH and HCR probe design pages using st.navigation.
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

pg = st.navigation([smfish_page, hcr_page, isoform_page], position="hidden")

# Top navigation bar
nav_cols = st.columns([1, 1, 2, 4])
with nav_cols[0]:
    if st.button("smFISH", use_container_width=True,
                 type="primary" if pg.title == "smFISH Probe Design" else "secondary"):
        st.switch_page(smfish_page)
with nav_cols[1]:
    if st.button("HCR", use_container_width=True,
                 type="primary" if pg.title == "HCR Probe Design" else "secondary"):
        st.switch_page(hcr_page)
with nav_cols[2]:
    if st.button("Isoform-explorer", use_container_width=True,
                 type="primary" if pg.title == "Isoform Explorer" else "secondary"):
        st.switch_page(isoform_page)
st.divider()

pg.run()
