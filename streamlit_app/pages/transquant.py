"""TransQuant — Streamlit page.

Embeds the TransQuant probe-weight tool (external package `transquant_w`,
https://github.com/JeffLeeLab/TransQuant-probe-weight) as a page of this
app. Computes the probe weight factor W, gene length L and probe
localisation profile N for an smFISH probe set.

Self-contained: no data handoff from the smFISH / HCR / Isoform Explorer
tabs. `transquant_w.ui.render()` takes no parameters, owns none of the page
config, and reads no session state shared with the rest of this app.
"""
from __future__ import annotations

import streamlit as st

try:
    from transquant_w.ui import render
except ImportError:
    st.error(
        "The `transquant_w` package is not installed. Install it into the "
        "probedesign environment with:\n\n"
        "```\n"
        "pip install \"git+https://github.com/JeffLeeLab/TransQuant-probe-weight.git@v0.2.0\"\n"
        "```\n\n"
        "Or use the standalone hosted app: "
        "https://jeffleelab.github.io/TransQuant-probe-weight/"
    )
else:
    render()
