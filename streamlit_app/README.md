# ProbeDesign — Web App

> For full installation, usage, and troubleshooting, see the [root README](../README.md).

Four tabs in the top navigation bar, one file each in `pages/`:

| Tab | File |
|-|-|
| smFISH probe design (single + batch) | `pages/smfish.py` |
| HCR split-initiator probe design (single + batch) | `pages/hcr.py` |
| Isoform Explorer | `pages/isoform_explorer.py` |
| TransQuant probe-weight tool | `pages/transquant.py` |

`app.py` only builds the navigation. Backend helpers are in `utils.py`.

## Launch

```bash
# Activate the probedesign environment
mamba activate probedesign   # or: conda activate probedesign

# From the repository root:
streamlit run streamlit_app/app.py
# Opens http://localhost:8501
```

## Stop

Press `Ctrl+C` in the terminal.

## Alternative port

```bash
streamlit run streamlit_app/app.py --server.port 8502
```
