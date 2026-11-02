# Legacy code (archived)

Code inherited from the upstream Raj Lab MATLAB repository. It is kept for reference only: nothing here is used by the Python package (`src/probedesign/`), the Streamlit app or `setup_all.sh`, and none of it is maintained.

| Folder | What it was |
|-|-|
| `probedesign/` | Original MATLAB probe designer. `findprobesLocal.m` is the algorithm the Python package was ported from |
| `panprobedesign/` | MATLAB pan-probe design across aligned sequences (`findpanprobesHD.m`) |
| `DesignServer/` | 2011 web design server: Python 2 scripts plus MATLAB callers for local bowtie |
| `maskprobes/` | Early stand-alone Python probe design/masking scripts |

The pseudogene FASTAs that used to sit in `probedesign/pseudogeneDBs/` are live data and now live in [`data/pseudogeneDBs/`](../data/pseudogeneDBs/).
