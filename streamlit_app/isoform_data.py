"""Data-access helpers and isoform structure plot for the Isoform Explorer page."""
from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import streamlit as st

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "isoform_explorer"

SPECIES_LABELS: dict[str, str] = {
    "Hsap": "Human (Hsap)",
    "Mmus": "Mouse (Mmus)",
    "Dmel": "Drosophila (Dmel)",
}

# TSL colour palette: 1 (best/darkest) → 5 (worst/lightest), sequential blue
_TSL_COLORS: dict[int, str] = {
    1: "#1f4e79",
    2: "#2e75b6",
    3: "#5ba3d9",
    4: "#9dc3e6",
    5: "#c9dff0",
}
_NA_COLOR        = "#888888"   # exon with no TSL annotation
_EXCL_ALPHA      = 0.25        # isoforms excluded by the TSL filter
_CONSTITUTIVE_COLOR = "#bb5a38"  # app accent / primaryColor
_BASELINE_COLOR  = "#cccccc"
_BG_COLOR        = "#f4f3ed"
_TEXT_COLOR      = "#3d3a2a"
_BORDER_COLOR    = "#d3d2ca"


# ── Parquet path helpers ──────────────────────────────────────────────────────

def _summary_path(species: str) -> Path:
    return DATA_DIR / species / f"{species}_summary.parquet"


def _intervals_path(species: str) -> Path:
    return DATA_DIR / species / f"{species}_intervals.parquet"


def safe_tsl_label(tsl_cutoff: str) -> str:
    """Return a FASTA-header-safe version of a tsl_cutoff string.

    Replaces operators and special characters with plain ASCII equivalents:
      'TSL<=2'              -> 'TSL_leq_2'
      'TSL==3'              -> 'TSL_eq_3'
      'all (no TSL)'        -> 'all_no_TSL'
      'all isoforms (no TSL)' -> 'all_isoforms_no_TSL'
    """
    s = tsl_cutoff
    s = s.replace("<=", "_leq_")
    s = s.replace("==", "_eq_")
    s = s.replace(" ", "_")
    s = s.replace("(", "").replace(")", "")
    return s.strip("_")


# ── Data-access helpers (cached) ──────────────────────────────────────────────

@st.cache_data(show_spinner="Loading gene list…")
def load_species_metadata(species: str) -> pd.DataFrame:
    """Light columnar read — only the columns needed for the search dropdown."""
    return pd.read_parquet(
        _summary_path(species),
        columns=["gene_id", "gene_name", "gene_biotype", "constitutive_length", "constitutive_exist"],
    )


def get_gene_record(species: str, gene_id: str) -> pd.Series:
    """Full row for one gene from the summary parquet (includes sequence)."""
    df = pd.read_parquet(
        _summary_path(species),
        filters=[("gene_id", "==", gene_id)],
    )
    if df.empty:
        raise KeyError(f"Gene {gene_id!r} not found in {species} summary parquet")
    return df.iloc[0]


@st.cache_data(show_spinner="Loading isoform intervals…")
def get_gene_intervals(species: str, gene_id: str) -> pd.DataFrame:
    """All interval rows for one gene (exons + constitutive track)."""
    return pd.read_parquet(
        _intervals_path(species),
        filters=[("gene_id", "==", gene_id)],
    )


# ── Isoform structure plot ────────────────────────────────────────────────────

def plot_isoforms(intervals_df: pd.DataFrame, record: pd.Series) -> plt.Figure:
    """
    Port of plot_isoforms() from isoform_explorer_script.qmd, extended with
    TSL-based exon colouring.

    - Mammals: each isoform coloured by transcript_support_level (1=darkest).
      Isoforms excluded by the adaptive TSL filter are drawn at reduced alpha.
    - Dmel / no-TSL genes: all exons uniform gray, no TSL legend.
    - CONSTITUTIVE track: accent colour #bb5a38 at top of the plot.
    """
    gene_id    = str(record["gene_id"])
    gene_name  = str(record["gene_name"])
    strand     = str(record["strand"])
    const_len  = int(record["constitutive_length"])

    exon_df  = intervals_df[intervals_df["feature"] == "exon"].copy()
    const_df = intervals_df[intervals_df["feature"] == "constitutive"].copy()

    # Transcript order: alphabetical; CONSTITUTIVE lane at top
    tx_ids = sorted(exon_df["transcript_id"].unique())
    tracks = (["CONSTITUTIVE"] if not const_df.empty else []) + tx_ids
    n_tracks = len(tracks)

    # Per-transcript TSL level and inclusion flag
    has_tsl = exon_df["transcript_support_level"].notna().any()
    tx_info: dict[str, dict] = {}
    for tx in tx_ids:
        rows = exon_df[exon_df["transcript_id"] == tx]
        tsl_series = rows["transcript_support_level"].dropna()
        tsl = int(tsl_series.iloc[0]) if len(tsl_series) > 0 else None
        included = bool(rows["tsl_included"].iloc[0]) if len(rows) > 0 else True
        tx_info[tx] = {"tsl": tsl, "included": included}

    # Genome span
    all_coords = pd.concat([
        intervals_df["start"],
        intervals_df["end"],
    ])
    g_min = int(all_coords.min())
    g_max = int(all_coords.max())

    fig_h = max(4.0, n_tracks * 0.45 + 1.8)
    fig, ax = plt.subplots(figsize=(12, fig_h))
    fig.patch.set_facecolor(_BG_COLOR)
    ax.set_facecolor(_BG_COLOR)

    # y-positions: track index 0 = bottom; reversed so CONSTITUTIVE is at top
    y_pos = {track: (n_tracks - 1 - i) for i, track in enumerate(tracks)}

    for track in tracks:
        y = y_pos[track]

        if track == "CONSTITUTIVE":
            ax.plot([g_min, g_max], [y, y], color=_BASELINE_COLOR, lw=0.8, zorder=1)
            for _, row in const_df.iterrows():
                ax.plot(
                    [int(row["start"]), int(row["end"])], [y, y],
                    color=_CONSTITUTIVE_COLOR, lw=14, solid_capstyle="butt", zorder=2,
                )
        else:
            info = tx_info[track]
            tsl = info["tsl"]
            included = info["included"]
            color = _TSL_COLORS.get(tsl, _NA_COLOR) if has_tsl else "#606060"
            alpha = 1.0 if included else _EXCL_ALPHA

            ax.plot([g_min, g_max], [y, y], color=_BASELINE_COLOR, lw=0.8, zorder=1, alpha=alpha)
            tx_rows = exon_df[exon_df["transcript_id"] == track]
            for _, row in tx_rows.iterrows():
                ax.plot(
                    [int(row["start"]), int(row["end"])], [y, y],
                    color=color, lw=6, solid_capstyle="butt", zorder=2, alpha=alpha,
                )

    # Axes
    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels(list(tracks), fontsize=9, color=_TEXT_COLOR)
    ax.set_ylim(-0.6, n_tracks - 0.4)
    ax.set_xlabel("Chromosomal coordinate", fontsize=11, color=_TEXT_COLOR)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.tick_params(colors=_TEXT_COLOR, labelsize=10)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(_BORDER_COLOR)
    ax.spines["left"].set_color(_BORDER_COLOR)

    # Title
    ax.set_title(
        f"{gene_id}  {gene_name}  ({strand} strand, constitutive {const_len:,} nt)",
        fontsize=12, color=_TEXT_COLOR, pad=8,
    )

    # Legend — placed outside the axes on the right to avoid overlapping exons
    legend_handles = [mpatches.Patch(color=_CONSTITUTIVE_COLOR, label="Constitutive")]
    if has_tsl:
        for level, color in _TSL_COLORS.items():
            legend_handles.append(mpatches.Patch(color=color, label=f"TSL {level}"))
        legend_handles.append(mpatches.Patch(color=_NA_COLOR, label="TSL NA"))
        legend_handles.append(
            mpatches.Patch(facecolor="#aaaaaa", alpha=_EXCL_ALPHA, label="Excluded by filter")
        )
    else:
        legend_handles.append(mpatches.Patch(color="#606060", label="Exon"))

    ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        borderaxespad=0,
        fontsize=8,
        framealpha=0.85,
        facecolor=_BG_COLOR,
        edgecolor=_BORDER_COLOR,
    )

    fig.tight_layout()
    fig.subplots_adjust(right=0.82)
    return fig
