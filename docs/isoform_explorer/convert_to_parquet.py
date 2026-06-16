#!/usr/bin/env python3
"""
convert_to_parquet.py

Reads the intermediate CSVs written by export_to_parquet.R and writes
zstd-compressed Parquet files to data/isoform_explorer/{Sp}/.

Run from the repo root with the probedesign env:
    /Users/jefflee/miniforge3/envs/probedesign/bin/python \
        docs/HCR_plan/isoform_explorer/convert_to_parquet.py
"""

import sys
import os
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

CSV_DIR  = Path("/tmp/isoform_explorer_csv")
DATA_DIR = Path("data/isoform_explorer")
SPECIES  = ["Hsap", "Mmus", "Dmel"]

SUMMARY_DTYPES = {
    "gene_id":                    "string",
    "gene_name":                  "string",
    "gene_biotype":               "string",
    "strand":                     "string",
    "isoform_count":              "Int32",
    "constitutive_exist":         "boolean",
    "constitutive_length":        "Int32",
    "constitutive_seq_Ncollapsed":"string",
    "tsl_cutoff":                 "string",
}

INTERVALS_DTYPES = {
    "gene_id":                   "string",
    "transcript_id":             "string",
    "feature":                   "string",
    "start":                     "Int32",
    "end":                       "Int32",
    "strand":                    "string",
    "transcript_support_level":  "Int8",
    "tsl_included":              "boolean",
}


def write_parquet(df: pd.DataFrame, out_path: Path, sort_by: str = "gene_id") -> None:
    df = df.sort_values(sort_by).reset_index(drop=True)
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(
        table,
        str(out_path),
        compression="zstd",
        compression_level=9,
        row_group_size=10_000,  # small groups help predicate-pushdown on gene_id
    )
    mb = out_path.stat().st_size / 1024 / 1024
    print(f"  -> {out_path}  ({mb:.1f} MB)")


def convert_species(sp: str) -> None:
    print(f"\n=== {sp} ===")
    out_dir = DATA_DIR / sp
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Summary ──────────────────────────────────────────────────────────────
    csv_path = CSV_DIR / f"{sp}_summary.csv"
    print(f"Reading {csv_path} ...")
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False, na_values=[""])
    # Apply proper dtypes
    for col, dtype in SUMMARY_DTYPES.items():
        if col in df.columns:
            if dtype == "boolean":
                df[col] = df[col].map({"TRUE": True, "FALSE": False, "True": True, "False": False}).astype("boolean")
            elif dtype in ("Int32", "Int8"):
                df[col] = pd.to_numeric(df[col], errors="coerce").astype(dtype)
            else:
                df[col] = df[col].astype("string")
    print(f"  {len(df):,} rows")
    write_parquet(df, out_dir / f"{sp}_summary.parquet")

    # ── Intervals ────────────────────────────────────────────────────────────
    csv_path = CSV_DIR / f"{sp}_intervals.csv"
    print(f"Reading {csv_path} ...")
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False, na_values=[""])
    for col, dtype in INTERVALS_DTYPES.items():
        if col in df.columns:
            if dtype == "boolean":
                df[col] = df[col].map({"TRUE": True, "FALSE": False, "True": True, "False": False}).astype("boolean")
            elif dtype in ("Int32", "Int8"):
                df[col] = pd.to_numeric(df[col], errors="coerce").astype(dtype)
            else:
                df[col] = df[col].astype("string")
    print(f"  {len(df):,} rows")
    write_parquet(df, out_dir / f"{sp}_intervals.parquet")


def verify_parquet(sp: str) -> None:
    print(f"\n--- Verify {sp} ---")
    out_dir = DATA_DIR / sp

    summary_path   = out_dir / f"{sp}_summary.parquet"
    intervals_path = out_dir / f"{sp}_intervals.parquet"

    s = pd.read_parquet(summary_path, columns=["gene_id", "gene_name", "constitutive_exist", "constitutive_length", "tsl_cutoff"])
    print(f"  summary: {len(s):,} rows, columns: {list(s.columns)}")
    print(f"  tsl_cutoff values: {s['tsl_cutoff'].value_counts().head(5).to_dict()}")

    i = pd.read_parquet(intervals_path, columns=["gene_id", "feature", "transcript_support_level", "tsl_included"])
    print(f"  intervals: {len(i):,} rows, feature counts: {i['feature'].value_counts().to_dict()}")

    # Test predicate-pushdown: read one gene
    sample_gene = s["gene_id"].iloc[0]
    filters = [("gene_id", "==", sample_gene)]
    g = pd.read_parquet(intervals_path, filters=filters)
    print(f"  predicate-pushdown test ({sample_gene}): {len(g)} interval rows")


if __name__ == "__main__":
    for sp in SPECIES:
        convert_species(sp)
    for sp in SPECIES:
        verify_parquet(sp)
    print("\n=== All done. Parquet files written to data/isoform_explorer/ ===")
    print("After verifying, remove source files:")
    print("  rm /tmp/isoform_explorer_csv/*.csv")
    for sp in SPECIES:
        print(f"  rm data/isoform_explorer/{sp}/*.qs data/isoform_explorer/{sp}/*_summary.csv")
