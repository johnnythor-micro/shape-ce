"""
shapece.report
==============
Export reactivities to the formats the field uses:

* ``.shape`` -- two columns (position, reactivity); -999 = no data.
* ``.map``   -- RNAstructure format (position, reactivity, std-error, nucleotide).
* ``.xlsx``  -- a tidy workbook designed to drop straight into GraphPad Prism
  (one row per nucleotide, one column per replicate, plus mean/SD), so users can
  regenerate publication graphs and run their own statistics in Prism.

Nothing here assumes a particular experiment; you pass in whatever samples and
conditions you have.
"""
from __future__ import annotations
import numpy as np
import pandas as pd


def reactivity_table(nucleotide, base, samples: dict) -> pd.DataFrame:
    """Build a tidy per-nucleotide table.

    Parameters
    ----------
    nucleotide : array-like
        Nucleotide numbers.
    base : array-like or str
        Per-nucleotide base identity (or the sequence string).
    samples : dict[str, array-like]
        Column label -> reactivity array (one entry per replicate/condition),
        e.g. ``{"cond1_rep1": r1, "cond1_rep2": r2, "cond2_rep1": r3}``.
    """
    if isinstance(base, str):
        base = list(base)
    df = pd.DataFrame({"nucleotide": np.asarray(nucleotide, int),
                       "base": list(base)[:len(nucleotide)]})
    for name, vals in samples.items():
        df[name] = np.asarray(vals, float)
    return df


def add_group_stats(df: pd.DataFrame, groups: dict) -> pd.DataFrame:
    """Append mean and SD columns for named groups of replicate columns.

    ``groups`` maps a group label to the list of replicate column names, e.g.
    ``{"cond1": ["cond1_rep1","cond1_rep2","cond1_rep3"]}``.
    """
    out = df.copy()
    for label, cols in groups.items():
        out[f"{label}_mean"] = out[cols].mean(axis=1)
        out[f"{label}_SD"] = out[cols].std(axis=1, ddof=1)
    return out


def write_shape(path: str, nucleotide, reactivity, seq_len: int | None = None):
    """Write a ``.shape`` file (1-based, -999 for no data / gaps)."""
    seq_len = int(np.nanmax(nucleotide)) if seq_len is None else seq_len
    arr = np.full(seq_len, -999.0)
    for nt, v in zip(nucleotide, reactivity):
        if np.isfinite(nt) and 1 <= int(nt) <= seq_len and np.isfinite(v):
            arr[int(nt) - 1] = v
    with open(path, "w") as fh:
        for i in range(seq_len):
            fh.write(f"{i+1}\t{arr[i]:.4f}\n")


def write_map(path: str, nucleotide, reactivity, seq: str, stderr=None):
    """Write an RNAstructure ``.map`` file (position, reactivity, stderr, base)."""
    seq_len = len(seq)
    react = np.full(seq_len, -999.0); err = np.zeros(seq_len)
    for k, nt in enumerate(nucleotide):
        if np.isfinite(nt) and 1 <= int(nt) <= seq_len:
            react[int(nt) - 1] = reactivity[k]
            if stderr is not None:
                err[int(nt) - 1] = stderr[k]
    with open(path, "w") as fh:
        for i in range(seq_len):
            fh.write(f"{i+1}\t{react[i]:.4f}\t{err[i]:.4f}\t{seq[i]}\n")


def write_excel_for_prism(path: str, df: pd.DataFrame, notes: str | None = None):
    """Write a Prism-friendly ``.xlsx``.

    Sheet 'reactivity' is a wide table (rows = nucleotides, columns = samples),
    which imports into a Prism 'Column' or 'XY' table directly. Sheet 'README'
    documents the columns and how to import.
    """
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        df.to_excel(xw, sheet_name="reactivity", index=False)
        readme = pd.DataFrame({"how_to_use": [
            "Sheet 'reactivity': one row per nucleotide.",
            "Column 'nucleotide' = position; 'base' = identity.",
            "Remaining columns = per-replicate normalized reactivity.",
            "",
            "GraphPad Prism import:",
            "1. New table > XY (or Grouped).",
            "2. Copy 'nucleotide' into X; copy replicate columns into Y groups.",
            "3. For mean+/-SD per condition, group replicate columns as sub-columns.",
            "",
            (notes or "")]})
        readme.to_excel(xw, sheet_name="README", index=False)
