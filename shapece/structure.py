"""
shapece.structure
=================
SHAPE-directed secondary-structure prediction with ViennaRNA, plus simple 2D
plots. Wraps the Deigan (pseudo-energy) and Zarringhalam (probability) SHAPE
methods. Requires the ``ViennaRNA`` package (``pip install ViennaRNA``).
"""
from __future__ import annotations
import numpy as np


def _to_vienna_shape(reactivity: np.ndarray, seq_len: int) -> list:
    """Build a 1-indexed reactivity list for ViennaRNA (index 0 unused; -999=no data)."""
    r = [-999.0] * (seq_len + 1)
    for i, v in enumerate(reactivity):
        if i + 1 <= seq_len and np.isfinite(v):
            r[i + 1] = float(v)
    return r


def fold(seq: str, reactivity=None, method: str = "deigan",
         slope: float = 1.8, intercept: float = -0.6, temperature: float = 37.0):
    """Fold ``seq`` (optionally SHAPE-directed). Returns (dot_bracket, mfe).

    Parameters
    ----------
    seq : str
        RNA sequence (U or T accepted).
    reactivity : array-like or None
        Per-nucleotide reactivity aligned to ``seq`` (1:1, use NaN where no data).
        ``None`` -> unconstrained thermodynamic fold.
    method : {"deigan", "zarringhalam"}
        SHAPE incorporation method. Deigan uses (slope, intercept) pseudo-energies.
    slope, intercept : float
        Deigan m and b (defaults 1.8, -0.6, the standard values).
    temperature : float
        Folding temperature in Celsius.
    """
    import RNA
    md = RNA.md(); md.temperature = temperature
    fc = RNA.fold_compound(seq.replace("T", "U"), md)
    if reactivity is not None:
        shp = _to_vienna_shape(np.asarray(reactivity, float), len(seq))
        if method == "deigan":
            fc.sc_add_SHAPE_deigan(shp, slope, intercept)
        elif method == "zarringhalam":
            fc.sc_add_SHAPE_zarringham(shp, 0.5, 0.8, "O")
        else:
            raise ValueError("method must be 'deigan' or 'zarringhalam'")
    ss, mfe = fc.mfe()
    return ss, float(mfe)


def base_pair_probabilities(seq: str, reactivity=None, **kw):
    """Return the base-pair probability matrix from the partition function."""
    import RNA
    md = RNA.md(); md.temperature = kw.get("temperature", 37.0)
    fc = RNA.fold_compound(seq.replace("T", "U"), md)
    if reactivity is not None:
        shp = _to_vienna_shape(np.asarray(reactivity, float), len(seq))
        fc.sc_add_SHAPE_deigan(shp, kw.get("slope", 1.8), kw.get("intercept", -0.6))
    fc.pf()
    return np.array(fc.bpp())


def pair_table(dot_bracket: str) -> list:
    """1-based pair table: pt[i]=j if i pairs j, else 0 (index 0 unused)."""
    pt = [0] * (len(dot_bracket) + 1)
    stack = []
    for i, c in enumerate(dot_bracket, 1):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop(); pt[i] = j; pt[j] = i
    return pt


def shape_color(r) -> str:
    """SHAPE reactivity color convention."""
    if r is None or not np.isfinite(r):
        return "#cccccc"
    if r >= 0.85:
        return "#d62728"    # high  (flexible / likely unpaired)
    if r >= 0.40:
        return "#ff7f0e"    # medium
    return "#1a1a1a"        # low   (likely paired)


def plot_arc(seq: str, dot_bracket: str, reactivity=None, ax=None, title=""):
    """Arc diagram colored by reactivity (no external layout dependency)."""
    import matplotlib.pyplot as plt
    pt = pair_table(dot_bracket)
    n = len(seq)
    if ax is None:
        _, ax = plt.subplots(figsize=(min(20, n / 8 + 3), 3.5))
    for i in range(1, n + 1):
        if pt[i] > i:
            j = pt[i]; c = (i + j) / 2; rad = (j - i) / 2
            th = np.linspace(0, np.pi, 50)
            ax.plot(c + rad * np.cos(th), rad * np.sin(th), color="#9ecae1", lw=0.6)
    cols = (["#cccccc"] * n if reactivity is None
            else [shape_color(v) for v in reactivity])
    ax.scatter(np.arange(1, n + 1), np.zeros(n), c=cols, s=12, marker="s")
    ax.set_ylim(-2, n / 2 + 2); ax.set_yticks([]); ax.set_xlabel("nucleotide")
    ax.set_title(title)
    return ax


def write_ct(path: str, seq: str, dot_bracket: str, title: str = "structure"):
    """Write a .ct connectivity file (RNAstructure/compatible)."""
    pt = pair_table(dot_bracket)
    n = len(seq)
    with open(path, "w") as fh:
        fh.write(f"{n}\t{title}\n")
        for i in range(1, n + 1):
            fh.write(f"{i}\t{seq[i-1]}\t{i-1}\t{i+1 if i<n else 0}\t{pt[i]}\t{i}\n")
