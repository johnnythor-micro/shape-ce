"""
shapece.plots
=============
Publication-oriented visualizations for reactivity profiles, differences, and
structures. Each function returns a matplotlib ``Axes`` so plots compose in a
notebook and stay editable for figure preparation.

Choosing a plot
---------------
============  ==========================================================
Plot          What it shows / when to use it
============  ==========================================================
skyline       Reactivity as a step function, one step per nucleotide.
              Steps span the full nucleotide width, so **two or more
              profiles overlay without occluding each other**. The
              standard way to compare conditions or replicates.
profile       A single condition as colored bars (SHAPE convention) with
              error bars. Use to **report one dataset**: reactive vs.
              protected nucleotides at a glance.
delta_shape   Smoothed difference between conditions with significant
              sites shaded, optionally with Z-factor / standard-score
              tracks. The signature deltaSHAPE output.
arc           Base pairs (or any pair list / probability matrix) drawn as
              arcs over a sequence axis; a second dataset can be mirrored
              below. Nested arcs = helices; crossing arcs = pseudoknots.
              Best for **short to medium RNAs**.
circle        Sequence wrapped on a circle, pairs drawn as chords. Arcs
              become unreadably tall on long RNAs; chords do not, so use
              this for **long RNAs and long-range contacts**.
linreg        Replicate-vs-replicate scatter with R². The reproducibility
              check that should **precede** any comparison.
heatmap       Replicates x nucleotides. Use to **spot a bad lane**.
============  ==========================================================

Color convention (shared with :mod:`shapece.structure`): reactivity >= 0.85 red
(flexible / likely unpaired), 0.40-0.85 orange (intermediate), < 0.40 black
(constrained / likely paired), no data gray.
"""
from __future__ import annotations

import numpy as np

__all__ = ["shape_colors", "skyline", "profile", "delta_shape_plot",
           "arc", "circle", "linreg", "heatmap"]

HIGH = "#d62728"
MED = "#ff7f0e"
LOW = "#1a1a1a"
NODATA = "#cccccc"


def _ax(ax, figsize):
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    return ax


def shape_colors(reactivity) -> list:
    """Map reactivities to the standard SHAPE color convention."""
    out = []
    for r in np.asarray(reactivity, float):
        if not np.isfinite(r):
            out.append(NODATA)
        elif r >= 0.85:
            out.append(HIGH)
        elif r >= 0.40:
            out.append(MED)
        else:
            out.append(LOW)
    return out


def _pairs_from_dotbracket(db: str):
    stack, pairs = [], []
    for i, c in enumerate(db, 1):
        if c == "(":
            stack.append(i)
        elif c == ")":
            pairs.append((stack.pop(), i))
    return pairs


def _as_pairs(structure, threshold: float = 0.3):
    """Accept a dot-bracket string, a list of (i, j) pairs, or a probability
    matrix, and return a list of ``(i, j, weight)`` with 1-based indices."""
    if isinstance(structure, str):
        return [(i, j, 1.0) for i, j in _pairs_from_dotbracket(structure)]
    arr = np.asarray(structure)
    if arr.ndim == 2 and arr.shape[0] == arr.shape[1] and arr.shape[0] > 2:
        out = []
        n = arr.shape[0]
        for i in range(n):
            for j in range(i + 1, n):
                w = arr[i, j]
                if np.isfinite(w) and w >= threshold:
                    out.append((i + 1, j + 1, float(w)))
        return out
    return [(int(p[0]), int(p[1]), float(p[2]) if len(p) > 2 else 1.0)
            for p in structure]


# ---------------------------------------------------------------------------
def skyline(profiles: dict, nucleotides=None, ax=None, colors=None,
            figsize=(14, 3.5), title: str = "", ylabel: str = "normalized reactivity"):
    """Overlay reactivity profiles as step ("skyline") plots.

    Parameters
    ----------
    profiles : dict[str, array-like]
        label -> per-nucleotide reactivity. Overlaid on one axis.
    nucleotides : array-like, optional
        Nucleotide numbers (defaults to 1..N).

    Why a skyline? Each nucleotide is drawn as a full-width step rather than a
    thin bar, so several conditions can be superimposed and compared directly
    without hiding one another.
    """
    ax = _ax(ax, figsize)
    for k, (label, prof) in enumerate(profiles.items()):
        y = np.asarray(prof, float)
        x = np.arange(1, len(y) + 1) if nucleotides is None else np.asarray(nucleotides)
        c = None if colors is None else colors[k % len(colors)]
        ax.step(x, y, where="mid", lw=1.0, label=label, color=c)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("nucleotide")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if len(profiles) > 1:
        ax.legend(fontsize=8)
    return ax


def profile(reactivity, error=None, nucleotides=None, ax=None,
            figsize=(14, 3.5), title: str = ""):
    """Single-condition bar plot, colored by the SHAPE convention.

    Error bars, when supplied, are typically the SEM across replicates.
    """
    ax = _ax(ax, figsize)
    y = np.asarray(reactivity, float)
    x = np.arange(1, len(y) + 1) if nucleotides is None else np.asarray(nucleotides)
    ax.bar(x, np.nan_to_num(y), color=shape_colors(y),
           yerr=None if error is None else np.asarray(error, float),
           ecolor="0.6", error_kw={"lw": 0.6})
    ax.axhline(1.0, color="k", lw=0.4, ls="--")
    ax.set_xlabel("nucleotide")
    ax.set_ylabel("normalized reactivity")
    if title:
        ax.set_title(title)
    return ax


def delta_shape_plot(result: dict, nucleotides=None, ax=None, figsize=(14, 4),
                     title: str = "", show_dots: bool = True,
                     labels=("condition 1", "condition 2")):
    """Plot a :func:`shapece.stats.delta_shape` result.

    Shows the smoothed difference (condition 1 - condition 2). Significant sites
    are shaded: **orange** where condition 1 is more reactive, **blue** where
    condition 2 is. Optional dots mark nucleotides passing the Z-factor filter.

    A positive region means the RNA became *more* flexible/accessible in
    condition 1 -- e.g. a helix melting or a protein leaving.
    """
    ax = _ax(ax, figsize)
    sd = np.asarray(result["smoothed_diff"], float)
    x = np.arange(1, len(sd) + 1) if nucleotides is None else np.asarray(nucleotides)
    ax.step(x, sd, where="mid", color="0.25", lw=0.9)
    ax.axhline(0, color="k", lw=0.6)

    sig = result["significant"]
    pos = sig & (sd >= 0)
    neg = sig & (sd < 0)
    ax.fill_between(x, 0, np.where(pos, sd, 0), step="mid", color="#e67e22", alpha=0.85,
                    label=f"up in {labels[0]}")
    ax.fill_between(x, 0, np.where(neg, sd, 0), step="mid", color="#2471a3", alpha=0.85,
                    label=f"up in {labels[1]}")
    if show_dots:
        zf = np.asarray(result["z_factors"], float)
        ok = np.isfinite(zf) & (zf > 0)
        yy = np.nanmax(np.abs(sd)) * 1.15 if np.isfinite(sd).any() else 1.0
        ax.plot(x[ok], np.full(ok.sum(), yy), ".", ms=3, color="#555",
                label="Z-factor > 0")
    ax.set_xlabel("nucleotide")
    ax.set_ylabel("Δ reactivity (smoothed)")
    ax.set_title(title or f"deltaSHAPE: {result['n_sites']} site(s) called")
    ax.legend(fontsize=7, loc="upper right")
    return ax


def arc(structure, reactivity=None, sequence=None, structure2=None, ax=None,
        figsize=(14, 5), title: str = "", threshold: float = 0.3,
        labels=("", "")):
    """Arc diagram of base pairs, optionally comparing two structures.

    ``structure`` may be a dot-bracket string, a list of ``(i, j)`` pairs, or a
    square probability matrix (pairs with weight >= ``threshold`` are drawn, with
    opacity scaled by weight).

    If ``structure2`` is given it is mirrored **below** the axis, so two models
    (e.g. SHAPE-directed vs thermodynamic) can be compared arc-for-arc. Shared
    pairs appear as mirror images; disagreements stand out immediately.

    Nested arcs correspond to helices; crossing arcs indicate pseudoknots.
    """
    ax = _ax(ax, figsize)
    pairs = _as_pairs(structure, threshold)
    n = (len(sequence) if sequence is not None
         else max([j for _, j, _ in pairs], default=1))

    def _draw(prs, sign, color):
        for i, j, w in prs:
            c = (i + j) / 2.0
            r = (j - i) / 2.0
            th = np.linspace(0, np.pi, 60)
            ax.plot(c + r * np.cos(th), sign * r * np.sin(th),
                    color=color, lw=0.8, alpha=min(1.0, 0.25 + 0.75 * w))

    _draw(pairs, 1, "#9ecae1")
    if structure2 is not None:
        _draw(_as_pairs(structure2, threshold), -1, "#fdae6b")

    if reactivity is not None:
        ax.scatter(np.arange(1, n + 1), np.zeros(n),
                   c=shape_colors(np.asarray(reactivity, float)[:n]), s=12, marker="s",
                   zorder=3)
    else:
        ax.plot([1, n], [0, 0], color="0.6", lw=1)

    top = max((j - i) / 2.0 for i, j, _ in pairs) if pairs else 1
    ax.set_ylim(-(top + 2) if structure2 is not None else -2, top + 2)
    ax.set_yticks([])
    ax.set_xlabel("nucleotide")
    if labels[0]:
        ax.text(0.01, 0.95, labels[0], transform=ax.transAxes, fontsize=8, va="top")
    if labels[1]:
        ax.text(0.01, 0.05, labels[1], transform=ax.transAxes, fontsize=8, va="bottom")
    if title:
        ax.set_title(title)
    return ax


def circle(structure, reactivity=None, sequence_length=None, ax=None,
           figsize=(7.5, 7.5), title: str = "", threshold: float = 0.3):
    """Circle (circos-style) plot: sequence on a circle, pairs as chords.

    On long RNAs an arc diagram becomes unreadable because arc height grows with
    pairing distance. A circle plot keeps every contact the same visual scale, so
    **long-range interactions are easy to see**. Chord opacity scales with pair
    weight (e.g. base-pair probability).

    Nucleotides are drawn around the rim, colored by reactivity if supplied.
    """
    ax = _ax(ax, figsize)
    pairs = _as_pairs(structure, threshold)
    n = (sequence_length if sequence_length is not None
         else max([j for _, j, _ in pairs], default=1))
    ang = lambda k: 2 * np.pi * (k - 1) / n + np.pi / 2  # noqa: E731
    xy = lambda k: (np.cos(ang(k)), np.sin(ang(k)))      # noqa: E731

    for i, j, w in pairs:
        x1, y1 = xy(i)
        x2, y2 = xy(j)
        # quadratic Bezier bowed toward the centre; deeper bow for longer range
        sep = min(abs(j - i), n - abs(j - i)) / (n / 2.0)
        cx, cy = (x1 + x2) / 2 * (1 - sep), (y1 + y2) / 2 * (1 - sep)
        t = np.linspace(0, 1, 60)
        bx = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * cx + t ** 2 * x2
        by = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * cy + t ** 2 * y2
        ax.plot(bx, by, color="#6baed6", lw=0.7, alpha=min(1.0, 0.25 + 0.75 * w))

    ks = np.arange(1, n + 1)
    px, py = np.cos(ang(ks)), np.sin(ang(ks))
    cols = (["#999999"] * n if reactivity is None
            else shape_colors(np.asarray(reactivity, float)[:n]))
    ax.scatter(px, py, c=cols, s=8, zorder=3)
    for k in range(0, n, max(1, n // 12)):
        ax.text(px[k] * 1.08, py[k] * 1.08, str(k + 1), ha="center", va="center", fontsize=6)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title)
    return ax


def linreg(rep1, rep2, ax=None, figsize=(4.5, 4.5), labels=("replicate 1", "replicate 2"),
           title: str = ""):
    """Replicate scatter with a fitted line and R^2 -- the reproducibility check.

    Run this **before** any between-condition comparison. Points should hug the
    diagonal; a low R^2 means the comparison that follows is not interpretable.
    """
    from scipy import stats as _st
    ax = _ax(ax, figsize)
    a = np.asarray(rep1, float)
    b = np.asarray(rep2, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    ax.scatter(a, b, s=10, alpha=0.7, color="#2471a3")
    if len(a) > 2:
        sl, ic, r, _, _ = _st.linregress(a, b)
        xs = np.linspace(a.min(), a.max(), 10)
        ax.plot(xs, sl * xs + ic, "k--", lw=0.8)
        ax.set_title(title or f"R² = {r**2:.3f}")
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    return ax


def heatmap(matrix, row_labels=None, nucleotides=None, ax=None, figsize=(14, 3),
            title: str = "", cmap: str = "magma", vmin=None, vmax=None):
    """Replicates x nucleotides heatmap -- use it to spot an outlier lane."""
    ax = _ax(ax, figsize)
    m = np.asarray(matrix, float)
    extent = None
    if nucleotides is not None:
        nt = np.asarray(nucleotides)
        extent = [nt.min(), nt.max(), m.shape[0] - 0.5, -0.5]
    im = ax.imshow(m, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, extent=extent)
    if row_labels is not None:
        ax.set_yticks(range(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xlabel("nucleotide")
    if title:
        ax.set_title(title)
    ax.figure.colorbar(im, ax=ax, fraction=0.025, label="reactivity")
    return ax
