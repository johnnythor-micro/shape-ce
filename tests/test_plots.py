"""Tests for `shapece.plots`: every plotter returns an Axes and tolerates no-data."""
import matplotlib
matplotlib.use("Agg")

import numpy as np
from matplotlib.axes import Axes

from shapece import plots as P, stats as S


def _profiles(n=60):
    rng = np.random.RandomState(1)
    a = np.abs(rng.normal(0.4, 0.3, n))
    b = a.copy(); b[25:31] += 1.0
    a[10] = np.nan; b[10] = np.nan          # a no-data nucleotide
    return a, b


DB = "".join("(" if 5 <= i < 15 else (")" if 25 <= i < 35 else ".") for i in range(60))


def test_shape_colors_convention():
    cols = P.shape_colors([0.1, 0.5, 1.2, np.nan])
    assert cols == [P.LOW, P.MED, P.HIGH, P.NODATA]


def test_skyline_overlays_multiple_profiles():
    a, b = _profiles()
    ax = P.skyline({"cond1": a, "cond2": b})
    assert isinstance(ax, Axes)
    assert len(ax.lines) >= 2          # one step line per profile (+ zero line)


def test_profile_with_errors():
    a, _ = _profiles()
    ax = P.profile(a, error=np.full(len(a), 0.05))
    assert isinstance(ax, Axes) and len(ax.patches) == len(a)


def test_delta_shape_plot():
    a, b = _profiles()
    e = np.full(len(a), 0.05)
    res = S.delta_shape(a, b, e, e.copy(), pad=1)
    ax = P.delta_shape_plot(res, labels=("A", "B"))
    assert isinstance(ax, Axes)


def test_arc_accepts_dotbracket_and_second_structure():
    a, _ = _profiles()
    ax = P.arc(DB, reactivity=a, sequence="N" * 60, structure2=DB)
    assert isinstance(ax, Axes)
    # arcs drawn above and below the axis
    ys = np.concatenate([ln.get_ydata() for ln in ax.lines])
    assert ys.max() > 0 and ys.min() < 0


def test_arc_accepts_pair_list():
    ax = P.arc([(5, 30), (6, 29)], sequence="N" * 60)
    assert isinstance(ax, Axes)


def test_circle_accepts_probability_matrix():
    n = 60
    probs = np.zeros((n, n))
    for i, j in zip(range(5, 15), range(34, 24, -1)):
        probs[i, j] = 0.9
    ax = P.circle(probs, sequence_length=n)
    assert isinstance(ax, Axes) and len(ax.lines) == 10   # one chord per pair


def test_circle_threshold_filters_weak_pairs():
    n = 40
    probs = np.zeros((n, n))
    probs[3, 30] = 0.9
    probs[4, 29] = 0.1                # below default threshold
    ax = P.circle(probs, sequence_length=n, threshold=0.3)
    assert len(ax.lines) == 1


def test_linreg_returns_axes():
    a, b = _profiles()
    ax = P.linreg(a, b)
    assert isinstance(ax, Axes) and "R²" in ax.get_title()


def test_heatmap_returns_axes():
    a, b = _profiles()
    ax = P.heatmap(np.vstack([a, b]), row_labels=["r1", "r2"])
    assert isinstance(ax, Axes)
