"""
shapece.reactivity
==================
Turn per-peak areas into per-nucleotide reactivity. Two models are provided.

* ``area_difference`` (**default**, QuShape-compatible): scale the (-) background
  to the (+) reagent on the low-reactivity peaks, subtract, then boxplot-normalize.
  This is the long-standing Weeks-lab CE-SHAPE quantity, for continuity with the
  field.

* ``stop_fraction`` (alternative): convert each lane's areas to a per-nucleotide
  RT-stop probability, ``-ln(1 - f)`` (Aviran/Weeks), which is loading-independent
  and depletion-corrected, then subtract (+)-(-). Useful for MaP-like or
  stop-based readouts.

Normalization is the model-free 2-8% boxplot rule (Low & Weeks).
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize_scalar

trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def scale_background(area_rx: np.ndarray, area_bg: np.ndarray, rate: float = 0.25) -> float:
    """Find the factor that scales BG onto RX using the lowest-reactivity peaks.

    QuShape's logic: at unreactive positions RX ~ BG, so the bottom ``rate``
    fraction of RX peaks (paired with their BG) fixes the loading ratio without
    being biased by genuinely reactive positions. Returns the scale factor.
    """
    n = len(area_rx)
    order = np.argsort(area_rx)[:max(3, int(n * rate))]
    A, B = area_rx[order], area_bg[order]
    res = minimize_scalar(lambda f: np.abs(A - f * B).sum(),
                          bounds=(0.1, 10.0), method="bounded")
    return float(res.x) if res.success else 1.0


def area_difference(area_rx: np.ndarray, area_bg: np.ndarray,
                    scale: bool = True, rate: float = 0.25,
                    clip: bool = False) -> np.ndarray:
    """Reactivity = area(RX) - factor*area(BG).

    ``clip=False`` (default) matches QuShape, which retains small negative values
    (they carry noise information and affect normalization). Set ``clip=True`` to
    floor at zero, which some downstream folding tools prefer.
    """
    factor = scale_background(area_rx, area_bg, rate) if scale else 1.0
    diff = area_rx - factor * area_bg
    return np.clip(diff, 0, None) if clip else diff


def stop_fraction(trace: np.ndarray, positions: np.ndarray, half_width: float = 6.0,
                  end: int | None = None) -> np.ndarray:
    """Per-lane RT-stop reactivity ``-ln(1 - f)`` from one aligned trace.

    ``f`` = (area at a position) / (area from that position to ``end``), i.e. the
    fraction of molecules reaching a nucleotide that stopped there. Loading-
    independent and depletion-corrected.
    """
    end = len(trace) if end is None else end
    react = np.empty(len(positions))
    for i, p in enumerate(positions):
        p = int(p)
        num = trapezoid(trace[max(0, p - int(half_width)):p + int(half_width) + 1])
        den = trapezoid(trace[max(0, p - int(half_width)):end])
        f = num / den if den > 0 else 0.0
        react[i] = -np.log(1.0 - min(f, 0.999))
    return react


def _find_p_outlier_box(data: np.ndarray) -> tuple[float, float]:
    """Verbatim port of QuShape's ``findPOutlierBox``.

    Outliers are values above Q3 + 1.5*IQR (traditional boxplot rule). For fewer
    than 100 data points the "top 10%" becomes the *top 10 values*.
    Operates on **all** values (including negatives), exactly as QuShape does.
    """
    n = len(data)
    if n < 50:
        return 2.0, 10.0
    ds = np.sort(data)
    q1 = ds[int(n * 0.25)]
    q3 = ds[int(n * 0.75)]
    qthres = 1.5 * (q3 - q1) + q3
    n_out = 0
    for i in range(n - 1, 0, -1):
        if ds[i] > qthres:
            n_out += 1
        else:
            break
    p_out = float(n_out) / float(n) * 100.0
    p_aver = 10.0
    if n < 100:
        p_aver = (10.0 / float(n)) * 100.0
    return p_out, p_aver


def _norm_simple(data: np.ndarray, p_outlier: float = 2.0,
                 p_aver: float = 10.0) -> tuple[np.ndarray, float]:
    """Verbatim port of QuShape's ``normSimple``: divide by the mean of the top
    ``p_aver`` percent of values, after discarding the top ``p_outlier`` percent."""
    n = len(data)
    n_out = int(float(n) * float(p_outlier) / 100.0)
    if n_out < 1:
        n_out = 1
    n_aver = int(float(n) * float(p_aver) / 100.0) + n_out
    ds = np.sort(data)
    aver = float(np.average(ds[-n_aver:-n_out]))
    return data / aver, aver


def boxplot_normalize(reactivity: np.ndarray) -> tuple[np.ndarray, float]:
    """Model-free boxplot / 2-8% normalization (Low & Weeks), QuShape-compatible.

    Faithful port of QuShape's ``normBox``: the boxplot outlier rule and the
    top-10%/top-10-values average are computed over **all** reactivity values
    (negatives included, no clipping), so results are directly comparable to
    QuShape output. Returns ``(normalized, factor)``.
    """
    r = np.asarray(reactivity, float)
    if len(r) < 5:
        return r, 1.0
    p_out, p_aver = _find_p_outlier_box(r)
    normed, factor = _norm_simple(r, p_out, p_aver)
    return (normed if factor > 0 else r), factor
