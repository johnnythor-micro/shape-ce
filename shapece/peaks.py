"""
shapece.peaks
=============
Peak detection and quantification. Two quantification modes are provided:

* ``"gaussian"`` (default, QuShape-style): iteratively fit a Gaussian to each
  peak and report area = amplitude x width. Best for overlapping peaks.
* ``"trapezoid"``: simple numerical integration in a fixed window. Fast and
  robust when peaks are well separated.

Peaks are detected on a chosen reference trace and then the *same* positions are
quantified in every lane, so lanes stay row-comparable.
"""
from __future__ import annotations
import numpy as np
from scipy.signal import find_peaks

trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def detect_peaks(reference: np.ndarray, lo: int = 0, hi: int | None = None,
                 min_spacing: float = 12.0, prominence_frac: float = 0.03) -> np.ndarray:
    """Detect consensus peak positions (sample indices) on a reference trace.

    Parameters
    ----------
    reference : np.ndarray
        Trace to detect peaks on (typically the mean of the (+) lanes, which has
        the richest peak structure).
    lo, hi : int
        Analysis window (sample indices). ``hi=None`` -> end of trace.
    min_spacing : float
        Minimum samples between peaks (~ one nucleotide). Prevents splitting.
    prominence_frac : float
        Peak prominence as a fraction of the in-window maximum.
    """
    hi = len(reference) if hi is None else hi
    win = np.zeros(len(reference), bool)
    win[lo:hi] = True
    # robust scale: use the 98th percentile so a saturated primer spike does not
    # set the threshold and suppress genuine peaks.
    robust_max = np.percentile(reference[win], 98)
    prom = prominence_frac * robust_max
    idx, _ = find_peaks(np.where(win, reference, 0.0),
                        prominence=prom, distance=int(round(min_spacing)))
    return idx


def _gaussian(x, pos, amp, wid):
    return amp * np.exp(-2 * (x - pos) ** 2 / wid ** 2)


def _fit_widths(trace, positions, init_wid):
    """Estimate a per-peak Gaussian width by minimizing local reconstruction error."""
    wid = np.full(len(positions), init_wid, float)
    for i in range(1, len(positions) - 1):
        lo, hi = positions[i - 1], positions[i + 1] + 1
        x = np.arange(lo, hi)
        y = trace[lo:hi]
        best_err, best_w = np.inf, init_wid
        for w in np.linspace(init_wid * 0.7, init_wid * 1.3, 9):
            model = (_gaussian(x, positions[i - 1], trace[positions[i - 1]], wid[i - 1])
                     + _gaussian(x, positions[i], trace[positions[i]], w)
                     + _gaussian(x, positions[i + 1], trace[positions[i + 1]], wid[i + 1]))
            err = np.abs(model - y).sum()
            if err < best_err:
                best_err, best_w = err, w
        wid[i] = best_w
    return wid


def quantify(trace: np.ndarray, positions: np.ndarray, mode: str = "gaussian",
             half_width: float = 6.0) -> dict:
    """Quantify each peak in ``trace`` at the given ``positions``.

    Returns a dict with arrays ``pos``, ``amp``, ``wid``, ``area`` (one per peak).
    ``area`` is the quantity reactivity is computed from.
    """
    positions = np.asarray(positions, int)
    amp = trace[positions].astype(float)
    if mode == "gaussian":
        spacing = np.median(np.diff(positions)) if len(positions) > 1 else 2 * half_width
        wid = _fit_widths(trace, positions, init_wid=0.45 * spacing)
        area = np.abs(amp * wid)
    elif mode == "trapezoid":
        wid = np.full(len(positions), half_width, float)
        area = np.empty(len(positions))
        for i, p in enumerate(positions):
            lo = max(0, int(p - half_width)); hi = min(len(trace), int(p + half_width) + 1)
            area[i] = trapezoid(trace[lo:hi], np.arange(lo, hi))
    else:
        raise ValueError("mode must be 'gaussian' or 'trapezoid'")
    return {"pos": positions.astype(float), "amp": amp, "wid": wid, "area": area}
