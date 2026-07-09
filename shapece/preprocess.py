"""
shapece.preprocess
==================
Signal conditioning for CE traces, modernizing the corrections QuShape applies
in its "Tool Data" stage. Every function is pure (returns a new array) and
individually optional, so a user can build exactly the preprocessing chain their
data needs.

Chain (typical order): saturation -> smoothing -> baseline -> (mobility shift) ->
(decay correction).
"""
from __future__ import annotations
import numpy as np
from scipy.signal import savgol_filter
from scipy.ndimage import grey_erosion, grey_dilation
from scipy.optimize import curve_fit


def correct_saturation(trace: np.ndarray, threshold: float | None = None) -> np.ndarray:
    """Repair off-scale (clipped) peaks by cubic interpolation across the plateau.

    CE detectors saturate at a hardware ceiling; the true peak apex is lost. We
    detect runs at/above ``threshold`` and interpolate across them.

    Parameters
    ----------
    trace : np.ndarray
    threshold : float, optional
        Saturation level. Defaults to 0.999 * observed max.
    """
    y = trace.astype(float).copy()
    if threshold is None:
        threshold = 0.999 * y.max()
    sat = y >= threshold
    if not sat.any():
        return y
    good = ~sat
    x = np.arange(len(y))
    y[sat] = np.interp(x[sat], x[good], y[good])   # linear fill is safe & monotone
    return y


def smooth(trace: np.ndarray, window: int = 11, polyorder: int = 3,
           method: str = "savgol") -> np.ndarray:
    """Smooth a trace while preserving peak area.

    ``method="savgol"`` (default) uses a Savitzky-Golay filter (recommended: it
    preserves peak height/area far better than a moving average).
    ``method="triangle"`` reproduces QuShape's triangular smoothing.
    """
    if method == "savgol":
        if window % 2 == 0:
            window += 1
        return savgol_filter(trace, window, polyorder)
    if method == "triangle":
        w = np.arange(1, window // 2 + 2)
        w = np.concatenate([w, w[-2::-1]]).astype(float)
        w /= w.sum()
        return np.convolve(trace, w, mode="same")
    raise ValueError("method must be 'savgol' or 'triangle'")


def baseline(trace: np.ndarray, window: int = 150) -> np.ndarray:
    """Estimate and subtract a slowly varying baseline via morphological opening.

    An opening (erosion then dilation) with a structuring element much wider than
    a peak follows the floor under the peaks without being pulled up by them —
    more robust than polynomial baselines against tall peaks. Returns the
    baseline-subtracted trace (clipped at 0).
    """
    w = max(3, int(window))
    base = grey_dilation(grey_erosion(trace, size=w), size=w)
    return np.clip(trace - base, 0, None)


def _exp_decay(x, a, k, c):
    return a * np.exp(-k * x) + c


def correct_decay(trace: np.ndarray, floor_percentile: float = 10.0) -> np.ndarray:
    """Flatten the exponential loss of signal across the read (RT processivity).

    Fits an exponential to the lower envelope of the trace and divides it out, so
    early and late nucleotides are on a comparable amplitude scale. This is the
    modern analogue of QuShape's ``decayCorrectionExp``.
    """
    y = trace.astype(float)
    x = np.arange(len(y))
    # lower envelope: rolling low percentile
    win = max(25, len(y) // 40)
    env = np.array([np.percentile(y[max(0, i - win):i + win + 1], floor_percentile)
                    for i in range(0, len(y), win)])
    ex = np.arange(0, len(y), win)
    try:
        p0 = (y.max(), 1.0 / len(y), np.median(y))
        popt, _ = curve_fit(_exp_decay, ex, env, p0=p0, maxfev=5000)
        trend = _exp_decay(x, *popt)
        trend = np.clip(trend, 1e-6, None)
        trend /= trend.max()
        return y / trend
    except Exception:
        return y   # if the fit fails, leave the data untouched (documented behavior)


def mobility_shift(reference: np.ndarray, moving: np.ndarray,
                   max_shift: int = 40) -> np.ndarray:
    """Correct a constant offset between two co-loaded dye channels.

    Different fluorophores migrate slightly differently; when a sequencing ladder
    and the reagent trace are co-loaded in one capillary, their peaks are offset.
    This finds the integer shift maximizing cross-correlation and applies it.
    (For separate-capillary designs use :mod:`shapece.align` instead.)
    """
    ref = reference - reference.mean()
    mov = moving - moving.mean()
    best_r, best_k = -np.inf, 0
    for k in range(-max_shift, max_shift + 1):
        r = np.corrcoef(np.roll(mov, k), ref)[0, 1]
        if r > best_r:
            best_r, best_k = r, k
    return np.roll(moving, best_k)
