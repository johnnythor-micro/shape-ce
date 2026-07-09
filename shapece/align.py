"""
shapece.align
=============
Co-register lanes so a given sample index means the same nucleotide everywhere.
Two backends are provided; pick based on your experiment:

* :func:`size_standard_align` -- **recommended when an internal size standard
  (e.g. GeneScan LIZ) was co-loaded in every lane.** Warps each lane onto a
  reference lane using the matched standard peaks. Because the standard runs in
  the same capillary as the data, this co-registers the data channels exactly.
  (This is the method that resolved a real cross-lane calibration artifact in our
  development data, where the instrument called inconsistent numbers of standard
  peaks per lane.)

* :func:`dtw_align` -- QuShape-style **banded dynamic time warping** between two
  traces, for designs with no internal size standard. Aligns on shared peak
  structure directly.
"""
from __future__ import annotations
import numpy as np
from scipy.signal import find_peaks
from scipy.interpolate import PchipInterpolator


def detect_standard_peaks(size_trace: np.ndarray, n_expected: int | None = None,
                          prominence_frac: float = 0.02, min_spacing: int = 25) -> np.ndarray:
    """Detect size-standard peak positions (sample indices), consistently.

    Detecting the standard directly (rather than trusting instrument peak calls,
    which can vary lane to lane) is what makes cross-lane warping reliable.
    """
    p, _ = find_peaks(size_trace, prominence=prominence_frac * size_trace.max(),
                      distance=min_spacing)
    p = np.sort(p)
    if n_expected is not None and len(p) > n_expected:
        # keep the n_expected most prominent, then re-sort by position
        proms = size_trace[p]
        keep = np.sort(np.argsort(proms)[::-1][:n_expected])
        p = p[keep]
    return p.astype(float)


def size_standard_align(traces: dict, size_channels: dict, reference: str,
                        roles=("RX", "BG", "S1", "S2")) -> dict:
    """Warp every lane's data channels onto the reference lane's sample axis.

    Parameters
    ----------
    traces : dict[str, dict[str, np.ndarray]]
        lane_name -> {role -> trace}. Must include the roles you want warped.
    size_channels : dict[str, np.ndarray]
        lane_name -> size-standard trace for that lane.
    reference : str
        Lane whose axis becomes the common ruler.
    roles : tuple
        Which roles to warp.

    Returns
    -------
    dict
        lane_name -> {role -> warped trace} (all on the reference axis).
    """
    ref_peaks = detect_standard_peaks(size_channels[reference])
    n = len(ref_peaks)
    ref_idx = np.arange(len(size_channels[reference]), dtype=float)
    out = {}
    for lane, chans in traces.items():
        lp = detect_standard_peaks(size_channels[lane])
        # match peak counts to the reference (keep the n most prominent, in order)
        if len(lp) != n:
            lp = detect_standard_peaks(size_channels[lane], n_expected=n)
        if len(lp) != n:
            # last resort: truncate/pad by order (documented as approximate)
            m = min(len(lp), n)
            lp, rp = lp[:m], ref_peaks[:m]
        else:
            rp = ref_peaks
        # map reference sample index -> this lane's sample index, then resample
        to_lane = PchipInterpolator(rp, lp, extrapolate=True)
        warped = {}
        for role in roles:
            if role in chans:
                src = chans[role]
                warped[role] = np.interp(to_lane(ref_idx), np.arange(len(src)), src)
        out[lane] = warped
    return out


def _banded_dtw(x: np.ndarray, y: np.ndarray, band_frac: float = 0.1):
    """Banded DTW returning the warping path (i, j) index arrays. O(n*band)."""
    n, m = len(x), len(y)
    r = int(max(band_frac * max(n, m), 5))
    INF = np.inf
    D = np.full((n + 1, m + 1), INF)
    D[0, 0] = 0.0
    for i in range(1, n + 1):
        jlo = max(1, int(i * m / n) - r)
        jhi = min(m, int(i * m / n) + r)
        for j in range(jlo, jhi + 1):
            cost = abs(x[i - 1] - y[j - 1])
            D[i, j] = cost + min(D[i - 1, j], D[i, j - 1], D[i - 1, j - 1])
    # backtrace
    i, j = n, m
    pi, pj = [i - 1], [j - 1]
    while i > 1 or j > 1:
        step = np.argmin([D[i - 1, j - 1], D[i - 1, j], D[i, j - 1]])
        if step == 0:
            i, j = i - 1, j - 1
        elif step == 1:
            i -= 1
        else:
            j -= 1
        pi.append(i - 1); pj.append(j - 1)
    return np.array(pi[::-1]), np.array(pj[::-1])


def dtw_align(traces: dict, reference: str, align_on: str = "S1",
              roles=("RX", "BG", "S1", "S2"), band_frac: float = 0.1) -> dict:
    """Align lanes to a reference by DTW on a shared channel (QuShape-style).

    ``align_on`` is the role used to compute the warp (a sequencing ladder shared
    by all lanes works well); the warp is then applied to all ``roles``.
    """
    ref = traces[reference][align_on]
    ref = (ref - ref.mean()) / (ref.std() + 1e-9)
    ref_idx = np.arange(len(ref), dtype=float)
    out = {}
    for lane, chans in traces.items():
        mov = chans[align_on]
        movn = (mov - mov.mean()) / (mov.std() + 1e-9)
        pi, pj = _banded_dtw(ref, movn, band_frac)
        # build reference-index -> lane-index map from the path
        to_lane = PchipInterpolator(*_dedup(pi, pj), extrapolate=True)
        warped = {}
        for role in roles:
            if role in chans:
                src = chans[role]
                warped[role] = np.interp(to_lane(ref_idx), np.arange(len(src)), src)
        out[lane] = warped
    return out


def _dedup(a, b):
    """Make x strictly increasing for PchipInterpolator (average duplicate x)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    ua, inv = np.unique(a, return_inverse=True)
    ub = np.array([b[inv == k].mean() for k in range(len(ua))])
    return ua, ub
