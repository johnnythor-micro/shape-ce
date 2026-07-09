"""
shapece.stats
=============
Between-condition comparison of normalized reactivity profiles, for
capillary-electrophoresis (CE) probing experiments with **replicates**.

Three complementary layers, from cheapest to most rigorous:

1. :func:`probe_signal_qc` -- a go/no-go gate. Did the probe actually write signal
   into the RNA? If not, nothing downstream is meaningful. **Run this first.**
2. :func:`delta_shape` -- a faithful port of deltaSHAPE (Smola et al. 2015):
   window-smoothed differences, Z-factors, standard scores, and windowed
   site-calling. Answers *"which regions changed?"*
3. :func:`ttest_fdr` / :func:`permutation_test` -- classical per-nucleotide testing
   with multiple-testing correction, and a global condition-separation test.
   Answers *"which individual nucleotides changed?"* and *"did anything change at
   all?"* respectively.

`delta_shape` and `ttest_fdr` are **complementary, not redundant**: the former
gains power by requiring several changed nucleotides to cluster in a window (a
structural rearrangement is local), the latter treats nucleotides independently
and controls the false-discovery rate across them. Neither rescues an experiment
that fails :func:`probe_signal_qc`.

Adaptation note for CE
----------------------
deltaSHAPE was designed for SHAPE-MaP, which yields a per-nucleotide standard
error natively. CE does not. Here the per-nucleotide error is the **standard error
of the mean across replicates** (:func:`sem_from_replicates`), which is the
appropriate and defensible substitute. At least two replicates per condition are
required; three or more are strongly recommended.

References
----------
Smola MJ, Calabrese JM, Weeks KM. Detection of RNA-protein interactions in living
cells with SHAPE. *Biochemistry* 54:6867-6875 (2015).  [deltaSHAPE]
"""
from __future__ import annotations

import warnings
from contextlib import contextmanager
from itertools import combinations, groupby

import numpy as np
from scipy import stats as _st


@contextmanager
def _quiet():
    """Silence all-NaN-slice warnings: no-data nucleotides are expected in SHAPE."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            from scipy.stats import SmallSampleWarning
            warnings.simplefilter("ignore", SmallSampleWarning)
        except ImportError:
            pass
        with np.errstate(invalid="ignore", divide="ignore"):
            yield

__all__ = [
    "sem_from_replicates", "replicate_stats", "probe_signal_qc",
    "smooth_profile", "z_factor", "standard_score", "delta_shape",
    "ttest_fdr", "permutation_test", "benjamini_hochberg",
]


# ---------------------------------------------------------------------------
# Replicate handling
# ---------------------------------------------------------------------------
def sem_from_replicates(matrix) -> tuple[np.ndarray, np.ndarray]:
    """Mean and standard error of the mean across replicates.

    Parameters
    ----------
    matrix : array-like, shape (n_replicates, n_nucleotides)
        Normalized reactivity, one row per replicate.

    Returns
    -------
    (mean, sem) : tuple of np.ndarray

    Raises
    ------
    ValueError
        If fewer than two replicates are supplied (no error can be estimated).
    """
    m = np.asarray(matrix, float)
    if m.ndim != 2 or m.shape[0] < 2:
        raise ValueError(
            "At least 2 replicates are required to estimate a per-nucleotide "
            "standard error for CE data (got shape %r)." % (m.shape,))
    n = m.shape[0]
    with _quiet():
        mean = np.nanmean(m, axis=0)
        sd = np.nanstd(m, axis=0, ddof=1)
    return mean, sd / np.sqrt(n)


def replicate_stats(matrix) -> dict:
    """Summarize a replicate matrix: mean, SD, SEM, and reproducibility.

    ``loo_r`` is the leave-one-out Pearson correlation of each replicate against
    the mean of the others -- the standard reproducibility check that should
    precede any between-condition comparison.
    """
    m = np.asarray(matrix, float)
    mean, sem = sem_from_replicates(m)
    n = m.shape[0]
    loo = []
    with _quiet():
        for i in range(n):
            others = np.nanmean(np.delete(m, i, axis=0), axis=0)
            ok = np.isfinite(m[i]) & np.isfinite(others)
            loo.append(_st.pearsonr(m[i][ok], others[ok])[0] if ok.sum() > 2 else np.nan)
        pair = [
            _st.pearsonr(*_finite_pair(m[i], m[j]))[0]
            for i, j in combinations(range(n), 2)
        ]
        sd_all = np.nanstd(m, axis=0, ddof=1)
    return {"mean": mean, "sd": sd_all, "sem": sem,
            "loo_r": np.array(loo), "pairwise_r": np.array(pair),
            "n_replicates": n}


def _finite_pair(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return a[ok], b[ok]


# ---------------------------------------------------------------------------
# QC gate -- run before any differential analysis
# ---------------------------------------------------------------------------
def probe_signal_qc(area_rx, area_bg, scale: float | None = None) -> dict:
    """Did the chemical probe produce signal above background?

    In a healthy probing experiment the (+) reagent lane departs substantially
    from its (-) background lane at reactive nucleotides. If the two lanes are
    nearly identical, background subtraction cannot isolate reactivity, and what
    survives is the natural-stop pattern rescaled -- which can masquerade as a
    beautifully reproducible profile.

    Returns a dict with the three diagnostics and a boolean ``passed``:

    ``rx_bg_correlation``
        Pearson r between (+) and (-) peak areas. Should be **< ~0.95**.
    ``mean_area_ratio``
        mean(RX)/mean(BG). Should be clearly **> 1**.
    ``background_carryover``
        Pearson r between the background-subtracted reactivity and the background
        areas. Should collapse toward **0**; a large positive value means the
        "reactivity" is still the background pattern.

    Notes
    -----
    Thresholds are deliberately permissive defaults for flagging, not hard
    physical constants; inspect the numbers rather than trusting the boolean.
    """
    rxa = np.asarray(area_rx, float)
    bga = np.asarray(area_bg, float)
    a, b = _finite_pair(rxa, bga)
    r_rxbg = _st.pearsonr(a, b)[0] if len(a) > 2 else np.nan

    if scale is None:
        # same low-reactivity scaling used by shapece.reactivity
        order = np.argsort(a)[: max(3, int(len(a) * 0.25))]
        A, B = a[order], b[order]
        scale = float(np.sum(A * B) / np.sum(B * B)) if np.sum(B * B) > 0 else 1.0
    diff = a - scale * b
    r_carry = _st.pearsonr(diff, b)[0] if len(a) > 2 else np.nan
    ratio = float(np.nanmean(rxa) / np.nanmean(bga)) if np.nanmean(bga) else np.nan

    reasons = []
    if not (r_rxbg < 0.95):
        reasons.append("(+) and (-) lanes are nearly identical (r >= 0.95)")
    if not (ratio > 1.05):
        reasons.append("mean (+)/(-) area ratio is not appreciably > 1")
    if not (abs(r_carry) < 0.5):
        reasons.append("reactivity still tracks the background (|r| >= 0.5)")

    return {"rx_bg_correlation": float(r_rxbg),
            "mean_area_ratio": ratio,
            "background_carryover": float(r_carry),
            "scale_factor": float(scale),
            "passed": len(reasons) == 0,
            "reasons": reasons}


# ---------------------------------------------------------------------------
# deltaSHAPE -- faithful port
# ---------------------------------------------------------------------------
def smooth_profile(data, err, pad: int = 1):
    """Centered-window smoothing of a profile and its errors (deltaSHAPE).

    Data are averaged over ``2*pad+1`` nucleotides, ignoring no-data positions.
    Errors are combined as ``sqrt(sum(err**2)) / window_length``.

    .. note::
       The error normalization divides by the *window length*, not by
       ``sqrt(n_valid)``. This is preserved verbatim from the published
       deltaSHAPE implementation so results remain comparable to it.

    Positions within ``pad`` of either end (where no centered window fits), and
    any position that was no-data in the input, are returned as NaN.
    """
    data = np.asarray(data, float).copy()
    err = np.asarray(err, float).copy()
    data[data == -999] = np.nan
    err[~np.isfinite(data)] = np.nan
    n = len(data)
    if pad == 0:
        return data, err

    mask = ~np.isfinite(data)
    out_d = np.full(n, np.nan)
    out_e = np.full(n, np.nan)
    win_len = 2 * pad + 1
    for i in range(pad, n - pad):
        seg_d = data[i - pad:i + pad + 1]
        seg_e = err[i - pad:i + pad + 1]
        if np.isfinite(seg_d).any():
            with _quiet():
                out_d[i] = np.nanmean(seg_d)
        fin = np.isfinite(seg_e)
        if fin.any():
            out_e[i] = np.sqrt(np.nansum(seg_e[fin] ** 2)) / win_len
    out_d[mask] = np.nan
    out_e[mask] = np.nan
    return out_d, out_e


def z_factor(data1, data2, err1, err2, coeff: float = 1.96) -> np.ndarray:
    """deltaSHAPE Z-factor: ``1 - coeff*(err1+err2)/|data2-data1|``.

    ``Z > 0`` means the two confidence intervals (at ``coeff`` standard errors,
    1.96 -> 95%) do **not** overlap. Undefined (NaN) where the difference is zero
    or either input is missing.
    """
    d1 = np.asarray(data1, float); d2 = np.asarray(data2, float)
    e1 = np.asarray(err1, float); e2 = np.asarray(err2, float)
    bot = np.abs(d2 - d1)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = 1.0 - (coeff * (e1 + e2) / bot)
    z[~np.isfinite(d1) | ~np.isfinite(d2) | (bot == 0)] = np.nan
    return z


def standard_score(diffs) -> np.ndarray:
    """Standard (Z) score of the smoothed differences: ``(x - mean) / sigma``."""
    d = np.asarray(diffs, float)
    with _quiet():
        mu = np.nanmean(d)
        sd = np.nanstd(d)
    if not np.isfinite(sd) or sd == 0:
        return np.full(len(d), np.nan)
    return (d - mu) / sd


def _group_consecutive(indices):
    out = []
    for _, g in groupby(enumerate(sorted(indices)), lambda t: t[0] - t[1]):
        out.append([t[1] for t in g])
    return out


def delta_shape(profile1, profile2, err1=None, err2=None, *,
                replicates1=None, replicates2=None,
                pad: int = 1, z_coeff: float = 1.96, z_thresh: float = 0.0,
                ss_thresh: float = 1.0, site_pad: int = 2, site_min: int = 3,
                mask5: int = 0, mask3: int = 0) -> dict:
    """deltaSHAPE differential analysis between two conditions.

    Supply either mean profiles **plus** per-nucleotide errors, or replicate
    matrices (in which case the mean and the SEM are computed for you -- the
    recommended route for CE data).

    Parameters
    ----------
    profile1, profile2 : array-like or None
        Mean normalized reactivity per nucleotide for each condition.
    err1, err2 : array-like or None
        Per-nucleotide standard errors.
    replicates1, replicates2 : array-like, shape (n_reps, n_nt), optional
        Replicate matrices; if given, override ``profile*``/``err*``.
    pad : int
        Smoothing half-window (window = ``2*pad+1``). Default 1, per deltaSHAPE.
        Set 0 to disable smoothing.
    z_coeff, z_thresh : float
        Z-factor coefficient (1.96 -> 95% CI) and the threshold a nucleotide must
        exceed. Default requires non-overlapping 95% CIs (``Z > 0``).
    ss_thresh : float
        Minimum |standard score| of the smoothed difference.
    site_pad, site_min : int
        A site requires ``site_min`` significant nucleotides within a
        ``2*site_pad+1``-nt window. Default: 3 within 5.
    mask5, mask3 : int
        Nucleotides to ignore at the 5'/3' ends (e.g. primer-binding sites).

    Returns
    -------
    dict with keys
        ``diff``            raw difference (cond1 - cond2)
        ``smoothed_diff``   window-smoothed difference
        ``z_factors``       per-nucleotide Z-factor
        ``z_scores``        standard score of the smoothed difference
        ``significant``     boolean mask of nucleotides passing both filters
                            *and* belonging to a called site
        ``sites``           list of dicts (start, end, direction, mean_diff),
                            1-based inclusive coordinates
        ``n_sites``         number of called sites

    Notes
    -----
    Sign convention follows deltaSHAPE: ``diff = profile1 - profile2``, so a
    positive site is *more reactive in condition 1*.
    """
    if site_min > 2 * site_pad + 1:
        raise ValueError("site_min cannot exceed the window size 2*site_pad+1.")

    if replicates1 is not None and replicates2 is not None:
        profile1, err1 = sem_from_replicates(replicates1)
        profile2, err2 = sem_from_replicates(replicates2)
    if profile1 is None or profile2 is None or err1 is None or err2 is None:
        raise ValueError("Provide replicate matrices, or profiles together with errors.")

    d1 = np.asarray(profile1, float).copy()
    d2 = np.asarray(profile2, float).copy()
    e1 = np.asarray(err1, float).copy()
    e2 = np.asarray(err2, float).copy()
    for arr in (d1, d2):
        arr[arr == -999] = np.nan
    n = len(d1)
    if not (len(d2) == len(e1) == len(e2) == n):
        raise ValueError("All inputs must have the same length.")

    # mask primer / user-excluded regions
    if mask5:
        d1[:mask5] = d2[:mask5] = np.nan
    if mask3:
        d1[n - mask3:] = d2[n - mask3:] = np.nan
    e1[~np.isfinite(d1)] = np.nan
    e2[~np.isfinite(d2)] = np.nan

    s_d1, s_e1 = smooth_profile(d1, e1, pad)
    s_d2, s_e2 = smooth_profile(d2, e2, pad)

    diff = d1 - d2
    s_diff, _ = smooth_profile(diff, e1, pad)

    zf = z_factor(s_d1, s_d2, s_e1, s_e2, z_coeff)
    zs = standard_score(s_diff)

    passes = (zf > z_thresh) & (np.abs(zs) >= ss_thresh)
    passes = np.where(np.isfinite(zf) & np.isfinite(zs), passes, False)

    # windowed site calling: site_min hits inside a (2*site_pad+1) window
    hits: set[int] = set()
    for i in range(site_pad, n - site_pad):
        win = range(i - site_pad, i + site_pad + 1)
        maybes = [j for j in win if passes[j]]
        if len(maybes) >= site_min:
            hits.update(maybes)

    significant = np.zeros(n, bool)
    significant[list(hits)] = True

    sites = []
    for grp in _group_consecutive(hits):
        md = float(np.nanmean(s_diff[grp]))
        sites.append({"start": grp[0] + 1, "end": grp[-1] + 1,
                      "direction": "increase" if md >= 0 else "decrease",
                      "mean_diff": md, "n_nt": len(grp)})
    sites.sort(key=lambda s: s["start"])

    return {"diff": diff, "smoothed_diff": s_diff,
            "z_factors": zf, "z_scores": zs,
            "significant": significant, "sites": sites, "n_sites": len(sites),
            "profile1": d1, "profile2": d2, "err1": e1, "err2": e2,
            "smoothed1": s_d1, "smoothed2": s_d2}


# ---------------------------------------------------------------------------
# Classical per-nucleotide and global tests
# ---------------------------------------------------------------------------
def benjamini_hochberg(pvals) -> np.ndarray:
    """Benjamini-Hochberg FDR q-values. NaN p-values are treated as 1.0."""
    p = np.asarray(pvals, float).copy()
    p[~np.isfinite(p)] = 1.0
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.clip(ranked, 0, 1)
    return q


def ttest_fdr(replicates1, replicates2, min_delta: float = 0.3,
              q_thresh: float = 0.10) -> dict:
    """Per-nucleotide Welch t-test with BH-FDR control and an effect-size gate.

    A nucleotide is called only if ``q < q_thresh`` **and** ``|delta| >= min_delta``.
    The effect-size gate prevents calling statistically detectable but
    biologically trivial differences.

    .. warning::
       With n=3 replicates and hundreds of nucleotides this test is typically
       **underpowered after FDR correction** -- it is common for zero nucleotides
       to survive even when a real global difference exists. Interpret a null
       result as "insufficient power", not "no difference". See also
       :func:`delta_shape`, which gains power by requiring changes to cluster.
    """
    a = np.asarray(replicates1, float)
    b = np.asarray(replicates2, float)
    with _quiet():
        m1, m2 = np.nanmean(a, axis=0), np.nanmean(b, axis=0)
        delta = m1 - m2
        t, p = _st.ttest_ind(a, b, axis=0, equal_var=False, nan_policy="omit")
    p = np.asarray(p, float)
    p[~np.isfinite(p)] = 1.0
    q = benjamini_hochberg(p)
    sig = (q < q_thresh) & (np.abs(delta) >= min_delta)
    return {"delta": delta, "mean1": m1, "mean2": m2,
            "t": np.asarray(t, float), "p": p, "q": q, "significant": sig,
            "n_significant": int(sig.sum())}


def permutation_test(replicates1, replicates2) -> dict:
    """Exact label-permutation test of global condition separation.

    The statistic is (mean within-condition correlation) - (mean between-condition
    correlation). All balanced relabelings are enumerated exactly.

    .. warning::
       With 3 vs 3 replicates there are only 10 distinct partitions, so the
       smallest attainable p-value is 0.10. A non-significant result under this
       design does not imply the conditions are identical.
    """
    a = np.asarray(replicates1, float)
    b = np.asarray(replicates2, float)
    allr = np.vstack([a, b])
    n1, n2 = a.shape[0], b.shape[0]
    idx = list(range(n1 + n2))

    def _stat(g1, g2):
        def mc(g):
            if len(g) < 2:
                return np.nan
            return np.nanmean([_st.pearsonr(*_finite_pair(allr[i], allr[j]))[0]
                               for i, j in combinations(g, 2)])
        within = np.nanmean([mc(g1), mc(g2)])
        between = np.nanmean([_st.pearsonr(*_finite_pair(allr[i], allr[j]))[0]
                              for i in g1 for j in g2])
        return within - between

    with _quiet():
        obs = _stat(list(range(n1)), list(range(n1, n1 + n2)))
    stats_all = []
    seen = set()
    for g1 in combinations(idx, n1):
        g2 = tuple(sorted(set(idx) - set(g1)))
        key = tuple(sorted([g1, g2]))
        if key in seen:
            continue
        seen.add(key)
        with _quiet():
            stats_all.append(_stat(list(g1), list(g2)))
    stats_all = np.array(stats_all)
    p = float(np.mean(stats_all >= obs))
    return {"observed": float(obs), "null": stats_all, "p_value": p,
            "n_partitions": len(stats_all),
            "min_attainable_p": 1.0 / len(stats_all)}
