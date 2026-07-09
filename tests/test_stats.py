"""Tests for `shapece.stats`.

Includes a fidelity check of the deltaSHAPE port against a direct transcription
of the published Python-2 implementation (Smola et al. 2015), on identical inputs.
"""
import numpy as np
import pytest

from shapece import stats as S


# ---------------------------------------------------------------------------
# Verbatim transcriptions of the published deltaSHAPE routines, for fidelity.
# ---------------------------------------------------------------------------
def _ref_smooth(data, err, pad):
    new_data, new_err = [], []
    mask = [i for i in range(len(data)) if (data[i] == -999 or np.isnan(data[i]))]
    for _ in range(pad):
        new_data.append(np.nan); new_err.append(np.nan)
    for i in range(pad, len(data) - pad):
        seg = [j for j in data[i - pad:i + pad + 1]]
        new_data.append(np.mean(np.ma.MaskedArray(seg, np.isnan(seg))))
        errs = np.array(err[i - pad:i + pad + 1])
        squerrs = np.power([j for j in errs if not np.isnan(j)], 2)
        new_err.append(np.sqrt(np.sum(squerrs)) / len(data[i - pad:i + pad + 1]))
    for _ in range(pad):
        new_data.append(np.nan); new_err.append(np.nan)
    for i in mask:
        new_data[i] = np.nan; new_err[i] = np.nan
    return np.array(new_data, float), np.array(new_err, float)


def _ref_zfactor(d1, d2, e1, e2, factor=1.96):
    out = []
    for i in range(len(d1)):
        if np.isnan(d1[i]) or np.isnan(d2[i]):
            out.append(np.nan); continue
        top = factor * (e2[i] + e1[i])
        bot = abs(d2[i] - d1[i])
        out.append(np.nan if bot == 0 else 1 - (top / bot))
    return np.array(out)


def _fixture():
    rng = np.random.RandomState(7)
    n = 60
    d1 = np.abs(rng.normal(0.4, 0.3, n))
    d2 = d1.copy()
    d2[20:25] += 1.2                      # planted 5-nt change
    e1 = np.abs(rng.normal(0.05, 0.01, n))
    e2 = np.abs(rng.normal(0.05, 0.01, n))
    d1[5] = np.nan; d2[5] = np.nan        # a no-data nucleotide
    e1[5] = np.nan; e2[5] = np.nan
    return d1, d2, e1, e2


def test_smoothing_matches_published_implementation():
    d1, d2, e1, e2 = _fixture()
    mine_d, mine_e = S.smooth_profile(d1, e1, pad=1)
    ref_d, ref_e = _ref_smooth(d1, e1, 1)
    assert np.allclose(mine_d, ref_d, equal_nan=True, atol=1e-12)
    assert np.allclose(mine_e, ref_e, equal_nan=True, atol=1e-12)


def test_z_factor_matches_published_implementation():
    d1, d2, e1, e2 = _fixture()
    s1, se1 = S.smooth_profile(d1, e1, 1)
    s2, se2 = S.smooth_profile(d2, e2, 1)
    assert np.allclose(S.z_factor(s1, s2, se1, se2),
                       _ref_zfactor(s1, s2, se1, se2), equal_nan=True, atol=1e-12)


def test_z_factor_sign_semantics():
    """Z > 0 exactly when the 1.96-SE intervals do not overlap."""
    # difference 1.0, errors 0.1 each -> 1.96*0.2 = 0.392 < 1.0 -> Z > 0
    assert S.z_factor([0.0], [1.0], [0.1], [0.1])[0] > 0
    # difference 0.1, errors 0.1 each -> 1.96*0.2 = 0.392 > 0.1 -> Z < 0
    assert S.z_factor([0.0], [0.1], [0.1], [0.1])[0] < 0


def test_delta_shape_recovers_planted_site():
    d1, d2, e1, e2 = _fixture()
    res = S.delta_shape(d1, d2, e1, e2, pad=1)
    assert res["n_sites"] == 1
    site = res["sites"][0]
    # planted at 0-based 20..24 -> 1-based 21..25; smoothing widens by +/-1
    assert site["start"] <= 21 and site["end"] >= 25
    assert site["direction"] == "decrease"      # cond2 higher -> diff negative


def test_delta_shape_no_false_positive_on_identical_profiles():
    rng = np.random.RandomState(3)
    d = np.abs(rng.normal(0.4, 0.3, 60))
    e = np.full(60, 0.05)
    res = S.delta_shape(d, d.copy(), e, e.copy(), pad=1)
    assert res["n_sites"] == 0


def test_delta_shape_from_replicates_uses_sem():
    rng = np.random.RandomState(11)
    n = 60
    base = np.abs(rng.normal(0.4, 0.2, n))
    r1 = np.array([base + rng.normal(0, 0.02, n) for _ in range(3)])
    hot = base.copy(); hot[30:36] += 1.0
    r2 = np.array([hot + rng.normal(0, 0.02, n) for _ in range(3)])
    res = S.delta_shape(None, None, replicates1=r1, replicates2=r2)
    assert res["n_sites"] >= 1
    assert any(s["direction"] == "decrease" for s in res["sites"])


def test_masking_excludes_ends():
    d1, d2, e1, e2 = _fixture()
    res = S.delta_shape(d1, d2, e1, e2, mask5=30, mask3=30)
    assert res["n_sites"] == 0          # the planted site (21-25) is masked out


def test_sem_requires_two_replicates():
    with pytest.raises(ValueError):
        S.sem_from_replicates(np.zeros((1, 10)))


def test_site_min_cannot_exceed_window():
    with pytest.raises(ValueError):
        S.delta_shape(np.zeros(10), np.zeros(10), np.ones(10), np.ones(10),
                      site_pad=1, site_min=5)


def test_benjamini_hochberg_monotone_and_bounded():
    p = np.array([0.001, 0.01, 0.02, 0.5, 0.9])
    q = S.benjamini_hochberg(p)
    assert np.all(q >= p - 1e-12) and np.all(q <= 1.0)
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)   # non-decreasing


def test_ttest_fdr_effect_size_gate():
    rng = np.random.RandomState(5)
    n = 50
    a = np.array([rng.normal(1.0, 0.001, n) for _ in range(3)])
    b = np.array([rng.normal(1.05, 0.001, n) for _ in range(3)])   # tiny but sig.
    res = S.ttest_fdr(a, b, min_delta=0.3)
    assert res["n_significant"] == 0        # gate rejects trivial effect sizes
    res2 = S.ttest_fdr(a, b, min_delta=0.01)
    assert res2["n_significant"] > 0


def test_permutation_test_reports_power_limit():
    rng = np.random.RandomState(2)
    a = np.array([rng.normal(0, 1, 40) for _ in range(3)])
    b = np.array([rng.normal(0, 1, 40) for _ in range(3)])
    res = S.permutation_test(a, b)
    assert res["n_partitions"] == 10
    assert np.isclose(res["min_attainable_p"], 0.1)


def test_probe_signal_qc_flags_degenerate_experiment():
    rng = np.random.RandomState(9)
    bg = np.abs(rng.normal(100, 30, 100))
    rx = bg * 1.02 + rng.normal(0, 0.5, 100)      # (+) ~= (-)
    qc = S.probe_signal_qc(rx, bg)
    assert qc["passed"] is False
    assert qc["rx_bg_correlation"] > 0.95

    good_bg = np.abs(rng.normal(100, 30, 100))
    good_rx = np.abs(rng.normal(400, 200, 100))
    assert S.probe_signal_qc(good_rx, good_bg)["passed"] is True


def test_replicate_stats_shapes():
    rng = np.random.RandomState(4)
    m = rng.normal(0.5, 0.1, (3, 40))
    rs = S.replicate_stats(m)
    assert rs["mean"].shape == (40,) and rs["sem"].shape == (40,)
    assert rs["loo_r"].shape == (3,) and rs["pairwise_r"].shape == (3,)
