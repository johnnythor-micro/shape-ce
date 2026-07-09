"""
Validation of `shapece` against QuShape's published TPP-riboswitch reference.

QuShape ships a *finished* analysis (`test1001_done.qushape`) alongside the raw
`.fsa` files. That project file is a Python-2 `shelve` (Berkeley DB) containing,
among other things, the published peak areas (`dPeakRX/BG`), the background-
subtracted reactivity (`areaDiff`), and the normalized reactivity (`normDiff`).

These tests reproduce that reference with the modern pipeline at four levels of
strictness. Requires QuShape's practice data; skipped if it is not present.

    git clone https://github.com/Weeks-UNC/QuShape.git
    pytest tests/test_validation_tpp.py -q

Reading the reference file needs a Berkeley DB binding and GUI-class stubs
(the pickles embed PyQt4/sip objects):

    apt-get install -y libdb-dev && pip install bsddb3
"""
import io as _io
import os
import pickle

import numpy as np
import pytest
from scipy.stats import pearsonr

import shapece as sc
from shapece.io import read_abif, get_channel
from shapece import preprocess as pp, align, peaks, reactivity as rx

QUSHAPE_DIR = os.environ.get("QUSHAPE_DIR", "QuShape/TPP_Practice_Data")
REF_PROJECT = os.path.join(QUSHAPE_DIR, "test1001_done.qushape")
RX_FSA = os.path.join(QUSHAPE_DIR, "TPP_+1M7.fsa")
BG_FSA = os.path.join(QUSHAPE_DIR, "TPP_DMSO.fsa")

# QuShape's region of interest, recovered by cross-correlating our preprocessed
# raw trace against the project's stored preprocessed trace.
ROI_OFFSET = 1329
RX_DATA_CHANNEL = 2   # per TPP README: column 2 = reagent
S1_DATA_CHANNEL = 3   # per TPP README: column 3 = ddC sequencing ladder


def _load_reference():
    """Unpickle QuShape's finished project, stubbing out its PyQt4/sip classes."""
    bsddb3 = pytest.importorskip("bsddb3", reason="pip install bsddb3 (needs libdb-dev)")

    class _Dummy:
        def __init__(self, *a, **k): pass
        def __setstate__(self, s): self.state = s

    class _SafeUnpickler(pickle.Unpickler):
        def find_class(self, module, name):
            if module.startswith(("sip", "PyQt4", "PyQt5")):
                return type("D", (_Dummy,), {})
            try:
                return super().find_class(module, name)
            except Exception:
                return type("D", (_Dummy,), {})

    db = bsddb3.db
    d = db.DB(); d.open(REF_PROJECT, None, db.DB_HASH, db.DB_RDONLY)
    proj = _SafeUnpickler(_io.BytesIO(d[b"dProject"]), encoding="latin1").load()
    d.close()
    return proj


needs_data = pytest.mark.skipif(
    not os.path.exists(REF_PROJECT), reason="QuShape TPP practice data not present")


@needs_data
def test_area_difference_is_exact_subtraction():
    """QuShape's areaDiff equals RX area - BG area exactly (BG pre-scaled)."""
    P = _load_reference()
    rx_a = np.asarray(P["dPeakRX"]["area"], float)
    bg_a = np.asarray(P["dPeakBG"]["area"], float)
    ad = np.asarray(P["areaDiff"], float)
    coef, *_ = np.linalg.lstsq(np.vstack([rx_a, bg_a]).T, ad, rcond=None)
    assert np.allclose(coef, [1.0, -1.0], atol=1e-6)


@needs_data
def test_normalization_matches_reference_up_to_a_constant():
    """shapece's boxplot normalization reproduces normDiff exactly in shape.

    Pearson r is 1.0 to 10 decimals. The normalization *divisor* differs by
    ~6.5% because the shipped reference project was produced by an earlier
    QuShape build whose outlier count on these data was 1 rather than the 2 that
    the current published source computes (`findPOutlierBox`). We therefore assert
    exact proportionality, not bit-identity.
    """
    P = _load_reference()
    ad = np.asarray(P["areaDiff"], float)
    nd = np.asarray(P["normDiff"], float)
    mine, factor = rx.boxplot_normalize(ad)
    assert pearsonr(mine, nd)[0] > 1 - 1e-9
    ratio = (ad / nd)[np.isfinite(ad / nd)]
    assert np.allclose(ratio, np.median(ratio), rtol=1e-6)   # single constant


@needs_data
def test_peak_areas_agree():
    """Gaussian deconvolution reproduces QuShape's peak areas (r > 0.99)."""
    P = _load_reference()
    for lane, dkey in [("dPeakRX", "RX"), ("dPeakBG", "BG")]:
        pos = np.asarray(P[lane]["pos"], int)
        ref_area = np.asarray(P[lane]["area"], float)
        trace = np.asarray(P["dData"][dkey], float)
        mine = peaks.quantify(trace, pos, mode="gaussian")["area"]
        assert pearsonr(mine, ref_area)[0] > 0.99


@needs_data
def test_end_to_end_from_raw_fsa():
    """Full pipeline from raw .fsa reproduces the published profile (r > 0.85)."""
    P = _load_reference()
    nd = np.asarray(P["normDiff"], float)
    ref_pos = np.asarray(P["dPeakRX"]["pos"], int)
    L = len(np.asarray(P["dData"]["RX"], float))

    def prep(path):
        e = read_abif(path)
        f = lambda c: pp.baseline(pp.smooth(pp.correct_saturation(get_channel(e, c))))
        return f(RX_DATA_CHANNEL), f(S1_DATA_CHANNEL)

    rxT, rxS = prep(RX_FSA)
    bgT, bgS = prep(BG_FSA)
    sl = slice(ROI_OFFSET, ROI_OFFSET + L)
    traces = {"RX": {"RX": rxT[sl], "S1": rxS[sl]},
              "BG": {"RX": bgT[sl], "S1": bgS[sl]}}
    # no internal size standard in this dataset -> DTW on the shared ddC ladder
    al = align.dtw_align(traces, reference="RX", align_on="S1",
                         roles=("RX", "S1"), band_frac=0.06)
    assert pearsonr(al["RX"]["S1"], al["BG"]["S1"])[0] > 0.95   # alignment worked

    q_rx = peaks.quantify(al["RX"]["RX"], ref_pos, mode="gaussian")
    q_bg = peaks.quantify(al["BG"]["RX"], ref_pos, mode="gaussian")
    norm, _ = rx.boxplot_normalize(
        rx.area_difference(q_rx["area"], q_bg["area"], scale=True))
    assert pearsonr(norm, nd)[0] > 0.85
