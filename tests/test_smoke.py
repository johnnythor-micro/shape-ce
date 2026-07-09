"""Smoke tests + a template for validating against QuShape's TPP practice data.

Run with:  pytest -q
"""
import numpy as np
import shapece as sc


def test_is_a_real_package_not_a_namespace_package():
    """Guard against __init__.py being lost during packaging.

    Without __init__.py Python silently treats shapece/ as a *namespace package*:
    `import shapece` still succeeds and `from shapece import stats` works, but no
    submodule is auto-imported and __version__ is absent. This test fails loudly
    in that case.
    """
    assert getattr(sc, "__file__", None) is not None, (
        "shapece imported as a namespace package -- __init__.py is missing "
        "from the installed distribution.")
    assert hasattr(sc, "__version__")


def test_imports():
    for m in ["io", "preprocess", "align", "peaks", "reactivity",
              "sequence", "structure", "report", "stats", "plots"]:
        assert hasattr(sc, m), f"shapece.{m} not exposed by __init__.py"


def test_reactivity_and_norm():
    rx = np.array([10, 100, 12, 90, 11, 8, 200, 9.0] * 8)
    bg = np.array([9, 10, 11, 12, 10, 8, 15, 9.0] * 8)
    r = sc.reactivity.area_difference(rx, bg)
    norm, factor = sc.reactivity.boxplot_normalize(r)
    assert factor > 0 and np.nanmax(norm) >= 1.0


def test_fold_shape_directed():
    seq = "GGGAAACUUCGGUUUCCCAAAGGGAAACUUCGGUUUCCC"
    react = np.array([1.5 if b in "AU" else 0.1 for b in seq])
    ss, mfe = sc.structure.fold(seq, react, method="deigan")
    assert len(ss) == len(seq) and mfe < 0


# --- TEMPLATE: validation against QuShape TPP practice data ---
# Point these at QuShape/TPP_Practice_Data/ and compare exported reactivities to the
# finished .qushape reference. Left as a template so the repo ships without large data.
def test_tpp_validation_template():
    pass
