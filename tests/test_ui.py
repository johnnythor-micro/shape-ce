"""Tests for `shapece.ui` (widget logic; no browser needed)."""
import pytest

import shapece as sc
from shapece.ui import clean_sequence

pytest.importorskip("ipywidgets", reason="pip install ipywidgets")

from shapece.ui import ChannelPicker, SamplePairer, SequenceInput  # noqa: E402


# ---------------------------------------------------------------- sequence
def test_clean_sequence_strips_fasta_whitespace_digits():
    assert clean_sequence(">my rna\nGGGA UCG\nTCAC 123") == "GGGAUCGUCAC"


def test_clean_sequence_uppercases_and_converts_t_to_u():
    assert clean_sequence("acgt") == "ACGU"


def test_clean_sequence_can_keep_dna():
    assert clean_sequence("acgt", to_rna=False) == "ACGT"


def test_clean_sequence_rejects_invalid_characters():
    with pytest.raises(ValueError):
        clean_sequence("ACGXZ")


def test_sequence_input_reports_length_and_gc():
    box = SequenceInput()
    box._area.value = ">t\nAUGC gg"          # noqa: SLF001
    assert box.value == "AUGCGG"
    assert "6 nt" in box._status.value       # noqa: SLF001


def test_sequence_input_flags_bad_characters():
    box = SequenceInput()
    box._area.value = "ACGXZ"                # noqa: SLF001
    assert box.value == ""
    assert "invalid" in box._status.value    # noqa: SLF001


# ---------------------------------------------------------------- pairing
def _pairer():
    paths = {"plus_r1.fsa": "/tmp/a", "minus_r1.fsa": "/tmp/b",
             "plus_r2.fsa": "/tmp/c", "minus_r2.fsa": "/tmp/d",
             "ddATP.fsa": "/tmp/e"}
    sp = SamplePairer(paths, conditions=("A", "B"))
    for name, (role, cond, rep) in sp.rows.items():
        if name.startswith("dd"):
            role.value = "ladder"
        else:
            role.value = "plus" if name.startswith("plus") else "minus"
            rep.value = int(name.split("_r")[1][0])
        cond.value = "A"
    return sp


def test_sample_pairer_matches_plus_and_minus_by_condition_and_replicate():
    cfg = _pairer().to_config()
    assert sorted(cfg["pairs"]) == [("plus_r1.fsa", "minus_r1.fsa"),
                                    ("plus_r2.fsa", "minus_r2.fsa")]
    assert list(cfg["ladders"]) == ["ddATP.fsa"]
    assert cfg["unmatched_plus"] == []
    assert "ddATP.fsa" not in cfg["samples"]      # ladders are not samples


def test_sample_pairer_reports_unmatched_plus_lanes():
    sp = SamplePairer({"p1.fsa": "/tmp/x", "p2.fsa": "/tmp/y"})
    for i, (_, (role, _c, rep)) in enumerate(sp.rows.items()):
        role.value = "plus"
        rep.value = i + 1
    assert len(sp.to_config()["unmatched_plus"]) == 2


def test_sample_pairer_ignore_role_excludes_file():
    sp = SamplePairer({"a.fsa": "/tmp/a", "b.fsa": "/tmp/b"})
    for role, _c, _r in sp.rows.values():
        role.value = "ignore"
    cfg = sp.to_config()
    assert cfg["samples"] == {} and cfg["pairs"] == []


# ---------------------------------------------------------------- channels
def test_channel_picker_options_come_from_the_file(tmp_path):
    """Options are read from the ABIF file, so users pick what actually exists."""
    import os
    qs = os.environ.get("QUSHAPE_DIR", "QuShape/TPP_Practice_Data")
    fsa = os.path.join(qs, "TPP_+1M7.fsa")
    if not os.path.exists(fsa):
        pytest.skip("QuShape practice data not present")
    cp = ChannelPicker(fsa)
    cfg = cp.to_config()
    assert cfg["data_channel"] == 2          # 4-dye file: sensible default
    assert cfg["size_channel"] is None       # no LIZ standard -> DTW alignment
    assert any("DATA1" in label for label, _ in cp.data.options)


def test_ui_is_exposed_on_the_package():
    assert hasattr(sc, "ui")
