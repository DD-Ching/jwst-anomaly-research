"""Offline checks of the SPT0615 image list and rule wrappers (PLAN.md rules 5 and 7)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import spt0615_images as sp  # noqa: E402


def test_q1_answer():
    lt, gl = sp.MODELS
    assert sp.q1_answer({lt: "lenstool_like", gl: "near_1"}) == "H_sys"
    assert sp.q1_answer({lt: "near_1", gl: "near_1"}) == "H_none"
    assert sp.q1_answer({lt: "intermediate", gl: "near_1"}) == "undecided"
    assert sp.q1_answer({lt: None, gl: "near_1"}) == "undecided"


def test_models_registered():
    for m in sp.MODELS:
        assert sp.lc.MODELS[m]["z_lens"] == 0.972


def test_amplitude_median_skips_3b():
    fams = {"1": {"scale": 0.6}, "3": {"scale": 0.8}, "3b": {"scale": 5.0}}
    assert abs(sp.amplitude_median(fams) - 0.7) < 1e-12
    assert sp.amplitude_median({"3b": {"scale": 1.0}}) is None


def test_image_list():
    t = sp.load()
    assert len(t) == 14
    assert sorted(set(t["system"])) == ["1", "2", "3", "3b"]
    assert int(np.sum(t["held_out"])) == 6
    assert list(t["image_id"][t["z_spec"]]) == ["1.1", "1.2", "1.3", "3.1", "3.2", "3.5"]
    assert set(sp.READ_A) <= set(t["image_id"][t["held_out"]])
    # 3.1 is 06 15 49.96, -57 46 53.5
    r = t[t["image_id"] == "3.1"][0]
    assert abs(r["ra"] - 93.958167) < 1e-5 and abs(r["dec"] + 57.781528) < 1e-5


def test_build_matches_committed_list():
    b, t = sp.build_images(), sp.load()
    assert list(b["image_id"]) == list(t["image_id"])
    assert np.allclose(b["ra"], t["ra"]) and np.allclose(b["dec"], t["dec"])
    assert list(b["z"]) == list(t["z"]) and list(b["use_spec"]) == list(t["z_spec"])
