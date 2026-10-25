"""Offline checks of the SLICE MACS J0358 rule wrapper (PLAN.md rule 7) and image list."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import slice_relics_macs0358 as sm  # noqa: E402


def test_q1_answer():
    lt, gl = sm.MODELS
    assert sm.q1_answer({lt: "lenstool_like", gl: "near_1"}) == "H_sys"
    assert sm.q1_answer({lt: "near_1", gl: "near_1"}) == "H_none"
    assert sm.q1_answer({lt: "near_1", gl: None}) == "undecided"
    assert sm.q1_answer({lt: None, gl: "near_1"}) == "undecided"
    assert sm.q1_answer({lt: "intermediate", gl: "near_1"}) == "undecided"
    assert sm.q1_answer({lt: "lenstool_like", gl: "intermediate"}) == "undecided"


def test_models_registered():
    for m in sm.MODELS:
        assert m in sm.lc.MODELS
        assert sm.lc.MODELS[m]["z_lens"] == 0.425


def test_image_list():
    t = sm.load()
    assert len(t) == 16
    assert sorted(set(t["system"])) == [str(i) for i in range(1, 7)]
    assert sorted(set(t["system"][t["held_out"]])) == ["4", "5", "6"]
    assert int(np.sum(t["held_out"])) == 9
    assert list(t["image_id"][t["z_spec"]]) == ["1.1", "1.2", "1.3", "2.1", "2.2", "3.1", "3.2"]
    z = np.asarray(t["z"])[t["z_spec"]]
    assert list(z) == [1.261] * 3 + [3.07] * 4
