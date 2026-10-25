"""Offline checks of the SLICE RX J0152 rule functions (PLAN.md rules 5, 7) and image list."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import slice_relics_rxj0152 as sr  # noqa: E402


def test_q1_answer():
    lt, gl, zi = sr.MODELS
    assert sr.q1_answer({lt: "lenstool_like", gl: "near_1", zi: "near_1"}) == "H_sys"
    assert sr.q1_answer({lt: "near_1", gl: "near_1", zi: "near_1"}) == "H_none"
    assert sr.q1_answer({lt: "near_1", gl: "near_1", zi: None}) == "H_none"
    assert sr.q1_answer({lt: "lenstool_like", gl: "intermediate", zi: "near_1"}) == "undecided"
    assert sr.q1_answer({lt: None, gl: "near_1", zi: "near_1"}) == "undecided"
    assert sr.q1_answer({lt: "lenstool_like", gl: None, zi: None}) == "undecided"
    assert sr.q1_answer({lt: "intermediate", gl: "near_1", zi: "near_1"}) == "undecided"


def test_per_model_misses():
    ids = np.array(["7.1", "7.2", "7.3"])
    d = np.array([2.9, 3.1, np.nan])
    assert sr.per_model_misses(ids, d) == ["7.2", "7.3"]


def test_models_registered():
    for m in sr.MODELS:
        assert m in sr.lc.MODELS
        assert sr.lc.MODELS[m]["z_lens"] == 0.8269


def test_image_list():
    t = sr.load()
    assert len(t) == 20
    assert sorted(set(t["system"])) == [str(i) for i in range(1, 8)]
    assert list(t["image_id"][t["held_out"]]) == ["7.1", "7.2"]
    assert list(t["image_id"][t["z_spec"]]) == ["1.1", "1.2", "1.3", "1.4"]
    assert np.all(np.asarray(t["z"])[t["z_spec"]] == 3.93)


def test_unmatched_predictions():
    pred = np.array([[0.0, 0.0], [1.0, 0.0], [10.0, 10.0]])
    obs = np.array([[0.1, 0.0], [0.9, 0.1]])
    assert list(sr.unmatched_predictions(pred, obs)) == [2]
