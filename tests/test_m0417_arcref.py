"""Frozen rules of scripts/m0417_arcref.py (counterpart, units, reference, class, gate C2)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np

_SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_SCRIPTS))
_spec = importlib.util.spec_from_file_location("m0417_arcref", _SCRIPTS / "m0417_arcref.py")
ma = importlib.util.module_from_spec(_spec)
sys.modules["m0417_arcref"] = ma
_spec.loader.exec_module(ma)


def _cat(**kw):
    base = {
        "ra": np.array([10.0, 10.0]),
        "dec": np.array([0.0, 0.0005]),
        "a": np.array([0.5, 0.1]),
        "b": np.array([0.1, 0.1]),
        "pa": np.array([0.0, 0.0]),  # major axis North-South
        "flux": np.array([10.0, 1.0]),
        "snr": np.array([50.0, 3.0]),
    }
    base.update(kw)
    return base


def test_ellipse_radius_axes():
    assert np.isclose(ma.ellipse_radius(0.0, 1.0, 1.0, 0.5, 0.0), 1.0)  # along the major axis
    assert np.isclose(ma.ellipse_radius(1.0, 0.0, 1.0, 0.5, 0.0), 2.0)  # along the minor axis
    assert np.isclose(ma.ellipse_radius(1.0, 0.0, 1.0, 0.5, 90.0), 1.0)  # major axis East-West


def test_counterpart_uses_the_elongated_ellipse():
    cat = _cat()
    # 1.2" North of source 0: inside its 3-sigma major axis (1.5"), outside along the minor axis
    assert ma.counterpart(10.0, 1.2 / 3600, cat) == 0
    assert ma.counterpart(10.0 + 1.2 / 3600, 0.0, cat) == -1
    # a large segment far away still counts when its ellipse contains the position (no distance cut)
    big = _cat(dec=np.array([0.0, 30 / 3600]), a=np.array([0.5, 12.0]), b=np.array([0.5, 1.0]))
    assert ma.counterpart(10.0, 1 / 3600, big) == 0  # both contain it; the smaller radius wins
    assert ma.counterpart(10.0, 3 / 3600, big) == 1  # only the far segment contains it
    assert ma.counterpart(10.0, 50 / 3600, big) == 1


def test_units_merge_shared_counterparts():
    ims = [
        {"id": "8.1", "cp": 3, "mu": 10.0},
        {"id": "8.2", "cp": 3, "mu": -7.0},
        {"id": "8.3", "cp": -1, "mu": 4.0},
    ]
    u = ma.units(ims)
    merged = next(x for x in u if x["cp"] == 3)
    assert merged["ids"] == ["8.1", "8.2"] and sum(merged["mus"]) == 17.0
    assert any(x["cp"] == -1 and x["ids"] == ["8.3"] for x in u)


def test_reference_rules():
    cat = _cat()
    assert ma.reference([{"cp": 1, "mus": [5.0]}], cat)[1] == "no_reference"  # S/N 3
    assert ma.reference([{"cp": 0, "mus": [60.0]}], cat)[1] == "no_reference"  # |mu| > 50
    lref, state, k = ma.reference([{"cp": 0, "mus": [5.0]}], cat)
    assert state == "" and lref == 2.0 and k == 0
    neg = _cat(flux=np.array([-10.0, 1.0]), snr=np.array([50.0, 3.0]))
    assert ma.reference([{"cp": 0, "mus": [5.0]}], neg)[1] == "no_reference"  # flux <= 0
    cat2 = _cat(snr=np.array([50.0, 50.0]))
    assert ma.reference([{"cp": 0, "mus": [1.0]}, {"cp": 1, "mus": [1.0]}], cat2)[1] == (
        "inconsistent_reference"
    )
    lref, state, k = ma.reference([{"cp": 0, "mus": [4.0]}, {"cp": 1, "mus": [0.5]}], cat2)
    assert state == "" and lref == 2.0 and k == 1  # smallest sum |mu| wins (L 2.5 vs 2.0 within 3x)


def test_classify():
    depth = 28.0
    f10 = 3631.0 * 10 ** (-0.4 * 26.0)  # m 26: brighter than the 10-sigma limit 27.25
    assert ma.classify(f10, depth, [0.5]) == "recovered"
    assert ma.classify(f10, depth, [5.0]) == "confused"
    assert ma.classify(f10, depth, [0.1]) == "absent"
    assert ma.classify(f10, depth, []) == "absent"
    f7 = 3631.0 * 10 ** (-0.4 * 27.6)  # between 10 and 5 sigma
    assert ma.classify(f7, depth, []) == "ambiguous"
    assert ma.classify(3631.0 * 10 ** (-0.4 * 28.5), depth, []) == "undetectable"


def test_gate_c2():
    ok = ["recovered"] * 8 + ["absent"]
    assert ma.gate_c2(ok, [True] * 3 + [False] * 6)["passes"]
    assert not ma.gate_c2(ok, [True] * 2 + [False] * 7)["passes"]  # too few arc references
    arc_bad = ["absent"] * 2 + ["recovered"] * 8
    assert not ma.gate_c2(arc_bad, [True] * 3 + [False] * 7)["passes"]  # arc fraction 1/3
    assert not ma.gate_c2(["recovered"] * 7, [True] * 7)["passes"]  # fewer than 8 testable
