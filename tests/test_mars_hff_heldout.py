import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import mars_hff_heldout as mo  # noqa: E402

DEC0 = -24.0


def _pos(dx_arcsec, dy_arcsec, ra0=64.0):
    dx, dy = np.asarray(dx_arcsec, float), np.asarray(dy_arcsec, float)
    return ra0 + dx / 3600.0 / np.cos(np.deg2rad(DEC0)), DEC0 + dy / 3600.0


def _mars(dx, dy):
    ra, dec = _pos(dx, dy)
    return Table({"ra": ra, "dec": dec})


def _canucs(ids, dx, dy, z):
    ra, dec = _pos(dx, dy)
    return Table({"image_id": ids, "ra": ra, "dec": dec, "z": np.asarray(z, float)})


def test_frame_offset_median_and_too_few_pairs():
    base = np.arange(6) * 10.0
    mars = _mars(base + 0.2, np.zeros(6) - 0.1)
    can = _canucs([f"{i}.1" for i in range(6)], base, np.zeros(6), [1.0] * 6)
    off = mo.frame_offset(can, mars)
    assert off["applied"] and off["n_pairs"] == 6
    assert np.isclose(off["dx_arcsec"], 0.2, atol=1e-6) and np.isclose(
        off["dy_arcsec"], -0.1, atol=1e-6
    )
    assert not mo.frame_offset(can, mars[:0])["applied"]
    off4 = mo.frame_offset(can[:4], mars)
    assert not off4["applied"] and off4["dx_arcsec"] == 0.0
    sh = mo.shift(can, off)
    _, _, d = mo._offsets_arcsec(sh["ra"], sh["dec"], mars["ra"], mars["dec"])
    assert np.allclose(d.min(axis=1), 0.0, atol=1e-6)


def test_split_seen_radius_and_z():
    mars = _mars([0.0, 10.0], [0.0, 0.0])
    can = _canucs(
        ["1.1", "1.2", "2.1", "3.1"], [0.4, 10.6, 20.0, 0.0], [0, 0, 0, 0], [1.0, 1.0, 2.0, 0.0]
    )
    t = mo.split_seen(can, mars)
    assert list(t["image_id"]) == ["1.1", "1.2", "2.1"]  # z = 0 dropped
    assert list(t["seen"]) == [True, False, False]
    assert list(t["system"]) == ["1", "1", "2"]


def test_anchors_types():
    t = Table(
        {
            "system": ["1", "1", "1", "2", "2", "3", "3", "4", "5"],
            "seen": [True, True, False, False, False, True, False, False, True],
        }
    )
    a = mo.anchors(t)
    assert set(a) == {"1", "2", "3"}  # 4: one held-out image alone; 5: no held-out image
    assert a["1"][0] == "E" and list(np.where(a["1"][1])[0]) == [0, 1]
    assert a["2"][0] == "N" and list(np.where(a["2"][1])[0]) == [3, 4]
    assert a["3"][0] == "M" and list(np.where(a["3"][1])[0]) == [5, 6]


def _r(d, inside=None):
    d = np.asarray(d, float)
    inside = np.ones(len(d), bool) if inside is None else np.asarray(inside, bool)
    return Table({"dtheta_arcsec": d, "inside": inside, "type": ["E"] * len(d)})


def test_gate_h_thresholds_and_no_prediction():
    assert mo.gate_h(_r([0.5] * 5))["passes"]
    g = mo.gate_h(_r([0.5, 0.5, 0.5, 0.5, np.inf]))  # 80 % within 3", rms of the finite ones
    assert g["passes"] and np.isclose(g["matched_fraction"], 0.8)
    assert not mo.gate_h(_r([0.5, 0.5, 0.5, np.inf, np.inf]))["passes"]
    assert not mo.gate_h(_r([1.2] * 5))["passes"]
    g = mo.gate_h(_r([0.5] * 5 + [np.nan], inside=[True] * 5 + [False]))
    assert g["n_tested"] == 5 and g["passes"]
    assert mo.gate_h(_r([]))["n_tested"] == 0


def test_outcome():
    assert mo.outcome({"n_tested": 4, "passes": True}) == "not testable"
    assert mo.outcome({"n_tested": 5, "passes": True}) == "reproduces"
    assert mo.outcome({"n_tested": 9, "passes": False}) == "does not reproduce"
