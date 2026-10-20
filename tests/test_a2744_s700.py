"""Offline tests for scripts/a2744_s700.py (pure helpers; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_s700 as s7  # noqa: E402

RA0, DEC0 = 3.5797, -30.3577
ARC = 1 / 3600.0


def _t(offsets):
    cd = np.cos(np.radians(DEC0))
    return Table(
        {
            "ra": [RA0 + dx * ARC / cd for dx, _ in offsets],
            "dec": [DEC0 + dy * ARC for _, dy in offsets],
            "magnification": [-3.0] * len(offsets),
        }
    )


def test_assign_is_one_to_one():
    cat = _t([(0, 0), (1.7, 0), (15, 5)])
    # one predicted image between a and b: nearest is shared, the assignment is not
    pred = _t([(0.85, 0), (15.2, 5), (40, 0)])
    idx, sep = s7.assign(pred, cat)
    assert len(set(idx)) == 3
    assert sep[0] < 1.0 and sep[2] < 1.0 and sep[1] > 20
    near = s7.nearest(pred, cat)
    assert np.all(near < 1.0)


def test_assign_too_few_images():
    idx, sep = s7.assign(_t([(0, 0)]), _t([(0, 0), (2, 0)]))
    assert idx == [] and np.all(np.isnan(sep))


def test_verdict_outcomes():
    rows = [
        {"z": 1.0, "rms_arcsec": 3.0, "reproduced": False},
        {"z": 1.1, "rms_arcsec": None, "reproduced": False},
    ]
    v = s7.verdict(rows)
    assert v["outcome"].startswith("B") and v["best_assignment"]["z"] == 1.0
    rows.append({"z": 1.2, "rms_arcsec": 0.4, "reproduced": True})
    v = s7.verdict(rows)
    assert v["outcome"].startswith("A") and v["z_reproduced"] == [1.2, 1.2]


def test_gated_match_picks_smallest_passing_radius():
    k = {1.0: {"fraction": 0.7}, 1.5: {"fraction": 0.85}}
    assert s7.gated_match(k) == 1.5
    k[1.0]["fraction"] = 0.8
    assert s7.gated_match(k) == 1.0
    assert s7.gated_match({1.0: {"fraction": 0.5}, 1.5: {"fraction": None}}) is None


def test_counterparts_nearest_with_redshift_metadata(tmp_path):
    cd = np.cos(np.radians(DEC0))
    raw = Table(
        {
            "ID": [78.21, 5.1],
            "RA": [RA0 + 0.02 * ARC / cd, RA0 + 30 * ARC / cd],
            "Dec": [DEC0, DEC0],
            "z": [2.3151, 1.5],
            "Redshift_type": ["spec", "phot"],
            "Redshift_reference": ["Price+24", "x"],
            "Redshift_instrument": ["JWST/NIRSpec", "x"],
        }
    )
    path = tmp_path / "images.fits"
    raw.write(path)
    cat = _t([(0, 0)])
    cat["image_id"] = ["700.1a"]
    (row,) = s7.counterparts(cat, path)
    assert row["uncover_id"] == "78.21" and row["redshift_type"] == "spec"
    assert row["sep_arcsec"] < 0.05


def test_failure_kind():
    assert s7.failure_kind(3, [0, 1], [0.1, 0.1, 0.1], 1.5, n_pred=2) == "missing"
    assert s7.failure_kind(3, [0, 0, 1], [0.5, 1.2, 0.3], 1.5, n_pred=3) == "shared"
    assert s7.failure_kind(3, [0, 1, 2], [0.5, 1.7, 0.3], 1.5, n_pred=3) == "offset"
    assert s7.failure_kind(2, [0, 0], [0.5, 1.7], 1.5, n_pred=3) == "offset"


def test_nearest_index_shares():
    cat = _t([(0, 0), (1.7, 0)])
    idx, sep = s7.nearest_index(_t([(0.85, 0), (20, 0)]), cat)
    assert idx == [0, 0] and np.all(sep < 1.0)
    idx, sep = s7.nearest_index(_t([]), cat)
    assert idx == [-1, -1] and np.all(np.isnan(sep))


def test_counterparts_beyond_radius_is_none(tmp_path):
    cd = np.cos(np.radians(DEC0))
    raw = Table(
        {
            "ID": [5.1],
            "RA": [RA0 + 30 * ARC / cd],
            "Dec": [DEC0],
            "z": [1.5],
            "Redshift_type": ["phot"],
            "Redshift_reference": ["x"],
            "Redshift_instrument": ["x"],
        }
    )
    raw.write(tmp_path / "i.fits")
    cat = _t([(0, 0)])
    cat["image_id"] = ["700.1a"]
    (row,) = s7.counterparts(cat, tmp_path / "i.fits")
    assert row["uncover_id"] is None and row["sep_arcsec"] > 29
    assert s7.counterpart_sources([row]) == frozenset()
    assert s7.counterpart_sources([{"uncover_id": "78.21"}]) == frozenset({"78"})


def test_pair_geometry_midpoints():
    t = _t([(0, 0), (2, 0), (0.5, 0), (1.5, 0)])
    t["system"] = ["s1", "s1", "s2", "s2"]
    t["image_id"] = ["s1a", "s1b", "s2a", "s2b"]
    g = s7.pair_geometry(t)
    assert [round(x["sep_arcsec"], 3) for x in g] == [2.0, 1.0]
    assert abs(g[0]["mid_dx_arcsec"] - g[1]["mid_dx_arcsec"]) < 1e-3
