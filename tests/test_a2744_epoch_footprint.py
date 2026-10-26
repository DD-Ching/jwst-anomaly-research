import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from a2744_epoch_footprint import (  # noqa: E402
    classify,
    gnomonic,
    grid_points,
    inside,
    outcome,
    parse_s_region,
)


def test_parse_s_region_frames_and_multiple_polygons():
    s = "POLYGON ICRS 1 2 3 4 5 6 POLYGON 7 8 9 10 11 12 13 14"
    polys = parse_s_region(s)
    assert [p.shape for p in polys] == [(3, 2), (4, 2)]
    assert polys[1][0].tolist() == [7.0, 8.0]
    with pytest.raises(ValueError):
        parse_s_region("CIRCLE 1 2 3")


def test_gnomonic_scale_near_centre():
    xy = gnomonic(np.array([3.5875]), np.array([-30.3975 + 1 / 60]), 3.5875, -30.3975)
    assert xy[0, 0] == pytest.approx(0.0, abs=1e-9)
    assert xy[0, 1] == pytest.approx(1.0, rel=1e-4)
    # one arcmin of true angle along RA at this Dec
    dra = (1 / 60) / np.cos(np.radians(-30.3975))
    xy = gnomonic(np.array([3.5875 + dra]), np.array([-30.3975]), 3.5875, -30.3975)
    assert xy[0, 0] == pytest.approx(1.0, rel=1e-3)


def test_grid_area_of_a_square():
    pts, cell = grid_points(5.0, 1.0)
    sq = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)
    assert inside(pts, [sq]).sum() * cell == pytest.approx(4.0, rel=0.01)
    disc = np.hypot(pts[:, 0], pts[:, 1]) <= 1.0
    assert disc.sum() * cell == pytest.approx(np.pi, rel=0.01)


def test_classify_and_outcome():
    half = 0.5 * np.pi
    assert classify({"F200W": {"a_ov": 9.0, "a_core": half + 0.01}}) == "core"
    assert classify({"F200W": {"a_ov": 1.0, "a_core": 0.1}}) == "overlap"
    assert classify({"F200W": {"a_ov": 0.9, "a_core": 0.0}}) == "none"
    assert classify({}) == "none"
    one = {"dt_days": 10, "per_filter": {"F200W": {"a_ov": 9, "a_core": 2.0}}, "class": "core"}
    two = {
        "dt_days": 5,
        "per_filter": {f: {"a_ov": 9, "a_core": 2.0} for f in ("F200W", "F277W")},
        "class": "core",
    }
    three = {**two, "dt_days": 50}
    assert outcome([one]) == ("O", None)
    assert outcome([one, two, three])[1] is three
    assert outcome([{"per_filter": {}, "class": "none", "dt_days": 1}]) == ("N", None)
