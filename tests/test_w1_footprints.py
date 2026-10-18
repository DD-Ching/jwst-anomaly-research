import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_footprints as wf  # noqa: E402


def _square(ra, dec, half_arcmin=1.0):
    d = half_arcmin / 60
    return (
        f"POLYGON ICRS {ra - d} {dec - d} {ra + d} {dec - d} {ra + d} {dec + d} {ra - d} {dec + d}"
    )


def test_polygons_skip_frame_words():
    (p,) = wf.polygons("POLYGON J2000 10.0 20.0 10.1 20.0 10.1 20.1")
    assert p.shape == (3, 2) and p[0].tolist() == [10.0, 20.0]


def test_footprints_overlap_and_ra_wrap():
    # two 2' squares on the equator straddling RA 0, offset by 1': overlap 1 x 2 arcmin²
    rows = [
        {"obs_id": "a_nircam", "s_region": _square(359.99, 0.0), "t_exptime": 1.0},
        {"obs_id": "b_nircam", "s_region": _square(359.99 + 1 / 60 - 360, 0.0), "t_exptime": 1.0},
    ]
    res = wf.footprints(rows)
    assert res["area_arcmin2"]["a"] == pytest.approx(4.0, abs=0.1)
    assert res["overlap_arcmin2"]["a|b"] == pytest.approx(2.0, abs=0.1)
    assert res["centroid"]["a"][0] == pytest.approx(359.99, abs=1e-3)
    assert res["dropped"] == ["a"]  # tie on the summed excess: name order


def test_footprint_outside_grid_raises():
    rows = [
        {"obs_id": "a_nircam", "s_region": _square(10.0, 0.0), "t_exptime": 1.0},
        {"obs_id": "b_nircam", "s_region": _square(12.0, 0.0), "t_exptime": 1.0},
    ]
    with pytest.raises(ValueError, match="grid"):
        wf.footprints(rows)


def test_drop_overlapping_removes_the_worst_pointing_only():
    ov = {"x|y": 6.2, "y|z": 2.7, "x|z": 0.3, "z|w": 0.58}
    assert wf.drop_overlapping(ov) == ["y"]
    assert wf.drop_overlapping({"x|y": 0.6}) == []


def test_tile_grid_assigns_disjoint_tiles_to_the_covering_pointing():
    # a 6' x 2' strip and a 2' square inside its east end: every tile goes to the strip (ties go
    # to name order, "a" < "b"), tiles do not overlap and cover the strip
    d = 1 / 60
    strip = f"POLYGON ICRS {10 - 3 * d} {-d} {10 + 3 * d} {-d} {10 + 3 * d} {d} {10 - 3 * d} {d}"
    rows = [
        {"obs_id": "a_nircam", "s_region": strip, "t_exptime": 1.0},
        {"obs_id": "b_nircam", "s_region": _square(10 + 2 * d, 0.0), "t_exptime": 1.0},
    ]
    tiles = wf.tile_grid(rows, ["a", "b"], half_arcsec=60.0, min_frac=0.3)
    assert {t["pointing"] for t in tiles} == {"a"}
    assert len(tiles) == 3
    ras = sorted(t["centre"][0] for t in tiles)
    assert np.diff(ras) == pytest.approx([2 * d, 2 * d], abs=0.03 * d)


def test_in_box_is_square_and_wraps_ra():
    import inject_radial as ir

    ra = np.array([359.99, 0.0 + 100 / 3600, 359.99 - 150 / 3600, 359.99])
    dec = np.array([0.0, 0.0, 0.0, 141 / 3600])
    assert ir.in_box(ra, dec, (359.99, 0.0), 140.0).tolist() == [True, True, False, False]
