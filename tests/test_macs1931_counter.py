"""Offline checks of the MACS J1931 footprint step (polygon parsing and containment)."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import macs1931_counter as mc  # noqa: E402

SQUARE = "POLYGON ICRS 10.0 -5.0 10.1 -5.0 10.1 -4.9 10.0 -4.9"
TWO = SQUARE + " POLYGON 20.0 1.0 20.2 1.0 20.2 1.2 20.0 1.2"


def test_parse_polygons_one_and_two():
    (p,) = mc.parse_polygons(SQUARE)
    assert p.shape == (4, 2)
    assert p[1].tolist() == [10.1, -5.0]
    assert len(mc.parse_polygons(TWO)) == 2


@pytest.mark.parametrize("bad", ["", "CIRCLE 1 2 3", "POLYGON 1 2 3"])
def test_parse_polygons_refuses_malformed(bad):
    with pytest.raises(ValueError):
        mc.parse_polygons(bad)


def test_inside_any_polygon():
    polys = mc.parse_polygons(TWO)
    hit = mc.inside([10.05, 20.1, 15.0, 10.2], [-4.95, 1.1, 0.0, -4.95], polys)
    assert hit.tolist() == [True, True, False, False]


def test_footprint_counts_used_and_spec_rows():
    images = Table(
        {
            "ra": [10.05, 10.06, 10.5, 10.07],
            "dec": [-4.95, -4.94, -4.95, -4.93],
            "used": [True, True, True, False],
            "z_spec": [True, False, True, True],
        }
    )
    obs = Table(
        {
            "obs_id": ["a"],
            "filters": ["F200W"],
            "target_name": ["T"],
            "t_min": [60000.0],
            "s_region": [SQUARE],
        }
    )
    (row,) = mc.footprint(images, obs)
    assert (row["n_used_inside"], row["n_used"]) == (2, 3)
    assert (row["n_spec_inside"], row["n_spec"]) == (1, 2)


def test_committed_summary_covers_every_used_image():
    path = Path(__file__).resolve().parents[1] / mc.OUT / "summary.json"
    s = __import__("json").loads(path.read_text())
    assert s["mosaics"], "no mosaic recorded"
    assert all(r["n_used_inside"] == r["n_used"] == 51 for r in s["mosaics"])
    assert np.isclose(s["f200w_i2d_header"]["pixel_arcsec"], 0.0312, atol=2e-4)
