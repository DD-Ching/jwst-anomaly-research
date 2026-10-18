import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_cal_trail as t  # noqa: E402


def test_classify():
    assert t.classify([]).startswith("undecided")
    assert t.classify([True, True]).startswith("in every exposure")
    assert t.classify([False, False]).startswith("in no exposure")
    assert t.classify([True, False, False]).startswith("in 1 of 3")


def test_contrast_and_track_on_synthetic_line():
    rng = np.random.default_rng(1)
    stamp = rng.normal(0.0, 1.0, (201, 600))
    rows = np.round(100 + 0.01 * (np.arange(600) - 300)).astype(int)
    stamp[rows, np.arange(600)] += 30.0
    prof = t.row_profile(stamp)
    c = t.contrast(prof, 100, scale=0.03)
    assert c["contrast"] > t.IN_SIGMA
    tr = t.line_track(stamp, 100)
    assert abs(tr["slope_px_per_px"] - 0.01) < 0.003
    flat = t.contrast(t.row_profile(rng.normal(0.0, 1.0, (201, 600))), 100, scale=0.03)
    assert flat["contrast"] < t.IN_SIGMA


def test_star_geometry_miss_distance():
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [10.0, 0.0]
    w.wcs.crpix = [1.0, 1.0]
    w.wcs.cdelt = [-0.03 / 3600, 0.03 / 3600]
    sx, sy = (float(v) for v in w.world_to_pixel_values(10.0, 0.0))
    xs = [100.0, 200.0, 300.0]
    stamps = [
        {
            "name": "e",
            "track_post_hoc": {
                "slope_px_per_px": 0.0,
                "rms_resid_px": 0.1,
                "chunk_x_px": xs,
                "chunk_row_px": [4.0, 4.0, 4.0],
            },
        }
    ]
    g = t.star_geometry(w, (10.0, 0.0), 0, 0, stamps)
    assert abs(g["stamp_x_px"] - sx) < 1e-9
    assert abs(g["tracks"]["e"]["miss_distance_px"] - abs(4.0 - sy)) < 1e-9
