import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import g09725_sw as g  # noqa: E402

PIX = 0.03


def test_measure_recovers_point_source_and_noise():
    rng = np.random.default_rng(1)
    shape = (201, 201)
    diff = rng.normal(0, 1.0, shape)
    yy, xx = np.indices(shape, dtype=float)
    psf = np.exp(-0.5 * ((xx - 100) ** 2 + (yy - 100) ** 2) / 1.5**2)
    diff += 500 * psf / psf.sum()
    err = np.full(shape, 1.0)
    r = g.measure(diff, err, (100.0, 100.0), PIX)
    assert abs(r["flux"] - 500) < 5 * r["sigma"]
    assert r["flux"] / r["sigma"] > 5
    # aperture of radius 0.2"/0.03" ~ 6.7 px: sqrt(pi r^2) ~ 11.8 for unit noise
    assert 8 < r["sigma_prop"] < 15
    assert r["sigma"] >= r["sigma_prop"]


def test_classify_rules():
    assert g.classify(6.0, 1.0, True) == "brighter"
    assert g.classify(-5.0, 1.0, True) == "fainter"
    assert g.classify(4.9, 1.0, True) == "none"
    assert g.classify(10.0, 1.0, False) == "unregistered"
    assert g.classify(np.nan, 1.0, True) == "unregistered"


def test_outcome_rules():
    lw = {"f277w": "brighter", "f356w": "brighter", "f444w": "none"}
    assert g.outcome({**lw, "f115w": "none", "f150w": "brighter", "f200w": "none"}) == "SW_detected"
    assert g.outcome({**lw, "f115w": "none", "f150w": "none", "f200w": "none"}) == "SW_none"
    assert g.outcome({**lw, "f115w": "fainter", "f150w": "none", "f200w": "none"}) == "SW_fainter"
    assert g.outcome({**lw, "f115w": "unregistered", "f150w": "none", "f200w": "none"}) == "U"
    bad = {"f277w": "brighter", "f356w": "none", "f444w": "none"}
    assert g.outcome({**bad, "f115w": "brighter", "f150w": "none", "f200w": "none"}) == "U"


def test_uri():
    u = g.uri(g.LATE, "f150w")
    assert u.endswith("jw02561/L3/t/o006/jw02561-o006_t007_nircam_clear-f150w_i2d.fits")
