"""Unit tests for scripts/hit_lightcurve.py and hit_lightcurve_controls.py (offline, synthetic)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table


def _load(name: str):
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


hl = _load("hit_lightcurve")
hc = _load("hit_lightcurve_controls")


def _gauss(n: int, x: float, y: float, sigma: float, flux: float) -> np.ndarray:
    yy, xx = np.mgrid[:n, :n]
    return flux * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sigma**2)) / (2 * np.pi * sigma**2)


def test_aperture_recovers_point_flux_over_constant_background():
    img = 3.0 + _gauss(81, 40.0, 40.0, 1.5, 1000.0)
    f = hl.aperture_photometry(img, 40.0, 40.0, 6.0, (8.0, 12.0))
    assert abs(f - 1000.0) / 1000.0 < 0.01


def test_aperture_off_image_or_blank_is_nan():
    img = np.zeros((40, 40))
    assert np.isnan(hl.aperture_photometry(img, 2.0, 20.0, 3.0, (5.0, 8.0)))
    assert np.isnan(hl.aperture_photometry(img, 20.0, 20.0, 3.0, (5.0, 8.0)))


def test_empirical_error_matches_white_noise():
    rng = np.random.default_rng(0)
    img = rng.normal(0.0, 1.0, (300, 300))
    r = 3.0
    err = hl.empirical_error(img, r, (6.0, 10.0), np.random.default_rng(1))
    n_ap = np.sum(np.hypot(*np.mgrid[-3:4, -3:4]) <= r)
    assert 0.8 * np.sqrt(n_ap) < err < 1.3 * np.sqrt(n_ap)


def test_groups_split_by_distance():
    g = hl.groups([("a", 3.55, -30.37), ("b", 3.5502, -30.3702), ("c", 3.60, -30.37)])
    assert [[p[0] for p in grp] for grp in g] == [["a", "b"], ["c"]]


def test_compare_calibrates_against_controls():
    rng = np.random.default_rng(2)
    rows = []
    for i in range(40):
        f1 = 50.0
        f2 = f1 * 0.9 + rng.normal(0, 1.0)  # a common 10 % program offset
        rows += [(f"c{i:02d}", "e1", f1, 1.0), (f"c{i:02d}", "e2", f2, 1.0)]
    rows += [("h1", "e1", 50.0, 1.0), ("h1", "e2", 45.0, 1.0)]  # same as controls
    rows += [("h2", "e1", 50.0, 1.0), ("h2", "e2", 80.0, 1.0)]  # a real brightening
    lc = Table(rows=rows, names=("name", "obs_id", "flux_njy", "err_njy"))
    t = hc.compare(lc, "e1", "e2")
    z = dict(zip(t["name"], t["z_cal"], strict=True))
    assert abs(z["h1"]) < 2
    assert z["h2"] > 10
    tail = dict(zip(t["name"], t["ctrl_tail_frac"], strict=True))
    assert tail["h2"] == 0.0


def test_aperture_local_cutout_matches_full_grid():
    rng = np.random.default_rng(3)
    img = rng.normal(0.0, 1.0, (120, 150)) + _gauss(150, 70.3, 55.7, 2.0, 500.0)[:120]
    r, ann = 4.0, (6.0, 9.5)
    yy, xx = np.ogrid[:120, :150]
    d2 = (xx - 70.3) ** 2 + (yy - 55.7) ** 2
    ap = d2 <= r**2
    bg = img[(d2 >= ann[0] ** 2) & (d2 <= ann[1] ** 2)]
    expect = img[ap].sum() - np.median(bg) * ap.sum()
    assert np.isclose(hl.aperture_photometry(img, 70.3, 55.7, r, ann), expect)


def test_compare_controls_by_name_pattern_and_per_pair_stats():
    rng = np.random.default_rng(4)
    rows = []
    for i in range(30):
        rows += [(f"c{i:02d}", "e1", 50.0, 1.0), (f"c{i:02d}", "e2", 50.0 + rng.normal(), 1.0)]
    rows += [("cand1", "e1", 50.0, 1.0), ("cand1", "e2", 90.0, 1.0)]  # a hit named 'c...'
    lc = Table(rows=rows, names=("name", "obs_id", "flux_njy", "err_njy"))
    t = hc.compare(lc, "e1", "e2")
    assert list(t["name"]) == ["cand1"]
    assert t["n_controls"][0] == 30 and t["ctrl_rsd"][0] > 0
