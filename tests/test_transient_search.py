"""Tests for scripts/transient_search.py (synthetic two-epoch catalogs)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("transient_search", _DIR / "transient_search.py")
ts = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ts)

RA0, DEC0 = 110.0, -73.0


def _catalog(n, rng, mag, err):
    cosd = np.cos(np.deg2rad(DEC0))
    return Table(
        {
            "label": np.arange(1, n + 1),
            "ra": RA0 + rng.uniform(-30, 30, n) / 3600 / cosd,
            "dec": DEC0 + rng.uniform(-30, 30, n) / 3600,
            "aper50_abmag": mag,
            "aper50_abmag_err": err,
        }
    )


def test_finds_variable_appeared_and_disappeared_sources():
    rng = np.random.default_rng(5)
    n = 300
    mag = rng.uniform(22, 27, n)
    mag[0] = 23.0  # the variable source is bright enough to test (S/N ~80)
    err = 2.5 * np.log10(1 + 1 / (200 * 10 ** (-0.4 * (mag - 22))))  # S/N 200 at mag 22
    cat1 = _catalog(n, rng, mag, err)
    cat2 = cat1.copy()
    cat2["aper50_abmag"] = mag + 0.05 + rng.normal(0, 0.005, n)  # zero point offset 0.05
    cat2["aper50_abmag"][0] = mag[0] + 0.05 - 1.0  # source 0 brightened by 1 mag
    cat2.remove_row(1)  # source 1 disappeared
    cosd = np.cos(np.deg2rad(DEC0))
    new = {
        "label": 9999,
        "ra": RA0 + 1.0 / 3600 / cosd,
        "dec": DEC0,
        "aper50_abmag": 22.5,
        "aper50_abmag_err": 0.01,
    }
    cat2.add_row(new)  # a new source inside the footprint
    cat1["aper50_abmag"][1] = 22.2
    cat1["aper50_abmag_err"][1] = 0.005
    res = ts.search(cat1, cat2, tie_radius_arcsec=100.0)
    kinds = {k: list(res["kind"]).count(k) for k in ("variable", "appeared", "disappeared")}
    assert kinds == {"variable": 1, "appeared": 1, "disappeared": 1}, list(res)
    var = res[res["kind"] == "variable"][0]
    assert abs(var["dmag"] + 1.0) < 0.05 and var["label1"] == 1  # zero point removed
    assert res[res["kind"] == "appeared"][0]["label2"] == 9999
    assert res[res["kind"] == "disappeared"][0]["label1"] == 2
    assert res.meta["thresholds"]["provenance"] == "assumption"


def test_quiet_epochs_give_no_candidates():
    rng = np.random.default_rng(6)
    n = 200
    mag = rng.uniform(22, 26, n)
    err = np.full(n, 0.02)
    cat1 = _catalog(n, rng, mag, err)
    cat2 = cat1.copy()
    cat2["aper50_abmag"] = mag + rng.normal(0, 0.02, n)
    assert len(ts.search(cat1, cat2, tie_radius_arcsec=100.0)) == 0


_spec2 = importlib.util.spec_from_file_location("transient_forced", _DIR / "transient_forced.py")
tf = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(tf)


def test_forced_aperture_flux_and_comparison():
    yy, xx = np.mgrid[0:41, 0:41]
    star = 100.0 * np.exp(-((xx - 20) ** 2 + (yy - 20) ** 2) / (2 * 1.5**2)) + 5.0  # pedestal 5
    f, e = tf.aperture_flux(star, None, 20, 20, 6, 12, 18)
    assert f > 0 and np.isfinite(e)
    assert abs(f - 100.0 * 2 * np.pi * 1.5**2) / f < 0.02  # background removed
    hole = star.copy()
    hole[20, 20] = np.nan
    assert np.isnan(tf.aperture_flux(hole, None, 20, 20, 6, 12, 18)[0])
    f1 = np.array([10.0, 10.0, 10.0, 10.0])
    f2 = np.array([10.0, 10.0, 10.0, 2.5])  # the last one faded by 1.5 mag
    err = np.full(4, 1e-3)
    dm, sig, sig_flux = tf.compare(f1, err, f2, err, sys_floor=0.05, min_zp_refs=3)
    assert abs(dm[-1] - 2.5 * np.log10(4)) < 1e-6 and abs(sig[-1]) > 25 and abs(sig[0]) < 1
    # a source that vanished (flux ~0 in epoch 2) has no magnitude but a clear flux change
    dm, sig, sig_flux = tf.compare(
        np.array([10.0]), np.array([0.1]), np.array([-0.05]), np.array([0.1])
    )
    assert np.isnan(dm[0]) and sig_flux[0] < -15  # floor-limited: 0.05 mag of 10 units
    # too few positions: no zero point is fitted, so one real variable is not absorbed
    dm, _, _ = tf.compare(np.array([10.0]), np.array([0.01]), np.array([5.0]), np.array([0.01]))
    assert abs(dm[0] - 2.5 * np.log10(2)) < 1e-6


def test_frame_offset_between_epochs_is_tied_before_matching():
    rng = np.random.default_rng(7)
    n = 300
    mag = rng.uniform(22, 25, n)
    cat1 = _catalog(n, rng, mag, np.full(n, 0.01))
    cat2 = cat1.copy()
    cosd = np.cos(np.deg2rad(DEC0))
    cat2["ra"] = cat1["ra"] + 0.25 / 3600 / cosd  # 0.25" offset, close to the match radius
    res = ts.search(cat1, cat2, tie_radius_arcsec=100.0)
    assert len(res) == 0  # no false appeared/disappeared pairs
    assert abs(res.meta["frame_shift_arcsec"][0] - 0.25) < 0.01
    assert res.meta["n_without_local_tie"] == 0


_spec3 = importlib.util.spec_from_file_location("transient_combine", _DIR / "transient_combine.py")
tc = importlib.util.module_from_spec(_spec3)
_spec3.loader.exec_module(tc)


def _cands(rows):
    return Table(rows=rows, names=("kind", "ra", "dec"), dtype=(str, float, float))


def test_coincident_keeps_same_kind_in_two_bands_once():
    d = 0.1 / 3600  # 0.1" apart
    a = _cands(
        [("variable", RA0, DEC0), ("appeared", RA0 + 0.01, DEC0), ("variable", RA0 + 0.02, DEC0)]
    )
    b = _cands([("variable", RA0 + d, DEC0), ("disappeared", RA0 + 0.01, DEC0)])
    c = _cands([("variable", RA0, DEC0 + d)])
    out = tc.coincident({"F090W": a, "F115W": b, "F277W": c})
    assert len(out) == 1  # different kinds at RA0+0.01 and the single-band one are dropped
    assert out[0]["n_bands"] == 3 and out[0]["bands"] == "F090W,F115W,F277W"
    assert len(tc.coincident({"F090W": a, "F115W": b}, min_bands=3)) == 0


def test_bright_gaia_exclusion_radius_grows_with_brightness():
    r = tc.exclusion_radius(np.array([22.0, 20.0, 17.5, 5.0, np.nan]))
    assert (
        np.allclose(r[:2], 1.5)
        and abs(r[2] - 1.5 * 10**0.5) < 1e-9
        and r[3] == 12.0
        and r[4] == 1.5
    )
    gaia = Table({"ra": [RA0], "dec": [DEC0], "gmag": [15.0]})  # 15" radius clipped to 12"
    cosd = np.cos(np.deg2rad(DEC0))
    ra = RA0 + np.array([5.0, 11.0, 13.0]) / 3600 / cosd
    assert list(tc.near_bright(ra, np.full(3, DEC0), gaia)) == [0, 0, -1]


def test_recentre_finds_an_offset_source():
    yy, xx = np.mgrid[0:41, 0:41]
    img = 50.0 * np.exp(-((xx - 23.0) ** 2 + (yy - 18.0) ** 2) / (2 * 1.5**2))
    cx, cy = tf.recentre(img, 20.0, 20.0, 5)
    assert abs(cx - 23.0) < 0.2 and abs(cy - 18.0) < 0.2
    assert tf.recentre(np.zeros((41, 41)), 20.0, 20.0, 5) == (20.0, 20.0)  # empty box: unchanged
    x, y = tf.recentre(img, np.nan, 20.0, 5)
    assert np.isnan(x) and y == 20.0


def test_coincident_merges_chains_into_one_group():
    cosd = np.cos(np.deg2rad(DEC0))
    step = 0.2 / 3600 / cosd  # A~B and B~C within 0.3", A and C 0.4" apart
    tables = {
        "F090W": _cands([("variable", RA0, DEC0)]),
        "F115W": _cands([("variable", RA0 + step, DEC0)]),
        "F277W": _cands([("variable", RA0 + 2 * step, DEC0)]),
    }
    out = tc.coincident(tables)
    assert len(out) == 1 and out[0]["n_bands"] == 3
    empty = tc.coincident({"F090W": _cands([]), "F115W": _cands([])})
    assert len(empty) == 0 and "n_bands" in empty.colnames and empty.meta["provenance"] == "derived"


def test_recentre_iterates_beyond_its_box_and_rejects_off_image_starts():
    yy, xx = np.mgrid[0:41, 0:41]
    img = 50.0 * np.exp(-((xx - 24.5) ** 2 + (yy - 20.0) ** 2) / (2 * 1.5**2)) + 3.0
    cx, cy = tf.recentre(img, 20.0, 20.0, 3)  # 4.5 px away, box half-width 3
    assert abs(cx - 24.5) < 0.2 and abs(cy - 20.0) < 0.2
    assert tf.recentre(img, -10.0, 20.0, 2) == (-10.0, 20.0)


def test_coincident_without_tables_is_an_error():
    import pytest

    with pytest.raises(ValueError):
        tc.coincident({})
