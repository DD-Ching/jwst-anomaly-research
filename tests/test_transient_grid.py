"""Tests for scripts/transient_grid.py (synthetic two-epoch mosaics)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("transient_grid", _DIR / "transient_grid.py")
tg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tg)

SCALE = 0.03  # arcsec per pixel
N = 400


def _wcs() -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [24.35, -8.45]
    w.wcs.crpix = [N / 2, N / 2]
    w.wcs.cdelt = [-SCALE / 3600, SCALE / 3600]
    return w


def _image(path: Path, xy: np.ndarray, flux: np.ndarray, rng: np.random.Generator) -> None:
    sci = rng.normal(0, 0.01, (N, N))
    for (x, y), f in zip(xy, flux, strict=True):
        sci += tg.gaussian((N, N), x, y, 1.5, f)
    hdus = [fits.PrimaryHDU()]
    hdus.append(fits.ImageHDU(sci.astype("f4"), header=_wcs().to_header(), name="SCI"))
    hdus.append(fits.ImageHDU(np.full((N, N), 0.01, "f4"), header=_wcs().to_header(), name="ERR"))
    fits.HDUList(hdus).writeto(path)


def _field(rng: np.random.Generator, n: int = 150) -> np.ndarray:
    axis = np.arange(30, N - 30, 26)
    g = np.stack(np.meshgrid(axis, axis), -1).reshape(-1, 2)
    return g[:n] + rng.uniform(-2, 2, (min(n, len(g)), 2))


def test_union_positions_merges_duplicates():
    a = Table({"ra": [24.0, 24.001], "dec": [-8.0, -8.0]})
    b = Table({"ra": [24.0 + 0.05 / 3600, 24.002], "dec": [-8.0, -8.0]})
    u = tg.union_positions([a, b], 0.2)
    assert len(u) == 3
    assert list(u["origin"]) == [0, 0, 1]


def test_hits_needs_same_sign_in_min_bands():
    s = np.array([[1, 1, 0, 0], [1, -1, 0, 0], [0, -1, -1, 0], [1, 0, 0, 0]])
    assert list(tg.hits(s, 2)) == [True, False, True, False]


def test_band_pass_with_fixed_zero_point_keeps_a_common_brightening():
    rng = np.random.default_rng(1)
    f1 = rng.uniform(5, 10, 200)
    m = {"f1": f1, "e1": np.full(200, 0.05), "f2": 2 * f1, "e2": np.full(200, 0.05)}
    # its own zero point absorbs a brightening shared by every position
    assert (tg.band_pass(m, 0.05, 0.3, 5.0)[2] == 0).all()
    # the real grid's zero point (0) and scale keep it
    assert (tg.band_pass(m, 0.05, 0.3, 5.0, scale=1.0, zp=0.0)[2] == 1).all()


def test_band_pass_ignores_a_difference_of_two_negative_fluxes():
    m = {"f1": np.full(200, -1.0), "e1": np.full(200, 0.01), "f2": np.full(200, -3.0)}
    m["e2"] = np.full(200, 0.01)
    assert (tg.band_pass(m, 0.05, 0.3, 5.0, scale=1.0, zp=0.0)[2] == 0).all()
    m["f2"] = np.full(200, 2.0)  # appeared: no dmag, but epoch 2 has flux
    assert (tg.band_pass(m, 0.05, 0.3, 5.0, scale=1.0, zp=0.0)[2] == 1).all()


def test_group_hits_links_neighbours():
    ra = np.array([24.0, 24.0 + 1 / 3600, 24.0 + 2 / 3600, 24.1, 24.2])
    dec = np.full(5, -8.0)
    g = tg.group_hits(ra, dec, np.array([True, True, True, True, False]), 1.5)
    assert list(g) == [0, 0, 0, 1, -1]


def test_main_finds_the_brightened_source_and_recovers_injections(tmp_path):
    rng = np.random.default_rng(0)
    xy = _field(rng)
    flux = rng.uniform(2, 6, len(xy))
    w = _wcs()
    ra, dec = w.pixel_to_world_values(xy[:, 0], xy[:, 1])
    flux2 = flux.copy()
    flux2[20] *= 3  # one source brightens in both bands
    args = []
    for band in ("A", "B"):
        p1, p2 = tmp_path / f"{band}1.fits", tmp_path / f"{band}2.fits"
        _image(p1, xy, flux, rng)
        _image(p2, xy, flux2, rng)
        args += ["--band", band, str(p1), str(p2)]
    Table({"ra": ra, "dec": dec}).write(tmp_path / "cat.ecsv", format="ascii.ecsv")
    out = tmp_path / "out"
    rc = tg.main(
        args
        + ["--catalog", str(tmp_path / "cat.ecsv"), "--epoch2-catalogs", "--out", str(out)]
        + ["--inject", "40", "--inject-frac", "1.0", "--procs", "1"]
    )
    assert rc == 0
    grid = Table.read(out / "grid.ecsv")
    assert list(np.flatnonzero(grid["hit"])) == [20]
    summary = json.loads((out / "summary.json").read_text())
    assert summary["n_hits"] == 1
    assert summary["n_hit_groups"] == 1
    assert summary["injection"]["runs"][0]["recovered_frac"] >= 0.9


def test_spike_zone_radius_scales_with_star_brightness():
    ra0, dec0 = 24.35, -8.45
    d = np.array([10.0, 30.0, 50.0]) / 3600  # arcsec north of the star
    ra, dec = np.full(3, ra0), dec0 + d
    bright = tg.spike_zone(ra, dec, np.array([ra0]), np.array([dec0]), np.array([17.3]), 45.0, 17.3)
    assert bright.tolist() == [True, True, False]
    faint = tg.spike_zone(ra, dec, np.array([ra0]), np.array([dec0]), np.array([19.3]), 45.0, 17.3)
    assert faint.tolist() == [True, False, False]  # R = 45 * 10**-0.4 = 17.9 arcsec


def test_bright_point_sources_skips_extended_and_faint():
    from astropy.coordinates import SkyCoord

    cat = Table(
        {
            "sky_centroid": SkyCoord([24.30, 24.31, 24.32], [-8.4, -8.4, -8.4], unit="deg"),
            "aper_total_abmag": [17.0, 17.0, 21.0],
            "is_extended": [False, True, False],
        }
    )
    ra, _, m = tg.bright_point_sources(cat, 19.5)
    assert np.allclose(ra, [24.30]) and m.tolist() == [17.0]
