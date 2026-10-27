"""Tests for scripts/transient_grid.py (synthetic two-epoch mosaics)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
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


def _con() -> np.ndarray:
    """Context planes: one exposure on the left half, three on the right (bit 33 on plane 1)."""
    con = np.zeros((2, N, N), "i4")
    con[0, :, : N // 2] = 0b1
    con[0, :, N // 2 :] = 0b11
    con[1, :, N // 2 :] = 0b1
    return con


def _image(
    path: Path, xy: np.ndarray, flux: np.ndarray, rng: np.random.Generator, con: bool = False
) -> None:
    sci = rng.normal(0, 0.01, (N, N))
    for (x, y), f in zip(xy, flux, strict=True):
        sci += tg.gaussian((N, N), x, y, 1.5, f)
    hdus = [fits.PrimaryHDU()]
    hdus.append(fits.ImageHDU(sci.astype("f4"), header=_wcs().to_header(), name="SCI"))
    hdus.append(fits.ImageHDU(np.full((N, N), 0.01, "f4"), header=_wcs().to_header(), name="ERR"))
    if con:
        hdus.append(fits.ImageHDU(_con(), name="CON"))
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


def test_hits_any_sign_counts_a_colour_change():
    s = np.array([[1, -1, 0, 0], [1, 1, 0, 0], [1, 0, 0, 0]])
    assert list(tg.hits(s, 2, any_sign=True)) == [True, True, False]
    assert list(tg.hits(s, 2)) == [False, True, False]


def test_catalog_bands_reads_whole_tokens():
    names = [
        "jw01176-o211_t009_nircam_clear-f115w_cat.ecsv",
        "x/jw01176-o212_t009_nircam_clear-f356w_cat.ecsv",
    ]
    assert tg.catalog_bands(names, ["F115W", "F356W"]) == ["F115W", "F356W"]
    # 'a' is inside 'cat' but not a token of it
    assert tg.catalog_bands(["cat_a.ecsv"], ["A", "B"]) == ["A"]
    with pytest.raises(SystemExit):
        tg.catalog_bands(["cat.ecsv"], ["A", "B"])


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
    Table({"ra": ra, "dec": dec}).write(tmp_path / "cat_a.ecsv", format="ascii.ecsv")
    out = tmp_path / "out"
    rc = tg.main(
        args
        + ["--catalog", str(tmp_path / "cat_a.ecsv"), "--epoch2-catalogs", "--out", str(out)]
        + ["--inject", "40", "--inject-frac", "1.0", "--procs", "1"]
    )
    assert rc == 0
    grid = Table.read(out / "grid.ecsv")
    assert list(np.flatnonzero(grid["hit"])) == [20]
    summary = json.loads((out / "summary.json").read_text())
    assert summary["n_hits"] == 1
    assert summary["n_hit_groups"] == 1
    assert summary["injection"]["runs"][0]["recovered_frac"] >= 0.9


def test_n_exposures_counts_bits_over_every_plane(tmp_path):
    p = tmp_path / "a.fits"
    _image(p, np.empty((0, 2)), np.empty(0), np.random.default_rng(0), con=True)
    w = _wcs()
    ra, dec = w.pixel_to_world_values([50.0, 350.0, -100.0], [200.0, 200.0, 200.0])
    assert list(tg.n_exposures(str(p), np.asarray(ra), np.asarray(dec))) == [1, 3, 0]


def test_main_writes_injections_and_splits_recovery_by_coverage(tmp_path):
    rng = np.random.default_rng(1)
    xy = _field(rng)
    flux = rng.uniform(2, 6, len(xy))
    w = _wcs()
    ra, dec = w.pixel_to_world_values(xy[:, 0], xy[:, 1])
    args = []
    for band in ("A", "B"):
        p1, p2 = tmp_path / f"{band}1.fits", tmp_path / f"{band}2.fits"
        _image(p1, xy, flux, rng, con=True)
        _image(p2, xy, flux, rng, con=True)
        args += ["--band", band, str(p1), str(p2)]
    Table({"ra": ra, "dec": dec}).write(tmp_path / "cat_a.ecsv", format="ascii.ecsv")
    out = tmp_path / "out"
    tg.main(
        args
        + ["--catalog", str(tmp_path / "cat_a.ecsv"), "--epoch2-catalogs", "--out", str(out)]
        + ["--inject", "40", "--inject-frac", "1.0", "0.5", "--procs", "1", "--coverage"]
    )
    grid = Table.read(out / "grid.ecsv")
    left = np.asarray(grid["ra"]) > w.pixel_to_world_values(N / 2, N / 2)[0]  # RA grows leftwards
    assert set(grid["A_e1_n_exp"][left]) == {1} and set(grid["B_e2_n_exp"][~left]) == {3}
    inj = Table.read(out / "injections.ecsv")
    assert len(inj) == 80 and set(inj["frac"]) == {1.0, 0.5}
    run = json.loads((out / "summary.json").read_text())["injection"]["runs"][0]
    first = inj[inj["frac"] == 1.0]
    assert run["recovered_frac"] == pytest.approx(np.mean(first["recovered"]))
    split = run["by_e1_n_exp_min"]
    assert split["0"]["n"] + split["1"]["n"] + split[">=2"]["n"] == 40
    assert split["1"]["n"] == int((first["e1_n_exp_min"] == 1).sum())
    assert np.all(first["A_inj_flux"] > 0) and set(first["B_pass"]) <= {-1, 0, 1}


def test_coverage_without_con_stops_before_measuring(tmp_path):
    p = tmp_path / "a.fits"
    _image(p, np.empty((0, 2)), np.empty(0), np.random.default_rng(0))
    Table({"ra": [24.35], "dec": [-8.45]}).write(tmp_path / "c.ecsv", format="ascii.ecsv")
    with pytest.raises(SystemExit, match="CON"):
        tg.main(
            ["--band", "A", str(p), str(p), "--catalog", str(tmp_path / "c.ecsv")]
            + ["--epoch2-catalogs", "--out", str(tmp_path / "o"), "--coverage", "--procs", "1"]
        )


def _colour_change_field(tmp_path, recentre_in):
    """Source 20 fades in band A and appears in band B, where a bright neighbour sits 8 px away."""
    tmp_path = tmp_path / recentre_in
    tmp_path.mkdir()
    rng = np.random.default_rng(3)
    xy = _field(rng)
    w = _wcs()
    ra, dec = w.pixel_to_world_values(xy[:, 0], xy[:, 1])
    flux = rng.uniform(2, 6, len(xy))
    nb = xy[20] + (8.0, 0.0)
    args = []
    for band in ("A", "B"):
        f1, f2 = flux.copy(), flux.copy()
        if band == "A":
            f1[20], f2[20] = 3.0, 1.0
        else:
            f1[20], f2[20] = 0.0, 3.0
        pts = np.vstack([xy, nb]) if band == "B" else xy
        p1, p2 = tmp_path / f"{band}1.fits", tmp_path / f"{band}2.fits"
        _image(p1, pts, np.append(f1, 30.0) if band == "B" else f1, rng)
        _image(p2, pts, np.append(f2, 30.0) if band == "B" else f2, rng)
        args += ["--band", band, str(p1), str(p2)]
    Table({"ra": ra, "dec": dec}).write(tmp_path / "cat_a.ecsv", format="ascii.ecsv")
    out = tmp_path / "out"
    tg.main(
        args
        + ["--catalog", str(tmp_path / "cat_a.ecsv"), "--epoch2-catalogs", "--out", str(out)]
        + ["--recentre-in", recentre_in, "--procs", "1"]
    )
    return Table.read(out / "grid.ecsv")


def test_origin_recentring_keeps_the_aperture_on_a_source_absent_in_one_band(tmp_path):
    grid = _colour_change_field(tmp_path, "origin")
    assert list(np.flatnonzero(grid["hit"])) == [20]
    assert grid["mixed"][20]
    assert grid["A_pass"][20] == -1 and grid["B_pass"][20] == 1
    # per-band recentring follows the neighbour in band B, epoch 1, and loses the source
    band = _colour_change_field(tmp_path, "band")
    assert not band["hit"][20]


def test_first_of_groups_keeps_one_per_close_pair():
    ra = np.array([24.0, 24.0 + 0.1 / 3600, 24.0 + 1 / 3600])
    assert list(tg.first_of_groups(ra, np.full(3, -8.0), 0.2)) == [True, False, True]
