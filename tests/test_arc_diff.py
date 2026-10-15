"""Unit tests for scripts/arc_diff.py on synthetic images (offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "arc_diff.py"
_spec = importlib.util.spec_from_file_location("arc_diff", _PATH)
ad = importlib.util.module_from_spec(_spec)
sys.modules["arc_diff"] = ad
_spec.loader.exec_module(ad)


def _scene(rng: np.random.Generator, n: int = 128) -> np.ndarray:
    img = np.zeros((n, n))
    for _ in range(40):
        img += ad.gaussian_star((n, n), *rng.uniform(10, n - 10, 2), rng.uniform(1.0, 3.0), 50.0)
    return img


def test_register_recovers_subpixel_shift():
    rng = np.random.default_rng(1)
    a = _scene(rng)
    b = ad.shift(a, -1.3, 0.6)  # b is a moved by (-1.3, +0.6)
    dy, dx = ad.register(a, b)
    assert abs(dy - 1.3) < 0.1 and abs(dx + 0.6) < 0.1


def test_psf_match_blurs_the_sharper_epoch():
    rng = np.random.default_rng(2)
    e1 = _scene(rng)
    e2 = ad._blur(e1, 1.0)  # epoch 2 broader
    s, k, m1, m2 = ad.psf_match(e1, 1.03 * e2, np.round(np.linspace(-1.5, 1.5, 13), 3))
    assert 0.75 <= s <= 1.25  # blur epoch 1
    assert abs(k - 1.03) < 0.01
    assert ad._robust_std(m2 - m1) < 0.2 * ad._robust_std(e2 - e1)
    assert ad._robust_std(m2 - m1) < 0.2 * ad._robust_std(1.03 * e2 - e1)


def test_point_source_detected_and_noise_calibrated():
    rng = np.random.default_rng(3)
    n = 160
    noise = 1.0
    e1 = rng.normal(0, noise, (n, n))
    e2 = rng.normal(0, noise, (n, n))
    sig = 1.5
    flux = 12 * noise * np.sqrt(4 * np.pi) * sig  # ~12 sigma matched-filter
    e2 += ad.gaussian_star((n, n), 80.3, 70.6, sig, flux)
    var = np.full((n, n), 2 * noise**2)
    s, d = ad.significance(e2 - e1, var, sig, 21)
    assert 0.8 < ad._robust_std(s) < 1.2
    pk = ad.peaks(s, 5.0, 3)
    assert any(abs(y - 80) <= 1 and abs(x - 71) <= 1 and v > 8 for y, x, v in pk)
    assert len(pk) <= 3  # pure noise elsewhere rarely reaches 5 sigma


def test_group_links_bands_within_radius():
    d = [
        {"band": "F150W", "ra": 64.0, "dec": -24.0, "sigma": 6.0},
        {"band": "F277W", "ra": 64.0 + 0.05 / 3600, "dec": -24.0, "sigma": 7.0},
        {"band": "F356W", "ra": 64.0, "dec": -24.0 + 1.0 / 3600, "sigma": -6.0},
    ]
    g = sorted(ad.group(d, 0.1), key=len)
    assert [len(x) for x in g] == [1, 2]


def test_mag_to_flux_roundtrip():
    pixar = 2.3e-14
    f = ad.mag_to_mjysr(25.0, pixar)
    mag = -2.5 * np.log10(f * 1e6 * pixar / 3631.0)
    assert abs(mag - 25.0) < 1e-9


def test_local_shift_removes_a_varying_offset():
    rng = np.random.default_rng(4)
    a = _scene(rng, 192)
    # epoch 2: left half moved by +0.6 px in x, right half by -0.4 px (a smooth-ish warp)
    left = ad.shift(a, 0.0, 0.6)
    right = ad.shift(a, 0.0, -0.4)
    w = np.clip((np.arange(192) - 64) / 64.0, 0, 1)[None, :]
    b = (1 - w) * left + w * right
    glob = ad.shift(b, *ad.register(a, b))
    loc, loc2 = ad.local_shift(a, glob, 48, glob)
    assert np.allclose(loc, loc2, equal_nan=True)
    inner = (slice(20, -20), slice(20, -20))
    assert ad._robust_std((loc - a)[inner]) < 0.6 * ad._robust_std((glob - a)[inner])


def test_kernel_match_removes_psf_and_scale_difference_but_keeps_a_transient():
    rng = np.random.default_rng(5)
    n = 160
    e1 = _scene(rng, n)
    # epoch 2: broader, shifted by a fraction of a pixel, 3 % brighter, plus a transient
    e2 = 1.03 * ad.shift(ad._blur(e1, 1.2), 0.3, -0.2)
    e2 += ad.gaussian_star((n, n), 81.0, 77.0, 1.6, 40.0)
    e1n = e1 + rng.normal(0, 0.01, (n, n))
    e2n = e2 + rng.normal(0, 0.01, (n, n))
    v = np.full((n, n), 0.01**2)
    m1, m2, var, k, scale, conv = ad.kernel_match(e1n, e2n, v, v, 1.0, 4)
    assert conv == 1
    assert abs(scale - 1.03) < 0.02
    d = m2 - m1
    inner = (slice(10, -10), slice(10, -10))
    assert ad._robust_std(d[inner]) < 3 * 0.01 * np.sqrt(2)
    assert d[81, 77] > 0.5 * ad.gaussian_star((n, n), 81.0, 77.0, 1.6, 40.0)[81, 77]
    assert np.all(var > 0)


def test_compact_stamps_cover_point_sources_and_skip_masked_cores():
    rng = np.random.default_rng(6)
    n = 160
    img = rng.normal(0, 0.01, (n, n))
    img += ad.gaussian_star((n, n), 40.0, 40.0, 1.2, 30.0)
    img += ad.gaussian_star((n, n), 120.0, 100.0, 1.2, 30.0)
    img[118:123, 98:103] = np.nan  # saturated core
    st = ad.compact_stamps(img, 4)
    assert st[40, 40] and not st[120, 100]


def test_kernel_match_with_stamps_keeps_a_transient_on_a_compact_source():
    rng = np.random.default_rng(7)
    n = 192
    e1 = _scene(rng, n)
    e2 = ad._blur(e1, 1.0)
    y, x = 100.0, 90.0
    e1 += ad.gaussian_star((n, n), y, x, 1.5, 60.0)  # a compact host in both epochs
    e2 += ad.gaussian_star((n, n), y, x, np.hypot(1.5, 1.0), 60.0)
    tr = ad.gaussian_star((n, n), y, x, np.hypot(1.5, 1.0), 30.0)  # it brightens by 50 %
    e1n = e1 + rng.normal(0, 0.01, (n, n))
    e2n = e2 + tr + rng.normal(0, 0.01, (n, n))
    v = np.full((n, n), 0.01**2)
    m1, m2, _, _, _, _ = ad.kernel_match(e1n, e2n, v, v, 1.0, 4, stamps=True)
    d = m2 - m1
    assert d[100, 90] > 0.7 * tr[100, 90]


def test_add_star_matches_gaussian_star_total_flux():
    img = np.zeros((60, 60))
    ad.add_star(img, 30.2, 29.7, 1.5, 100.0)
    assert np.isclose(img.sum(), 100.0, rtol=1e-6)
    ref = ad.gaussian_star(img.shape, 30.2, 29.7, 1.5, 100.0)
    assert np.allclose(img, ref, rtol=0, atol=1e-5)


def test_grid_positions_bright_separated_and_inside_margin():
    from astropy.wcs import WCS

    rng = np.random.default_rng(0)
    img = _scene(rng, 200)
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [64.0, -24.0]
    w.wcs.crpix = [100, 100]
    w.wcs.cdelt = [-0.04 / 3600, 0.04 / 3600]
    pos = ad.grid_positions(img, w, 0.04, 30, 0.6, 0.8, 90.0, 1)
    assert 5 <= len(pos) <= 30
    xy = np.array([w.all_world2pix(ra, dec, 0) for ra, dec in pos])
    thr = np.percentile(img, 90)
    for x, y in xy:
        assert 20 - 1 <= x <= 180 and 20 - 1 <= y <= 180
        assert img[int(round(y)), int(round(x))] > thr * 0.5
    d = np.hypot(*(xy[:, None, :] - xy[None, :, :]).transpose(2, 0, 1))
    assert d[np.triu_indices(len(xy), 1)].min() * 0.04 >= 0.6 - 0.04


def _i2d(path, sci, err, scale):
    from astropy.io import fits
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [64.0, -24.0]
    w.wcs.crpix = [sci.shape[1] / 2 + 0.5, sci.shape[0] / 2 + 0.5]
    w.wcs.cdelt = [-scale / 3600, scale / 3600]
    h = w.to_header()
    h["PIXAR_SR"] = (scale / 206264.806) ** 2
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(sci.astype(np.float32), h, name="SCI"),
            fits.ImageHDU(err.astype(np.float32), h, name="ERR"),
        ]
    ).writeto(path)


def test_inject_grid_runs_through_the_whole_chain(tmp_path):
    rng = np.random.default_rng(5)
    n, scale = 400, 0.063
    args = []
    for band in ("F277W", "F356W"):
        sky = _scene(rng, n) * 0.05
        paths = []
        for ep in (1, 2):
            p = tmp_path / f"{band}_{ep}.fits"
            _i2d(p, sky + rng.normal(0, 0.002, (n, n)), np.full((n, n), 0.002), scale)
            paths.append(str(p))
        args += ["--band", band, *paths]
    out = tmp_path / "out"
    rc = ad.main(
        args
        + ["--ra", "64.0", "--dec", "-24.0", "--half-arcsec", "10", "--tile-arcsec", "0"]
        + ["--inject-grid", "12", "--inject-grid-mags", "24.0", "--procs", "1"]
        + ["--out", str(out), "--maps", str(tmp_path / "m.npz")]
    )
    assert rc == 0
    import csv
    import json

    rows = list(csv.DictReader(open(out / "injected.csv")))
    assert len(rows) >= 6
    assert all(r["hit"] for r in rows)  # bright injections are all recovered as hits
    summ = json.loads((out / "summary.json").read_text())
    assert summ["inject_grid"]["n_hit"] == len(rows)
