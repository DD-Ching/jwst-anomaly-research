"""Tests for scripts/fade_calib.py and arc_diff pre-matching injection (synthetic, offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("fade_calib", _DIR / "fade_calib.py")
fc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fc)


def test_class_mags_and_cycle():
    assert fc.class_mags("vanish", 26.0) == (26.0, None)
    assert fc.class_mags("appear", 26.0) == (None, 26.0)
    assert fc.class_mags("partial", 26.0) == (26.0, 26.75)
    rows = fc.sample_rows([(1.0, 2.0)] * 12, [26.0, 27.0, 27.5])
    cells = {(r["cls"], r["mag"]) for r in rows}
    assert len(cells) == 12  # 4 classes x 3 mags each appear once in 12 rows


def test_sign_ok():
    assert fc.sign_ok("vanish", "-", {})
    assert fc.sign_ok("vanish", "mixed", {"F150W": -8.0, "F444W": 3.0})
    assert not fc.sign_ok("vanish", "mixed", {"F150W": -3.0, "F444W": 8.0})
    assert not fc.sign_ok("vanish", "+", {})
    assert fc.sign_ok("appear", "+", {})


def test_summarize_and_pass_rule():
    rows = []
    for i in range(10):
        ok = i < 8
        rows.append(
            {
                "cls": "vanish",
                "mag": "26.00",
                "hit": ok,
                "unflagged": ok,
                "morph_pass": ok,
                "sign_ok": ok,
                "review": ok,
                "review_sign": ok,
            }
        )
        rows.append(
            {
                "cls": "unchanged",
                "mag": "26.00",
                "hit": i == 0,
                "unflagged": False,
                "morph_pass": False,
                "sign_ok": False,
                "review": False,
                "review_sign": False,
            }
        )
    s = fc.summarize(rows)
    assert s["vanish"]["26.00"]["review_sign_frac"] == 0.8
    assert s["unchanged_false_hits"]["frac"] == 0.1
    assert fc.pass_rule(s)["passed"]
    lo, hi = fc.wilson(8, 10)
    assert lo < 0.8 < hi


def test_add_star_flux_and_mag_scale():
    import arc_diff as ad

    img = np.zeros((41, 41))
    ad.add_star(img, 20.0, 20.0, 1.5, 3.0)
    assert np.isclose(img.sum(), 3.0, rtol=1e-6)
    f1 = ad.mag_to_mjysr(26.0, 2.3e-14)
    f2 = ad.mag_to_mjysr(26.75, 2.3e-14)
    assert np.isclose(f2 / f1, 10 ** (-0.4 * 0.75))


def _ad():
    import arc_diff

    return arc_diff


def test_register_median_ignores_minority_tiles():
    ad = _ad()
    rng = np.random.default_rng(1)
    a = ndimage_gauss(rng.normal(size=(128, 128)))
    b = np.roll(a, (1, -1), axis=(0, 1)).copy()
    # one tile (of 16) dominated by a bright feature shifted differently
    b[:32, :32] = np.roll(a[:32, :32], (0, 2), axis=(0, 1)) * 50
    a2 = a.copy()
    a2[:32, :32] *= 50
    dy, dx, n = ad.register_median(a2, b, 32)
    assert n == 16
    assert abs(dy + 1) < 0.3 and abs(dx - 1) < 0.3


def ndimage_gauss(x):
    from scipy import ndimage

    return ndimage.gaussian_filter(x, 1.5)


def test_add_stamp_conserves_flux_and_centres():
    ad = _ad()
    psf = ad.gaussian_star((15, 15), 7, 7, 1.2, 1.0)
    img = np.zeros((40, 40))
    ad.add_stamp(img, 20.3, 18.6, psf, 5.0)
    assert np.isclose(img.sum(), 5.0, rtol=2e-3)
    yy, xx = np.mgrid[0:40, 0:40]
    assert abs((img * yy).sum() / img.sum() - 20.3) < 0.02
    assert abs((img * xx).sum() / img.sum() - 18.6) < 0.02
    edge = np.zeros((40, 40))
    ad.add_stamp(edge, 1.0, 1.0, psf, 1.0)  # clipped at the border, no error
    assert 0 < edge.sum() < 1.0


def test_empirical_psfs_per_epoch_width():
    ad = _ad()
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [10.0, -40.0]
    w.wcs.crpix = [100.0, 100.0]
    w.wcs.cdelt = [-0.03 / 3600, 0.03 / 3600]
    rng = np.random.default_rng(2)
    boxes = []
    for sig in (1.0, 1.6):
        sci = rng.normal(0, 0.01, (200, 200))
        for y, x in ((40, 40), (40, 150), (150, 60), (120, 160), (90, 100)):
            ad.add_star(sci, y + 0.2, x - 0.3, sig, 100.0)
        boxes.append({"sci": sci, "err": np.full_like(sci, 0.01), "wcs": w, "scale": 0.03})
    p1, p2, used = ad.empirical_psfs(boxes[0], boxes[1], 0.05, 0.3, n=5)
    assert len(used) == 5 and np.isclose(p1.sum(), 1) and np.isclose(p2.sum(), 1)
    assert p1.max() > 1.5 * p2.max()  # epoch 1 sharper, as injected


def test_cached_returns_fresh_copies(tmp_path):
    ad = _ad()
    calls = []

    def boxfn(n):
        calls.append(n)
        return {"sci": np.zeros(n)}

    a = ad.cached(str(tmp_path), boxfn, 3)
    a["sci"] += 1  # an injection on the returned copy must not reach the cache
    b = ad.cached(str(tmp_path), boxfn, 3)
    assert calls == [3] and b["sci"].sum() == 0
    assert ad.cached(str(tmp_path), boxfn, 4)["sci"].size == 4  # every argument is in the key
    assert ad.cached(None, boxfn, 2)["sci"].size == 2


def test_register_median_drops_empty_tiles():
    ad = _ad()
    rng = np.random.default_rng(3)
    a = ndimage_gauss(rng.normal(size=(128, 128)))
    b = np.roll(a, (1, -1), axis=(0, 1)).copy()
    b[:, :64] = np.nan  # half the box outside epoch 2's footprint
    dy, dx, n = ad.register_median(a, b, 32)
    assert n == 8
    assert abs(dy + 1) < 0.3 and abs(dx - 1) < 0.3


def _tan(crpix, crval=(15.73, -49.253), rot=0.0):
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = list(crval)
    w.wcs.crpix = list(crpix)
    c, s = np.cos(np.radians(rot)), np.sin(np.radians(rot))
    d = 0.03 / 3600
    w.wcs.cd = np.array([[-d * c, d * s], [d * s, d * c]])
    return w


def test_epoch2_pixel_follows_registration():
    ad = _ad()
    w1 = _tan((100.0, 100.0))
    w2 = _tan((110.0, 95.0), rot=0.84)
    ra, dec = w1.all_pix2world(60.0, 70.0, 0)
    y2, x2 = ad.epoch2_pixel(w1, w2, float(ra), float(dec), 0.0, 0.0)
    xr, yr = w2.all_world2pix(ra, dec, 0)
    assert np.allclose([y2, x2], [yr, xr], atol=1e-6)  # zero shift = WCS position
    y2, x2 = ad.epoch2_pixel(w1, w1, float(ra), float(dec), 0.7, -0.57)
    assert np.allclose([y2, x2], [70.0 - 0.7, 60.0 + 0.57], atol=1e-6)
    # after the chain's shift by +(dy, dx), the epoch-2 source sits back on the epoch-1 pixel
    img = np.zeros((140, 140))
    ad.add_star(img, y2, x2, 1.2, 1.0)
    moved = ad.shift(img, 0.7, -0.57)
    yy, xx = np.mgrid[0:140, 0:140]
    m = np.nan_to_num(moved)
    assert abs((m * yy).sum() / m.sum() - 70.0) < 0.05
    assert abs((m * xx).sum() / m.sum() - 60.0) < 0.05
