"""Tests for scripts/w1_streak_mask.py (synthetic inputs only; no acstools needed)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_streak_mask as wsm  # noqa: E402


def _wcs(scale_arcsec: float = 0.25) -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [150.0, 2.0]
    w.wcs.crpix = [1.0, 1.0]
    w.wcs.cdelt = [-scale_arcsec / 3600.0, scale_arcsec / 3600.0]
    return w


def test_i2d_uri_from_catalogue_name():
    stem = "jw01208-o010_t005_nircam_clear-f200w"
    cat = f"cache/mast/{stem}/{stem}_cat.ecsv"
    assert wsm.i2d_uri(cat) == (
        "s3://stpubdata/jwst/public/jw01208/L3/t/o010/jw01208-o010_t005_nircam_clear-f200w_i2d.fits"
    )


def test_block_median_ignores_zero_fill():
    a = np.ones((8, 8), np.float32)
    a[:4, :4] = 0.0  # off-mosaic block
    a[4:, 4:] = 3.0
    out = wsm.block_median(a, 4)
    assert np.isnan(out[0, 0]) and out[0, 1] == 1.0 and out[1, 1] == 3.0


def test_normalise_clips_zeroes_edges_and_rebins():
    rng = np.random.default_rng(0)
    d = rng.normal(10.0, 1.0, (80, 80))
    d[40, 40] = 1e4  # a bright pixel: clipped, not removed
    d[:4, :] = np.nan  # off mosaic
    z = wsm.normalise(d)
    assert z.shape == (40, 40)
    assert 0 < z[20, 20] <= wsm.CLIP_NSIG / 4 + 1  # one clipped pixel in a 2 x 2 mean
    assert (z[: (4 + wsm.EDGE_ERODE) // 2, :] == 0).all()  # eroded footprint edge
    assert abs(np.median(z[20:30, 20:30])) < 0.5


def test_veto_rows_within_tolerance_of_mask():
    w = _wcs(0.25)
    mask = np.zeros((100, 100), bool)
    mask[:, 50] = True  # a vertical streak one pixel wide
    on = w.pixel_to_world_values(50, 30)
    near = w.pixel_to_world_values(53, 30)  # 0.75" away: vetoed at 1"
    far = w.pixel_to_world_values(55, 30)  # 1.25" away: kept
    outside = w.pixel_to_world_values(500, 30)  # off the mask image: kept
    ra = [on[0], near[0], far[0], outside[0]]
    dec = [on[1], near[1], far[1], outside[1]]
    assert wsm.veto_rows(ra, dec, mask, w, 1.0).tolist() == [True, True, False, False]


def test_veto_rows_empty_mask_vetoes_nothing():
    w = _wcs()
    assert not wsm.veto_rows([150.0], [2.0], np.zeros((10, 10), bool), w).any()


def test_mask_round_trip(tmp_path):
    w = _wcs()
    m = np.zeros((12, 9), bool)
    m[3, :] = True
    p = tmp_path / "x_streak_mask.fits.gz"
    wsm.write_mask(p, m, w)
    m2, w2 = wsm.read_mask(p)
    assert (m2 == m).all()
    assert np.allclose(w2.wcs.cdelt, w.wcs.cdelt)


def test_known_case_and_adoption_rule():
    assert wsm.adopt(wsm.known_case(True, False))
    assert not wsm.adopt(wsm.known_case(False, False))  # K1 fails
    assert not wsm.adopt(wsm.known_case(True, True))  # MACS1423 still passes: K2 fails


def test_evaluate_reads_k1_k2(tmp_path):
    w = _wcs()
    m = np.zeros((40, 40), bool)
    m[:, 20] = True
    wsm.write_mask(tmp_path / "macs1423_streak_mask.fits.gz", m, w)
    ra, dec = w.pixel_to_world_values(20, 10)
    hit = tmp_path / "hit.json"
    hit.write_text(json.dumps({"screen": {"s_max_radec": [float(ra), float(dec)]}}))
    new = tmp_path / "new.json"
    scr = {"s_max": 3.0, "p_random_max": 0.3, "b_max": 3.5, "passes_eb_rule": False}
    new.write_text(
        json.dumps({"field": "macs1423", "assumptions": {"streak_veto": {}}, "screen": scr})
    )
    res = wsm.evaluate(tmp_path, new, hit)
    assert res["K1"] and res["K2"] and res["adopted"]


def test_rewriting_an_identical_mask_gives_the_same_file(tmp_path):
    import time

    mask = np.zeros((10, 10), bool)
    mask[3, :] = True
    a, b = tmp_path / "a_streak_mask.fits.gz", tmp_path / "b_streak_mask.fits.gz"
    wsm.write_mask(a, mask, _wcs())
    time.sleep(1.1)  # a gzip timestamp would differ by now
    wsm.write_mask(b, mask, _wcs())
    ia, ib = wsm.mask_identity(a), wsm.mask_identity(b)
    assert ia["sha256"] == ib["sha256"] and ia["data_sha256"] == ib["data_sha256"]
    back, _ = wsm.read_mask(a)
    np.testing.assert_array_equal(back, mask)
    assert ia["mask_pixels"] == 10


def test_find_trails_passes_the_threshold(monkeypatch):
    """--mrt-threshold reaches TrailFinder (results/w1_streak_rowdiag/PLAN.md); stubbed."""
    import types

    seen = {}

    class Stub:
        def __init__(self, z, **kw):
            seen.update(kw)
            self.mask, self.source_list = None, None

        def run_all(self):
            pass

    mod = types.ModuleType("acstools.findsat_mrt")
    mod.TrailFinder = Stub
    monkeypatch.setitem(sys.modules, "acstools", types.ModuleType("acstools"))
    monkeypatch.setitem(sys.modules, "acstools.findsat_mrt", mod)
    mask, cat = wsm.find_trails(np.zeros((4, 4)), processes=1, threshold=2.5)
    assert seen["threshold"] == 2.5 and not mask.any() and cat is None
    wsm.find_trails(np.zeros((4, 4)), processes=1)
    assert seen["threshold"] == wsm.MRT_THRESHOLD
