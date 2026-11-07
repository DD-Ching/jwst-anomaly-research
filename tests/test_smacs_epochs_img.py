"""Frozen outcome rules and helpers of scripts/smacs_epochs_img.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import smacs_epochs_img as m  # noqa: E402


def _band(res, med=5.0, p95=20.0):
    return {
        "target_residual_mas": list(res),
        "target_residual_total_mas": float(np.hypot(*res)),
        "peer_residual_median_mas": med,
        "peer_residual_p95_mas": p95,
    }


def test_decide_moves_needs_both_bands_and_direction():
    assert m.decide({"f150w": _band((25, 0)), "f444w": _band((24, 5))}) == "moves: star"
    assert m.decide({"f150w": _band((25, 0)), "f444w": _band((0, 25))}) == "undetermined"
    assert m.decide({"f150w": _band((25, 0)), "f444w": _band((15, 0))}) == "undetermined"


def test_decide_no_motion_and_missing_band():
    assert m.decide({"f150w": _band((1, 1)), "f444w": _band((2, -2))}) == "no motion detected"
    assert m.decide({"f150w": _band((1, 1)), "f444w": _band((6, 0))}) == "undetermined"
    assert m.decide({"f150w": _band((1, 1)), "f444w": None}) == "undetermined"


def test_band_stats_min_peers_and_tie():
    peers = np.tile([10.0, 20.0], (m.MIN_PEERS, 1))
    st = m.band_stats(np.array([13.0, 24.0]), peers)
    assert st["tie_mas"] == [10.0, 20.0]
    assert np.isclose(st["target_residual_total_mas"], 5.0)
    assert m.band_stats(np.array([0.0, 0.0]), peers[:-1]) is None
    assert m.band_stats(None, peers) is None


def _wcs():
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [110.0, -73.0]
    w.wcs.crpix = [16.0, 16.0]
    w.wcs.cdelt = [-0.03 / 3600, 0.03 / 3600]
    return w


def test_centroid_recovers_injected_shift():
    w = _wcs()
    yy, xx = np.mgrid[:31, :31]
    img = 100 * np.exp(-((xx - 15.2) ** 2 + (yy - 14.9) ** 2) / (2 * 1.5**2))
    c0 = m.centroid_sky(img, w)
    c1 = m.centroid_sky(m.shifted(img, w, c0, (20.0, 0.0)), w)
    dra, ddec = c0.spherical_offsets_to(c1)
    assert abs(dra.to_value("mas") - 20.0) < 0.5
    assert abs(ddec.to_value("mas")) < 0.5
