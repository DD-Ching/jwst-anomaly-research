"""Unit tests for scripts/arc_diff_peers.py on synthetic maps (offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table
from astropy.wcs import WCS

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arc_diff_peers", _DIR / "arc_diff_peers.py")
adp = importlib.util.module_from_spec(_spec)
sys.modules["arc_diff_peers"] = adp
_spec.loader.exec_module(adp)

N = 201
SCALE = 0.063


def _wcs() -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [3.5, -30.4]
    w.wcs.crpix = [N // 2 + 1, N // 2 + 1]
    w.wcs.cdelt = [-SCALE / 3600, SCALE / 3600]
    return w


def _scene(peer_s: float, hit_s: float):
    """Epoch 1: a lattice of equal compact sources; S: peer_s at each, hit_s at one of them."""
    yy, xx = np.mgrid[:N, :N]
    e1 = np.full((N, N), 1e-3)
    s = np.zeros((N, N))
    rng = np.random.default_rng(1)
    centres = [(y, x) for y in range(20, N - 20, 16) for x in range(20, N - 20, 16)]
    for y, x in centres:
        e1 += 0.05 * np.exp(-((yy - y) ** 2 + (xx - x) ** 2) / (2 * 1.2**2))
        s[y, x] = peer_s + rng.normal(0, 1.0)
    hy, hx = centres[len(centres) // 2]
    s[hy, hx] = hit_s
    w = _wcs()
    ra, dec = w.all_pix2world(hx, hy, 0)
    hits = Table(rows=[{"hit_id": "h000", "ra": float(ra), "dec": float(dec)}])
    return s, e1, w, hits


def test_peer_offset_absorbs_systematic():
    # Every compact source has S ≈ −6 (program PSF mismatch); a hit at −6 is ordinary.
    s, e1, w, hits = _scene(-6.0, -6.2)
    t = adp.peer_table(s, e1, w, hits)
    assert t["n_peers"][0] > 50
    assert abs(t["peer_med_s"][0] + 6.0) < 0.5
    assert not t["peer_pass"][0]


def test_real_change_passes():
    s, e1, w, hits = _scene(-6.0, 8.0)
    t = adp.peer_table(s, e1, w, hits)
    assert t["z_peers"][0] > 3
    assert t["peer_pass"][0]


def test_no_peers_is_no_test():
    s, e1, w, hits = _scene(-6.0, -6.0)
    t = adp.peer_table(s, e1, w, hits, dex=1e-9, self_arcsec=10.0)
    assert t["n_peers"][0] == 0
    assert t["peer_pass"][0]


def test_stored_z_uses_nearest_tile_rows():
    t = Table(
        rows=[
            {"hit": "h0", "s_hit": 8.0, "log_peak": -1.0},
            {"hit": "h1", "s_hit": -5.0, "log_peak": -2.0},
        ]
    )
    ref = Table(
        rows=[
            {"tile": "a", "log_peak": -1.0, "peer_med_s": 0.0, "peer_mad_s": 2.0},
            {"tile": "a", "log_peak": -2.0, "peer_med_s": -4.0, "peer_mad_s": 1.0},
            {"tile": "b", "log_peak": -2.0, "peer_med_s": 50.0, "peer_mad_s": 1.0},
        ]
    )
    out = adp.stored_z(t, ref, "a", k=1)
    assert out["z_stored"][0] == 4.0 and out["stored_pass"][0]
    assert out["z_stored"][1] == -1.0 and not out["stored_pass"][1]


def test_stored_z_unknown_tile_raises():
    t = Table(rows=[{"hit": "h0", "s_hit": 8.0, "log_peak": -1.0}])
    ref = Table(rows=[{"tile": "a", "log_peak": -1.0, "peer_med_s": 0.0, "peer_mad_s": 2.0}])

    with pytest.raises(ValueError):
        adp.stored_z(t, ref, "zz")


def test_window_outside_is_empty_and_extreme_skips_nan():
    a = np.arange(25.0).reshape(5, 5)
    assert adp.window(a, -10, 2, 2).size == 0
    assert adp.window(a, 0, 0, 1).shape == (2, 2)
    w = np.array([[np.nan, 0.0], [-2.0, 1.0]])
    assert adp.signed_extreme(w) == -2.0
    assert np.isnan(adp.signed_extreme(np.full((2, 2), np.nan)))


def test_peaks_edge_px_zero_still_finds_peaks():
    e1 = np.full((40, 40), 1e-3)
    e1[20, 20] = 1.0
    y, x, _ = adp.peaks(e1, 3, 0)
    assert (20, 20) in set(zip(y.tolist(), x.tolist(), strict=True))
