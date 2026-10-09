"""Offline tests for the E1 GW sky-map reduction and the directional GW pair classes."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

pytest.importorskip("astropy_healpix")

import e1_gw_directional as D  # noqa: E402
import e1_gw_skymaps as S  # noqa: E402
from astropy.table import Table  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402


def test_flat_nested_map_is_degraded_and_normalised():
    nside = 64
    prob = np.zeros(12 * nside**2)
    prob[4 * 7 : 4 * 7 + 4] = 0.25  # the four nside-64 children of nside-32 pixel 7
    t = Table({"PROB": prob})
    out = S.to_nested32(t.as_array(), {"ORDERING": "NESTED"})
    assert out.shape == (S.NPIX,)
    assert out[7] == pytest.approx(1.0)


def test_ring_map_is_reordered():
    from astropy_healpix.core import nested_to_ring

    prob = np.zeros(S.NPIX)
    prob[nested_to_ring(np.array([123]), S.NSIDE)[0]] = 1.0
    out = S.to_nested32(Table({"PROB": prob}).as_array(), {"ORDERING": "RING"})
    assert np.argmax(out) == 123


def test_multi_order_map_finer_and_coarser_pixels():
    # one order-6 pixel (child of order-5 pixel 10) and one order-4 pixel (parent of 4*3..4*3+3)
    o6, o4 = 6, 4
    uniq = np.array([4 * 4**o6 + 4 * 10, 4 * 4**o4 + 3])
    area = 4 * np.pi / (12 * 4.0 ** np.array([o6, o4]))
    dens = np.array([0.5, 0.5]) / area
    out = S.to_nested32(Table({"UNIQ": uniq, "PROBDENSITY": dens}).as_array(), {})
    assert out[10] == pytest.approx(0.5)
    assert out[12:16] == pytest.approx(np.full(4, 0.125))


def test_credible_level_and_region_distance():
    vec = D.pix_vectors()
    p = D.gaussian_map(100.0, 20.0, 3.0, vec)
    cl = D.credible_level(p)
    assert cl.min() == pytest.approx(p.max())
    assert cl.max() == pytest.approx(1.0)
    d = D.region_distance(cl <= 0.9, vec)
    assert d[cl <= 0.9].max() == pytest.approx(0.0, abs=0.05)
    assert d.max() > 150


def _maps(tmp_path, centres, sigma=2.0):
    vec = D.pix_vectors()
    prob = np.array([D.gaussian_map(r, d, sigma, vec) for r, d in centres], dtype=np.float32)
    names = np.array([f"GW{k:06d}_000000" for k in range(len(centres))])
    path = tmp_path / "m.npz"
    np.savez(path, name=names, prob=prob)
    return D.GWMaps(list(names), path), names


def test_classify_x_same_wide_antipodal(tmp_path):
    m, _ = _maps(tmp_path, [(100.0, 20.0)])
    rot = np.zeros(1)
    gi = np.zeros(3, int)
    ra = np.array([100.0, 160.0, 280.0])
    dec = np.array([20.0, -40.0, -20.0])
    c = D.classify_x(m, gi, rot, ra, dec, np.full(3, 0.5))
    assert c[0] & 1 and not c[0] & 2  # same
    assert c[1] & 2 and not c[1] & 1  # wide
    assert c[2] & 4 and c[2] & 2  # antipodal (also outside the 99 % region)


def test_rotation_follows_hour_angle(tmp_path):
    # Shift the GW event by 6 sidereal hours: a partner at RA + 90 deg is now "same".
    m, _ = _maps(tmp_path, [(100.0, 20.0)])
    gw = en.Sample(
        "GW",
        np.array([60000.0]),
        np.array([np.nan]),
        np.array([np.nan]),
        np.array([np.nan]),
        np.array([2023]),
    )
    orig = gw.mjd.copy()
    moved = gw.with_times(gw.mjd + en.SIDEREAL_DAY / 4 / en.DAY)
    rot = D.rotation(moved, orig)
    c = D.classify_x(m, np.array([0]), rot, np.array([190.0]), np.array([20.0]), np.array([0.5]))
    assert c[0] & 1


def test_gw_gw_overlap_and_disjoint(tmp_path):
    m, _ = _maps(tmp_path, [(100.0, 20.0), (101.0, 21.0), (250.0, -60.0), (280.0, -20.0)])
    rot = np.zeros(4)
    c = D.classify_gw_gw(m, np.array([0, 0, 0]), np.array([1, 2, 3]), rot)
    assert c[0] & 1
    assert c[1] & 2 and not c[1] & 1
    assert c[2] & 4


def test_injected_same_pairs_are_counted(tmp_path):
    m, _ = _maps(tmp_path, [(100.0, 20.0)])
    gw = en.Sample(
        "GW",
        np.array([60000.0]),
        np.array([np.nan]),
        np.array([np.nan]),
        np.array([np.nan]),
        np.array([2023]),
    )
    rng = np.random.default_rng(1)
    b = en.Sample(
        "ICECAT",
        60000.0 + rng.uniform(-20, 20, 50),
        rng.uniform(0, 360, 50),
        rng.uniform(-60, 60, 50),
        np.full(50, 0.5),
        np.full(50, 2023),
    )
    base = D.count_channel(m, gw, gw.mjd, np.zeros(1), b)[0, 0]
    hits = 0
    for _ in range(20):
        sb = D.inject(m, gw, b, 1, 0.0, 10.0, "same", rng)
        hits += D.count_channel(m, gw, gw.mjd, np.zeros(1), sb)[0, 0] == base + 1
    # a counterpart drawn from the map lies outside the 90 % region about 10 % of the time
    assert hits >= 15
