"""Offline tests for the signed-lag GW cells with sky-map classes."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

pytest.importorskip("astropy_healpix")

import e1_gw_directional as D  # noqa: E402
import e1_gw_signed as G  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402


def _setup(tmp_path):
    vec = D.pix_vectors()
    prob = np.array([D.gaussian_map(100.0, 20.0, 2.0, vec)], dtype=np.float32)
    path = tmp_path / "m.npz"
    np.savez(path, name=np.array(["GW000000_000000"]), prob=prob)
    m = D.GWMaps(["GW000000_000000"], path)
    gw = en.Sample(
        "GW",
        np.array([60000.0]),
        np.array([np.nan]),
        np.array([np.nan]),
        np.array([np.nan]),
        np.array([2023]),
    )
    return m, gw


def test_sign_counts_after_minus_before(tmp_path):
    m, gw = _setup(tmp_path)
    b = en.Sample(
        "GBM",
        60000.0 + np.array([2.0, 5.0, -3.0, 50.0]) / en.DAY,
        np.array([100.0, 100.5, 100.0, 280.0]),
        np.array([20.0, 20.0, 20.0, -20.0]),
        np.full(4, 0.5),
        np.full(4, 2023),
    )
    d = G.signed_channel(m, gw, gw.mjd, b)
    assert d[0, 0] == 1  # same, 0-10 s: two after, one before
    assert d[1, 2] == 1  # antipodal, 10-100 s: one after
    assert G.signed_channel(m, gw, gw.mjd, b, skip=[(0, 0)])[0, 0] == 0


def test_injection_after_gives_positive_d(tmp_path):
    m, gw = _setup(tmp_path)
    rng = np.random.default_rng(3)
    b = en.Sample(
        "ICECAT",
        60000.0 + rng.uniform(-20, 20, 30),
        rng.uniform(0, 360, 30),
        rng.uniform(-60, 60, 30),
        np.full(30, 0.5),
        np.full(30, 2023),
    )
    base = G.signed_channel(m, gw, gw.mjd, b)[0, 0]
    sb = D.inject(m, gw, b, 1, 0.0, 10.0, "same", rng, sign=1)
    assert G.signed_channel(m, gw, gw.mjd, sb)[0, 0] >= base
