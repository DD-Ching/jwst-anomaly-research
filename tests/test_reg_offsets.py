"""Unit tests for scripts/reg_offsets.py (offline, synthetic catalogues)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table


def _load(name: str):
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ro = _load("reg_offsets")


def _cat(ra, dec, mag, ext=None):
    t = Table()
    t["label"] = np.arange(len(ra))
    t["sky_centroid"] = SkyCoord(ra * u.deg, dec * u.deg)
    t["aper_total_abmag"] = mag
    t["is_extended"] = np.zeros(len(ra), bool) if ext is None else ext
    return t


def test_global_shift_recovered_and_held_out_passes():
    rng = np.random.default_rng(3)
    n = 600
    dec0 = -30.0
    x = rng.uniform(-90, 90, n)  # arcsec
    y = rng.uniform(-90, 90, n)
    ra = 3.5 + x / 3600 / np.cos(np.radians(dec0))
    dec = dec0 + y / 3600
    mag = rng.uniform(21.5, 25.0, n)
    ref = _cat(ra, dec, mag)
    sx, sy = 0.025, -0.012  # arcsec, the shift to recover
    noise = 0.008
    ra2 = ra + (sx + rng.normal(0, noise, n)) / 3600 / np.cos(np.radians(dec))
    dec2 = dec + (sy + rng.normal(0, noise, n)) / 3600
    tab = ro.controls(_cat(ra2, dec2, mag), ref)
    assert len(tab) > 100 and set(tab["set"]) == {"A", "B"}
    f = ro.fit(tab, "M0")
    assert abs(f["px"][0] - sx) < 0.002 and abs(f["py"][0] - sy) < 0.002
    b = ro.score(tab, f)
    assert max(abs(m) for m in b["median"]) < 0.003
    assert all(abs(r - noise) < 0.003 for r in b["rms"])


def test_controls_drop_extended_bright_and_crowded():
    ra = np.array([3.5, 3.5 + 0.0001, 3.51, 3.52, 3.515])
    dec = np.array([-30.0, -30.0, -30.01, -30.02, -30.005])
    mag = np.array([23.0, 23.0, 19.0, 23.0, 23.0])
    ext = np.array([False, False, False, True, False])
    # pad with a footprint ring far away so the hull edge test keeps the inner points
    ring = np.linspace(0, 2 * np.pi, 40, endpoint=False)
    ra = np.r_[ra, 3.51 + 0.03 * np.cos(ring)]
    dec = np.r_[dec, -30.01 + 0.03 * np.sin(ring)]
    mag = np.r_[mag, np.full(40, 30.0)]
    ext = np.r_[ext, np.zeros(40, bool)]
    tab = ro.controls(_cat(ra, dec, mag, ext), _cat(ra, dec, mag, ext))
    # 0/1 are a crowded pair, 2 is too bright, 3 is extended, the ring is too faint: only 4 stays
    assert len(tab) == 1 and abs(tab["ra"][0] - 3.515) < 1e-9
