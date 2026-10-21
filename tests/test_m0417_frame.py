"""Frozen rules of scripts/m0417_frame.py (offset measurement, gates F1 and F2)."""

import importlib.util
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

_SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_SCRIPTS))
_spec = importlib.util.spec_from_file_location("m0417_frame", _SCRIPTS / "m0417_frame.py")
mf = importlib.util.module_from_spec(_spec)
sys.modules["m0417_frame"] = mf
_spec.loader.exec_module(mf)


def _synthetic(n=200, shift=(-0.09, -0.27), noise=0.02, seed=1):
    rng = np.random.default_rng(seed)
    ra = 64.39 + rng.uniform(-0.02, 0.02, n)
    dec = -11.90 + rng.uniform(-0.02, 0.02, n)
    relics = Table({"id": np.arange(n), "ra": ra, "dec": dec})
    c = SkyCoord(ra * u.deg, dec * u.deg).spherical_offsets_by(
        (shift[0] + rng.normal(0, noise, n)) * u.arcsec,
        (shift[1] + rng.normal(0, noise, n)) * u.arcsec,
    )
    return relics, c


def test_measure_recovers_injected_offset():
    relics, jwst = _synthetic()
    m = mf.measure(relics, jwst)
    assert abs(m["offset"]["dra"] + 0.09) < 0.01
    assert abs(m["offset"]["ddec"] + 0.27) < 0.01
    assert m["offset"]["n"] == 200
    assert m["F1"]


def test_gate_f1_needs_count_precision_and_halves():
    ok = {"n": 30, "dra_err": 0.02, "ddec_err": 0.01}
    assert mf.gate_f1(ok, 0.03)
    assert not mf.gate_f1({**ok, "n": 29}, 0.0)
    assert not mf.gate_f1({**ok, "ddec_err": 0.021}, 0.0)
    assert not mf.gate_f1(ok, 0.031)


def test_gate_f2_threshold():
    assert mf.gate_f2(0.18)
    assert not mf.gate_f2(0.181)
