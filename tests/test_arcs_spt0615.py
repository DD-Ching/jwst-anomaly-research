"""Tests for scripts/arcs_spt0615.py: RELICS BPZ conversion, catalogue reading and the astrometric
shift rule (synthetic tables; no network)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arcs_spt0615", _DIR / "arcs_spt0615.py")
asp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(asp)


def _bpz() -> Table:
    return Table(
        {
            "id": [1, 2, 3],
            "RA": [94.0, 94.001, 94.002],
            "Dec": [-57.78, -57.78, -57.78],
            "stel": [0.02, 0.95, 0.79],
            "zb": [2.0, 0.1, 3.5],
            "zbmin": [1.6, 0.0, 3.0],
            "zbmax": [2.4, 0.3, 4.0],
        }
    )


def test_bpz_zout_maps_columns_and_drops_stars():
    z = asp.bpz_zout(_bpz())
    assert list(z["id"]) == [1, 3]  # stel 0.95 dropped, 0.79 kept
    assert list(z["z_phot"]) == [2.0, 3.5]
    assert list(z["z160"]) == [1.6, 3.0]  # 95 % lower bound stands in for z160
    assert list(z["z840"]) == [2.4, 4.0]


def test_read_bpz_names_from_header(tmp_path):
    p = tmp_path / "cat.txt"
    p.write_text(
        "## comment\n#  1 id   Object ID\n#  2 RA   Right ascension\n#  3 zb   BPZ\n"
        "1 94.0 2.0\n2 94.1 0.5\n"
    )
    t = asp.read_bpz(p)
    assert t.colnames == ["id", "RA", "zb"]
    assert list(t["zb"]) == [2.0, 0.5]


def test_shift_only_above_threshold():
    z = asp.bpz_zout(_bpz())
    ra0, dec0 = np.array(z["ra"]), np.array(z["dec"])
    assert not asp.shift(z, {"dra_arcsec": 0.1, "ddec_arcsec": -0.05})
    assert np.array_equal(z["ra"], ra0)
    assert asp.shift(z, {"dra_arcsec": 0.0, "ddec_arcsec": 0.36})
    assert np.allclose((z["dec"] - dec0) * 3600, 0.36)
    assert np.allclose(z["ra"], ra0)


def test_shift_ra_is_on_sky_arcsec():
    z = asp.bpz_zout(_bpz())
    ra0 = np.array(z["ra"])
    asp.shift(z, {"dra_arcsec": 0.3, "ddec_arcsec": 0.0})
    dra = (z["ra"] - ra0) * np.cos(np.radians(z["dec"])) * 3600
    assert np.allclose(dra, 0.3)


def test_shift_refuses_undefined_offset():
    z = asp.bpz_zout(_bpz())
    with pytest.raises(SystemExit):
        asp.shift(z, {"n_pairs": 0, "dra_arcsec": float("nan"), "ddec_arcsec": float("nan")})
