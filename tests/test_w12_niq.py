"""Tests for scripts/w12_niq.py helpers (tiny synthetic tables, offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w12_niq", _DIR / "w12_niq.py")
niq = importlib.util.module_from_spec(_spec)
sys.modules["w12_niq"] = niq
_spec.loader.exec_module(niq)


def test_sdss_name_radec():
    ra, dec = niq.sdss_name_radec("J132236.41+105239.4")
    assert ra == pytest.approx(15 * (13 + 22 / 60 + 36.41 / 3600))
    assert dec == pytest.approx(10 + 52 / 60 + 39.4 / 3600)
    assert niq.sdss_name_radec("J212956.56-005152.4")[1] < 0
    assert np.isnan(niq.sdss_name_radec("not a name")[0])


@pytest.mark.parametrize(
    "comment, group",
    [
        ("No lensing object", "rejected"),
        ("No lens object, binary QSO", "rejected"),
        ("QSO pair, no lensing object", "rejected"),
        ("SDSS lens", "control"),
        ("Known lens", "control"),
        ("QSO+star", ""),
        ("", ""),
    ],
)
def test_sqls_group(comment, group):
    assert niq.sqls_group(comment) == group


def test_read_tsv_skips_units_and_dashes(tmp_path):
    f = tmp_path / "t.tsv"
    f.write_text("#comment\nName\tSep\n \tarcsec\n----\t---\nJ1\t 2.1\nJ2\n")
    t = niq.read_tsv(f)
    assert list(t["Name"]) == ["J1", "J2"] and list(t["Sep"]) == ["2.1", ""]


def _lemon(rows):
    return Table(
        rows=rows, names=("Name", "RAJ2000", "DEJ2000", "z", "Sep", "Class"), dtype=[str] * 6
    )


def _sqls(rows):
    return Table(rows=rows, names=("SDSS", "z", "theta", "Com"), dtype=[str] * 4)


def test_build_sample_groups_and_dedup():
    lemon = _lemon(
        [
            ("J1000+0100", "150.0", "1.0", "1.5", "2.2", "UQP"),
            ("J1100+0100", "165.0", "1.0", "2.0", "2.5", "lens"),
            ("J1200+0100", "180.0", "1.0", "1.0", "3.0", "QSO + star"),
        ]
    )
    # SQLS: one duplicate of the Lemon UQP (listed as a lens there), one new rejection
    sqls = _sqls(
        [
            ("J100000.00+010000.0", "1.5", "2.2", "SDSS lens"),
            ("J130000.00+020000.0", "", "2.4", "No lensing object"),
        ]
    )
    s = niq.build_sample({"J/MNRAS/520/3305/table1": lemon, "J/AJ/143/119/table4": sqls})
    assert len(s) == 3  # the Lemon "QSO + star" is dropped, the SQLS duplicate merged
    row = s[s["name"] == "J1000+0100"][0]
    assert row["group"] == "control" and "listed as lens elsewhere" in row["comment"]
    new = s[s["name"] == "J130000.00+020000.0"][0]
    assert new["group"] == "rejected" and new["ra"] == pytest.approx(195.0)
    assert np.allclose(s["theta_e"], s["sep_cat"] / 2)
    assert set(s["selection"]) == {"quasar"}


def test_pair_colour_difference():
    src = Table({"mag_g": [20.0, 21.0, 19.0], "mag_z": [19.5, 19.5, np.nan]})
    images = Table({"img1": [0, 0, 0], "img2": [1, 2, -1]})
    d = niq.pair_colour_difference(images, src)
    assert d[0] == pytest.approx(0.5 - 1.5) and np.isnan(d[1]) and np.isnan(d[2])


def test_binary_match():
    s = Table({"ra": [150.0, 10.0], "dec": [1.0, -5.0]})
    binq = Table({"_RA": ["150.0003", "30.0"], "_DE": ["1.0", "2.0"]}, dtype=[str, str])
    assert list(niq.binary_match(s, binq)) == [True, False]


def test_binary_match_sexagesimal_and_unparseable():
    s = Table({"ra": [150.0], "dec": [1.0]})
    binq = Table({"RA1": ["10:00:00.05"], "DE1": ["+01:00:00.5"]}, dtype=[str, str])
    assert list(niq.binary_match(s, binq)) == [True]
    with pytest.raises(ValueError):
        niq.binary_match(s, Table({"RA1": ["x"], "DE1": ["y"]}, dtype=[str, str]))
