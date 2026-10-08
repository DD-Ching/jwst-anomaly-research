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
        ("SDSS lens; QSO pair nearby", "control"),
        ("QSO pair", "rejected"),
        ("Binary QSO (z=0.799, 0.799)", "rejected"),
        ("QSO+star", "nonpair"),
        ("QSO pair (different SED)", "nonpair"),
        ("QSO + unknown lens candidate", "nonpair"),
        ("not a binary", ""),
        ("not a known lens", ""),
        ("known lensed quasar", "control"),
        ("Different SED, not QSO", "nonpair"),
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
        rows=[r + ("", "") for r in rows],
        names=("Name", "RAJ2000", "DEJ2000", "z", "Sep", "Class", "z2", "f_z"),
        dtype=[str] * 8,
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
    # the system's kept row is its lens entry, with the Lemon UQP noted in the comment
    row = s[s["name"] == "J100000.00+010000.0"][0]
    assert row["group"] == "control" and "Lemon2023: UQP" in row["comment"]
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


def test_build_sample_drops_wide_and_component_rows():
    sqls = _sqls(
        [
            ("J100000.00+010000.0", "1.5", "19.3", "No lensing object"),
            ("J103000.00+010000.0", "1.5", "", "No lensing object"),
            ("J110000.00+010000.0", "1.5", "2.3", "SDSS Lens (component A)"),
            ("J120000.00+010000.0", "1.5", "2.3", "No lensing object"),
        ]
    )
    s = niq.build_sample({"J/AJ/143/119/table4": sqls})
    assert list(s["name"]) == ["J120000.00+010000.0"]
    assert s.meta["dropped"]["sep_cat > 3.0 arcsec"] == 1
    assert s.meta["dropped"]["no catalogued separation"] == 1


def test_dedup_is_transitive():
    # A (rejected) - B (control, 2" from A) - C (rejected, 2" from B, 4" from A)
    lemon = _lemon([("A", "150.0", "1.0", "1.5", "2.2", "UQP")])
    sqls = _sqls(
        [
            ("J100000.26+010000.0", "1.5", "2.2", "No lensing object"),  # C: 3.9" from A
            ("J100000.13+010000.0", "1.5", "2.2", "SDSS lens"),  # B: 1.95" east of A, after C
        ]
    )
    s = niq.build_sample({"J/MNRAS/520/3305/table1": lemon, "J/AJ/143/119/table4": sqls})
    groups = dict(zip(s["name"], s["group"], strict=True))
    # transitive: A, B and C are one system, which holds a catalogued lens
    assert list(groups.items()) == [("J100000.13+010000.0", "control")]


def test_data_sha256_ignores_header_time():
    a = b"#Date: 2026-10-08T21:56:54\nName\tSep\nJ1\t2.1\n"
    b = b"#Date: 2026-10-09T08:00:00\nName\tSep\nJ1\t2.1\n"
    assert niq.data_sha256(a) == niq.data_sha256(b) != niq.data_sha256(b + b"J2\t3\n")


def test_check_pin_refuses_mismatch_and_missing():
    niq.check_pin("t", "abc", None)  # --repin
    niq.check_pin("t", "abc", {"t": "abc"})
    with pytest.raises(RuntimeError):
        niq.check_pin("t", "abc", {"t": "def"})
    with pytest.raises(RuntimeError):
        niq.check_pin("t", "abc", {})


def test_redshift_pair_and_difference():
    assert niq.comment_z_pair("QSO pair (z=1.686, 1.600)") == (1.686, 1.6)
    assert niq.comment_z_pair("Binary QSO (z = 0.799, 0.799)") == (0.799, 0.799)
    assert np.isnan(niq.comment_z_pair("No lensing object")[0])
    d = niq.different_redshift([1.686, 0.827, 1.0], [1.600, 0.824, np.nan])
    assert list(d) == [True, False, False]


def test_lemon_z2_used_only_for_a_second_quasar():
    t = Table(
        rows=[
            ("Q", "150.0", "1.0", "1.5", "2.2", "QSO pair", "1.40", "", " ", "  "),  # real format
            ("L", "160.0", "1.0", "1.5", "2.2", "lens", "0.40", "z_lens=", " ", "  "),
            ("F", "170.0", "1.0", "1.5", "2.2", "UQP", "1.45", "", "?", "  "),  # flagged
        ],
        names=("Name", "RAJ2000", "DEJ2000", "z", "Sep", "Class", "z2", "n_z2", "f_z2", "f_z"),
        dtype=[str] * 10,
    )
    s = niq.build_sample({"J/MNRAS/520/3305/table1": t})
    z2 = dict(zip(s["name"], s["z2"], strict=True))
    assert z2["Q"] == pytest.approx(1.4) and np.isnan(z2["L"]) and np.isnan(z2["F"])


def test_wide_lens_keeps_its_system_out_of_the_rejected_sample():
    lemon = _lemon([("A", "150.0", "1.0", "1.5", "2.2", "UQP")])
    sqls = _sqls(
        [
            ("J100000.13+010000.0", "1.5", "3.4", "SDSS lens"),  # wide lens 1.95" from A
            ("J110000.00+010000.0", "1.5", "2.4", "QSO pair (z=1.686, 1.600)"),
            ("J110000.13+010000.0", "1.5", "", "No lensing object"),
        ]
    )
    s = niq.build_sample({"J/MNRAS/520/3305/table1": lemon, "J/AJ/143/119/table4": sqls})
    # A's system is a lens (dropped with the wide lens row), never a rejection
    assert "A" not in s["name"]
    row = s[s["name"] == "J110000.00+010000.0"][0]
    assert row["different_z"] and "No lensing object" in row["comment"]


def test_check_vizier_refuses_errors_and_empty(tmp_path):
    ok = b"#INFO x\nName\n \n----\nJ1\n"
    niq.check_vizier(ok, "t", tmp_path)
    with pytest.raises(RuntimeError):
        niq.check_vizier(b"#INFO QUERY_STATUS=ERROR\nName\n \n----\nJ1\n", "t", tmp_path)
    with pytest.raises(RuntimeError):
        niq.check_vizier(b"Name\n \n----\n", "t", tmp_path)


def test_nonpair_elsewhere_vetoes_a_rejection_and_sqls_z_pair_is_used():
    lemon = _lemon([("S", "150.0", "1.0", "1.5", "2.2", "QSO + star")])
    sqls = _sqls(
        [
            ("J100000.01+010000.0", "1.5", "2.2", "No lens object"),  # same system as S
            ("J110000.00+010000.0", "1.600", "2.4", "QSO pair (z=1.686, 1.600)"),
        ]
    )
    s = niq.build_sample({"J/MNRAS/520/3305/table1": lemon, "J/AJ/143/119/table4": sqls})
    assert list(s["name"]) == ["J110000.00+010000.0"]
    assert s.meta["dropped"]["non-pair classifications (incl. vetoed rejections)"] == 1
    assert s["z_source"][0] == pytest.approx(1.686) and s["z2"][0] == pytest.approx(1.6)


def test_check_vizier_ignores_column_descriptions(tmp_path):
    ok = b"#Column\te_RA\t(F5.2)\tMean error on RA [ucd=stat.error]\nName\n \n----\nJ1\n"
    niq.check_vizier(ok, "t", tmp_path)


def test_sqls_pair_format_redshifts():
    t = Table(
        rows=[
            ("J004757.25+144741.9", " ", "1.612", "", ""),
            ("J004757.87+144744.7", " ", "2.790", "9.42", "QSO pair"),
            ("J074013.44+292648.4", " ", "0.980", "", ""),
            ("J074013.42+292645.8", "(", "0.978", "2.64", "QSO pair"),
        ],
        names=("SDSS", "f_z", "z", "theta", "Com"),
        dtype=[str] * 5,
    )
    assert niq.sqls_redshifts(t, t[1]) == (1.612, 2.79)
    z, z2 = niq.sqls_redshifts(t, t[3])
    assert z == 0.98 and np.isnan(z2)  # flagged companion z is not used
    assert niq.comment_z_pair("QSO pair (z=1.686, 1.600.)") == (1.686, 1.6)


def test_read_tsv_refuses_a_second_resource(tmp_path):
    f = tmp_path / "t.tsv"
    f.write_text("Name\tSep\n \tarcsec\n----\t---\nJ1\t2.1\t\nName\tSep\n")
    with pytest.raises(ValueError):
        niq.read_tsv(f)


def test_sqls_nonpair_companion_does_not_veto_and_second_companion_keeps_primary_z():
    t = Table(
        rows=[
            ("J120000.00+010000.0", " ", "1.500", "", ""),
            ("J120000.10+010000.0", " ", "", "1.50", "QSO pair"),
            ("J120000.00+010002.5", " ", "", "2.50", "QSO+star"),
        ],
        names=("SDSS", "f_z", "z", "theta", "Com"),
        dtype=[str] * 5,
    )
    assert niq.sqls_redshifts(t, t[2])[0] == 1.5
    s = niq.build_sample({"J/AJ/143/119/table4": t})
    assert list(s["group"]) == ["rejected"]


def test_coords_refuses_mixed_rows():
    with pytest.raises(ValueError):
        niq._coords(["150.0003"], ["+01:00:00"])
    ra, de = niq._coords(["12:00:00"], ["+05"])  # whole-degree Dec with a sexagesimal RA
    assert ra[0] == pytest.approx(180.0) and de[0] == pytest.approx(5.0)


def test_sqls_newer_nonpair_row_does_not_veto_older_rejection():
    dr7 = _sqls([("J100000.00+010000.0", "1.5", "2.2", "QSO+star")])
    dr5 = Table(
        rows=[("J100000.01+010000.0", "1.5", "2.2", "No lens object")],
        names=("SDSS", "z", "theta", "Com"),
        dtype=[str] * 4,
    )
    s = niq.build_sample({"J/AJ/143/119/table4": dr7, "J/AJ/140/403/table3": dr5})
    assert list(s["group"]) == ["rejected"] and s["catalogue"][0] == "SQLS-DR5"
