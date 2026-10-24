"""Rule functions S1-S4 of results/lens_consistency/macs0647_specz/PLAN.md."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import macs0647_specz as m  # noqa: E402


def test_image_z_grade3_only():
    rows = [{"grade": "2", "z_best": "3.1"}, {"grade": "3", "z_best": "2.0"}]
    assert m.image_z(rows) == ("spec_z", 2.0)
    assert m.image_z([{"grade": "2", "z_best": "3.1"}]) == ("none", None)
    assert m.image_z([]) == ("none", None)


def test_image_z_ambiguous():
    rows = [{"grade": "3", "z_best": "2.0"}, {"grade": "3", "z_best": "3.0"}]
    assert m.image_z(rows) == ("ambiguous", None)
    close = [{"grade": "3", "z_best": "2.000"}, {"grade": "3", "z_best": "2.02"}]
    assert m.image_z(close)[0] == "spec_z"


def test_agree_uses_full_spread():
    # 2.00 vs 2.05: spread 0.0165 (1 + z) > 0.01, though each is 0.0083 from the median
    rows = [{"grade": "3", "z_best": "2.00"}, {"grade": "3", "z_best": "2.05"}]
    assert m.image_z(rows) == ("ambiguous", None)
    assert m.family_z({"a": 2.00, "b": 2.05}) == ("inconsistent", None)


def test_match_rows_and_shared_sources(tmp_path):
    from astropy.table import Table

    images = Table(
        {"image_id": ["a", "b", "c"], "ra": [10.0, 10.0, 11.0], "dec": [0.0, 0.0001, 0.0]}
    )
    csv_path = tmp_path / "t.csv"
    csv_path.write_text(
        "file,srcid,ra,dec,grating,root,sn50,zline,z_best,z_phot,grade\n"
        "f1,1,10.00002,0.00005,P,r,1,2,2.0,2,3\n"
        "f2,2,bad,0.0,P,r,1,2,2.0,2,3\n"
        "f3,3,11.001,0.0,P,r,1,2,2.0,2,3\n"
    )
    matches, n = m.match_rows(csv_path, images)
    assert n == 3
    assert set(matches) == {"a", "b"}  # f1 is 0.19 / 0.20 arcsec from both; f3 is 3.6 arcsec off
    assert m.shared_sources(matches, {"a", "b"}) == {"1": ["a", "b"]}
    assert m.shared_sources(matches, {"a"}) == {}


def test_family_z():
    assert m.family_z({}) == ("none", None)
    assert m.family_z({"1.1": 2.0, "1.2": 2.01})[0] == "spec_z"
    assert m.family_z({"1.1": 2.0, "1.2": 2.5}) == ("inconsistent", None)


def test_gate_s3():
    assert m.gate_s3(np.array([0.5, 1.0, 2.0, 2.9]))["passes"]
    assert not m.gate_s3(np.array([0.5, 1.0, 2.0]))["passes"]  # < 4 images
    assert not m.gate_s3(np.array([0.5, 1.0, np.nan, 3.5, 0.1]))["passes"]  # 60 %
    assert m.gate_s3(np.array([0.5, 1.0, 0.2, 3.5, 0.1]))["passes"]  # 80 %


def test_readings():
    d = np.array([1.0, 3.5, np.nan])
    assert m.readings(d, True) == ["no_flag", "flag", "flag"]
    assert m.readings(d, False) == ["not_decided"] * 3
