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
