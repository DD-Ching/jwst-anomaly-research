"""Rules 3-4 of results/lens_consistency/elgordo_third_specz/PLAN.md."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import elgordo_third_specz as m  # noqa: E402


def test_agrees_scales_with_one_plus_z():
    assert m.agrees(4.23, 4.1879)  # 0.008 (1 + z)
    assert not m.agrees(4.25, 4.1879)  # 0.012 (1 + z)


def test_gate_z_counts_and_fraction():
    ok = [(2.56, 2.5635), (2.83, 2.8254), (4.19, 4.1879)]
    assert m.gate_z(ok)["passes"]
    assert not m.gate_z(ok[:2])["passes"]  # fewer than 3 images
    bad = ok + [(1.0, 2.5635), (0.5, 2.8254)]  # 3 / 5 agree
    assert not m.gate_z(bad)["passes"]
    assert not m.gate_z([])["passes"]


def test_reading():
    assert m.reading("spec_z", 4.19, 4.1879, True) == "confirmed"
    assert m.reading("spec_z", 0.87, 4.1879, True) == "z_differs"
    assert m.reading("none", None, 4.1879, True) == "no_spectrum"
    assert m.reading("ambiguous", None, 4.1879, True) == "ambiguous"
    assert m.reading("spec_z", 4.19, 4.1879, False) == "not_decided"


def test_family_parses_ids():
    assert m._family("1.1a") == 1
    assert m._family("20c") == 20
    assert m._family("6a") == 6
