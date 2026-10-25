"""Frozen rule functions of scripts/rxj0152_fam4_photometry.py (PLAN.md rules 5-8 and outcomes)."""

import importlib.util
from pathlib import Path

import pytest

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "rxj0152_fam4_photometry.py"
_spec = importlib.util.spec_from_file_location("rxj0152_fam4_photometry", _PATH)
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

SED = {"f115w": 1.0, "f150w": 2.0, "f277w": 4.0, "f444w": 5.0}


def test_flux_window():
    lo, hi = m.flux_window(453.0)
    assert lo == pytest.approx(453.0 / (3 * 15.10))
    assert hi == pytest.approx(3 * 453.0 / 10.05)


def test_colour_match_scale_free_and_tolerance():
    assert m.colour_match({b: v / 13 for b, v in SED.items()}, SED)
    redder = dict(SED, f444w=5.0 * 10**0.31)
    assert not m.colour_match(redder, SED)
    assert m.colour_match(dict(SED, f444w=5.0 * 10**0.29), SED)
    assert not m.colour_match(dict(SED, f115w=0.0), SED)
    assert not m.colour_match(dict(SED, f115w=float("nan")), SED)


def test_outcomes():
    assert m.outcome(True, True, True, 1) == "O1_counterpart"
    assert m.outcome(True, True, False, 1) == "O1_counterpart"
    assert m.outcome(True, True, True, 0) == "O2_missing"
    assert m.outcome(True, True, False, 0) == "O3_undetermined"
    assert m.outcome(False, True, True, 1) == "O3_undetermined"
    assert m.outcome(True, False, True, 0) == "O3_undetermined"


def test_known_case_gate():
    assert m.known_case_pass(6)
    assert not m.known_case_pass(5)


def test_region():
    assert m.in_region(*m.PREDICTED["glafic"])
    assert not m.in_region(28.183443, -13.956815)  # image 4.1
