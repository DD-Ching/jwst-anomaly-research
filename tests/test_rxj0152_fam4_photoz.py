"""Frozen rule functions of scripts/rxj0152_fam4_photoz.py (PLAN.md steps 2-4 and the outcome)."""

import importlib.util
from pathlib import Path

import pytest
from astropy.table import Table

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "rxj0152_fam4_photoz.py"
_spec = importlib.util.spec_from_file_location("rxj0152_fam4_photoz", _PATH)
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)


def _row(zbmin, zbmax, odds):
    return {"zbmin": zbmin, "zbmax": zbmax, "odds": odds}


def test_classify():
    assert m.classify(None) == "undetermined"
    assert m.classify(_row(1.5, 1.8, 0.2)) == "consistent"
    assert m.classify(_row(1.69, 1.69, 0.99)) == "consistent"
    assert m.classify(_row(0.2, 0.5, 0.95)) == "excluded"
    assert m.classify(_row(0.2, 0.5, 0.9)) == "excluded"
    assert m.classify(_row(0.2, 0.5, 0.89)) == "undetermined"


def test_match_radius_and_stars():
    dec0 = -13.95
    cat = Table(
        {
            "id": [1, 2, 3],
            "RA": [28.18, 28.18, 28.18],
            "Dec": [dec0 + 0.1 / 3600, dec0 + 0.5 / 3600, dec0 + 0.7 / 3600],
            "stel": [0.9, 0.1, 0.1],
            "zb": [1.0, 2.0, 3.0],
            "zbmin": [0.5, 1.5, 2.5],
            "zbmax": [1.5, 2.5, 3.5],
            "odds": [0.5, 0.95, 0.95],
        }
    )
    got = m.match(28.18, dec0, cat)
    assert got["relics_id"] == 2  # the nearer row is a star
    assert got["sep_arcsec"] == pytest.approx(0.5, abs=1e-3)
    assert m.match(28.18, dec0, cat[[0, 2]]) is None  # star only, or 0.7" > 0.6"


def test_gate_p_and_outcome():
    ok = {"zbmin": 1.0, "zbmax": 2.0, "odds": 1.0}
    rows = (
        [{"z": 1.5, "match": ok}] * 5 + [{"z": 3.0, "match": ok}] * 3 + [{"z": 1.5, "match": None}]
    )
    g = m.gate_p(rows)
    assert (g["n_matched"], g["n_contain"], g["pass"]) == (8, 5, True)
    g = m.gate_p(rows[:7])
    assert not g["pass"]  # 7 matched < 8
    g = m.gate_p([{"z": 1.5, "match": ok}] * 4 + [{"z": 3.0, "match": ok}] * 4)
    assert g["frac"] == 0.5 and not g["pass"]
    assert m.outcome(False, ["consistent"]) == "U"
    assert m.outcome(True, ["excluded", "consistent", "undetermined"]) == "C"
    assert m.outcome(True, ["excluded"] * 3) == "E"
    assert m.outcome(True, ["excluded", "excluded", "undetermined"]) == "U"
