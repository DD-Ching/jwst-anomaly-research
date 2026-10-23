"""Frozen rule of scripts/a2744_vet78_cats.py (PLAN_vet78.md)."""

import sys
from pathlib import Path

from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_vet78_cats as v  # noqa: E402


def _knots(own=(0.5, 0.5), t_cat=2.0, t_unc=0.3, t_ra=1.0):
    k = {f"78.{i + 1}": {"c_ra": 1.0, "c_to_own_catalogued_arcsec": d} for i, d in enumerate(own)}
    k["78.3"] = {"c_ra": t_ra, "c_to_catalogued_arcsec": t_cat, "c_to_uncover_78_3c_arcsec": t_unc}
    return k


def test_decide_branches():
    assert v.decide(_knots(own=(0.5, 1.6))) == "tool not validated"
    assert v.decide(_knots(t_ra=None)) == "tool not validated"
    assert v.decide(_knots(t_cat=0.9)).startswith("explained")
    assert v.decide(_knots(t_cat=2.0, t_unc=1.2)).startswith("untestable")
    assert v.decide(_knots(t_cat=2.0, t_unc=0.3)) == "flag survives test 1"


def test_decide_missing_c_of_control_fails_tool():
    k = _knots()
    k["78.1"]["c_ra"] = None
    assert v.decide(k) == "tool not validated"


def test_c_type_is_farthest_from_ab_midpoint():
    ab = Table({"ra": [10.0, 10.0002], "dec": [0.0, 0.0]})
    pred = Table({"ra": [10.0, 10.0002, 10.002], "dec": [0.0, 0.0, 0.001]})
    assert v.c_type(pred, ab) == 2
    assert v.c_type(pred[:2], ab) is None


def _knots_1b(own=(1.2, 1.1), t_cat=1.5):
    k = {f"78.{i + 1}": {"c_ra": 1.0, "c_to_own_catalogued_arcsec": d} for i, d in enumerate(own)}
    k["78.3"] = {"c_ra": 1.0, "c_to_catalogued_arcsec": t_cat}
    return k


def test_decide_1b_branches():
    k = _knots_1b()
    k["78.3"]["c_ra"] = None
    assert v.decide_1b(k) == "tool not validated"
    assert v.decide_1b(_knots_1b(own=(1.2, 1.1))).startswith("explained: this model's c positions")
    assert v.decide_1b(_knots_1b(own=(0.8, 1.1), t_cat=0.9)).startswith("explained: model position")
    assert v.decide_1b(_knots_1b(own=(0.8, 1.1), t_cat=1.5)) == "flag survives test 1b"
