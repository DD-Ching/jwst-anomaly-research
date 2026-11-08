"""Frozen decision rule of scripts/athena_hst18470.py (offline)."""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location(
    "athena_hst18470", ROOT / "scripts" / "athena_hst18470.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


@pytest.mark.parametrize(
    ("delta", "sigma", "seen", "expected"),
    [
        (250.0, 40.0, True, "appeared"),
        (200.0, 40.0, True, "appeared"),
        (200.0, 40.0, False, "undetermined"),
        (199.0, 40.0, True, "undetermined"),
        (120.0, 40.0, True, "undetermined"),
        (119.0, 40.0, False, "not seen"),
        (-80.0, 40.0, False, "not seen"),
        (100.0, 0.0, True, "undetermined"),
        (float("nan"), 40.0, True, "undetermined"),
    ],
)
def test_decide(delta, sigma, seen, expected):
    assert mod.decide(delta, sigma, seen) == expected


@pytest.mark.parametrize(
    ("ds", "expected"),
    [
        (["appeared", "not seen"], "appeared"),
        (["not seen", "not seen"], "not seen"),
        (["not seen", "undetermined"], "undetermined"),
        ([], "undetermined"),
    ],
)
def test_combine(ds, expected):
    assert mod.combine(ds) == expected


def test_evaluate_missing_visit():
    res = [{"visit": "a", "error": "not found"}, {"visit": "b", "error": "not found"}]
    assert mod.evaluate(res, mod.SEEN)["outcome"] == "undetermined"


def test_evaluate_not_seen():
    row = {name: {"flux_njy": 100.0, "err_njy": 30.0} for name in mod.POS}
    res = [{"visit": "a", "mjd": 1.0, **row}, {"visit": "b", "mjd": 2.0, **row}]
    out = mod.evaluate(res, mod.SEEN)
    assert out["outcome"] == "not seen"
    assert out["positions"]["athena2_12781"]["z"] == 0.0
