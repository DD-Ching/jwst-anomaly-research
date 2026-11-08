"""Frozen decision rule of scripts/macs0308_ares.py (offline)."""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location(
    "macs0308_ares", ROOT / "scripts" / "macs0308_ares.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


@pytest.mark.parametrize(
    ("path_ok", "sep", "delays", "expected"),
    [
        (False, 0.0, [70.0, 75.0], "stopped: path"),
        (True, 0.5, [70.0, 75.0], "stopped: no SN image"),
        (True, float("nan"), [70.0, 75.0], "stopped: no SN image"),
        (True, 0.1, [70.0, 75.0], "agrees"),
        (True, 0.1, [10.0, 75.0], "agrees"),
        (True, 0.1, [9.9, 75.0], "undetermined"),
        (True, 0.1, [75.0], "undetermined"),
        (True, 0.1, [-1.0, 75.0], "conflicts"),
        (True, 0.1, [-0.99, 75.0], "undetermined"),
        (True, 0.3, [-30.0], "conflicts"),
    ],
)
def test_decide(path_ok, sep, delays, expected):
    assert mod.decide(path_ok, sep, delays) == expected
