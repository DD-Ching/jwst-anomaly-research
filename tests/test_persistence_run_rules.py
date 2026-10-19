"""Frozen decision rules of results/dimming/persistence_run/PLAN.md (summarise.py)."""

import importlib.util
from pathlib import Path

from astropy.table import Table

_P = Path(__file__).resolve().parents[1] / "results/dimming/persistence_run/summarise.py"
_spec = importlib.util.spec_from_file_location("persistence_run_summarise", _P)
sm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sm)


def _mon(pairs):
    return Table(rows=pairs, names=["mosaic_snr", "verdict"])


def test_efficiency_bins_and_rule():
    rows = [(25.0, "on_sky")] * 9 + [(25.0, "inconclusive")] + [(5.0, "undetected")] * 4
    eff = sm.efficiency(_mon(rows))
    assert [b["n_pairs"] for b in eff] == [4, 0, 10, 0]
    assert eff[2]["on_sky_fraction"] == 0.9
    assert sm.usable(eff, 0.0)
    assert not sm.usable(eff, 0.03)  # false-afterimage rate above 2 %


def test_rule_needs_a_populated_high_snr_bin():
    eff = sm.efficiency(_mon([(25.0, "on_sky")] * 9))  # 9 < MIN_PAIRS
    assert not sm.usable(eff, 0.0)
    eff = sm.efficiency(_mon([(60.0, "on_sky")] * 7 + [(60.0, "persistence")] * 3))
    assert not sm.usable(eff, 0.0)  # 0.7 < 0.8
