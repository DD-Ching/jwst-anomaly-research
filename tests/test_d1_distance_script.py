"""Offline tests for scripts/d1_distance.py FRB helpers (ISM-model comparison, YMW16 lookup)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("d1_distance", _DIR / "d1_distance.py")
d1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(d1)


def _table(model: str, ism: list[float], z_low: list[float], flag_low: list[bool]) -> Table:
    n = len(ism)
    return Table(
        {
            "frb": [f"F{i}" for i in range(n)],
            "z": np.full(n, 0.1),
            "dm_obs": np.full(n, 300.0),
            f"dm_ism_{model}": ism,
            "z_low": z_low,
            "z_high": np.zeros(n),
            "flag_low": flag_low,
            "flag_high": np.zeros(n, bool),
        }
    )


def test_compare_ism_flags_either_model(tmp_path, monkeypatch):
    monkeypatch.setattr(d1, "OUT", tmp_path)
    ne = _table("ne2001", [100.0, 50.0], [1.0, 6.0], [False, True])
    ym = _table("ymw16", [150.0, 40.0], [2.5, 4.0], [False, False])
    s = d1.compare_ism(ym, ne, "ymw16", 5.8)
    assert s["n_joined"] == 2
    assert s["flag_either"] == ["F1"]  # flagged under NE2001 only still counts
    assert s["ism_ratio_range"] == [0.8, 1.5]
    assert s["max_abs_pull_change"] == ["F1", 2.0]
    assert (tmp_path / "frb_ism_compare.ecsv").exists()


def test_ymw16_known_sightline():
    pytest.importorskip("ymw16")
    # FRB 20121102A (l = 174.95, b = -0.22): YMW16 gives ~287 pc cm^-3 through the whole disc
    dm = float(d1.ymw16_dm_ism(82.99458, 33.14793)[0])
    assert 250.0 < dm < 320.0


def test_td_lens_uses_a_final_j1206_chain_as_given():
    td = {
        n: {"ddt_model": np.full(4, 5000.0), "kappa": np.full(4, 0.1), "n_chain": 4}
        for n in ("J1206", "RXJ1131")
    }
    plain = d1._td_lens(td, "kext")
    assert plain["J1206"]["ddt"] == pytest.approx(np.full(4, 5000.0 / 0.9))
    final = np.array([5700.0, 5710.0, 5720.0, 5730.0])
    swapped = d1._td_lens(td, "kext", final)
    assert swapped["J1206"]["ddt"] is final  # kappa_ext already folded in: not applied twice
    assert swapped["RXJ1131"]["ddt"] == pytest.approx(plain["RXJ1131"]["ddt"])
    assert set(d1.J1206_FINAL) == {"final_power_law", "final_composite"}
