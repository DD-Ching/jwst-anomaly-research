"""Unit tests for scripts/cal_dither_phot.py on synthetic per-exposure tables (offline)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("cal_dither_phot", _DIR / "cal_dither_phot.py")
cdp = importlib.util.module_from_spec(_spec)
sys.modules["cal_dither_phot"] = cdp
_spec.loader.exec_module(cdp)


def _table(target_gain: float, seed: int = 1) -> tuple[Table, list[str]]:
    """Two epochs of 6 exposures; a common per-exposure factor; target brighter in o212."""
    rng = np.random.default_rng(seed)
    knots = [f"k{i:02d}" for i in range(8)]
    base = {u: 100.0 + 20 * i for i, u in enumerate(["tgt", *knots])}
    rows = []
    for epoch in ("o211", "o212"):
        for j in range(6):
            name = f"jw01176{epoch[1:]}001_0210{j}_00001_nrcb1_cal.fits"
            common = rng.uniform(0.9, 1.1)
            for u, f in base.items():
                gain = 1 + target_gain if (u == "tgt" and epoch == "o212") else 1.0
                rows.append(
                    {
                        "uid": u,
                        "file": name,
                        "epoch": cdp.epoch_of(name),
                        "f0": f * common * gain * (1 + rng.normal(0, 0.003)),
                    }
                )
    return Table(rows=rows), knots


def test_epoch_of():
    assert cdp.epoch_of("jw01176211001_02101_00001_nrcb1_cal.fits") == "o211"
    assert cdp.epoch_of("s3://b/jw01208004001_x_cal.fits") == "o004"
    assert cdp.epoch_of("other.fits") == "?"


def test_target_excess_recovered_and_knots_null():
    tab, knots = _table(0.05)
    t = cdp.normalise(tab, ["tgt", *knots], knots, "f0")
    fe = dict(zip(tab["file"], tab["epoch"], strict=True))
    ep = np.array([fe[f] for f in t["_files"]])
    ex, err = cdp.epoch_excess(t["tgt"], ep, "o211", "o212")
    assert abs(ex - 0.05) < 0.005 and 0 < err < 0.01
    kex = np.array([cdp.epoch_excess(t[k], ep, "o211", "o212")[0] for k in knots])
    assert np.all(np.abs(kex) < 0.01)
    assert cdp.robust_z(ex, kex)[2] > 5


def test_no_excess_gives_small_z():
    tab, knots = _table(0.0, seed=3)
    t = cdp.normalise(tab, ["tgt", *knots], knots, "f0")
    fe = dict(zip(tab["file"], tab["epoch"], strict=True))
    ep = np.array([fe[f] for f in t["_files"]])
    ex, _ = cdp.epoch_excess(t["tgt"], ep, "o211", "o212")
    assert abs(ex) < 0.01


def test_load_knots_range_exclusion_and_dedup(tmp_path):
    d = {
        "bands": {
            "F150W": {
                "knots": [
                    {"ra": 64.0, "dec": -24.0, "e1_njy": 200.0},  # kept
                    {"ra": 64.01, "dec": -24.0, "e1_njy": 100.0},  # too faint
                    {"ra": 64.02, "dec": -24.0, "e1_njy": 300.0},  # next to the target
                ]
            },
            "F200W": {"knots": [{"ra": 64.0000100, "dec": -24.0, "e1_njy": 220.0}]},  # duplicate
        }
    }
    p = tmp_path / "k.json"
    p.write_text(json.dumps(d))
    k = cdp.load_knots(p, {"t": (64.02, -24.00001)}, ["F150W", "F200W"])
    assert list(k.values()) == [(64.0, -24.0)]
