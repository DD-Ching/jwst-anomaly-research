"""Offline tests for scripts/m0417_images.py (frozen rules; no network)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import m0417_images as mi  # noqa: E402

OUT = ROOT / "results" / "lens_consistency" / "macs0417_images"


def _ip(ids, d):
    return Table({"image_id": ids, "system": [i.split(".")[0] for i in ids], "dtheta_arcsec": d})


def test_gate_k_rms_and_fraction():
    ok = mi.gate_k(_ip(["1.1", "1.2", "1.3", "1.4", "1.5"], [0.2, 0.3, 0.4, 0.5, 0.6]))
    assert ok["passes"] and ok["matched_fraction"] == 1.0
    # one 5" outlier pushes the rms over 1.0" although 4 / 5 are matched
    bad = mi.gate_k(_ip(["1.1", "1.2", "1.3", "1.4", "1.5"], [0.2, 0.3, 0.4, 0.5, 5.0]))
    assert not bad["passes"] and bad["matched_fraction"] == 0.8
    # an unsolved image counts as unmatched
    nan = mi.gate_k(_ip(["1.1", "1.2"], [0.2, np.nan]))
    assert nan["matched_fraction"] == 0.5 and not nan["passes"]


def test_best_z_prefers_full_match_over_lower_rms():
    d = [np.array([0.1, 3.5]), np.array([0.9, 1.0]), np.array([np.nan, np.nan])]
    assert mi.best_z(d) == 1  # rms of z0 is lower, but only z1 matches every image
    assert mi.best_z([np.array([4.0, 5.0]), np.array([3.5, 3.6])]) == 1  # none full: min rms


def test_readings_flag_only_when_every_passing_model_misses():
    images = Table({"image_id": ["1.1", "1.2"], "system": ["1", "1"], "z_spec": [True, False]})
    tables = {"a": _ip(["1.1", "1.2"], [4.0, 0.5]), "b": _ip(["1.1", "1.2"], [5.0, 6.0])}
    r = mi.readings(tables, ["a", "b"], images)
    assert list(r["reading"]) == ["lensing_residual_input", "no_flag"]
    assert list(mi.readings(tables, [], images)["reading"]) == ["no_flag", "no_flag"]


def test_committed_summary_matches_frozen_outcome():
    s = json.loads((OUT / "summary.json").read_text())
    for field in ("macs0417", "macs1423"):
        assert s[field]["passing_models"] == []
        assert s[field]["flags"] == []
    assert s["rules"]["K_MAX_RMS_ARCSEC"] == mi.K_MAX_RMS_ARCSEC
