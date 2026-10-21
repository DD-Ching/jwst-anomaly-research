"""Offline tests for scripts/s1063_d035.py (frozen-rule helpers; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import s1063_d035 as d  # noqa: E402


def test_whole_list_gate_rms_and_unsolved():
    g = d.whole_list_gate([0.3, 0.4, 0.5], 0.48, 1.5)
    assert g["passes"] and np.isclose(g["rms_arcsec"], np.sqrt((0.09 + 0.16 + 0.25) / 3))
    assert not d.whole_list_gate([0.3, 2.0], 0.48, 1.5)["passes"]  # rms 1.43 > 0.72
    g = d.whole_list_gate([0.1, np.nan], 0.48, 1.5)
    assert not g["passes"] and g["n_solved"] == 1  # an unsolved image fails R1
    assert not d.whole_list_gate([], 0.48, 1.5)["passes"]


def test_outcome_order():
    assert d.outcome({"passes": True}, ["1"]) == "open: all z_spec systems"
    assert d.outcome({"passes": False}, ["1", "5"]) == "open: z_spec systems 1, 5"
    assert d.outcome({"passes": False}, []) == "map only"


def test_chan_images_drops_dup_rows():
    t = d.chan_images()
    assert not any(i in ("26.1", "27.1", "27.2", "27.3") for i in t["image_id"])
    assert len(t) == 117 and set(t["system"]) >= {"1", "2"}
    assert np.all(t["a"] == 1)
