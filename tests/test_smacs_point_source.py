"""Offline tests of the SMACS point-source screen rule (scripts/smacs_point_source.py)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import smacs_point_source as ps  # noqa: E402


def test_locus_cut_median_plus_three_mad():
    r = np.array([2.0, 2.1, 1.9, 2.0, 2.2, 1.8, np.nan])
    med, s, cut = ps.locus_cut(r)
    assert med == 2.0
    assert np.isclose(s, 1.4826 * 0.1)
    assert np.isclose(cut, 2.0 + 3 * 1.4826 * 0.1)


def test_locus_mask_uses_fallback_window_then_undetermined():
    r = np.full(12, 2.0)
    mag = np.array([18.5] * 5 + [20.0] * 7)  # only 7 in [19, 24], 12 in [18, 24]
    m, w = ps.locus_mask(mag, r)
    assert w == [18.0, 24.0] and m.sum() == 12
    m, w = ps.locus_mask(np.full(12, 17.0), r)
    assert w is None and not m.any()


def test_held_out_recovery_both_directions():
    ids = np.arange(20)
    r = np.where(ids == 3, 9.0, 2.0 + 0.01 * ids)  # one odd outlier
    rec = ps.held_out_recovery(ids, r)
    assert len(rec) == 2
    assert rec[0] == 0.9  # train even, test odd: 9 / 10 recovered
    assert rec[1] == 1.0


def test_decide_gates():
    assert ps.decide([0.95, 0.92], 30, 0.05) == "validated"
    assert ps.decide([0.95, 0.85], 30, 0.05) == "not validated (V1 fails)"
    assert ps.decide([0.95, 0.95], 9, 0.0) == "not validated (V2 fails)"
    assert ps.decide([0.95, 0.95], 30, 0.2) == "not validated (V2 fails)"
    assert ps.decide([], 30, None) == "not validated (V1 fails) (V2 fails)"
