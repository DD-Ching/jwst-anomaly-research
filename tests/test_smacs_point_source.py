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
    assert ps.decide([], 30, None) == "undetermined (V1 could not be split)"


def test_held_out_by_group_and_spatial_halves():
    ra = np.arange(20, dtype=float)
    g = ps.spatial_groups(ra)
    assert g.sum() == 10 and g[-1] and not g[0]
    r = np.where(ra == 2.0, 9.0, 2.0 + 0.01 * (ra % 2))  # flat locus, one west outlier
    rec = ps.held_out_by_group(g, r)
    assert rec == [0.9, 1.0]  # east cut tested on west: 9 / 10
    assert ps.held_out_by_group(np.ones(5, bool), np.ones(5)) == []


def test_decide_v2_gates_and_reach():
    ok = [0.95, 0.92]
    assert ps.decide_v2(ok, ok, 19, 0.0, 12, 0.85) == "validated to ref_mag 27"
    assert ps.decide_v2(ok, ok, 19, 0.0, 12, 0.5) == "validated to ref_mag 24 only"
    assert ps.decide_v2(ok, ok, 19, 0.0, 9, 1.0) == "validated to ref_mag 24 only"  # V3 descriptive
    assert ps.decide_v2(ok, [0.95, 0.85], 19, 0.0, 12, 0.9) == "not validated (V1b fail)"
    assert ps.decide_v2([0.8, 0.9], ok, 5, 0.0, 12, 0.9) == "not validated (V1a, V2 fail)"
    assert ps.decide_v2([], ok, 19, 0.0, 12, 0.9) == "undetermined (V1 could not be split)"
