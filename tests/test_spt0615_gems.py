"""Rules of scripts/spt0615_gems.py (frozen plan results/lens_consistency/spt0615_gems/PLAN.md)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import spt0615_gems as g  # noqa: E402


def test_sign_changes_interpolates():
    t = np.linspace(0, 1, 5)
    c = g.sign_changes(np.array([1.0, 0.5, -0.5, -1.0, -2.0]), t)
    assert np.allclose(c, [0.375])
    assert len(g.sign_changes(np.ones(5), t)) == 0
    # a sample exactly on zero still yields one crossing
    assert len(g.sign_changes(np.array([1.0, 0.0, -1.0, -1.0, -1.0]), t)) == 1


def test_classify_crossing():
    assert g.classify_crossing(np.array([0.5]), (0.4, 0.6)) == "through"
    assert g.classify_crossing(np.array([0.2]), (0.4, 0.6)) == "between"
    assert g.classify_crossing(np.array([]), (0.4, 0.6)) == "miss"
    assert g.classify_crossing(np.array([0.1, 0.9]), (0.4, 0.6)) == "miss"
    assert g.classify_crossing(np.array([0.5]), (0.6, 0.4)) == "through"  # order-free
    # three crossings, one inside: odd inside wins
    assert g.classify_crossing(np.array([0.1, 0.5, 0.9]), (0.4, 0.6)) == "through"


def test_parity_split():
    assert g.parity_split(50.0, -40.0)
    assert not g.parity_split(5.0, 4.0)


def test_classify_counter():
    assert g.classify_counter(1.0) == "consistent"
    assert g.classify_counter(2.2) == "consistent"
    assert g.classify_counter(3.0) == "offset"
    assert g.classify_counter(5.0) == "offset"
    assert g.classify_counter(5.1) == "none"
    assert g.classify_counter(None) == "none"


def test_joint_flag():
    miss = {"a": {"10.0": "miss", "10.2": "miss"}, "b": {"10.0": "miss", "10.2": "miss"}}
    ok = {"a": {"10.0": "consistent"}, "b": {"10.0": "consistent"}}
    assert g.joint_flag(miss, ok)
    mixed = {"a": {"10.0": "miss", "10.2": "through"}, "b": {"10.0": "miss", "10.2": "miss"}}
    assert not g.joint_flag(mixed, ok)
    none = {"a": {"10.0": "none"}, "b": {"10.0": "none"}}
    assert g.joint_flag(mixed, none)
    assert not g.joint_flag({}, {})


def test_project():
    assert np.isclose(g.project(np.array(1.0), np.array(0.0), (0.0, 0.0), (2.0, 0.0)), 0.5)
