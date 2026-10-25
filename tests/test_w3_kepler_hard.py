"""Offline tests of the hard-control generators in scripts/w3_kepler_hard.py."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w3_kepler as wk  # noqa: E402
import w3_kepler_hard as wh  # noqa: E402


def _star(n_q=3, n=500):
    t = np.concatenate([q * (n * wk.CAD + 0.5) + np.arange(n) * wk.CAD for q in range(n_q)])
    return {
        "TIME": t,
        "CADENCENO": np.arange(len(t), dtype=float),
        "SAP_FLUX": np.full(len(t), 1000.0),
        "SAP_FLUX_ERR": np.full(len(t), 0.1),
        "SAP_QUALITY": np.zeros(len(t)),
        "QUARTER": np.repeat(np.arange(1, n_q + 1, dtype=float), n),
    }


def test_single_point_adds_one_outlier_per_side():
    st = _star()
    good = wk.good_mask(st)
    for j in range(20):
        f, _ = wh.hard_flux("single_point", st, good, np.random.default_rng([1, j]))
        assert 1 <= np.sum(np.isclose(f, 1.2)) <= 2 and (f < 0.4).sum() >= 1


def test_quarter_edge_dip_touches_a_boundary():
    st = _star()
    good = wk.good_mask(st)
    for j in range(20):
        f, _ = wh.hard_flux("quarter_edge", st, good, np.random.default_rng([2, j]))
        dim = np.flatnonzero(f < 0.5)
        q = st["QUARTER"]
        ends = [np.flatnonzero(q == k).max() for k in (1, 2)]
        starts = [np.flatnonzero(q == k).min() for k in (2, 3)]
        near = min(min(abs(dim.max() - e) for e in ends), min(abs(dim.min() - s) for s in starts))
        spans = any(dim.min() <= e < dim.max() for e in ends)
        assert len(dim) and (near <= 6 or spans)


def test_synth_bright_eclipse_has_plateau_and_dip():
    st = _star()
    good = wk.good_mask(st)
    hits = 0
    for j in range(10):
        f, _ = wh.hard_flux("synth_bright_eclipse", st, good, np.random.default_rng([3, j]))
        assert f.max() <= 1.2001
        hits += f.max() >= 1.099 and f.min() < 0.3
    assert hits >= 3


def test_does_not_mutate_pilot_settings_on_import():
    assert wk.OUT.name == "w3_kepler" and wk.LEDGER_CAP == 2_000_000_000


def _fams(passed):
    return {"a": {"passed": passed, "n": 100}}


def _eb(windows=100, win_pass=0, stars=6, star_pass=0):
    return {
        "flagged_windows": windows,
        "windows_passing_single_window": win_pass,
        "stars_with_flagged_window": stars,
        "stars_with_single_window_pass": star_pass,
    }


def test_hard_gate_rules():
    assert wh.hard_gate(_fams(10), _eb(), 4)["pass"]
    assert not wh.hard_gate(_fams(11), _eb(), 4)["pass"]
    assert not wh.hard_gate(_fams(0), _eb(win_pass=11), 4)["pass"]
    assert not wh.hard_gate(_fams(0), _eb(stars=4), 4)["pass"]  # fewer than 5 EB stars
    assert not wh.hard_gate(_fams(0), _eb(), 2)["pass"]  # fewer than 3 brightening sources
    assert not wh.hard_gate(_fams(0), _eb(windows=0, stars=0), 4)["pass"]
