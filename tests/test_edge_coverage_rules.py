"""Frozen rules of results/dimming/edge_coverage/PLAN.md (edge.edge_pair, gate_v, decision)."""

import importlib.util
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

_P = Path(__file__).resolve().parents[1] / "results/dimming/edge_coverage/edge.py"
_spec = importlib.util.spec_from_file_location("edge_coverage_edge", _P)
ed = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ed)


def test_edge_pair_needs_coverage_and_low_weight():
    cov = np.array([True, True, True, False, True])
    rel = np.array([0.49, 0.5, np.nan, 0.0, 1.2])
    assert ed.edge_pair(cov, rel).tolist() == [True, False, True, False, False]


def test_gate_v():
    assert ed.gate_v(9) and ed.gate_v(11)
    assert not ed.gate_v(8)


def test_decision_any_t_e_over_ten_percent():
    u = np.array([10.0, 20.0, 30.0])
    assert ed.decision(np.array([9.5, 19.0, 27.5]), u) == "gain_stands"
    assert ed.decision(np.array([10.0, 17.9, 30.0]), u) == "over_counts"


def test_wht_rel_takes_the_3x3_minimum_and_zero_off_grid():
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [10.0, 0.0]
    w.wcs.crpix = [1.0, 1.0]
    w.wcs.cdelt = [-1e-5, 1e-5]
    grid = np.ones((10, 10), np.float32)
    grid[5, 6] = np.nan  # neighbour of node (5, 5)
    step = 4
    ra, dec = w.all_pix2world([5 * step, 2 * step, -100 * step], [5 * step, 2 * step, 2 * step], 0)
    rel = ed.wht_rel(grid, w, ra, dec, step)
    assert rel[0] == 0.0  # NaN neighbour counts as 0
    assert rel[1] == 1.0
    assert rel[2] == 0.0  # off the grid
