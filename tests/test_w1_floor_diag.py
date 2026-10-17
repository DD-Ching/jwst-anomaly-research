"""Tests for the W1 B-mode floor diagnostic helpers (scripts/w1_floor_diag.py)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w1_floor_diag", _DIR / "w1_floor_diag.py")
fd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fd)


def test_bright_neighbour_reach_and_dmag():
    # a bright source at the origin (sigma 1"), faint ones at 2", 4" and 6"; one bright one at 2"
    x = np.array([0.0, 2.0, 4.0, 6.0, 0.0])
    y = np.array([0.0, 0.0, 0.0, 0.0, 2.0])
    mag = np.array([20.0, 24.0, 24.0, 24.0, 21.0])
    a = np.array([1.0, 0.1, 0.1, 0.1, 0.1])
    out = fd.bright_neighbour(x, y, mag, a, dmag=2.0, k=5.0)
    assert out.tolist() == [False, True, True, False, False]
    assert not fd.bright_neighbour(x, y, mag, a, dmag=2.0, k=1.0).any()


def test_chain_neighbours_counts_collinear_aligned_only():
    # three segments on a north-south line with PA 0, plus one off-line and one misaligned
    x = np.array([0.0, 0.0, 0.0, 5.0, 0.0])
    y = np.array([0.0, 3.0, 6.0, 0.0, -3.0])
    pa = np.array([0.0, 2.0, 178.0, 0.0, 60.0])
    out = fd.chain_neighbours(x, y, pa, radius=10.0, pa_tol=10.0, line_tol=10.0)
    assert out.tolist()[:3] == [2, 2, 2]
    assert out[3] == 0 and out[4] == 0
