import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import radial_axial as rax  # noqa: E402


def test_n_axes_merges_azimuths_mod_180():
    assert rax.n_axes([10.0, 190.0, 15.0]) == 1
    assert rax.n_axes([0.0, 60.0, 120.0]) == 3
    assert rax.n_axes([175.0, 2.0]) == 1  # wraps at 180
    assert rax.n_axes([]) == 0


def test_on_line_counts_collinear_arcs():
    x = np.array([0.0, 1.0, 2.0, 3.0, 0.0])
    y = np.array([0.0, 1.0, 2.1, 3.0, 3.0])
    assert rax.on_line(x, y) == 4
    assert rax.on_line(x[:1], y[:1]) == 1


def test_gate_k_uses_the_strongest_centre():
    peaks = [{"n_lines": 6, "axial_a": 0.73}, {"n_lines": 3, "axial_a": 0.99}]
    assert not rax.gate_k(peaks)
    peaks[0]["axial_a"] = 0.95
    assert rax.gate_k(peaks)
