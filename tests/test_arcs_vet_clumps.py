"""Tests for scripts/arcs_vet_clumps.py (closed-form values; no network)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("arcs_vet_clumps", _DIR / "arcs_vet_clumps.py")
vc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vc)


def test_clump_anti_fraction_closed_form():
    # |tan φ'| = |tan φ| (1 − g) / (1 + g) for |g| < 1, so anti ⇔ |φ| ≥ atan(√3 (1 + g) / (1 − g))
    for g in (0.0, 0.2, 0.44):
        expect = 1 - np.degrees(np.arctan(np.sqrt(3) * (1 + g) / (1 - g))) / 90
        assert abs(vc.clump_anti_fraction(g) - expect) < 1e-3
    assert abs(vc.clump_anti_fraction(0.0) - 1 / 3) < 1e-3
