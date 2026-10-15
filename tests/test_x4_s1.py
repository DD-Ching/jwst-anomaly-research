"""Offline tests for scripts/x4_s1_degeneracy.py (X4 S1). Waveform-dependent parts need lalsuite."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_P = Path(__file__).resolve().parents[1] / "scripts" / "x4_s1_degeneracy.py"
_spec = importlib.util.spec_from_file_location("x4_s1_degeneracy", _P)
x4 = importlib.util.module_from_spec(_spec)
sys.modules.setdefault("x4_s1_degeneracy", x4)
_spec.loader.exec_module(x4)


def test_chirp_mass_round_trip():
    for mtot, q in ((30.0, 1.5), (200.0, 4.0)):
        m1, m2 = x4.masses(x4.chirp(mtot, q), q)
        assert abs(m1 + m2 - mtot) < 1e-9 and abs(m1 / m2 - q) < 1e-9


def test_rho_req_uses_exact_residual():
    assert abs(x4.rho_req(0.0) - 3.0 / np.sqrt(1e-24)) / (3e12) < 1e-6
    mm = 0.01
    assert abs(x4.rho_req(mm) - 3.0 / np.sqrt(1 - 0.99**2)) < 1e-12


def test_frozen_case_list():
    c = x4.cases()
    assert len(c) == 36 and len(set(c)) == 36
    assert {m for m, _, _ in c} == {30.0, 60.0, 200.0}


def test_overlap_invariances():
    pytest.importorskip("lalsimulation")
    h = x4.waveform(x4.chirp(60.0, 1.5), 1.5, 0.2, 0.2)
    shifted = 3.7 * h * np.exp(1j * (0.4 + 2 * np.pi * x4._FREQS * 0.01))
    assert x4.overlap(h, shifted) > 0.999


def test_geometric_fallback_matches_exact_at_large_w():
    for s in (1, -1):
        for y in (1.5, 3.0, 5.0):
            e = x4.gl.amp_factor(150.0, y, s)
            g = x4.amp_geometric(150.0, y, s)
            assert abs(e - g) <= 2e-3 * max(abs(e), 1e-3)
