"""Offline tests for scripts/wh_support_scales.py (arithmetic from published formulas)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import wh_support_scales as ws  # noqa: E402


def test_mmp_binding_energy_identity():
    # (2.3) and (5.31) evaluated separately give |E_min| r_e = g^2 hbar c / (256 pi) for any q
    for q, g in ((1e6, 0.3), (3e12, 0.36)):
        r_e, e_ev = ws.mmp_binding_from_q(q, g)
        assert math.isclose(
            e_ev * r_e / ws.HBAR_C_EV_M, ws.mmp_binding_times_radius(g), rel_tol=1e-12
        )


def test_mm_reproduces_published_example():
    # Maldacena-Milekhin (3.27): c2 ~ 7e72, l ~ 3e3 ly, E_bin ~ -5e9 kg for r_e = 0.05 s, R5 = 50 um
    p = ws.mm_parameters(1.5e7, 50e-6)
    assert 6e72 < p["c2"] < 8e72
    assert 2.5e3 < p["ell_ly"] < 3.5e3
    assert -6e9 < p["E_bin_kg"] < -4e9


def test_scales_entries_carry_sources():
    res = ws.scales()
    for entry in res["entries"].values():
        assert entry["source"] and entry["outputs"]
