"""Offline tests for scripts/s1_energy.py (sign convention, quality rule, threshold matching)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import s1_energy as se  # noqa: E402


def test_zval_convention_physical_vs_program_stretch():
    # b is a's source with t -> s_phys t, E -> E / s_phys; the program stretches b by 1 / s_phys
    lna, s_phys = math.log(300.0), 1.6
    lnb = lna - math.log(s_phys)
    assert abs(se.zval(lna, 0.1, lnb, 0.1, 1.0 / s_phys, 0.02)) < 1e-12
    # using the physical factor instead of the program's would be wrong by 2 ln s_phys
    assert abs(se.zval(lna, 0.1, lnb, 0.1, s_phys, 0.02)) > 2.0


def test_quality_rule_treats_bad_values_as_missing():
    lne = np.log(np.array([200.0, 5.0, 2e4, 300.0, np.nan]))
    sig = np.array([0.1, 0.1, 0.1, 0.9, 0.1])
    assert se.quality(lne, sig).tolist() == [True, False, False, False, False]


def test_thr_for_matches_count():
    m2 = np.array([0.99, 0.98, 0.97, 0.96, 0.95, 0.94])
    zp = np.array([False, True, True, False, True, True])
    thr = se.thr_for(m2, zp, 2)
    assert ((m2 > thr) & zp).sum() == 2
    assert thr == 0.95
