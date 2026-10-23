"""Frozen decision rules of the STScI PSF capability check (w2_j2244_stpsf/PLAN.md section 5).

Pure functions, tested in the project environment (tests/test_w2_j2244_stpsf_rules.py). Every
threshold is an ASSUMPTION named in Params.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Params:
    hookup_pass_px: float = 0.1  # case-1 sky-position rms (px at 0.128"): pass at or below
    hookup_stop_px: float = 0.3  # above this: stop (mapping problem)
    star_redchi_max: float = 2.0  # case 1: chi2/nu within r <= 5 px, every exposure
    star_z_max: float = 3.0  # case 1: |z| in every annulus, every exposure
    pair_z_max: float = 3.0  # case 2: |z_P| below this
    inj_z_min: float = 5.0  # case 3: z_P at or above this counts as recovered
    inj_r_min: float = 0.5  # case 3: recovered fraction at 20 mag at or above this


P = Params()


def hookup(rms_px: float, p: Params = P) -> str:
    if not math.isfinite(rms_px) or rms_px > p.hookup_stop_px:
        return "stop"
    return "pass" if rms_px <= p.hookup_pass_px else "pass (flagged)"


def star_ok(redchi: list[float], z_annuli: list[list[float]], p: Params = P) -> bool:
    """redchi: one value per exposure; z_annuli: per exposure, one z per annulus."""
    return all(r <= p.star_redchi_max for r in redchi) and all(
        abs(z) <= p.star_z_max for zs in z_annuli for z in zs
    )


def pair_ok(z_p: float, p: Params = P) -> bool:
    return math.isfinite(z_p) and abs(z_p) < p.pair_z_max


def injection_ok(z_p: float, recovered: float, p: Params = P) -> bool:
    return z_p >= p.inj_z_min and recovered >= p.inj_r_min


def worth_validation(h: str, s: bool, pr: bool, i20: bool, z21: float, p: Params = P) -> bool:
    return h.startswith("pass") and s and pr and i20 and z21 >= p.inj_z_min


def combine(signals: list[float], noises: list[float]) -> float:
    """z over exposures: sum of signals / quadrature sum of noises."""
    return sum(signals) / math.sqrt(sum(n * n for n in noises))
