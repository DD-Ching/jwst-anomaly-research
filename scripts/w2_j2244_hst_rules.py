"""Frozen rules of the J2244-0550 HST test (results/w2_j2244_hst/PLAN.md), as pure functions.

Imported by scripts/w2_j2244_hst.py (run in an isolated environment with lenstronomy) and tested in
the project environment (tests/test_w2_j2244_hst_rules.py). Every threshold is an ASSUMPTION named
in Params.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Params:
    t_min: float = (
        25.0  # floor of the detection threshold on delta chi^2 (about 5 sigma for one parameter)
    )
    t_factor: float = 2.0  # T = max(t_min, t_factor * largest development-negative delta chi^2)
    centre_max: float = 0.6  # extended-light centre within this of the image midpoint (arcsec)
    redchi_max: float = 10.0  # v2 (from development): only catastrophic fits are "fit unreliable"
    r_min: dict = field(
        default_factory=lambda: {"F160W": 0.15, "F814W": 0.075}
    )  # v2: resolved Sersic
    pos_tol: float = 0.15  # fitted point sources within this of the initial peaks (arcsec)
    gate_min_neg: int = 5  # validation negatives needed
    gate_min_pos: int = 3  # validation positives needed
    gate_pos_frac: float = 0.8  # fraction of validation positives detected
    absorb_min: float = 0.5  # recovered / injected flux at the brightest injection
    sep_range: tuple = (0.7, 1.8)  # pool geometry: two point sources this far apart (arcsec)
    inj_reff: tuple = (0.2, 0.5)  # injected Sersic r_eff (arcsec)
    inj_n: float = 4.0
    inj_q: float = 0.8
    inj_frac: float = 0.4  # injection on the image line, this fraction from the fainter image
    inj_mags: dict = field(
        default_factory=lambda: {
            "F160W": (19.0, 20.0, 21.0, 22.0, 23.0),
            "F814W": (20.0, 21.0, 22.0, 23.0, 24.0),
        }
    )


P = Params()


def threshold(dev_neg_dchi2: list[float], p: Params = P) -> float:
    """Detection threshold of a band, fixed from the development negatives before validation."""
    finite = [d for d in dev_neg_dchi2 if d is not None and math.isfinite(d)]
    return max(p.t_min, p.t_factor * max(finite)) if finite else p.t_min


def status(
    dchi2: float | None,
    t: float,
    centre_offset: float | None,
    redchi: float | None,
    pos_ok: bool,
    p: Params = P,
    sersic_r: float | None = None,
    sersic_flux: float | None = None,
    band: str | None = None,
) -> str:
    """Outcome: fit unreliable, extended light (offset), unresolved or none."""
    if dchi2 is None or redchi is None or not math.isfinite(dchi2) or not pos_ok:
        return "fit unreliable"
    if redchi > p.redchi_max:
        return "fit unreliable"
    if dchi2 < t or (sersic_flux is not None and sersic_flux <= 0):
        return "none"
    if band is not None and sersic_r is not None and sersic_r < p.r_min[band]:
        return "unresolved"
    if centre_offset is None or centre_offset > p.centre_max:
        return "extended light offset"
    return "extended light"


def detected(st: str) -> bool:
    return st == "extended light"


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 1.0)
    f = k / n
    d = 1 + z * z / n
    c = (f + z * z / (2 * n)) / d
    h = z * math.sqrt(f * (1 - f) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def gate(neg_status: list[str], pos_status: list[str], absorb_ok: bool, p: Params = P) -> dict:
    """Per-band gate on the validation controls (never on development controls)."""
    neg_dec = [s for s in neg_status if s != "fit unreliable"]
    pos_dec = [s for s in pos_status if s != "fit unreliable"]
    false_det = sum(detected(s) for s in neg_dec)
    pos_det = sum(detected(s) for s in pos_dec)
    enough = len(neg_dec) >= p.gate_min_neg and len(pos_dec) >= p.gate_min_pos
    pos_ok = enough and pos_det >= math.ceil(p.gate_pos_frac * len(pos_dec) - 1e-9)
    return {
        "n_neg": len(neg_dec),
        "neg_false_detections": false_det,
        "n_pos": len(pos_dec),
        "pos_detected": pos_det,
        "pos_wilson": wilson(pos_det, len(pos_dec)),
        "neg_false_wilson": wilson(false_det, len(neg_dec)),
        "unreliable": len(neg_status) - len(neg_dec) + len(pos_status) - len(pos_dec),
        "enough_controls": enough,
        "absorb_ok": absorb_ok,
        "passed": bool(enough and false_det == 0 and pos_ok and absorb_ok),
    }


def contiguous_limit(mags: list[float], recovered: list[bool]) -> float | None:
    """Faintest magnitude of the unbroken run of recoveries from the bright end (None if the
    brightest fails)."""
    lim = None
    for m, r in sorted(zip(mags, recovered, strict=True)):
        if not r:
            break
        lim = m
    return lim
