"""Order-of-magnitude scales from published wormhole-support formulas (support screen, 2026-10-29).

Arithmetic only: no solver, no evolution, no scan. Each entry records the published formula (paper
and equation), its inputs with units, the assumptions and the outputs. Formulas printed in natural
units (hbar = c = 1, l_p = sqrt(G_N)) are converted with hbar*c and l_p below.
Record: docs/plans/2026-10-29-wh-support-screen.md.

  python scripts/wh_support_scales.py   # -> results/wh_support_screen/scales.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "results" / "wh_support_screen"

# CODATA 2018 / IAU values
C = 299_792_458.0  # m s^-1
G = 6.67430e-11  # m^3 kg^-1 s^-2
HBAR_C_EV_M = 1.973269804e-7  # eV m
L_P = 1.616255e-35  # m
LY = 9.4607304725808e15  # m
YR = 3.15576e7  # s (Julian year)
M_SUN = 1.98841e30  # kg


def mmp_binding_times_radius(g: float) -> float:
    """|E_min| r_e / (hbar c) for MMP (arXiv:1807.04726).

    (2.3) r_e^2 = pi q^2 l_p^2 / g^2 and (5.31) E_min = -G_N q^2 / (256 r_e^3) combine to
    g^2 / (256 pi), independent of q.
    """
    return g**2 / (256.0 * math.pi)


def mmp_binding_from_q(q: float, g: float) -> tuple[float, float]:
    """(r_e in m, |E_min| in eV) from MMP (2.3) and (5.31) evaluated separately (unit test)."""
    r_e = math.sqrt(math.pi) * q * L_P / g
    e_min = L_P**2 * q**2 / (256.0 * r_e**3)  # 1/m, natural units
    return r_e, e_min * HBAR_C_EV_M


def mm_parameters(r_e: float, r5: float) -> dict:
    """Maldacena-Milekhin (arXiv:2008.06618), formulas after their (3.26), SI inputs.

    c2 = 4 pi R5 r_e / (G4 sqrt(L)), l = (4/pi) (r_e^2/R5) sqrt(L),
    E_bin = -(pi^2/16) R5^2 / (G4 r_e L), with L = ln(r_e/R5).
    """
    lg = math.log(r_e / r5)
    c2 = 4.0 * math.pi * r5 * r_e / (L_P**2 * math.sqrt(lg))
    ell = (4.0 / math.pi) * r_e**2 / r5 * math.sqrt(lg)
    e_bin_kg = -(math.pi**2 / 16.0) * r5**2 * C**2 / (G * r_e * lg)
    return {"c2": c2, "ell_m": ell, "ell_ly": ell / LY, "E_bin_kg": e_bin_kg}


def scales() -> dict:
    out = {"units": "SI unless stated; energies in eV or kg (E / c^2)", "entries": {}}
    e = out["entries"]

    e["mmp_max_signal_energy"] = {
        "source": (
            "Maldacena, Milekhin & Popov, arXiv:1807.04726: (2.3), (5.31); section 5.4: a "
            "signal above the throat binding energy makes a near-extremal black hole"
        ),
        "formula": (
            "E_max ~ |E_min| = G_N q^2/(256 r_e^3) with r_e^2 = pi q^2 l_p^2/g^2, "
            "so E_max r_e/(hbar c) = g^2/(256 pi) [der]"
        ),
        "assumptions": (
            "U(1) coupling g = 0.30 (electromagnetic) or 0.36 (hypercharge near the weak "
            "scale) [assumption]; MMP give the threshold only to order of magnitude"
        ),
        "outputs": {
            f"E_max_r_e_over_hbar_c_g{g}": mmp_binding_times_radius(g) for g in (0.30, 0.36)
        },
        "reading": (
            "one signal may carry about 1e-4 of the energy of a single quantum whose "
            "wavelength is the mouth radius r_e"
        ),
    }

    d_max = HBAR_C_EV_M / 1e12  # 1/TeV in m
    e["mmp_standard_model_size"] = {
        "source": (
            "arXiv:1807.04726 section 8: with Standard Model fields the distance d between "
            "the mouths must be below the electroweak scale (say 1/TeV), r_e much smaller"
        ),
        "formula": "d_max = hbar c / (1 TeV) [der: unit conversion]",
        "outputs": {"d_max_m": d_max},
        "reading": (
            "with Standard Model matter alone the wormhole is sub-attometre; larger ones "
            "need new massless charged fields"
        ),
    }

    r_e, r5 = 1.5e7, 50e-6  # MM (3.26) minimal r_e for a human traveller; R5 = 50 micron
    mm = mm_parameters(r_e, r5)
    d_min_ly = 0.5  # MM (B.32) for the (3.27) parameters
    transit_yr = math.pi * mm["ell_m"] / C / YR
    e["mm_humanly_traversable"] = {
        "source": (
            "Maldacena & Milekhin, arXiv:2008.06618: formulas after (3.26), numbers (3.27), "
            "orbit condition (B.32) d > 0.5 ly, information <= ~c2 quanta per l (section 2)"
        ),
        "inputs": {"r_e_m": r_e, "R5_m": r5, "d_min_ly": d_min_ly},
        "outputs": {
            **mm,
            "outside_transit_pi_l_yr": transit_yr,
            "exterior_light_time_at_d_min_yr": d_min_ly,
            "transit_over_exterior_at_d_min": transit_yr / d_min_ly,
            "E_bin_J": abs(mm["E_bin_kg"]) * C**2,
            "mouth_mass_r_e_c2_over_G_kg": r_e * C**2 / G,
            "mouth_mass_solar": r_e * C**2 / G / M_SUN,
        },
        "assumptions": (
            "a Randall-Sundrum II dark sector with a gauged dark U(1) (the authors call it "
            "science fiction); no formation mechanism; mouths kept apart by orbiting"
        ),
    }

    garattini_lp, hps_lp = 1.0166, (0.02, 67.0)
    e["casimir_static_throats"] = {
        "source": (
            "Garattini, arXiv:1907.03623: (29), (41)-(42), footnote 6 (omega = 3 gives "
            "r0 = 1.0166 l_P); Hochberg, Popov & Sushkov, gr-qc/9701064: numerical throat "
            "r(0) = 0.02 l_P and an analytic example r(0) = 67 l_P from their Eq. (8)"
        ),
        "formula": "r0 [m] = (r0 / l_P) * l_P [der: unit conversion]",
        "outputs": {
            "garattini_r0_m": garattini_lp * L_P,
            "hps_r0_m_range": [x * L_P for x in hps_lp],
        },
        "reading": (
            "the static Casimir- or vacuum-polarisation-supported throats in these papers "
            "are Planck scale; neither paper treats dynamics or a signal"
        ),
    }
    return out


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1) + "\n")


def main() -> int:
    res = scales()
    write_json(OUT / "scales.json", res)
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
