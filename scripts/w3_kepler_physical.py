"""W3 / Kepler: physical mapping of the tested grid and the population requirement (#559).

Plan: ``results/w3_kepler/physical/PLAN.md``. No light curve, image or injection: the stars'
radii and distances come from the official Kepler DR25 stellar table, the efficiencies from
``results/w3_kepler/s2000/combined.json`` (``round_plus_546``).

Stages (``python scripts/w3_kepler_physical.py <stage>``):

- ``stellar``: DR25 radius and distance (with errors) of the 2,000 searched stars (ledger);
- ``map``: for every star and W3 cell, theta_*, theta_E, mu_rel and t_* = rho t_E; for the
  scenario speeds, the lens distance fraction x = v t_* / R_*, D_L, R_E and |M| (``mapping.json``);
- ``population``: the uniform lens mass density needed for N = 1 or 3 expected detections per
  tested t_E, now and with +500 / +5,000 stars (``population.json``).

Labels: radii and distances are catalogue (model) values; everything else is derived. Nothing here
is an abundance, a physical event-rate limit or a discovery probability.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import sys
import urllib.parse
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import net_ledger as nl  # noqa: E402
import w3_kepler as wk  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "w3_kepler" / "physical"
CAP = 150_000_000  # owner limit for this round, bytes
TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
COLS = "kepid,radius,radius_err1,radius_err2,dist,dist_err1,dist_err2,mass,teff,logg,kepmag"

AU_KM = 1.495978707e8
R_SUN_AU = 6.957e5 / AU_KM  # 0.00465 AU
PC_AU = 206264.806
SCHW4_AU = 4 * 6.674e-11 * 1.98892e30 / 2.99792458e8**2 / 1.495978707e11  # 4 G M_sun / c^2 in AU
FOUR_PI_G_C2 = 4 * math.pi * 6.674e-11 * 1.98892e30 / 2.99792458e8**2 / 3.0857e16  # pc / M_sun
DAY_YR = 1 / 365.25
SPEEDS = {"S1_halo_like_220": 220.0, "S2_disk_like_30": 30.0}  # km/s; S2 is an ASSUMPTION
REF = {
    "rho_dm_extrapolated_read2014": [0.005, 0.015],  # M_sun / pc^3 (Read 2014, conclusions)
    "rho_baryon_flynn2006_via_read2014": [0.0914, 0.009],  # value, error (M_sun / pc^3)
}


def searched() -> list[int]:
    """The 2,000 stars searched under one rule: #542 / #546 and #554 batches 1-3."""
    k = ROOT / "results" / "w3_kepler"
    a = json.loads((k / "expand" / "sample_expand.json").read_text())["stars"]
    s = json.loads((k / "s2000" / "sample_2000.json").read_text())["stars"]
    return list(a) + [x["kepid"] for x in s if x["batch"] <= 3]


# ----------------------------------------------------------------------------------------- stellar


def stage_stellar(_a) -> None:
    lg = nl.Ledger(OUT / "net_ledger.jsonl", CAP)
    kids = searched()
    rows = []
    for i in range(0, len(kids), 500):
        part = ",".join(str(k) for k in kids[i : i + 500])
        q = (
            f"select {COLS} from keplerstellar where st_delivname='q1_q17_dr25_stellar'"
            f" and kepid in ({part})"
        )
        b = nl.fetch(TAP + "?query=" + urllib.parse.quote(q) + "&format=csv", lg, 200_000)
        rows += list(csv.DictReader(io.StringIO(b.decode())))
    missing = set(kids) - {int(r["kepid"]) for r in rows}
    if missing:  # a TAP error page or skipped kepids must not pass as a complete table
        raise SystemExit(f"{len(missing)} searched stars missing from the DR25 response")
    rows.sort(key=lambda r: int(r["kepid"]))
    (OUT / "inputs").mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / "inputs" / "stellar_dr25.csv.gz", "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS.split(","))
        w.writeheader()
        w.writerows(rows)
    raw = gzip.open(OUT / "inputs" / "stellar_dr25.csv.gz", "rb").read()
    print(len(kids), "asked;", len(rows), "rows;", hashlib.sha256(raw).hexdigest(), lg.used())


def load_stellar(complete: bool = True) -> list[dict]:
    """Catalogue rows as floats (empty -> NaN); complete=True keeps rows with radius, distance
    and both errors of each."""
    with gzip.open(OUT / "inputs" / "stellar_dr25.csv.gz", "rt") as fh:
        rows = [
            {k: float(v) if v else float("nan") for k, v in r.items()} for r in csv.DictReader(fh)
        ]
    need = ("radius", "dist", "radius_err1", "radius_err2", "dist_err1", "dist_err2")
    return [r for r in rows if not complete or all(np.isfinite(r[k]) for k in need)]


# --------------------------------------------------------------------------------------- mapping


def theta_star_uas(r_sun: float, d_pc: float) -> float:
    """Angular source radius in micro-arcseconds."""
    return r_sun * R_SUN_AU / (d_pc / 1000.0) * 1000.0  # 1 AU at 1 kpc = 1 mas


def cell_map(r_sun: float, d_pc: float, rho: float, te_d: float, v_kms: float) -> dict:
    """PLAN section 1 for one star, cell and speed. x = v t_* / R_* (lens distance fraction);
    reachable only for 0 < x < 1; then D_L = x D_S, R_E = v t_E and |M| = R_E^2 / (4GM_sun/c^2
    D_L (1 - x))."""
    th_s = theta_star_uas(r_sun, d_pc)
    th_e = th_s / rho
    t_star = rho * te_d
    mu = th_e / 1000.0 / (te_d * DAY_YR)  # mas / yr
    v_au_d = v_kms * 86400.0 / AU_KM
    x = v_au_d * t_star / (r_sun * R_SUN_AU)
    out = {
        "theta_star_uas": th_s,
        "theta_E_uas": th_e,
        "mu_rel_mas_yr": mu,
        "t_star_d": t_star,
        "x": x,
    }
    if 0 < x < 1:
        d_l_au = x * d_pc * PC_AU
        r_e = v_au_d * te_d
        out |= {"D_L_pc": x * d_pc, "R_E_au": r_e, "M_msun": r_e**2 / (SCHW4_AU * d_l_au * (1 - x))}
    return out


def theta_star_corners(s: dict) -> tuple[tuple[float, float], tuple[float, float]]:
    """(R, D) giving the largest and the smallest theta_*: R high + D low, and the reverse.
    Magnitudes of the errors are used: in DR25, dist_err2 is stored positive (all 1,971 stars)."""
    r, d = s["radius"], s["dist"]
    return (
        (r + abs(s["radius_err1"]), d - abs(s["dist_err2"])),
        (r - abs(s["radius_err2"]), d + abs(s["dist_err1"])),
    )


def _pct(a) -> list | None:
    a = np.asarray([v for v in a if v is not None and np.isfinite(v)])
    if not len(a):
        return None
    return [float(f"{v:.4g}") for v in np.percentile(a, [5, 50, 95])]


def stage_map(_a) -> None:
    st = load_stellar()
    cells = []
    for c, (rho, blend, te) in enumerate(wk.CELLS):
        rec = {"cell": c, "rho": rho, "blend": blend, "tE_d": te, "t_star_d": rho * te}
        base = [cell_map(s["radius"], s["dist"], rho, te, 220.0) for s in st]
        rec["theta_star_uas_p5_50_95"] = _pct([b["theta_star_uas"] for b in base])
        rec["theta_E_uas_p5_50_95"] = _pct([b["theta_E_uas"] for b in base])
        rec["mu_rel_mas_yr_p5_50_95"] = _pct([b["mu_rel_mas_yr"] for b in base])
        for name, v in SPEEDS.items():
            m = [cell_map(s["radius"], s["dist"], rho, te, v) for s in st]
            # catalogue-error range of theta_* moves x and |M| (theta_star_corners)
            corners = [theta_star_corners(s) for s in st]
            hi = [cell_map(*c[0], rho, te, v) for c in corners]
            lo = [cell_map(*c[1], rho, te, v) for c in corners]
            rec[name] = {
                "share_reachable_x_lt_1": round(sum("M_msun" in r for r in m) / len(m), 4),
                "x_p5_50_95": _pct([r["x"] for r in m]),
                "D_L_pc_p5_50_95": _pct([r.get("D_L_pc") for r in m]),
                "M_msun_p5_50_95": _pct([r.get("M_msun") for r in m]),
                "M_msun_p5_50_95_theta_star_high": _pct([r.get("M_msun") for r in hi]),
                "M_msun_p5_50_95_theta_star_low": _pct([r.get("M_msun") for r in lo]),
            }
        cells.append(rec)
    out = {
        "label": "derived from catalogue radii and distances (DR25 model values); no light curve",
        "plan": "results/w3_kepler/physical/PLAN.md",
        "stars": len(st),
        "stars_without_complete_radius_distance": len(load_stellar(False)) - len(st),
        "radius_rsun_p5_50_95": _pct([s["radius"] for s in st]),
        "dist_pc_p5_50_95": _pct([s["dist"] for s in st]),
        "speeds_kms": SPEEDS,
        "cells": cells,
    }
    (OUT / "mapping.json").write_text(json.dumps(out, indent=1) + "\n")
    for r in cells:
        s1, s2 = r["S1_halo_like_220"], r["S2_disk_like_30"]
        print(
            r["cell"],
            r["rho"],
            r["tE_d"],
            r["t_star_d"],
            s1["share_reachable_x_lt_1"],
            s1["M_msun_p5_50_95"],
            s2["share_reachable_x_lt_1"],
            s2["M_msun_p5_50_95"],
        )


# ------------------------------------------------------------------------------------ population


def tau_per_density(d_pc: np.ndarray) -> np.ndarray:
    """Einstein-radius optical depth per unit lens mass density (M_sun / pc^3) for a uniform density
    between observer and source: (4 pi G / c^2) D_S^2 / 6."""
    return FOUR_PI_G_C2 * np.asarray(d_pc) ** 2 / 6.0


def detections_per_density(d_pc, t_yr, eps: float, te_d: float) -> float:
    """N per unit |rho_neg| = sum_i eps (4/pi) tau_i T_i / t_E (umbra crossings, u0 < 2)."""
    return float(
        eps * 4 / math.pi * (tau_per_density(d_pc) * np.asarray(t_yr)).sum() / (te_d * DAY_YR)
    )


def tested_rho_share(r_sun: np.ndarray, te_d: float, v_kms: float, lo=0.01, hi=0.1) -> np.ndarray:
    """Share of each star's optical depth from lenses whose rho = x R_* / (v t_E) lies in [lo, hi]
    (single speed v): the x(1 - x) weight integrated over that x interval, over its full-sightline
    value 1/6."""
    r_e_scale = v_kms * 86400.0 / AU_KM * te_d / (np.asarray(r_sun) * R_SUN_AU)  # x per unit rho
    x1 = np.clip(lo * r_e_scale, 0, 1)
    x2 = np.clip(hi * r_e_scale, 0, 1)

    def prim(x):
        return x**2 / 2 - x**3 / 3

    return (prim(x2) - prim(x1)) * 6.0


POOL_UNREAD_AFTER_554 = 17740 - 1500  # sample_2000.json candidates minus the 1,500 read in #554


def _inv(n: float, per: float) -> float | None:
    """Density giving n expected detections, or None when the rate is zero (eps 0 at some cell)."""
    return float(f"{n / per:.4g}") if per > 0 else None


def stage_population(_a) -> None:
    st = load_stellar()
    d = np.array([s["dist"] for s in st])
    cmb = json.loads((ROOT / "results" / "w3_kepler" / "s2000" / "combined.json").read_text())
    cells = cmb["round_plus_546"]["cells"]
    win = {c["sampling_star_years"] for c in cells}
    assert len(win) == 1, win  # every searched star has the same window (1459.49 d)
    t = np.full(len(d), win.pop() / cmb["round_plus_546"]["stars"])
    rows = []
    for te in sorted({c["tE_d"] for c in cells}):
        eps = [c["recovery_rate"] for c in cells if c["tE_d"] == te]
        rec = {"tE_d": te, "eps_min": min(eps), "eps_max": max(eps)}
        for n_star, scale in (
            ("now_2000", 1.0),
            ("plus_500", 2500 / 2000),
            ("plus_5000", 7000 / 2000),
        ):
            per_lo = detections_per_density(d, t, min(eps), te) * scale
            per_hi = detections_per_density(d, t, max(eps), te) * scale
            rec[n_star] = {
                "N_per_unit_density": [float(f"{per_lo:.4g}"), float(f"{per_hi:.4g}")],
                "density_for_N1": [_inv(1.0, per_hi), _inv(1.0, per_lo)],
                "density_for_N3": [_inv(2.996, per_hi), _inv(2.996, per_lo)],
            }
        best = detections_per_density(d, t, max(eps), te)
        rec["star_factor_needed_vs_now"] = {
            f"N{label}_at_{name}": float(f"{n / best / ref:.4g}")
            for label, n in ((1, 1.0), (3, 2.996))
            for name, ref in (("rho_dm_high_0.015", 0.015), ("rho_baryon_0.0914", 0.0914))
        }
        # the floor that holds whatever happens outside the tested rho range: every umbra
        # crossing along the whole sight line detected (eps = 1)
        one = detections_per_density(d, t, 1.0, te)
        pool = (2000 + POOL_UNREAD_AFTER_554) / 2000
        rec["eps1_floor"] = {
            "density_for_N1_now": float(f"{1 / one:.4g}"),
            "density_for_N1_whole_pool": float(f"{1 / one / pool:.4g}"),
            "density_for_N1_whole_pool_at_eps_max": float(f"{1 / best / pool:.4g}"),
        }
        # lenses restricted to the tested rho range (0.01-0.1), measured eps, single speed
        r = np.array([x["radius"] for x in st])
        rec["tested_rho_only"] = {}
        for name, v in SPEEDS.items():
            share = tested_rho_share(r, te, v)
            w = (tau_per_density(d) * t * share).sum()
            n_lo = min(eps) * 4 / math.pi * w / (te * DAY_YR)
            n_hi = max(eps) * 4 / math.pi * w / (te * DAY_YR)
            rec["tested_rho_only"][name] = {
                "share_of_tau": float(f"{w / (tau_per_density(d) * t).sum():.4g}"),
                "density_for_N1": [_inv(1.0, n_hi), _inv(1.0, n_lo)] if n_hi > 0 else None,
            }
        rows.append(rec)
    out = {
        "label": "derived (scenario requirement, not an abundance or a limit)",
        "plan": "results/w3_kepler/physical/PLAN.md",
        "assumptions": [
            "uniform lens mass density |rho_neg| (M_sun/pc^3) between observer and each star",
            "every event has the stated t_E, u0 < 2; rows: (a) the measured eps applied to lenses"
            " at every distance (the rho of most of them is outside the tested 0.01-0.1, where"
            " eps is not measured), (b) eps = 1 at every distance (a floor that holds whatever"
            " eps is outside the tested range), (c) only lenses whose rho is in 0.01-0.1, one"
            " speed",
            "efficiency between the lowest and highest measured cell at that t_E (2,000 stars)",
            "SAP-level efficiency applies to real events (sky-to-SAP response not verified)",
            "+500 / +5,000: zero result, similar efficiency and the same distance distribution",
        ],
        "stars": len(st),
        "stars_left_out_no_distance": len(load_stellar(False)) - len(st),
        "pool_unread_after_554": POOL_UNREAD_AFTER_554,
        # N per density scales as 1 / t_E: with eps = 1 and these 2,000 stars, the t_E at which one
        # expected detection needs only the reference density, and the umbra length 4 t_E (u0 = 0)
        # in 29.4-min cadences (efficiency below t_E 0.1 d is not measured)
        "eps1_tE_needed_for_N1_at_reference": {
            name: {
                "tE_d": float(f"{0.1 * ref * detections_per_density(d, t, 1.0, 0.1):.4g}"),
                "umbra_max_cadences": float(
                    f"{4 * 0.1 * ref * detections_per_density(d, t, 1.0, 0.1) / wk.CAD:.3g}"
                ),
            }
            for name, ref in (("rho_baryon_0.0914", 0.0914), ("rho_dm_high_0.015", 0.015))
        },
        "sum_D_S2_T_pc2_yr": float(f"{(d**2 * t).sum():.4g}"),
        "references_msun_pc3": REF,
        "rows": rows,
    }
    (OUT / "population.json").write_text(json.dumps(out, indent=1) + "\n")
    for r in rows:
        print(
            r["tE_d"],
            r["eps_min"],
            r["eps_max"],
            r["now_2000"]["density_for_N1"],
            r["now_2000"]["density_for_N3"],
            r["plus_5000"]["density_for_N3"],
        )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)
    for name, fn in (
        ("stellar", stage_stellar),
        ("map", stage_map),
        ("population", stage_population),
    ):
        p = sub.add_parser(name)
        p.set_defaults(fn=fn)
    a = ap.parse_args(argv)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
