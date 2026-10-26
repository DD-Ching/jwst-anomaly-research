"""Cosmic Gems arc (SPT-CL J0615-5746) vs the HST-era RELICS Lenstool v1 and GLAFIC v2 maps.

Frozen plan: ``results/lens_consistency/spt0615_gems/PLAN.md`` (D-034 addendum). Per map and
redshift: does the tangential critical curve cross the arc between its mirror clumps (rule 1,
cross-checked by the parity split, rule 2), and where is the third image predicted relative to
the JWST counter-image candidate (rule 3); joint flag (rule 4). Gates P and frame F are read
from the #537 run (``results/lens_consistency/spt0615_images/summary.json``).

    python scripts/spt0615_gems.py [--grid-step 0.25]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/spt0615_gems")
GATES = Path("results/lens_consistency/spt0615_images/summary.json")
MODELS = ("spt0615-relics-lenstool", "spt0615-relics-glafic")
Z_GRID = (10.0, 10.2, 10.4)  # arc photo-z 10.2 +- 0.2 (95 %), Bradley et al.
N_SAMPLES = 401  # points on the A.1 -> A.2 segment
COUNTER_MIN_ARCSEC = (
    5.0  # ASSUMPTION: a predicted image farther than this from the arc is a counter-image
)
CONSISTENT_ARCSEC = (
    2.2  # ASSUMPTION: largest JWST-era model offset from the candidate (Bradley et al.)
)
OFFSET_ARCSEC = 5.0  # ASSUMPTION: beyond this, no counterpart
TCC_HALF_ARCSEC = 15.0  # nearest-critical-curve search box half-width around the symmetry point
TCC_STEP_ARCSEC = 0.1

# Adamo et al. 2024 (arXiv:2401.03224v2) Extended Data Table 3 (deg, J2000); observed
CLUMPS = {
    "A.1": (93.979828, -57.772475),
    "B.1": (93.979698, -57.772395),
    "C.1": (93.979660, -57.772373),
    "D.1": (93.979635, -57.772355),
    "E.1": (93.979599, -57.772335),
    "A.2": (93.979316, -57.772180),
    "B.2": (93.979416, -57.772235),
    "C.2+D.2": (93.979460, -57.772256),
}
# Bradley et al. 2025 (arXiv:2404.10770v2), text; observed (candidate), model_prediction (others)
CANDIDATE = (93.95000245, -57.77021846)
JWST_PREDICTIONS = {
    "lenstool_a": (93.9490607, -57.7701814),
    "glafic": (93.9504865, -57.7696559),
}


def sign_changes(v: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Positions (linear interpolation in ``t``) where ``v`` changes sign between samples."""
    v = np.asarray(v, float)
    t = np.asarray(t, float)
    s = np.where(v >= 0, 1, -1)  # an exact zero counts as positive: one crossing, not none
    k = np.nonzero(s[:-1] != s[1:])[0]
    return t[k] + (t[k + 1] - t[k]) * v[k] / (v[k] - v[k + 1])


def classify_crossing(crossings: np.ndarray, t_inner: tuple[float, float]) -> str:
    """Rule 1: ``through`` (odd count inside ``t_inner``), ``between`` (odd on the whole
    segment only) or ``miss`` (even count on the whole segment)."""
    c = np.asarray(crossings, float)
    lo, hi = sorted(t_inner)
    if int(np.sum((c >= lo) & (c <= hi))) % 2 == 1:
        return "through"
    return "between" if len(c) % 2 == 1 else "miss"


def parity_split(mu_a1: float, mu_a2: float) -> bool:
    """Rule 2: the two outer clumps have opposite parity."""
    return bool(np.sign(mu_a1) * np.sign(mu_a2) < 0)


def classify_counter(d_arcsec: float | None) -> str:
    """Rule 3 on the smaller candidate distance of the two source points (None: no
    counter-image)."""
    if d_arcsec is None or not np.isfinite(d_arcsec) or d_arcsec > OFFSET_ARCSEC:
        return "none"
    return "consistent" if d_arcsec <= CONSISTENT_ARCSEC else "offset"


def joint_flag(crossing: dict[str, dict[str, str]], counter: dict[str, dict[str, str]]) -> bool:
    """Rule 4: every usable map ``miss`` at every z, or every usable map ``none`` at every z."""
    if not crossing:
        return False
    all_miss = all(all(c == "miss" for c in per_z.values()) for per_z in crossing.values())
    all_none = all(all(c == "none" for c in per_z.values()) for per_z in counter.values())
    return all_miss or all_none


def project(x: np.ndarray, y: np.ndarray, p0: tuple, p1: tuple) -> np.ndarray:
    """Fractional position of points on the line p0 -> p1 (0 at p0, 1 at p1)."""
    d = np.array(p1) - np.array(p0)
    return ((np.asarray(x) - p0[0]) * d[0] + (np.asarray(y) - p0[1]) * d[1]) / float(d @ d)


def lam_t(model: lensmodel.LensModel, x: np.ndarray, y: np.ndarray, z: float) -> np.ndarray:
    ra, dec = model.to_sky(x, y)
    e = model.evaluate(ra, dec, z)
    return 1.0 - np.asarray(e["kappa"], float) - np.asarray(e["gamma"], float)


def arc_table() -> Table:
    t = Table(
        rows=[(k, *v) for k, v in CLUMPS.items()],
        names=("image_id", "ra", "dec"),
        dtype=(str, float, float),
    )
    t.meta["provenance"] = "observed"
    t.meta["source"] = "Adamo et al. 2024 (arXiv:2401.03224v2) Extended Data Table 3, transcribed"
    return t


def run_model(name: str, grid_step: float) -> dict:
    model, _, _ = lc.load_model(name)
    frame = lc.apply_frame_offset(name, model)
    arc = arc_table()
    pts = Table(
        {
            "ra": np.r_[arc["ra"], CANDIDATE[0], [v[0] for v in JWST_PREDICTIONS.values()]],
            "dec": np.r_[arc["dec"], CANDIDATE[1], [v[1] for v in JWST_PREDICTIONS.values()]],
        }
    )
    half = lc.grid_half_width(model, pts)
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, half, grid_step), half, grid_step
    )
    fx, fy = model.to_frame(np.asarray(arc["ra"]), np.asarray(arc["dec"]))
    pos = {k: (float(a), float(b)) for k, a, b in zip(arc["image_id"], fx, fy, strict=True)}
    p0, p1 = pos["A.1"], pos["A.2"]
    t = np.linspace(0.0, 1.0, N_SAMPLES)
    sx, sy = p0[0] + t * (p1[0] - p0[0]), p0[1] + t * (p1[1] - p0[1])
    t_e = float(project(*pos["E.1"], p0, p1))
    t_cd = float(project(*pos["C.2+D.2"], p0, p1))
    t_sym = 0.5 * (t_e + t_cd)
    seg_len = float(np.hypot(p1[0] - p0[0], p1[1] - p0[1]))
    sym = (p0[0] + t_sym * (p1[0] - p0[0]), p0[1] + t_sym * (p1[1] - p0[1]))
    cand = model.to_frame(np.array([CANDIDATE[0]]), np.array([CANDIDATE[1]]))
    cand = (float(cand[0][0]), float(cand[1][0]))
    jw = {
        k: tuple(float(c[0]) for c in model.to_frame(np.array([v[0]]), np.array([v[1]])))
        for k, v in JWST_PREDICTIONS.items()
    }
    g = np.arange(-TCC_HALF_ARCSEC, TCC_HALF_ARCSEC + 1e-9, TCC_STEP_ARCSEC)
    gx, gy = np.meshgrid(sym[0] + g, sym[1] + g)
    arc_mid = (0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]))
    out: dict = {
        "frame_offset_arcsec": list(frame),
        "segment_arcsec": round(seg_len, 3),
        "per_z": {},
    }
    for z in Z_GRID:
        lt = lam_t(model, sx, sy, z)
        cr = sign_changes(lt, t)
        cls = classify_crossing(cr, (t_e, t_cd))
        ra_a, dec_a = model.to_sky(np.array([p0[0], p1[0]]), np.array([p0[1], p1[1]]))
        mu = np.asarray(model.evaluate(ra_a, dec_a, z)["magnification"], float)
        par = parity_split(mu[0], mu[1])
        if par != (len(cr) % 2 == 1):
            raise RuntimeError(
                f"{name} z {z}: parity split {par} disagrees with {len(cr)} crossings"
            )
        lg = lam_t(model, gx.ravel(), gy.ravel(), z).reshape(gx.shape)
        flip = np.zeros(lg.shape, bool)
        flip[:, :-1] |= np.sign(lg[:, :-1]) != np.sign(lg[:, 1:])
        flip[:-1, :] |= np.sign(lg[:-1, :]) != np.sign(lg[1:, :])
        d_tcc = (
            float(np.min(np.hypot(gx[flip] - sym[0], gy[flip] - sym[1]))) if flip.any() else None
        )
        lm = model.lens_map(np.array([p0[0], p1[0]]), np.array([p0[1], p1[1]]), z)
        sources = {}
        best = None
        for k, bx, by in zip(("A.1", "A.2"), lm["beta_x"], lm["beta_y"], strict=True):
            im = lensmodel.find_images(model, grid, float(bx), float(by), z)
            ix, iy = np.asarray(im["x"], float), np.asarray(im["y"], float)
            far = np.hypot(ix - arc_mid[0], iy - arc_mid[1]) > COUNTER_MIN_ARCSEC
            rec = {"n_images": len(im), "n_counter": int(far.sum())}
            if far.any():
                d = np.hypot(ix[far] - cand[0], iy[far] - cand[1])
                j = int(np.argmin(d))
                cx, cy = ix[far][j], iy[far][j]
                cra, cdec = model.to_sky(np.array([cx]), np.array([cy]))
                rec.update(
                    {
                        "d_candidate_arcsec": round(float(d[j]), 3),
                        "counter_ra": float(cra[0]),
                        "counter_dec": float(cdec[0]),
                        "counter_mu": round(
                            float(np.asarray(im["magnification"], float)[far][j]), 3
                        ),
                        **{
                            f"d_{n}_arcsec": round(float(np.hypot(cx - v[0], cy - v[1])), 3)
                            for n, v in jw.items()
                        },
                    }
                )
                best = (
                    rec["d_candidate_arcsec"]
                    if best is None
                    else min(best, rec["d_candidate_arcsec"])
                )
            sources[k] = rec
        out["per_z"][f"{z:.1f}"] = {
            "n_crossings_A1_A2": len(cr),
            "crossing_offset_from_symmetry_arcsec": [
                round(float((c - t_sym) * seg_len), 3) for c in cr
            ],
            "crossing": cls,
            "mu_A1": round(float(mu[0]), 2),
            "mu_A2": round(float(mu[1]), 2),
            "parity_split": par,
            "nearest_tcc_from_symmetry_arcsec": None if d_tcc is None else round(d_tcc, 2),
            "sources": sources,
            "d_candidate_min_arcsec": best,
            "counter": classify_counter(best),
        }
    out["inner_t"] = [round(t_e, 4), round(t_cd, 4)]
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid-step", type=float, default=0.25)
    args = ap.parse_args(argv)
    gates = json.loads(GATES.read_text())["models"]
    summary: dict = {
        "plan": str(OUT / "PLAN.md"),
        "provenance": "positions observed; crossings, magnifications, predicted images "
        "model_prediction; classes derived; thresholds ASSUMPTION",
        "rules": {
            "Z_GRID": Z_GRID,
            "N_SAMPLES": N_SAMPLES,
            "COUNTER_MIN_ARCSEC": COUNTER_MIN_ARCSEC,
            "CONSISTENT_ARCSEC": CONSISTENT_ARCSEC,
            "OFFSET_ARCSEC": OFFSET_ARCSEC,
            "TCC_HALF_ARCSEC": TCC_HALF_ARCSEC,
            "TCC_STEP_ARCSEC": TCC_STEP_ARCSEC,
            "grid_step_arcsec": args.grid_step,
        },
        "candidate": CANDIDATE,
        "models": {},
    }
    crossing, counter = {}, {}
    for name in MODELS:
        usable = bool(gates[name]["gate_P"]["passes"] and gates[name]["frame_F"]["passes"])
        rec = {"usable": usable, "gates_from": str(GATES)}
        if usable:
            rec.update(run_model(name, args.grid_step))
            crossing[name] = {z: v["crossing"] for z, v in rec["per_z"].items()}
            counter[name] = {z: v["counter"] for z, v in rec["per_z"].items()}
        summary["models"][name] = rec
    summary["joint_flag"] = joint_flag(crossing, counter)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    sep = SkyCoord(*CANDIDATE, unit="deg").separation(SkyCoord(*CLUMPS["A.1"], unit="deg")).arcsec
    print(json.dumps(summary, indent=1))
    print(f"candidate - A.1 separation {sep:.1f} arcsec (not stored)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
