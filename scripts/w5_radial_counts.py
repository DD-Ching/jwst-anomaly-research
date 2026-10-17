"""W5 count deficit stacked at the ``radial`` screen's line-convergence centres.

Plan and decision rule: ``results/w5_radial_counts/PLAN.md``. Fields load as in
``w1_stack_radial.py`` (``inject_radial.load_field`` on ``inject_shear.field_spec``). Counted rows:
background (z_phot > z_l + ``DZ`` where a photo-z exists), magnitude brighter than the field's
``MAG_PCT`` percentile of background rows, and not a ``radial`` arc row. N_obs(R) at the centres
is compared with control positions in the centre's cluster-centric annulus; the power forecast
scales the control sums by ``exotic_sim.count_ratio`` averaged over r <= R.

    python scripts/w5_radial_counts.py --out results/w5_radial_counts [--n-draws 2000]

All outputs are ``derived`` (forecast: ``model_prediction``); thresholds are ASSUMPTIONs (PLAN.md).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))

import inject_shear as ish  # noqa: E402
import w1_stack_radial as wsr  # noqa: E402

from jwst_anomaly import exotic_sim, paths  # noqa: E402

ir = ish.ir
RADII = (3.0, 6.0, 9.0)  # arcsec
MASSES = (1e12, 4e12, 1e13)  # M_sun
DZ = 0.1  # background margin above z_l (ASSUMPTION)
MAG_PCT = 80.0  # counted rows: brighter than this percentile of background mags (ASSUMPTION)
ANNULUS = 10.0  # arcsec half-width of the control annulus in cluster-centric radius
FOOTPRINT = 20.0  # arcsec: a control needs >= 1 counted row within this (footprint proxy)
Z_SRC = 2.0


def counted_rows(shapes, z_l: float, arcs: np.ndarray) -> tuple[np.ndarray, float]:
    """Bool per row (background, brighter than the limit, not a radial arc) and the mag limit."""
    mag = np.asarray(shapes["mag"], float)
    bg = np.ones(len(shapes), bool)
    if "z_phot" in shapes.colnames:
        z = np.asarray(shapes["z_phot"], float)
        bg = ~np.isfinite(z) | (z > z_l + DZ)
    ok = bg & np.isfinite(mag)
    m_lim = float(np.percentile(mag[ok], MAG_PCT))
    return ok & (mag <= m_lim) & ~arcs, m_lim


def mean_ratio(radius: float, theta_e: float, mags: np.ndarray, m_lim: float) -> float:
    """Area mean of ``count_ratio`` over r <= radius (``model_prediction``), with N(>S) the
    empirical counts of ``mags`` (fluxes 10^(-0.4 m)); zero beyond the brightest row."""
    s_sorted = np.sort(10 ** (-0.4 * mags))

    def cum(s):
        return len(s_sorted) - np.searchsorted(s_sorted, np.asarray(s, float), side="left")

    r = np.linspace(radius / 400, radius, 400)
    q = exotic_sim.count_ratio(r / theta_e, cum, 10 ** (-0.4 * m_lim), n=1.0, sign=-1)
    return float(np.sum(q * r) / np.sum(r))


def control_positions(
    rc: float, centres: np.ndarray, tree, radius: float, n: int, rng
) -> np.ndarray:
    """``n`` control positions (x, y) in the annulus rc ± ANNULUS, inside the footprint proxy and
    >= 2 radius from every centre; drawn by rejection (NaN rows if the annulus is exhausted)."""
    out = np.full((n, 2), np.nan)
    got = 0
    for _ in range(200):
        k = 4 * n
        r = np.sqrt(rng.uniform(max(rc - ANNULUS, 0) ** 2, (rc + ANNULUS) ** 2, k))
        phi = rng.uniform(0, 2 * np.pi, k)
        p = np.c_[r * np.cos(phi), r * np.sin(phi)]
        far = cKDTree(centres).query(p)[0] >= 2 * radius
        cover = np.array([len(v) > 0 for v in tree.query_ball_point(p, FOOTPRINT)])
        p = p[far & cover]
        take = min(len(p), n - got)
        out[got : got + take] = p[:take]
        got += take
        if got == n:
            break
    return out


def field_counts(name: str, out: Path, n_draws: int, rng, radii=RADII, near: float = 0.0) -> dict:
    spec = ish.field_spec(name)
    (out / name).mkdir(parents=True, exist_ok=True)
    model, shapes, extra, _, _ = ir.load_field(name, out / name, spec)
    arcs = wsr.radial_arc_rows(model, shapes, spec, extra)
    if near > 0:  # post hoc diagnostic: the arcs' neighbours (e.g. deblended fragments) too
        arcs = wsr.near_rows(model, shapes, arcs, near)
    use, m_lim = counted_rows(shapes, model.z_lens, arcs)
    x, y = model.to_frame(np.asarray(shapes["ra"]), np.asarray(shapes["dec"]))
    xy = np.c_[np.asarray(x)[use], np.asarray(y)[use]]
    tree = cKDTree(xy)
    peaks = json.loads(wsr.RADIAL[name].read_text())["peaks"]
    cx, cy = model.to_frame(np.array([p["ra"] for p in peaks]), np.array([p["dec"] for p in peaks]))
    centres = np.c_[cx, cy]
    in_fp = np.array([len(v) > 0 for v in tree.query_ball_point(centres, FOOTPRINT)])
    res = {
        "field": name,
        "model": spec["model"],
        "z_lens": float(model.z_lens),
        "m_lim": m_lim,
        "n_counted": int(use.sum()),
        "n_arc_rows": int(arcs.sum()),
        "n_centres": len(peaks),
        "n_centres_in_footprint": int(in_fp.sum()),
        "radii": {},
    }
    mags = np.asarray(shapes["mag"], float)[use]
    for rad in radii:
        obs = np.array([len(v) for v in tree.query_ball_point(centres[in_fp], rad)])
        ctrl = np.zeros((int(in_fp.sum()), n_draws))
        for j, (a, b) in enumerate(centres[in_fp]):
            pos = control_positions(np.hypot(a, b), centres, tree, rad, n_draws, rng)
            ok = np.isfinite(pos[:, 0])
            if not ok.any():  # a NaN control sum would read as p_low = 0 in combine
                raise ValueError(f"{name}: no control position for centre {j} at R = {rad}")
            cnt = np.full(n_draws, np.nan)
            cnt[ok] = [len(v) for v in tree.query_ball_point(pos[ok], rad)]
            ctrl[j] = np.where(ok, cnt, np.nanmean(cnt))  # exhausted annulus: mean (rare)
        theta = {m: float(ir.theta_e_arcsec(m, model.z_lens, np.array([Z_SRC]))[0]) for m in MASSES}
        q_h = {m: mean_ratio(rad, t, mags, m_lim) for m, t in theta.items()}
        res["radii"][str(rad)] = {
            "n_obs": int(obs.sum()),
            "per_centre": obs.tolist(),
            "ctrl_sum": ctrl.sum(axis=0).tolist(),
            "theta_e_arcsec": {f"{m:.0e}": t for m, t in theta.items()},
            "q_hypothesis": {f"{m:.0e}": q for m, q in q_h.items()},
        }
    return res


def combine(fields: list[dict]) -> dict:
    """Pooled Q, one-sided p and power per radius (PLAN.md)."""
    out = {}
    for rad in fields[0]["radii"]:
        n_obs = sum(f["radii"][rad]["n_obs"] for f in fields)
        ctrl = np.sum([f["radii"][rad]["ctrl_sum"] for f in fields], axis=0)
        p5 = float(np.percentile(ctrl, 5))
        power = {}
        for m in (f"{v:.0e}" for v in MASSES):
            scaled = np.sum(
                [np.asarray(f["radii"][rad]["ctrl_sum"]) * f["radii"][rad]["q_hypothesis"][m]
                 for f in fields], axis=0)  # fmt: skip
            power[m] = float(np.mean(scaled <= p5))
        out[rad] = {
            "n_obs": n_obs,
            "ctrl_mean": float(ctrl.mean()),
            "ctrl_p5": p5,
            "ctrl_p95": float(np.percentile(ctrl, 95)),
            "q": float(n_obs / ctrl.mean()) if ctrl.mean() > 0 else None,
            "p_low": float(np.mean(ctrl <= n_obs)),
            "p_high": float(np.mean(ctrl >= n_obs)),
            "power": power,
        }
    p_min = min(v["p_low"] for v in out.values())
    return {"radii": out, "p_low_min": p_min, "rule_outcome": decide(p_min)}


def decide(p_low_min: float) -> str:
    """PLAN.md decision rule on the smallest one-sided deficit p over the radii (ASSUMPTIONs)."""
    if p_low_min < 0.01:
        return "vet the centres"
    return "exploratory note" if p_low_min < 0.05 else "null"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fields", nargs="+", default=list(wsr.RADIAL), choices=list(wsr.RADIAL))
    ap.add_argument("--n-draws", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "w5_radial_counts")
    ap.add_argument("--work", type=Path, default=None, help="Gaia / zout cache (default OUT)")
    ap.add_argument("--radii", nargs="+", type=float, default=list(RADII), help="arcsec")
    ap.add_argument(
        "--exclude-near-arcsec", type=float, default=0.0,
        help="also drop rows within this of a radial arc (post hoc diagnostic)",
    )  # fmt: skip
    args = ap.parse_args(argv)
    work = args.work or args.out
    rng = np.random.default_rng(args.seed)
    fields = [
        field_counts(n, work, args.n_draws, rng, tuple(args.radii), args.exclude_near_arcsec)
        for n in args.fields
    ]
    summary = {"plan": "results/w5_radial_counts/PLAN.md", "provenance": "derived",
               "argv": sys.argv[1:] if argv is None else argv,
               "fields": fields, "stack": combine(fields)}  # fmt: skip
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["stack"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
