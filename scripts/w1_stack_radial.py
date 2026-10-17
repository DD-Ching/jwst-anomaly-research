"""W1 shear statistic S stacked at the ``radial`` screen's line-convergence centres.

Plan and decision rule: ``results/w1_stack_radial/PLAN.md``. Each field is loaded as in
``inject_shear.run_field`` (``inject_radial.load_field`` on ``inject_shear.field_spec``); the used
ellipticities come from ``exotic_screens.shear_sources`` with the shear-screen defaults, and S, S_×
are evaluated with ``exotic_screens.ApertureMass`` at the centres of ``radial.json``. The null
rotates every used ellipticity (positions and |e| kept); one rotation per draw and field.
Sign check: S with the cluster shear kept (uncorrected ε), averaged over centres within
``SIGN_RADIUS`` of the model origin, must be negative (tangential cluster shear).

    python scripts/w1_stack_radial.py --out results/w1_stack_radial [--n-random 2000]

All outputs are ``derived``; thresholds are ASSUMPTIONs (PLAN.md).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import exotic_screens as es  # noqa: E402
import inject_shear as ish  # noqa: E402

from jwst_anomaly import paths  # noqa: E402

ir = ish.ir
REPO = Path(__file__).resolve().parents[1]
RADIAL = {
    "abells1063": REPO / "results/radial_more/abells1063_reproduce/radial.json",
    "macs0417": REPO / "results/radial_more/macs0417/radial.json",
    "macs1423": REPO / "results/radial_more/macs1423/radial.json",
    "whl0137": REPO / "results/radial_more/whl0137/radial.json",
}
SIGN_RADIUS = 30.0  # arcsec, centres used for the sign check (ASSUMPTION)


def stack_null(
    ap: es.ApertureMass, e: np.ndarray, n_draws: int, rng
) -> tuple[np.ndarray, np.ndarray]:
    """Per-draw sum and count of finite S over ``ap``'s centres for random rotations of ``e``."""
    rot = np.exp(1j * rng.uniform(0, 2 * np.pi, (len(e), n_draws)))
    s, _ = ap.snr(e[:, None] * rot)
    return np.nansum(s, axis=0), np.isfinite(s).sum(axis=0)


def radial_arc_rows(model, shapes, spec: dict, extra) -> np.ndarray:
    """Bool per shape row: the ``radial`` screen's selected arcs
    (``exotic_screens.radial_candidates``
    with the ``radial_more`` settings), matched by catalogue label."""
    rargs = argparse.Namespace(**es.radial_defaults())
    rargs.max_radius = spec["max_radius"]
    cand, _ = es.radial_candidates(model, shapes, rargs, extra)
    return np.isin(np.asarray(shapes["label"]), np.asarray(cand["label"]))


def near_rows(model, shapes, rows: np.ndarray, radius: float) -> np.ndarray:
    """Bool per shape row: within ``radius`` arcsec of any row in ``rows`` (``rows`` included)."""
    from scipy.spatial import cKDTree

    x, y = model.to_frame(np.asarray(shapes["ra"]), np.asarray(shapes["dec"]))
    if not rows.any():
        return rows.copy()
    d, _ = cKDTree(np.c_[x[rows], y[rows]]).query(np.c_[x, y])
    return d <= radius


def neighbour_alignment(model, shapes, arcs: np.ndarray, nb: np.ndarray) -> dict:
    """Mean cos 2(PA_neighbour - PA_nearest arc) over used neighbour rows (1 parallel, 0 random)."""
    from scipy.spatial import cKDTree

    if not (arcs.any() and nb.any()):
        return {"n": int(nb.sum()), "mean_cos2dpa": None}
    x, y = model.to_frame(np.asarray(shapes["ra"]), np.asarray(shapes["dec"]))
    pa = np.radians(np.asarray(shapes["pa_obs"], float))
    _, k = cKDTree(np.c_[x[arcs], y[arcs]]).query(np.c_[x[nb], y[nb]])
    c = np.cos(2 * (pa[nb] - pa[arcs][k]))
    return {
        "n": int(nb.sum()),
        "mean_cos2dpa": float(np.mean(c)),
        "sem": float(np.std(c) / np.sqrt(len(c))),
    }


def field_stack(
    name: str,
    out: Path,
    n_random: int,
    seed,
    exclude_radial: bool = False,
    near: float = 0.0,
    shear_model: str | None = None,
) -> dict:
    spec = ish.field_spec(name)
    radial_model = spec["model"]
    if shear_model is not None:  # a second model removes the cluster shear; centres stay fixed
        spec = {**spec, "model": shear_model}
    (out / name).mkdir(parents=True, exist_ok=True)
    model, shapes, extra, _, _ = ir.load_field(name, out / name, spec)
    sargs = argparse.Namespace(**ish.SHEAR_DEFAULTS)
    psf = es.psf_sigma_px(shapes)
    e, counts = es.shear_sources(model, shapes, psf, sargs.min_snr, sargs.max_g, None, extra)
    n_radial_used = None
    if exclude_radial:  # vetting step: the arcs that define the centres leave the shear sample
        rmodel = model
        if radial_model != spec["model"]:  # the arc selection is the radial screen's own model
            rmodel, _, _ = ir.lc.load_model(radial_model)
            ir.lc.apply_frame_offset(radial_model, rmodel)
        arcs = radial_arc_rows(rmodel, shapes, spec, extra)
        if near > 0:  # also the arcs' neighbours, e.g. deblended fragments of the same arc
            nb = near_rows(model, shapes, arcs, near) & ~arcs
            align = neighbour_alignment(model, shapes, arcs, nb & np.isfinite(e))
            arcs = arcs | nb
        n_radial_used = int((arcs & np.isfinite(e)).sum())
        e = np.where(arcs, np.nan + 0j, e)
    use = np.isfinite(e)
    x, y = model.to_frame(np.asarray(shapes["ra"])[use], np.asarray(shapes["dec"])[use])
    peaks = json.loads(RADIAL[name].read_text())["peaks"]
    if json.loads(RADIAL[name].read_text())["model"] != radial_model:
        raise ValueError(f"{name}: radial.json model differs from the shear field spec")
    cx, cy = model.to_frame(np.array([p["ra"] for p in peaks]), np.array([p["dec"] for p in peaks]))
    kw = dict(radius=sargs.aperture_arcsec, kind=sargs.filter, r_min=sargs.r_min_arcsec)
    ap = es.ApertureMass(x, y, cx, cy, min_n=sargs.min_sources, **kw)
    s, sx = ap.snr(e[use])
    null_sum, null_n = stack_null(ap, e[use], n_random, np.random.default_rng(seed))
    # sign check: cluster shear kept
    gx, gy = es.radial_grid(SIGN_RADIUS, sargs.grid_arcsec)
    inner = np.hypot(gx, gy) <= SIGN_RADIUS
    ap0 = es.ApertureMass(x, y, gx[inner], gy[inner], min_n=sargs.min_sources, **kw)
    s_keep, _ = ap0.snr(counts["eps"][use])
    s_corr, _ = ap0.snr(e[use])
    return {
        "field": name,
        "model": spec["model"],
        "n_used": int(use.sum()),
        "n_rows_removed": n_radial_used,  # radial arcs (+ neighbours with --exclude-near-arcsec)
        "neighbour_alignment": align if exclude_radial and near > 0 else None,
        "responsivity": counts["responsivity"],
        "centres": [
            {"ra": p["ra"], "dec": p["dec"], "n_lines": p["n_lines"], "s": _f(a), "s_cross": _f(b)}
            for p, a, b in zip(peaks, s, sx, strict=True)
        ],
        "null_sum": null_sum,
        "null_n": null_n,
        "sign_check": {
            "n_centres": int(np.isfinite(s_keep).sum()),
            "mean_s_cluster_kept": float(np.nanmean(s_keep)),
            "mean_s_cluster_removed": float(np.nanmean(s_corr)),
        },
    }


def decide(p: float, sign_ok: bool) -> str:
    """PLAN.md decision rule on the pooled one-sided p (ASSUMPTION thresholds); without a passed
    sign check the stack gives no result."""
    if not sign_ok:
        return "no result (sign check failed)"
    if p < 0.01:
        return "vet the centres"
    if p < 0.05:
        return "exploratory note"
    return "null"


def _f(v) -> float | None:
    return float(v) if np.isfinite(v) else None


def combine(fields: list[dict]) -> dict:
    s = np.array([c["s"] for f in fields for c in f["centres"] if c["s"] is not None])
    sx = np.array([c["s_cross"] for f in fields for c in f["centres"] if c["s_cross"] is not None])
    null_mean = sum(f["null_sum"] for f in fields) / np.maximum(sum(f["null_n"] for f in fields), 1)
    obs = float(s.mean()) if len(s) else float("nan")
    return {
        "n_centres": int(len(s)),
        "mean_s": obs,
        "mean_s_cross": float(sx.mean()) if len(sx) else float("nan"),
        "p_one_sided": float((null_mean >= obs).mean()),
        "null_mean_p50_p95": [
            float(np.percentile(null_mean, 50)),
            float(np.percentile(null_mean, 95)),
        ],
        "n_random": int(len(null_mean)),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fields", nargs="+", default=list(RADIAL), choices=list(RADIAL))
    ap.add_argument("--n-random", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=314)
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "w1_stack_radial")
    ap.add_argument(
        "--exclude-radial-arcs",
        action="store_true",
        help="vetting step: drop the radial screen's selected arcs from the shear sample",
    )
    ap.add_argument(
        "--exclude-near-arcsec",
        type=float,
        default=0.0,
        help="with --exclude-radial-arcs, also drop sources within this distance of those arcs",
    )
    ap.add_argument(
        "--shear-model",
        default=None,
        help="lens model that removes the cluster shear (one field only; centres stay fixed)",
    )
    ap.add_argument("--work", type=Path, default=None, help="Gaia star / zout cache (default OUT)")
    args = ap.parse_args(argv)
    if args.shear_model is not None and len(args.fields) != 1:
        ap.error("--shear-model needs exactly one field")
    work = args.work or args.out
    seeds = np.random.SeedSequence(args.seed).spawn(len(args.fields))
    fields = [
        field_stack(
            n,
            work,
            args.n_random,
            sd,
            args.exclude_radial_arcs,
            args.exclude_near_arcsec,
            args.shear_model,
        )
        for n, sd in zip(args.fields, seeds, strict=True)
    ]
    summary = {
        "provenance": "derived",
        "plan": "results/w1_stack_radial/PLAN.md",
        "seed": args.seed,
        "assumptions": {
            **{
                k: ish.SHEAR_DEFAULTS[k]
                for k in (
                    "aperture_arcsec",
                    "filter",
                    "r_min_arcsec",
                    "min_sources",
                    "min_snr",
                    "max_g",
                )
            },
            "sign_radius_arcsec": SIGN_RADIUS,
            "streak_veto": False,
            "exclude_radial_arcs": args.exclude_radial_arcs,
            "exclude_near_arcsec": args.exclude_near_arcsec,
            "shear_model": args.shear_model,
        },
        "stack": combine(fields),
        "fields": [{k: v for k, v in f.items() if k not in ("null_sum", "null_n")} for f in fields],
    }
    sign_ok = all(f["sign_check"]["mean_s_cluster_kept"] < 0 for f in fields)
    summary["sign_check_passed"] = sign_ok
    summary["outcome"] = decide(summary["stack"]["p_one_sided"], sign_ok)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["stack"], indent=1), summary["outcome"])
    for f in summary["fields"]:
        print(f["field"], f["n_used"], f["sign_check"], [c["s"] for c in f["centres"]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
