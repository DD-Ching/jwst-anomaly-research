"""W1 lenses injected at the ``radial`` centres through the W1 shear stack.

Plan and rules: ``results/w1_stack_radial/inject/PLAN.md``. The sample is the stack's vet step 1
(``w1_stack_radial.py --exclude-radial-arcs``): the ``radial`` screen's selected arcs leave the
shear sample and are not painted. One lens of mass M is painted at each valid centre in turn
(``inject_shear.ShearInjector.lensed_sample``: the D-050/D-053 chain) and S is evaluated at that
centre only. A fraction f of the centres then takes its injected S, the others keep the real S,
and the pooled stacked mean is compared with the observed mean and with the real rotation null.

    python scripts/w1_stack_inject.py --out results/w1_stack_radial/inject

Injected shapes are ``simulated``; S values, stacks and rates are ``derived``; thresholds are
ASSUMPTIONs (PLAN.md).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import inject_shear as ish  # noqa: E402
import w1_stack_radial as wsr  # noqa: E402

from jwst_anomaly import paths  # noqa: E402

ir = ish.ir
es = ish.es
FIELDS = ("abells1063", "macs0417")  # the fields with valid centres in the stack (PLAN.md)
FRACTIONS = (0.1, 0.25, 0.5, 1.0)
EXCLUDE_LEVEL = 0.95  # share of draws above the observed mean that excludes (M, f) (ASSUMPTION)
NULL_PCT = 95.0  # the real rotation null's percentile for the efficiency (ASSUMPTION)
GATE_Z_TOL = 1e-9


def excluded(stacks: np.ndarray, observed: float, level: float = EXCLUDE_LEVEL) -> bool:
    """PLAN rule: (M, f) is excluded when at least ``level`` of the draws' stacked means exceed
    the observed stacked mean."""
    return bool(np.mean(np.asarray(stacks) > observed) >= level)


def efficiency(stacks: np.ndarray, null_p95: float) -> float:
    """Share of draws whose stacked mean exceeds the real null's 95th percentile."""
    return float(np.mean(np.asarray(stacks) > null_p95))


def centre_null(ap: es.ApertureMass, e: np.ndarray, n_draws: int, rng) -> np.ndarray:
    """S per centre (rows) and random rotation of ``e`` (columns); the draws of
    ``w1_stack_radial.stack_null``, kept per centre."""
    rot = np.exp(1j * rng.uniform(0, 2 * np.pi, (len(e), n_draws)))
    return ap.snr(e[:, None] * rot)[0]


def noise_stacks(s_null: np.ndarray, delta: np.ndarray, f: float, rng) -> np.ndarray:
    """Post hoc (RESULT.md): stacked means when ``max(1, round(f N))`` random centres add their
    injected change ``delta`` = S_inj - S_obs to one rotation-null realisation ``s_null[:, i]``
    (one subset per draw). A centre whose injected S is undefined (NaN delta) drops out of that
    draw's mean, as in the real stack."""
    n, n_draws = s_null.shape
    k = max(1, int(round(f * n)))
    out = np.empty(n_draws)
    for i in range(n_draws):
        s = s_null[:, i].copy()
        pick = np.arange(n) if k >= n else rng.choice(n, k, replace=False)
        s[pick] = s[pick] + delta[pick]
        out[i] = np.nanmean(s)
    return out


def fraction_stacks(
    s_obs: np.ndarray, s_inj: np.ndarray, f: float, n_draws: int, rng
) -> np.ndarray:
    """Stacked means when ``max(1, round(f N))`` random centres take ``s_inj``, the rest
    ``s_obs`` (one draw when all centres are taken). A centre whose injected S is undefined (NaN:
    too few sources left in its aperture) drops out of the mean, as in the real stack (post hoc
    detail; PLAN.md does not say)."""
    n = len(s_obs)
    k = max(1, int(round(f * n)))
    if k >= n:
        return np.array([float(np.nanmean(s_inj))])
    out = np.empty(n_draws)
    for i in range(n_draws):
        pick = rng.choice(n, k, replace=False)
        s = s_obs.copy()
        s[pick] = s_inj[pick]
        out[i] = np.nanmean(s)
    return out


def field_injector(name: str, out: Path):
    """The vet-step-1 shear sample of one field as a ``ShearInjector`` plus its valid centres."""
    spec = ish.field_spec(name)
    (out / name).mkdir(parents=True, exist_ok=True)
    model, shapes, extra, cat, _ = ir.load_field(name, out / name, spec)
    sargs = SimpleNamespace(**ish.SHEAR_DEFAULTS)
    sargs.max_radius = spec["max_radius"]
    sargs.recover_tol = 0.0
    sargs.min_coverage = 0.0  # the stack applies no coverage cut
    psf = es.psf_sigma_px(shapes)
    inj = ish.ShearInjector(model, shapes, sargs, psf, extra)
    arcs = wsr.radial_arc_rows(model, shapes, spec, extra)
    inj.e = np.where(arcs, np.nan + 0j, inj.e)
    inj.use = inj.use & ~arcs
    inj.background = inj.background & ~arcs
    peaks = json.loads(wsr.RADIAL[name].read_text())["peaks"]
    cx, cy = model.to_frame(np.array([p["ra"] for p in peaks]), np.array([p["dec"] for p in peaks]))
    return inj, np.atleast_1d(cx), np.atleast_1d(cy), ir.pixel_scale_arcsec(cat), shapes


def s_at(inj, x, y, e, cx: float, cy: float) -> float:
    a = inj.args
    ap = es.ApertureMass(
        x, y, np.array([cx]), np.array([cy]), a.aperture_arcsec, a.filter, a.r_min_arcsec,
        a.min_sources,
    )  # fmt: skip
    s, _ = ap.snr(e)
    return float(s[0])


def field_run(name: str, out: Path, masses, n_random: int, seed) -> dict:
    inj, cx, cy, scale, shapes = field_injector(name, out)
    xs, ys, e = inj.xs[inj.use], inj.ys[inj.use], inj.e[inj.use]
    s_obs = np.array([s_at(inj, xs, ys, e, a, b) for a, b in zip(cx, cy, strict=True)])
    valid = np.isfinite(s_obs)
    a = inj.args
    ap = es.ApertureMass(xs, ys, cx, cy, a.aperture_arcsec, a.filter, a.r_min_arcsec, a.min_sources)
    s_null = centre_null(ap, e, n_random, np.random.default_rng(seed))
    zero = np.zeros(len(shapes))  # gate Z: theta_E = 0 paints nothing
    s_zero = np.array(
        [
            s_at(inj, *inj.lensed_sample(x, y, zero, scale)[:3], x, y)
            for x, y in zip(cx, cy, strict=True)
        ]
    )
    z_src = ir.source_redshift(shapes)
    rows = []
    for m in masses:
        theta = np.where(inj.background, ir.theta_e_arcsec(m, inj.model.z_lens, z_src), np.nan)
        for i in np.flatnonzero(valid):
            x, y, ee, info, n_img = inj.lensed_sample(cx[i], cy[i], theta, scale)
            rows.append(
                {"mass_msun": m, "centre": int(i), "s_inj": s_at(inj, x, y, ee, cx[i], cy[i]),
                 "n_lensed": info["n_lensed"], "n_image_used": n_img}
            )  # fmt: skip
    return {
        "field": name,
        "model": inj.model.name if hasattr(inj.model, "name") else ish.field_spec(name)["model"],
        "n_used": int(inj.use.sum()),
        "s_obs": [wsr._f(v) for v in s_obs],
        "gate_z_max_abs_diff": float(np.nanmax(np.abs(s_zero[valid] - s_obs[valid]))),
        "gate_z_valid_same": bool(np.array_equal(np.isfinite(s_zero), valid)),
        "injections": rows,
        "s_null": s_null,
    }


def combine(fields: list[dict], masses, n_draws: int, rng) -> dict:
    s_obs = np.array([v for f in fields for v in f["s_obs"] if v is not None])
    obs = float(s_obs.mean())
    s_null = np.vstack([f["s_null"][np.array([v is not None for v in f["s_obs"]])] for f in fields])
    null_mean = s_null.mean(axis=0)
    p95 = float(np.percentile(null_mean, NULL_PCT))
    gate_z = all(f["gate_z_max_abs_diff"] <= GATE_Z_TOL and f["gate_z_valid_same"] for f in fields)
    cells = []
    for m in masses:
        s_inj = np.array(
            [r["s_inj"] for f in fields for r in f["injections"] if r["mass_msun"] == m]
        )
        delta = s_inj - s_obs
        for fr in FRACTIONS:
            st = fraction_stacks(s_obs, s_inj, fr, n_draws, rng)
            nz = noise_stacks(s_null, delta, fr, rng)
            cells.append(
                {
                    "mass_msun": m,
                    "fraction": fr,
                    "n_draws": len(st),
                    "median_stack": float(np.nanmedian(st)),
                    "share_above_observed": float(np.mean(st > obs)),
                    "excluded": excluded(st, obs),
                    "efficiency": efficiency(st, p95),
                    "n_centres_undefined": int(np.isnan(s_inj).sum()),
                    "noise_aware": {  # post hoc (RESULT.md)
                        "n_draws": len(nz),
                        "median_stack": float(np.median(nz)),
                        "share_above_observed": float(np.mean(nz > obs)),
                        "excluded": excluded(nz, obs),
                        "efficiency": efficiency(nz, p95),
                    },
                }
            )
    top = [c for c in cells if c["mass_msun"] == max(masses) and c["fraction"] == 1.0]
    gate_r = bool(top and top[0]["median_stack"] > p95)
    return {
        "n_centres": int(len(s_obs)),
        "observed_mean_s": obs,
        "null_mean_p50_p95": [float(np.percentile(null_mean, 50)), p95],
        "p_one_sided": float(np.mean(null_mean >= obs)),
        "gate_z_passed": gate_z,
        "gate_r_passed": gate_r,
        "outcome": "limit"
        if gate_z and gate_r
        else ("no limit (gate Z)" if not gate_z else "no reach"),
        "cells": cells,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fields", nargs="+", default=list(FIELDS), choices=list(wsr.RADIAL))
    ap.add_argument("--mass", nargs="+", type=float, default=list(ir.MASSES))
    ap.add_argument("--n-random", type=int, default=2000)
    ap.add_argument("--n-draws", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=316)
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "w1_stack_inject")
    ap.add_argument("--work", type=Path, default=None, help="Gaia star / zout cache (default OUT)")
    args = ap.parse_args(argv)
    work = args.work or args.out
    ss = np.random.SeedSequence(args.seed)
    seeds = ss.spawn(len(args.fields) + 1)
    fields = [
        field_run(n, work, args.mass, args.n_random, sd)
        for n, sd in zip(args.fields, seeds[:-1], strict=True)
    ]
    stack = combine(fields, args.mass, args.n_draws, np.random.default_rng(seeds[-1]))
    summary = {
        "provenance": "derived (injected shapes simulated)",
        "plan": "results/w1_stack_radial/inject/PLAN.md",
        "seed": args.seed,
        "assumptions": {
            "exclude_level": EXCLUDE_LEVEL,
            "null_percentile": NULL_PCT,
            "fractions": list(FRACTIONS),
            "one_lens_changes_only_its_centre": True,
            "exclude_radial_arcs": True,
            "screen": {
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
        },  # fmt: skip
        "stack": stack,
        "fields": [{k: v for k, v in f.items() if k != "s_null"} for f in fields],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(es.json_safe(summary), indent=1) + "\n")
    print(json.dumps({k: v for k, v in stack.items() if k != "cells"}, indent=1))
    for c in stack["cells"]:
        print(c)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
