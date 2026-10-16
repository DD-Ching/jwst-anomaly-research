"""Vetting step 4 for Abell S1063 ``arcs``: clump-pair census of the 24 strong-shear rows.

Frozen plan: ``results/arcs/abells1063-cats/vetting/PLAN_census.md`` (written before the run;
the sample test is post hoc, so its result is exploratory).

Each row of ``bgsub_shapes.ecsv`` gets the step-3 deblend (``arcs_vet_deblend.fit_models`` and
``deblend_class``) and its own-flux gate. ``two-knot`` rows get the clump-orientation anti
probability, the others the screen's shape-noise probability; the Poisson-binomial P(≥ k) of the
observed anti count replaces the screen's shape-noise P.

    python scripts/arcs_vet_census.py [--workers 4]

Output: ``results/arcs/abells1063-cats/vetting/census.{json,ecsv}`` (derived; pixels observed,
injections simulations, probabilities model_prediction).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")  # pool workers: no BLAS oversubscription

import argparse
import importlib.util
import json
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy import stats
from scipy.ndimage import gaussian_filter

from jwst_anomaly import cutouts, paths


def _load_module(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


vd = _load_module("arcs_vet_deblend", "arcs_vet_deblend.py")
vc = _load_module("arcs_vet_clumps", "arcs_vet_clumps.py")
vn = _load_module("arcs_vet_noise", "arcs_vet_noise.py")
lc = _load_module("lens_consistency", "lens_consistency.py")


# ---- frozen rules ---------------------------------------------------------------------------
def gate_separation(sep_arcsec: float) -> float:
    """Pair separation for the gate: the observed one clipped to 0.10–0.30″, rounded to 0.01″."""
    return round(float(np.clip(sep_arcsec, 0.10, 0.30)), 2)


def row_p_anti(cls: str, gate_passed: bool, p_cp: float, p_sn: float, f_hat: float) -> float:
    """Per-row anti probability (PLAN_census.md)."""
    if not gate_passed:
        return f_hat * p_cp + (1.0 - f_hat) * p_sn
    return p_cp if cls == "two-knot" else p_sn


def p_tail(probs: dict, labels: list, k: int, drop: tuple = ()) -> float:
    """Poisson-binomial P(≥ k) over ``labels``; a row removed from the count is also removed from
    the null (PLAN_fixes.md fix 4)."""
    return vn.excess_p([probs[lab] for lab in labels if lab not in drop], k)


def clopper_pearson_upper(k: int, n: int, cl: float = 0.95) -> float:
    return 1.0 if k >= n else float(stats.beta.ppf(1 - (1 - cl) / 2, k + 1, n - k))


def outcome(p_primary: float) -> str:
    if not np.isfinite(p_primary):
        return "undetermined (fewer than 12 gate-passed rows)"
    if p_primary >= 0.05:
        return "anti excess explained by clump orientation (exploratory)"
    return "excess remains under the clump-aware null"


# ---- per row --------------------------------------------------------------------------------
def deblend_row(path: str, ra: float, dec: float, fit_arcsec: float, n_inject: int, seed: int):
    sci, sub, err, wcs, px = vd.load(path)
    xt, yt = wcs.world_to_pixel_values(ra, dec)
    r_px = fit_arcsec / px
    yy, xx = np.indices(sub.shape, dtype=float)
    near = np.hypot(xx - xt, yy - yt) <= 3
    j, i = np.unravel_index(
        np.argmax(np.where(near, gaussian_filter(np.nan_to_num(sub), 1.0), -np.inf)), sub.shape
    )
    obs = vd.fit_models(sub, err, float(i), float(j), r_px)
    m1, m2 = obs["m1"], obs["m2"]
    sep = obs["sep_px"] * px
    cls = vd.deblend_class(obs["dbic"], obs["snr1"], obs["snr2"], sep)
    sep_g = gate_separation(sep)
    sizes = {"single": (m1[2], m1[3], m1[4]), "two": (m2[2], m2[3], m2[6], m2[7])}
    g = vd.gate(
        sci, err, xt, yt, px, r_px, sizes, n_inject, np.random.default_rng(seed), sep_arcsec=sep_g
    )
    pa_sep = vd.sky_pa(wcs, m2[0], m2[1], m2[4], m2[5])
    return {
        "dbic": obs["dbic"],
        "snr1": obs["snr1"],
        "snr2": obs["snr2"],
        "sep_arcsec": sep,
        "pa_sep": pa_sep,
        "m1_sigma_arcsec": m1[3] * px,
        "m1_q": m1[4],
        "class": cls,
        "gate_sep_arcsec": sep_g,
        "gate": g,
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model", default="abells1063-cats")
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    p.add_argument("--fit-arcsec", type=float, default=0.4, help="ASSUMPTION")
    p.add_argument("--n-inject", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261016)
    p.add_argument("--sigma-eps", type=float, default=0.26, help="ASSUMPTION (screen value)")
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--cut-dir", type=Path, default=None, help="reuse cutouts (else a temp dir)")
    a = p.parse_args(argv)
    d = a.arcs_dir / a.model
    rows = Table.read(d / "bgsub_shapes.ecsv")
    full = Table.read(d / "arcs.ecsv")
    gmin = {str(r["label"]): float(r["g_min"]) for r in full}
    screen_anti = {str(r["label"]) for r in rows if r["orientation_class"] == "anti"}
    bgsub_anti = {str(r["label"]) for r in rows if r["class_bgsub"] == "anti"}
    tg = Table(
        {"source_uid": [str(v) for v in rows["label"]], "ra": rows["ra"], "dec": rows["dec"]}
    )
    labels = list(tg["source_uid"])

    with tempfile.TemporaryDirectory() as tmp:
        base = a.cut_dir or Path(tmp)
        c = cutouts.make_cutouts(
            vd.URI.format(b="f200w"),
            tg,
            size_arcsec=8.0,
            out_dir=base / "f200w",
            extensions=["ERR", "WHT"],
        )
        path = {str(u): str(pth) for u, pth in zip(c["source_uid"], c["path"], strict=True)}
        outside = sorted(lab for lab, pth in path.items() if not pth)
        if outside:
            raise SystemExit(f"error: rows outside the F200W image: {outside}")
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = {
                lab: ex.submit(
                    deblend_row,
                    path[lab],
                    float(r["ra"]),
                    float(r["dec"]),
                    a.fit_arcsec,
                    a.n_inject,
                    a.seed + k,
                )
                for k, (lab, r) in enumerate(zip(labels, rows, strict=True))
            }
            res = {lab: f.result() for lab, f in futs.items()}

    passed = [lab for lab in labels if res[lab]["gate"]["passed"]]
    n_two = sum(res[lab]["class"] == "two-knot" for lab in passed)
    f_hat = n_two / len(passed) if passed else float("nan")
    f_up = clopper_pearson_upper(n_two, len(passed)) if passed else float("nan")
    p_cp = {lab: vc.clump_anti_fraction(gmin[lab]) for lab in labels}
    p_sn = {lab: lc.anti_probability(gmin[lab], a.sigma_eps, 0.5, 60.0) for lab in labels}
    p_row = {
        lab: row_p_anti(res[lab]["class"], res[lab]["gate"]["passed"], p_cp[lab], p_sn[lab], f_hat)
        for lab in labels
    }

    def tail(pr: dict, k: int, drop: tuple[str, ...] = ()) -> float:
        return p_tail(pr, labels, k, drop) if len(passed) >= 12 else float("nan")

    mix = {f: {lab: f * p_cp[lab] + (1 - f) * p_sn[lab] for lab in labels} for f in (f_hat, f_up)}
    # 2184: noise-explicable at core flux (step 2); counts and nulls are derived from the same sets
    k_counts = {"k4_screen": len(screen_anti), "k3_bgsub": len(bgsub_anti)}
    drops = {}
    if "2184" in bgsub_anti:  # the core-flux count exists only when the step-2 exclusion applies
        drops["k2_core_flux"] = ("2184",)
        k_counts["k2_core_flux"] = len(bgsub_anti - {"2184"})
    p_primary = tail(p_row, 4)
    sample = {
        "n_rows": len(labels),
        "n_gate_passed": len(passed),
        "n_two_knot_gate_passed": n_two,
        "f_two_knot": round(f_hat, 3),
        "f_two_knot_95_upper": round(f_up, 3),
        "expected_anti_shape_noise_only": round(sum(p_sn.values()), 3),
        "expected_anti_clump_aware": round(sum(p_row.values()), 3),
        "p_ge_k_clump_aware": {n: tail(p_row, k, drops.get(n, ())) for n, k in k_counts.items()},
        "p_ge_k_mixture_f_hat": {
            n: tail(mix[f_hat], k, drops.get(n, ())) for n, k in k_counts.items()
        },
        "p_ge_k_mixture_f_upper": {
            n: tail(mix[f_up], k, drops.get(n, ())) for n, k in k_counts.items()
        },
        "p_ge_k_shape_noise_only_check": {
            n: p_tail(p_sn, labels, k, drops.get(n, ())) for n, k in k_counts.items()
        },
        "k_counts": k_counts,
        "rows_dropped_from_count_and_null": {n: list(v) for n, v in drops.items()},
        "two_knot_among_screen_aligned": sum(
            res[lab]["class"] == "two-knot" and res[lab]["gate"]["passed"]
            for lab in labels
            if lab not in screen_anti
        ),
        "two_knot_among_screen_anti": sorted(
            lab for lab in screen_anti if res[lab]["class"] == "two-knot" and lab in passed
        ),
        "outcome": outcome(p_primary),
    }
    t = Table(
        {
            "label": labels,
            "g_min": [gmin[lab] for lab in labels],
            "screen_class": [str(r["orientation_class"]) for r in rows],
            "bgsub_class": [str(r["class_bgsub"]) for r in rows],
            "mag": rows["mag"],
            "dbic": [res[lab]["dbic"] for lab in labels],
            "snr_min": [min(res[lab]["snr1"], res[lab]["snr2"]) for lab in labels],
            "sep_arcsec": [res[lab]["sep_arcsec"] for lab in labels],
            "deblend_class": [res[lab]["class"] for lab in labels],
            "gate_single_ok": [res[lab]["gate"]["frac_single_not_two_knot"] for lab in labels],
            "gate_two_ok": [res[lab]["gate"]["frac_two_called_two_knot"] for lab in labels],
            "gate_passed": [res[lab]["gate"]["passed"] for lab in labels],
            "p_cp": [p_cp[lab] for lab in labels],
            "p_sn": [p_sn[lab] for lab in labels],
            "p_anti_row": [p_row[lab] for lab in labels],
        },
        meta={
            "provenance": "derived (exploratory, post hoc sample test; PLAN_census.md)",
            "image": vd.URI.format(b="f200w"),
            "screen_anti": sorted(screen_anti),
            "bgsub_anti": sorted(bgsub_anti),
        },
    )
    out_dir = d / "vetting"
    t.write(out_dir / "census.ecsv", overwrite=True)
    out = {
        "plan": "PLAN_census.md",
        "label": "exploratory; fits derived, injections simulated, probabilities model_prediction",
        "assumptions": {k: getattr(a, k) for k in ("fit_arcsec", "n_inject", "seed", "sigma_eps")}
        | {"filter_arcsec": 1.0, "min_ellipticity": 0.5, "anti_deg": 60.0},
        "sample": sample,
        "rows": {lab: res[lab] for lab in labels},
    }
    (out_dir / "census.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(sample, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
