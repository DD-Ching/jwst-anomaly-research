"""MACS1423 ``arcs`` vetting step 2: clump-pair census of the 55 strong-shear background rows.

Frozen plan: ``results/arcs/macs1423-canucs/vetting/PLAN_census.md`` (written before the run).
Reuses #273's ``arcs_vet_census.py`` unchanged at commit ``43fbbbd``, extracted to a directory:

    mkdir -p S && for f in arcs_vet_census arcs_vet_deblend arcs_vet_clumps arcs_vet_noise \\
        lens_consistency; do git show 43fbbbd:scripts/$f.py > S/$f.py; done
    python scripts/arcs_macs1423_census.py --s1063-scripts S [--workers 4]

Only the image and the anti sets (photo-z basis) differ from the S1063 run. Output:
``results/arcs/macs1423-canucs/vetting/census.{json,ecsv}`` (derived; pixels observed, injections
simulated, probabilities model_prediction).
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

from astropy.table import Table

from jwst_anomaly import cutouts, paths

IMAGE = "s3://stpubdata/jwst/public/jw01208/L3/t/o010/jw01208-o010_t005_nircam_clear-f200w_i2d.fits"
MIN_GATE_PASS_FRAC = 0.5  # ASSUMPTION (PLAN_census.md: below it the census is undetermined here)


def load_census(s1063_dir: Path):
    """Import #273's census module from ``s1063_dir`` (its siblings load from the same place)."""
    sys.path.insert(0, str(s1063_dir))  # spawned pool workers re-import by module name
    spec = importlib.util.spec_from_file_location(
        "arcs_vet_census", s1063_dir / "arcs_vet_census.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["arcs_vet_census"] = mod
    spec.loader.exec_module(mod)
    return mod


def anti_sets(rows: Table) -> dict[str, list[str]]:
    """Primary: screen anti (photo-z basis); secondary: those still anti after ICL subtraction."""
    screen = sorted(str(r["label"]) for r in rows if r["orientation_class"] == "anti")
    staying = sorted(
        str(r["label"])
        for r in rows
        if r["orientation_class"] == "anti" and r["class_bgsub"] == "anti"
    )
    return {"screen_anti": screen, "screen_anti_staying_bgsub": staying}


def final_outcome(census_outcome: str, n_passed: int, n_rows: int) -> str:
    """#273's outcome, unless the gate fails on most rows of this image (PLAN_census.md)."""
    if n_rows == 0 or n_passed / n_rows < MIN_GATE_PASS_FRAC:
        return "undetermined (deblend gate fails on most rows of this image)"
    return census_outcome


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--s1063-scripts", type=Path, required=True)
    p.add_argument("--model", default="macs1423-canucs-maps")
    p.add_argument(
        "--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs" / "macs1423-canucs"
    )
    p.add_argument("--fit-arcsec", type=float, default=0.4, help="ASSUMPTION")
    p.add_argument("--n-inject", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261016)
    p.add_argument("--sigma-eps", type=float, default=0.26, help="ASSUMPTION (screen value)")
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--cut-dir", type=Path, default=None, help="reuse cutouts (else a temp dir)")
    a = p.parse_args(argv)
    vc = load_census(a.s1063_scripts.resolve())
    d = a.arcs_dir / a.model
    rows = Table.read(d / "bgsub_shapes.ecsv")
    full = Table.read(d / "arcs.ecsv")
    gmin = {str(r["label"]): float(r["g_min"]) for r in full}
    sets = anti_sets(rows)
    labels = [str(v) for v in rows["label"]]
    tg = Table({"source_uid": labels, "ra": rows["ra"], "dec": rows["dec"]})

    with tempfile.TemporaryDirectory() as tmp:
        base = a.cut_dir or Path(tmp)
        c = cutouts.make_cutouts(
            IMAGE, tg, size_arcsec=8.0, out_dir=base / "f200w", extensions=["ERR", "WHT"]
        )
        path = {str(u): str(pth) for u, pth in zip(c["source_uid"], c["path"], strict=True)}
        outside = sorted(lab for lab, pth in path.items() if not pth)
        if outside:
            raise SystemExit(f"error: rows outside the F200W image: {outside}")
        px_scale = float(vc.vd.load(path[labels[0]])[4])
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = {
                lab: ex.submit(
                    vc.deblend_row,
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
    f_up = vc.clopper_pearson_upper(n_two, len(passed)) if passed else float("nan")
    p_cp = {lab: vc.vc.clump_anti_fraction(gmin[lab]) for lab in labels}
    p_sn = {lab: vc.lc.anti_probability(gmin[lab], a.sigma_eps, 0.5, 60.0) for lab in labels}
    p_row = {
        lab: vc.row_p_anti(
            res[lab]["class"], res[lab]["gate"]["passed"], p_cp[lab], p_sn[lab], f_hat
        )
        for lab in labels
    }

    def tail(pr: dict, k: int) -> float:
        return vc.vn.excess_p([pr[lab] for lab in labels], k) if len(passed) >= 12 else float("nan")

    mix = {f: {lab: f * p_cp[lab] + (1 - f) * p_sn[lab] for lab in labels} for f in (f_hat, f_up)}
    k_counts = {
        "k4_screen": len(sets["screen_anti"]),
        "k3_staying_bgsub": len(sets["screen_anti_staying_bgsub"]),
    }
    p_primary = tail(p_row, k_counts["k4_screen"])
    sample = {
        "n_rows": len(labels),
        "n_gate_passed": len(passed),
        "n_two_knot_gate_passed": n_two,
        "f_two_knot": round(f_hat, 3),
        "f_two_knot_95_upper": round(f_up, 3),
        "expected_anti_shape_noise_only": round(sum(p_sn.values()), 3),
        "expected_anti_clump_aware": round(sum(p_row.values()), 3),
        "p_ge_k_clump_aware": {n: tail(p_row, k) for n, k in k_counts.items()},
        "p_ge_k_mixture_f_hat": {n: tail(mix[f_hat], k) for n, k in k_counts.items()},
        "p_ge_k_mixture_f_upper": {n: tail(mix[f_up], k) for n, k in k_counts.items()},
        "p_ge_k_shape_noise_only_check": {
            n: vc.vn.excess_p(list(p_sn.values()), k) for n, k in k_counts.items()
        },
        "two_knot_among_screen_aligned": sum(
            res[lab]["class"] == "two-knot" and res[lab]["gate"]["passed"]
            for lab in labels
            if lab not in sets["screen_anti"]
        ),
        "classes_of_screen_anti": {
            lab: (res[lab]["class"], bool(res[lab]["gate"]["passed"]))
            for lab in sets["screen_anti"]
        },
        "census_outcome_273": vc.outcome(p_primary),
        "outcome": final_outcome(vc.outcome(p_primary), len(passed), len(labels)),
    }
    t = Table(
        {
            "label": labels,
            "g_min": [gmin[lab] for lab in labels],
            "screen_class": [str(r["orientation_class"]) for r in rows],
            "bgsub_class_z2": [str(r["class_bgsub"]) for r in rows],
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
            "provenance": "derived (exploratory; PLAN_census.md)",
            "image": IMAGE,
            "code": "#273 arcs_vet_census.py at 43fbbbd (unchanged)",
        }
        | sets,
    )
    out_dir = a.arcs_dir / "vetting"
    t.write(out_dir / "census.ecsv", overwrite=True)
    out = {
        "plan": "PLAN_census.md",
        "label": "exploratory; fits derived, injections simulated, probabilities model_prediction",
        "assumptions": {k: getattr(a, k) for k in ("fit_arcsec", "n_inject", "seed", "sigma_eps")}
        | {
            "filter_arcsec": 1.0,
            "min_ellipticity": 0.5,
            "anti_deg": 60.0,
            "min_gate_pass_frac": MIN_GATE_PASS_FRAC,
        },
        "pixel_scale_arcsec": px_scale,
        "anti_sets": sets,
        "sample": sample,
        "rows": {lab: res[lab] for lab in labels},
    }
    (out_dir / "census.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(sample, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
