"""Clump-pair census of the pooled `arcs` anti rows (variant B fields; exploratory).

Plan: ``results/arcs/pooled-census/PLAN.md`` (written before the run). Reuses main's
``arcs_vet_census.py`` (#273's frozen rules) unchanged and ``arcs_pooled.background`` for the rows;
F200W pixels come as 8″ S3 byte-range cutouts and are not stored.

    python scripts/arcs_pooled_census.py [--workers 4] [--fields "RXC J0018"]

Output: ``results/arcs/pooled-census/{census.ecsv,summary.json}`` (pixels observed, fits derived,
injections simulated, probabilities model_prediction).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")  # pool workers: no BLAS oversubscription

import argparse
import json
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

sys.path.insert(0, str(Path(__file__).resolve().parent))  # spawned workers re-import by name
import arcs_pooled as ap  # noqa: E402
import arcs_vet_census as vc  # noqa: E402

from jwst_anomaly import cutouts  # noqa: E402

OUT = ap.ROOT / "pooled-census"
S3 = "s3://stpubdata/jwst/public/{p}/L3/t/{o}/{p}-{o}_{t}_nircam_clear-f200w_i2d.fits"
# field -> F200W mosaic of its screen catalogue (each field's arcs.json / RESULT.md)
IMAGES = {
    "MACS0417": ("jw01208", "o006", "t003"),
    "WHL0137": ("jw02282", "o010", "t001"),
    "RXC J0600": ("jw06796", "o010", "t001"),
    "MACS0308": ("jw05293", "o005", "t001"),
    "RXC J0018": ("jw06882", "o001", "t001"),
}
FIELDS_B = ap.VARIANTS["B"]
# anti rows that stay anti under the second published models (../anti-model-swap/summary.json)
PERSIST = {("WHL0137", "2194"), ("RXC J0600", "581"), ("MACS0308", "2157"), ("MACS0417", "4071")}
MIN_GATE_PASS_FRAC = 0.5  # ASSUMPTION (plan)
MIN_GATE_PASSED = 12  # #273's outcome rule
POOLED_B_EXPECTED = 2.273  # ../pooled-tally/summary.json, variant B at σ_ε 0.26 (reproduced)


def image_url(field: str) -> str:
    p, o, t = IMAGES[field]
    return S3.format(p=p, o=o, t=t)


def load_rows(fields=FIELDS_B) -> Table:
    """Strong-shear background rows of each field's deciding run, with a field column."""
    parts = []
    for f in fields:
        t = ap.background(ap.ROOT / ap.FIELDS[f])
        t = t["label", "ra", "dec", "g_min", "orientation_class", "mag", "z_phot"]
        t["field"] = f
        t["label"] = [str(v) for v in t["label"]]
        parts.append(t)
    return vstack(parts)


def row_probs(gmin, cls, passed, f_hat, sigma_eps=0.26, e_min=0.5, anti_deg=60.0):
    """Per-row (P_cp, P_sn, P_anti) under #273's rules; g_min rounded as arcs_pooled."""
    g = round(float(gmin), 2)
    p_cp = vc.vc.clump_anti_fraction(g, anti_deg)
    p_sn = ap._p_anti(g, sigma_eps, e_min, anti_deg)
    return p_cp, p_sn, vc.row_p_anti(cls, passed, p_cp, p_sn, f_hat)


def decide(n_passed: int, n_rows: int, p_primary: float) -> str:
    """Frozen outcome (PLAN.md): gate checks first, then #273's outcome()."""
    if n_rows == 0 or n_passed < MIN_GATE_PASSED or n_passed / n_rows < MIN_GATE_PASS_FRAC:
        return "undetermined (deblend gate fails on too many rows)"
    return vc.outcome(p_primary)


def tail(probs, k: int) -> float:
    return float(vc.vn.excess_p(list(probs), k))


def posthoc(t: Table) -> dict:
    """Post hoc (not in the plan): anti observed vs expected within each gate-passed class."""
    out = {}
    passed = np.asarray(t["gate_passed"], bool)
    anti = np.asarray(t["screen_class"]) == "anti"
    for name, sel in (
        ("two_knot_gate_passed", passed & (np.asarray(t["deblend_class"]) == "two-knot")),
        ("other_gate_passed", passed & (np.asarray(t["deblend_class"]) != "two-knot")),
        ("gate_failed", ~passed),
    ):
        probs = [float(v) for v in t["p_anti_row"][sel]]
        k = int(np.sum(anti & sel))
        out[name] = {
            "n": int(np.sum(sel)),
            "anti": k,
            "expected": round(sum(probs), 3),
            "p_ge_k": tail(probs, k) if probs else None,
            "p_le_k": float(1.0 - tail(probs, k + 1)) if probs else None,
        }
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--fields", nargs="*", default=list(FIELDS_B))
    p.add_argument("--fit-arcsec", type=float, default=0.4, help="ASSUMPTION")
    p.add_argument("--n-inject", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261108)
    p.add_argument("--sigma-eps", type=float, default=0.26, help="ASSUMPTION (screen value)")
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--posthoc", action="store_true", help="only re-read census.ecsv, post hoc")
    a = p.parse_args(argv)
    if a.posthoc:
        ph = {"label": "post hoc (not in PLAN.md); derived from census.ecsv"}
        ph |= posthoc(Table.read(a.out / "census.ecsv"))
        rows = load_rows()
        ph["anti_rows_inputs"] = {
            f"{r['field']} {r['label']}": {
                "z_phot": round(float(r["z_phot"]), 2),
                "mag": round(float(r["mag"]), 1),
                "g_min": round(float(r["g_min"]), 2),
            }
            for r in rows[rows["orientation_class"] == "anti"]
        }
        (a.out / "posthoc.json").write_text(json.dumps(ph, indent=1) + "\n")
        print(json.dumps(ph, indent=1))
        return 0
    rows = load_rows(tuple(a.fields))
    full = tuple(a.fields) == tuple(FIELDS_B)
    lam_check = sum(row_probs(g, "single", True, 0.0, a.sigma_eps)[1] for g in rows["g_min"])
    if full and round(lam_check, 3) != POOLED_B_EXPECTED:
        raise SystemExit(
            f"error: shape-noise sum {lam_check:.3f} != pooled tally {POOLED_B_EXPECTED}"
        )

    res, px_scale = {}, {}
    with tempfile.TemporaryDirectory() as tmp, ProcessPoolExecutor(a.workers) as ex:
        futs = {}
        gi = 0
        for f in a.fields:
            sel = rows[rows["field"] == f]
            tg = Table({"source_uid": list(sel["label"]), "ra": sel["ra"], "dec": sel["dec"]})
            c = cutouts.make_cutouts(
                image_url(f),
                tg,
                size_arcsec=8.0,
                out_dir=Path(tmp) / f.replace(" ", "_"),
                extensions=["ERR", "WHT"],
            )
            path = {str(u): str(q) for u, q in zip(c["source_uid"], c["path"], strict=True)}
            outside = sorted(lab for lab, q in path.items() if not q)
            if outside:
                raise SystemExit(f"error: {f} rows outside the F200W image: {outside}")
            px_scale[f] = float(vc.vd.load(path[str(sel["label"][0])])[4])
            for r in sel:
                key = (f, str(r["label"]))
                futs[key] = ex.submit(
                    vc.deblend_row,
                    path[key[1]],
                    float(r["ra"]),
                    float(r["dec"]),
                    a.fit_arcsec,
                    a.n_inject,
                    a.seed + gi,
                )
                gi += 1
        res = {k: fu.result() for k, fu in futs.items()}

    keys = [(str(r["field"]), str(r["label"])) for r in rows]
    passed = [k for k in keys if res[k]["gate"]["passed"]]
    n_two = sum(res[k]["class"] == "two-knot" for k in passed)
    f_hat = n_two / len(passed) if passed else float("nan")
    f_up = vc.clopper_pearson_upper(n_two, len(passed)) if passed else float("nan")
    gmin = dict(zip(keys, (float(g) for g in rows["g_min"]), strict=True))
    pr = {
        k: row_probs(gmin[k], res[k]["class"], res[k]["gate"]["passed"], f_hat, a.sigma_eps)
        for k in keys
    }
    anti = [k for k, c in zip(keys, rows["orientation_class"], strict=True) if c == "anti"]
    k_b = len(anti)
    p_primary = tail((pr[k][2] for k in keys), k_b) if len(passed) >= MIN_GATE_PASSED else np.nan

    per_field = {}
    for f in a.fields:
        fk = [k for k in keys if k[0] == f]
        fp = [k for k in fk if res[k]["gate"]["passed"]]
        nt = sum(res[k]["class"] == "two-knot" for k in fp)
        per_field[f] = {
            "image": image_url(f),
            "pixel_scale_arcsec": px_scale[f],
            "n_rows": len(fk),
            "n_gate_passed": len(fp),
            "n_two_knot_gate_passed": nt,
            "f_two_knot": round(nt / len(fp), 3) if fp else None,
            "anti": sum(k in anti for k in fk),
            "expected_shape_noise": round(sum(pr[k][1] for k in fk), 3),
            "expected_clump_aware": round(sum(pr[k][2] for k in fk), 3),
        }
    sens = {}
    for name, f in (("f_hat", f_hat), ("f_95_upper", f_up)):
        mix = [f * pr[k][0] + (1 - f) * pr[k][1] for k in keys]
        sens[name] = {"f": round(f, 3), "p_ge_k": tail(mix, k_b)}
    persist_keys = [k for k in anti if k in PERSIST]
    secondary = {
        "p_ge_persisting": tail((pr[k][2] for k in keys), len(persist_keys))
        if full and len(passed) >= MIN_GATE_PASSED
        else None,
        "k_persisting": len(persist_keys),
    }
    m1423 = ap.ROOT / "macs1423-canucs" / "vetting" / "census.ecsv"
    if full and m1423.exists():
        t = Table.read(m1423)
        probs_a = [pr[k][2] for k in keys] + [float(v) for v in t["p_anti_row"]]
        k_a = k_b + int(np.sum(t["screen_class"] == "anti"))
        secondary["variant_A_with_macs1423_stored"] = {"k": k_a, "p_ge_k": tail(probs_a, k_a)}
    sample = {
        "n_rows": len(keys),
        "n_gate_passed": len(passed),
        "n_two_knot_gate_passed": n_two,
        "f_two_knot_pooled": round(f_hat, 3),
        "f_two_knot_95_upper": round(f_up, 3),
        "k_anti": k_b,
        "expected_shape_noise": round(sum(pr[k][1] for k in keys), 3),
        "expected_clump_aware": round(sum(pr[k][2] for k in keys), 3),
        "p_ge_k_shape_noise_check": tail((pr[k][1] for k in keys), k_b),
        "p_ge_k_clump_aware": p_primary,
        "classes_of_anti": {
            f"{k[0]} {k[1]}": [res[k]["class"], bool(res[k]["gate"]["passed"])] for k in anti
        },
        "two_knot_among_non_anti": sum(
            res[k]["class"] == "two-knot" and res[k]["gate"]["passed"]
            for k in keys
            if k not in anti
        ),
        "outcome": decide(len(passed), len(keys), p_primary),
    }
    t = Table(
        {
            "field": [k[0] for k in keys],
            "label": [k[1] for k in keys],
            "g_min": [gmin[k] for k in keys],
            "screen_class": [str(v) for v in rows["orientation_class"]],
            "mag": rows["mag"],
            "dbic": [res[k]["dbic"] for k in keys],
            "snr_min": [min(res[k]["snr1"], res[k]["snr2"]) for k in keys],
            "sep_arcsec": [res[k]["sep_arcsec"] for k in keys],
            "pa_sep": [res[k]["pa_sep"] for k in keys],
            "deblend_class": [res[k]["class"] for k in keys],
            "gate_single_ok": [res[k]["gate"]["frac_single_not_two_knot"] for k in keys],
            "gate_two_ok": [res[k]["gate"]["frac_two_called_two_knot"] for k in keys],
            "gate_passed": [res[k]["gate"]["passed"] for k in keys],
            "p_cp": [pr[k][0] for k in keys],
            "p_sn": [pr[k][1] for k in keys],
            "p_anti_row": [pr[k][2] for k in keys],
        },
        meta={"provenance": "derived (exploratory; results/arcs/pooled-census/PLAN.md)"},
    )
    a.out.mkdir(parents=True, exist_ok=True)
    t.write(a.out / "census.ecsv", overwrite=True)
    out = {
        "plan": "PLAN.md",
        "label": "exploratory; pixels observed, fits derived, injections simulated, "
        "P model_prediction",
        "fields": list(a.fields),
        "assumptions": {k: getattr(a, k) for k in ("fit_arcsec", "n_inject", "seed", "sigma_eps")}
        | {"min_gate_pass_frac": MIN_GATE_PASS_FRAC, "min_gate_passed": MIN_GATE_PASSED},
        "sample": sample,
        "per_field": per_field,
        "sensitivity_mixture": sens,
        "secondary": secondary,
    }
    (a.out / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(sample, indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
