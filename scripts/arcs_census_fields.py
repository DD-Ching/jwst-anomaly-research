"""Clump-pair census on five other ``arcs`` fields, to test the S1063 single-field clump rate.

Frozen plan: ``results/arcs/census-fields/PLAN.md`` (written before the run). The per-row deblend
and own-flux gate are the unchanged S1063 census code (``arcs_vet_census.deblend_row``). The pooled
two-knot fraction then replaces the S1063-only f̂ in the S1063 mixture null.

    python scripts/arcs_census_fields.py [--workers 4] [--cut-dir DIR] [--sheet PNG]

Output: ``results/arcs/census-fields/census_fields.{json,ecsv}`` (derived; pixels observed,
injections simulations, f̂ and P model_prediction).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")  # pool workers: no BLAS oversubscription

import argparse
import importlib.util
import json
import re
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

from jwst_anomaly import cutouts, paths

FIELDS = (
    "abell2744-bergamini23",
    "abell370-canucs",
    "elgordo-caminha23",
    "macs0416-canucs",
    "macs1149-cats",
)
S1063 = "abells1063-cats"
K_S1063 = 4  # the S1063 screen anti count (census "k4_screen")


def _load_module(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / file)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # pool workers pickle its functions by module name
    spec.loader.exec_module(mod)
    return mod


census = _load_module("arcs_vet_census", "arcs_vet_census.py")


def i2d_uri(catalog: str) -> str:
    """MAST S3 URI of the F200W i2d mosaic a pipeline ``_cat.ecsv`` was made from."""
    m = re.search(r"jw(\d{5})-o(\d{3})_(t\d{3})_nircam_clear-f200w_cat\.ecsv$", catalog)
    if not m:
        raise ValueError(f"not a NIRCam F200W L3 catalogue name: {catalog}")
    prog, obs, tgt = m.groups()
    return (
        f"s3://stpubdata/jwst/public/jw{prog}/L3/t/o{obs}/"
        f"jw{prog}-o{obs}_{tgt}_nircam_clear-f200w_i2d.fits"
    )


def strong_background(t: Table) -> Table:
    """Rows the ``arcs`` sample test uses: strong shear, background population."""
    return t[np.asarray(t["strong_shear"], bool) & (np.asarray(t["population"]) == "background")]


def mixture_p(p_cp, p_sn, f: float, k: int) -> float:
    """P(≥ k anti) when each row is a clump pair with probability f (S1063 census mixture)."""
    p = f * np.asarray(p_cp, float) + (1.0 - f) * np.asarray(p_sn, float)
    return census.vn.excess_p(p, k)


def plan_outcome(p_pooled: float) -> str:
    if not np.isfinite(p_pooled):
        return "undetermined (no gate-passed rows)"
    if p_pooled >= 0.05:
        return "S1063 excess explained by clump pairs at the five-field rate (exploratory)"
    return "S1063 excess remains at the five-field clump rate (exploratory)"


def field_rows(arcs_dir: Path, model: str) -> tuple[Table, str]:
    meta = json.loads((arcs_dir / model / "arcs.json").read_text())
    rows = strong_background(Table.read(arcs_dir / model / "arcs.ecsv"))
    return rows, i2d_uri(Path(meta["catalog"]).name)


def run_field(model, rows, uri, base: Path, ex, a) -> Table:
    tg = Table(
        {"source_uid": [str(v) for v in rows["label"]], "ra": rows["ra"], "dec": rows["dec"]}
    )
    c = cutouts.make_cutouts(
        uri, tg, size_arcsec=8.0, out_dir=base / model, extensions=["ERR", "WHT"]
    )
    path = {str(u): str(p) for u, p in zip(c["source_uid"], c["path"], strict=True)}
    keep = [k for k, lab in enumerate(tg["source_uid"]) if path.get(lab)]
    futs = [
        ex.submit(
            census.deblend_row,
            path[tg["source_uid"][k]],
            float(rows["ra"][k]),
            float(rows["dec"][k]),
            a.fit_arcsec,
            a.n_inject,
            a.seed + k,
        )
        for k in keep
    ]
    res = [f.result() for f in futs]
    return Table(
        {
            "field": [model] * len(keep),
            "label": [str(rows["label"][k]) for k in keep],
            "ra": [float(rows["ra"][k]) for k in keep],
            "dec": [float(rows["dec"][k]) for k in keep],
            "mag": [float(rows["mag"][k]) for k in keep],
            "g_min": [float(rows["g_min"][k]) for k in keep],
            "screen_class": [str(rows["orientation_class"][k]) for k in keep],
            "dbic": [r["dbic"] for r in res],
            "snr_min": [min(r["snr1"], r["snr2"]) for r in res],
            "sep_arcsec": [r["sep_arcsec"] for r in res],
            "deblend_class": [r["class"] for r in res],
            "gate_single_ok": [r["gate"]["frac_single_not_two_knot"] for r in res],
            "gate_two_ok": [r["gate"]["frac_two_called_two_knot"] for r in res],
            "gate_passed": [bool(r["gate"]["passed"]) for r in res],
            "cutout": [path[tg["source_uid"][k]] for k in keep],
        },
        meta={"n_rows": len(rows), "n_outside": len(rows) - len(keep), "image": uri},
    )


def summarize(t: Table) -> dict:
    passed = np.asarray(t["gate_passed"], bool)
    n_pass = int(passed.sum())
    n_two = int((passed & (np.asarray(t["deblend_class"]) == "two-knot")).sum())
    return {
        "n_rows": len(t),
        "n_gate_passed": n_pass,
        "n_two_knot_gate_passed": n_two,
        "f_two_knot": n_two / n_pass if n_pass else float("nan"),
        "f_two_knot_95_upper": census.clopper_pearson_upper(n_two, n_pass)
        if n_pass
        else float("nan"),
        "interpretable_alone": n_pass >= 12,
    }


def contact_sheet(t: Table, png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sel = t[np.asarray(t["deblend_class"]) == "two-knot"]
    if not len(sel):
        return
    n = len(sel)
    ncol = min(6, n)
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.2 * ncol, 2.4 * nrow), squeeze=False)
    for ax in axes.flat:
        ax.axis("off")
    for ax, r in zip(axes.flat, sel, strict=False):
        sci, sub, err, wcs, px = census.vd.load(r["cutout"])
        x, y = wcs.world_to_pixel_values(r["ra"], r["dec"])
        h = int(round(0.6 / px))
        xi, yi = int(round(float(x))), int(round(float(y)))
        st = sub[max(yi - h, 0) : yi + h + 1, max(xi - h, 0) : xi + h + 1]
        v = np.nanpercentile(st, [5, 99.5])
        ax.imshow(st, origin="lower", cmap="gray", vmin=v[0], vmax=v[1])
        ax.set_title(
            f"{r['field'].split('-')[0]} {r['label']}\n{r['screen_class']} "
            f'sep {r["sep_arcsec"]:.2f}"{"" if r["gate_passed"] else " gate-fail"}',
            fontsize=7,
        )
    fig.suptitle('two-knot rows, F200W 1.2" stamps, 1"-filtered (north not up)', fontsize=8)
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    p.add_argument("--fields", nargs="+", default=list(FIELDS))
    p.add_argument("--fit-arcsec", type=float, default=0.4, help="ASSUMPTION (S1063 census)")
    p.add_argument("--n-inject", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261016)
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--cut-dir", type=Path, default=None, help="reuse cutouts (else a temp dir)")
    p.add_argument("--sheet", type=Path, default=None, help="contact sheet of two-knot rows")
    a = p.parse_args(argv)
    out_dir = a.arcs_dir / "census-fields"

    with tempfile.TemporaryDirectory() as tmp, ProcessPoolExecutor(a.workers) as ex:
        base = a.cut_dir or Path(tmp)
        tabs = []
        for model in a.fields:
            rows, uri = field_rows(a.arcs_dir, model)
            tabs.append(run_field(model, rows, uri, base, ex, a))
            print(model, json.dumps(summarize(tabs[-1]), default=float), flush=True)
        t = vstack(tabs, metadata_conflicts="silent")
        if a.sheet:
            contact_sheet(t, a.sheet)
    per_field = {
        m: summarize(tt) | {"n_outside": tt.meta["n_outside"], "image": tt.meta["image"]}
        for m, tt in zip(a.fields, tabs, strict=True)
    }
    pooled = summarize(t)

    s = Table.read(a.arcs_dir / S1063 / "vetting" / "census.ecsv")
    p_s1063 = {
        name: mixture_p(s["p_cp"], s["p_sn"], f, K_S1063)
        for name, f in (
            ("pooled_f_hat", pooled["f_two_knot"]),
            ("pooled_f_upper", pooled["f_two_knot_95_upper"]),
        )
    }
    out = {
        "plan": "PLAN.md",
        "label": "exploratory; fits derived, injections simulated, f and P model_prediction",
        "assumptions": {k: getattr(a, k) for k in ("fit_arcsec", "n_inject", "seed")},
        "per_field": per_field,
        "pooled_five_fields": pooled,
        "s1063_p_ge_4_anti": p_s1063,
        "outcome": plan_outcome(p_s1063["pooled_f_hat"]),
    }
    t.remove_column("cutout")
    t.meta = {"provenance": "derived (exploratory; results/arcs/census-fields/PLAN.md)"}
    t.write(out_dir / "census_fields.ecsv", overwrite=True)
    (out_dir / "census_fields.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
