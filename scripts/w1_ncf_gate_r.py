"""Gate R of results/w1_shear_ncf/PLAN.md: the shear responsivity of F150W pipeline-catalogue
moments, fitted on the CANUCS cluster pointings of program 1208 (MACS0417, MACS1423) with the
models and photo-z of their F200W runs, then compared with the stored F200W values.

    python scripts/w1_ncf_gate_r.py --out results/w1_shear_ncf/gate_r.json

Writes the fitted R per field and band and the gate outcome (``gate_r``). Fitted R ``derived``;
every threshold an ASSUMPTION (PLAN.md).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))

import inject_shear as ish  # noqa: E402

from jwst_anomaly import paths, schema  # noqa: E402

ir, es = ish.ir, ish.es
CLUSTERS = {"macs0417": "o006_t003", "macs1423": "o010_t005"}
F200W_SUMMARY = "results/w1_volume/summaries/z0.3_{0}.json"
MAX_DIFF = 0.15  # |R_F150W - R_F200W| allowed per field (ASSUMPTION, PLAN.md)


def gate_r(r150: dict, r200: dict, max_diff: float = MAX_DIFF) -> dict:
    """Gate R as frozen: every field's F150W fit usable (``shear_responsivity`` raised otherwise,
    recorded as None) *and* within ``max_diff`` of its F200W value, in every field. Adopted R is
    the smallest F150W value."""
    ok = {f: r150.get(f) is not None and abs(r150[f] - r200[f]) <= max_diff for f in sorted(r200)}
    passed = bool(ok) and all(ok.values()) and set(r150) == set(r200)
    return {
        "per_field_pass": ok,
        "passed": passed,
        "adopted_r": min(r150.values()) if passed else None,
    }


def fit_r(field: str, obs: str, out: Path) -> dict:
    spec = ish.field_spec(field) | {
        "catalog": ish.ncf_catalog(obs),
    }  # Gaia spike stars as in the F200W runs (spec["spike_stars"])
    model, shapes, extra, cat, _ = ir.load_field(field, out, spec)
    psf = es.psf_sigma_px(shapes)
    a = SimpleNamespace(**ish.SHEAR_DEFAULTS)
    try:
        _, counts = es.shear_sources(model, shapes, psf, a.min_snr, a.max_g, extra=extra)
        return {
            "catalog": str(cat),
            "psf_sigma_px": psf,
            "r": counts["responsivity"],
            "r_err": counts["responsivity_err"],
            "n_weak_shear": counts["n_weak_shear"],
            "n_spike_stars": None if extra is None else len(extra),
        }
    except ValueError as exc:  # R not usable: the gate fails for this field
        return {"catalog": str(cat), "psf_sigma_px": psf, "r": None, "error": str(exc)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    work = paths.outputs_dir() / "w1_ncf_gate_r"
    fits, r200 = {}, {}
    for field, obs in CLUSTERS.items():
        (work / field).mkdir(parents=True, exist_ok=True)
        fits[field] = fit_r(field, obs, work / field)
        stored = json.loads(Path(F200W_SUMMARY.format(field)).read_text())["screen"]
        r200[field] = stored["responsivity"]
        fits[field]["r_f200w_stored"] = r200[field]
        fits[field]["r_f200w_err_stored"] = stored["responsivity_err"]
    result = {
        "provenance": schema.Provenance.DERIVED.value,
        "plan": "results/w1_shear_ncf/PLAN.md (gate R)",
        "max_diff": MAX_DIFF,
        "fields": fits,
        "gate_r": gate_r({f: v["r"] for f, v in fits.items()}, r200),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(es.json_safe(result), indent=1) + "\n")
    print(json.dumps(result["gate_r"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
