"""`arcs` on MACS J1423.8+2404: a further independent test of the Abell S1063 anti excess.

Plan (frozen before the run): ``results/arcs/macs1423-canucs/PLAN.md``. This script converts the
CANUCS DR1 photometry to zout columns, runs ``lens_consistency.py validate`` and ``arcs`` (D-024
defaults) with the CANUCS DR1 maps on the program 1208 F200W catalogue, and applies the frozen
outcome rule (:func:`outcome`). Shapes are ``observed``, predicted axes ``model_prediction``,
classes, counts and P values ``derived``; thresholds are ASSUMPTIONs.

    python scripts/arcs_macs1423.py --catalog <jw01208-o010_t005 F200W _cat.ecsv> \
        --photometry <macs1423-clu photometry-cat.fits.gz>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import paths

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arcs_macs0647 import P_REPLICATE, _background_by_gmin, _run, control_passes  # noqa: E402
from inject_radial import canucs_zout  # noqa: E402

MODEL = "macs1423-canucs-maps"
OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "macs1423-canucs"


def outcome(run: dict) -> dict:
    """The frozen rule of PLAN.md for one ``arcs.json`` summary: no result if the control
    fails; else replicates if P(>= n_anti) < 0.05 under the shape-noise expectation."""
    if not control_passes(run["strong_shear_background"]):
        return {"outcome": "no result", "p_ge_observed": None}
    p = float(run["anti_expected_shape_noise_background"]["p_ge_observed"])
    return {
        "outcome": "replicates" if p < P_REPLICATE else "does not replicate",
        "p_ge_observed": p,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--photometry", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    zout = paths.cache_dir() / "macs1423" / "macs1423-clu_zout.fits"  # gitignored
    zout.parent.mkdir(parents=True, exist_ok=True)
    canucs_zout(args.photometry).write(zout, overwrite=True)
    run = _run(MODEL, args.out, args.catalog, zout)
    t = Table.read(args.out / MODEL / "arcs.ecsv")
    anti = t[
        t["strong_shear"] & (t["population"] == "background") & (t["orientation_class"] == "anti")
    ]
    summary = {
        "plan": "results/arcs/macs1423-canucs/PLAN.md",
        "inputs": {"catalog": args.catalog.name, "photometry": args.photometry.name},
        "model": MODEL,
        "strong_shear_background": run["strong_shear_background"],
        "anti_expected_shape_noise_background": run["anti_expected_shape_noise_background"],
        "control_passes": control_passes(run["strong_shear_background"]),
        "background_by_gmin": _background_by_gmin(t),
        **outcome(run),
        "anti_background": [
            {
                "label": int(r["label"]),
                "ra": round(float(r["ra"]), 6),
                "dec": round(float(r["dec"]), 6),
                "mag": round(float(r["mag"]), 2),
                "g_min": round(float(r["g_min"]), 3),
                "z_phot": round(float(r["z_phot"]), 2),
                "ellipticity": round(float(r["ellipticity"]), 3),
            }
            for r in anti
        ],
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("outcome", "p_ge_observed", "control_passes")}))
    return 0


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    raise SystemExit(main())
