"""`arcs` on MACS J0647.7+7015: a second independent test of the Abell S1063 anti excess.

Plan (frozen before the run): ``results/arcs/macs0647-clash/PLAN.md``. This script runs
``lens_consistency.py validate`` and ``arcs`` (D-024 defaults) with the two CLASH Zitrin map models
on the JWST program 1433 epoch 1 F200W catalogue, applies the frozen outcome rule
(:func:`outcome`), and reports the pre-registered secondary checks (epoch 2 class, other-model
class) for the anti rows of the deciding run. Shapes are ``observed``, predicted axes
``model_prediction``, classes, counts and P values ``derived``; thresholds are ASSUMPTIONs.

    python scripts/arcs_macs0647.py --epoch1 <o010 F200W _cat.ecsv> --epoch2 <o020 ...> \
        --photoz <macs0647-grizli-v7.0-fix.eazypy.zout.fits>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

PRIMARY = "macs0647-clash-nfw"
SECOND = "macs0647-clash-ltm"
CONTROL_MIN_ALIGNED = 0.70  # ASSUMPTION (as for MACS0417): background aligned fraction
P_REPLICATE = 0.05  # ASSUMPTION: shape-noise Poisson P(>= n_anti) below this replicates
MATCH_ARCSEC = 0.3  # ASSUMPTION: epoch-2 match radius
OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "macs0647-clash"


def control_passes(background: dict) -> bool:
    """Whether a run's strong-shear background rows are >= 70 % aligned (n = 0 fails)."""
    n = int(background["n"])
    return n > 0 and int(background["aligned"]) / n >= CONTROL_MIN_ALIGNED


def outcome(primary: dict, second: dict) -> dict:
    """The frozen rule of PLAN.md. ``primary``/``second`` are ``arcs.json`` summaries.

    The primary run decides if its control passes; otherwise the second run decides if its
    control passes; otherwise there is no result. The deciding run replicates the S1063 pattern
    if P(>= n_anti) < 0.05 under the shape-noise expectation, else it does not replicate.
    """
    for run in (primary, second):
        if control_passes(run["strong_shear_background"]):
            p = float(run["anti_expected_shape_noise_background"]["p_ge_observed"])
            return {
                "deciding_model": run["model"],
                "outcome": "replicates" if p < P_REPLICATE else "does not replicate",
                "p_ge_observed": p,
            }
    return {"deciding_model": None, "outcome": "no result", "p_ge_observed": None}


def match_class(rows: Table, other: Table, radius_arcsec: float | None) -> list[str]:
    """Class of each row in ``other``: by label (``radius_arcsec`` None, same catalogue) or by
    nearest position within ``radius_arcsec``; 'not selected' where there is no match."""
    out = []
    if radius_arcsec is None:
        lookup = {int(r["label"]): str(r["orientation_class"]) for r in other}
        return [lookup.get(int(v), "not selected") for v in rows["label"]]
    if not len(other):
        return ["not selected"] * len(rows)
    c = SkyCoord(rows["ra"], rows["dec"], unit="deg")
    co = SkyCoord(other["ra"], other["dec"], unit="deg")
    idx, sep, _ = c.match_to_catalog_sky(co)
    for i, s in zip(idx, sep.arcsec, strict=True):
        out.append(str(other["orientation_class"][i]) if s <= radius_arcsec else "not selected")
    return out


def _run(model: str, out: Path, catalog: Path, photoz: Path) -> dict:
    lc.main(["--model", model, "--out", str(out), "validate"])
    lc.main(
        ["--model", model, "--out", str(out), "arcs", "--catalog", str(catalog)]
        + ["--photoz", str(photoz)]
    )
    return json.loads((out / model / "arcs.json").read_text())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--epoch1", type=Path, required=True)
    ap.add_argument("--epoch2", type=Path, required=True)
    ap.add_argument("--photoz", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    runs = {m: _run(m, args.out / "epoch1", args.epoch1, args.photoz) for m in (PRIMARY, SECOND)}
    lc.main(
        ["--model", PRIMARY, "--out", str(args.out / "epoch2"), "arcs"]
        + ["--catalog", str(args.epoch2), "--photoz", str(args.photoz)]
    )
    res = outcome(runs[PRIMARY], runs[SECOND])
    summary = {
        "plan": "results/arcs/macs0647-clash/PLAN.md",
        "inputs": {
            "epoch1": args.epoch1.name,
            "epoch2": args.epoch2.name,
            "photoz": args.photoz.name,
        },
        "runs": {
            m: {
                "strong_shear_background": r["strong_shear_background"],
                "anti_expected_shape_noise_background": r["anti_expected_shape_noise_background"],
                "control_passes": control_passes(r["strong_shear_background"]),
            }
            for m, r in runs.items()
        },
        **res,
    }
    if res["deciding_model"]:
        dec = res["deciding_model"]
        other_model = SECOND if dec == PRIMARY else PRIMARY
        t = Table.read(args.out / "epoch1" / dec / "arcs.ecsv")
        anti = t[
            t["strong_shear"]
            & (t["population"] == "background")
            & (t["orientation_class"] == "anti")
        ]
        other = Table.read(args.out / "epoch1" / other_model / "arcs.ecsv")
        ep2 = Table.read(args.out / "epoch2" / PRIMARY / "arcs.ecsv")
        summary["anti_background"] = [
            {
                "label": int(r["label"]),
                "ra": round(float(r["ra"]), 6),
                "dec": round(float(r["dec"]), 6),
                "mag": round(float(r["mag"]), 2),
                "g_min": round(float(r["g_min"]), 3),
                "z_phot": round(float(r["z_phot"]), 2),
                "ellipticity": round(float(r["ellipticity"]), 3),
                "class_other_model": c_other,
                "class_epoch2_primary": c_ep2,
            }
            for r, c_other, c_ep2 in zip(
                anti,
                match_class(anti, other, None),
                match_class(anti, ep2, MATCH_ARCSEC),
                strict=True,
            )
        ]
        summary["epoch2_note"] = (
            "epoch-2 classes use the primary model regardless of the deciding run"
        )
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("outcome", "deciding_model", "p_ge_observed")}))
    return 0


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    raise SystemExit(main())
