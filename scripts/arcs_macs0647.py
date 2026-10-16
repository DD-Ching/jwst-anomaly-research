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
from epoch_compare import mutual_matches  # noqa: E402

PRIMARY = "macs0647-clash-nfw"
SECOND = "macs0647-clash-ltm"
CONTROL_MIN_ALIGNED = 0.70  # ASSUMPTION (as for MACS0417): background aligned fraction
P_REPLICATE = 0.05  # ASSUMPTION: shape-noise Poisson P(>= n_anti) below this replicates
MATCH_ARCSEC = 0.3  # ASSUMPTION: epoch-2 match radius
VALIDATE_MAX_REL_DIFF = 1e-3  # ASSUMPTION: kappa from the deflections must reproduce the map
OFFSET_ARCSEC = 0.3  # ASSUMPTION (post hoc): frame-offset sensitivity step, cf. D-034's 0.1-0.3"
REPO = Path(__file__).resolve().parents[1]
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
    """Class in ``other`` of each row of ``rows``: by label (``radius_arcsec`` None, same
    catalogue) or by mutual nearest neighbour within ``radius_arcsec``; 'not selected' where
    there is no match."""
    if radius_arcsec is None:
        lookup = {int(r["label"]): str(r["orientation_class"]) for r in other}
        return [lookup.get(int(v), "not selected") for v in rows["label"]]
    out = ["not selected"] * len(rows)
    if not len(other) or not len(rows):
        return out
    c = SkyCoord(rows["ra"], rows["dec"], unit="deg")
    co = SkyCoord(other["ra"], other["dec"], unit="deg")
    for i, j in zip(*mutual_matches(c, co, radius_arcsec), strict=True):
        out[i] = str(other["orientation_class"][j])
    return out


def validate_passes(summary: dict) -> bool:
    """Whether ``validate``'s kappa check reproduces the published map (median rel. diff.)."""
    return float(summary["kappa_map"]["median_rel_diff"]) <= VALIDATE_MAX_REL_DIFF


def _relative(value):
    """``value`` with every absolute path under the repository made repo-relative."""
    if isinstance(value, dict):
        return {k: _relative(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_relative(v) for v in value]
    if isinstance(value, str) and value and Path(value).is_absolute():
        try:
            return Path(value).resolve().relative_to(REPO).as_posix()
        except (ValueError, OSError):
            return value
    return value


def _relative_paths(path: Path) -> None:
    """Store repository-relative paths in a written JSON file (no session paths in git)."""
    path.write_text(json.dumps(_relative(json.loads(path.read_text())), indent=1) + "\n")


def _arcs(model: str, out: Path, catalog: Path, photoz: Path) -> dict:
    lc.main(
        ["--model", model, "--out", str(out), "arcs", "--catalog", str(catalog)]
        + ["--photoz", str(photoz)]
    )
    _relative_paths(out / model / "arcs.json")
    return json.loads((out / model / "arcs.json").read_text())


def _run(model: str, out: Path, catalog: Path, photoz: Path) -> dict:
    lc.main(["--model", model, "--out", str(out), "validate"])
    _relative_paths(out / model / "validate.json")
    if not validate_passes(json.loads((out / model / "validate.json").read_text())):
        raise SystemExit(f"{model}: validate fails; no arcs run")
    return _arcs(model, out, catalog, photoz)


def _background_by_gmin(t: Table) -> dict:
    b = t[t["strong_shear"] & (t["population"] == "background")]
    d = {}
    for name, sel in (("g_min<1", b["g_min"] < 1), ("g_min>=1", b["g_min"] >= 1)):
        x = b[sel]
        d[name] = {
            "n": len(x),
            "aligned": int((x["orientation_class"] == "aligned").sum()),
            "anti": int((x["orientation_class"] == "anti").sum()),
        }
    return d


def posthoc(args) -> dict:
    """Post-hoc diagnostics (not in the plan; they never change the frozen outcome): background
    classes by g_min, the epoch-1 anti rows' nearest epoch-2 rows, and the background alignment
    with the model frame shifted by +-0.3" along RA and Dec."""
    out: dict = {
        "note": "post hoc diagnostic (not in the plan); does not change the frozen outcome"
    }
    for ep in ("epoch1", "epoch2"):
        for m in (PRIMARY, SECOND):
            out[f"{ep}/{m}"] = _background_by_gmin(Table.read(args.out / ep / m / "arcs.ecsv"))
    a = Table.read(args.out / "epoch1" / PRIMARY / "arcs.ecsv")
    b = Table.read(args.out / "epoch2" / PRIMARY / "arcs.ecsv")
    anti = a[
        a["strong_shear"] & (a["population"] == "background") & (a["orientation_class"] == "anti")
    ]
    cb = SkyCoord(b["ra"], b["dec"], unit="deg")
    nearest = {}
    for r in anti:
        sep = SkyCoord(r["ra"], r["dec"], unit="deg").separation(cb).arcsec
        j = int(np.argmin(sep))
        nearest[str(int(r["label"]))] = {
            "g_min": round(float(r["g_min"]), 3),
            "mag": round(float(r["mag"]), 2),
            "epoch2_label": int(b["label"][j]),
            "sep_arcsec": round(float(sep[j]), 3),
            "epoch2_class": str(b["orientation_class"][j]),
        }
    out["epoch1_primary_anti_nearest_epoch2_row"] = nearest
    shifts = {}
    for m in (PRIMARY, SECOND):
        for dra, ddec in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            off = (dra * OFFSET_ARCSEC, ddec * OFFSET_ARCSEC)
            lc.MODELS[m]["frame_offset_arcsec"] = off
            try:
                tmp = args.out / "offset_scratch"
                r = _arcs(m, tmp, args.epoch1, args.photoz)["strong_shear_background"]
            finally:
                lc.MODELS[m].pop("frame_offset_arcsec", None)
            shifts[f"{m} {off[0]:+.1f},{off[1]:+.1f}"] = {
                k: r[k] for k in ("n", "aligned", "anti")
            } | {"aligned_fraction": round(r["aligned"] / r["n"], 3) if r["n"] else None}
    out["frame_offset_sensitivity_epoch1"] = shifts
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--epoch1", type=Path, required=True)
    ap.add_argument("--epoch2", type=Path, required=True)
    ap.add_argument("--photoz", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--posthoc", action="store_true", help="only the post-hoc diagnostics")
    args = ap.parse_args(argv)
    if args.posthoc:
        (args.out / "posthoc.json").write_text(json.dumps(posthoc(args), indent=1) + "\n")
        return 0

    runs = {m: _run(m, args.out / "epoch1", args.epoch1, args.photoz) for m in (PRIMARY, SECOND)}
    for m in (PRIMARY, SECOND):
        _arcs(m, args.out / "epoch2", args.epoch2, args.photoz)
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
        ep2 = Table.read(args.out / "epoch2" / dec / "arcs.ecsv")
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
                "class_epoch2": c_ep2,
            }
            for r, c_other, c_ep2 in zip(
                anti,
                match_class(anti, other, None),
                match_class(anti, ep2, MATCH_ARCSEC),
                strict=True,
            )
        ]
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("outcome", "deciding_model", "p_ge_observed")}))
    return 0


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    raise SystemExit(main())
