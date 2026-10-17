"""W1 shear screen: the aperture-coverage cut on grid centres (``results/w1_edge_cut/PLAN.md``).

``diagnose`` re-runs each field's real screen (``inject_shear`` inputs, cut off) and reports the
aperture coverage of its S_max and |S_×| max centres and the share of centres a cut would drop.
``decide`` applies the plan's frozen rules to two ``inject_shear.py`` output directories (cut off
and on): the c = 0 reproduction check, the D-053 rule per field and the adoption rule. Every
output is ``derived``.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

_DIR = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("inject_shear", _DIR / "inject_shear.py")
ish = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ish)
es, ir = ish.es, ish.ir

LIMIT_MASSES = ("2e+12", "8e+12", "2e+13")  # the mass cells with a limit at c = 0 (PLAN rule 2)
COMMITTED = {  # the committed c = 0 runs the reproduction check compares with
    "abell2744": "results/w1_shear_s1063/rerun_4fields/abell2744_summary.json",
    "abell370": "results/w1_shear_s1063/rerun_4fields/abell370_summary.json",
    "macs0416": "results/w1_shear_s1063/rerun_4fields/macs0416_summary.json",
    "macs1149": "results/w1_shear_s1063/rerun_4fields/macs1149_summary.json",
    "abells1063": "results/w1_shear_s1063/inject_summary_abells1063.json",
}


def recoveries(summary: dict) -> dict[str, int]:
    return {k: v["n_recovered"] for k, v in summary["efficiency"].items()}


def adopt(limits_off: dict, limits_on: dict, masses=LIMIT_MASSES) -> bool:
    """PLAN rule 2: the cut is adopted only if the combined 95 % upper limit is lower with it at
    *each* of ``masses``; a missing (None) limit with the cut counts as not lower."""
    for m in masses:
        off, on = limits_off[m]["upper_limit_deg2"], limits_on[m]["upper_limit_deg2"]
        if on is None or off is None or not on < off:
            return False
    return True


def eb_pass(screen: dict) -> bool:
    """PLAN rule 1 (D-053) on one field's real screen: p_rot < 0.05 and S_max > b_max."""
    return bool(screen["p_random_max"] < ir.P_RECOVER and screen["s_max"] > screen["b_max"])


def cmd_decide(args) -> dict:
    out = {"fields": {}, "reproduces_committed": True}
    for f in ish.DEFAULT_FIELDS:
        off = json.loads((args.off / f / "summary.json").read_text())
        on = json.loads((args.on / f / "summary.json").read_text())
        ref = json.loads(Path(COMMITTED[f]).read_text())
        same = recoveries(off) == recoveries(ref)
        out["reproduces_committed"] &= same
        out["fields"][f] = {
            "recovered_c0": recoveries(off),
            "recovered_committed": recoveries(ref),
            "c0_equals_committed": same,
            "recovered_cut": recoveries(on),
            "b_max_c0": off["screen"]["b_max"],
            "b_max_cut": on["screen"]["b_max"],
            "s_max_cut": on["screen"]["s_max"],
            "p_random_max_cut": on["screen"]["p_random_max"],
            "n_centres_valid_c0": off["screen"]["n_centres_valid"],
            "n_centres_valid_cut": on["screen"]["n_centres_valid"],
            "passes_eb_rule_cut": eb_pass(on["screen"]),
        }
    lim_off = json.loads((args.off / "combined.json").read_text())["limits_95"]
    lim_on = json.loads((args.on / "combined.json").read_text())["limits_95"]
    out["limits_c0"] = {m: lim_off[m]["upper_limit_deg2"] for m in lim_off}
    out["limits_cut"] = {m: lim_on[m]["upper_limit_deg2"] for m in lim_on}
    out["any_real_pass_cut"] = any(v["passes_eb_rule_cut"] for v in out["fields"].values())
    out["adopt_cut"] = bool(out["reproduces_committed"] and adopt(lim_off, lim_on))
    out["provenance"] = "derived"
    args.dest.write_text(json.dumps(es.json_safe(out), indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return out


def cmd_diagnose(args) -> dict:
    out = {"min_coverage_tested": args.min_coverage, "provenance": "derived", "fields": {}}
    for f in ish.DEFAULT_FIELDS:
        spec = ish.field_spec(f)
        model, shapes, extra, _, _ = ir.load_field(f, args.work / f, spec)
        sargs = SimpleNamespace(**ish.SHEAR_DEFAULTS)
        sargs.max_radius = spec["max_radius"]
        res = es.shear_screen(model, shapes, sargs, es.psf_sigma_px(shapes), extra)
        xa, ya = model.to_frame(np.asarray(shapes["ra"]), np.asarray(shapes["dec"]))
        cov = es.aperture_coverage(
            xa, ya, res["gx"], res["gy"], sargs.aperture_arcsec, sargs.r_min_arcsec
        )
        s, sx = res["s"], res["s_cross"]
        fin = np.isfinite(s)
        row = {"n_centres_valid": int(fin.sum())}
        for name, m in (("s_max", s), ("abs_s_cross_max", np.abs(sx))):
            iy, ix = np.unravel_index(np.nanargmax(m), m.shape)
            row[name] = {
                "value": float(m[iy, ix]),
                "xy": [float(res["gx"][iy, ix]), float(res["gy"][iy, ix])],
                "coverage": float(cov[iy, ix]),
            }
        row["frac_valid_below_cut"] = float(np.mean(cov[fin] < args.min_coverage))
        out["fields"][f] = row
    args.dest.write_text(json.dumps(es.json_safe(out), indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("diagnose")
    d.add_argument("--work", type=Path, required=True, help="scratch dir (Gaia star tables)")
    d.add_argument("--min-coverage", type=float, default=0.5)
    d.add_argument("--dest", type=Path, required=True)
    c = sub.add_parser("decide")
    c.add_argument("--off", type=Path, required=True, help="inject_shear.py --out, cut off")
    c.add_argument("--on", type=Path, required=True, help="inject_shear.py --out, cut on")
    c.add_argument("--dest", type=Path, required=True)
    args = ap.parse_args(argv)
    {"diagnose": cmd_diagnose, "decide": cmd_decide}[args.cmd](args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
