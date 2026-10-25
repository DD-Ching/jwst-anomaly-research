"""SLICE JWST-era image lists vs HST-era RELICS GLAFIC v2 maps (MACS J0553, RXC J2211).

Frozen plan: ``results/lens_consistency/slice_relics_glafic/PLAN.md`` (D-034 addendum bd). Same
image list, families and rules 1-6 as ``slice_relics_heldout`` (addendum bc), with the GLAFIC
model in place of Lenstool v1; adds the frozen amplitude factor (rule 7), the product check gate P
(rule 8) and the two-model statement for MACS J0553 11.1 (rule 9). Writes ``summary.json`` and
per-model tables.

    python scripts/lens_consistency.py validate --model macs0553-relics-glafic \
        --out results/lens_consistency/slice_relics_glafic/validate   # gate P inputs, per model
    python scripts/slice_relics_glafic.py [--grid-step 0.25]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import macs0647_clash_heldout as ho
import macs2129_clash_heldout as h29
import numpy as np
import slice_relics_heldout as srh
from astropy.table import Table

OUT = Path("results/lens_consistency/slice_relics_glafic")
LENSTOOL = Path("results/lens_consistency/slice_relics_heldout/summary.json")
# cluster -> GLAFIC model; families and Table 8 blocks come from srh.CLUSTERS
MODELS = {"macs0553": "macs0553-relics-glafic", "rxc2211": "rxc2211-relics-glafic"}
NEAR_ONE = (0.85, 1.15)  # rule 7, ASSUMPTION
LENSTOOL_LIKE_MAX = 0.75  # rule 7, ASSUMPTION
P_RATIO = (0.95, 1.05)  # rule 8, ASSUMPTION
FLAG_11 = "11.1"  # rule 9


def amplitude_class(median: float) -> str:
    """Rule 7: classify the median D_LS/D_S factor over families."""
    if NEAR_ONE[0] <= median <= NEAR_ONE[1]:
        return "near_1"
    if median <= LENSTOOL_LIKE_MAX:
        return "lenstool_like"
    return "intermediate"


def gate_p(validate: dict) -> bool:
    """Rule 8: κ and μ(z = 2) median ratios (deflection-derived / published) within P_RATIO."""
    ratios = [validate["kappa_map"]["median_ratio"], validate["magnification_z2"]["median_ratio"]]
    return all(P_RATIO[0] <= r <= P_RATIO[1] for r in ratios)


def two_model_11(usable_glafic: bool, glafic_dtheta: float | None, lenstool_flagged: bool) -> str:
    """Rule 9: the 11.1 statement from the GLAFIC run and the frozen Lenstool flag."""
    if not usable_glafic:
        return "glafic_not_usable"
    glafic_miss = glafic_dtheta is not None and glafic_dtheta > ho.MATCH_ARCSEC
    if glafic_miss:
        return "missed_by_both" if lenstool_flagged else "glafic_only"
    return "lenstool_only" if lenstool_flagged else "neither"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid-step", type=float, default=0.25)
    args = ap.parse_args(argv)
    images = Table.read(srh.IMAGES, format="ascii.ecsv")
    lenstool = json.loads(LENSTOOL.read_text())
    summary: dict = {
        "plan": str(OUT / "PLAN.md"),
        "rules": {
            "K0_MAX_RMS_ARCSEC": ho.K0_MAX_RMS_ARCSEC,
            "K0_MIN_MATCHED": ho.K0_MIN_MATCHED,
            "MATCH_ARCSEC": ho.MATCH_ARCSEC,
            "FRAME_MAX_ARCSEC": ho.FRAME_MAX_ARCSEC,
            "z_grid": [ho.Z_STEPS, ho.Z_MIN_OFFSET, ho.Z_MAX],
            "grid_step": args.grid_step,
            "near_one": NEAR_ONE,
            "lenstool_like_max": LENSTOOL_LIKE_MAX,
            "p_ratio": P_RATIO,
        },
    }
    for cluster, model in MODELS.items():
        t = srh.load_cluster(images, cluster)
        held = np.asarray(t["held_out"], bool)
        spec = np.asarray(t["z_spec"], bool)
        validate = json.loads((OUT / "validate" / model / "validate.json").read_text())
        p_ok = gate_p(validate)
        entry, ip, ds = h29.run_model(model, t, args.grid_step)
        usable = entry["usable"] and p_ok
        names = [model] if usable else []
        ra = ho.readings({model: np.asarray(ip["dtheta_arcsec"], float)[held]}, names)
        rb = h29.reading_b({model: ds}, spec, names) if spec.any() else np.array([], str)
        # rule 7 uses the addendum-bc factor definition, on this model's fitted z
        amp = srh.posthoc_scale({cluster: entry}, images, {cluster: model})[cluster]
        factors = [r["scale"] for r in amp["families"].values()]
        med = float(np.median(factors)) if factors else None
        entry.update(
            {
                "model": model,
                "gate_P": {"passes": p_ok, "validate": str(OUT / "validate" / model)},
                "usable_with_P": usable,
                "images": len(t),
                "held_out_images": int(held.sum()),
                "spec_z_images": int(spec.sum()),
                "flags_A": [
                    str(i) for i, v in zip(t["image_id"][held], ra, strict=True) if v != "no_flag"
                ],
                "flags_B": [
                    str(i) for i, v in zip(t["image_id"][spec], rb, strict=True) if v != "no_flag"
                ],
                "amplitude": {
                    **amp,
                    "median": med,
                    "min": min(factors, default=None),
                    "max": max(factors, default=None),
                    "class": amplitude_class(med) if usable and med is not None else None,
                },
            }
        )
        if cluster == "macs0553":
            ids = np.asarray(ip["image_id"], str)
            d11 = float(np.asarray(ip["dtheta_arcsec"], float)[ids == FLAG_11][0])
            entry["image_11_1"] = {
                "glafic_dtheta_arcsec": d11,
                "lenstool_flagged": FLAG_11 in lenstool[cluster]["flags_A"],
                "statement": two_model_11(usable, d11, FLAG_11 in lenstool[cluster]["flags_A"]),
            }
        if not spec.any():
            entry["spec_z_B"] = entry["spec_fitted_z"] = None
        summary[cluster] = entry
        ip.meta["model"] = model
        lc._write(ip, OUT / f"imageplane_{cluster}.ecsv")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    for c in MODELS:
        e = summary[c]
        print(c, e["usable_with_P"], e["gate_K0"], e["frame_F"], e["gate_P"], sep="\n  ")
        print("  A", e["held_out_A"], "\n  B", e["spec_z_B"])
        print("  flags", e["flags_A"], e["flags_B"])
        print("  amplitude", {k: e["amplitude"][k] for k in ("median", "min", "max", "class")})
        if "image_11_1" in e:
            print("  11.1", e["image_11_1"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
