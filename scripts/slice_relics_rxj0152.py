"""SLICE JWST-era images of RX J0152.7-1357 vs three HST-era RELICS maps (Lenstool, GLAFIC, LTM).

Frozen plan: ``results/lens_consistency/slice_relics_rxj0152/PLAN.md`` (D-034 addendum). Image
list: SLICE Table 8 block (clump 1; observed), parsed with ``slice_relics_heldout.parse_block``.
Per model: gate P, K0 and F (rules 2-4), reading A on family 7 and reading B on family 1 at z_spec
(rule 5, joint over usable models), multiplicity (rule 6) and the amplitude class (rule 7, Q1).

    python scripts/slice_relics_rxj0152.py --build-images main.tex   # once, from the e-print
    python scripts/lens_consistency.py --model cl0152-relics-lenstool \
        --out results/lens_consistency/slice_relics_rxj0152/validate validate
    python scripts/slice_relics_rxj0152.py [--grid-step 0.25]
    python scripts/slice_relics_rxj0152.py --posthoc   # post hoc: unmatched predicted images
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import lens_consistency as lc
import m0417_images as mi
import macs0647_clash_heldout as ho
import macs2129_clash_heldout as h29
import numpy as np
import slice_relics_glafic as sg
import slice_relics_heldout as srh
from astropy.table import Table

OUT = Path("results/lens_consistency/slice_relics_rxj0152")
IMAGES = OUT / "slice_images.ecsv"
CLUSTER = "rxj0152"
SIDEHEAD = "RX\\,J0152.7$-$1357"
HELD = ("7",)  # rule 5 reading A: new in the JWST imaging (SLICE text)
SPEC = ("1",)  # rule 5 reading B: z_spec 3.93 (Acebron et al. 2019)
LENSTOOL = "cl0152-relics-lenstool"
CENTRE = (28.17764473, -13.95911996)  # SLICE Table 1 (deg)
MODELS = (LENSTOOL, "cl0152-relics-glafic", "cl0152-relics-zitrin")


def q1_answer(classes: dict[str, str | None]) -> str:
    """Rule 7: H_sys / H_none / undecided from the amplitude class of each model (None = not
    usable)."""
    others = [c for m, c in classes.items() if m != LENSTOOL and c is not None]
    lt = classes.get(LENSTOOL)
    if lt is None or not others:
        return "undecided"
    if lt == "lenstool_like" and all(c == "near_1" for c in others):
        return "H_sys"
    if lt == "near_1" and all(c == "near_1" for c in others):
        return "H_none"
    return "undecided"


def per_model_misses(ids: np.ndarray, d: np.ndarray) -> list[str]:
    """Image ids whose nearest predicted image is > MATCH_ARCSEC away or missing."""
    d = np.asarray(d, float)
    miss = ~(np.isfinite(d) & (d <= ho.MATCH_ARCSEC))
    return [str(i) for i in np.asarray(ids, str)[miss]]


def unmatched_predictions(pred_xy: np.ndarray, obs_xy: np.ndarray) -> np.ndarray:
    """Indices of predicted images that are no catalogued image's nearest prediction."""
    d = np.hypot(*(obs_xy[:, None, :] - pred_xy[None, :, :]).transpose(2, 0, 1))
    return np.setdiff1d(np.arange(len(pred_xy)), np.argmin(d, axis=1))


def posthoc(grid_step: float) -> dict:
    """Post hoc (not in PLAN.md): sky positions of the predicted images that match no catalogued
    image, for every family with extra predictions, at the run's fitted z."""
    from jwst_anomaly import lensmodel

    summary = json.loads((OUT / "summary.json").read_text())
    t = load()
    fam = np.asarray(t["system"], str)
    out: dict = {"note": "post hoc, not in PLAN.md; provenance model_prediction", "models": {}}
    for model in MODELS:
        e = summary["models"][model]
        solver = mi.Solver(model, t, grid_step)
        rows = {}
        for f in e["extra_predicted"]:
            sub = t[fam == f]
            sub["z"] = e["z_fit"][f]
            bt = lensmodel.backtrace_images(solver.model, sub, {})
            pred = lensmodel.find_images(
                solver.model,
                solver.grid,
                float(np.mean(bt["beta_x"])),
                float(np.mean(bt["beta_y"])),
                float(bt["z_used"][0]),
            )
            ox, oy = solver.model.to_frame(
                np.asarray(sub["ra"], float), np.asarray(sub["dec"], float)
            )
            pxy = np.column_stack([pred["x"], pred["y"]])
            k = unmatched_predictions(pxy, np.column_stack([ox, oy]))
            ra, dec = solver.model.to_sky(pxy[k, 0], pxy[k, 1])
            mu = np.asarray(pred["magnification"], float)
            near = np.argmin(np.hypot(ox[:, None] - pxy[:, 0], oy[:, None] - pxy[:, 1]), axis=1)
            mu_cat = np.abs(mu[near])
            rows[f] = {
                "z_fit": e["z_fit"][f],
                "n_catalogued": len(sub),
                "n_predicted": len(pred),
                "catalogued_abs_magnification": [float(m) for m in mu_cat],
                "unmatched": [
                    {
                        "ra": float(r),
                        "dec": float(d),
                        "magnification": float(m),
                        "arcsec_to_nearest_catalogued": float(
                            np.min(np.hypot(ox - pxy[j, 0], oy - pxy[j, 1]))
                        ),
                        "arcsec_to_cluster_centre": _sep(r, d, *CENTRE),
                        "abs_mu_ratio_catalogued_over_this": [float(c / abs(m)) for c in mu_cat],
                    }
                    for r, d, m, j in zip(ra, dec, mu[k], k, strict=True)
                ],
            }
        out["models"][model] = rows
    # spread across models of family 4's brightest positive-parity unmatched prediction
    pts = {}
    for m, r in out["models"].items():
        pos = [u for u in r.get("4", {}).get("unmatched", []) if u["magnification"] > 0]
        if pos:
            u = max(pos, key=lambda u: u["magnification"])
            pts[m] = (u["ra"], u["dec"])
    names = list(pts)
    out["family_4_brightest_positive_unmatched_separations_arcsec"] = {
        f"{a}|{b}": _sep(*pts[a], *pts[b]) for i, a in enumerate(names) for b in names[i + 1 :]
    }
    return out


def _sep(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    from astropy.coordinates import SkyCoord

    return float(SkyCoord(ra1, dec1, unit="deg").separation(SkyCoord(ra2, dec2, unit="deg")).arcsec)


def build_images(tex_path: Path) -> Table:
    raw = tex_path.read_bytes()
    t = srh.parse_block(raw.decode(), SIDEHEAD)
    fam = np.asarray(t["system"], str)
    t["cluster"] = CLUSTER
    t["new_in_jwst"] = np.isin(fam, HELD)
    t["use_spec"] = np.asarray(t["z_spec"], bool) & np.isin(fam, SPEC)
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        "Cerny et al. 2026, ApJ 1001, 60 (SLICE; arXiv:2503.17498v2), Table 8 block of "
        f"RX J0152.7-1357 (main.tex sha256 {hashlib.sha256(raw).hexdigest()}; e-print tarball "
        f"sha256 {srh.EPRINT_SHA256}); clump 1 of each source only; z_spec True: spectroscopic, "
        "False: SLICE model-optimized; new_in_jwst and use_spec per the paper's text (PLAN.md)"
    )
    return t


def load() -> Table:
    t = Table.read(IMAGES, format="ascii.ecsv")
    t["image_id"] = np.asarray(t["image_id"], str)
    t["system"] = np.asarray(t["system"], str)
    t["held_out"] = np.asarray(t["new_in_jwst"], bool)
    if set(t["system"][t["held_out"]]) != set(HELD):
        raise ValueError("stale image list (new_in_jwst vs HELD); rebuild it")
    t["z_spec"] = np.asarray(t["use_spec"], bool)
    t["a"] = np.ones(len(t))
    return t


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid-step", type=float, default=0.25)
    ap.add_argument(
        "--build-images", type=Path, help="SLICE e-print main.tex: write the image list"
    )
    ap.add_argument("--posthoc", action="store_true", help="write posthoc.json from summary.json")
    args = ap.parse_args(argv)
    if args.posthoc:
        ph = posthoc(args.grid_step)
        (OUT / "posthoc.json").write_text(json.dumps(ph, indent=1) + "\n")
        print(json.dumps(ph, indent=1))
        return 0
    if args.build_images:
        t = build_images(args.build_images)
        t.write(IMAGES, format="ascii.ecsv", overwrite=True)
        print(len(t), "images ->", IMAGES)
        return 0
    t = load()
    images = Table.read(IMAGES, format="ascii.ecsv")
    held = np.asarray(t["held_out"], bool)
    spec = np.asarray(t["z_spec"], bool)
    ids = np.asarray(t["image_id"], str)
    summary: dict = {
        "plan": str(OUT / "PLAN.md"),
        "rules": {
            "K0_MAX_RMS_ARCSEC": ho.K0_MAX_RMS_ARCSEC,
            "K0_MIN_MATCHED": ho.K0_MIN_MATCHED,
            "MATCH_ARCSEC": ho.MATCH_ARCSEC,
            "FRAME_MAX_ARCSEC": ho.FRAME_MAX_ARCSEC,
            "z_grid": [ho.Z_STEPS, ho.Z_MIN_OFFSET, ho.Z_MAX],
            "grid_step": args.grid_step,
            "near_one": sg.NEAR_ONE,
            "lenstool_like_max": sg.LENSTOOL_LIKE_MAX,
            "p_ratio": sg.P_RATIO,
        },
        "models": {},
    }
    d_held, d_spec, usable, classes = {}, {}, [], {}
    for model in MODELS:
        validate = json.loads((OUT / "validate" / model / "validate.json").read_text())
        p_ok = sg.gate_p(validate)
        entry, ip, ds = h29.run_model(model, t, args.grid_step)
        ok = entry["usable"] and p_ok
        d = np.asarray(ip["dtheta_arcsec"], float)
        d_held[model], d_spec[model] = d[held], ds[spec]
        if ok:
            usable.append(model)
        amp = srh.posthoc_scale({CLUSTER: entry}, images, {CLUSTER: model})[CLUSTER]
        factors = [r["scale"] for r in amp["families"].values()]
        med = float(np.median(factors)) if factors else None
        classes[model] = sg.amplitude_class(med) if ok and med is not None else None
        entry.update(
            {
                "gate_P": {"passes": p_ok, "validate": str(OUT / "validate" / model)},
                "usable_with_P": ok,
                "misses_A": per_model_misses(ids[held], d[held]),
                "misses_B": per_model_misses(ids[spec], ds[spec]),
                "amplitude": {
                    **amp,
                    "median": med,
                    "min": min(factors, default=None),
                    "max": max(factors, default=None),
                    "class": classes[model],
                },
            }
        )
        summary["models"][model] = entry
        ip.meta["model"] = model
        lc._write(ip, OUT / f"imageplane_{model}.ecsv")
    ra = ho.readings(d_held, usable)
    rb = h29.label_b(ho.readings(d_spec, usable))
    summary.update(
        {
            "images": len(t),
            "held_out_images": int(held.sum()),
            "spec_z_images": int(spec.sum()),
            "usable": usable,
            "flags_A": [str(i) for i, v in zip(ids[held], ra, strict=True) if v != "no_flag"],
            "flags_B": [str(i) for i, v in zip(ids[spec], rb, strict=True) if v != "no_flag"],
            "amplitude_classes": classes,
            "q1_answer": q1_answer(classes),
        }
    )
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    for m, e in summary["models"].items():
        print(m, e["usable_with_P"], e["gate_P"]["passes"], e["gate_K0"], e["frame_F"], sep="\n  ")
        print("  A", e["held_out_A"], e["misses_A"], "\n  B", e["spec_z_B"], e["misses_B"])
        print("  extra", e["extra_predicted"], "fewer", e["fewer_predicted"])
        print("  amplitude", {k: e["amplitude"][k] for k in ("median", "min", "max", "class")})
    print("joint", summary["flags_A"], summary["flags_B"], summary["q1_answer"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
