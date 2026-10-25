"""MACS J2129: the HST-era CLASH Zitrin v2 maps on the SLICE JWST-era image list.

Frozen plan: ``results/lens_consistency/macs2129_clash_heldout/PLAN.md``. Images: SLICE Table 8
(``slice_images.ecsv``, observed). Per model: each family at its fitted z (rule 1), gate K0 and
frame check F on families 1-14 (rules 2-3), held-out reading A on the SLICE-new families 16-18
(rule 4) and reading B on the spectroscopic families at their spectroscopic z (rule 5). Rule
functions are those of ``macs0647_clash_heldout`` (same frozen rules). Writes ``summary.json``
and per-model image-plane tables.

    python scripts/macs2129_clash_heldout.py [--grid-step 0.25]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import m0417_images as mi
import macs0647_clash_heldout as ho
import numpy as np
from astropy.table import Table

OUT = Path("results/lens_consistency/macs2129_clash_heldout")
IMAGES = OUT / "slice_images.ecsv"
MODELS = ("macs2129-clash-nfw", "macs2129-clash-ltm")
NEW_FAMILIES = ("16", "17", "18")  # SLICE sources first identified in JWST imaging (paper text)


def load_images(path: Path = IMAGES) -> Table:
    """All 53 rows; ``held_out`` = families 16-18 (rule 4); ``z_spec`` rows feed reading B."""
    t = Table.read(path, format="ascii.ecsv")
    t["image_id"] = np.asarray(t["image_id"], str)
    t["system"] = np.asarray(t["system"], str)
    t["held_out"] = np.asarray(t["new_in_jwst"], bool)
    if set(t["system"][t["held_out"]]) != set(NEW_FAMILIES):
        raise ValueError("new_in_jwst disagrees with the frozen held-out families 16-18")
    t["a"] = np.ones(len(t))
    return t


def label_b(r: np.ndarray) -> np.ndarray:
    """Rule 5 label: the rule-4 decision under its own name."""
    return np.where(np.asarray(r) == "lensing_residual_input", "redshift_residual_input", r)


def reading_b(d_spec: dict[str, np.ndarray], spec: np.ndarray, usable: list[str]) -> np.ndarray:
    """Rule 5: the rule-4 decision on the z_spec offsets of the spectroscopic rows only."""
    spec = np.asarray(spec, bool)
    return label_b(ho.readings({n: np.asarray(d, float)[spec] for n, d in d_spec.items()}, usable))


def run_model(name: str, images: Table, grid_step: float) -> tuple[dict, Table, np.ndarray]:
    mi.Z_STEPS, mi.Z_MIN_OFFSET, mi.Z_MAX, mi.MATCH_ARCSEC = (
        ho.Z_STEPS,
        ho.Z_MIN_OFFSET,
        ho.Z_MAX,
        ho.MATCH_ARCSEC,
    )
    solver = mi.Solver(name, images, grid_step)
    ip = mi.fitted_z(solver, images)
    order = {i: k for k, i in enumerate(np.asarray(ip["image_id"], str))}
    ip = ip[[order[i] for i in np.asarray(images["image_id"], str)]]
    held = np.asarray(images["held_out"], bool)
    d = np.asarray(ip["dtheta_arcsec"], float)
    dra, ddec = ho.signed_offsets(solver, ip, images)
    ip["dra_cosdec_arcsec"] = dra
    ip["ddec_arcsec"] = ddec
    ip["held_out"] = held
    fam = np.asarray(images["system"], str)
    z_fit = {f: float(ip["z_fit"][fam == f][0]) for f in dict.fromkeys(fam)}
    counts = {f: (int(np.sum(fam == f)), int(np.max(ip["n_predicted"][fam == f]))) for f in z_fit}
    # rule 5: spectroscopic families at their spectroscopic z, no fit
    spec = np.asarray(images["z_spec"], bool)
    ts = images[spec]
    ips = solver.image_plane(ts)
    ds = np.full(len(images), np.nan)
    ds[spec] = np.asarray(ips["dtheta_arcsec"], float)
    ip["dtheta_specz_arcsec"] = ds
    entry = {
        "model_sha256": solver.model.sha256,
        "gate_K0": ho.gate_k0(d[~held]),
        "frame_F": ho.frame_check(dra[~held], ddec[~held]),
        "held_out_A": ho.gate_k0(d[held]),  # statistics only; "passes" is not a rule here
        "spec_z_B": ho.gate_k0(ds[spec]),
        "spec_fitted_z": ho.gate_k0(d[spec]),
        "z_fit": z_fit,
        "z_spec": {f: float(images["z"][fam == f][0]) for f in dict.fromkeys(fam[spec])},
        "extra_predicted": {f: c for f, c in counts.items() if c[1] > c[0]},
        "fewer_predicted": {f: c for f, c in counts.items() if c[1] < c[0]},
    }
    for k in ("held_out_A", "spec_z_B", "spec_fitted_z"):
        entry[k].pop("passes")
    entry["usable"] = entry["gate_K0"]["passes"] and entry["frame_F"]["passes"]
    return entry, ip, ds


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid-step", type=float, default=0.25)
    args = ap.parse_args(argv)
    images = load_images()
    held = np.asarray(images["held_out"], bool)
    spec = np.asarray(images["z_spec"], bool)
    summary: dict = {
        "plan": str(OUT / "PLAN.md"),
        "images": len(images),
        "held_out_images": int(held.sum()),
        "spec_z_images": int(spec.sum()),
        "rules": {
            "K0_MAX_RMS_ARCSEC": ho.K0_MAX_RMS_ARCSEC,
            "K0_MIN_MATCHED": ho.K0_MIN_MATCHED,
            "MATCH_ARCSEC": ho.MATCH_ARCSEC,
            "FRAME_MAX_ARCSEC": ho.FRAME_MAX_ARCSEC,
            "z_grid": [ho.Z_STEPS, ho.Z_MIN_OFFSET, ho.Z_MAX],
            "grid_step": args.grid_step,
        },
    }
    d_fit, d_spec = {}, {}
    for name in MODELS:
        entry, ip, ds = run_model(name, images, args.grid_step)
        summary[name] = entry
        d_fit[name] = np.asarray(ip["dtheta_arcsec"], float)
        d_spec[name] = ds
        ip.meta["model"] = name
        lc._write(ip, OUT / f"imageplane_{name}.ecsv")
    usable = [n for n in MODELS if summary[n]["usable"]]
    ra = ho.readings({n: d[held] for n, d in d_fit.items()}, usable)
    rb = reading_b(d_spec, spec, usable)
    summary["usable_models"] = usable
    summary["flags_A"] = [
        str(i) for i, v in zip(images["image_id"][held], ra, strict=True) if v != "no_flag"
    ]
    summary["flags_B"] = [
        str(i) for i, v in zip(images["image_id"][spec], rb, strict=True) if v != "no_flag"
    ]
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("usable_models", "flags_A", "flags_B")}, indent=1))
    for n in MODELS:
        e = summary[n]
        print(n, e["gate_K0"], e["frame_F"], e["held_out_A"], e["spec_z_B"], sep="\n  ")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
