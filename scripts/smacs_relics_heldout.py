"""RELICS HST-era SMACS 0723 models vs the JWST-era Mahler+22 ICLv2 images (held-out stage).

Frozen plan: ``results/lens_consistency/smacs_relics_heldout/PLAN.md``. Rules 1-5 of the
``mars_hff_heldout`` stage, unchanged, per RELICS model (glafic v2, zitrin-ltm-gauss v2,
lenstool v1), with the Golubchik+22 HST-era images as the "seen" list. No flag, flux or
extra-image reading.

    python scripts/smacs_relics_heldout.py [--out DIR] [--keep-maps]

The maps are deleted after use unless ``--keep-maps`` (cloud disk rule).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import m0417_images as mi
import mars_hff_heldout as mo
import mars_hff_inputs as mh
import numpy as np
import smacs_s17_vet as sv
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel
from jwst_anomaly.photometry import fetch_catalog

OUT = Path("results/lens_consistency/smacs_relics_heldout")
SEEN = OUT / "golubchik22_images.ecsv"
ICLV2 = "smacs0723-iclv2"
# name -> (HLSP subdirectory, sha256 x, sha256 y); arcsec deflections at D_LS/D_S = 1
MODELS = {
    **sv.RELICS,
    "lenstool_v1": (  # sha256 as downloaded 2026-10-25
        "lenstool/",
        "e3ef2551e11c218e6760bf826bfc6172b2e1e1aeac091d06e68c710c397959dc",
        "c7a51d2af371657ed4f610c61188948903576487d54e70ab7b519ed520d95b9a",
    ),
}


def heldout_images() -> Table:
    """The ICLv2 images with z (PLAN.md Inputs): spectroscopic, else best.par's optimised z."""
    _, files, par = lc.load_model(ICLV2)
    images, z = lc.image_list(ICLV2, files, par)
    t = images.copy()
    fam = [lensmodel.image_family(i) for i in t["image_id"]]
    listed = np.asarray(t["z"], float)
    t["z"] = [zl if zl > 0 else float(z.get(f, 0.0)) for zl, f in zip(listed, fam, strict=True)]
    return t


def fetch_maps(name: str) -> tuple[Path, Path]:
    """The pinned x / y arcsec deflection maps of one RELICS model (checksummed downloads)."""
    sub, sha_x, sha_y = MODELS[name]
    base = f"{sv.RELICS_URL}{sub}hlsp_relics_model_model_smacs0723m73_{name}"
    px = fetch_catalog(f"{base}_x-arcsec-deflect.fits", sha_x)
    try:
        return px, fetch_catalog(f"{base}_y-arcsec-deflect.fits", sha_y)
    except Exception:
        px.unlink(missing_ok=True)
        raise


def load_relics(name: str, centre: tuple[float, float], keep: bool = False):
    """One RELICS map model (as ``smacs_s17_vet.relics``, without its frame shift: rule 1 moves
    the image list into the RELICS frame instead)."""
    px, py = fetch_maps(name)
    try:
        h0, om = sv.RELICS_COSMO
        return lensmodel.MapLensModel.from_fits(
            px, py, sv.RELICS_Z_LENS, FlatLambdaCDM(H0=h0, Om0=om), source=name, centre=centre
        )
    finally:
        if not keep:
            px.unlink(missing_ok=True)
            py.unlink(missing_ok=True)


def seen_rms(model, grid, t: Table) -> dict:
    """Reported, not a gate (mars_hff_heldout rule 4's training comparison): residuals of the
    seen images of each family with >= 2 of them, from the mean source of those images."""
    seen = np.asarray(t["seen"], bool)
    s = t[seen]
    bt = lensmodel.backtrace_images(model, s, {})
    ox, oy = model.to_frame(np.asarray(s["ra"], float), np.asarray(s["dec"], float))
    sysv = np.asarray(s["system"], str)
    d = []
    for fam in dict.fromkeys(sysv):
        m = sysv == fam
        bx, by = np.asarray(bt["beta_x"])[m], np.asarray(bt["beta_y"])[m]
        if m.sum() < 2 or not (np.isfinite(bx).all() and np.isfinite(by).all()):
            continue
        z = float(s["z"][np.where(m)[0][0]])
        pred = lensmodel.find_images(model, grid, float(bx.mean()), float(by.mean()), z)
        for i in np.where(m)[0]:
            d.append(np.hypot(ox[i] - pred["x"], oy[i] - pred["y"]).min() if len(pred) else np.inf)
    d = np.asarray(d, float)
    return {
        "n": len(d),
        "rms_arcsec": mi.rms_of(d),
        "matched_fraction": float(np.mean(d <= mi.MATCH_ARCSEC)) if len(d) else 0.0,
    }


def model_block(name: str, images: Table, seen: Table, out: Path, keep: bool) -> dict:
    off = mo.frame_offset(images, seen)
    t = mo.split_seen(mo.shift(images, off), seen)
    centre = (float(np.mean(seen["ra"])), float(np.mean(seen["dec"])))
    model = load_relics(name, centre, keep)
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, t), mh.GRID_STEP)
    r = mo.heldout_residuals(model, grid, t)
    r["dtheta_arcsec"].info.format = ".3f"
    r.write(out / f"heldout_{name}.ecsv", format="ascii.ecsv", overwrite=True)
    g = mo.gate_h(r)
    held = ~np.asarray(t["seen"], bool)
    tested = set(np.asarray(r["image_id"], str)) if len(r) else set()
    d = np.asarray(r["dtheta_arcsec"], float) if len(r) else np.array([])
    typ = np.asarray(r["type"], str) if len(r) else np.array([], str)
    return {
        "frame_offset": off,
        "seen": int((~held).sum()),
        "seen_systems": sorted(set(np.asarray(t["system"], str)[~held]), key=float),
        "held_out": int(held.sum()),
        "held_out_untested": [str(i) for i in t["image_id"][held] if str(i) not in tested],
        "types": {k: int(np.sum(typ == k)) for k in "ENM"},
        "outside_map": int(np.sum(~np.asarray(r["inside"], bool))) if len(r) else 0,
        "gate_H": g,
        "gate_H_type_E": mo.gate_h(r[typ == "E"]) if len(r) else None,
        "outcome": mo.outcome(g),
        "seen_rms_reported": seen_rms(model, grid, t),
        "over_3arcsec": [str(i) for i in r["image_id"][~(d <= 3.0)]] if len(r) else [],
    }


def run(out: Path = OUT, keep: bool = False) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    images = heldout_images()
    seen = Table.read(SEEN)
    summary = {
        "plan": str(out / "PLAN.md"),
        "held_out_list": f"{ICLV2} arcs.dat ({len(images)} images)",
        "seen_list": f"{SEEN} ({len(seen)} rows)",
        "z_used": {
            lensmodel.image_family(i): float(z)
            for i, z in zip(images["image_id"], images["z"], strict=True)
        },
        "provenance": "positions: observed (Mahler+22, Golubchik+22); z: observed (spec) or "
        "model value (ICLv2 best.par); residuals: model_prediction; thresholds: ASSUMPTION "
        "(mars_hff_heldout/PLAN.md)",
        "models": {n: model_block(n, images, seen, out, keep) for n in MODELS},
    }
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--keep-maps", action="store_true")
    a = ap.parse_args(argv)
    s = run(a.out, a.keep_maps)
    print(json.dumps(_clean(s["models"]), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
