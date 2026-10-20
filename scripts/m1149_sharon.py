"""MACS1149: known-case gates for the HFF Sharon v4cor map model (``macs1149-sharon``).

G1 compares the model's κ and z = 2 |μ| (derived from the deflection maps) with the published maps;
G2 back-traces the model's own gold spectroscopic images (readme Table 5) and re-solves them.
The same images go through ``macs1149-cats`` for comparison. Frozen plan:
``results/lens_consistency/macs1149_sharon/PLAN.md``.

    python scripts/m1149_sharon.py [--out results/lens_consistency/macs1149_sharon]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import numpy as np
from astropy.table import Table

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/macs1149_sharon")
# Frozen gates (PLAN.md; ASSUMPTIONs)
G1_MAX_MEDIAN_REL = 0.05  # median relative difference of κ and of |μ(z = 2)| vs the published maps
G2_MAX_RMS_ARCSEC = 1.0  # image-plane rms of the gold images
G2_MAX_MATCH_ARCSEC = 3.0  # every catalogued image matched by a predicted image within this
# Post hoc (deviation from PLAN.md): system 13 is free in CATS v4.1 (params.txt z_m_limit 13: 1.301)
CATS_FITTED_Z = {"13": 1.301}


def gold_images(path: Path) -> Table:
    """The transcribed gold images as a :func:`lensmodel.load_lenstool_images`-like table."""
    t = Table.read(path, format="ascii.ecsv")
    t["image_id"] = np.asarray(t["image_id"], str)
    t["system"] = [lensmodel.image_family(i) for i in t["image_id"]]
    t["a"] = np.ones(len(t))
    t.meta.setdefault("source", str(path))
    return t


def footprint(name: str, images: Table) -> np.ndarray:
    """Whether each image lies inside the model's deflection maps (finite deflection)."""
    model, _, _ = lc.load_model(name)
    x, y = model.to_frame(np.asarray(images["ra"], float), np.asarray(images["dec"], float))
    lm = model.lens_map(x, y, np.asarray(images["z"], float))
    return np.isfinite(np.asarray(lm["beta_x"], float))


def gates(maps: dict, ip: Table) -> dict:
    """Apply G1 and G2 (frozen) to the map checks and the image-plane table."""
    g1 = {k: v["median_rel_diff"] for k, v in maps.items()}
    g1_pass = all(np.isfinite(v) and v <= G1_MAX_MEDIAN_REL for v in g1.values())
    d = np.asarray(ip["dtheta_arcsec"], float)
    systems = np.asarray(ip["system"], str)
    ok = np.isfinite(d)
    rms = float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else float("nan")
    per_sys = {}
    for s in dict.fromkeys(systems):
        m = systems == s
        per_sys[s] = bool(np.all(np.isfinite(d[m])) and np.all(d[m] <= G2_MAX_MATCH_ARCSEC))
    g2_pass = bool(ok.all() and rms <= G2_MAX_RMS_ARCSEC and all(per_sys.values()))
    return {
        "G1": {"median_rel_diff": g1, "max": G1_MAX_MEDIAN_REL, "passes": bool(g1_pass)},
        "G2": {
            "rms_dtheta_arcsec": rms,
            "max_rms": G2_MAX_RMS_ARCSEC,
            "systems_reproduced": per_sys,
            "max_match_arcsec": G2_MAX_MATCH_ARCSEC,
            "passes": g2_pass,
        },
        "pinned_for_use": bool(g1_pass and g2_pass),
    }


def image_plane(name: str, images: Table, grid_step: float) -> Table:
    model, _, _ = lc.load_model(name)
    lc.apply_frame_offset(name, model)
    bt = lensmodel.backtrace_images(model, images, {})
    half = lc.grid_half_width(model, images)
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, half, grid_step), half, grid_step
    )
    ip = lensmodel.imageplane_residuals(model, grid, bt)
    ip.meta["model"] = name
    return ip


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--grid-step", type=float, default=0.25)
    ap.add_argument("--map-step", type=int, default=37)
    args = ap.parse_args(argv)
    images = gold_images(args.out / "sharon_gold_images.ecsv")
    name = "macs1149-sharon"
    model, _, _ = lc.load_model(name)
    files = lc.model_files(name)
    maps = {
        "kappa": lc.map_check(model, files["kappa_map"], "kappa", None, args.map_step),
        "mu_z2": lc.map_check(model, files["mag_z2"], "mu", 2.0, args.map_step),
    }
    for m in maps.values():
        m["map"] = Path(m["map"]).name  # cache file name (sha256 prefix), not a local path
    summary = {"model": name, "model_sha256": model.sha256, "maps": maps}
    tables = {m: image_plane(m, images, args.grid_step) for m in (name, "macs1149-cats")}
    summary["gates"] = gates(maps, tables[name])
    # post hoc: which images the published maps cover, and G2 on the covered systems only
    inside = footprint(name, images)
    summary["post_hoc"] = {
        "images_outside_maps": [str(i) for i in images["image_id"][~inside]],
        "systems_outside_maps": sorted({str(s) for s in images["system"][~inside]}),
    }
    covered = ~np.isin(
        np.asarray(tables[name]["system"], str), summary["post_hoc"]["systems_outside_maps"]
    )
    summary["post_hoc"]["G2_inside_footprint"] = gates(maps, tables[name][covered])["G2"]
    refit = images.copy()
    for s_id, z in CATS_FITTED_Z.items():
        refit["z"][np.asarray(refit["system"], str) == s_id] = z
    tables["macs1149-cats-fittedz"] = image_plane("macs1149-cats", refit, args.grid_step)
    for m, ip in tables.items():
        lc._write(ip, args.out / f"imageplane_{m}.ecsv")
        d = np.asarray(ip["dtheta_arcsec"], float)
        ok = np.isfinite(d)
        summary[f"rms_{m}"] = float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else None
        summary[f"n_solved_{m}"] = int(ok.sum())
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["gates"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
