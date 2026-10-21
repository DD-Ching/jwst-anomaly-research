"""MACS0717: known-case gates for the HFF Sharon v4cor and Keeton v4 map models.

Same gates as ``m1149_sharon`` (G1: κ and z = 2 |μ| derived from the deflection maps against the
published maps; G2: the readme Table 7 spectroscopic images back-traced and re-solved), applied to
``macs0717-sharon`` and ``macs0717-keeton``; the same images go through ``macs0717-cats`` for
comparison. Frozen plan: ``results/lens_consistency/macs0717_models/PLAN.md``.

    python scripts/m0717_models.py [--out results/lens_consistency/macs0717_models]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import lens_consistency as lc
import m1149_sharon as ms
import numpy as np

OUT = Path("results/lens_consistency/macs0717_models")
GATED = ("macs0717-sharon", "macs0717-keeton")
COMPARISON = "macs0717-cats"


def rms(ip) -> tuple[float | None, int]:
    d = np.asarray(ip["dtheta_arcsec"], float)
    ok = np.isfinite(d)
    return (float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else None), int(ok.sum())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--grid-step", type=float, default=0.25)
    ap.add_argument("--map-step", type=int, default=37)
    args = ap.parse_args(argv)
    images = ms.gold_images(args.out / "sharon_gold_images.ecsv")
    summary: dict = {"images": len(images)}
    for name in (*GATED, COMPARISON):
        model, _, _ = lc.load_model(name)
        entry: dict = {"model_sha256": model.sha256}
        ip = ms.image_plane(name, images, args.grid_step)
        lc._write(ip, args.out / f"imageplane_{name}.ecsv")
        entry["rms_arcsec"], entry["n_solved"] = rms(ip)
        inside = ms.footprint(name, images)
        entry["images_outside_maps"] = [str(i) for i in images["image_id"][~inside]]
        if name in GATED:
            files = lc.model_files(name)
            maps = {
                "kappa": lc.map_check(model, files["kappa_map"], "kappa", None, args.map_step),
                "mu_z2": lc.map_check(model, files["mag_z2"], "mu", 2.0, args.map_step),
            }
            for m in maps.values():
                m["map"] = Path(m["map"]).name  # cache file name (sha256 prefix), not a local path
            entry["maps"] = maps
            entry["gates"] = ms.gates(maps, ip)
            # post hoc (deviation from PLAN.md): G2 on the systems with every image inside the maps
            outside = sorted({str(s) for s in images["system"][~inside]})
            covered = ~np.isin(np.asarray(ip["system"], str), outside)
            entry["post_hoc"] = {
                "systems_outside_maps": outside,
                "G2_inside_footprint": ms.gates(maps, ip[covered])["G2"],
            }
        summary[name] = entry
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({n: summary[n]["gates"] for n in GATED}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
