"""Post hoc checks for D-034 addendum bg (not in the frozen PLAN.md): gate V misses, gate C colours.

Reads the cutouts saved by ``rxj0152_fam4_photometry_v2.py --cutout-dir DIR`` (not stored), writes
``posthoc.json``; ``--sheet PNG`` also writes a contact sheet (not stored).

    python scripts/rxj0152_fam4_photometry_v2_posthoc.py DIR [--sheet PNG]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

_P = Path(__file__).resolve().with_name("rxj0152_fam4_photometry_v2.py")
_spec = importlib.util.spec_from_file_location("rxj0152_fam4_photometry_v2", _P)
v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v2)


def load(cdir: Path) -> dict:
    out = {}
    for b in v2.BANDS:
        (f,) = (cdir / f"jw06882-o008_t008_nircam_clear-{b}_i2d").glob("*.fits")
        with fits.open(f) as h:
            out[b] = (v2.v1.background(np.array(h["SCI"].data, float)), WCS(h["SCI"].header))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cutout_dir", type=Path)
    ap.add_argument("--sheet", type=Path, default=None)
    args = ap.parse_args()
    images = load(args.cutout_dir)
    cat = Table.read(v2.v1.IMAGES)
    pos = {str(r["image_id"]): (float(r["ra"]), float(r["dec"])) for r in cat}
    v2._POS41 = pos["4.1"]
    res = v2.measure(images, pos)
    frozen = json.loads((v2.OUT / "result.json").read_text())
    misses = {}
    for d in frozen["gate_V_detail"]:
        if d["recovered"]:
            continue
        t = (d["ra"], d["dec"])
        r2 = v2.measure(
            v2.inject(images, res["segmaps"], res["near41"], t, 1.0 / v2.RATIO41[d["model"]]),
            pos,
            res["sigma"],
        )
        near = [
            {
                "label": s["label"],
                "sep_arcsec": v2._sep((s["ra"], s["dec"]), t),
                "f277w": s["fluxes"]["f277w"],
                "expected_f277w": res["f41"]["f277w"] / v2.RATIO41[d["model"]],
                "in_flux_window": s["in_flux_window"],
                "colour_match": s["colour_match"],
            }
            for s in r2["region"]
            if v2._sep((s["ra"], s["dec"]), t) < 1.5
        ]
        misses[f"{d['model']}_{d['k']}"] = near
    # gate C: colour offsets per band for each catalogued image vs its siblings
    ids = [i for i in pos if i.split(".")[0] in ("2", "3", "5", "6")]
    fc = dict(zip(ids, v2.forced(images, [pos[i] for i in ids]), strict=True))
    worst = {}
    for i in ids:
        sib = [j for j in ids if j != i and j.split(".")[0] == i.split(".")[0]]
        best = None
        for j in sib:
            dmax = max(
                abs(
                    np.log10(max(fc[i][b], 1e-9) / fc[i]["f277w"])
                    - np.log10(fc[j][b] / fc[j]["f277w"])
                )
                for b in v2.BANDS
                if b != "f277w" and fc[j][b] > 0
            )
            best = dmax if best is None else min(best, dmax)
        worst[i] = {"f277w": fc[i]["f277w"], "min_over_siblings_max_abs_dex": best}
    out = {
        "note": "post hoc, not in PLAN.md",
        "gate_V_misses_segments_within_1.5arcsec": misses,
        "gate_C_colour_offsets": worst,
    }
    (v2.OUT / "posthoc.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))
    if args.sheet:
        sheet(images, res, pos, args.sheet)


def sheet(images: dict, res: dict, pos: dict, path: Path) -> None:
    import astropy.units as u
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.coordinates import SkyCoord
    from astropy.nddata import Cutout2D

    centre = SkyCoord(28.1810, -13.9545, unit="deg")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
    for ax, b in zip(axes, ("f150w", "f277w", "f444w"), strict=True):
        sci, wcs = images[b]
        c = Cutout2D(sci, centre, (20 * u.arcsec, 26 * u.arcsec), wcs=wcs)
        lim = np.nanpercentile(c.data, [5, 99.5])
        ax.imshow(np.arcsinh(c.data / lim[1] * 10), origin="lower", cmap="gray")
        for name, p in v2.PREDICTED.items():
            x, y = c.wcs.world_to_pixel(SkyCoord(*p, unit="deg"))
            r = 3.0 / c.wcs.proj_plane_pixel_scales()[0].to("arcsec").value
            ax.add_patch(plt.Circle((x, y), r, fill=False, color="c", lw=0.8))
            ax.text(x, y, name[0].upper(), color="c", ha="center", va="center")
        for s in res["region"]:
            x, y = c.wcs.world_to_pixel(SkyCoord(s["ra"], s["dec"], unit="deg"))
            ax.text(x + 4, y + 4, str(s["label"]), color="y" if s["passes"] else "r", fontsize=8)
        for i in ("4.1", "4.2"):
            x, y = c.wcs.world_to_pixel(SkyCoord(*pos[i], unit="deg"))
            ax.text(x, y, i, color="m", fontsize=8)
        ax.set_title(b)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=110)


if __name__ == "__main__":
    main()
