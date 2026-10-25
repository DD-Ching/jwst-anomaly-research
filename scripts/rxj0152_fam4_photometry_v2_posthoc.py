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

import astropy.units as u
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
    res = v2.measure(images, pos)
    frozen = json.loads((v2.OUT / "result.json").read_text())
    pre = [(c["ra"], c["dec"]) for c in res["counterparts"]]
    misses, credited = {}, {}
    for d in frozen["gate_V_detail"]:
        t = (d["ra"], d["dec"])
        scale = 1.0 / v2.RATIO41[d["model"]]
        r2 = v2.measure(
            v2.inject(images, res["segmaps"], res["near41"], pos["4.1"], t, scale),
            pos,
            res["sigma"],
        )
        if d["recovered"]:
            # was the recovery the injected copy, or a counterpart already there without injection?
            hits = [c for c in r2["counterparts"] if v2._sep((c["ra"], c["dec"]), t) < 1.0]
            new = [c for c in hits if min(v2._sep((c["ra"], c["dec"]), q) for q in pre) > 0.3]
            credited[f"{d['model']}_{d['k']}"] = {
                "n_hits": len(hits),
                "n_hits_not_preexisting": len(new),
            }
            continue
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
    n_pre_only = sum(v["n_hits_not_preexisting"] == 0 for v in credited.values())
    # gate C: per sibling, rule 6 detectable bands, non-positive fluxes, largest offset
    ids = [i for i in pos if i.split(".")[0] in ("2", "3", "5", "6")]
    fc = dict(zip(ids, v2.forced(images, [pos[i] for i in ids]), strict=True))
    sig = res["sigma"]
    worst = {}
    for i in ids:
        rows = {}
        for j in (j for j in ids if j != i and j.split(".")[0] == i.split(".")[0]):
            used = [
                b
                for b in v2.BANDS
                if b != "f277w"
                and fc[j][b] > 0
                and fc[j][b] * fc[i]["f277w"] / fc[j]["f277w"] >= v2.DETECTABLE_SIGMA * sig[b]
            ]
            nonpos = [b for b in used if fc[i][b] <= 0]
            offs = [
                abs(np.log10(fc[i][b] / fc[i]["f277w"]) - np.log10(fc[j][b] / fc[j]["f277w"]))
                for b in used
                if b not in nonpos
            ]
            rows[j] = {
                "n_detectable_bands": len(used),
                "nonpositive_bands": nonpos,
                "max_abs_dex": max(offs) if offs else None,
            }
        worst[i] = {"f277w": fc[i]["f277w"], "vs_siblings": rows}
    # NaN pixels inside the 0.3" apertures (forced() treats them as 0)
    from astropy.coordinates import SkyCoord
    from photutils.aperture import SkyCircularAperture

    pts = {"4.1": pos["4.1"], **{str(s["label"]): (s["ra"], s["dec"]) for s in res["region"]}}
    pts.update({i: pos[i] for i in ids})
    nan_hits = []
    for b in v2.BANDS:
        sci, wcs = images[b]
        for k, p in pts.items():
            ap = SkyCircularAperture(SkyCoord(*p, unit="deg"), r=v2.APER_RADIUS * u.arcsec)
            mask = ap.to_pixel(wcs).to_mask(method="exact")
            cut = mask.cutout(sci, fill_value=np.nan)
            if cut is None or np.any(~np.isfinite(cut[mask.data > 0])):
                nan_hits.append([k, b])
    out = {
        "note": "post hoc, not in PLAN.md",
        "gate_V_misses_segments_within_1.5arcsec": misses,
        "gate_V_recovered_credit": credited,
        "gate_V_recovered_only_by_preexisting_counterpart": n_pre_only,
        "aperture_nonfinite_pixels": nan_hits,
        "gate_C_colour_offsets": worst,
    }
    (v2.OUT / "posthoc.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))
    if args.sheet:
        sheet(images, res, pos, args.sheet)


def sheet(images: dict, res: dict, pos: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.coordinates import SkyCoord
    from astropy.nddata import Cutout2D

    centre = SkyCoord(*v2.v1.CENTRE, unit="deg")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
    for ax, b in zip(axes, ("f150w", "f277w", "f444w"), strict=True):
        sci, wcs = images[b]
        c = Cutout2D(sci, centre, (20 * u.arcsec, 26 * u.arcsec), wcs=wcs)
        lim = np.nanpercentile(c.data, [5, 99.5])
        ax.imshow(np.arcsinh(c.data / lim[1] * 10), origin="lower", cmap="gray")
        for name, p in v2.PREDICTED.items():
            x, y = c.wcs.world_to_pixel(SkyCoord(*p, unit="deg"))
            r = v2.REGION_RADIUS / c.wcs.proj_plane_pixel_scales()[0].to("arcsec").value
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
