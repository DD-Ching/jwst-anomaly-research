"""Full-chain check of arc difference imaging: detection -> morphology -> arc mask (2026-10-16).

Implements ``results/transients/arc_chain/PLAN.md`` (written before the run). The arc mask is a
reproducible, image-based definition of "on the arc": the epoch-1 image of the reddest band of an
``arc_diff.py --maps`` file, minus a 3" median filter, smoothed (sigma 0.06"), thresholded at 3x
the robust sky sigma; the 8-connected components that contain a catalogued multiple-image position
(within 0.3"), dilated by 0.3". Every number is ``derived``; the parameters are ASSUMPTIONs fixed in
the plan.

Subcommands:
- ``known``: per known signal, the stage reached (hit, morph pass, in mask).
- ``inject``: for an ``--inject-grid`` run, injected stars inside the mask that survive each stage.
- ``review`` (owner decision 2026-10-16, ``results/transients/review_list/PLAN.md``): position is a
  label, not a stage; every unflagged morph ``pass`` is listed with position labels and a priority.
- ``inject-review``: whole-box injection efficiency of that chain, by mask and arc-level stratum.

    python scripts/arc_chain.py known --maps m.npz --hits hits.ecsv --morph morph.json \\
        --known known.csv --field smacs --out chain.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parent))

MEDIAN_ARCSEC = 3.0
SMOOTH_ARCSEC = 0.06
THRESH_SIGMA = 3.0
SEED_ARCSEC = 0.3
DILATE_ARCSEC = 0.3
MAX_AREA_FRAC = 0.25
MATCH_ARCSEC = 0.15  # a known signal is "a hit" if a hit lies within this (arc_diff default)
BAND_ORDER = ("F090W", "F115W", "F150W", "F200W", "F277W", "F356W", "F410M", "F444W")


def image_positions(field: str) -> tuple[np.ndarray, np.ndarray]:
    """Catalogued multiple-image positions (RA, Dec in degrees) for a field (PLAN.md lists)."""
    if field == "sunrise":
        from sunrise_scofield import IMAGES

        v = np.array([IMAGES[k][:2] for k in IMAGES], float)
        return v[:, 0], v[:, 1]
    import lens_consistency as lc

    if field == "smacs":
        import smacs_arcdiff_images as sai

        from jwst_anomaly.photometry import fetch_catalog

        url, sha = lc.MODELS[sai.MODEL]["files"]["arcs.dat"]
        t = lc.shift_images(sai.MODEL, sai.load_lenstool_images(fetch_catalog(url, sha)))
    elif field == "macs0416":
        t, _ = lc.image_list("macs0416-canucs", {}, None)
        lc.shift_images("macs0416-canucs", t)
    elif field == "elgordo":
        _, files, par = lc.load_model("elgordo-caminha23")
        t, _ = lc.image_list("elgordo-caminha23", files, par)
        lc.shift_images("elgordo-caminha23", t)
    else:
        raise ValueError(field)
    return np.asarray(t["ra"], float), np.asarray(t["dec"], float)


def filtered_snr(img: np.ndarray, scale: float) -> np.ndarray:
    """Median-subtracted, smoothed image in units of its robust sky sigma (PLAN.md steps 1-3)."""
    from scipy import ndimage

    im = np.where(np.isfinite(img), img, np.nanmedian(img))
    k = max(int(round(MEDIAN_ARCSEC / scale)) | 1, 3)
    sm = ndimage.gaussian_filter(im - ndimage.median_filter(im, size=k), SMOOTH_ARCSEC / scale)
    med = np.median(sm)
    mad = np.median(np.abs(sm - med))
    if not mad > 0:
        raise ValueError("robust sigma is zero: the box is mostly constant or blank")
    return (sm - med) / (1.4826 * mad)


def arc_mask(img: np.ndarray, scale: float, seeds_px: list[tuple[float, float]]) -> np.ndarray:
    """PLAN.md arc mask on an epoch-1 image (pixel scale ``scale`` arcsec, seeds in pixels)."""
    from scipy import ndimage

    im = img
    det = filtered_snr(img, scale) > THRESH_SIGMA
    lab, _ = ndimage.label(det, structure=np.ones((3, 3), bool))
    r = SEED_ARCSEC / scale
    keep = set()
    yy, xx = np.indices(im.shape)
    for x, y in seeds_px:
        if not (0 <= x < im.shape[1] and 0 <= y < im.shape[0]):
            continue
        near = ((xx - x) ** 2 + (yy - y) ** 2 <= r * r) & (lab > 0)
        keep.update(np.unique(lab[near]).tolist())
    mask = np.isin(lab, sorted(keep)) if keep else np.zeros(im.shape, bool)
    d = int(round(DILATE_ARCSEC / scale))
    if d > 0 and mask.any():
        disk = (np.add.outer(np.arange(-d, d + 1) ** 2, np.arange(-d, d + 1) ** 2)) <= d * d
        mask = ndimage.binary_dilation(mask, structure=disk)
    return mask


def load_mask(maps: Path, field: str, with_snr: bool = False):
    with np.load(maps) as m:
        bands = sorted({k.split("_")[0] for k in m.files}, key=BAND_ORDER.index)
        b = bands[-1]
        img = m[f"{b}_e1img"]
        w = WCS(fits.Header.fromstring(str(m[f"{b}_wcs"])))
    scale = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600
    ra, dec = image_positions(field)
    x, y = w.all_world2pix(ra, dec, 0)
    mask = arc_mask(img, scale, list(zip(x, y, strict=True)))
    if with_snr:
        return mask, w, b, scale, filtered_snr(img, scale)
    return mask, w, b, scale


def _pixel(shape, w: WCS, ra: float, dec: float) -> tuple[int, int] | None:
    """Nearest pixel (x, y) of a sky position, or None when it is off the array."""
    x, y = (float(v) for v in w.all_world2pix(ra, dec, 0))
    if not (np.isfinite(x) and np.isfinite(y)):
        return None
    x, y = int(round(x)), int(round(y))
    return (x, y) if 0 <= x < shape[1] and 0 <= y < shape[0] else None


def in_mask(mask: np.ndarray, w: WCS, ra: float, dec: float) -> bool:
    p = _pixel(mask.shape, w, ra, dec)
    return bool(p is not None and mask[p[1], p[0]])


def _sep(ra1, dec1, ra2, dec2):
    return 3600 * np.hypot((ra1 - ra2) * np.cos(np.radians(dec2)), dec1 - dec2)


def known_chain(maps, hits, morph, known, field, offset=(0.0, 0.0)) -> dict:
    from scipy import ndimage

    mask, w, band, scale, snr = load_mask(maps, field, with_snr=True)
    edge = ndimage.distance_transform_edt(~mask) * scale  # arcsec to the mask, 0 inside
    ira, idec = image_positions(field)
    h = Table.read(hits)
    verdict = {t["id"]: t["verdict"] for t in json.loads(Path(morph).read_text())["targets"]}
    from arc_diff import load_known  # same frame shift as the detection run

    ks = load_known(known, tuple(offset))
    rows = []
    for k in ks:
        ra, dec = float(k["ra"]), float(k["dec"])
        d = _sep(np.asarray(h["ra"], float), np.asarray(h["dec"], float), ra, dec)
        j = int(np.argmin(d)) if len(d) else -1
        hit = str(h["hit_id"][j]) if j >= 0 and d[j] <= MATCH_ARCSEC else None
        p = _pixel(mask.shape, w, ra, dec)
        rows.append(
            {
                "name": k["name"],
                "hit": hit,
                "near_bright_core": bool(h["near_bright_core"][j]) if hit else None,
                "morph": verdict.get(hit) if hit else None,
                "in_mask": in_mask(mask, w, ra, dec),
                "mask_distance_arcsec": round(float(edge[p[1], p[0]]), 2) if p else None,
                "nearest_image_arcsec": round(float(_sep(ira, idec, ra, dec).min()), 2),
                "arc_level_sigma": round(float(snr[p[1], p[0]]), 2) if p else None,
                "full_chain": bool(
                    hit and verdict.get(hit) == "pass" and in_mask(mask, w, ra, dec)
                ),
            }
        )
    hits_in = [
        {
            "id": str(r["hit_id"]),
            "verdict": verdict.get(str(r["hit_id"])),
            "near_bright_core": bool(r["near_bright_core"]),
        }
        for r in h
        if in_mask(mask, w, float(r["ra"]), float(r["dec"]))
    ]
    return {
        "mask_band": band,
        "mask_area_frac": float(mask.mean()),
        "mask_area_ok": bool(mask.mean() <= MAX_AREA_FRAC),
        "known": rows,
        "hits_in_mask": hits_in,
        "n_hits": len(h),
    }


def inject_chain(runs, field) -> dict:
    """Pool ``(maps, injected, morph)`` runs of one box (e.g. several ``--seed`` values)."""
    per: dict[str, dict] = {}
    masks = []
    for maps, injected, morph in runs:
        mask, w, band, _ = load_mask(maps, field)
        masks.append({"maps": Path(maps).name, "band": band, "area_frac": float(mask.mean())})
        verdict = {t["id"]: t["verdict"] for t in json.loads(Path(morph).read_text())["targets"]}
        with open(injected, newline="") as f:
            for r in csv.DictReader(f):
                if not in_mask(mask, w, float(r["ra"]), float(r["dec"])):
                    continue
                c = per.setdefault(
                    r["mag"],
                    {"n_in_mask": 0, "hit": 0, "morph_pass": 0, "chain": 0, "chain_unflagged": 0},
                )
                c["n_in_mask"] += 1
                if not r["hit"]:
                    continue
                c["hit"] += 1
                if verdict.get(r["hit"]) == "pass":
                    c["morph_pass"] += 1
                    c["chain"] += 1  # position: the star is inside the mask by selection
                    # the default search (no --all-hits) also drops near_bright_core hits
                    c["chain_unflagged"] += r["near_bright_core"] not in ("True", "true", "1")
    for c in per.values():
        n = c["n_in_mask"]
        for k in ("hit", "chain", "chain_unflagged"):
            c[f"eff_{k}"] = round(c[k] / n, 3) if n else None
    return {"masks": masks, "n_runs": len(runs), "per_mag": per}


OLD_RULE_ARCSEC = 1.0  # frozen position rule of #274 / #275, now only a label
PRIORITY_B_ARCSEC = 2.0  # review-list priority B (results/transients/review_list/PLAN.md)


def _truthy(v) -> bool:
    return str(v) in ("True", "true", "1")


def _labels(mask, w, scale, snr, edge, ira, idec, ra, dec) -> dict:
    """Position labels of one sky position (review-list PLAN.md); they never veto."""
    p = _pixel(mask.shape, w, ra, dec)
    near = float(_sep(ira, idec, ra, dec).min())
    inm = bool(p is not None and mask[p[1], p[0]])
    return {
        "in_mask": inm,
        "mask_distance_arcsec": round(float(edge[p[1], p[0]]), 2) if p else None,
        "arc_level_sigma": round(float(snr[p[1], p[0]]), 2) if p else None,
        "nearest_image_arcsec": round(near, 2),
        "old_rule": bool(near <= OLD_RULE_ARCSEC),
        "priority": "A" if inm else ("B" if near <= PRIORITY_B_ARCSEC else "C"),
    }


def review_list(maps, hits, morph, field, known=None, offset=(0.0, 0.0)) -> dict:
    """Every unflagged morph ``pass`` hit, with position labels; position is not a stage."""
    from scipy import ndimage

    mask, w, band, scale, snr = load_mask(maps, field, with_snr=True)
    edge = ndimage.distance_transform_edt(~mask) * scale
    ira, idec = image_positions(field)
    h = Table.read(hits)
    verdict = {t["id"]: t["verdict"] for t in json.loads(Path(morph).read_text())["targets"]}
    ks = []
    if known:
        from arc_diff import load_known

        ks = load_known(known, tuple(offset))
    rows, n_unflagged = [], 0
    for r in h:
        if _truthy(r["near_bright_core"]):
            continue
        n_unflagged += 1
        hid = str(r["hit_id"])
        if hid not in verdict:
            raise ValueError(f"{hid} has no morphology verdict; run arc_diff_morph.py on all hits")
        if verdict[hid] != "pass":
            continue
        ra, dec = float(r["ra"]), float(r["dec"])
        names = [
            k["name"] for k in ks if _sep(float(k["ra"]), float(k["dec"]), ra, dec) <= MATCH_ARCSEC
        ]
        rows.append(
            {
                "id": hid,
                "ra": ra,
                "dec": dec,
                "bands": str(r["bands"]) if "bands" in h.colnames else None,
                "max_abs_sigma": float(r["max_abs_sigma"])
                if "max_abs_sigma" in h.colnames
                else None,
                "known": names[0] if names else None,
                **_labels(mask, w, scale, snr, edge, ira, idec, ra, dec),
            }
        )
    rows.sort(key=lambda x: (x["priority"], -(x["max_abs_sigma"] or 0)))
    return {
        "mask_band": band,
        "mask_area_frac": float(mask.mean()),
        "n_hits": len(h),
        "n_unflagged": n_unflagged,
        "n_review": len(rows),
        "n_added": sum(not x["old_rule"] for x in rows),
        "by_priority": {k: sum(x["priority"] == k for x in rows) for k in "ABC"},
        "review": rows,
    }


def inject_review(runs, field) -> dict:
    """Whole-box efficiency of the review-list chain (hit, unflagged, morph pass), by stratum."""
    per: dict[str, dict] = {}
    masks = []
    for maps, injected, morph in runs:
        mask, w, band, scale, snr = load_mask(maps, field, with_snr=True)
        masks.append({"maps": Path(maps).name, "band": band, "area_frac": float(mask.mean())})
        verdict = {t["id"]: t["verdict"] for t in json.loads(Path(morph).read_text())["targets"]}
        with open(injected, newline="") as f:
            for r in csv.DictReader(f):
                p = _pixel(mask.shape, w, float(r["ra"]), float(r["dec"]))
                if p is None:
                    raise ValueError(
                        f"{r['name']} lies off the maps; every injected star must count"
                    )
                if r["hit"] and r["hit"] not in verdict:
                    raise ValueError(f"{r['hit']} has no morphology verdict; run --all-hits")
                inm = bool(mask[p[1], p[0]])
                level = "level_ge3" if snr[p[1], p[0]] >= THRESH_SIGMA else "level_lt3"
                ok_hit = bool(r["hit"])
                ok = (
                    ok_hit
                    and not _truthy(r["near_bright_core"])
                    and verdict.get(r["hit"]) == "pass"
                )
                for stratum in ("all", "in_mask" if inm else "out_mask", level):
                    c = per.setdefault(r["mag"], {}).setdefault(
                        stratum, {"n": 0, "hit": 0, "review": 0, "hit_outcomes": {}}
                    )
                    c["n"] += 1
                    c["hit"] += ok_hit
                    c["review"] += ok
                    if ok_hit:  # where a recovered star is lost: bright-core flag or morph verdict
                        o = (
                            "near_bright_core"
                            if _truthy(r["near_bright_core"])
                            else verdict.get(r["hit"])
                        )
                        c["hit_outcomes"][str(o)] = c["hit_outcomes"].get(str(o), 0) + 1
    for strata in per.values():
        for c in strata.values():
            c["eff_hit"] = round(c["hit"] / c["n"], 3) if c["n"] else None
            c["eff_review"] = round(c["review"] / c["n"], 3) if c["n"] else None
    return {"masks": masks, "n_runs": len(runs), "per_mag": per}


def mask_figure(maps, field, known, offset, out) -> None:
    """Mask-band epoch-1 image, mask contour, catalogued images (+) and known signals (o)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from arc_diff import load_known

    mask, w, band, _ = load_mask(maps, field)
    with np.load(maps) as m:
        img = m[f"{band}_e1img"]
    x, y = w.all_world2pix(*image_positions(field), 0)
    ks = load_known(known, tuple(offset))
    kx, ky = w.all_world2pix([float(k["ra"]) for k in ks], [float(k["dec"]) for k in ks], 0)
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    lo, hi = np.nanpercentile(img, [5, 99.5])
    for a in ax:
        a.imshow(img, origin="lower", cmap="gray_r", vmin=lo, vmax=hi)
        a.plot(x, y, "r+", ms=10)
        a.plot(kx, ky, "o", mfc="none", mec="m", ms=14)
        a.set_xlim(0, img.shape[1])
        a.set_ylim(0, img.shape[0])
    ax[1].contour(mask, [0.5], colors="c")
    ax[0].set_title(f"{field} {band} epoch 1; + catalogued images, o known")
    ax[1].set_title(f"arc mask ({mask.mean():.3f} of box)")
    fig.savefig(out, dpi=90, bbox_inches="tight")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("known", "inject", "review", "inject-review"):
        p = sub.add_parser(name)
        nargs = "+" if name.startswith("inject") else None
        p.add_argument("--maps", type=Path, required=True, nargs=nargs)
        p.add_argument("--morph", type=Path, required=True, nargs=nargs)
        p.add_argument(
            "--field", choices=("smacs", "macs0416", "sunrise", "elgordo"), required=True
        )
        p.add_argument("--out", type=Path, required=True)
        if name in ("known", "review"):
            p.add_argument("--hits", type=Path, required=True)
            p.add_argument("--known", type=Path, required=name == "known")
            p.add_argument("--frame-offset", type=float, nargs=2, default=[0.0, 0.0])
            if name == "known":
                p.add_argument("--figure", type=Path, help="also draw the mask figure (PNG)")
        else:
            p.add_argument("--injected", type=Path, required=True, nargs="+")
    a = ap.parse_args(argv)
    if a.cmd == "known":
        res = known_chain(a.maps, a.hits, a.morph, a.known, a.field, a.frame_offset)
        if a.figure:
            mask_figure(a.maps, a.field, a.known, a.frame_offset, a.figure)
    elif a.cmd == "review":
        res = review_list(a.maps, a.hits, a.morph, a.field, a.known, a.frame_offset)
    else:
        if not len(a.maps) == len(a.injected) == len(a.morph):
            ap.error("--maps, --injected and --morph need one entry per run")
        runs = list(zip(a.maps, a.injected, a.morph, strict=True))
        res = inject_chain(runs, a.field) if a.cmd == "inject" else inject_review(runs, a.field)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    print(
        json.dumps({k: v for k, v in res.items() if k not in ("hits_in_mask", "review")}, indent=1)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
