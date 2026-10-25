"""MACS0416: forced photometry at the MARS-only extra images (follow-up of D-034 addendum ar)

Frozen plan: ``results/lens_consistency/m0416_mars_only_phot/PLAN.md``. Solves the 40 families of
``m0416_mars_extra`` with the MARS maps, pairs MARS images with observed images (``match_pairs``),
takes the MARS images paired with no observed image and far from every CANUCS extra prediction
(``mars_only``) as targets, and runs the ``m0417_arcref`` photometry with the CANUCS DR1 catalogue
at MARS positions and MARS |mu| (controls: the families' units), then the local-noise rule.

    python scripts/m0416_mars_only_phot.py --zip DIR --cat <photometry-cat.fits.gz> \
        --image <f200w bgsub-sci.fits.gz> [--png OUT.png]
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import lens_consistency as lc
import m0416_canucsref as mc
import m0416_mars_extra as me
import m0417_arcref as ar
import m0417_canucsref as cr
import m1423_canucsref as m14
import mars_hff_heldout as mo
import mars_hff_inputs as mh
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel
from jwst_anomaly.provenance import file_sha256

OUT = Path("results/lens_consistency/m0416_mars_only_phot")
MATCH_ARCSEC = me.MATCH_ARCSEC  # 3.0", PLAN.md rules 1-2 (ASSUMPTION)
SEARCH_ARCSEC = 1.5  # ASSUMPTION (PLAN.md rule 3: one MARS map pixel)


def match_pairs(px, py, ox, oy, r: float = MATCH_ARCSEC) -> np.ndarray:
    """PLAN.md rule 1: for every observed image, the index of its paired predicted image (-1 if
    none); one-to-one, greedy by distance, within ``r`` (same order as ``me.match_observed``)."""
    px, py, ox, oy = (np.asarray(v, float) for v in (px, py, ox, oy))
    pair = np.full(len(ox), -1, int)
    if not len(px) or not len(ox):
        return pair
    d = np.hypot(ox[:, None] - px[None, :], oy[:, None] - py[None, :])
    used = np.zeros(len(px), bool)
    for flat in np.argsort(d, axis=None):
        i, j = np.unravel_index(flat, d.shape)
        if d[i, j] > r:
            break
        if pair[i] < 0 and not used[j]:
            pair[i], used[j] = j, True
    return pair


def mars_only(px, py, pair, ex, ey, r: float = MATCH_ARCSEC) -> np.ndarray:
    """PLAN.md rule 2: mask of predicted images paired with no observed image and farther than
    ``r`` from every extra prediction (``ex``, ``ey``) of the family."""
    px, py = np.asarray(px, float), np.asarray(py, float)
    free = np.ones(len(px), bool)
    free[np.asarray(pair, int)[np.asarray(pair, int) >= 0]] = False
    ex, ey = np.asarray(ex, float), np.asarray(ey, float)
    if len(ex):
        d = np.hypot(px[:, None] - ex[None, :], py[:, None] - ey[None, :]).min(axis=1)
        free &= d > r
    return free


def reading(class_local: str, c2_passes: bool) -> str:
    """PLAN.md rule 6."""
    if class_local == "absent" and c2_passes:
        return "lensing_residual_input"
    if class_local == "recovered":
        return "seen_where_only_mars_predicts"
    return "no_reading"


def unshift(model, x, y, off: dict) -> tuple[np.ndarray, np.ndarray]:
    """MARS-frame positions (arcsec) to CANUCS-frame RA / Dec (PLAN.md rule 2)."""
    ra, dec = model.to_sky(x, y)
    dec = np.asarray(dec, float)
    ra = np.asarray(ra, float) - off["dx_arcsec"] / 3600.0 / np.cos(np.deg2rad(dec))
    return ra, dec - off["dy_arcsec"] / 3600.0


def solve(model, grid, t: Table, pred: Table, off: dict, systems=None) -> tuple[dict, list[dict]]:
    """Per family (``systems``; default: the families with an extra prediction): observed images
    (CANUCS frame) with their paired MARS position / |mu|, and the MARS-only targets (CANUCS
    frame)."""
    bt = lensmodel.backtrace_images(model, t, {})
    sysv = np.asarray(t["system"], str)
    ox, oy = model.to_frame(np.asarray(t["ra"], float), np.asarray(t["dec"], float))
    ex, ey = model.to_frame(np.asarray(pred["ra"], float), np.asarray(pred["dec"], float))
    psys = np.asarray(pred["system"], str)
    fams, targets = {}, []
    for s in dict.fromkeys(psys) if systems is None else systems:
        fam = np.flatnonzero(sysv == s)
        bx, by = np.asarray(bt["beta_x"])[fam], np.asarray(bt["beta_y"])[fam]
        if not len(fam) or not (np.isfinite(bx).all() and np.isfinite(by).all()):
            continue
        z = float(np.asarray(t["z"])[fam][0])
        p = lensmodel.find_images(model, grid, float(bx.mean()), float(by.mean()), z)
        px, py = np.asarray(p["x"], float), np.asarray(p["y"], float)
        mu = np.abs(np.asarray(p["magnification"], float))
        pair = match_pairs(px, py, ox[fam], oy[fam])
        pra, pdec = unshift(model, px, py, off)
        ims = []
        for i, j in zip(fam, pair, strict=True):
            if j < 0:
                continue
            ra0, dec0 = unshift(model, ox[i], oy[i], off)  # the image's own CANUCS position
            ims.append(
                {
                    "id": str(t["image_id"][i]),
                    "ra": float(ra0),
                    "dec": float(dec0),
                    "mars_ra": float(pra[j]),
                    "mars_dec": float(pdec[j]),
                    "mu": float(mu[j]),
                }
            )
        fams[s] = ims
        k = psys == s
        for j in np.flatnonzero(mars_only(px, py, pair, ex[k], ey[k])):
            targets.append(
                {"system": s, "ra": float(pra[j]), "dec": float(pdec[j]), "mu": float(mu[j])}
            )
    return fams, targets


def photometry(fams: dict, targets: list[dict], cat: dict, depth: float) -> tuple[list, list]:
    """PLAN.md rules 3-4 (``m0417_arcref`` reference / classify, MARS positions and |mu|)."""

    def test(ra, dec, mu_t, cands) -> dict:
        lref, state, rcp = ar.reference(cands, cat)
        row = {"L_ref": lref, "ref_counterpart": -1, "arc_ref": False, "f_pred": float("nan")}
        if state:
            return {**row, "class": state, "n_sources": 0, "max_ratio": float("nan")}
        f_pred = lref * mu_t
        idx = ar.compared_sources(ra, dec, SEARCH_ARCSEC, cat, rcp, contain=True)
        ratios = [float(cat["flux"][k] / f_pred) for k in idx]
        return {
            **row,
            "ref_counterpart": int(cat["label"][rcp]),
            "arc_ref": ar.is_arc(rcp, cat),
            "f_pred": f_pred,
            "class": ar.classify(f_pred, depth, ratios),
            "n_sources": len(ratios),
            "max_ratio": max(ratios) if ratios else float("nan"),
        }

    units = {}
    for s, ims in fams.items():
        for im in ims:
            im["cp"] = ar.counterpart(im["ra"], im["dec"], cat)
        units[s] = ar.units(ims)
    controls = []
    for s, us in units.items():
        for u in us:
            first = next(im for im in fams[s] if im["id"] == u["ids"][0])
            res = test(
                first["mars_ra"], first["mars_dec"], sum(u["mus"]), [v for v in us if v is not u]
            )
            controls.append(
                {
                    "system": s,
                    "images": "+".join(u["ids"]),
                    "ra": first["mars_ra"],
                    "dec": first["mars_dec"],
                    "sum_abs_mu": sum(u["mus"]),
                    **res,
                }
            )
    rows = []
    for tg in targets:
        res = test(tg["ra"], tg["dec"], tg["mu"], units.get(tg["system"], []))
        rows.append(
            {"system": tg["system"], "ra": tg["ra"], "dec": tg["dec"], "abs_mu": tg["mu"], **res}
        )
    return controls, rows


def cutout_png(d, w, rows: list[dict], fams: dict, path: Path, half: float = 4.0) -> None:
    """Contact sheet: each target and its family's observed images, 8" boxes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = []
    for r in rows:
        panels.append((f"{r['system']} MARS-only |mu| {r['abs_mu']:.1f}", r["ra"], r["dec"]))
        for im in fams.get(r["system"], []):
            panels.append((f"{im['id']} observed", im["ra"], im["dec"]))
    if not panels:
        return
    fig, axs = plt.subplots(1, len(panels), figsize=(3 * len(panels), 3.2), squeeze=False)
    h = half / m14.PIX
    for ax, (lab, ra, dec) in zip(axs[0], panels, strict=True):
        x, y = (float(v) for v in w.world_to_pixel_values(ra, dec))
        x0, y0 = int(round(x)), int(round(y))
        sub = d[max(0, int(y0 - h)) : int(y0 + h), max(0, int(x0 - h)) : int(x0 + h)]
        v = np.nanpercentile(sub, [5, 99.5]) if np.isfinite(sub).any() else (0, 1)
        ax.imshow(sub, origin="lower", cmap="gray_r", vmin=v[0], vmax=v[1])
        c = plt.Circle((x - max(0, int(x0 - h)), y - max(0, int(y0 - h))), SEARCH_ARCSEC / m14.PIX)
        c.set(fill=False, color="r", lw=0.8)
        ax.add_patch(c)
        ax.set_title(lab, fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def load(zdir: Path):
    """The MARS model and the shifted CANUCS image / extra-prediction tables (addendum as)."""
    meta = mh.fetch_zip(requests.Session(), zdir)
    folder, tag, z_lens = mh.FIELDS[me.FIELD]
    base = f"{mh.ROOT}{folder}/deflection_angle_{tag}_"
    with zipfile.ZipFile(zdir / mh.ZIP) as zf:
        mcat = mh.parse_catalogue(zf.read(mh.CATALOGUE).decode())
        (ax, hx), (ay, _) = (mh._fits(zf, base + f"{c}_w_header.fits") for c in "xy")
    model = mh.build_model(
        ax, ay, hx, False, z_lens, "see mars_hff_inputs", f"mars-{tag} (A_arcsec)"
    )
    mars = mcat[np.asarray(mcat["field"]) == me.FIELD]
    canucs = Table.read(mo.FIELDS[me.FIELD][0])
    off = mo.frame_offset(canucs, mars)
    t = mo.split_seen(mo.shift(canucs, off), mars)
    pred = mo.shift(Table.read(me.PREDICTED), off)
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, t), mh.GRID_STEP)
    return meta, model, t, pred, off, grid


def run(zdir: Path, cat_path: Path, image: Path, out: Path = OUT) -> dict:
    if file_sha256(cat_path) != mc.CAT_SHA or file_sha256(image) != mc.IMAGE_SHA:
        raise ValueError("--cat / --image are not the pinned CANUCS DR1 MACS0416 files")
    meta, model, t, pred, off, grid = load(zdir)
    fams, targets = solve(model, grid, t, pred, off)

    tab = Table.read(cat_path)
    cat = cr.to_cat(tab)
    depth = cr.depth_mag(cat["flux"], cat["snr"])
    controls, rows = photometry(fams, targets, cat, depth)
    c2 = ar.gate_c2([c["class"] for c in controls], [c["arc_ref"] for c in controls])
    with fits.open(image) as hl:
        w, d = WCS(hl[0].header), np.asarray(hl[0].data, float)
    src_x, src_y = (
        np.asarray(v, float) for v in w.world_to_pixel(SkyCoord(tab["RA"], tab["DEC"], unit="deg"))
    )
    rng = np.random.default_rng(m14.SEED)
    ctab = Table(rows=controls)
    m14.apply_local(ctab, d, w, src_x, src_y, rng)
    ttab = Table(rows=rows) if rows else Table()
    if rows:
        m14.apply_local(ttab, d, w, src_x, src_y, rng)
        ttab["reading"] = [reading(str(c), c2["passes"]) for c in ttab["class_local"]]
    out.mkdir(parents=True, exist_ok=True)
    for tb, name in ((ctab, "controls.ecsv"), (ttab, "targets.ecsv")):
        tb.meta["provenance"] = "derived"
        tb.meta["source"] = f"{out / 'PLAN.md'}; MARS maps (model_prediction) + CANUCS DR1 F200W"
        tb.write(out / name, format="ascii.ecsv", overwrite=True)
    cutout_png(d, w, [dict(r) for r in rows], fams, out / "targets.png")
    summary = {
        "plan": str(out / "PLAN.md"),
        "zip": meta,
        "canucs_cat_sha256": mc.CAT_SHA,
        "image_sha256": mc.IMAGE_SHA,
        "frame_offset": off,
        "depth_mag": depth,
        "search_arcsec": SEARCH_ARCSEC,
        "n_families": len(fams),
        "gate_C2": c2,
        "control_classes": m14._count(ctab["class"]),
        "control_classes_local": m14._count(ctab["class_local"]),
        "targets": [
            {
                k: (r[k].item() if hasattr(r[k], "item") else r[k])
                for k in (
                    "system",
                    "ra",
                    "dec",
                    "abs_mu",
                    "class",
                    "f_pred",
                    "n_sources",
                    "max_ratio",
                    "sigma_loc_njy",
                    "pred_over_sigma",
                    "class_local",
                    "reading",
                )
            }
            for r in ttab
        ]
        if rows
        else [],
    }
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zip", type=Path, required=True, help="directory for the MARS zip")
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--image", type=Path, required=True)
    ap.add_argument("--posthoc", action="store_true", help="only the post hoc table (after a run)")
    a = ap.parse_args(argv)
    if a.posthoc:
        print(json.dumps(posthoc(a.cat), indent=1))
        return 0
    try:
        s = run(a.zip, a.cat, a.image)
    finally:
        (a.zip / mh.ZIP).unlink(missing_ok=True)
    print(json.dumps(_clean({k: v for k, v in s.items() if k != "zip"}), indent=1))
    return 0


def posthoc(cat_path: Path, out: Path = OUT) -> dict:
    """Post hoc (not a frozen rule): distance from each target to its source's observed images
    (D115a / D115b share one source) and the compared sources with their CANUCS counterparts."""
    cat = cr.to_cat(Table.read(cat_path))
    images = Table.read(mo.FIELDS[me.FIELD][0])
    ctl = Table.read(out / "controls.ecsv")
    paired = {i for s in ctl["images"] for i in str(s).split("+")}
    res = []
    for tg in Table.read(out / "targets.ecsv"):
        ra, dec = float(tg["ra"]), float(tg["dec"])
        stem = str(tg["system"]).rstrip("abcdefgh")
        near = []
        for r in images:
            if lensmodel.image_family(str(r["image_id"])).rstrip("abcdefgh") != stem:
                continue
            k = ar.counterpart(float(r["ra"]), float(r["dec"]), cat)
            near.append(
                {
                    "image_id": str(r["image_id"]),
                    "d_arcsec": float(
                        np.hypot((r["ra"] - ra) * np.cos(np.radians(dec)), r["dec"] - dec) * 3600
                    ),
                    "counterpart": int(cat["label"][k]) if k >= 0 else -1,
                    "paired_to_mars": str(r["image_id"]) in paired,
                }
            )
        k_ref = np.flatnonzero(cat["label"] == tg["ref_counterpart"])
        rcp = int(k_ref[0]) if len(k_ref) else -1  # -1: no reference (no exclusion)
        srcs = [
            {
                "label": int(cat["label"][k]),
                "d_arcsec": float(
                    np.hypot((cat["ra"][k] - ra) * np.cos(np.radians(dec)), cat["dec"][k] - dec)
                    * 3600
                ),
                "f_over_fpred": float(cat["flux"][k] / tg["f_pred"]),
            }
            for k in ar.compared_sources(ra, dec, SEARCH_ARCSEC, cat, rcp, contain=True)
        ]
        res.append({"system": str(tg["system"]), "observed_images": near, "compared_sources": srcs})
    p = {"label": "post hoc, not a frozen rule", "targets": res}
    (out / "posthoc.json").write_text(json.dumps(_clean(p), indent=1, allow_nan=False) + "\n")
    return p


if __name__ == "__main__":
    raise SystemExit(main())
