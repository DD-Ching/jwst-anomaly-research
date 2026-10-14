"""Abell 2744: the UNCOVER v2.0 lens model (Furtak+2023, Price+2025) as an independent second model.

Step 1 of the UNCOVER cross-check (TASKS.md "Supporting" item 1): validate the published v2.0
best-fit deflection maps (0.1"/pix, D_LS/D_S = 1) before any counter-image test uses them.

* **Known cases (own list):** UNCOVER's own multiple-image table (``A2744_multiple-images_v2.0``,
  ids ``S.KI``: source ``S.K``, image ``I``), every image at its catalogued redshift, traced
  through ``lensmodel.MapLensModel`` and re-imaged from one source per family. The release quotes
  a lens-plane rms of 0.60"; ours must come within ``lens_consistency.IMAGE_LIST_RMS_FACTOR`` of
  it (ASSUMPTION, as for the HFF maps, D-035) before the model is used. The source is the
  magnification-weighted least-squares source (``wls_source``): the mean back-traced source
  fails the gate (rms 4.8") because it falls outside the caustic for fold pairs (34.1, 80.1);
  the solver grid step (0.25" vs 0.1") changes nothing.
* **Cross-model (Bergamini+23 list):** the Bergamini+23 ``arcs.dat`` images with their
  ``best.par`` redshifts through the UNCOVER maps: per family the image-plane residuals and the
  number of predicted images against the number catalogued. Families with predicted images that
  are not catalogued (``n_extra``) are where the two models disagree on multiplicity; those are
  the inputs of the D-041 cross-model rule in the next step.

Everything computed is ``model_prediction``; the image positions are ``observed``.

Data: the maps are fetched once from the UNCOVER Google Drive release (sha256-pinned; 216 MB
compressed in three files, ~250 MB on disk while the run lasts; the maps are deleted at the end
unless ``--keep-maps``).

* **4.2c (D-041 second-model test):** the Bergamini+23 ``absent`` flag 4.2c
  (docs/fields/abell2744.md). ``check_42c`` solves family 4.2 in the UNCOVER maps and measures
  1"-high-pass F150W/F277W/F444W peak fluxes (r = 0.2") at the UNCOVER 4.2c and at 4.2a/b, plus
  injections at 4.2c of the most magnified catalogued image's knot, scaled by the mu ratio.

Usage: python scripts/a2744_uncover.py [--out results/lens_consistency/abell2744-uncover2]
       [--check-42c]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

_DRIVE = "https://drive.usercontent.google.com/download?export=download&confirm=t&id="
# UNCOVER DR4 lensing maps v2.0 (2024-08-07), Best-model_low-resolution_100mas; accessed 2026-10-14
UNCOVER_V2 = {
    "alpha_x": (
        _DRIVE + "1BHZ63x9nYAd5O0s1K_mO1UndyaekeJjz",
        "7ba107377e1585327d9ad95f8d847941dbcb7864611c442af6b712b455df4803",
    ),
    "alpha_y": (
        _DRIVE + "1Ht5KmqeXifqI8Q0tSRcVzPvEz6iAaH3w",
        "23469d560df560f6c7136cf273695aebe0fedd09853f256599b55c30335d156a",
    ),
    "images": (
        _DRIVE + "1uQQQBA1W_zd-y8SILZfBYKpVGAGoTxMG",
        "0d2649b3991129df0b775965caadf10109714d32cd804dcc76afa17c9c789119",
    ),
}
Z_LENS = 0.308
COSMOLOGY = (70.0, 0.3)  # ASSUMPTION: the Zitrin-code default; maps are at D_LS/D_S = 1 anyway
QUOTED_RMS = 0.60  # README v2.0, lens-plane rms (arcsec)
MATCH_ARCSEC = 1.5  # a predicted image within this of a catalogued one is that image (ASSUMPTION)


def wls_source(beta: np.ndarray, jac: np.ndarray) -> np.ndarray:
    """Magnification-weighted least-squares source of one family.

    ``beta`` (n, 2) are the back-traced sources, ``jac`` (n, 2, 2) the lens Jacobians A at the
    images. Minimises sum |A_i^-1 (beta - beta_i)|^2, the image-plane offsets to first order
    (the source-plane chi^2 of Keeton 2001 / Oguri 2010 without position errors)."""
    minv = np.linalg.inv(jac)
    w = np.einsum("nki,nkj->nij", minv, minv)
    return np.linalg.solve(w.sum(0), np.einsum("nij,nj->i", w, beta))


def family_source(model, obs: Table) -> tuple[float, float]:
    """``wls_source`` of a family's back-traced rows (``backtrace_images`` output)."""
    z = float(obs["z_used"][0])
    ox, oy = model.to_frame(np.asarray(obs["ra"], float), np.asarray(obs["dec"], float))
    m = model.lens_map(ox, oy, z)
    jac = np.stack([np.stack([m["a11"], m["a12"]], -1), np.stack([m["a21"], m["a22"]], -1)], -2)
    beta = np.stack([np.asarray(obs["beta_x"], float), np.asarray(obs["beta_y"], float)], -1)
    bx, by = wls_source(beta, jac)
    return float(bx), float(by)


def imageplane(model, grid, bt: Table) -> Table:
    """Per catalogued image: offset to the nearest predicted image of its family's
    ``family_source`` (``lensmodel.imageplane_residuals`` with the WLS source)."""
    systems = np.asarray(bt["system"]).astype(str)
    dtheta = np.full(len(bt), np.nan)
    mu = np.full(len(bt), np.nan)
    n_pred = np.zeros(len(bt), int)
    shared = np.zeros(len(bt), bool)
    for fam in dict.fromkeys(systems):
        sel = np.flatnonzero((systems == fam) & np.isfinite(bt["beta_x"]))
        if len(sel) < 2:
            continue
        obs = bt[sel]
        pred = lensmodel.find_images(
            model, grid, *family_source(model, obs), float(obs["z_used"][0])
        )
        n_pred[sel] = len(pred)
        if not len(pred):
            continue
        ox, oy = model.to_frame(np.asarray(obs["ra"], float), np.asarray(obs["dec"], float))
        d = np.hypot(ox[:, None] - pred["x"][None, :], oy[:, None] - pred["y"][None, :])
        k = np.argmin(d, axis=1)
        dtheta[sel] = d[np.arange(len(sel)), k]
        mu[sel] = np.asarray(pred["magnification"])[k]
        shared[sel] = np.bincount(k, minlength=len(pred))[k] > 1
    return Table(
        {
            "image_id": bt["image_id"],
            "system": systems,
            "z_used": bt["z_used"],
            "dtheta_arcsec": dtheta,
            "magnification": mu,
            "n_predicted": n_pred,
            "shared_match": shared,
        }
    )


def uncover_system(image_id: str) -> str:
    """Source of an UNCOVER image id: "S.KI" (two decimals) -> "S.K"; "S.I" (one decimal,
    single-knot systems such as 83.1-83.3) -> "S". Ids are stored as floats, so "S.K0" cannot
    occur (image numbers start at 1)."""
    whole, frac = image_id.split(".")
    return f"{whole}.{frac[0]}" if len(frac) == 2 else whole


def uncover_images(path: Path) -> Table:
    """UNCOVER's multiple-image table as a ``load_lenstool_images``-like table."""
    t = Table.read(path)
    ids = [repr(float(i)) for i in t["ID"]]
    out = Table(
        {
            "image_id": ids,
            "system": [uncover_system(i) for i in ids],
            "ra": np.asarray(t["RA"], float),
            "dec": np.asarray(t["Dec"], float),
            "a": np.full(len(t), np.nan),
            "z": np.asarray(t["z"], float),
        }
    )
    out.meta.update(provenance="observed", source=str(path))
    return out


def multiplicity(model, grid, bt: Table) -> Table:
    """Per family: catalogued and predicted image counts, and predicted images more than
    ``MATCH_ARCSEC`` from every catalogued image (``n_extra``; their |mu| listed)."""
    rows = []
    for fam in np.unique(np.asarray(bt["system"]).astype(str)):
        sel = bt[(np.asarray(bt["system"]).astype(str) == fam) & np.isfinite(bt["beta_x"])]
        if len(sel) < 2:
            continue
        z_s = float(sel["z_used"][0])
        pred = lensmodel.find_images(model, grid, *family_source(model, sel), z_s)
        ox, oy = model.to_frame(np.asarray(sel["ra"], float), np.asarray(sel["dec"], float))
        if len(pred):
            d = np.hypot(pred["x"][:, None] - ox[None, :], pred["y"][:, None] - oy[None, :])
            extra = np.flatnonzero(d.min(axis=1) > MATCH_ARCSEC)
        else:
            extra = np.array([], int)
        rows.append(
            (
                fam,
                z_s,
                len(sel),
                len(pred),
                len(extra),
                ",".join(f"{abs(float(pred['magnification'][k])):.1f}" for k in extra),
                ",".join(f"{float(pred['ra'][k]):.5f}:{float(pred['dec'][k]):.5f}" for k in extra),
            )
        )
    return Table(
        rows=rows or None,
        names=(
            "system",
            "z_used",
            "n_catalogued",
            "n_predicted",
            "n_extra",
            "extra_mu",
            "extra_radec",
        ),
        dtype=(str, float, int, int, int, str, str),
    )


S3_I2D = "s3://stpubdata/jwst/public/jw02561/L3/t/o001/jw02561-o001_t003_nircam_clear-{}_i2d.fits"
BERG_42C = (3.57972, -30.40835)  # Bergamini+23 prediction of 4.2c (docs/fields/abell2744.md)
HIGHPASS_ARCSEC = 1.0  # median-filter width removing BCG/ICL light (ASSUMPTION, as halo_phot)
PEAK_RADIUS = {"ref": 0.4, "target": 1.0}  # peak search radius, arcsec (ASSUMPTION)
INJECT_SCALES = (1.0, 0.5, 0.25)


def _highpass(sci, xx, yy):
    return lc._median_subtract(sci, xx, yy, HIGHPASS_ARCSEC)


def _peak_flux(hp, err, xx, yy, radius: float, r_ap: float = 0.2) -> dict:
    inner = np.where((np.hypot(xx, yy) <= radius) & np.isfinite(hp), hp, -np.inf)
    k = int(np.argmax(inner))
    cx, cy = float(xx.flat[k]), float(yy.flat[k])
    ap = np.hypot(xx - cx, yy - cy) <= r_ap
    if not np.isfinite(inner).any():
        return {"dx_west": np.nan, "dy_north": np.nan, "flux": np.nan, "snr": np.nan}
    f, e = float(np.nansum(hp[ap])), float(np.sqrt(np.nansum(err[ap] ** 2)))
    return {"dx_west": cx, "dy_north": cy, "flux": f, "snr": f / e if e > 0 else np.nan}


def check_42c(model, grid, bt: Table) -> dict:
    """UNCOVER's 4.2 solution, high-pass peak photometry at its 4.2c and at 4.2a/b, and
    injections of the most magnified catalogued image's knot (scaled by the mu ratio) at 4.2c."""
    sel = bt[(np.asarray(bt["system"]).astype(str) == "4.2") & np.isfinite(bt["beta_x"])]
    if len(sel) < 2:
        return {"error": "family 4.2 has fewer than two traced images"}
    z = float(sel["z_used"][0])
    # the grid covers the whole cluster (4.2c lies ~44" from 4.2a/b)
    pred = lensmodel.find_images(model, grid, *family_source(model, sel), z)
    if len(pred) <= len(sel):
        return {"error": f"no extra image predicted ({len(pred)} images)"}
    ox, oy = model.to_frame(np.asarray(sel["ra"], float), np.asarray(sel["dec"], float))
    d = np.hypot(pred["x"][:, None] - ox[None, :], pred["y"][:, None] - oy[None, :])
    k_new = int(np.argmax(d.min(axis=1)))  # the predicted image farthest from a and b
    match = np.argmin(d, axis=0)  # predicted image of each catalogued one
    c = pred[k_new]
    bx, by = model.to_frame(*BERG_42C)
    out = {
        "z": z,
        "predicted": [
            {"ra": float(r["ra"]), "dec": float(r["dec"]), "mu": float(r["magnification"])}
            for r in pred
        ],
        "c_ra_dec": [float(c["ra"]), float(c["dec"])],
        "c_mu": float(c["magnification"]),
        "c_offset_from_bergamini_arcsec": float(np.hypot(c["x"] - bx, c["y"] - by)),
        "bands": {},
    }
    pos = {
        str(sel["image_id"][i]): (float(sel["ra"][i]), float(sel["dec"][i]))
        for i in range(len(sel))
    }
    mus = {str(sel["image_id"][i]): float(pred["magnification"][match[i]]) for i in range(len(sel))}
    out["mu_catalogued"] = mus
    ref_id = max(mus, key=lambda i: abs(mus[i]))
    for band in ("f150w", "f277w", "f444w"):
        stamp, close = lc.image_stamper(S3_I2D.format(band), 4.0)
        try:
            res = {}
            for iid, (ra, dec) in pos.items():
                sci, err, xx, yy = stamp(ra, dec)
                res[iid] = _peak_flux(_highpass(sci, xx, yy), err, xx, yy, PEAK_RADIUS["ref"])
            sci, err, xx, yy = stamp(float(c["ra"]), float(c["dec"]))
            res["4.2c_uncover"] = _peak_flux(
                _highpass(sci, xx, yy), err, xx, yy, PEAK_RADIUS["target"]
            )
            # injection: the reference image's positive high-pass knot (1.2" box) scaled to mu(c)
            rsci, _, rx, ry = stamp(*pos[ref_id])
            rhp = _highpass(rsci, rx, ry)
            j = np.unravel_index(
                np.argmax(np.where(np.hypot(rx, ry) < PEAK_RADIUS["ref"], rhp, -np.inf)), rhp.shape
            )
            n = int(round(0.6 / abs(float(rx[0, 1] - rx[0, 0]))))
            knot = np.clip(rhp[j[0] - n : j[0] + n + 1, j[1] - n : j[1] + n + 1], 0, None)
            ci, cj = np.array(sci.shape) // 2
            inj = {}
            for s in INJECT_SCALES:
                im = sci.copy()
                im[ci - n : ci + n + 1, cj - n : cj + n + 1] += (
                    s * knot * abs(out["c_mu"] / mus[ref_id])
                )
                inj[str(s)] = _peak_flux(_highpass(im, xx, yy), err, xx, yy, PEAK_RADIUS["target"])
            res["injections"] = inj
            res["expected_c_flux"] = res[ref_id]["flux"] * abs(out["c_mu"] / mus[ref_id])
            out["bands"][band] = res
        finally:
            close()
    out["reference"] = ref_id
    return out


def run(model, images: Table, zml: dict, label: str, out: Path) -> tuple[dict, object, Table]:
    bt = lensmodel.backtrace_images(model, images, zml)
    traced = images[np.isfinite(bt["beta_x"])]
    half = lc.grid_half_width(model, traced)
    grid = lensmodel.DeflectionGrid.compute(model, half, 0.25)
    ip = imageplane(model, grid, bt)
    d = np.asarray(ip["dtheta_arcsec"], float)
    ok = np.isfinite(d)
    mult = multiplicity(model, grid, bt)
    fam = np.asarray(bt["system"]).astype(str)
    traced = np.isfinite(bt["beta_x"])
    multi = np.array([np.sum(traced & (fam == f)) >= 2 for f in fam])
    for t in (ip, mult):
        t.meta.update(provenance="model_prediction", source=f"abell2744-uncover2 vs {label}")
    ip.write(out / f"imageplane_{label}.ecsv", overwrite=True)
    mult.write(out / f"multiplicity_{label}.ecsv", overwrite=True)
    return (
        {
            "image_list": label,
            "n_images": len(images),
            "n_traced": int(np.isfinite(bt["beta_x"]).sum()),
            "n_reimaged": int(ok.sum()),
            # traced images of families with >= 2 traced images that got no prediction
            "n_unsolved": int(((ip["n_predicted"] == 0) & np.isfinite(bt["beta_x"]) & multi).sum()),
            "n_shared_match": int(np.asarray(ip["shared_match"]).sum()),
            "rms_dtheta_arcsec": float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else None,
            "median_dtheta_arcsec": float(np.median(d[ok])) if ok.any() else None,
            "p90_dtheta_arcsec": float(np.percentile(d[ok], 90)) if ok.any() else None,
            "max_dtheta_arcsec": float(d[ok].max()) if ok.any() else None,
            "n_over_2arcsec": int((d[ok] > 2.0).sum()),
            "families": len(mult),
            "families_with_extra": int((mult["n_extra"] > 0).sum()) if len(mult) else 0,
            "families_short": int((mult["n_predicted"] < mult["n_catalogued"]).sum())
            if len(mult)
            else 0,
        },
        grid,
        bt,
    )


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path("results/lens_consistency/abell2744-uncover2"))
    ap.add_argument("--check-42c", action="store_true", help="4.2c photometry (S3 byte ranges)")
    ap.add_argument("--keep-maps", action="store_true", help="keep the cached deflection maps")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    files = {k: fetch_catalog(url, sha) for k, (url, sha) in UNCOVER_V2.items()}
    h0, om = COSMOLOGY
    model = lensmodel.MapLensModel.from_fits(
        files["alpha_x"],
        files["alpha_y"],
        Z_LENS,
        FlatLambdaCDM(H0=h0, Om0=om),
        source="abell2744-uncover2",
    )
    summary = {
        "model": "abell2744-uncover2",
        "model_sha256": model.sha256,
        "quoted_rms_arcsec": QUOTED_RMS,
    }
    own, _, _ = run(model, uncover_images(files["images"]), {}, "uncover2", args.out)
    # every traced multi-image family must re-image; the rms of the rest is not enough
    own["gate_passes"] = bool(
        own["rms_dtheta_arcsec"] is not None
        and own["n_unsolved"] == 0
        and own["rms_dtheta_arcsec"] <= lc.IMAGE_LIST_RMS_FACTOR * QUOTED_RMS
    )
    summary["own_list"] = own
    _, bfiles, par = lc.load_model("abell2744-bergamini23")
    berg = lensmodel.load_lenstool_images(bfiles["arcs.dat"])
    summary["bergamini23_list"], grid, bt = run(
        model, berg, par["z_m_limit"], "bergamini23", args.out
    )
    if args.check_42c:
        summary["check_42c"] = check_42c(model, grid, bt)
    if not args.keep_maps:  # cloud-disk rule: ~250 MB of maps, deleted after use
        for key in ("alpha_x", "alpha_y"):
            files[key].unlink(missing_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return summary


if __name__ == "__main__":
    main()
