"""WHL0137 (Sunrise): every image the four RELICS models predict for systems 1 and 2.

``sunrise_relics.py`` tests the catalogued images (Scofield+2025 Table 3) against the RELICS v1
models (Glafic, Lenstool, WSLAP+, Zitrin LTM) by local patches. Here each model's deflection
maps are read over a box around the cluster (HTTP byte ranges, held in memory, nothing stored),
wrapped as a ``lensmodel.MapLensModel`` and solved with ``lensmodel.find_images``:

* each knot (1.1, 1.2, 1.7, 2) is traced to the source plane at its redshift (system 1: z = 6.2,
  spectroscopic; system 2: the model's best z from ``results/sunrise/relics_positions.json``,
  a model_prediction) and the mean source position is re-imaged;
* predicted images within ``--self-arcsec`` of a catalogued image of the same system are that
  system's known images; the rest are counter-images (model_prediction);
* each counter-image with |mu| >= ``--min-mu`` gets forced photometry on the JWST F200W mosaic
  (``lens_consistency`` aperture, ERR scale and classes, D-029). Its predicted aperture flux is
  the median over the knot's catalogued images with |mu| < ``--max-ref-mu`` of their aperture
  flux times |mu_c / mu_ref| (point-source scaling; ASSUMPTION: it over-predicts extended
  images, so ``absent`` is conservative only for compact knots). With no such reference the
  class is ``no_reference``.

The search radius per counter-image is the larger of ``--search-arcsec`` and the spread of the
same counter-image across the system's redshift range (system 2: best z +- ``--dz``), capped at
``lens_consistency.MAX_SEARCH_ARCSEC``.

Rules (D-041, thresholds are ASSUMPTIONs):
* model copy: not applied as a class. With the rms of the re-predicted images it cannot fire
  (every counter-image is > ``--self-arcsec`` away), and with the missed images included it
  widens itself. Each counter-image records ``copy_distance_arcsec`` to the nearest catalogued
  image of its knot that the model does not re-predict;
* cross-model: ``n_models`` counts the models predicting a counter-image of the same knot within
  ``--cross-arcsec`` (demagnified images of other models do not count); an image
  predicted by only some models is model-dependent;
* injection: a compact image of the predicted flux is added at each tested position and the
  verdict repeated (``inject_class``, simulated).
A flag is an ``absent`` counter-image whose injection reads ``recovered`` and that every one of
at least two models predicts. Positions without mosaic coverage are ``off_image``.

Known cases: the catalogued images re-predicted within ``--self-arcsec``
(``n_catalogued_recovered``) and, as the positive control of the flux chain, each catalogued
image with |mu| < ``--max-ref-mu`` classified from the other references of its knot
(``control``; must come out ``recovered``).

Usage: python scripts/sunrise_counter.py [--out results/sunrise/counter_images.json]
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests
from astropy.cosmology import FlatLambdaCDM
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402
from sunrise_relics import BASE, DTYPES, KEYS, MODELS, _get, frame_check, read_header  # noqa: E402
from sunrise_scofield import IMAGES, system_of  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

Z_LENS = 0.566
COSMO = FlatLambdaCDM(H0=70.0, Om0=0.3)
CENTRE = (24.3555, -8.4570)  # model-frame origin near the cluster core (RA, Dec; ASSUMPTION)
F200W = "s3://stpubdata/jwst/public/jw02282/L3/t/o010/jw02282-o010_t001_nircam_clear-f200w_i2d.fits"
CHUNK = 8 * 2**20  # bytes per range request
MAP_THREADS, MODEL_THREADS = 4, 2  # <= 8 concurrent range requests to MAST (scripts/CLAUDE.md)
INJECT_FWHM_ARCSEC = lc.HALO_INJECT_FWHM_ARCSEC  # compact injected image (~F200W PSF; ASSUMPTION)
FLAG_CLASSES = ("absent",)
NOT_INDEPENDENT = ("demagnified",)  # not counted as another model's prediction of an image


def cross_model(results: list[dict], radius: float) -> None:
    """Add ``n_models``: how many models predict a counter-image of the same knot within
    ``radius`` of each counter-image (D-041: an image predicted by only some models is
    model-dependent and untestable). Other models' demagnified images do not count. Flags are
    ``absent`` images, recovered when injected, that every one of at least two models predicts."""
    pts: dict = {}
    for r in results:
        for q in r["counter_images"]:
            if q["class"] not in NOT_INDEPENDENT:
                pts.setdefault((r["model"], q["knot"]), []).append([q["ra"], q["dec"]])
    pts = {k: np.array(v) for k, v in pts.items()}
    for r in results:
        for c in r["counter_images"]:
            n = 0
            for other in results:
                p = pts.get((other["model"], c["knot"]))
                if p is None or not len(p):
                    continue
                dx = (p[:, 0] - c["ra"]) * np.cos(np.deg2rad(c["dec"])) * 3600.0
                dy = (p[:, 1] - c["dec"]) * 3600.0
                n += bool(np.hypot(dx, dy).min() <= radius)
            c["n_models"] = n
            c["flag"] = bool(
                c["class"] in FLAG_CLASSES
                and c.get("inject_class") == "recovered"
                and len(results) >= 2
                and n == len(results)
            )


def read_block(session, url: str, wcs_hdr, offset: int, rows: tuple[int, int]) -> np.ndarray:
    """Full-width rows [rows[0], rows[1]) of one map, read as parallel byte ranges (the
    columns are cut afterwards: a row segment per request would mean thousands of requests)."""
    nx, bpp = wcs_hdr["NAXIS1"], abs(wcs_hdr["BITPIX"]) // 8
    start, stop = offset + rows[0] * nx * bpp, offset + rows[1] * nx * bpp
    spans = [(s, min(s + CHUNK, stop)) for s in range(start, stop, CHUNK)]
    with ThreadPoolExecutor(MAP_THREADS) as ex:
        raw = b"".join(ex.map(lambda sp: _get(session, url, *sp), spans))
    return np.frombuffer(raw, DTYPES[wcs_hdr["BITPIX"]]).reshape(rows[1] - rows[0], nx)


def load_model(name: str, half: float) -> lensmodel.MapLensModel:
    """The model's deflection maps over a (2 half + margin) box around ``CENTRE``."""
    session = requests.Session()
    out, grid_keys = {}, ("NAXIS1", "NAXIS2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2")
    for k in ("defx", "defy"):
        url = BASE.format(m=name, k=KEYS[k])
        hdr, off = read_header(session, url)
        grid = tuple(hdr.get(c) for c in grid_keys)
        if out.setdefault("grid", grid) != grid:
            raise RuntimeError(f"{name}: defy grid differs from defx")
        if hdr.get("BSCALE", 1) != 1 or hdr.get("BZERO", 0) != 0:
            raise RuntimeError(f"{name}: scaled map pixels are not supported")
        wcs = WCS(hdr)
        pix = frame_check(wcs)  # north-up, east-left, square: checked before any data is read
        ic, jc = (float(v) for v in wcs.world_to_pixel_values(*CENTRE))
        m = int(np.ceil((half + 2.0) / pix))  # 2" margin for the finite differences
        i0, i1 = max(0, int(ic) - m), min(hdr["NAXIS1"], int(ic) + m + 1)
        j0, j1 = max(0, int(jc) - m), min(hdr["NAXIS2"], int(jc) + m + 1)
        out[k] = read_block(session, url, hdr, off, (j0, j1))[:, i0:i1]
        out["wcs"] = wcs[j0:j1, i0:i1]
    return lensmodel.MapLensModel(
        out["defx"],
        out["defy"],
        out["wcs"],
        Z_LENS,
        COSMO,
        source=BASE.format(m=name, k="{defx,defy}"),
        centre=CENTRE,
    )


def knots() -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for name in IMAGES:
        out.setdefault(system_of(name), []).append(name)
    return out


def predict(model, grid, knot: str, z: float) -> dict:
    """Catalogued images of one knot traced at z, mean source re-imaged (model frame)."""
    names = knots()[knot]
    ra = np.array([IMAGES[n][0] for n in names])
    dec = np.array([IMAGES[n][1] for n in names])
    ev = model.evaluate(ra, dec, z)
    bx = np.asarray(ev["x"] - ev["alpha_x"], float)
    by = np.asarray(ev["y"] - ev["alpha_y"], float)
    if not (np.isfinite(bx).all() and np.isfinite(by).all()):
        raise RuntimeError(f"knot {knot}: a catalogued image lies outside the loaded maps")
    imgs = lensmodel.find_images(model, grid, float(bx.mean()), float(by.mean()), z)
    return {
        "names": names,
        "ref_mu": np.asarray(ev["magnification"], float),
        "ref_xy": np.column_stack([ev["x"], ev["y"]]),
        "src_scatter": float(np.sqrt(np.mean((bx - bx.mean()) ** 2 + (by - by.mean()) ** 2))),
        "images": imgs,
    }


def system_xy(model, system: str) -> np.ndarray:
    names = [n for n in IMAGES if system_of(n).split(".")[0] == system]
    x, y = model.to_frame(
        np.array([IMAGES[n][0] for n in names]), np.array([IMAGES[n][1] for n in names])
    )
    return np.column_stack([x, y])


def split_images(pred: dict, sys_xy: np.ndarray, self_arcsec: float):
    """(known, counter, recovered mask, image-plane rms) of one knot's predicted images.

    Known images lie within ``self_arcsec`` of a catalogued image of the system (any knot: the
    knots are clumps of one galaxy, so another knot's image marks the same arc segment),
    counter-images do not. ``recovered`` marks the knot's catalogued images with a prediction
    within ``self_arcsec``; the rms is over those images' distances to their nearest prediction
    (the model's image-plane rms for this knot, derived; NaN when none is re-predicted)."""
    imgs = pred["images"]
    xy = np.column_stack([imgs["x"], imgs["y"]])
    if not len(xy):
        n = len(pred["ref_xy"])
        return imgs, imgs, np.zeros(n, bool), float("nan")
    d = np.hypot(*(xy[:, None, :] - sys_xy[None, :, :]).transpose(2, 0, 1))
    near = d.min(axis=1) <= self_arcsec
    dk = np.hypot(*(pred["ref_xy"][:, None, :] - xy[None, :, :]).transpose(2, 0, 1)).min(axis=1)
    ok = dk <= self_arcsec
    # rms over the re-predicted images only (diagnostic; the misses are listed separately)
    rms = float(np.sqrt(np.mean(dk[ok] ** 2))) if ok.any() else float("nan")
    return imgs[near], imgs[~near], ok, rms


def classify(stamp, ra: float, dec: float, pred: float, search: float) -> dict:
    """D-029 forced photometry of one predicted image, and the same verdict with a compact
    image of the predicted flux injected at the position (simulated; must read ``recovered``
    for an ``absent`` verdict to count)."""
    sci, err, xx, yy = stamp(ra, dec)
    _, e0 = lc.aperture_snr(sci, err, xx, yy, 0.0, 0.0)
    if not np.isfinite(e0) or e0 <= 0:
        return {"class": "off_image"}
    snr, fbest, dx, dy = lc.best_within(sci, err, xx, yy, search)
    out = {"search_arcsec": search, "best_snr": float(snr), "best_flux": float(fbest)}
    out["best_offset_arcsec"] = float(np.hypot(dx, dy))
    if not np.isfinite(pred) or pred <= 0:
        out["class"] = "no_reference"
        return out
    pred_snr, ratio = pred / e0, fbest / pred
    out |= {"pred_flux": pred, "pred_snr": pred_snr, "flux_ratio": ratio}
    out["class"] = lc.forced_class(pred_snr, float(snr), ratio)
    s = INJECT_FWHM_ARCSEC / 2.3548
    g = np.exp(-(xx**2 + yy**2) / (2 * s**2))
    g *= pred / g.sum()
    isnr, ifl, _, _ = lc.best_within(sci + g, err, xx, yy, search)
    out["inject_class"] = lc.forced_class(pred_snr, float(isnr), ifl / pred)
    return out


def match_spread(main_xy: np.ndarray, others: list[np.ndarray]) -> np.ndarray:
    """For each main counter-image, the largest distance to the nearest counter-image of each
    other redshift (inf when one has none)."""
    out = np.zeros(len(main_xy))
    for o in others:
        if not len(o):
            out[:] = np.inf
            continue
        d = np.hypot(*(main_xy[:, None, :] - o[None, :, :]).transpose(2, 0, 1)).min(axis=1)
        out = np.maximum(out, d)
    return out


def run_model(name: str, model, args, zbest: dict, stamp) -> dict:
    grid = lensmodel.DeflectionGrid.compute(model, args.half_width, args.step)
    # solver nodes outside the loaded maps cannot hold images (reported, not hidden)
    nan_frac = float(np.mean(~np.isfinite(grid.alpha_x)))
    sys_xy = {s: system_xy(model, s) for s in ("1", "2")}
    rows, checks, controls = [], [], []
    for knot in sorted(knots()):
        system = knot.split(".")[0]
        zs = [6.2] if system == "1" else [zbest[name]]
        if system == "2":
            zs += [zbest[name] - args.dz, zbest[name] + args.dz]
        preds = [predict(model, grid, knot, z) for z in zs]
        main = preds[0]
        _, counter, recovered, rms = split_images(main, sys_xy[system], args.self_arcsec)
        # match_spread pairs nearest counter-images across z (approximate; system 2 only)
        others = [split_images(p, sys_xy[system], args.self_arcsec)[1] for p in preds[1:]]
        unrecovered = main["ref_xy"][~recovered]
        cxy = np.column_stack([counter["x"], counter["y"]]) if len(counter) else np.zeros((0, 2))
        spread = match_spread(cxy, [np.column_stack([o["x"], o["y"]]) for o in others])
        checks.append(
            {
                "knot": knot,
                "z": zs[0],
                "n_catalogued": len(main["names"]),
                "n_catalogued_recovered": int(recovered.sum()),
                "unrecovered": [
                    n for n, ok in zip(main["names"], recovered, strict=True) if not ok
                ],
                "image_plane_rms_arcsec": rms,
                "n_predicted": len(main["images"]),
                "n_counter": len(counter),
                "source_scatter_arcsec": main["src_scatter"],
            }
        )
        # reference fluxes: catalogued images of this knot with moderate |mu|
        refs = []
        for n, mu in zip(main["names"], main["ref_mu"], strict=True):
            sci, err, xx, yy = stamp(*IMAGES[n])
            f, e = lc.aperture_snr(sci, err, xx, yy, 0.0, 0.0)
            refs.append({"image": n, "mu": float(mu), "flux": f, "err": e})
        good = [
            r
            for r in refs
            if abs(r["mu"]) < args.max_ref_mu and np.isfinite(r["flux"]) and r["flux"] > 0
        ]
        for r in good:  # positive control: leave-one-out, through the counter-image chain
            rest = [q for q in good if q is not r]
            if not rest:
                continue
            pred = float(np.median([q["flux"] * abs(r["mu"] / q["mu"]) for q in rest]))
            v = classify(stamp, *IMAGES[r["image"]], pred, args.search_arcsec)
            controls.append({"image": r["image"], "mu": r["mu"]} | v)
        for c, sp in zip(counter, spread, strict=True):
            mu_c = float(c["magnification"])
            row = {
                "model": name,
                "knot": knot,
                "z": zs[0],
                "ra": float(c["ra"]),
                "dec": float(c["dec"]),
                "x": float(c["x"]),
                "y": float(c["y"]),
                "mu": mu_c,
                "z_spread_arcsec": float(sp),
                "references": [r["image"] for r in good],
            }
            if len(unrecovered):  # distance to a catalogued image of this knot the model misses
                row["copy_distance_arcsec"] = float(
                    np.hypot(*(unrecovered - [row["x"], row["y"]]).T).min()
                )
            if abs(mu_c) < args.min_mu:
                row["class"] = "demagnified"
                rows.append(row)
                continue
            search = float(min(max(args.search_arcsec, sp), lc.MAX_SEARCH_ARCSEC))
            pred = (
                float(np.median([r["flux"] * abs(mu_c / r["mu"]) for r in good]))
                if good
                else float("nan")
            )
            row.update(classify(stamp, row["ra"], row["dec"], pred, search))
            rows.append(row)
    return {
        "model": name,
        "grid_nan_fraction": nan_frac,
        "checks": checks,
        "controls": controls,
        "counter_images": rows,
    }


def _finite(o):
    """NaN and inf as null, so the output is strict JSON."""
    if isinstance(o, dict):
        return {k: _finite(v) for k, v in o.items()}
    if isinstance(o, list | tuple):
        return [_finite(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=Path("results/sunrise/counter_images.json"))
    ap.add_argument("--models", nargs="+", default=list(MODELS), choices=MODELS)
    ap.add_argument("--half-width", type=float, default=50.0, help="solver box half-width, arcsec")
    ap.add_argument("--step", type=float, default=0.1, help="solver grid step, arcsec")
    ap.add_argument("--self-arcsec", type=float, default=2.0)
    ap.add_argument("--min-mu", type=float, default=0.5)
    ap.add_argument("--max-ref-mu", type=float, default=50.0)
    ap.add_argument("--search-arcsec", type=float, default=1.0)
    ap.add_argument(
        "--cross-arcsec",
        type=float,
        default=lc.MAX_SEARCH_ARCSEC,
        help="two models' counter-images of one knot within this are the same image",
    )
    ap.add_argument("--dz", type=float, default=0.3, help="system 2 redshift half-range")
    ap.add_argument("--relics", type=Path, default=Path("results/sunrise/relics_positions.json"))
    ap.add_argument("--image", default=F200W)
    args = ap.parse_args(argv)
    zbest = {s["model"]: s["sys2_best_z"] for s in json.loads(args.relics.read_text())["summary"]}
    raw_stamp, close = lc.image_stamper(
        args.image, lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.1
    )
    cache: dict = {}

    def stamp(ra: float, dec: float):  # catalogued images are measured once, not per model
        key = (round(ra, 7), round(dec, 7))
        if key not in cache:
            cache[key] = raw_stamp(ra, dec)
        return cache[key]

    try:
        # fetch every model's maps concurrently (I/O), then solve one model at a time
        with ThreadPoolExecutor(MODEL_THREADS) as ex:
            models = list(ex.map(lambda m: load_model(m, args.half_width), args.models))
        results = [
            run_model(m, model, args, zbest, stamp)
            for m, model in zip(args.models, models, strict=True)
        ]
    finally:
        close()
    cross_model(results, args.cross_arcsec)
    out = {
        "source": "RELICS v1 lens models, HLSP DOI 10.17909/T9SP45, accessed 2026-10-14",
        "access": "HTTP byte ranges of full-width map rows over the solver box; nothing stored",
        "positions": "Scofield+2025 Table 3 (observed, JWST frame)",
        "image": args.image,
        "provenance": (
            "image positions and mu: model_prediction; fluxes and classes: derived; "
            "thresholds are ASSUMPTIONs"
        ),
        "assumptions": {
            k: getattr(args, k)
            for k in (
                "half_width",
                "step",
                "self_arcsec",
                "min_mu",
                "max_ref_mu",
                "search_arcsec",
                "cross_arcsec",
                "dz",
            )
        }
        | {"centre": CENTRE},
        "results": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(_finite(out), indent=1) + "\n")
    print("flags:", sum(c["flag"] for r in results for c in r["counter_images"]))
    for r in results:
        for c in r["checks"] + r["controls"]:
            print(r["model"], json.dumps(c))
        for c in r["counter_images"]:
            print(
                r["model"],
                c["knot"],
                f"mu={c['mu']:.2f}",
                c["class"],
                f"n_models={c['n_models']}",
                f"inject={c.get('inject_class', '-')}",
                f"pred_snr={c.get('pred_snr', float('nan')):.1f}",
                f"best_snr={c.get('best_snr', float('nan')):.1f}",
                f"({c['ra']:.6f}, {c['dec']:.6f})",
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
