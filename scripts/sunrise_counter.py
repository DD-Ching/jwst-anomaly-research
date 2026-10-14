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
* model copy: a counter-image within ``COPY_FACTOR`` x the knot's image-plane rms of a catalogued
  image the model does not re-predict is the model's displaced copy of it, not a missing image;
* cross-model: ``n_models`` counts the models predicting a counter-image of the same knot within
  ``--cross-arcsec``; an image predicted by only some models is model-dependent. A flag is an
  ``absent`` counter-image that every model predicts.

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
from sunrise_relics import BASE, DTYPES, KEYS, MODELS, _get, read_header
from sunrise_scofield import IMAGES, system_of

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

Z_LENS = 0.566
COSMO = FlatLambdaCDM(H0=70.0, Om0=0.3)
CENTRE = (24.3555, -8.4570)  # model-frame origin near the cluster core (RA, Dec; ASSUMPTION)
F200W = "s3://stpubdata/jwst/public/jw02282/L3/t/o010/jw02282-o010_t001_nircam_clear-f200w_i2d.fits"
CHUNK = 8 * 2**20  # bytes per range request
COPY_FACTOR = 1.75  # model copy within this x the knot's image-plane rms (D-041; ASSUMPTION)
FLAG_CLASSES = ("absent",)


def cross_model(results: list[dict], radius: float) -> None:
    """Add ``n_models``: how many models predict a counter-image of the same knot within
    ``radius`` of each counter-image (D-041: an image predicted by only some models is
    model-dependent and untestable). Flags are ``absent`` images that every model predicts."""
    pts = {
        (r["model"], c["knot"]): np.array(
            [[q["ra"], q["dec"]] for q in r["counter_images"] if q["knot"] == c["knot"]]
        )
        for r in results
        for c in r["counter_images"]
    }
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
            c["flag"] = c["class"] in FLAG_CLASSES and n == len(results)


def read_block(session, url: str, wcs_hdr, offset: int, rows: tuple[int, int]) -> np.ndarray:
    """Full-width rows [rows[0], rows[1]) of one map, read as parallel byte ranges."""
    nx, bpp = wcs_hdr["NAXIS1"], abs(wcs_hdr["BITPIX"]) // 8
    start, stop = offset + rows[0] * nx * bpp, offset + rows[1] * nx * bpp
    spans = [(s, min(s + CHUNK, stop)) for s in range(start, stop, CHUNK)]
    with ThreadPoolExecutor(8) as ex:
        raw = b"".join(ex.map(lambda sp: _get(session, url, *sp), spans))
    return np.frombuffer(raw, DTYPES[wcs_hdr["BITPIX"]]).reshape(rows[1] - rows[0], nx)


def load_model(name: str, half: float) -> lensmodel.MapLensModel:
    """The model's deflection maps over a (2 half + margin) box around ``CENTRE``."""
    session = requests.Session()
    out = {}
    for k in ("defx", "defy"):
        url = BASE.format(m=name, k=KEYS[k])
        hdr, off = read_header(session, url)
        wcs = WCS(hdr)
        pix = float(abs(wcs.pixel_scale_matrix[1, 1]) * 3600.0)
        ic, jc = (float(v) for v in wcs.world_to_pixel_values(*CENTRE))
        m = int(np.ceil((half + 2.0) / pix))  # 2" margin for the finite differences
        i0, i1 = max(0, int(ic) - m), min(hdr["NAXIS1"], int(ic) + m + 1)
        j0, j1 = max(0, int(jc) - m), min(hdr["NAXIS2"], int(jc) + m + 1)
        out[k] = read_block(session, url, hdr, off, (j0, j1))[:, i0:i1]
        out["wcs"] = wcs[j0:j1, i0:i1]
        out["box"] = (i0, i1, j0, j1, pix)
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

    Known images lie within ``self_arcsec`` of a catalogued image of the system, counter-images
    do not. ``recovered`` marks the knot's catalogued images with a prediction within
    ``self_arcsec``; the rms is that of each catalogued image's distance to its nearest
    prediction (the model's image-plane rms for this knot, derived)."""
    imgs = pred["images"]
    xy = np.column_stack([imgs["x"], imgs["y"]])
    if not len(xy):
        n = len(pred["ref_xy"])
        return imgs, imgs, np.zeros(n, bool), float("inf")
    d = np.hypot(*(xy[:, None, :] - sys_xy[None, :, :]).transpose(2, 0, 1))
    near = d.min(axis=1) <= self_arcsec
    dk = np.hypot(*(pred["ref_xy"][:, None, :] - xy[None, :, :]).transpose(2, 0, 1)).min(axis=1)
    return imgs[near], imgs[~near], dk <= self_arcsec, float(np.sqrt(np.mean(dk**2)))


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


def run_model(name: str, args, zbest: dict, stamp) -> dict:
    model = load_model(name, args.half_width)
    grid = lensmodel.DeflectionGrid.compute(model, args.half_width, args.step)
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
        for r in good:  # positive control: leave-one-out prediction of a catalogued image
            rest = [q for q in good if q is not r]
            if not rest:
                continue
            pred = float(np.median([q["flux"] * abs(r["mu"] / q["mu"]) for q in rest]))
            ratio = r["flux"] / pred
            controls.append(
                {
                    "image": r["image"],
                    "mu": r["mu"],
                    "pred_snr": pred / r["err"],
                    "best_snr": r["flux"] / r["err"],
                    "flux_ratio": ratio,
                    "class": lc.forced_class(pred / r["err"], r["flux"] / r["err"], ratio),
                }
            )
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
            # D-041 model copy: near a catalogued image of this knot that the model misses
            if len(unrecovered):
                dcopy = float(np.hypot(*(unrecovered - [row["x"], row["y"]]).T).min())
                row["copy_distance_arcsec"] = dcopy
                if dcopy <= COPY_FACTOR * rms:
                    row["class"] = "model_copy"
                    rows.append(row)
                    continue
            if abs(mu_c) < args.min_mu:
                row["class"] = "demagnified"
                rows.append(row)
                continue
            search = float(min(max(args.search_arcsec, sp), lc.MAX_SEARCH_ARCSEC))
            sci, err, xx, yy = stamp(row["ra"], row["dec"])
            f0, e0 = lc.aperture_snr(sci, err, xx, yy, 0.0, 0.0)
            snr, fbest, dx, dy = lc.best_within(sci, err, xx, yy, search)
            pred = (
                float(np.median([r["flux"] * abs(mu_c / r["mu"]) for r in good]))
                if good
                else float("nan")
            )
            pred_snr = pred / e0 if np.isfinite(e0) and e0 > 0 else float("nan")
            ratio = fbest / pred if np.isfinite(pred) and pred > 0 else float("nan")
            row.update(
                {
                    "search_arcsec": search,
                    "pred_flux": pred,
                    "pred_snr": pred_snr,
                    "best_snr": float(snr),
                    "best_flux": float(fbest),
                    "best_offset_arcsec": float(np.hypot(dx, dy)),
                    "flux_ratio": ratio,
                    "class": "no_reference"
                    if not good
                    else lc.forced_class(
                        pred_snr, float(snr), ratio if np.isfinite(ratio) else 0.0
                    ),
                }
            )
            rows.append(row)
    return {"model": name, "checks": checks, "controls": controls, "counter_images": rows}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=Path("results/sunrise/counter_images.json"))
    ap.add_argument("--models", nargs="+", default=list(MODELS))
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
    stamp, close = lc.image_stamper(args.image, lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.1)
    try:
        results = [run_model(m, args, zbest, stamp) for m in args.models]
    finally:
        close()
    cross_model(results, args.cross_arcsec)
    out = {
        "source": "RELICS v1 lens models, HLSP DOI 10.17909/T9SP45, accessed 2026-10-14",
        "access": "HTTP byte ranges of the deflection maps over the solver box; nothing stored",
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
        | {"centre": CENTRE, "copy_factor": COPY_FACTOR},
        "results": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
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
                f"pred_snr={c.get('pred_snr', float('nan')):.1f}",
                f"best_snr={c.get('best_snr', float('nan')):.1f}",
                f"({c['ra']:.6f}, {c['dec']:.6f})",
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
