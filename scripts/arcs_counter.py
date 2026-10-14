"""Counter-images of uncatalogued bright arcs: does every arc have the images its model predicts?

``lens_consistency.py images`` (D-029) tests the predicted images of *catalogued* systems. Here
the sources are the strongly sheared background arcs of a stored ``arcs`` run that match no
catalogued multiple image. Each arc's catalogue centre is traced to the source plane at its
photo-z (z16, z50, z84 of the stored row) and ``lensmodel.find_images`` predicts every image.
Images within ``--self-arcsec`` of the arc are the arc itself (a merging pair across a critical
curve); the rest are counter-images (model_prediction). A counter-image's predicted total
flux is the arc's isophotal flux (observed) times |mu_c / mu_arc| (model).

Each counter-image predicted at z50 with |mu| >= ``--min-mu`` is then checked (derived):
- ``catalog``: a catalogue source within the search radius no fainter than the predicted
  magnitude + ``--mag-slack`` (the D-029 1/3 flux bound is 1.19 mag);
- forced photometry otherwise (D-029 aperture, ERR scale, classes): the predicted aperture
  flux is the total times ``--ee`` (encircled fraction, ASSUMPTION: 0.5 is low for compact
  counter-images, so ``absent`` is conservative). The search radius is the largest offset of
  the matching z16 / z84 counter-image from the z50 one, at least ``--search-arcsec`` and at most
  ``lens_consistency.MAX_SEARCH_ARCSEC``.

A counter-image is ``robust`` when z16 and z84 both predict one within that capped radius. Flags
are robust, absent images (predicted >= 10 sigma) with nothing detected: best S/N below
``--max-best-snr``. The isophotal-to-aperture flux prediction over-predicts bright extended
arcs by up to x4 (Abell 2744 control 7247: a catalogued image detected at 86 sigma is
``absent`` by the D-029 1/3 flux bound), so a detected but faint source is not a flag. Arcs
within ``--match-cat-arcsec`` (2") of a catalogued image are excluded as knots of that system:
Abell 2744 arc 2528 lies 1.1" from 26.1a, and its photo-z (5.16 against z_spec 3.054) moved the
predicted counter-images by up to 8". Arcs with |mu_arc| > ``--max-arc-mu``
are skipped: near a critical curve mu_arc, and so every flux prediction, is unreliable.
``--control`` runs the same route on the arcs that are catalogued multiple images, whose
counter-images are known to exist (positive control; outputs under ``control/``).
Outputs ``results/arcs_counter/<model>/counter.ecsv`` and ``counter.json``.

    python scripts/arcs_counter.py --model abell2744-bergamini23
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack

from jwst_anomaly import lensmodel, paths, schema
from jwst_anomaly.photometry import fetch_catalog

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ab = _load("arcs_bands")
ac = _load("arcs_curvature")
lc = ab.lc

ASSUMPTIONS = (
    "self_arcsec",
    "min_mu",
    "max_arc_mu",
    "match_cat_arcsec",
    "search_arcsec",
    "mag_slack",
    "ee",
    "max_best_snr",
)


def catalogued_mask(arcs: Table, images: Table, radius: float) -> np.ndarray:
    """True for arcs within ``radius`` arcsec of a catalogued multiple image."""
    if not len(images):
        return np.zeros(len(arcs), bool)
    ca = SkyCoord(arcs["ra"], arcs["dec"], unit="deg")
    ci = SkyCoord(images["ra"], images["dec"], unit="deg")
    _, sep, _ = ca.match_to_catalog_sky(ci)
    return sep.arcsec <= radius


def known_images(name: str, files: dict, par) -> Table:
    """Catalogued multiple images of the cluster in the JWST frame, for excluding arcs that are
    knots of known systems: the model's own list even where D-035/D-044 gate it off for position
    tests (identification needs only the observed positions), plus the HFF CATS list of the same
    cluster when ``lens_consistency.MODELS`` has one. Abell 370 3881 is CATS image 7.5 (system 7,
    z_spec 2.751) but matches no CANUCS image within 7", and its z_phot 2.99 moved the predicted
    counter-image 2.6" off the detected image 7.4."""
    names = [name] + [c for c in (f"{name.split('-')[0]}-cats",) if c != name and c in lc.MODELS]
    parts = []
    for n in names:
        if "arcs.dat" not in lc.MODELS[n]["files"]:
            continue
        fs = files if n == name else {}
        imgs, _ = lc.image_list(n, fs, par if n == name else None)
        parts.append(lc.shift_images(n, imgs)["ra", "dec"])
    return vstack(parts, metadata_conflicts="silent") if parts else Table({"ra": [], "dec": []})


def counter_images(model, grid, ra: float, dec: float, z: float, self_arcsec: float) -> Table:
    """Images of the arc's source at ``z`` farther than ``self_arcsec`` from the arc, with
    ``mu_ratio`` = |mu_image / mu_arc| (mu_arc from the model at the arc centre)."""
    x0, y0 = (float(v) for v in model.to_frame(ra, dec))
    bx, by = model.source_points(np.array([x0]), np.array([y0]), z)
    pred = lensmodel.find_images(model, grid, float(bx[0]), float(by[0]), z)
    mu_arc = float(model.magnification(ra, dec, z)["magnification"][0])
    if not len(pred):
        return pred, mu_arc
    far = np.hypot(np.asarray(pred["x"]) - x0, np.asarray(pred["y"]) - y0) > self_arcsec
    out = pred[far]
    out["mu_ratio"] = np.abs(np.asarray(out["magnification"], float)) / abs(mu_arc)
    return out, mu_arc


def nearest_offset(x: float, y: float, other: Table) -> float:
    if other is None or not len(other):
        return float("inf")
    return float(np.min(np.hypot(np.asarray(other["x"]) - x, np.asarray(other["y"]) - y)))


def catalog_match(shapes: Table, ra: float, dec: float, radius: float, mag_max: float):
    """Brightest catalogue source within ``radius`` arcsec, and whether it is bright enough."""
    c = SkyCoord(shapes["ra"], shapes["dec"], unit="deg")
    sep = c.separation(SkyCoord(ra, dec, unit="deg")).arcsec
    near = sep <= radius
    if not near.any():
        return -1, np.nan, np.nan, False
    mags = np.where(near, np.asarray(shapes["mag"], float), np.inf)
    k = int(np.argmin(mags))
    return int(shapes["label"][k]), float(mags[k]), float(sep[k]), bool(mags[k] <= mag_max)


def mag_to_image_units(mag: float, pixar_sr: float) -> float:
    """AB magnitude to a sum of MJy/sr pixels (the i2d SCI unit)."""
    return 3631.0 * 10 ** (-0.4 * mag) / 1e6 / pixar_sr


def pixar_sr(uri: str) -> float:
    import fsspec

    opts = {"anon": True} if uri.startswith("s3://") else {}
    fs, path = fsspec.core.url_to_fs(uri, **opts)
    with fs.open(path, "rb", block_size=2**18) as fo, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with fits.open(fo, lazy_load_hdus=True) as h:
            return float(h["SCI"].header["PIXAR_SR"])


def nearest_image(model, ra: float, dec: float, z: float, x: float, y: float, half: float, step):
    """Offset (arcsec) and signed mu of the image of the source at (ra, dec) nearest to model
    position (x, y) within ``half`` arcsec, from a local grid refined by Newton steps; NaN when
    no image lies there (source-plane misfit > 0.01")."""
    x0, y0 = (float(v) for v in model.to_frame(ra, dec))
    b = model.source_points(np.array([x0]), np.array([y0]), z)
    g = np.arange(-half, half + step / 2, step)
    gx, gy = np.meshgrid(x + g, y + g)
    bx, by = model.source_points(gx, gy, z)
    k = int(np.argmin(np.hypot(bx - b[0][0], by - b[1][0])))
    tx, ty = float(gx.flat[k]), float(gy.flat[k])
    for _ in range(12):
        m = model.lens_map(np.array([tx]), np.array([ty]), z)
        r = np.array([m["beta_x"][0] - b[0][0], m["beta_y"][0] - b[1][0]])
        a = np.array([[m["a11"][0], m["a12"][0]], [m["a21"][0], m["a22"][0]]])
        try:
            d = np.linalg.solve(a, r)
        except np.linalg.LinAlgError:
            break
        tx, ty = tx - float(d[0]), ty - float(d[1])
    m = model.lens_map(np.array([tx]), np.array([ty]), z)
    miss = float(np.hypot(m["beta_x"][0] - b[0][0], m["beta_y"][0] - b[1][0]))
    off = float(np.hypot(tx - x, ty - y))
    if miss > 0.01 or off > half * np.sqrt(2):
        return float("nan"), float("nan")
    det = m["a11"][0] * m["a22"][0] - m["a12"][0] * m["a21"][0]
    return off, float(1.0 / det)


def chain_check(name: str, par: dict, flags: Table, n_samples: int, half: float = 3.0) -> list:
    """For each flagged counter-image: over ``n_samples`` random rows of the model's
    ``bayes.dat``, the fraction that still predicts an image within ``half`` arcsec of it, and
    those images' offsets and |mu| (model_prediction; D-045 chain)."""
    bayes, _, _, reference = lc.load_chain(name, par)
    rows = np.random.default_rng(1).choice(len(bayes), n_samples, replace=False)
    res = {i: [] for i in range(len(flags))}
    for row in rows:
        m = lensmodel.LensModel.from_par(lensmodel.posterior_par(par, bayes, int(row), reference))
        lc.apply_frame_offset(name, m)
        for i, f in enumerate(flags):
            res[i].append(
                nearest_image(m, f["arc_ra"], f["arc_dec"], f["z_phot"], f["x"], f["y"], half, 0.05)
            )
    out = []
    for i, f in enumerate(flags):
        a = np.array(res[i], float)
        ok = np.isfinite(a[:, 0])
        out.append(
            {
                "label": int(f["label"]),
                "n_samples": n_samples,
                "frac_with_image": float(ok.mean()),
                "offset_arcsec_p16_50_84": lc._pct(a[ok, 0]) if ok.any() else None,
                "abs_mu_p16_50_84": lc._pct(np.abs(a[ok, 1])) if ok.any() else None,
            }
        )
    return out


def run(args) -> dict:
    arcs = Table.read(args.arcs_dir / args.model / "arcs.ecsv")
    model, files, par = lc.load_model(args.model)
    lc.apply_frame_offset(args.model, model)
    images = known_images(args.model, files, par)

    bg = (arcs["population"] == "background") & np.asarray(arcs["strong_shear"], bool)
    bg &= np.isfinite(np.asarray(arcs["z_phot"], float)) & np.isfinite(
        np.asarray(arcs["mag"], float)
    )
    cat_img = catalogued_mask(arcs, images, args.match_cat_arcsec)
    sel = arcs[bg & (cat_img if args.control else ~cat_img)]

    xs, ys = model.to_frame(np.asarray(sel["ra"]), np.asarray(sel["dec"]))
    half = float(np.max(np.hypot(xs, ys))) + 25.0 if len(sel) else 60.0
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, half, args.step), half, args.step
    )
    cat_name = ac.IMAGES[args.model].replace("_i2d.fits", "_cat.ecsv")
    manifest = ab.FIELDS[args.model][0]
    shapes = lc.load_shapes(fetch_catalog(ab.s3_url(cat_name), ab.manifest_sha(manifest, cat_name)))
    uri = ac.image_uri(ac.IMAGES[args.model])
    pix_sr = pixar_sr(uri)
    stamp, close = lc.image_stamper(uri, lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.1)

    rows, skipped = [], {"arc_mu": 0, "no_counter": 0}
    z_lens = float(arcs.meta.get("z_lens", 0.3))
    try:
        for a in sel:
            zs = {
                "z50": float(a["z_phot"]),
                "z16": float(a["z160"]),
                "z84": float(a["z840"]),
            }
            zs = {k: max(v, z_lens + 0.05) for k, v in zs.items() if np.isfinite(v)}
            preds = {
                k: counter_images(model, grid, a["ra"], a["dec"], z, args.self_arcsec)
                for k, z in zs.items()
            }
            pred, mu_arc = preds["z50"]
            if not np.isfinite(mu_arc) or abs(mu_arc) > args.max_arc_mu:
                skipped["arc_mu"] += 1
                continue
            pred = pred[np.abs(np.asarray(pred["magnification"], float)) >= args.min_mu]
            if not len(pred):
                skipped["no_counter"] += 1
                continue
            for p in pred:
                x, y = float(p["x"]), float(p["y"])
                offs = [nearest_offset(x, y, preds[k][0]) for k in ("z16", "z84") if k in preds]
                spread = max(offs) if offs else float("inf")
                robust = bool(offs) and spread <= lc.MAX_SEARCH_ARCSEC
                radius = float(np.clip(spread, args.search_arcsec, lc.MAX_SEARCH_ARCSEC))
                m_pred = float(a["mag"]) - 2.5 * np.log10(float(p["mu_ratio"]))
                lab, m_near, sep_near, cat_ok = catalog_match(
                    shapes, float(p["ra"]), float(p["dec"]), radius, m_pred + args.mag_slack
                )
                sci, err, xx, yy = stamp(float(p["ra"]), float(p["dec"]))
                _, e0 = lc.aperture_snr(sci, err, xx, yy, 0.0, 0.0)
                pred_snr = best_snr = flux_ratio = np.nan
                if cat_ok:
                    cls = "catalog"
                elif not (np.isfinite(e0) and e0 > 0):
                    cls = "off_image"
                else:
                    f_pred = args.ee * mag_to_image_units(m_pred, pix_sr)
                    pred_snr = f_pred / e0
                    best_snr, flux, _, _ = lc.best_within(sci, err, xx, yy, radius)
                    flux_ratio = flux / f_pred
                    cls = lc.forced_class(pred_snr, best_snr, flux_ratio)
                rows.append(
                    {
                        "label": int(a["label"]),
                        "arc_ra": float(a["ra"]),
                        "arc_dec": float(a["dec"]),
                        "arc_mag": float(a["mag"]),
                        "z_phot": zs["z50"],
                        "mu_arc": mu_arc,
                        "ra": float(p["ra"]),
                        "dec": float(p["dec"]),
                        "x": x,
                        "y": y,
                        "mu": float(p["magnification"]),
                        "mag_pred": m_pred,
                        "z_spread_arcsec": spread,
                        "robust": robust,
                        "search_arcsec": radius,
                        "nearest_label": lab,
                        "nearest_mag": m_near,
                        "nearest_sep_arcsec": sep_near,
                        "pred_snr": pred_snr,
                        "best_snr": best_snr,
                        "flux_ratio": flux_ratio,
                        "counter_class": cls,
                    }
                )
    finally:
        close()

    t = Table(rows=rows) if rows else Table()
    t.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        model=args.model,
        image=uri,
        catalog=cat_name,
        assumptions={k: getattr(args, k) for k in ASSUMPTIONS},
        forced={
            "r_ap_arcsec": lc.FORCED_R_AP,
            "annulus_arcsec": list(lc.FORCED_ANNULUS),
            "err_scale": lc.ERR_SCALE,
        },
    )
    classes = {}
    if len(t):
        for c in np.unique(t["counter_class"]):
            classes[str(c)] = int(np.sum(t["counter_class"] == c))
    flags = (
        t[(t["counter_class"] == "absent") & t["robust"] & (t["best_snr"] < args.max_best_snr)]
        if len(t)
        else t
    )
    summary = {
        "model": args.model,
        "assumptions": {k: getattr(args, k) for k in ASSUMPTIONS},
        "control": bool(args.control),
        "n_background_strong_shear": int(np.sum(bg)),
        "n_catalogued_excluded": int(np.sum(bg & cat_img)),
        "n_arcs_tested": int(len(set(t["label"]))) if len(t) else 0,
        "skipped": skipped,
        "n_counter_images": len(t),
        "n_robust": int(np.sum(t["robust"])) if len(t) else 0,
        "classes": classes,
        "flags_robust_absent": [
            {
                "label": int(r["label"]),
                "xy": [round(float(r["x"]), 2), round(float(r["y"]), 2)],
                "mu": round(float(r["mu"]), 2),
                "mu_arc": round(float(r["mu_arc"]), 2),
                "mag_pred": round(float(r["mag_pred"]), 2),
                "pred_snr": round(float(r["pred_snr"]), 1),
                "best_snr": round(float(r["best_snr"]), 1),
            }
            for r in flags
        ],
    }
    if args.chain_samples and len(flags):
        summary["chain"] = chain_check(args.model, par, flags, args.chain_samples)
    out = args.out / args.model / ("control" if args.control else "")
    out.mkdir(parents=True, exist_ok=True)
    if len(t):
        t.write(out / "counter.ecsv", overwrite=True)
    (out / "counter.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--model", required=True, choices=sorted(ac.IMAGES))
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    p.add_argument("--out", type=Path, default=paths.repo_root() / "results" / "arcs_counter")
    p.add_argument(
        "--control",
        action="store_true",
        help="positive control: test the arcs that ARE catalogued multiple images instead",
    )
    p.add_argument(
        "--chain-samples",
        type=int,
        default=0,
        help="bayes.dat rows to re-test each robust absent flag on (models with a chain)",
    )
    p.add_argument("--step", type=float, default=0.1, help="deflection grid step, arcsec")
    p.add_argument("--self-arcsec", type=float, default=2.0)
    p.add_argument("--min-mu", type=float, default=0.5)
    p.add_argument("--max-arc-mu", type=float, default=30.0)
    p.add_argument("--match-cat-arcsec", type=float, default=2.0)
    p.add_argument("--search-arcsec", type=float, default=1.0)
    p.add_argument("--mag-slack", type=float, default=1.19)
    p.add_argument("--max-best-snr", type=float, default=5.0)
    p.add_argument("--ee", type=float, default=0.5)
    return p


def main(argv=None) -> int:
    summary = run(build_parser().parse_args(argv))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
