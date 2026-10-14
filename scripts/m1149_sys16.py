"""MACS1149 system 16: forced photometry along the predicted images' redshift track.

System 16 (CATS ``arcs.txt``: 16.1, 16.2, 0.85" apart) has no spectroscopic redshift; CATS fits
z 4.419, the CANUCS DR1 photo-z of 16.2 is 2.25 (68 % 1.90-2.28). Both CATS v4.1 and Sharon v4cor
predict a third, "far" image ~10" south-west of the pair whose position moves several arcsec with
the source redshift (docs/fields/macs1149.md, Limits). This script scans the source redshift and,
for each model and z:

- solves the images of the magnification-weighted least-squares source of 16.1/16.2
  (``m1149_refs.second_model_images``) and assigns the two predicted images nearest to 16.1 and
  16.2 (distinct; ``pair_rms_arcsec`` is their image-plane rms, an observed-vs-model check);
- every other image is an "extra" (``far`` when more than ``Params.far_arcsec`` from the pair);
- predicts its flux from the measured pair flux, ``(f1 + f2) * |mu_x| / (|mu_1| + |mu_2|)``,
  with the model's |mu| of the two assigned images, and its S/N against the aperture error at
  the predicted position, in F150W, F277W and F444W on 1"-median-filtered stamps;
- measures the fixed-position S/N and the best aperture within ``Params.search_arcsec`` and applies
  ``lens_consistency.forced_class``.

A redshift is ``excluded`` for a model when an extra is ``absent`` in a band, or the pair is not
reproduced (``pair_rms_arcsec`` > ``Params.pair_rms_max``). A track position whose F277W best
aperture has S/N >= 5, with F150W and F444W measured at that same offset (S/N >=
``Params.colour_min_snr``) and colours matching the pair's, is listed as a possible counter-image
(``colour_matches``).

Provenance: positions and mu ``model_prediction``; 16.1/16.2 positions ``observed``; fluxes
``derived`` from S3 byte-range stamps of CANUCS program 1208. Thresholds are ASSUMPTIONs.

Usage: python scripts/m1149_sys16.py --sharon-x <fits> --sharon-y <fits> [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a2744_uncover as au  # noqa: E402
import lens_consistency as lc  # noqa: E402
import m1149_refs as mr  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

SYSTEM = "16"
Z_GRID = np.round(np.arange(2.0, 6.0001, 0.05), 3)  # source redshifts scanned


@dataclass(frozen=True)
class Params:
    """Thresholds (all ASSUMPTIONs)."""

    far_arcsec: float = 3.0  # an extra this far from both pair images is the "far" image
    search_arcsec: float = 0.7  # best-aperture radius around a predicted extra (~CATS rms 0.63")
    pair_rms_max: float = 1.0  # pair image-plane rms above this: z does not reproduce 16.1/16.2
    colour_tol_mag: float = 0.5  # |colour - pair colour| for a colour match, per colour
    colour_min_snr: float = 3.0  # every band at the common offset needs this S/N for a colour


P = Params()
STAMP_HALF_ARCSEC = 2.0
GRID_STEP = 0.25
DEFAULT_OUT = mr.DEFAULT_OUT
_r = au._r


def assign_pair(pred: Table, cat: Table) -> tuple[list[int], float]:
    """Indices of the predicted images assigned to the catalogued images (distinct, nearest
    first by the cheaper of the two assignments) and their image-plane rms (arcsec).
    ``cat`` must hold exactly the two pair images."""
    if len(cat) != 2:
        raise ValueError(f"assign_pair needs 2 catalogued images, got {len(cat)}")
    if len(pred) < len(cat):
        return [], float("inf")
    d = np.array(
        [au._sep_arcsec(pred["ra"], pred["dec"], float(c["ra"]), float(c["dec"])) for c in cat]
    )
    best, best_cost = None, np.inf
    n = len(pred)
    for i in range(n):  # two catalogued images: brute force over ordered pairs
        for j in range(n):
            if i == j:
                continue
            cost = d[0, i] ** 2 + d[1, j] ** 2
            if cost < best_cost:
                best, best_cost = [i, j], cost
    return best, float(np.sqrt(best_cost / 2))


def predicted_flux(pair_flux: float, mu_pair: tuple[float, float], mu_x: float) -> float:
    """Flux of an extra image from the summed pair flux and the model's |mu| of all three."""
    s = abs(mu_pair[0]) + abs(mu_pair[1])
    return pair_flux * abs(mu_x) / s if s > 0 and np.isfinite(s) else float("nan")


def colours(f: dict[str, float]) -> tuple[float, float]:
    """(F150W - F277W, F277W - F444W) in mag from positive fluxes, else NaN."""

    def c(a, b):
        return -2.5 * np.log10(f[a] / f[b]) if f[a] > 0 and f[b] > 0 else float("nan")

    return c("f150w", "f277w"), c("f277w", "f444w")


def z_verdict(pair_rms: float, classes: list[str]) -> str:
    if not pair_rms <= P.pair_rms_max:
        return "excluded: pair not reproduced"
    if "absent" in classes:
        return "excluded: extra image absent"
    return "allowed"


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sharon-x", type=Path, required=True)
    ap.add_argument("--sharon-y", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)
    cats, files, par = lc.load_model("macs1149-cats")
    arcs, _ = lc.image_list("macs1149-cats", files, par)
    cat = arcs[np.asarray(arcs["system"]).astype(str) == SYSTEM]
    models = {"cats": cats, "sharon": mr.sharon_model(args.sharon_x, args.sharon_y)}
    opened = {b: lc.image_stamper(mr.S3_I2D.format(b), STAMP_HALF_ARCSEC) for b in mr.BANDS}
    stamps = {b: lc.halo_stamper(o[0], au.HALO_FILTER_ARCSEC) for b, o in opened.items()}
    try:
        pair, pair_flux = {}, {}
        for b, st in stamps.items():
            sib, tot = {}, 0.0
            for c in cat:
                img, err, xx, yy = st(float(c["ra"]), float(c["dec"]))
                f, e = lc.peak_flux(img, err, xx, yy)
                tot += f  # NaN propagates: no prediction in that band
                sib[str(c["image_id"])] = {
                    "flux": _r(f, 4),
                    "snr": _r(f / e if e > 0 else np.nan, 1),
                }
            pair[b], pair_flux[b] = sib, tot
        pair_col = colours(pair_flux)
        rows = []
        for name, model in models.items():
            grid = lensmodel.DeflectionGrid.compute(
                model, lc.grid_half_width(model, cat), GRID_STEP
            )
            for z in Z_GRID:
                pred = mr.second_model_images(model, grid, cat, float(z))
                idx, rms = assign_pair(pred, cat)
                row = {
                    "model": name,
                    "z": float(z),
                    "n_images": len(pred),
                    "pair_rms_arcsec": _r(rms),
                }
                mu_pair = (
                    tuple(float(pred["magnification"][k]) for k in idx) if idx else (np.nan, np.nan)
                )
                row["mu_pair"] = [_r(m, 1) for m in mu_pair]
                extras = []
                for k in range(len(pred)):
                    if k in idx:
                        continue
                    ra, dec = float(pred["ra"][k]), float(pred["dec"][k])
                    mu_x = float(pred["magnification"][k])
                    sep = float(np.min(au._sep_arcsec(cat["ra"], cat["dec"], ra, dec)))
                    ex = {
                        "radec": [round(ra, 6), round(dec, 6)],
                        "mu": _r(mu_x, 1),
                        "sep_pair_arcsec": _r(sep),
                        "kind": "far" if sep > P.far_arcsec else "near",
                    }
                    with ThreadPoolExecutor(len(stamps)) as pool:  # one range reader per band
                        cut = dict(
                            zip(
                                stamps,
                                pool.map(lambda st, ra=ra, dec=dec: st(ra, dec), stamps.values()),
                                strict=True,
                            )
                        )
                    # colours at one offset: the F277W best aperture (as the pair's peak apertures)
                    snr_c, _, cx, cy = lc.best_within(*cut["f277w"], P.search_arcsec)
                    common = {}
                    for b, (img, err, xx, yy) in cut.items():
                        f0, e0 = lc.aperture_snr(img, err, xx, yy, 0.0, 0.0)
                        snr, flux, dx, dy = lc.best_within(img, err, xx, yy, P.search_arcsec)
                        pf = predicted_flux(pair_flux[b], mu_pair, mu_x)
                        psnr = pf / e0 if e0 > 0 else np.nan
                        ratio = flux / pf if pf > 0 else np.nan
                        fc, ec = lc.aperture_snr(img, err, xx, yy, cx, cy)
                        ok = snr_c >= 5 and ec > 0 and fc / ec >= P.colour_min_snr
                        common[b] = fc if ok else np.nan
                        ex[b] = {
                            "pred_snr": _r(psnr, 1),
                            "fixed_snr": _r(f0 / e0 if e0 > 0 else np.nan, 1),
                            "best_snr": _r(snr, 1),
                            "best_offset": [_r(dx, 1), _r(dy, 1)],
                            "flux_over_pred": _r(ratio),
                            "forced_class": lc.forced_class(psnr, snr, ratio),
                        }
                    col = colours(
                        {b: common[b] if np.isfinite(common[b]) else -1 for b in mr.BANDS}
                    )
                    ex["colour_offset"] = [_r(cx, 1), _r(cy, 1)]
                    ex["best_colours"] = [_r(c) for c in col]
                    ex["colour_match"] = bool(
                        all(np.isfinite(col))
                        and all(
                            abs(a - p) <= P.colour_tol_mag
                            for a, p in zip(col, pair_col, strict=True)
                        )
                    )
                    extras.append(ex)
                row["extras"] = extras
                classes = [e[b]["forced_class"] for e in extras for b in mr.BANDS]
                row["verdict"] = z_verdict(rms, classes)
                rows.append(row)
    finally:
        for _, close in opened.values():
            close()
    allowed = {
        n: [r["z"] for r in rows if r["model"] == n and r["verdict"] == "allowed"] for n in models
    }
    matches = [
        {"model": r["model"], "z": r["z"], **{k: e[k] for k in ("radec", "kind", "best_colours")}}
        for r in rows
        for e in r["extras"]
        if e["colour_match"]
    ]
    summary = {
        "provenance": "positions and mu: model_prediction (CATS v4.1, Sharon v4cor); 16.1/16.2 "
        "positions: observed; fluxes: derived",
        "images": {b: mr.S3_I2D.format(b) for b in mr.BANDS},
        "assumptions": {
            **asdict(P),
            "halo_filter_arcsec": au.HALO_FILTER_ARCSEC,
            "source": "magnification-weighted least-squares source of 16.1/16.2",
            "predicted_flux": "(f1 + f2) |mu_x| / (|mu_1| + |mu_2|), model mu of the assigned "
            "images",
        },
        "pair": pair,
        "pair_colours": [_r(c) for c in pair_col],
        "allowed_z": allowed,
        "colour_matches": matches,
        "rows": rows,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "sys16_track.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


if __name__ == "__main__":
    s = main()
    print("pair", s["pair"], "colours", s["pair_colours"])
    for r in s["rows"]:
        print(
            r["model"],
            r["z"],
            r["pair_rms_arcsec"],
            r["mu_pair"],
            r["verdict"],
            *(
                f"{e['kind']} mu={e['mu']} "
                + " ".join(
                    f"{b}:{e[b]['forced_class']}({e[b]['pred_snr']}/{e[b]['best_snr']})"
                    for b in mr.BANDS
                )
                + f" cm={e['colour_match']}"
                for e in r["extras"]
            ),
        )
    print("allowed", s["allowed_z"])
    print("colour matches", len(s["colour_matches"]))
