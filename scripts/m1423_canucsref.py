"""MACS1423: extra-image test with CANUCS DR1 reference fluxes and a local-noise rule.

Frozen plan: ``results/lens_consistency/macs1423_canucsref/PLAN.md``. Runs ``m0417_arcref.run``
on MACS1423 with the CANUCS DR1 catalogue mapped by ``m0417_canucsref.to_cat`` (rules unchanged),
then applies the plan's local-noise rule (rule 2) with apertures on the F200W BCG-subtracted
mosaic, and re-evaluates gate C2 and the flags on the updated classes.

    python scripts/m1423_canucsref.py --cat <photometry-cat.fits.gz> --image <f200w bgsub-sci>
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import m0417_arcref as ar
import m0417_canucsref as cr
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from m0417_canucsref_posthoc import _clean, aperture_sum

from jwst_anomaly.provenance import file_sha256

FIELD = "macs1423"
FRAME = Path("results/lens_consistency/macs0417_extra")
OUT = Path("results/lens_consistency/macs1423_canucsref")
BASE = "https://archive.stsci.edu/hlsps/canucs/dr1/macs1423/clu/"
CAT_NAME = "hlsp_canucs_jwst-hst_multi_macs1423-clu_multi_v1_photometry-cat.fits.gz"
IMAGE_NAME = "hlsp_canucs_jwst_nircam_macs1423-clu-40mas_f200w_v1_bgsub-sci.fits.gz"
# Frozen local-noise rule (PLAN.md rule 2; ASSUMPTIONs)
AP_ARCSEC, PIX = 0.3, 0.04
N_RAND, R_RAND = 200, (1.5, 4.0)
SRC_EXCL = 0.5
SEED = 20261022
MIN_SIGMA, ABSENT_SIGMA = 5.0, 10.0


def local_noise(d, x, y, src_x, src_y, rng) -> float:
    """sigma-clipped std of ``N_RAND`` random aperture sums 1.5-4" from (x, y) (pixels),
    redrawn within ``SRC_EXCL`` of a source or on NaN pixels."""
    r = AP_ARCSEC / PIX
    vals, tries = [], 0
    while len(vals) < N_RAND and tries < 50 * N_RAND:
        tries += 1
        a, s = rng.uniform(0, 2 * np.pi), rng.uniform(*R_RAND) / PIX
        rx, ry = x + s * np.cos(a), y + s * np.sin(a)
        if len(src_x) and np.hypot(src_x - rx, src_y - ry).min() * PIX < SRC_EXCL:
            continue
        v = aperture_sum(d, rx, ry, r)
        if np.isfinite(v):
            vals.append(v)
    if len(vals) < N_RAND:
        raise RuntimeError(f"only {len(vals)} of {N_RAND} random apertures accepted")
    return float(sigma_clipped_stats(vals)[2])


def local_class(cls: str, f_pred_njy: float, sigma_njy: float) -> str:
    """Plan rule 2 on a rule-1 class."""
    if cls not in ar.TESTABLE:
        return cls
    if not f_pred_njy >= MIN_SIGMA * sigma_njy:
        return "undetectable_local"
    if cls == "absent" and f_pred_njy < ABSENT_SIGMA * sigma_njy:
        return "ambiguous"
    return cls


def apply_local(tab: Table, d, w, src_x, src_y, rng) -> None:
    """Add ``sigma_loc_njy``, ``pred_over_sigma`` and ``class_local`` columns in place."""
    sig, cls = [], []
    for r in tab:
        f = float(r["f_pred"]) * 1e9
        if str(r["class"]) in ar.TESTABLE:
            x, y = (float(v) for v in w.world_to_pixel_values(r["ra"], r["dec"]))
            s = local_noise(d, x, y, src_x, src_y, rng)
        else:
            s = float("nan")
        sig.append(s)
        cls.append(local_class(str(r["class"]), f, s))
    tab["sigma_loc_njy"] = sig
    with np.errstate(divide="ignore", invalid="ignore"):
        tab["pred_over_sigma"] = np.asarray(tab["f_pred"], float) * 1e9 / np.asarray(sig)
    tab["class_local"] = cls


def _count(v) -> dict:
    u, n = np.unique(np.asarray(v, str), return_counts=True)
    return {str(a): int(b) for a, b in zip(u, n, strict=True)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--image", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.cat.name != CAT_NAME or args.image.name != IMAGE_NAME:
        raise ValueError("--cat / --image are not the CANUCS DR1 MACS1423 files named in the plan")
    with fits.open(args.image) as hl:
        if hl[0].header.get("BUNIT") != "nJy":
            raise ValueError(f"unexpected BUNIT {hl[0].header.get('BUNIT')!r}")
        w, d = WCS(hl[0].header), np.asarray(hl[0].data, float)
    if abs(proj_plane_pixel_scales(w.celestial).mean() * 3600.0 - PIX) > 1e-4:
        raise ValueError("image pixel scale is not 0.04 arcsec")
    t = Table.read(args.cat)
    cat = cr.to_cat(t)
    inputs = {
        "canucs_cat": args.cat.name,
        "canucs_url": BASE + CAT_NAME,
        "canucs_cat_sha256": file_sha256(args.cat),
        "n_rows_finite_flux": int(len(cat["flux"])),
        "flux_column": f"FLUX_KRON_TOTAL_{cr.BAND}",
        "image": args.image.name,
        "image_url": BASE + IMAGE_NAME,
        "image_sha256": file_sha256(args.image),
    }
    depth = cr.depth_mag(cat["flux"], cat["snr"])
    s = ar.run(
        Path("."),
        cat=cat,
        depth=depth,
        cat_inputs=inputs,
        out=OUT,
        field=FIELD,
        frame=FRAME,
        offset=(0.0, 0.0),
    )
    src_x, src_y = (
        np.asarray(v, float) for v in w.world_to_pixel(SkyCoord(t["RA"], t["DEC"], unit="deg"))
    )
    rng = np.random.default_rng(SEED)
    ctrl = Table.read(OUT / f"control_{FIELD}.ecsv")
    pred = Table.read(OUT / f"predicted_{FIELD}.ecsv")
    apply_local(ctrl, d, w, src_x, src_y, rng)
    if len(pred):
        apply_local(pred, d, w, src_x, src_y, rng)
    c2 = ar.gate_c2(list(ctrl["class_local"]), list(ctrl["arc_ref"]))
    flags = [
        f"{r['system']}@{r['ra']:.6f},{r['dec']:.6f}"
        for r in pred
        if c2["passes"] and str(r["class_local"]) == "absent"
    ]
    ctrl.write(OUT / f"control_{FIELD}.ecsv", overwrite=True)
    pred.write(OUT / f"predicted_{FIELD}.ecsv", overwrite=True)
    summary = {
        "frozen": {
            "gate_C2": c2,
            "control_classes": _count(ctrl["class_local"]),
            "target_classes": _count(pred["class_local"]) if len(pred) else {},
            "flags": flags,
        },
        "without_local_noise_rule": {
            k: s[k] for k in ("gate_C2", "control_classes", "target_classes", "flags")
        },
        "inputs": s["inputs"],
        "rules": {
            **s["rules"],
            "local_noise": {
                "aperture_radius_arcsec": AP_ARCSEC,
                "n_random": N_RAND,
                "annulus_arcsec": list(R_RAND),
                "source_exclusion_arcsec": SRC_EXCL,
                "seed": SEED,
                "min_pred_over_sigma": MIN_SIGMA,
                "absent_pred_over_sigma": ABSENT_SIGMA,
            },
        },
        "posthoc_absent_controls": s["posthoc_absent_controls"],
    }
    summary = _clean(summary)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
