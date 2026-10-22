"""MACS0417 CANUCS reference run: post hoc checks of the family-10 flags (not frozen rules).

For each family-10 prediction: 0.3" forced aperture sum on the CANUCS DR1 F200W BCG-subtracted
40 mas image, the local aperture noise (sigma-clipped std of 200 random apertures 1.5-4" away)
and the three nearest CANUCS DR1 sources with their photo-z. Output ``posthoc.json``.

    python scripts/m0417_canucsref_posthoc.py --cat <photometry-cat> --image <f200w bgsub-sci>
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from jwst_anomaly.provenance import file_sha256

OUT = Path("results/lens_consistency/macs0417_canucsref")
AP_ARCSEC, PIX = 0.3, 0.04
N_RAND, R_RAND = 200, (1.5, 4.0)
SRC_EXCL = 0.5  # random apertures centred within this of a CANUCS source are redrawn (arcsec)
URL = (
    "https://archive.stsci.edu/hlsps/canucs/dr1/macs0417/clu/"
    "hlsp_canucs_jwst_nircam_macs0417-clu-40mas_f200w_v1_bgsub-sci.fits.gz"
)


def aperture_sum(d: np.ndarray, x: float, y: float, r: float) -> float:
    """Sum inside radius ``r`` (pixels) at (x, y); NaN if the aperture leaves the image or holds
    any NaN pixel (no partial coverage)."""
    h = int(np.ceil(r)) + 1
    y0, x0 = int(round(y)), int(round(x))
    if y0 - h < 0 or x0 - h < 0 or y0 + h + 1 > d.shape[0] or x0 + h + 1 > d.shape[1]:
        return float("nan")
    yy, xx = np.mgrid[y0 - h : y0 + h + 1, x0 - h : x0 + h + 1]
    m = (xx - x) ** 2 + (yy - y) ** 2 <= r * r
    v = d[y0 - h : y0 + h + 1, x0 - h : x0 + h + 1][m]
    return float(v.sum()) if np.isfinite(v).all() else float("nan")


def _clean(v):
    """JSON-safe: NaN / inf -> None (``allow_nan=False`` below)."""
    if isinstance(v, float) and not np.isfinite(v):
        return None
    if isinstance(v, dict):
        return {k: _clean(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_clean(x) for x in v]
    return v


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--image", type=Path, required=True)
    args = ap.parse_args(argv)
    with fits.open(args.image) as hl:
        if hl[0].header.get("BUNIT") != "nJy":
            raise ValueError(f"unexpected BUNIT {hl[0].header.get('BUNIT')!r}")
        w, d = WCS(hl[0].header), np.asarray(hl[0].data, float)
    if args.image.name not in URL:
        raise ValueError("--image is not the F200W 40 mas BCG-subtracted mosaic named in URL")
    if abs(proj_plane_pixel_scales(w.celestial).mean() * 3600.0 - PIX) > 1e-4:
        raise ValueError("image pixel scale is not 0.04 arcsec")
    cat = Table.read(args.cat)
    cs = SkyCoord(cat["RA"], cat["DEC"], unit="deg")
    pred = Table.read(OUT / "predicted_macs0417.ecsv")
    rng = np.random.default_rng(20261022)
    cpx, cpy = (np.asarray(v, float) for v in w.world_to_pixel(cs))
    pairs = pred[pred["system"] == "10"]
    r = AP_ARCSEC / PIX
    rows = []
    for p in pairs:
        x, y = (float(v) for v in w.world_to_pixel_values(p["ra"], p["dec"]))
        rand, tries = [], 0
        while len(rand) < N_RAND and tries < 50 * N_RAND:
            tries += 1
            a, s = rng.uniform(0, 2 * np.pi), rng.uniform(*R_RAND) / PIX
            rx, ry = x + s * np.cos(a), y + s * np.sin(a)
            if np.hypot(cpx - rx, cpy - ry).min() * PIX < SRC_EXCL:
                continue
            v = aperture_sum(d, rx, ry, r)
            if np.isfinite(v):
                rand.append(v)
        if len(rand) < N_RAND:
            raise RuntimeError(f"only {len(rand)} of {N_RAND} random apertures accepted")
        _, med, std = sigma_clipped_stats(rand)
        sep = cs.separation(SkyCoord(p["ra"], p["dec"], unit="deg")).arcsec
        near = [
            {
                "source": int(cat["SOURCE"][k]),
                "sep_arcsec": float(sep[k]),
                "flux_f200w_njy": float(cat["FLUX_KRON_TOTAL_F200W"][k]),
                "flux_over_pred": float(cat["FLUX_KRON_TOTAL_F200W"][k] * 1e-9 / p["f_pred"]),
                "z_ml": float(cat["Z_ML"][k]),
                "z160": float(cat["Z160"][k]),
                "z840": float(cat["Z840"][k]),
                "flag_high_bcgcontami": bool(cat["FLAG_HIGH_BCGCONTAMI"][k]),
            }
            for k in np.argsort(sep)[:3]
        ]
        ap_sum = aperture_sum(d, x, y, r)
        rows.append(
            {
                "image_class": str(p["image_class"]),
                "class": str(p["class"]),
                "ra": float(p["ra"]),
                "dec": float(p["dec"]),
                "f_pred_njy": float(p["f_pred"] * 1e9),
                "ap03_njy": ap_sum,
                "local_median_njy": float(med),
                "local_std_njy": float(std),
                "pred_over_local_std": float(p["f_pred"] * 1e9 / std),
                "ap_minus_median_over_std": float((ap_sum - med) / std),
                "n_random_used": len(rand),
                "nearest_sources": near,
            }
        )
    ctrl = Table.read(OUT / "control_macs0417.ecsv")
    idx = {int(v): k for k, v in enumerate(cat["SOURCE"])}
    catalogued = [
        {
            "images": str(c["images"]),
            "counterpart": int(c["counterpart"]),
            "flux_f200w_njy": float(cat["FLUX_KRON_TOTAL_F200W"][idx[int(c["counterpart"])]]),
            **{
                k.lower(): float(cat[k][idx[int(c["counterpart"])]])
                for k in ("Z_ML", "Z160", "Z840")
            },
        }
        for c in ctrl[(ctrl["system"] == "10") & (ctrl["counterpart"] >= 0)]
    ]
    miss = [q for q in pairs if q["image_class"] == "missing"]
    m, m2 = (miss[0], miss[1]) if len(miss) == 2 else (None, None)
    res = {
        "label": "post hoc (not frozen rules); observed pixels, derived sums",
        "image": args.image.name,
        "image_url": URL,
        "image_sha256": file_sha256(args.image),
        "aperture_radius_arcsec": AP_ARCSEC,
        "random_apertures": {
            "n": N_RAND,
            "annulus_arcsec": list(R_RAND),
            "source_exclusion_arcsec": SRC_EXCL,
            "seed": 20261022,
        },
        "separation_between_missing_arcsec": (
            float(
                SkyCoord(m["ra"], m["dec"], unit="deg")
                .separation(SkyCoord(m2["ra"], m2["dec"], unit="deg"))
                .arcsec
            )
            if m is not None
            else None
        ),
        "note": "target apertures as they lie; random apertures avoid catalogue sources",
        "family10": rows,
        "family10_catalogued": catalogued,
    }
    res = _clean(res)
    (OUT / "posthoc.json").write_text(json.dumps(res, indent=1, allow_nan=False) + "\n")
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
