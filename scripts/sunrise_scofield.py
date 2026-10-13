"""WHL0137 (Sunrise): image-position consistency of the Scofield+2025 MrMARTIAN M5-G1 model.

The published maps (Zenodo 10.5281/zenodo.15110933, CC-BY-4.0) are 11792² float64 at 0.02",
1.1 GB each, scaled to z = 6.2. This script reads only small pixel patches around each
catalogued image (Scofield+2025 Table 3) with HTTP byte ranges, so nothing large is stored.

Per image it derives (provenance: model_prediction from the published maps; image positions are
observed as published):
* the deflection, and the source position beta = theta - alpha;
* kappa and mu from finite differences of the deflection, checked against the published kappa
  and magnification maps at the same pixel (validation of units, sign and scaling);
* per system, the source-plane offset of each image from the mean of the model's own constraint
  images, and its image-plane equivalent A^-1 delta-beta (linear; approximate near critical curves).

Images used as constraints by M5-G1 (1.1c,d; 1.2; 1.7a,c,d) are the known case: their
source-plane rms must be near the published 0.06". The others (1.1a,b,e; 1.7b; system 2) are
held out from this model.

Usage: python scripts/sunrise_scofield.py [--out results/sunrise/scofield25_positions.json]
"""

from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests

RECORD = "https://zenodo.org/api/records/15110933/files/{}/content"
FILES = {  # name -> (Zenodo key, md5 from the record API, 2026-10-13)
    "defx": ("whl0137-08_mrmartian_M5-G1_defx_250325.fits", "20f8cfbe8529dfac995a9b2685cd8392"),
    "defy": ("whl0137-08_mrmartian_M5-G1_defy_250325.fits", "de9d9a92bccac79fa4579992f475a541"),
    "kappa": ("whl0137-08_mrmartian_M5-G1_kappa_250325.fits", "5ef8d2ae623fb698f8cb3edcd31f14dc"),
    "mag": ("whl0137-08_mrmartian_M5-G1_mag_250325.fits", "c7c84e93a4f47c72bc70f89565c02b51"),
}
HEADER_BYTES = 2880  # one FITS block; checked against the header read at run time
NAXIS = 11792
BYTES_PER_PIX = 8  # BITPIX -64, big-endian

# Scofield+2025 Table 3 (observed positions, JWST frame, deg)
IMAGES = {
    "1.1a": (24.346591, -8.464253),
    "1.1b": (24.347093, -8.464749),
    "1.1c": (24.349662, -8.466417),
    "1.1d": (24.347503, -8.465067),
    "1.1e": (24.346822, -8.464511),
    "1.2a": (24.346396, -8.464045),
    "1.2b": (24.347901, -8.465367),
    "1.2c": (24.349429, -8.466311),
    "1.7a": (24.346351, -8.463995),
    "1.7b": (24.348570, -8.465825),
    "1.7c": (24.349359, -8.466278),
    "1.7d": (24.348005, -8.465436),
    "2a": (24.345958, -8.457367),
    "2b": (24.346171, -8.458089),
    "2c": (24.346522, -8.459075),
}
# M5-G1 constraints (Scofield+2025 Table 4): 1.1(c, d); 1.2; 1.7(a, c, d)
CONSTRAINTS = {"1.1c", "1.1d", "1.2a", "1.2b", "1.2c", "1.7a", "1.7c", "1.7d"}
PUBLISHED_SOURCE_RMS = 0.06  # arcsec, Scofield+2025 Table 4 (M5-G1)
HALF = 3  # patch half-width in pixels (7x7 patch)


def _get(session: requests.Session, url: str, start: int, stop: int) -> bytes:
    """One byte range, retried with back-off on errors, 429/5xx or short reads."""
    for attempt in range(8):
        try:
            r = session.get(url, headers={"Range": f"bytes={start}-{stop - 1}"}, timeout=60)
            if r.status_code == 206 and len(r.content) == stop - start:
                return r.content
        except requests.RequestException:
            pass
        time.sleep(min(2.0**attempt, 30.0))
    raise RuntimeError(f"range {start}-{stop} of {url} failed")


def read_header(session: requests.Session, key: str):
    from astropy.io import fits
    from astropy.wcs import WCS

    raw = _get(session, RECORD.format(key), 0, HEADER_BYTES)
    hdr = fits.Header.fromstring(raw.decode("ascii"))
    if (hdr["NAXIS1"], hdr["NAXIS2"], hdr["BITPIX"]) != (NAXIS, NAXIS, -64):
        raise RuntimeError(f"{key}: unexpected layout")
    return hdr, WCS(hdr)


def read_patch(session, key: str, i0: int, j0: int, half: int = HALF) -> np.ndarray:
    """Pixels [j0-half, j0+half] x [i0-half, i0+half] (0-based) of one map.

    Patches are cached as .npy under $SUNRISE_SCOFIELD_CACHE when set (scratch, not the repo)."""
    cache = os.environ.get("SUNRISE_SCOFIELD_CACHE")
    path = Path(cache) / f"{key}_{i0}_{j0}_{half}.npy" if cache else None
    if path is not None and path.exists():
        return np.load(path)
    url = RECORD.format(key)
    rows = []
    for j in range(j0 - half, j0 + half + 1):
        start = HEADER_BYTES + (j * NAXIS + i0 - half) * BYTES_PER_PIX
        raw = _get(session, url, start, start + (2 * half + 1) * BYTES_PER_PIX)
        rows.append(np.frombuffer(raw, ">f8"))
    patch = np.array(rows, float)
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.save(path, patch)
    return patch


def measure(session, wcs, name: str, ra: float, dec: float) -> dict:
    i, j = (float(v) for v in wcs.world_to_pixel_values(ra, dec))
    i0, j0 = int(round(i)), int(round(j))
    p = {k: read_patch(session, FILES[k][0], i0, j0) for k in FILES}
    step = abs(wcs.wcs.cdelt[1]) * 3600.0
    di, dj = i - i0, j - j0

    def at(a):  # bilinear at the image position
        c, r = HALF + di, HALF + dj
        c0, r0 = int(np.floor(c)), int(np.floor(r))
        fc, fr = c - c0, r - r0
        return float(
            a[r0, c0] * (1 - fc) * (1 - fr)
            + a[r0, c0 + 1] * fc * (1 - fr)
            + a[r0 + 1, c0] * (1 - fc) * fr
            + a[r0 + 1, c0 + 1] * fc * fr
        )

    # +i runs West (= -RA) on this north-up, east-left grid; alpha_x along +i, alpha_y along +j
    dax_dj, dax_di = np.gradient(p["defx"], step)
    day_dj, day_di = np.gradient(p["defy"], step)
    pxx, pyy, pxy = at(dax_di), at(day_dj), 0.5 * (at(dax_dj) + at(day_di))
    kappa = 0.5 * (pxx + pyy)
    det = (1 - pxx) * (1 - pyy) - pxy**2
    x = (ra - IMAGES["1.1e"][0]) * np.cos(np.deg2rad(dec)) * -3600.0  # West, arcsec
    y = (dec - IMAGES["1.1e"][1]) * 3600.0
    ax, ay = at(p["defx"]), at(p["defy"])
    return {
        "image": name,
        "x_arcsec": x,
        "y_arcsec": y,
        "alpha_x": ax,
        "alpha_y": ay,
        "beta_x": x - ax,
        "beta_y": y - ay,
        "psi": [pxx, pxy, pyy],
        "kappa_fd": kappa,
        "kappa_map": at(p["kappa"]),
        "mu_fd": 1.0 / det,
        "mu_map": at(p["mag"]),
        "constraint": name in CONSTRAINTS,
    }


def system_of(name: str) -> str:
    return name[:-1]


def summarise(rows: list[dict]) -> dict:
    out = {}
    for sys_ in sorted({system_of(r["image"]) for r in rows}):
        members = [r for r in rows if system_of(r["image"]) == sys_]
        ref = [r for r in members if r["constraint"]] or members
        bx = np.mean([r["beta_x"] for r in ref])
        by = np.mean([r["beta_y"] for r in ref])
        per = {}
        for r in members:
            d = np.array([r["beta_x"] - bx, r["beta_y"] - by])
            pxx, pxy, pyy = r["psi"]
            a = np.array([[1 - pxx, -pxy], [-pxy, 1 - pyy]])
            dtheta = np.linalg.solve(a, d)
            per[r["image"]] = {
                "constraint": r["constraint"],
                "dbeta_arcsec": float(np.hypot(*d)),
                "dtheta_linear_arcsec": float(np.hypot(*dtheta)),
                "mu_fd": r["mu_fd"],
            }
        out[sys_] = {
            "reference": "constraint images"
            if any(m["constraint"] for m in members)
            else "all images (held-out system)",
            "images": per,
        }
    cons = [r for r in rows if r["constraint"]]
    resid = []
    for r in cons:
        s = system_of(r["image"])
        resid.append(out[s]["images"][r["image"]]["dbeta_arcsec"])
    out["constraint_source_rms_arcsec"] = float(np.sqrt(np.mean(np.square(resid))))
    out["published_source_rms_arcsec"] = PUBLISHED_SOURCE_RMS
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="results/sunrise/scofield25_positions.json")
    args = ap.parse_args(argv)
    session = requests.Session()
    hdr, wcs = read_header(session, FILES["defx"][0])
    for k in ("defy", "kappa", "mag"):
        h2, _ = read_header(session, FILES[k][0])
        if (h2["CRVAL1"], h2["CRVAL2"], h2["CDELT1"]) != (
            hdr["CRVAL1"],
            hdr["CRVAL2"],
            hdr["CDELT1"],
        ):
            raise RuntimeError(f"{k}: grid differs from defx")
    with ThreadPoolExecutor(6) as pool:
        rows = list(
            pool.map(lambda kv: measure(requests.Session(), wcs, kv[0], *kv[1]), IMAGES.items())
        )
    result = {
        "provenance": "model_prediction (Scofield+2025 M5-G1 maps, z = 6.2) at observed positions",
        "source": "https://doi.org/10.5281/zenodo.15110933",
        "files_md5": {k: v[1] for k, v in FILES.items()},
        "access": "HTTP byte ranges, 7x7-pixel patches; md5 not verifiable (no full download)",
        "rows": rows,
        "summary": summarise(rows),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1))
    s = result["summary"]
    print(f"constraint source rms {s['constraint_source_rms_arcsec']:.3f} (published 0.06)")
    for r in rows:
        print(
            f"{r['image']:5s} kappa fd/map {r['kappa_fd']:.4f}/{r['kappa_map']:.4f} "
            f"mu fd/map {r['mu_fd']:9.2f}/{r['mu_map']:9.2f}"
        )
    for v in s.values():
        if isinstance(v, dict):
            for im, d in v["images"].items():
                print(
                    f"{im:5s} {'C' if d['constraint'] else 'H'} dbeta {d['dbeta_arcsec']:.3f} "
                    f"dtheta~{d['dtheta_linear_arcsec']:.2f} mu {d['mu_fd']:.1f}"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
