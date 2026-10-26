"""Post hoc checks of the MACS J1931 two-epoch grid hits (not in the frozen plan).

1. Mosaic level: for every hit and LW band, the epoch-2 ``_i2d`` peak fraction (brightest pixel
   over the background-subtracted 5 × 5 pixel sum; a PSF-sampled point source sits near the
   reference sources' value, a single hot pixel near 1) and the number of contributing exposures
   (bits set in ``CON``). Reference: the same peak fraction for every catalogued grid position in
   the band's brightest 20 % of epoch-2 fluxes, hit positions excluded.
2. Detector level: each hit mapped into the twelve 12299 o001 NRCBLONG ``_cal`` exposures (FITS
   SIP WCS from the header, pixels read by S3 byte range): detector pixel, DQ and peak fraction.
   The same detector pixel in all three filters marks a detector-fixed artefact.
3. Known case: distance of each SN Eos image to its nearest grid row, and that row's significances.

All numbers are ``derived``.

    python scripts/m1931_grid_posthoc.py --data data/m1931 \\
        --out results/transients/macs1931_grid/posthoc.json
"""

from __future__ import annotations

import argparse
import json
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

E1 = "jw06882-o038_t041_nircam_clear-"
E2 = "jw12299-o001_t001_nircam_clear-"
BANDS = ("f277w", "f356w", "f444w")
CAL = [
    f"jw12299001001_031{v}_0000{d}_nrcblong_cal.fits"
    for v in ("01", "03", "05")
    for d in range(1, 5)
]
SN_EOS = Path("results/transients/sneos/summary.json")


def peak_fraction(stamp: np.ndarray, background: float) -> float | None:
    """Brightest pixel over the background-subtracted sum of a 5 × 5 stamp; None if undefined."""
    if stamp.shape != (5, 5) or not np.isfinite(stamp).all():
        return None
    c = stamp.astype(float) - background
    total = float(c.sum())
    return round(float(c.max()) / total, 3) if total > 0 else None


def stamp_at(img, x: int, y: int, half: int = 2) -> np.ndarray | None:
    """(2 half + 1)² pixels centred on (x, y), or None if any part is off the image."""
    ny, nx = img.shape[-2:]
    if x - half < 0 or y - half < 0 or x + half >= nx or y + half >= ny:
        return None
    return np.asarray(img[y - half : y + half + 1, x - half : x + half + 1], dtype=float)


def n_exposures(con: np.ndarray, x: int, y: int) -> int | None:
    """Exposures contributing at (x, y): bits set across all context planes; None off the image."""
    ny, nx = con.shape[-2:]
    if not (0 <= x < nx and 0 <= y < ny):
        return None
    planes = con[:, y, x] if con.ndim == 3 else np.atleast_1d(con[y, x])
    return int(sum(bin(int(v) & 0xFFFFFFFF).count("1") for v in planes))


def _pix(wcs: WCS, ra: float, dec: float) -> tuple[int, int]:
    x, y = wcs.world_to_pixel_values(ra, dec)
    return int(round(float(x))), int(round(float(y)))


def _local_pf(img, x: int, y: int) -> float | None:
    st, ring = stamp_at(img, x, y), stamp_at(img, x, y, half=10)
    if st is None or ring is None:
        return None
    return peak_fraction(st, float(np.nanmedian(ring)))


def mosaic_level(data: Path, hits: Table, grid: Table) -> dict:
    hit_ids = set(str(p) for p in hits["pos_id"])
    out = {}
    for b in BANDS:
        with fits.open(data / f"{E2}{b}_i2d.fits", memmap=True) as h:
            sci, con, wcs = h["SCI"].data, h["CON"].data, WCS(h["SCI"].header)
            rows = []
            for r in hits:
                x, y = _pix(wcs, r["ra"], r["dec"])
                rows.append(
                    {
                        "group": int(r["group"]),
                        "pos_id": str(r["pos_id"]),
                        "peak_fraction": _local_pf(sci, x, y),
                        "n_exposures_e2": n_exposures(con, x, y),
                    }
                )
            f2 = np.asarray(grid[f"{b.upper()}_flux2"], dtype=float)
            ok = np.isfinite(f2) & (f2 > np.nanpercentile(f2, 80))
            ok &= np.array([str(p) not in hit_ids for p in grid["pos_id"]])
            ref = [_local_pf(sci, *_pix(wcs, r["ra"], r["dec"])) for r in grid[ok]]
            ref = np.array([v for v in ref if v is not None])
        out[b.upper()] = {
            "hits": rows,
            "reference": (
                {
                    "n": int(ref.size),
                    "p50": round(float(np.median(ref)), 3),
                    "p99": round(float(np.percentile(ref, 99)), 3),
                    "max": round(float(ref.max()), 3),
                }
                if ref.size
                else {"n": 0}
            ),
        }
    return out


def detector_level(hits: Table) -> list[dict]:
    from jwst_anomaly.query import resolve_cloud_uris

    uris = resolve_cloud_uris([f"mast:JWST/product/{n}" for n in CAL])

    def one(item: tuple[str, str]) -> dict:
        name, uri = item
        kw = {"use_fsspec": True, "fsspec_kwargs": {"anon": True}, "lazy_load_hdus": True}
        with fits.open(uri, **kw) as h:
            sci, dq, wcs = h["SCI"], h["DQ"], WCS(h["SCI"].header)
            rows = []
            for r in hits:
                x, y = _pix(wcs, r["ra"], r["dec"])
                st, ring = stamp_at(sci.section, x, y), stamp_at(sci.section, x, y, half=10)
                if st is None or ring is None:
                    continue
                bg = float(np.nanmedian(ring))
                rows.append(
                    {
                        "group": int(r["group"]),
                        "pos_id": str(r["pos_id"]),
                        "x": x,
                        "y": y,
                        "dq": int(dq.section[y, x]),
                        "peak": round(float(np.nanmax(st) - bg), 3),
                        "peak_fraction": peak_fraction(st, bg),
                    }
                )
            return {"file": name, "filter": h[0].header["FILTER"], "rows": rows}

    with ThreadPoolExecutor(len(CAL)) as ex:
        return list(ex.map(one, zip(CAL, uris, strict=True)))


def known_case(grid: Table) -> dict:
    images = json.loads(SN_EOS.read_text())["images"]
    out = {}
    for im in images:
        if not im["catalogued"]:
            continue
        ra, dec = im["ra"], im["dec"]
        d = np.hypot((grid["ra"] - ra) * math.cos(math.radians(dec)), grid["dec"] - dec) * 3600
        r = grid[int(np.argmin(d))]
        out[im["catalogued"]] = {
            "nearest_pos_id": str(r["pos_id"]),
            "distance_arcsec": round(float(np.min(d)), 3),
            "hit": bool(r["hit"]),
            **{f"{b}_sigma": _num(r[f"{b}_sigma"]) for b in ("F277W", "F356W", "F444W")},
            **{f"{b}_dmag": _num(r[f"{b}_dmag"]) for b in ("F277W", "F356W", "F444W")},
        }
    return out


def _num(v) -> float | None:
    v = float(v)
    return round(v, 3) if math.isfinite(v) else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    hits = Table.read(a.data / "out" / "hits.ecsv")
    grid = Table.read(a.data / "out" / "grid.ecsv")
    sizes = {
        p.name: p.stat().st_size
        for e in (E1, E2)
        for b in BANDS
        for p in [a.data / f"{e}{b}_i2d.fits"]
    }
    out = {
        "label": "derived, post hoc (not in PLAN.md)",
        "mosaic_bytes": sizes,
        "grid_ecsv_bytes": (a.data / "out" / "grid.ecsv").stat().st_size,
        "known_case": known_case(grid),
        "mosaic_level": mosaic_level(a.data, hits, grid),
        "detector_level": detector_level(hits),
    }
    a.out.write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
