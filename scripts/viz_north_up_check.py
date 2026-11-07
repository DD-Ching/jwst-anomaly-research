"""Real-data check of ``viz.north_up`` on SMACS 2736 o001 F200W (pipeline catalog neighbours).

Picks catalogued 20-23 mag sources whose nearest catalogued neighbour (< 24 mag) is 0.8-2.2"
away, cuts 5" cutouts from the S3 i2d, resamples them north-up and finds the brightest pixel
within +-5 px of the neighbour position predicted from its east/north sky offset (a peak on the
window edge does not count; ``PARAMS``). The control
(``stored_ratio``) is the same window in the stored, unrotated cutout. Also renders contact
sheets (stored; north-up with one shared scale) for looking.

Usage: ``python -I scripts/viz_north_up_check.py OUT_DIR``. Result: results/viz_north_up/RESULT.md.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import fsspec
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from jwst_anomaly import cutouts, viz

BASE = "s3://stpubdata/jwst/public/jw02736/L3/t/o001/jw02736-o001_t001_nircam_clear-f200w_"
# ASSUMPTIONs (selection and matching of the check; recorded in the output):
PARAMS = {
    "target_mag": (20.0, 23.0),  # AB, aper_total_abmag; bright enough for a clear peak
    "neighbour_mag_max": 24.0,
    "sep_arcsec": (0.8, 2.2),  # inside a 5" cutout, clear of the target core
    "n_max_pairs": 12,  # unordered pairs: a mutual pair is tested once
    "half_px": 5,  # search half-width around the predicted neighbour position
    "match_px": 4.0,  # a non-edge peak within this distance counts as a match
}


def main(out_dir: Path) -> list[dict]:
    with fsspec.open(BASE + "cat.ecsv", "rb", anon=True) as f:
        cat = Table.read(f.read().decode(), format="ascii.ecsv")
    sky = cat["sky_centroid"]
    mag = np.asarray(cat["aper_total_abmag"], float)
    idx_n, d2, _ = sky.match_to_catalog_sky(sky, nthneighbor=2)
    sep = d2.arcsec
    lo, hi = PARAMS["target_mag"]
    smin, smax = PARAMS["sep_arcsec"]
    pick = np.isfinite(mag) & (mag > lo) & (mag < hi) & (sep > smin) & (sep < smax)
    pick &= mag[idx_n] < PARAMS["neighbour_mag_max"]
    seen: set[frozenset[int]] = set()
    cands = []
    for i in np.where(pick)[0]:
        pair = frozenset((int(i), int(idx_n[i])))
        if pair not in seen:
            seen.add(pair)
            cands.append(int(i))
        if len(cands) == PARAMS["n_max_pairs"]:
            break
    tgt = Table(
        {"source_uid": [f"c{i}" for i in cands], "ra": sky.ra.deg[cands], "dec": sky.dec.deg[cands]}
    )
    out = cutouts.make_cutouts(BASE + "i2d.fits", tgt, size_arcsec=5.0, out_dir=out_dir)
    half = PARAMS["half_px"]
    res = []
    for row, i in zip(out, cands, strict=True):
        if not row["path"] or row["quality_flag"] != "ok":
            continue
        data, head = fits.getdata(row["path"], extname="SCI", header=True)
        nu = viz.north_up(data, head)

        s = float(np.sqrt(np.prod(proj_plane_pixel_scales(WCS(head).celestial)))) * 3600.0
        east, north = sky[i].spherical_offsets_to(sky[idx_n[i]])
        ny, nx = nu.shape
        pr = (ny - 1) / 2 + north.to_value("arcsec") / s
        pc = (nx - 1) / 2 - east.to_value("arcsec") / s
        r0, c0 = int(round(pr)), int(round(pc))
        win = nu[r0 - half : r0 + half + 1, c0 - half : c0 + half + 1]
        dr, dc = np.unravel_index(np.nanargmax(win), win.shape)
        win0 = data[r0 - half : r0 + half + 1, c0 - half : c0 + half + 1]
        on_edge = int(dr) in (0, win.shape[0] - 1) or int(dc) in (0, win.shape[1] - 1)
        offset = float(np.hypot(r0 - half + dr - pr, c0 - half + dc - pc))
        res.append(
            {
                "uid": str(row["source_uid"]),
                "sep_arcsec": round(float(sep[i]), 2),
                "pred_rc": [round(pr, 1), round(pc, 1)],
                "peak_rc": [r0 - half + int(dr), c0 - half + int(dc)],
                "offset_px": round(offset, 2),
                "peak_on_window_edge": bool(on_edge),
                "match": bool(not on_edge and offset <= PARAMS["match_px"]),
                "stored_ratio": round(float(np.nanmax(win0) / np.nanmax(win)), 3),
                "pc1_1": float(head.get("PC1_1", np.nan)),
                "pc2_1": float(head.get("PC2_1", np.nan)),
                "rotation_deg": round(
                    float(np.degrees(np.arctan2(head.get("PC2_1", 0.0), head.get("PC1_1", 1.0)))), 2
                ),
            }
        )
    shown = out[out["path"] != ""]
    viz.contact_sheet(shown, out_dir / "sheet_stored.png")
    viz.contact_sheet(shown, out_dir / "sheet_northup.png", north_up=True, scale="all")
    return res


if __name__ == "__main__":
    rows = main(Path(sys.argv[1]))
    summary = {
        "provenance": "derived (observed pixels and pipeline catalogue positions)",
        "source": "jw02736-o001_t001_nircam_clear-f200w i2d + cat.ecsv (s3://stpubdata)",
        "params_assumption": PARAMS,
        "n_pairs": len(rows),
        "n_match": sum(r["match"] for r in rows),
        "offset_px_median": float(np.median([r["offset_px"] for r in rows])),
        "offset_px_max_matched": max((r["offset_px"] for r in rows if r["match"]), default=None),
        "stored_ratio_max": max(r["stored_ratio"] for r in rows),
        "rows": rows,
    }
    print(json.dumps(summary, indent=1))
