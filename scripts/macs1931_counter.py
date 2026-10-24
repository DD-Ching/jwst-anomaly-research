"""MACS J1931.8-2635: imaging footprint for the counter-image route (no model prediction read).

The MACS J1931 gate K passes on the VENUS AstroLensPy v2.1 maps (``../macs1931_gatek``); six spec-z
systems have more predicted than catalogued images. A counter-image test needs imaging that covers
the field. This step lists the VENUS (GO 6882) level-3 NIRCam mosaics around the catalogued images
(MAST CAOM metadata only), checks which catalogued images fall inside each mosaic's ``s_region``,
and records the F200W level-3 products (sizes, S3 URIs) the test would read by byte range. Only
catalogued (observed) positions are used; no map value and no predicted position is read.
Writes ``summary.json``.

    python scripts/macs1931_counter.py footprint
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
from astropy.table import Table
from matplotlib.path import Path as MplPath

OUT = Path("results/lens_consistency/macs1931_counter")
IMAGES = Path("results/lens_consistency/macs1931_inputs/m1931_images.ecsv")
PROGRAM = "6882"
BOX_DEG = 0.1  # ASSUMPTION: search box half-width around the mean catalogued position
PRODUCT_TYPES = ("I2D", "CAT", "SEGM")
_NUM = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def parse_polygons(s_region: str) -> list[np.ndarray]:
    """Vertices (N x 2, RA and Dec in degrees) of each ``POLYGON`` in an STC-S ``s_region``."""
    out = []
    for chunk in re.split(r"POLYGON", s_region, flags=re.IGNORECASE)[1:]:
        v = np.array([float(x) for x in _NUM.findall(chunk)])
        if v.size < 6 or v.size % 2:
            raise ValueError(f"malformed polygon in s_region: {chunk[:60]!r}")
        out.append(v.reshape(-1, 2))
    if not out:
        raise ValueError("no POLYGON in s_region")
    return out


def inside(ra, dec, polygons: list[np.ndarray]) -> np.ndarray:
    """True where (ra, dec) lies inside any polygon (planar test; mosaics span a few arcmin)."""
    pts = np.column_stack([np.asarray(ra, float), np.asarray(dec, float)])
    hit = np.zeros(len(pts), bool)
    for v in polygons:
        hit |= MplPath(v).contains_points(pts)
    return hit


def footprint(images: Table, observations: Table) -> list[dict]:
    """Per mosaic: band, obs_id and how many used / spec-z images its footprint contains."""
    used = images[np.asarray(images["used"], bool)]
    spec = used[np.asarray(used["z_spec"], bool)]
    rows = []
    for o in observations:
        polys = parse_polygons(str(o["s_region"]))
        rows.append(
            {
                "obs_id": str(o["obs_id"]),
                "filter": str(o["filters"]),
                "target_name": str(o["target_name"]),
                "t_min_mjd": float(o["t_min"]),
                "n_used_inside": int(inside(used["ra"], used["dec"], polys).sum()),
                "n_used": len(used),
                "n_spec_inside": int(inside(spec["ra"], spec["dec"], polys).sum()),
                "n_spec": len(spec),
            }
        )
    return sorted(rows, key=lambda r: (r["obs_id"], r["filter"]))


def i2d_header(uri: str) -> dict:
    """SCI-header facts of a remote mosaic, read by byte range (no pixel read)."""
    from jwst_anomaly import cutouts

    with cutouts._open_image(cutouts._resolve_uri(uri), None, 2**16) as (h, _fileobj):
        hdr = h["SCI"].header
        return {
            "uri": uri,
            "naxis": [int(hdr["NAXIS1"]), int(hdr["NAXIS2"])],
            "bunit": str(hdr.get("BUNIT", "")),
            "pixel_arcsec": float(np.sqrt(float(hdr["PIXAR_A2"]))),
        }


def run_footprint() -> dict:
    from astroquery.mast import Observations

    from jwst_anomaly import query

    images = Table.read(IMAGES)
    ra0, dec0 = float(np.mean(images["ra"])), float(np.mean(images["dec"]))
    box_ra = BOX_DEG / np.cos(np.radians(dec0))  # same on-sky half-width in RA
    obs = Observations.query_criteria(
        proposal_id=PROGRAM,
        instrument_name="NIRCAM/IMAGE",
        calib_level=3,
        s_ra=[ra0 - box_ra, ra0 + box_ra],
        s_dec=[dec0 - BOX_DEG, dec0 + BOX_DEG],
    )
    rows = footprint(images, obs)
    f200 = obs[obs["filters"] == "F200W"]
    products = Observations.get_product_list(f200)
    keep = products[
        np.isin(products["productSubGroupDescription"], PRODUCT_TYPES)
        & np.char.startswith(np.asarray(products["productFilename"], str), "jw06882-")
    ]
    uris = query.resolve_cloud_uris(list(keep["dataURI"])) if len(keep) else []
    prods = [
        {"file": str(p["productFilename"]), "bytes": int(p["size"]), "s3": u}
        for p, u in zip(keep, uris, strict=True)
    ]
    i2d = next((p["s3"] for p in prods if p["file"].endswith("_i2d.fits") and p["s3"]), "")
    summary = {
        "program": PROGRAM,
        "f200w_i2d_header": i2d_header(i2d) if i2d else {},
        "centre_deg": [ra0, dec0],
        "box_deg": BOX_DEG,
        "mosaics": rows,
        "f200w_products": sorted(prods, key=lambda r: r["file"]),
        "provenance": "observed (MAST CAOM metadata); positions: VENUS image list (observed)",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("footprint", help="VENUS mosaics around the catalogued images (MAST)")
    args = ap.parse_args(argv)
    if args.cmd == "footprint":
        s = run_footprint()
        for r in s["mosaics"]:
            print(r["obs_id"], r["filter"], r["n_used_inside"], "/", r["n_used"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
