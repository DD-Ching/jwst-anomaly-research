"""Cut a level-3 ``_i2d.fits`` to the sky footprint of another mosaic (SCI and ERR only).

``transient_grid.py`` memory-maps whole local mosaics. Some survey mosaics are far larger than
the region two epochs share (UNCOVER ``jw02561-o001_t003`` SW bands: 5.2-7.7 GB, against a
second-epoch footprint of a few arcmin), so streaming them whole wastes disk and transfer.
This script reads only the rows of SCI and ERR that cover the reference mosaic's footprint (one
byte range per extension when remote, via ``jwst_anomaly.cutouts``), slices the columns, and
writes a local file with the primary header, a WCS shifted to the cut and the provenance
keywords ``CUTFROM``, ``CUTX0``, ``CUTY0``. Pixels are copied unchanged (``observed``).
Full-width rows are read (one byte range per extension), so memory peaks at one row strip of
the source mosaic.

    python scripts/i2d_region.py --ref e2_f150w_i2d.fits \\
        --src s3://stpubdata/jwst/public/jw02561/L3/t/o001/<UNCOVER o001 f150w>_i2d.fits \\
        --out data/transients/a2744/jw02561-o001_t003_nircam_clear-f150w_i2d.fits
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from jwst_anomaly import cutouts

PAD_PX = 50  # ASSUMPTION: margin around the footprint so recentring apertures stay inside


def footprint_corners(ref: str) -> np.ndarray:
    """Sky (RA, Dec) of the four corners of the reference mosaic's SCI array, degrees."""
    uri = cutouts._resolve_uri(ref)
    with cutouts._open_image(uri, None, 2**16) as (h, _fileobj):
        hdr = h["SCI"].header
        wcs = cutouts._read_wcs(hdr, uri)
        return np.asarray(wcs.calc_footprint(axes=(int(hdr["NAXIS1"]), int(hdr["NAXIS2"]))))


def pixel_box(
    wcs: WCS, corners: np.ndarray, shape: tuple[int, int], pad: int
) -> tuple[int, int, int, int]:
    """Pixel box (x0, x1, y0, y1; end-exclusive) of the sky corners, padded and clipped."""
    xy = np.asarray(wcs.all_world2pix(corners, 0))
    ny, nx = shape
    x0 = max(int(np.floor(xy[:, 0].min())) - pad, 0)
    x1 = min(int(np.ceil(xy[:, 0].max())) + pad + 1, nx)
    y0 = max(int(np.floor(xy[:, 1].min())) - pad, 0)
    y1 = min(int(np.ceil(xy[:, 1].max())) + pad + 1, ny)
    if x0 >= x1 or y0 >= y1:
        raise ValueError("the reference footprint does not overlap the source mosaic")
    return x0, x1, y0, y1


def cut(src: str, corners: np.ndarray, out: Path, pad: int = PAD_PX) -> tuple[int, ...]:
    """Write SCI and ERR of ``src`` over the corners' pixel box to ``out``; return the box."""
    uri = cutouts._resolve_uri(src)
    with cutouts._open_image(uri, None, 2**20) as (hdul, _fileobj):
        sci = hdul["SCI"]
        wcs = cutouts._read_wcs(sci.header, uri)
        shape = (int(sci.header["NAXIS2"]), int(sci.header["NAXIS1"]))
        x0, x1, y0, y1 = pixel_box(wcs, corners, shape, pad)
        prim = fits.PrimaryHDU(header=hdul[0].header.copy())
        prim.header["CUTFROM"] = (cutouts._image_stem(uri), "mosaic this region was cut from")
        prim.header["CUTX0"] = (x0, "0-based column of pixel (0, 0) in CUTFROM")
        prim.header["CUTY0"] = (y0, "0-based row of pixel (0, 0) in CUTFROM")
        hdus = [prim]
        for name in ("SCI", "ERR"):
            hdu = cutouts._get_hdu(hdul, name)
            if hdu is None:
                raise ValueError(f"{uri}: no {name} extension")
            sec = cutouts._image_section(hdu, shape, name, uri)
            data = np.asarray(sec[y0:y1, x0:x1], dtype=np.float32)
            hdr = hdu.header.copy()
            # keys that describe the full mosaic: S_REGION is its footprint; CD would override PC
            for k in ("S_REGION", "CD1_1", "CD1_2", "CD2_1", "CD2_2"):
                hdr.pop(k, None)
            hdr.update(wcs[y0:y1, x0:x1].to_header(relax=True))
            for k in ("NAXIS1", "NAXIS2"):
                hdr.pop(k, None)
            hdus.append(fits.ImageHDU(data, header=hdr, name=name))
    out.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(hdus).writeto(out, overwrite=True)
    return x0, x1, y0, y1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ref", required=True, help="mosaic whose footprint defines the region")
    ap.add_argument("--src", required=True, help="mosaic to cut (local path, s3://, mast:)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--pad-px", type=int, default=PAD_PX)
    a = ap.parse_args(argv)
    t = time.time()
    box = cut(a.src, footprint_corners(a.ref), a.out, a.pad_px)
    mb = a.out.stat().st_size / 1e6
    print(
        f"{a.out.name}: box x {box[0]}:{box[1]} y {box[2]}:{box[3]}, {mb:.0f} MB, "
        f"{time.time() - t:.0f} s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
