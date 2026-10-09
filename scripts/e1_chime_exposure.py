"""E1: reduce the CHIME/FRB Catalog 2 exposure file to a small tracked declination profile.

The file (doi:10.11570/25.0066, ``data/exposure/chimefrbcat2_exposure.h5``, 216 MB) holds two
HEALPix (nside 4096, RING, ICRS) maps of *time-integrated* exposure (s) over 2018-09-04 to
2023-09-15, for the upper and lower transits. It has no time axis, so it cannot give a per-day
uptime series. Reason for the download (review of PR #112): it was the only published CHIME
exposure model that might separate day-level uptime from event dependence. This script records what
it contains and keeps a < 1 MB reduction; the raw file is deleted afterwards by the caller.

Needs ``h5py`` (not a project dependency; run it from a scratch environment):

  python scripts/e1_chime_exposure.py <path to chimefrbcat2_exposure.h5>

Writes ``results/e1_events/chime_exposure_dec_profile.ecsv``: per 0.1 deg declination bin, the
pixel-mean exposure and its RA scatter for each transit.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "e1_events" / "chime_exposure_dec_profile.ecsv"
DEC_BIN = 0.1  # deg (ASSUMPTION: fine enough for a transit instrument's Dec dependence)


def ring_of_pixel(pix: np.ndarray, nside: int) -> np.ndarray:
    """HEALPix RING ordering: ring index (1 .. 4 nside - 1) of each pixel."""
    npix = 12 * nside * nside
    ncap = 2 * nside * (nside - 1)
    pix = np.asarray(pix, dtype=np.int64)
    ring = np.empty_like(pix)
    north = pix < ncap
    ring[north] = (1 + np.floor(np.sqrt(1 + 2 * pix[north]))).astype(np.int64) // 2
    south = pix >= npix - ncap
    ps = npix - pix[south]
    ring[south] = 4 * nside - (1 + np.floor(np.sqrt(2 * ps - 1))).astype(np.int64) // 2
    eq = ~(north | south)
    ring[eq] = (pix[eq] - ncap) // (4 * nside) + nside
    return ring


def ring_dec(ring: np.ndarray, nside: int) -> np.ndarray:
    """Declination (deg) of a RING-ordering ring index."""
    ring = np.asarray(ring, dtype=float)
    n = float(nside)
    z = np.where(
        ring < n,
        1 - ring**2 / (3 * n * n),
        np.where(
            ring <= 3 * n, 4 / 3 - 2 * ring / (3 * n), -(1 - (4 * n - ring) ** 2 / (3 * n * n))
        ),
    )
    return np.degrees(np.arcsin(z))


def reduce(path: Path, chunk: int = 1 << 22) -> Table:
    import h5py

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    with h5py.File(path, "r") as f:
        nside = int(f["upper"].attrs["nside"])
        nring = 4 * nside - 1
        sums = {k: np.zeros(nring + 1) for k in ("upper", "lower")}
        sq = {k: np.zeros(nring + 1) for k in ("upper", "lower")}
        cnt = np.zeros(nring + 1)
        npix = f["upper"].shape[0]
        for start in range(0, npix, chunk):
            stop = min(start + chunk, npix)
            r = ring_of_pixel(np.arange(start, stop), nside)
            cnt += np.bincount(r, minlength=nring + 1)
            for k in ("upper", "lower"):
                v = f[k][start:stop]
                sums[k] += np.bincount(r, weights=v, minlength=nring + 1)
                sq[k] += np.bincount(r, weights=v * v, minlength=nring + 1)
        meta_attrs = {k: str(v) for k, v in f.attrs.items()}
    rings = np.arange(1, nring + 1)
    dec = ring_dec(rings, nside)
    edges = np.arange(-90, 90 + DEC_BIN / 2, DEC_BIN)
    b = np.clip(np.digitize(dec, edges) - 1, 0, len(edges) - 2)
    out = {"dec_lo": np.round(edges[:-1], 1), "dec_hi": np.round(edges[1:], 1)}
    n = np.bincount(b, weights=cnt[1:], minlength=len(edges) - 1)
    for k in ("upper", "lower"):
        s = np.bincount(b, weights=sums[k][1:], minlength=len(edges) - 1)
        s2 = np.bincount(b, weights=sq[k][1:], minlength=len(edges) - 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = s / n
            out[f"{k}_mean_s"] = np.round(mean, 1)
            out[f"{k}_ra_rms_s"] = np.round(np.sqrt(np.maximum(s2 / n - mean**2, 0)), 1)
    out["n_pix"] = n.astype(np.int64)
    t = Table(out)
    t = t[t["n_pix"] > 0]
    t.meta = {
        "provenance": "derived",
        "source": f"CHIME/FRB Catalog 2 exposure (doi:10.11570/25.0066) {path.name}, sha256 "
        f"{h.hexdigest()}, nside {nside} RING; reduced by scripts/e1_chime_exposure.py",
        "file_attrs": meta_attrs,
        "note": "time-integrated exposure only (no time axis): it cannot model per-day uptime",
    }
    return t


if __name__ == "__main__":
    t = reduce(Path(sys.argv[1]))
    t.write(OUT, format="ascii.ecsv", overwrite=True)
    print(t.meta["source"])
    print(len(t), "rows; total exposure-weighted check:", float(np.nansum(t["upper_mean_s"])))
