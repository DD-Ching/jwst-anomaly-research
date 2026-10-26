"""MACS0647 1433 g85 (post-hoc vetting): per-exposure, same-visit F480M, mosaic-depth and aperture
evidence for the F200W brightening between 1433 o010 and o020.

Not part of the frozen plan (``results/transients/macs0647_1433/PLAN.md``). Reads the local
``_i2d`` / ``_cal`` files under ``$JWST_ANOMALY_DATA/m0647_1433`` (streamed from MAST, deleted
after use). The ``_cal`` positions use the FITS SIP approximation in their SCI headers. Every
number is derived from observed pixels; thresholds are ASSUMPTIONs.

    python scripts/macs0647_g85_vet.py
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

from jwst_anomaly import paths

RA, DEC = 101.99665, 70.250721  # arc_diff h1023 (results/transients/macs0647_1433/core/hits.ecsv)
APER_ARCSEC, ANNULUS_ARCSEC = 0.15, (0.35, 0.60)  # ASSUMPTION
OUT = Path("results/transients/macs0647_1433/g85_vet.json")
CROSSMATCH = {  # observed; astroquery calls of 2026-10-26T06:28:03Z, copied from the session
    "queried_utc": "2026-10-26T06:28:03Z",
    "calls": "Simbad.query_region / Ned.query_region (5 arcsec), Gaia.cone_search (3 arcsec)",
    "simbad_5arcsec": [],
    "gaia_dr3_3arcsec": 0,
    "ned_5arcsec": [
        {"name": "[CDK2017] 07064", "type": "G", "redshift": None, "sep_arcmin": 0.02},
        {"name": "[CDK2017] 07127", "type": "G", "redshift": None, "sep_arcmin": 0.082},
    ],
    "provenance": "observed (catalogue query)",
}
CAL = (
    "jw01433020001_02101_00001_nrca3_cal.fits",
    "jw01433020001_02101_00002_nrca4_cal.fits",
    "jw01433020001_02101_00003_nrca4_cal.fits",
    "jw01433020001_02101_00004_nrca3_cal.fits",
)


def robust_rms(a: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    return float(1.4826 * np.median(np.abs(a - np.median(a))))


def peak_snr(sci: np.ndarray, x: float, y: float, box: int = 3, bg: int = 40) -> dict:
    """Peak within ±box px of (x, y) over the robust rms of a ±bg px stamp (median removed)."""
    xi, yi = int(round(x)), int(round(y))
    ny, nx = sci.shape
    big = np.asarray(sci[max(yi - bg, 0) : yi + bg + 1, max(xi - bg, 0) : xi + bg + 1], float)
    med, rms = float(np.nanmedian(big)), robust_rms(big)
    small = np.asarray(sci[max(yi - box, 0) : yi + box + 1, max(xi - box, 0) : xi + box + 1], float)
    n_ok = int(np.isfinite(small).sum())
    if n_ok == 0 or not (0 <= xi < nx and 0 <= yi < ny):
        return {"x": x, "y": y, "n_finite_box": n_ok, "peak_snr": None, "peak_offset_px": None}
    k = np.unravel_index(np.nanargmax(small), small.shape)
    oy, ox = k[0] + max(yi - box, 0) - yi, k[1] + max(xi - box, 0) - xi
    edge = min(xi, yi, nx - 1 - xi, ny - 1 - yi)
    return {
        "x": x,
        "y": y,
        "edge_distance_px": int(edge),
        "n_finite_box": n_ok,
        "peak_snr": float((small[k] - med) / rms),
        "peak_offset_px": [int(ox), int(oy)],
    }


def aperture_njy(path: Path) -> dict:
    with fits.open(path, memmap=True) as h:
        w = WCS(h["SCI"].header)
        x, y = (float(v) for v in w.world_to_pixel_values(RA, DEC))
        pix = float(np.sqrt(h["SCI"].header["PIXAR_A2"]))
        n = int(ANNULUS_ARCSEC[1] / pix) + 2
        yy, xx = np.mgrid[int(y) - n : int(y) + n + 1, int(x) - n : int(x) + n + 1]
        s = np.asarray(h["SCI"].data[yy, xx], float)
        rr = np.hypot(xx - x, yy - y) * pix
        bg = float(np.nanmedian(s[(rr > ANNULUS_ARCSEC[0]) & (rr < ANNULUS_ARCSEC[1])]))
        flux = float(np.nansum(s[rr < APER_ARCSEC] - bg) * h["SCI"].header["PIXAR_SR"] * 1e15)
        con = np.asarray(h["CON"].data[:, int(round(y)), int(round(x))]).astype(np.int64)
        nexp = int(sum(bin(int(v)).count("1") for v in con))
        ps = peak_snr(h["SCI"].data, x, y, box=2, bg=25)
        pa = float(h["SCI"].header["PA_APER"])
    return {"flux_njy_aperture": flux, "n_exposures": nexp, "pa_aper_deg": pa, **ps}


def main() -> int:
    warnings.simplefilter("ignore", FITSFixedWarning)
    root = paths.data_root() / "m0647_1433"
    out: dict = {"provenance": "derived (post hoc; not in PLAN.md)", "position": [RA, DEC]}
    out["aperture_arcsec"], out["annulus_arcsec"] = APER_ARCSEC, list(ANNULUS_ARCSEC)
    for e in ("o010", "o020"):
        out[f"f200w_{e}"] = aperture_njy(root / f"jw01433-{e}_t001_nircam_clear-f200w_i2d.fits")
    d = out["f200w_o020"]["flux_njy_aperture"] - out["f200w_o010"]["flux_njy_aperture"]
    out["f200w_difference_njy_aperture"] = d
    out["f200w_difference_abmag_aperture"] = float(-2.5 * np.log10(d * 1e-9 / 3631.0))
    out["roll_difference_deg"] = abs(
        out["f200w_o010"]["pa_aper_deg"] - out["f200w_o020"]["pa_aper_deg"]
    )
    out["f480m_o020"] = aperture_njy(root / "cal" / "jw01433-o020_t001_nircam_clear-f480m_i2d.fits")
    rows = []
    for name in CAL:
        with fits.open(root / "cal" / name, memmap=True) as h:
            w = WCS(h["SCI"].header)
            x, y = (float(v) for v in w.world_to_pixel_values(RA, DEC))
            t0 = float(h[0].header["EXPSTART"])
            rows.append({"file": name, "expstart_mjd": t0, **peak_snr(h["SCI"].data, x, y)})
    out["cal_exposures"] = rows
    out["crossmatch"] = CROSSMATCH
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
