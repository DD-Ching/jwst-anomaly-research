"""Is a line through a W1 flag in one level-2 exposure only (a trail) or in all of them (sky)?

Frozen test: results/w1_shear_primer/vet/x3y0_cal_plan.md. Every ``_cal`` exposure that covers the
peak is resampled (bilinear, through each exposure's FITS WCS) onto the mosaic's pixel grid around
the peak, so the line has the same row in every stamp when it is on the sky and a different row
(or none) when it belongs to one exposure or to the detector. Per stamp, the band rows around the
line are compared with two control bands ``CONTROL_ARCSEC`` away along the mosaic's columns:

    contrast = (median of band rows − median of control rows) / robust σ of control row medians

    python scripts/w1_cal_trail.py --ra 34.50675170212336 --dec -5.259761111111111 \\
        --mosaic s3://stpubdata/jwst/public/jw01837/L3/t/c1002/jw01837-c1002_t000_nircam_clear-f200w_i2d.fits \\
        --cal s3://stpubdata/jwst/public/jw01837/jw01837001001/jw01837001001_08201_00001_nrca4_cal.fits ... \\
        --out results/w1_shear_primer/vet

Pixels ``observed``; stamps, profiles and contrasts ``derived``; every threshold an ASSUMPTION.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from persistence_check import _open  # noqa: E402

HALF_ARCSEC = 20.0  # stamp half-width along the line (the 40″ visit sheet; ASSUMPTION)
HALF_ROWS_ARCSEC = 5.0  # stamp half-height across the line (ASSUMPTION)
BAND_PX = 3  # band = line row ± this many mosaic pixels (plan; ASSUMPTION)
CONTROL_ARCSEC = 3.0  # control bands this far above and below (plan; ASSUMPTION)
IN_SIGMA = 5.0  # a line is "in" a stamp above this contrast (plan; ASSUMPTION)
DQ_DO_NOT_USE = 1


def mosaic_grid(
    uri: str, ra: float, dec: float, half_arcsec: float = HALF_ARCSEC
) -> tuple[WCS, float, float, float, np.ndarray]:
    """Mosaic WCS, the peak's pixel, the pixel scale (arcsec) and the mosaic stamp."""
    with _open(uri) as hdul, warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        hdu = hdul["SCI"]
        wcs = WCS(hdu.header, naxis=2)
        scale = float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix))) * 3600.0)
        x0, y0 = (float(v) for v in wcs.world_to_pixel_values(ra, dec))
        hx, hy = int(half_arcsec / scale), int(HALF_ROWS_ARCSEC / scale)
        xi, yi = int(round(x0)), int(round(y0))
        stamp = np.asarray(hdu.section[yi - hy : yi + hy + 1, xi - hx : xi + hx + 1], dtype=float)
    return wcs, x0, y0, scale, stamp


def resample_cal(uri: str, sky_x: np.ndarray, sky_y: np.ndarray) -> dict:
    """The exposure's SCI (DQ DO_NOT_USE as NaN) sampled at the given sky positions."""
    with _open(uri) as hdul, warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        h0, hdu = hdul[0].header, hdul["SCI"]
        wcs = WCS(hdu.header, naxis=2)
        px, py = wcs.world_to_pixel_values(sky_x, sky_y)
        ny, nx = hdu.header["NAXIS2"], hdu.header["NAXIS1"]
        ok = np.isfinite(px) & np.isfinite(py)
        x_lo = int(max(0, np.floor(px[ok].min()) - 2))
        x_hi = int(min(nx, np.ceil(px[ok].max()) + 3))
        y_lo = int(max(0, np.floor(py[ok].min()) - 2))
        y_hi = int(min(ny, np.ceil(py[ok].max()) + 3))
        if x_lo >= x_hi or y_lo >= y_hi:
            return {"file": Path(uri).name, "stamp": np.full(px.shape, np.nan)}
        sci = np.asarray(hdu.section[y_lo:y_hi, x_lo:x_hi], dtype=float)
        dq = np.asarray(hdul["DQ"].section[y_lo:y_hi, x_lo:x_hi])
        sci[(dq & DQ_DO_NOT_USE) != 0] = np.nan
        info = {
            "file": Path(uri).name,
            "detector": h0.get("DETECTOR"),
            "expstart_mjd": h0.get("EXPSTART"),
            "peak_detector_px": None,
        }
    bad = ~np.isfinite(sci)
    filled = np.where(bad, 0.0, sci)
    coords = [py - y_lo, px - x_lo]
    stamp = ndimage.map_coordinates(filled, coords, order=1, cval=np.nan)
    badmap = ndimage.map_coordinates(bad.astype(float), coords, order=1, cval=1.0)
    outside = (px < 0) | (px > nx - 1) | (py < 0) | (py > ny - 1)
    stamp[(badmap > 0.01) | outside | ~ok] = np.nan
    info["stamp"] = stamp
    return info


def row_profile(stamp: np.ndarray, min_cover: float = 0.3) -> np.ndarray:
    """Median per stamp row (NaN where under ``min_cover`` of the row is finite)."""
    finite = np.isfinite(stamp)
    prof = np.full(stamp.shape[0], np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(stamp, axis=1)
    keep = finite.mean(axis=1) >= min_cover
    prof[keep] = med[keep]
    return prof


def contrast(prof: np.ndarray, centre: int, scale: float) -> dict:
    """Band versus control contrast of one row profile (rows indexed from the stamp's top)."""
    off = int(round(CONTROL_ARCSEC / scale))
    band = prof[centre - BAND_PX : centre + BAND_PX + 1]
    ctrl = np.concatenate(
        [
            prof[centre - off - BAND_PX : centre - off + BAND_PX + 1],
            prof[centre + off - BAND_PX : centre + off + BAND_PX + 1],
        ]
    )
    ctrl = ctrl[np.isfinite(ctrl)]
    band = band[np.isfinite(band)]
    if band.size < 2 or ctrl.size < 4:
        return {"contrast": None, "n_band": int(band.size), "n_control": int(ctrl.size)}
    sigma = 1.4826 * np.median(np.abs(ctrl - np.median(ctrl)))
    sigma = max(sigma, 1e-12)
    peak = float(np.max(band) - np.median(ctrl))
    return {
        "contrast": float((np.median(band) - np.median(ctrl)) / sigma),
        "peak_row_contrast": float(peak / sigma),
        "band_median": float(np.median(band)),
        "control_median": float(np.median(ctrl)),
        "control_sigma": float(sigma),
        "n_band": int(band.size),
        "n_control": int(ctrl.size),
    }


def line_track(stamp: np.ndarray, centre: int, half: int = 20, chunk: int = 50) -> dict:
    """Post hoc: the line's row per column chunk (brightest median row within ``centre ± half``)
    and a straight-line fit, so the line's tilt can be compared between stamps (a line fixed on
    the sky has the same row and tilt in every stamp)."""
    xs, rows, amps = [], [], []
    for x0 in range(0, stamp.shape[1] - chunk + 1, chunk):
        sub = stamp[centre - half : centre + half + 1, x0 : x0 + chunk]
        if np.isfinite(sub).mean() < 0.5:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            prof = np.nanmedian(sub, axis=1)
        if not np.isfinite(prof).any():
            continue
        xs.append(x0 + chunk / 2)
        rows.append(float(np.nanargmax(prof) - half))
        bg = np.nanmedian(prof)
        amps.append(float(np.nanmax(prof) - bg))
    if len(xs) < 3:
        return {"n_chunks": len(xs)}
    coef, cov = np.polyfit(xs, rows, 1, cov=True)
    resid = np.asarray(rows) - np.polyval(coef, xs)
    return {
        "n_chunks": len(xs),
        "slope_px_per_px": float(coef[0]),
        "slope_err": float(np.sqrt(cov[0, 0])),
        "row_at_stamp_centre_px": float(np.polyval(coef, stamp.shape[1] / 2)),
        "rms_resid_px": float(np.sqrt(np.mean(resid**2))),
        "chunk_x_px": [float(x) for x in xs],
        "chunk_row_px": rows,
        "chunk_peak_above_median": amps,
    }


def classify(line_in: list[bool]) -> str:
    """The plan's frozen outcome from the per-exposure calls (``contrast > IN_SIGMA``)."""
    if not line_in:
        return "undecided: no exposure covers the band and controls"
    n_in = sum(line_in)
    if n_in == len(line_in):
        return "in every exposure: not a single-exposure trail (stays a flag)"
    if n_in == 0:
        return "in no exposure above the threshold: the mosaic line is not reproduced"
    return f"in {n_in} of {len(line_in)} exposures: not a feature of every exposure (trail or exposure-level artefact)"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--mosaic", required=True)
    ap.add_argument("--cal", nargs="+", required=True)
    ap.add_argument("--line-row-offset", type=int, default=1, help="line row minus peak row (mosaic px)")
    ap.add_argument("--name", default="x3y0")
    ap.add_argument("--half-arcsec", type=float, default=HALF_ARCSEC, help="stamp half-width along the line")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    wcs, x0, y0, scale, mstamp = mosaic_grid(args.mosaic, args.ra, args.dec, args.half_arcsec)
    hy, hx = (s // 2 for s in mstamp.shape)
    xi, yi = int(round(x0)), int(round(y0))
    gy, gx = np.mgrid[yi - hy : yi + hy + 1, xi - hx : xi + hx + 1]
    sky_x, sky_y = wcs.pixel_to_world_values(gx, gy)
    with ThreadPoolExecutor(len(args.cal)) as ex:
        cals = list(ex.map(lambda u: resample_cal(u, sky_x, sky_y), args.cal))

    # rows are indexed from the stamp's first (lowest-y) row; the line is at the peak row + offset
    centre = hy + args.line_row_offset
    out = {
        "question": "is the line through the peak present in one exposure only (trail) or in all (sky)?",
        "plan": "results/w1_shear_primer/vet/x3y0_cal_plan.md",
        "peak_radec": [args.ra, args.dec],
        "mosaic": Path(args.mosaic).name,
        "mosaic_scale_arcsec": scale,
        "half_arcsec": args.half_arcsec,
        "line_row_offset_px": args.line_row_offset,
        "band_px": BAND_PX,
        "control_arcsec": CONTROL_ARCSEC,
        "in_sigma": IN_SIGMA,
        "stamps": [],
    }
    profiles = {"mosaic": row_profile(mstamp)}
    for name, stamp in [("mosaic", mstamp)] + [(c["file"], c["stamp"]) for c in cals]:
        prof = profiles.get(name) if name == "mosaic" else row_profile(stamp)
        profiles[name] = prof
        meta = next((c for c in cals if c["file"] == name), {})
        rec = {
            "name": name,
            "detector": meta.get("detector"),
            "expstart_mjd": meta.get("expstart_mjd"),
            "finite_fraction": float(np.isfinite(stamp).mean()),
            "brightest_row_offset_px": (
                int(np.nanargmax(prof[centre - 20 : centre + 21]) - 20 + args.line_row_offset)
                if np.isfinite(prof[centre - 20 : centre + 21]).any()
                else None
            ),
            **contrast(prof, centre, scale),
        }
        c = rec.get("contrast")
        rec["line_in"] = None if c is None else bool(c > IN_SIGMA)
        rec["track_post_hoc"] = line_track(stamp, centre)
        out["stamps"].append(rec)
    exp = [s for s in out["stamps"] if s["name"] != "mosaic" and s["line_in"] is not None]
    out["n_exposures_tested"] = len(exp)
    out["n_exposures_with_line"] = sum(s["line_in"] for s in exp)
    out["outcome"] = classify([s["line_in"] for s in exp])
    out["provenance"] = {"pixels": "observed", "stamps_profiles_contrasts": "derived", "thresholds": "assumption"}

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / f"{args.name}_cal_check.json").write_text(json.dumps(out, indent=1) + "\n")
    _sheet(args.out / f"{args.name}_cal_sheet.png", mstamp, cals, profiles, centre, scale, out)
    print(json.dumps({k: out[k] for k in ("n_exposures_tested", "n_exposures_with_line", "outcome")}))
    for s in out["stamps"]:
        print(s["name"], s.get("contrast"), s.get("line_in"), s.get("brightest_row_offset_px"), s["track_post_hoc"])
    return 0


def _sheet(path, mstamp, cals, profiles, centre, scale, out) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stamps = [("mosaic", mstamp)] + [(c["file"], c["stamp"]) for c in cals]
    fig, axes = plt.subplots(len(stamps), 2, figsize=(13, 2.2 * len(stamps)), width_ratios=[4, 1])
    for (name, st), (ax, axp), rec in zip(stamps, axes, out["stamps"], strict=True):
        fin = st[np.isfinite(st)]
        lo, hi = (np.percentile(fin, [5, 99.5]) if fin.size else (0, 1))
        ax.imshow(st, origin="lower", cmap="gray_r", vmin=lo, vmax=hi, aspect="equal")
        ax.axhline(centre, color="tab:red", lw=0.5, alpha=0.5)
        off = int(round(CONTROL_ARCSEC / scale))
        for c in (centre - off, centre + off):
            ax.axhline(c, color="tab:blue", lw=0.5, alpha=0.5)
        ax.set_title(f"{name}  contrast {rec.get('contrast')!s:.6}", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        prof = profiles[name]
        axp.plot(prof, np.arange(prof.size), lw=0.7)
        axp.axhline(centre, color="tab:red", lw=0.5)
        axp.set_ylim(0, prof.size)
        axp.tick_params(labelsize=6)
    fig.suptitle(f"{out['mosaic']} vs _cal exposures; rows = mosaic rows; {out['outcome']}", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
