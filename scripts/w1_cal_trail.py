"""Is a line through a W1 flag in one level-2 exposure only (a trail) or in all of them (sky)?

Frozen test: results/w1_shear_primer/vet/x3y0_cal_plan.md. Every ``_cal`` exposure that covers the
peak is resampled (bilinear, through each exposure's FITS WCS) onto the mosaic's pixel grid around
the peak, so the line has the same row in every stamp when it is on the sky and a different row
(or none) when it belongs to one exposure or to the detector. Per stamp, the band rows around the
line are compared with two control bands ``CONTROL_ARCSEC`` away along the mosaic's columns:

    contrast = (median of band rows − median of control rows) / robust σ of control row medians

    S=s3://stpubdata/jwst/public/jw01837
    python scripts/w1_cal_trail.py --ra 34.50675170212336 --dec -5.259761111111111 \\
        --mosaic $S/L3/t/c1002/jw01837-c1002_t000_nircam_clear-f200w_i2d.fits \\
        --search-dirs $S/jw01837001001 ... $S/jw01837001022 --out results/w1_shear_primer/vet

``--search-dirs`` finds the covering exposures (or pass ``--cal``). The post hoc run adds
``--half-arcsec 70 --star RA DEC --star-mag G --star-id ID`` (results/w1_shear_primer/RESULT.md).

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
from astropy.wcs import WCS, FITSFixedWarning
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from persistence_check import _open  # noqa: E402

HALF_ARCSEC = 10.0  # stamp half-width along the line (plan: ±10″; ASSUMPTION)
HALF_ROWS_ARCSEC = 5.0  # stamp half-height across the line (ASSUMPTION)
BAND_PX = 3  # band = line row ± this many mosaic pixels (plan; ASSUMPTION)
CONTROL_ARCSEC = 3.0  # control bands this far above and below (plan; ASSUMPTION)
IN_SIGMA = 5.0  # a line is "in" a stamp above this contrast (plan; ASSUMPTION)
DQ_DO_NOT_USE = 1
SIP_MARGIN_PX = 100  # full (SIP) inversion only within this margin of the detector (ASSUMPTION)
COVER_MARGIN_PX = 50  # an exposure "covers" the peak within this margin (ASSUMPTION)
SPIKE_LAW = (3.0, 20.0)  # spike length 3″ × 10^(0.2 (20 − m)): exotic_screens.spike_radius


def dither_of(name: str) -> str:
    """``00002`` from ``jw01837001011_08201_00002_nrcb3_cal.fits`` (exposure number in a visit)."""
    parts = Path(name).name.split("_")
    return parts[2] if len(parts) > 3 else name


def find_cal(dirs: list[str], ra: float, dec: float, band: str) -> list[dict]:
    """Every SW ``_cal`` exposure in ``dirs`` taken in ``band`` whose detector covers (ra, dec)."""
    import re

    import fsspec
    from astropy.io import fits
    from persistence_check import list_cal

    with ThreadPoolExecutor(16) as ex:
        files = [f for fs in ex.map(list_cal, dirs) for f in fs]
    files = [f for f in files if re.search(r"_nrc[ab][1-4]_cal\.fits$", f)]

    def head(uri: str) -> tuple[fits.Header, fits.Header]:
        fs, path = fsspec.core.url_to_fs(uri, **({"anon": True} if uri.startswith("s3://") else {}))
        with fs.open(path, block_size=2**16) as f, fits.open(f) as hdul:
            return hdul[0].header.copy(), hdul["SCI"].header.copy()

    groups: dict[str, str] = {}
    for f in files:  # one exposure group = one filter: read one header per group
        groups.setdefault(Path(f).name.rsplit("_", 3)[0], f)
    with ThreadPoolExecutor(16) as ex:
        heads = dict(zip(groups, ex.map(lambda u: head(u)[0], groups.values()), strict=True))
    keep = {k for k, h in heads.items() if band in (h.get("FILTER"), h.get("PUPIL"))}
    cand = [f for f in files if Path(f).name.rsplit("_", 3)[0] in keep]

    def cover(uri: str) -> dict | None:
        h0, h = head(uri)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            x, y = (float(v) for v in WCS(h, naxis=2).wcs_world2pix(ra, dec, 0))
        m = COVER_MARGIN_PX
        if not (-m < x < h["NAXIS1"] + m and -m < y < h["NAXIS2"] + m):
            return None
        return {
            "uri": uri,
            "file": Path(uri).name,
            "detector": h0.get("DETECTOR"),
            "filter": h0.get("FILTER"),
            "pupil": h0.get("PUPIL"),
            "expstart_mjd": h0.get("EXPSTART"),
            "pa_aper": h.get("PA_APER"),
            "x_linear_wcs": x,
            "y_linear_wcs": y,
        }

    with ThreadPoolExecutor(16) as ex:
        found = [c for c in ex.map(cover, cand) if c is not None]
    return sorted(found, key=lambda c: c["expstart_mjd"])


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
        ny, nx = hdu.header["NAXIS2"], hdu.header["NAXIS1"]
        if xi - hx < 0 or yi - hy < 0 or xi + hx + 1 > nx or yi + hy + 1 > ny:
            raise ValueError(f"stamp of ±{half_arcsec}″ runs off the mosaic ({nx} × {ny})")
        stamp = np.asarray(hdu.section[yi - hy : yi + hy + 1, xi - hx : xi + hx + 1], dtype=float)
    return wcs, x0, y0, scale, stamp


def resample_cal(uri: str, sky_x: np.ndarray, sky_y: np.ndarray) -> dict:
    """The exposure's SCI (DQ DO_NOT_USE as NaN) sampled at the given sky positions."""
    with _open(uri) as hdul, warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        h0, hdu = hdul[0].header, hdul["SCI"]
        wcs = WCS(hdu.header, naxis=2)
        ny, nx = hdu.header["NAXIS2"], hdu.header["NAXIS1"]
        # SIP diverges far off the detector: linear WCS first, full inversion only near it
        lx, ly = wcs.wcs_world2pix(sky_x, sky_y, 0)
        m = SIP_MARGIN_PX
        near = (lx > -m) & (lx < nx + m) & (ly > -m) & (ly < ny + m)
        px = np.full(sky_x.shape, np.nan)
        py = np.full(sky_x.shape, np.nan)
        if near.any():
            fx, fy = wcs.all_world2pix(sky_x[near], sky_y[near], 0, quiet=True)
            px[near], py[near] = fx, fy
        ok = np.isfinite(px) & np.isfinite(py)
        if not ok.any():
            return {"file": Path(uri).name, "stamp": np.full(px.shape, np.nan)}
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
            "pa_aper": hdu.header.get("PA_APER"),
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
    if centre - off - BAND_PX < 0 or centre + off + BAND_PX + 1 > prof.size:
        raise ValueError("control bands run off the stamp: raise HALF_ROWS_ARCSEC")
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
        "coef": [float(c) for c in coef],
        "chunk_x_px": [float(x) for x in xs],
        "chunk_row_px": rows,
        "chunk_peak_above_median": amps,
    }


def star_geometry(
    wcs: WCS, star, x_origin: int, row_origin: int, stamps: list[dict], x_centre: float = 0.0
) -> dict:
    """Post hoc: the star's stamp position and how far each stamp's fitted line passes from it
    (rows at the star's column), so a line radial from the star can be told apart."""
    sx, sy = (float(v) for v in wcs.world_to_pixel_values(*star))
    sx, srow = sx - x_origin, sy - row_origin
    out = {
        "radec": list(star),
        "stamp_x_px": sx,
        "row_rel_line_centre_px": srow,
        "x_from_stamp_centre_px": x_centre - sx,
        "tracks": {},
    }
    for s in stamps:
        t = s["track_post_hoc"]
        if "coef" not in t or t["rms_resid_px"] > 2.0:
            continue
        coef = t["coef"]
        out["tracks"][s["name"]] = {
            "line_row_at_star_x_px": float(np.polyval(coef, sx)),
            "miss_distance_px": float(abs(np.polyval(coef, sx) - srow)),
        }
    return out


def classify(line_in: list[bool], dithers: list[str] | None = None) -> str:
    """The plan's frozen outcome from the per-exposure calls (``contrast > IN_SIGMA``) and each
    exposure's dither (exposure number within its visit). The plan's third branch ("fixed in
    detector pixels") has no computable threshold in the plan and is not decided here."""
    if not line_in:
        return "undecided: no exposure covers the band and controls"
    n_in = sum(line_in)
    if n_in == len(line_in):
        return "in every exposure: not a trail (stays a flag)"
    if n_in == 0:
        return "undecided: in no exposure above the threshold"
    if n_in == 1:
        return "trail: in exactly one exposure (x3y0 explained)"
    if dithers is not None and len({d for d, i in zip(dithers, line_in, strict=True) if i}) == 1:
        return f"trail: only in exposures of one dither position ({n_in} of {len(line_in)})"
    return f"in {n_in} of {len(line_in)} exposures from several dither positions: no frozen branch"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--mosaic", required=True)
    ap.add_argument("--cal", nargs="+", help="_cal URIs (else --search-dirs)")
    ap.add_argument("--search-dirs", nargs="+", help="visit directories searched for _cal files")
    ap.add_argument("--band", default="F200W")
    ap.add_argument(
        "--line-row-offset", type=int, default=0, help="line row minus peak row (mosaic px)"
    )
    ap.add_argument("--name", default="x3y0")
    ap.add_argument(
        "--star", type=float, nargs=2, metavar=("RA", "DEC"), help="post hoc: a suspected source"
    )
    ap.add_argument("--star-mag", type=float, help="the star's magnitude (e.g. Gaia G), recorded")
    ap.add_argument("--star-id", help="the star's catalogue id, recorded")
    ap.add_argument(
        "--half-arcsec", type=float, default=HALF_ARCSEC, help="stamp half-width along the line"
    )
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    found = None
    if args.cal is None:
        if not args.search_dirs:
            ap.error("give --cal or --search-dirs")
        found = find_cal(args.search_dirs, args.ra, args.dec, args.band)
        args.cal = [c["uri"] for c in found]
        if not args.cal:
            print("no exposure covers the position")
            return 1
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
        "question": (
            "is the line through the peak present in one exposure only (trail) or in all (sky)?"
        ),
        "plan": "results/w1_shear_primer/vet/x3y0_cal_plan.md",
        "peak_radec": [args.ra, args.dec],
        "mosaic": Path(args.mosaic).name,
        "mosaic_scale_arcsec": scale,
        "half_arcsec": args.half_arcsec,
        "line_row_offset_px": args.line_row_offset,
        "search_dirs": args.search_dirs,
        "exposures_found": found,
        "band_px": BAND_PX,
        "control_arcsec": CONTROL_ARCSEC,
        "in_sigma": IN_SIGMA,
        "stamps": [],
    }
    profiles = []
    for i, (name, stamp) in enumerate(
        [("mosaic", mstamp)] + [(c["file"], c["stamp"]) for c in cals]
    ):
        prof = row_profile(stamp)
        profiles.append(prof)
        meta = cals[i - 1] if i else {}
        rec = {
            "name": name,
            "uri": args.cal[i - 1] if i else args.mosaic,
            "dither": dither_of(name) if i else None,
            "detector": meta.get("detector"),
            "expstart_mjd": meta.get("expstart_mjd"),
            "pa_aper": meta.get("pa_aper"),
            "finite_fraction": float(np.isfinite(stamp).mean()),
            "brightest_row_offset_px": (
                int(np.nanargmax(prof[centre - 20 : centre + 21]) - 20 + args.line_row_offset)
                if centre >= 20 and np.isfinite(prof[centre - 20 : centre + 21]).any()
                else None
            ),
            **contrast(prof, centre, scale),
        }
        c = rec.get("contrast")
        rec["line_in"] = None if c is None else bool(c > IN_SIGMA)
        rec["track_post_hoc"] = line_track(stamp, centre)
        out["stamps"].append(rec)
    if args.star:
        geo = star_geometry(wcs, args.star, xi - hx, centre + yi - hy, out["stamps"], hx)
        sep = float(
            np.hypot(
                (args.star[0] - args.ra) * np.cos(np.radians(args.dec)), args.star[1] - args.dec
            )
            * 3600.0
        )
        geo |= {"id": args.star_id, "mag": args.star_mag, "sep_from_peak_arcsec": sep}
        if args.star_mag is not None:
            a, m0 = SPIKE_LAW
            geo["spike_law_uncapped_arcsec"] = a * 10 ** (0.2 * (m0 - args.star_mag))
        out["star_post_hoc"] = geo
    exp = [s for s in out["stamps"][1:] if s["line_in"] is not None]
    out["n_exposures_tested"] = len(exp)
    out["n_exposures_with_line"] = sum(s["line_in"] for s in exp)
    out["outcome"] = classify([s["line_in"] for s in exp], [s["dither"] for s in exp])
    out["provenance"] = {
        "pixels": "observed",
        "stamps_profiles_contrasts": "derived",
        "thresholds": "assumption",
    }

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / f"{args.name}_cal_check.json").write_text(json.dumps(out, indent=1) + "\n")
    _sheet(args.out / f"{args.name}_cal_sheet.png", mstamp, cals, profiles, centre, scale, out)
    print(
        json.dumps({k: out[k] for k in ("n_exposures_tested", "n_exposures_with_line", "outcome")})
    )
    for s in out["stamps"]:
        print(s["name"], s.get("contrast"), s.get("line_in"), s.get("brightest_row_offset_px"))
    return 0


def _sheet(path, mstamp, cals, profiles, centre, scale, out) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stamps = [("mosaic", mstamp)] + [(c["file"], c["stamp"]) for c in cals]
    fig, axes = plt.subplots(len(stamps), 2, figsize=(13, 2.2 * len(stamps)), width_ratios=[4, 1])
    for (name, st), (ax, axp), rec, prof in zip(stamps, axes, out["stamps"], profiles, strict=True):
        fin = st[np.isfinite(st)]
        lo, hi = np.percentile(fin, [5, 99.5]) if fin.size else (0, 1)
        ax.imshow(st, origin="lower", cmap="gray_r", vmin=lo, vmax=hi, aspect="equal")
        ax.axhline(centre, color="tab:red", lw=0.5, alpha=0.5)
        off = int(round(CONTROL_ARCSEC / scale))
        for c in (centre - off, centre + off):
            ax.axhline(c, color="tab:blue", lw=0.5, alpha=0.5)
        ax.set_title(f"{name}  contrast {rec.get('contrast')!s:.6}", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        axp.plot(prof, np.arange(prof.size), lw=0.7)
        axp.axhline(centre, color="tab:red", lw=0.5)
        axp.set_ylim(0, prof.size)
        axp.tick_params(labelsize=6)
    fig.suptitle(
        f"{out['mosaic']} vs _cal exposures; rows = mosaic rows; {out['outcome']}", fontsize=8
    )
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
