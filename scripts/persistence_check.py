"""Detector-persistence test for transient candidates on level-2 ``_cal`` exposures (D-039).

A source that appears in one epoch only can be an afterimage: a bright source illuminates a
detector pixel in one exposure, and the pixel still shows ~0.05–3 % of that signal in the next
exposures (Sunrise ``n0022`` and ``n0150``, docs/fields/sunrise.md). The level-3 mosaic maps
that pixel to a different sky position after each dither, so the afterimage looks like a faint,
compact source at a fixed or wandering sky position.

For each position this script lists every ``_cal`` exposure that covers it (headers read from S3
or local paths, no full downloads), measures the aperture flux there, and measures the same
*detector pixel* in every earlier exposure on the same detector within ``--lookback-hours``. A
detection (S/N ≥ ``--min-snr``) is ``suspect`` when an earlier exposure put
≥ ``--ratio`` × its flux, or any saturated pixel, on that pixel. Per position:

- ``persistence``: every detection is suspect;
- ``on_sky``: at least two detections are not suspect;
- ``inconclusive``: one clean detection; ``undetected``: none.

Every number is ``derived``; thresholds are ASSUMPTIONs (D-039).

    python scripts/persistence_check.py \\
        --position n0022 24.364244 -8.433747 --position n0150 24.339250 -8.442280 \\
        --cal-dir s3://stpubdata/jwst/public/jw02282/jw02282010001 --out outputs/persistence
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transient_forced import aperture_flux  # noqa: E402

RATIO = 20.0  # ASSUMPTION: earlier flux ≥ 20x the detection (observed afterimages: 0.05–3 %)
MIN_SNR = 5.0  # ASSUMPTION: ERR-based S/N of a detection
LOOKBACK_HOURS = 3.0  # ASSUMPTION: NIRCam afterimages above noise fade within hours
SATURATED = 2  # JWST DQ bit
HALF = 15  # cutout half-size, pixels
_PAD = 5  # accept positions up to this many pixels outside the array (cutout still useful)


def classify(snr: np.ndarray, suspect: np.ndarray, min_snr: float = MIN_SNR) -> str:
    """Verdict for one position from its per-exposure S/N and suspect flags."""
    snr = np.asarray(snr, dtype=float)
    det = np.isfinite(snr) & (snr >= min_snr)
    if not det.any():
        return "undetected"
    clean = int((det & ~np.asarray(suspect, dtype=bool)).sum())
    if clean == 0:
        return "persistence"
    return "on_sky" if clean >= 2 else "inconclusive"


def is_suspect(flux: float, prior_flux: float, prior_sat: int, ratio: float = RATIO) -> bool:
    """True when an earlier exposure put far more signal (or saturation) on this pixel."""
    if prior_sat > 0:
        return True
    return bool(
        np.isfinite(flux) and flux > 0 and np.isfinite(prior_flux) and prior_flux >= ratio * flux
    )


def _fs(uri: str):
    import fsspec

    opts = {"anon": True} if uri.startswith("s3://") else {}
    return fsspec.core.url_to_fs(uri, **opts)


def list_cal(cal_dir: str) -> list[str]:
    fs, path = _fs(cal_dir)
    prefix = "s3://" if cal_dir.startswith("s3://") else ""
    return sorted(prefix + p for p in fs.ls(path, detail=False) if p.endswith("_cal.fits"))


def _open(uri: str):
    fs, path = _fs(uri)
    return fs.open(path, "rb", block_size=2**20, cache_type="readahead")


def header_info(uri: str, positions: dict[str, tuple[float, float]]) -> dict[str, Any] | None:
    """Detector, filter, start time and the pixel of each position on this exposure.

    None for files that are not 2-D imaging exposures (no ``EXPSTART`` or a non-2-D ``SCI``)."""
    with _open(uri) as fo, fits.open(fo, lazy_load_hdus=True) as hdul, warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        primary, sci = hdul[0].header, hdul["SCI"].header
        if "EXPSTART" not in primary or sci.get("NAXIS") != 2:
            return None
        wcs = WCS(sci, naxis=2)
        nx, ny = sci["NAXIS1"], sci["NAXIS2"]
        pix = {}
        for uid, (ra, dec) in positions.items():
            x, y = (
                float(v) for v in wcs.wcs_world2pix(ra, dec, 0)
            )  # linear first: SIP diverges far off
            if not (-100 < x < nx + 100 and -100 < y < ny + 100):
                continue
            x, y = (float(v) for v in wcs.all_world2pix(ra, dec, 0))
            if -_PAD < x < nx + _PAD and -_PAD < y < ny + _PAD:
                pix[uid] = (x, y)
    filt = str(primary.get("FILTER", "")).strip()
    pupil = str(primary.get("PUPIL", "")).strip()
    return {
        "file": uri,
        "detector": str(primary.get("DETECTOR", "")).strip(),
        "filter": pupil if pupil.startswith("F") else filt,
        "mjd": float(primary["EXPSTART"]),
        "pix": pix,
    }


def measure_pixel(uri: str, x: float, y: float, long: bool) -> dict[str, float]:
    """Aperture flux, ERR-based error and saturated-pixel count at detector pixel (x, y)."""
    xi, yi = int(round(x)), int(round(y))
    with _open(uri) as fo, fits.open(fo, lazy_load_hdus=True) as hdul:
        ny, nx = hdul["SCI"].header["NAXIS2"], hdul["SCI"].header["NAXIS1"]
        y0, y1, x0, x1 = (
            max(0, yi - HALF),
            min(ny, yi + HALF + 1),
            max(0, xi - HALF),
            min(nx, xi + HALF + 1),
        )
        sci = np.asarray(hdul["SCI"].section[y0:y1, x0:x1], dtype=float)
        err = np.asarray(hdul["ERR"].section[y0:y1, x0:x1], dtype=float)
        dq = np.asarray(hdul["DQ"].section[y0:y1, x0:x1])
    r_px = 2.0 if long else 3.0  # ~0.1" in both channels
    flux, e = aperture_flux(sci, err, x - x0, y - y0, r_px, 2 * r_px + 1, 4 * r_px)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    core = np.hypot(xx - x, yy - y) <= r_px
    return {"flux": flux, "err": e, "n_sat": int(((dq[core] & SATURATED) > 0).sum())}


def run(
    positions: dict[str, tuple[float, float]],
    cal_files: Sequence[str],
    ratio: float = RATIO,
    min_snr: float = MIN_SNR,
    lookback_hours: float = LOOKBACK_HOURS,
    workers: int = 12,
) -> tuple[Table, Table]:
    with ThreadPoolExecutor(workers) as ex:
        info = [i for i in ex.map(lambda f: header_info(f, positions), cal_files) if i is not None]
    info.sort(key=lambda r: r["mjd"])
    jobs = []
    for i, cur in enumerate(info):
        for uid, (x, y) in cur["pix"].items():
            prior = [
                p
                for p in info[:i]
                if p["detector"] == cur["detector"]
                and 0 < cur["mjd"] - p["mjd"] <= lookback_hours / 24
            ]
            jobs.append((uid, cur, prior, x, y))

    def one(job):
        uid, cur, prior, x, y = job
        long = "LONG" in cur["detector"].upper()
        now = measure_pixel(cur["file"], x, y, long)
        before = [(p, measure_pixel(p["file"], x, y, long)) for p in prior]
        worst = max(
            before,
            key=lambda b: (b[1]["n_sat"], np.nan_to_num(b[1]["flux"], nan=-np.inf)),
            default=None,
        )
        prior_flux = worst[1]["flux"] if worst else np.nan
        prior_sat = sum(b[1]["n_sat"] for b in before)
        snr = now["flux"] / now["err"] if now["err"] > 0 else np.nan
        return {
            "uid": uid,
            "file": cur["file"].rsplit("/", 1)[-1],
            "detector": cur["detector"],
            "filter": cur["filter"],
            "mjd": cur["mjd"],
            "x": x,
            "y": y,
            "flux": now["flux"],
            "snr": snr,
            "prior_file": worst[0]["file"].rsplit("/", 1)[-1] if worst else "",
            "prior_flux": prior_flux,
            "prior_sat": prior_sat,
            "suspect": bool(
                np.isfinite(snr)
                and snr >= min_snr
                and is_suspect(now["flux"], prior_flux, prior_sat, ratio)
            ),
        }

    with ThreadPoolExecutor(workers) as ex:
        rows = list(ex.map(one, jobs))
    per_exp = (
        Table(rows=rows)
        if rows
        else Table(names=["uid", "snr", "suspect"], dtype=[str, float, bool])
    )
    summary = []
    for uid in positions:
        sel = per_exp[per_exp["uid"] == uid] if len(per_exp) else per_exp
        snr = np.asarray(sel["snr"], dtype=float)
        sus = np.asarray(sel["suspect"], dtype=bool)
        det = np.isfinite(snr) & (snr >= min_snr)
        summary.append(
            {
                "uid": uid,
                "n_exposures": len(sel),
                "n_detected": int(det.sum()),
                "n_suspect": int((det & sus).sum()),
                "verdict": classify(snr, sus, min_snr),
            }
        )
    meta = {
        "provenance": "derived",
        "source": "level-2 _cal exposures: "
        + ", ".join(sorted({f.rsplit("/", 1)[0] for f in cal_files})),
        "ratio": ratio,
        "min_snr": min_snr,
        "lookback_hours": lookback_hours,
    }
    per_exp.meta.update(meta)
    summ = Table(rows=summary)
    summ.meta.update(meta)
    return per_exp, summ


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--position", nargs=3, action="append", metavar=("UID", "RA", "DEC"), required=True
    )
    ap.add_argument(
        "--cal-dir", action="append", required=True, help="directory of _cal.fits (s3:// or local)"
    )
    ap.add_argument("--ratio", type=float, default=RATIO)
    ap.add_argument("--min-snr", type=float, default=MIN_SNR)
    ap.add_argument("--lookback-hours", type=float, default=LOOKBACK_HOURS)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    positions = {u: (float(r), float(d)) for u, r, d in args.position}
    files = [f for d in args.cal_dir for f in list_cal(d)]
    per_exp, summ = run(positions, files, args.ratio, args.min_snr, args.lookback_hours)
    args.out.mkdir(parents=True, exist_ok=True)
    per_exp.write(args.out / "persistence_exposures.ecsv", overwrite=True)
    summ.write(args.out / "persistence_summary.ecsv", overwrite=True)
    summ.pprint_all()
    print(json.dumps({"n_files": len(files), "out": str(args.out)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
