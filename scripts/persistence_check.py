"""Detector-persistence test for transient candidates on level-2 ``_cal`` exposures (D-039).

A source that appears in one epoch only can be an afterimage: a bright source illuminates a
detector pixel in one exposure, and the pixel still shows ~0.05–3 % of that signal in the next
exposures (Sunrise ``n0022`` and ``n0150``, docs/fields/sunrise.md). The level-3 mosaic maps
that pixel to a different sky position after each dither, so the afterimage looks like a faint,
compact source at a fixed or wandering sky position.

For each position this script lists every ``_cal`` exposure that covers it (headers read from S3
or local paths, no full downloads), measures the aperture flux there, and measures the same
*detector pixel* in every earlier exposure on the same detector within ``--lookback-hours``.
Earlier exposures that put the position itself on (nearly) the same pixel are skipped: the source
would be its own "illuminator". A detection (S/N ≥ ``--min-snr``) is ``suspect`` when an earlier
exposure put ≥ ``--ratio`` × its flux, or any saturated pixel, on that pixel. It is ``clean`` when
it is not suspect and at least one earlier exposure was checked (``n_prior`` > 0; the first
exposure of a visit cannot be cleared). Per position:

- ``persistence``: detections, none of them clean, at least one suspect;
- ``on_sky``: at least two clean detections;
- ``inconclusive``: otherwise; ``undetected``: no detection.

Only imaging exposures (``EXP_TYPE`` ``*_IMAGE`` with a celestial WCS) are used. Apertures that
leave the array give NaN. A file that cannot be read is listed in ``meta["failed"]``.

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
from collections.abc import Iterator, Sequence
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

from jwst_anomaly.cutouts import _band_from_header, _is_remote, _open_image, _resolve_uri

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transient_forced import aperture_flux  # noqa: E402

RATIO = 20.0  # ASSUMPTION: earlier flux ≥ 20x the detection (observed afterimages: 0.05–3 %)
MIN_SNR = 5.0  # ASSUMPTION: ERR-based S/N of a detection
LOOKBACK_HOURS = 3.0  # ASSUMPTION: NIRCam afterimages above noise fade within hours
SATURATED = 2  # JWST DQ bit
HALF = 15  # cutout half-size, pixels
_PAD = 5  # accept positions up to this many pixels outside the array (cutout still useful)
SAME_PIXEL = 2.0  # earlier exposures with the position within this many pixels are skipped
BLOCK = 2**20


def classify(
    snr: np.ndarray, suspect: np.ndarray, n_prior: np.ndarray, min_snr: float = MIN_SNR
) -> str:
    """Verdict for one position from its per-exposure S/N, suspect flags and prior counts."""
    snr = np.asarray(snr, dtype=float)
    sus = np.asarray(suspect, dtype=bool)
    det = np.isfinite(snr) & (snr >= min_snr)
    if not det.any():
        return "undetected"
    clean = int((det & ~sus & (np.asarray(n_prior) > 0)).sum())
    if clean >= 2:
        return "on_sky"
    if clean == 0 and (det & sus).any():
        return "persistence"
    return "inconclusive"


def is_suspect(flux: float, prior_flux: float, prior_sat: int, ratio: float = RATIO) -> bool:
    """True when an earlier exposure put far more signal (or saturation) on this pixel."""
    if prior_sat > 0:
        return True
    return bool(
        np.isfinite(flux) and flux > 0 and np.isfinite(prior_flux) and prior_flux >= ratio * flux
    )


def list_cal(cal_dir: str) -> list[str]:
    import fsspec

    opts = {"anon": True} if cal_dir.startswith("s3://") else {}
    fs, path = fsspec.core.url_to_fs(cal_dir, **opts)
    prefix = "s3://" if cal_dir.startswith("s3://") else ""
    return sorted(prefix + p for p in fs.ls(path, detail=False) if p.endswith("_cal.fits"))


@contextmanager
def _open(uri: str) -> Iterator[fits.HDUList]:
    """Remote: ``cutouts._open_image`` (fsspec read-ahead). Local: no memmap, because the
    scaled uint32 DQ extension of ``_cal`` files cannot be memory-mapped; ``.section`` still
    reads only the requested rows."""
    uri = _resolve_uri(uri)
    if _is_remote(uri):
        with _open_image(uri, None, BLOCK) as (hdul, _):
            yield hdul
    else:
        with fits.open(uri, memmap=False, lazy_load_hdus=True) as hdul:
            yield hdul


def header_info(uri: str, positions: dict[str, tuple[float, float]]) -> dict[str, Any] | None:
    """Detector, filter, start time and the pixel of each position on this exposure.

    None for files that are not 2-D imaging exposures (``EXP_TYPE`` not ``*_IMAGE``, no
    ``EXPSTART``, a non-2-D ``SCI`` or no celestial WCS)."""
    with (
        _open(uri) as hdul,
        warnings.catch_warnings(),
    ):
        warnings.simplefilter("ignore", FITSFixedWarning)
        primary, sci = hdul[0].header, hdul["SCI"].header
        exp_type = str(primary.get("EXP_TYPE", "")).strip().upper()
        if not exp_type.endswith("_IMAGE") or "EXPSTART" not in primary or sci.get("NAXIS") != 2:
            return None
        wcs = WCS(sci, naxis=2)
        if not wcs.has_celestial:
            return None
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
        band = _band_from_header(primary, sci, uri)
    return {
        "file": uri,
        "detector": str(primary.get("DETECTOR", "")).strip(),
        "filter": band,
        "mjd": float(primary["EXPSTART"]),
        "pix": pix,
    }


def measure_pixel(uri: str, x: float, y: float, long: bool) -> dict[str, float]:
    """Aperture flux, ERR-based error and saturated-pixel count at detector pixel (x, y).

    Flux and error are NaN when the aperture leaves the array."""
    xi, yi = int(round(x)), int(round(y))
    r_px = 2.0 if long else 3.0  # ~0.1" in both channels
    with _open(uri) as hdul:
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
    flux, e = aperture_flux(sci, err, x - x0, y - y0, r_px, 2 * r_px + 1, 4 * r_px)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    core = np.hypot(xx - x, yy - y) <= r_px
    if not (r_px <= x <= nx - 1 - r_px and r_px <= y <= ny - 1 - r_px):
        flux = e = float("nan")
    return {"flux": flux, "err": e, "n_sat": int(((dq[core] & SATURATED) > 0).sum())}


def _safe(fn, *args):
    try:
        return fn(*args), None
    except Exception as exc:  # noqa: BLE001 - one bad file must not lose the run
        return None, f"{args[0]}: {type(exc).__name__}: {exc}"


def run(
    positions: dict[str, tuple[float, float]],
    cal_files: Sequence[str],
    ratio: float = RATIO,
    min_snr: float = MIN_SNR,
    lookback_hours: float = LOOKBACK_HOURS,
    workers: int = 12,
) -> tuple[Table, Table]:
    failed: list[str] = []
    with ThreadPoolExecutor(workers) as ex:
        heads = list(ex.map(lambda f: _safe(header_info, f, positions), cal_files))
    info = [h for h, _ in heads if h is not None]
    failed += [err for _, err in heads if err]
    info.sort(key=lambda r: r["mjd"])
    jobs = []
    for i, cur in enumerate(info):
        for uid, (x, y) in cur["pix"].items():
            prior = []
            for p in info[:i]:
                if p["detector"] != cur["detector"]:
                    continue
                if not 0 < cur["mjd"] - p["mjd"] <= lookback_hours / 24:
                    continue
                own = p["pix"].get(uid)
                if own is not None and np.hypot(own[0] - x, own[1] - y) < SAME_PIXEL:
                    continue  # the position itself sat on this pixel then
                prior.append(p)
            jobs.append((uid, cur, prior, x, y))

    def one(job):
        uid, cur, prior, x, y = job
        long = "LONG" in cur["detector"].upper()
        now, err = _safe(measure_pixel, cur["file"], x, y, long)
        if now is None:
            failed.append(err)
            now = {"flux": np.nan, "err": np.nan, "n_sat": 0}
        snr = now["flux"] / now["err"] if now["err"] > 0 else np.nan
        before = []
        if np.isfinite(snr) and snr >= min_snr:  # priors only matter for detections
            for p in prior:
                m, err = _safe(measure_pixel, p["file"], x, y, long)
                if m is None:
                    failed.append(err)
                else:
                    before.append((p, m))
        worst = max(
            before,
            key=lambda b: (b[1]["n_sat"], np.nan_to_num(b[1]["flux"], nan=-np.inf)),
            default=None,
        )
        prior_flux = worst[1]["flux"] if worst else np.nan
        prior_sat = sum(b[1]["n_sat"] for b in before)
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
            "n_prior": len(before),
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
        else Table(names=["uid", "snr", "suspect", "n_prior"], dtype=[str, float, bool, int])
    )
    summary = []
    for uid in positions:
        sel = per_exp[per_exp["uid"] == uid] if len(per_exp) else per_exp
        snr = np.asarray(sel["snr"], dtype=float)
        sus = np.asarray(sel["suspect"], dtype=bool)
        npr = np.asarray(sel["n_prior"], dtype=int)
        det = np.isfinite(snr) & (snr >= min_snr)
        summary.append(
            {
                "uid": uid,
                "n_exposures": len(sel),
                "n_detected": int(det.sum()),
                "n_suspect": int((det & sus).sum()),
                "n_clean": int((det & ~sus & (npr > 0)).sum()),
                "verdict": classify(snr, sus, npr, min_snr),
            }
        )
    meta = {
        "provenance": "derived",
        "source": "level-2 _cal exposures: "
        + ", ".join(sorted({f.rsplit("/", 1)[0] for f in cal_files})),
        "ratio": ratio,
        "min_snr": min_snr,
        "lookback_hours": lookback_hours,
        "failed": failed,
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
