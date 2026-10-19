"""Exposure-level check of dimming-screen forced-photometry points on level-2 ``_cal`` frames.

A point that appears in one band and one epoch of the resampled ``_i2d`` mosaic (and nowhere
else) is the signature of a cosmic ray, snowball or detector residual that survived outlier
rejection, not of a sky source. A sky source sits at the same sky position in every dither of that
epoch; a residual sits in one exposure (or in the exposures of one detector). This script measures,
for each (source, epoch) pair, a fixed-position aperture flux (no recentring, so an absent source
reads as zero) in every ``_cal`` exposure of that epoch and filter that covers the position, plus
the number of aperture pixels flagged JUMP_DET (DQ bit 4) or DO_NOT_USE (bit 1).

Fluxes are in the ``_cal`` units (MJy/sr) summed over pixels and converted to µJy with the pixel
area from the header (``PIXAR_SR``). Every number is ``derived``; the aperture radius is the
forced stage's (ASSUMPTION).

    python scripts/dimming_exposure_check.py --field nexus \\
        --table results/dimming/nexus_calib/survivors_and_flagged_controls.ecsv \\
        --pair nexus-d53823:0 --pair nexus-d79956:6 --out outputs/dimming/nexus_vet
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import yaml
from astropy.table import Table
from astropy.wcs import FITSFixedWarning

sys.path.insert(0, str(Path(__file__).resolve().parent))
from persistence_check import _open, header_info, list_cal  # noqa: E402
from transient_forced import aperture_flux  # noqa: E402

CONFIG = Path(__file__).resolve().parents[1] / "configs" / "dimming_screen.yaml"
APER_ARCSEC = 0.15  # ASSUMPTION: the forced stage's aperture (transient_forced default)
EXTRA_APER_ARCSEC = (
    0.3,
    0.5,
)  # aperture-size dependence (extended-nucleus test); columns flux_30, flux_50
SCALE = {"F200W": 0.031, "F444W": 0.063}  # arcsec per pixel (SW / LW)
JUMP_DET, DO_NOT_USE = 4, 1
HALF = 30  # stamp half-size in pixels (fits the 0.5" aperture's annulus on SW pixels)


# ASSUMPTION: no astrophysical spectrum gives F_nu(F200W) / F_nu(F444W) above this; the
# Rayleigh-Jeans limit is (4.44 / 1.99)^2 = 5, and the factor 4 margin covers strong emission
# lines in F200W
MAX_SW_LW_RATIO = 20.0
NSIG = 3.0


def chromatic_veto(
    f_sw: float, f_lw: float, e_lw: float, max_ratio: float = MAX_SW_LW_RATIO
) -> bool:
    """True when the simultaneous F444W flux is too low for any sky source with this F200W flux:
    ``f_lw + NSIG * e_lw < f_sw / max_ratio``. F200W and F444W are taken in the same exposure
    (dichroic), so a sky source, lensed or not, appears in both; a lensing spike is achromatic.
    False when any input is not finite or ``f_sw <= 0`` (the test does not apply)."""
    if not all(np.isfinite(v) for v in (f_sw, f_lw, e_lw)) or f_sw <= 0:
        return False
    return bool(f_lw + NSIG * e_lw < f_sw / max_ratio)


def cal_dirs(obs: str) -> list[str]:
    """Every visit directory of an observation: ``jw05105-o001_t001`` ->
    ``s3://stpubdata/jwst/public/jw05105/jw05105001VVV`` for each visit VVV on S3."""
    import fsspec

    m = re.match(r"jw(\d{5})-o(\d{3})_", obs)
    if not m:
        raise ValueError(f"unrecognised observation id {obs!r}")
    prog, o = m.groups()
    fs = fsspec.filesystem("s3", anon=True)
    pre = f"jw{prog}{o}"
    names = [d.rsplit("/", 1)[-1] for d in fs.ls(f"stpubdata/jwst/public/jw{prog}", detail=False)]
    return [
        f"s3://stpubdata/jwst/public/jw{prog}/{n}"
        for n in sorted(names)
        if re.fullmatch(pre + r"\d{3}", n)
    ]


def measure(uri: str, x: float, y: float, band: str) -> dict | None:
    """Fixed-position aperture flux (µJy), its error, flagged-pixel counts and a stamp."""
    r_px = APER_ARCSEC / SCALE[band]
    with _open(uri) as hdul, warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        hdr = hdul["SCI"].header
        ny, nx = hdr["NAXIS2"], hdr["NAXIS1"]
        xi, yi = int(round(x)), int(round(y))
        if not (HALF <= xi < nx - HALF and HALF <= yi < ny - HALF):
            return None
        sl = np.s_[yi - HALF : yi + HALF + 1, xi - HALF : xi + HALF + 1]
        sci = np.asarray(hdul["SCI"].section[sl], dtype=float)
        err = np.asarray(hdul["ERR"].section[sl], dtype=float)
        dq = np.asarray(hdul["DQ"].section[sl]).astype(np.int64)
        pixar = float(hdr.get("PIXAR_SR", hdul[0].header.get("PIXAR_SR", np.nan)))
    lx, ly = x - (xi - HALF), y - (yi - HALF)
    yy, xx = np.mgrid[0 : sci.shape[0], 0 : sci.shape[1]]
    ap = np.hypot(xx - lx, yy - ly) <= r_px
    n_jump = int(((dq[ap] & JUMP_DET) > 0).sum())
    n_dnu = int(((dq[ap] & DO_NOT_USE) > 0).sum())
    clean = sci.copy()
    clean[(dq & DO_NOT_USE) > 0] = np.nan
    f, e = aperture_flux(clean, err, lx, ly, r_px, r_px + 3, r_px + 7)
    to_ujy = pixar * 1e12  # MJy/sr * sr -> MJy; * 1e12 -> µJy
    extra = {}
    for r in EXTRA_APER_ARCSEC:
        rp = r / SCALE[band]
        extra[f"flux_{round(r * 100):02d}"] = (
            aperture_flux(clean, None, lx, ly, rp, rp + 3, rp + 7)[0] * to_ujy
        )
    return {
        **extra,
        "flux_ujy": f * to_ujy,
        "err_ujy": e * to_ujy,
        "flux_raw_ujy": float(np.nansum(sci[ap] - np.nanmedian(sci[~ap]))) * to_ujy,
        "n_ap": int(ap.sum()),
        "n_jump": n_jump,
        "n_dnu": n_dnu,
        "stamp": sci,
    }


def sheet(stamps: dict, tab: Table, path: Path, ncol: int = 6) -> None:
    """One panel per (source, epoch, band, exposure) stamp, titled with flux and JUMP count."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keys = sorted(stamps)
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.2 * ncol, 2.4 * nrow), squeeze=False)
    for a in axes.flat:
        a.axis("off")
    for a, key in zip(axes.flat, keys, strict=False):
        uid, k, band, name = key
        s = stamps[key]
        lo, hi = np.nanpercentile(s, [5, 99.5])
        a.imshow(s, origin="lower", cmap="gray", vmin=lo, vmax=max(hi, lo + 1e-6))
        r = tab[
            (tab["uid"] == uid)
            & (tab["epoch"] == k)
            & (tab["band"] == band)
            & (tab["exposure"] == name)
        ][0]
        a.set_title(
            f"{uid.replace('nexus-', '')} e{k} {band}\n{name[-16:]} "
            f"f={r['flux_ujy']:.2f} J={r['n_jump']}",
            fontsize=6,
        )
    fig.tight_layout()
    fig.savefig(path, dpi=60)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--field", required=True)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--pair", action="append", required=True, help="uid:epoch_index")
    ap.add_argument("--bands", nargs="+", default=["F200W", "F444W"])
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--hit-band", default="F200W")
    ap.add_argument("--lw-band", default="F444W")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    epochs = yaml.safe_load(CONFIG.read_text())["fields"][args.field]["epochs"]
    tab = Table.read(args.table)
    pos = {str(r["uid"]): (float(r["ra"]), float(r["dec"])) for r in tab}
    pairs = [(p.rsplit(":", 1)[0], int(p.rsplit(":", 1)[1])) for p in args.pair]
    by_epoch: dict[int, dict[str, tuple[float, float]]] = {}
    for uid, k in pairs:
        by_epoch.setdefault(k, {})[uid] = pos[uid]

    rows = []
    stamps = {}
    with ThreadPoolExecutor(args.workers) as pool:
        for k, targets in sorted(by_epoch.items()):
            files = [f for d in cal_dirs(epochs[k]["obs"]) for f in list_cal(d) if "_nrc" in f]
            heads = [h for h in pool.map(header_info, files, [targets] * len(files)) if h]
            jobs = [
                (h, uid, xy)
                for h in heads
                if h["filter"] in args.bands
                for uid, xy in h["pix"].items()
            ]
            res = pool.map(lambda j: measure(j[0]["file"], *j[2], j[0]["filter"]), jobs)
            for (h, uid, _), m in zip(jobs, res, strict=True):
                if m is None:
                    continue
                name = Path(h["file"]).name.replace("_cal.fits", "")
                stamps[(uid, k, h["filter"], name)] = m.pop("stamp")
                rows.append(
                    {
                        "uid": uid,
                        "epoch": k,
                        "band": h["filter"],
                        "exposure": name,
                        "detector": h["detector"],
                        "mjd": h["mjd"],
                        **m,
                    }
                )
    args.out.mkdir(parents=True, exist_ok=True)
    out = Table(rows=rows) if rows else Table()
    out.write(args.out / "exposure_fluxes.ecsv", overwrite=True)
    summary = {}
    rows_by_uid = {str(r["uid"]): r for r in tab}
    for uid, k in pairs:
        r = rows_by_uid[uid]
        cols = [
            f"{args.hit_band}_e{k}_fflux",
            f"{args.lw_band}_e{k}_fflux",
            f"{args.lw_band}_e{k}_ferr",
        ]
        if all(c in tab.colnames for c in cols):
            f_sw, f_lw, e_lw = (float(r[c]) for c in cols)
            summary[f"{uid}|e{k}|forced"] = {
                "f_sw": f_sw,
                "f_lw": f_lw,
                "e_lw": e_lw,
                "chromatic_veto": chromatic_veto(f_sw, f_lw, e_lw),
            }
        for band in args.bands:
            sel = [r for r in rows if r["uid"] == uid and r["epoch"] == k and r["band"] == band]
            if not sel:
                continue
            f = np.array([r["flux_ujy"] for r in sel])
            summary[f"{uid}|e{k}|{band}"] = {
                "n_exposures": len(sel),
                "flux_ujy": [round(float(v), 4) for v in f],
                "n_jump": [r["n_jump"] for r in sel],
                "median_ujy": float(np.nanmedian(f)),
                "max_ujy": float(np.nanmax(f)) if np.isfinite(f).any() else None,
            }
    (args.out / "exposure_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    if stamps:
        sheet(stamps, out, args.out / "exposure_sheet.png")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
