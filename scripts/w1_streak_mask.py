"""Field-wide linear-streak masks for the W1 shear screen (frozen rule:
results/w1_streak_veto/PLAN.md with addendum A).

    python scripts/w1_streak_mask.py --fields macs1423 abell2744 ... \
        --out results/w1_streak_veto/masks

For each field, streams the SCI extension of the F200W ``_i2d.fits`` behind the screen's
catalogue (S3, not stored), block-medians it, masks sources and runs STScI's median-Radon-transform
trail finder (``acstools.findsat_mrt``; Stark et al. 2022). Trails that pass its persistence check
(status 2) form the mask, written as ``<field>_streak_mask.fits.gz`` (uint8 with the binned WCS)
plus a per-field JSON. ``inject_shear.py --streak-masks`` drops shape rows within ``VETO_ARCSEC`` of
the mask. Pixels ``observed``; masks ``derived``; every number here is an ASSUMPTION.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy import ndimage

BLOCK = 4  # block median of native pixels (ASSUMPTION)
REBIN = 2  # further mean rebin before the trail finder (ASSUMPTION)
CLIP_NSIG = 10.0  # normalised pixels clipped to ±this (ASSUMPTION; PLAN addendum A)
EDGE_ERODE = 12  # footprint eroded by this many block pixels, outside set to 0 (ASSUMPTION; A)
MRT_THRESHOLD = 5.0  # TrailFinder threshold (ASSUMPTION; other settings its defaults)
KEEP_STATUS = 2  # TrailFinder status of a persistent trail (ASSUMPTION)
VETO_ARCSEC = 1.0  # shape rows within this of a masked pixel are vetoed (ASSUMPTION)
ROWS_PER_READ = 256  # native rows per range read (multiple of BLOCK)


def i2d_uri(catalog: str | Path) -> str:
    """The level-3 F200W mosaic behind a MAST pipeline ``_cat.ecsv`` (public S3 bucket)."""
    stem = Path(str(catalog)).name.removesuffix("_cat.ecsv")
    prog, obs = stem.split("-")[0], stem.split("-")[1].split("_")[0]
    return f"s3://stpubdata/jwst/public/{prog}/L3/t/{obs}/{stem}_i2d.fits"


def block_median(a: np.ndarray, b: int) -> np.ndarray:
    """``b`` × ``b`` block median; zero (off-mosaic fill) and NaN are ignored."""
    ny, nx = (a.shape[0] // b) * b, (a.shape[1] // b) * b
    a = np.where(a[:ny, :nx] == 0, np.nan, a[:ny, :nx].astype(np.float32))
    with np.errstate(all="ignore"):
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmedian(a.reshape(ny // b, b, nx // b, b), axis=(1, 3))


def normalise(d: np.ndarray) -> np.ndarray:
    """Background-normalised image clipped to ±``CLIP_NSIG``, zero outside the eroded footprint,
    mean-rebinned. Sources are not removed: the median Radon transform is robust to compact
    sources, and zeroing bright pixels also removed bright streaks (PLAN addendum A)."""
    ok = np.isfinite(d)
    med = np.nanmedian(d)
    mad = 1.4826 * np.nanmedian(np.abs(d[ok] - med))
    z = np.zeros_like(d, dtype=float)
    if mad > 0:
        z[ok] = np.clip((d[ok] - med) / mad, -CLIP_NSIG, CLIP_NSIG)
    z[~ndimage.binary_erosion(ok, iterations=EDGE_ERODE)] = 0.0
    ny, nx = (z.shape[0] // REBIN) * REBIN, (z.shape[1] // REBIN) * REBIN
    return z[:ny, :nx].reshape(ny // REBIN, REBIN, nx // REBIN, REBIN).mean(axis=(1, 3))


def stream_binned(uri: str) -> tuple[np.ndarray, WCS]:
    """Block-median SCI image and its WCS, read in row strips (no full copy kept)."""
    import s3fs

    fs = s3fs.S3FileSystem(anon=True)
    path = uri.removeprefix("s3://")
    with fs.open(path, "rb", block_size=2**24) as f, fits.open(f, lazy_load_hdus=True) as h:
        hdu = h["SCI"]
        hd = hdu.header
        if hd["BITPIX"] != -32:
            raise ValueError(f"{uri}: SCI BITPIX {hd['BITPIX']}, expected -32")
        nx, ny = hd["NAXIS1"], hd["NAXIS2"]
        start = hdu.fileinfo()["datLoc"]
        out = np.full((ny // BLOCK, nx // BLOCK), np.nan, np.float32)
        for y0 in range(0, (ny // BLOCK) * BLOCK, ROWS_PER_READ):
            nr = min(ROWS_PER_READ, (ny // BLOCK) * BLOCK - y0)
            f.seek(start + y0 * nx * 4)
            a = np.frombuffer(f.read(nr * nx * 4), ">f4").reshape(nr, nx)
            out[y0 // BLOCK : (y0 + nr) // BLOCK] = block_median(a, BLOCK)
        wcs = WCS(hd)
    return out, wcs[::BLOCK, ::BLOCK]


def find_trails(z: np.ndarray, processes: int = 4):
    """Persistent-trail mask and trail table from ``acstools.findsat_mrt``."""
    from acstools.findsat_mrt import TrailFinder

    tf = TrailFinder(
        z,
        processes=processes,
        threshold=MRT_THRESHOLD,
        save_catalog=False,
        save_diagnostic=False,
        mask_include_status=[KEEP_STATUS],
    )
    tf.run_all()
    cat = tf.source_list
    mask = np.zeros(z.shape, bool) if tf.mask is None else np.asarray(tf.mask, bool)
    return mask, cat


def write_mask(path: Path, mask: np.ndarray, wcs: WCS) -> None:
    fits.PrimaryHDU(mask.astype(np.uint8), header=wcs.to_header()).writeto(path, overwrite=True)


def read_mask(path: Path) -> tuple[np.ndarray, WCS]:
    with fits.open(path) as h:
        return np.asarray(h[0].data, bool), WCS(h[0].header)


def veto_rows(ra, dec, mask: np.ndarray, wcs: WCS, veto_arcsec: float = VETO_ARCSEC) -> np.ndarray:
    """True within ``veto_arcsec`` of a masked pixel (positions off the mask image: False)."""
    ra, dec = np.atleast_1d(np.asarray(ra, float)), np.atleast_1d(np.asarray(dec, float))
    if not mask.any():
        return np.zeros(len(ra), bool)
    scale = proj_plane_pixel_scales(wcs)[0] * 3600.0
    dist = ndimage.distance_transform_edt(~mask) * scale
    x, y = wcs.world_to_pixel_values(ra, dec)
    ix, iy = np.round(np.asarray(x)).astype(int), np.round(np.asarray(y)).astype(int)
    inside = (ix >= 0) & (iy >= 0) & (ix < mask.shape[1]) & (iy < mask.shape[0])
    out = np.zeros(len(ra), bool)
    out[inside] = dist[iy[inside], ix[inside]] <= veto_arcsec
    return out


def known_case(mask_covers_hit: bool, macs1423_passes_with_veto: bool) -> dict:
    """K1 and K2 of the plan."""
    return {"K1": bool(mask_covers_hit), "K2": not macs1423_passes_with_veto}


def adopt(k: dict) -> bool:
    """Frozen adoption rule: the veto is adopted only if K1 and K2 both hold."""
    return bool(k["K1"] and k["K2"])


HIT_SUMMARY = "results/w1_shear_more/inject_summary_macs1423.json"  # the 2026-10-17 D-053 hit


def evaluate(mask_dir: Path, veto_summary: Path, hit_summary: Path) -> dict:
    """K1, K2 and the adoption decision (PLAN) from the MACS1423 mask and veto-run summary."""
    hit = json.loads(hit_summary.read_text())["screen"]
    mask, wcs = read_mask(mask_dir / "macs1423_streak_mask.fits.gz")
    covers = bool(veto_rows([hit["s_max_radec"][0]], [hit["s_max_radec"][1]], mask, wcs)[0])
    new = json.loads(veto_summary.read_text())
    if new["field"] != "macs1423" or "streak_veto" not in new["assumptions"]:
        raise ValueError(f"{veto_summary} is not a MACS1423 summary run with the streak veto")
    k = known_case(covers, bool(new["screen"]["passes_eb_rule"]))
    return {
        "hit_radec": hit["s_max_radec"],
        "mask_covers_hit": covers,
        "macs1423_with_veto": {
            key: new["screen"][key] for key in ("s_max", "p_random_max", "b_max", "passes_eb_rule")
        },
        **k,
        "adopted": adopt(k),
        "rule": "results/w1_streak_veto/PLAN.md (K1, K2, adoption rule)",
    }


def main(argv: list[str] | None = None) -> int:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import inject_shear as ish

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fields", nargs="+", default=[], choices=ish.FIELD_NAMES)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--processes", type=int, default=4)
    ap.add_argument("--png", action="store_true", help="also write a contact sheet per field")
    ap.add_argument(
        "--evaluate",
        type=Path,
        metavar="SUMMARY",
        help="MACS1423 inject_shear summary run with the veto: write <out>/adoption.json",
    )
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.evaluate is not None:
        res = evaluate(args.out, args.evaluate, Path(HIT_SUMMARY))
        (args.out / "adoption.json").write_text(json.dumps(res, indent=1) + "\n")
        print(json.dumps(res))
        return 0
    for name in args.fields:
        t0 = time.time()
        uri = i2d_uri(ish.field_spec(name)["catalog"])
        img, wcs = stream_binned(uri)
        t_read = time.time() - t0
        z = normalise(img)
        wz = wcs[::REBIN, ::REBIN]
        mask, cat = find_trails(z, args.processes)
        write_mask(args.out / f"{name}_streak_mask.fits.gz", mask, wz)
        trails = []
        if cat is not None:
            for r in cat[cat["status"] == KEEP_STATUS]:
                trails.append({"theta_deg": float(r["theta"]), "rho_px": float(r["rho"])})
        info = {
            "field": name,
            "i2d": uri,
            "binned_shape": list(z.shape),
            "pixel_scale_arcsec": float(proj_plane_pixel_scales(wz)[0] * 3600.0),
            "n_trails_status2": len(trails),
            "n_candidates_all_status": 0 if cat is None else len(cat),
            "trails": trails,
            "mask_pixels": int(mask.sum()),
            "settings": {
                "block": BLOCK,
                "rebin": REBIN,
                "clip_nsig": CLIP_NSIG,
                "edge_erode": EDGE_ERODE,
                "mrt_threshold": MRT_THRESHOLD,
                "keep_status": KEEP_STATUS,
                "veto_arcsec": VETO_ARCSEC,
            },
            "wall_time_s": {"read": round(t_read, 1), "total": round(time.time() - t0, 1)},
            "provenance": {"pixels": "observed", "mask": "derived", "settings": "assumption"},
        }
        (args.out / f"{name}_streak_mask.json").write_text(json.dumps(info, indent=1) + "\n")
        if args.png:
            _sheet(args.out / f"{name}_streak_mask.png", img, mask, name)
        print(name, info["n_trails_status2"], info["mask_pixels"], info["wall_time_s"], flush=True)
    return 0


def _sheet(path: Path, img: np.ndarray, mask: np.ndarray, name: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    v = np.nanpercentile(img, [5, 99])
    fig, ax = plt.subplots(figsize=(10, 10 * img.shape[0] / img.shape[1] + 0.6))
    ax.imshow(img, origin="lower", cmap="gray", vmin=v[0], vmax=v[0] + 0.1 * (v[1] - v[0]))
    big = np.kron(mask, np.ones((REBIN, REBIN), bool))
    ov = np.zeros((*img.shape, 4))
    ov[: big.shape[0], : big.shape[1]][big] = (1, 0, 0, 0.35)
    ax.imshow(ov, origin="lower")
    ax.set_title(f"{name}: F200W (block {BLOCK}) with the status-2 trail mask (red)")
    fig.savefig(path, dpi=60, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
