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
import gzip
import hashlib
import io
import json
import os
import sys
import time
import warnings
from pathlib import Path

# TrailFinder runs a process pool sized to the cores: one BLAS thread each (scripts/CLAUDE.md)
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.wcs import WCS  # noqa: E402
from astropy.wcs.utils import proj_plane_pixel_scales  # noqa: E402
from scipy import ndimage  # noqa: E402

BLOCK = 4  # block median of native pixels (ASSUMPTION)
REBIN = 2  # further mean rebin before the trail finder (ASSUMPTION)
CLIP_NSIG = 10.0  # normalised pixels clipped to ±this (ASSUMPTION; PLAN addendum A)
EDGE_ERODE = 12  # footprint eroded by this many block pixels, outside set to 0 (ASSUMPTION; A)
MRT_THRESHOLD = 5.0  # TrailFinder threshold (ASSUMPTION; other settings its defaults)
KEEP_STATUS = 2  # TrailFinder status of a persistent trail (ASSUMPTION)
VETO_ARCSEC = 1.0  # shape rows within this of a masked pixel are vetoed (ASSUMPTION)
ROWS_PER_READ = 256  # native rows per range read (multiple of BLOCK)
# row/column-segment stripe veto (results/w1_row_veto/PLAN.md; every number an ASSUMPTION)
SEG_MIN_VALID = 0.5  # finite share a segment needs
SEG_DETREND = 41  # running-median window along the segment column (block pixels)
SEG_NSIG = 5.0  # |residual| > this × 1.4826 MAD is a stripe
KEEP_FRAC = 0.75  # cost check: rows kept with the veto ≥ this × rows kept without it


def i2d_uri(catalog: str | Path) -> str:
    """The level-3 F200W mosaic behind a MAST pipeline ``_cat.ecsv`` (public S3 bucket)."""
    stem = Path(str(catalog)).name.removesuffix("_cat.ecsv")
    prog, obs = stem.split("-")[0], stem.split("-")[1].split("_")[0]
    return f"s3://stpubdata/jwst/public/{prog}/L3/t/{obs}/{stem}_i2d.fits"


def block_median(a: np.ndarray, b: int) -> np.ndarray:
    """``b`` × ``b`` block median; zero (off-mosaic fill) and NaN are ignored."""
    ny, nx = (a.shape[0] // b) * b, (a.shape[1] // b) * b
    a = np.where(a[:ny, :nx] == 0, np.nan, a[:ny, :nx].astype(np.float32))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN (off-mosaic) blocks
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


def stream_binned(uri: str, workers: int = 8) -> tuple[np.ndarray, WCS]:
    """Block-median SCI image and its WCS, read in row strips (no full copy kept); ``workers``
    concurrent range reads (the PRIMER SCI extensions are 2-3 GB)."""
    from concurrent.futures import ThreadPoolExecutor

    import s3fs

    fs = s3fs.S3FileSystem(anon=True)
    path = uri.removeprefix("s3://")
    with fs.open(path, "rb", block_size=2**20) as f, fits.open(f, lazy_load_hdus=True) as h:
        hdu = h["SCI"]
        hd = hdu.header
        if hd["BITPIX"] != -32:
            raise ValueError(f"{uri}: SCI BITPIX {hd['BITPIX']}, expected -32")
        nx, ny = hd["NAXIS1"], hd["NAXIS2"]
        start = hdu.fileinfo()["datLoc"]
        wcs = WCS(hd)
    out = np.full((ny // BLOCK, nx // BLOCK), np.nan, np.float32)

    def strip(y0: int) -> None:
        nr = min(ROWS_PER_READ, (ny // BLOCK) * BLOCK - y0)
        b = fs.cat_file(path, start=start + y0 * nx * 4, end=start + (y0 + nr) * nx * 4)
        a = np.frombuffer(b, ">f4").reshape(nr, nx)
        out[y0 // BLOCK : (y0 + nr) // BLOCK] = block_median(a, BLOCK)

    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(strip, range(0, (ny // BLOCK) * BLOCK, ROWS_PER_READ)))
    return out, wcs[::BLOCK, ::BLOCK]


def tile_slices(wcs: WCS, shape, centre, half_arcsec: float) -> tuple[slice, slice]:
    """Pixel slices (y, x) of the square tile of half-width ``half_arcsec`` around ``centre`` on an
    image with ``wcs`` (clipped to the image; results/w1_shear_primer/PLAN.md). The tile's rows
    are the ones ``inject_radial.in_box`` keeps; a one-pixel margin covers the rounding."""
    scale = proj_plane_pixel_scales(wcs)[0] * 3600.0
    x, y = wcs.world_to_pixel_values(*centre)
    h = half_arcsec / scale + 1
    ny, nx = shape
    ys = slice(max(0, int(np.floor(y - h))), min(ny, int(np.ceil(y + h)) + 1))
    xs = slice(max(0, int(np.floor(x - h))), min(nx, int(np.ceil(x + h)) + 1))
    return ys, xs


def _row_stripes(img: np.ndarray, seg_len: int) -> np.ndarray:
    """Row half of :func:`segment_stripes`: residual of each row segment's median against a
    running median over ``SEG_DETREND`` rows of the same segment column (NaN if undefined)."""
    ny, nx = img.shape
    ns = nx // seg_len
    seg = img[:, : ns * seg_len].reshape(ny, ns, seg_len)
    valid = np.isfinite(seg).mean(axis=2) >= SEG_MIN_VALID
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        val = np.where(valid, np.nanmedian(seg, axis=2), np.nan)
        trend = ndimage.generic_filter(val, np.nanmedian, size=(SEG_DETREND, 1), mode="nearest")
    return val - trend


def segment_stripes(img: np.ndarray, seg_len: int) -> np.ndarray:
    """Stripe pixels of a block-median image (PLAN steps 1–4): row segments of ``seg_len``
    pixels (last partial one dropped) with ≥ ``SEG_MIN_VALID`` finite pixels, median-valued,
    detrended along y, flagged if |residual| > ``SEG_NSIG`` σ (σ from the MAD of all residuals
    of that axis); the same on columns. Returns a bool image of ``img``'s shape."""
    img = np.asarray(img, float)
    out = np.zeros(img.shape, bool)
    for axis_img, put in ((img, lambda m: m), (img.T, lambda m: m.T)):
        res = _row_stripes(axis_img, seg_len)
        ok = np.isfinite(res)
        if not ok.any():
            continue
        sig = 1.4826 * np.median(np.abs(res[ok] - np.median(res[ok])))
        if not sig > 0:
            continue
        flag = np.zeros(res.shape, bool)
        flag[ok] = np.abs(res[ok]) > SEG_NSIG * sig
        full = np.zeros(axis_img.shape, bool)
        full[:, : flag.shape[1] * seg_len] = np.repeat(flag, seg_len, axis=1)
        out |= put(full)
    return out


def to_mask_grid(stripes: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """OR of each ``REBIN`` × ``REBIN`` block onto the trail mask's grid of ``shape``."""
    ny, nx = shape
    s = stripes[: ny * REBIN, : nx * REBIN]
    return s.reshape(ny, REBIN, nx, REBIN).any(axis=(1, 3))


def kept_rows(summary: dict) -> int:
    """Rows the screen used after the streak veto (PLAN cost check)."""
    sc = summary["screen"]
    return int(sc["n_used_before_veto"] - sc["n_veto_used"])


def cost_ok(new: dict, old: dict, frac: float = KEEP_FRAC) -> bool:
    """Cost check: the union-mask screen keeps ≥ ``frac`` of the rows the D-097 screen kept."""
    return kept_rows(new) >= frac * kept_rows(old)


def find_trails(z: np.ndarray, processes: int = 4, threshold: float = MRT_THRESHOLD):
    """Persistent-trail mask and trail table from ``acstools.findsat_mrt``."""
    from acstools.findsat_mrt import TrailFinder

    tf = TrailFinder(
        z,
        processes=processes,
        threshold=threshold,
        save_catalog=False,
        save_diagnostic=False,
        mask_include_status=[KEEP_STATUS],
    )
    tf.run_all()
    cat = tf.source_list
    mask = np.zeros(z.shape, bool) if tf.mask is None else np.asarray(tf.mask, bool)
    return mask, cat


def write_mask(path: Path, mask: np.ndarray, wcs: WCS) -> None:
    """Gzipped FITS with a zero gzip timestamp, so an identical mask gives an identical file."""
    buf = io.BytesIO()
    fits.PrimaryHDU(mask.astype(np.uint8), header=wcs.to_header()).writeto(buf)
    with open(path, "wb") as f, gzip.GzipFile(filename="", fileobj=f, mode="wb", mtime=0) as g:
        g.write(buf.getvalue())


def read_mask(path: Path) -> tuple[np.ndarray, WCS]:
    with fits.open(path) as h:
        return np.asarray(h[0].data, bool), WCS(h[0].header)


def veto_function(mask: np.ndarray, wcs: WCS, veto_arcsec: float = VETO_ARCSEC):
    """``(ra, dec) -> bool``: True within ``veto_arcsec`` of a masked pixel, measured from the
    nearest mask-pixel centre (so the tolerance is exact only to about half a mask pixel, 0.18″
    at 0.25″ pixels); positions off the mask image are False. The distance map is built once."""
    scale = proj_plane_pixel_scales(wcs)[0] * 3600.0
    dist = ndimage.distance_transform_edt(~mask) * scale if mask.any() else None

    def veto(ra, dec) -> np.ndarray:
        ra, dec = np.atleast_1d(np.asarray(ra, float)), np.atleast_1d(np.asarray(dec, float))
        out = np.zeros(len(ra), bool)
        if dist is None:
            return out
        x, y = wcs.world_to_pixel_values(ra, dec)
        ix, iy = np.round(np.asarray(x)).astype(int), np.round(np.asarray(y)).astype(int)
        inside = (ix >= 0) & (iy >= 0) & (ix < mask.shape[1]) & (iy < mask.shape[0])
        out[inside] = dist[iy[inside], ix[inside]] <= veto_arcsec
        return out

    return veto


def veto_rows(ra, dec, mask: np.ndarray, wcs: WCS, veto_arcsec: float = VETO_ARCSEC) -> np.ndarray:
    """One-off form of :func:`veto_function`."""
    return veto_function(mask, wcs, veto_arcsec)(ra, dec)


def mask_identity(path: Path) -> dict:
    """sha256 of the file and of the mask pixels (shape + uint8 bytes), and the masked-pixel
    count, recorded in each veto-run summary. Files written before 2026-10-17T14:30Z carry a gzip
    timestamp, so only ``data_sha256`` and ``mask_pixels`` reproduce for them."""
    mask, _ = read_mask(path)
    data = hashlib.sha256(str(mask.shape).encode() + mask.astype(np.uint8).tobytes())
    return {
        "file": path.name,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "data_sha256": data.hexdigest(),
        "mask_pixels": int(mask.sum()),
    }


def known_case(mask_covers_hit: bool, macs1423_passes_with_veto: bool) -> dict:
    """K1 and K2 of the plan."""
    return {"K1": bool(mask_covers_hit), "K2": not macs1423_passes_with_veto}


def adopt(k: dict) -> bool:
    """Frozen adoption rule: the veto is adopted only if K1 and K2 both hold."""
    return bool(k["K1"] and k["K2"])


# the 2026-10-17 D-053 hit
HIT_SUMMARY = (
    Path(__file__).resolve().parents[1] / "results/w1_shear_more/inject_summary_macs1423.json"
)


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
    ap.add_argument("--processes", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--png", action="store_true", help="also write a contact sheet per field")
    ap.add_argument(
        "--mrt-threshold",
        type=float,
        default=MRT_THRESHOLD,
        help="TrailFinder threshold (results/w1_streak_rowdiag/PLAN.md: faint-streak variants)",
    )
    ap.add_argument(
        "--seg-len",
        type=int,
        default=None,
        help="add the row/column-segment stripe mask with this segment length in block pixels "
        "(results/w1_row_veto/PLAN.md); default: trail mask only",
    )
    ap.add_argument(
        "--evaluate",
        type=Path,
        metavar="SUMMARY",
        help="MACS1423 inject_shear summary run with the veto: write <out>/adoption.json",
    )
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.evaluate is not None:
        if args.fields:
            ap.error("--evaluate makes no masks; run --fields separately")
        res = evaluate(args.out, args.evaluate, HIT_SUMMARY)
        (args.out / "adoption.json").write_text(json.dumps(res, indent=1) + "\n")
        print(json.dumps(res))
        return 0
    mosaics = {}  # PRIMER tiles: one streamed and normalised mosaic per i2d, cut per tile
    for name in args.fields:
        t0 = time.time()
        spec = ish.field_spec(name)
        uri = i2d_uri(spec["catalog"])
        if uri not in mosaics:
            mosaics.clear()  # one mosaic in memory at a time (fields come grouped by mosaic)
            img, wcs = stream_binned(uri)
            mosaics[uri] = (img, wcs, normalise(img), wcs[::REBIN, ::REBIN])
        img, wcs, z, wz = mosaics[uri]
        t_read = time.time() - t0
        tile = None
        if spec.get("box_half") is not None:  # normalised over the whole mosaic, then cut
            ys, xs = tile_slices(wz, z.shape, spec["centre"], spec["box_half"])
            z, wz = z[ys, xs], wz[ys, xs]
            img = img[ys.start * REBIN : ys.stop * REBIN, xs.start * REBIN : xs.stop * REBIN]
            tile = {"y": [ys.start, ys.stop], "x": [xs.start, xs.stop]}
        mask, cat = find_trails(z, args.processes, args.mrt_threshold)
        n_trail_px = int(mask.sum())
        seg = None
        if args.seg_len is not None:
            stripes = segment_stripes(img, args.seg_len)
            mask = mask | to_mask_grid(stripes, mask.shape)
            seg = {
                "seg_len": args.seg_len,
                "min_valid": SEG_MIN_VALID,
                "detrend_rows": SEG_DETREND,
                "nsig": SEG_NSIG,
                "stripe_block_px": int(stripes.sum()),
                "trail_mask_pixels": n_trail_px,
            }
        write_mask(args.out / f"{name}_streak_mask.fits.gz", mask, wz)
        trails = []
        if cat is not None:
            for r in cat[cat["status"] == KEEP_STATUS]:
                trails.append({"theta_deg": float(r["theta"]), "rho_px": float(r["rho"])})
        info = {
            "field": name,
            "i2d": uri,
            "binned_shape": list(z.shape),
            "tile_slice_rebinned_px": tile,
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
                "mrt_threshold": args.mrt_threshold,
                "keep_status": KEEP_STATUS,
                "veto_arcsec": VETO_ARCSEC,
            },
            "segment_stripes": seg,
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
