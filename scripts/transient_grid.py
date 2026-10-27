"""Forced two-epoch photometry of every catalogued source in a field (D-027 grid).

The catalog stage of the two-epoch search (``transient_search.py`` then ``transient_combine.py``)
only measures positions that some catalog flags as changed in two bands. It misses a source that
changes in one module only (a blue transient seen in F090W and F115W but not in the LW bands) or
that the deblender hides. This script measures **every** position of the union of the given
pipeline catalogs (both epochs, all bands; duplicates within ``--dedup-arcsec`` merged) in every
band pair, with the aperture, background annulus, recentring and zero point of
``transient_forced.py`` (whose functions it reuses). It reads whole local ``_i2d.fits`` mosaics
(memory-mapped), since a field-wide grid touches every row of the image; byte-range cutouts would
read more than the file.

Noise: per band, the robust std of the flux-difference significance over all measured positions
is the noise scale (most sources do not vary; never below 1). A band **passes** when the
calibrated |significance| ≥ ``--min-sigma`` and |Δm| ≥ ``--min-dmag`` (or Δm is undefined and the
other epoch has positive flux: a source absent in one epoch). A position is a **hit** when
≥ ``--min-bands`` bands pass (with either sign by default: see Recentring below); hits are
grouped within ``--group-arcsec``. Bands run in parallel processes (``--procs``). Every number
is ``derived``; thresholds are ASSUMPTIONs.

Recentring (``--recentre-in origin``, the default): each position is centroided once, within
``--recentre-arcsec``, in the image of the catalog it came from (that catalog's band and epoch;
the band is read from the catalog's file name), and every band is measured at that sky position.
Recentring each band separately (``--recentre-in band``, the 2026-10-15 Sunrise and SMACS runs)
lets a band in which the source is absent pull the aperture onto a neighbour: MACS0416 SN01 was
missed that way. A hit also counts bands that pass with **opposite** signs (column ``mixed``; a
supernova that fades in the blue and brightens in the red); ``--same-sign`` restores the old rule.

``--inject N`` adds, in epoch 2 only, a circular Gaussian (FWHM ``--inject-fwhm-arcsec``) at N
random positions per band, with a flux equal to ``--inject-frac`` times the position's epoch-1
aperture flux (or ``--inject-floor`` × the band's median aperture error when that is larger), and
reports the fraction recovered as hits, overall and where the fraction (not the floor) set the
flux. Injections are measured in a separate pass with the real grid's zero point and noise scale.
Recovery counts the automatic thresholds only (no by-eye vetting). Every injection is written to
``injections.ecsv`` (position, fraction, injected flux and pass sign per band, recovered), so a post
hoc split needs no re-run.

``--coverage`` stores the number of exposures (bits set in the mosaic ``CON`` context extension) at
each position per band and epoch (``<band>_e1_n_exp``, ``<band>_e2_n_exp``) and splits the injection
recovery by the smallest epoch-1 count over bands (1 or >= 2).

    python scripts/transient_grid.py \\
        --band F277W e1_f277w_i2d.fits e2_f277w_i2d.fits --band F356W ... \\
        --catalog e1_f277w_cat.ecsv --catalog e2_f277w_cat.ecsv ... --epoch2-catalogs 1 \\
        --out outputs/grid
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import warnings
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_scales

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transient_forced import aperture_flux, compare, recentre, robust_std  # noqa: E402

MIN_SCALE_N = 100  # ASSUMPTION: positions needed per band for a stable noise scale


def first_of_groups(ra: np.ndarray, dec: np.ndarray, link_arcsec: float) -> np.ndarray:
    """Mask keeping the first row of any pair closer than ``link_arcsec``."""
    keep = np.ones(len(ra), bool)
    c = SkyCoord(ra * u.deg, dec * u.deg)
    i1, i2, _, _ = c.search_around_sky(c, link_arcsec * u.arcsec)
    for a, b in zip(i1, i2, strict=True):
        if a < b and keep[a]:
            keep[b] = False
    return keep


def union_positions(cats: list[Table], dedup_arcsec: float) -> Table:
    """RA/Dec of the union of pipeline catalogs, keeping the first of any group closer than
    ``dedup_arcsec``. ``origin`` is the index of the catalog a position came from."""
    rows = []
    for k, c in enumerate(cats):
        if "ra" in c.colnames:
            ra, dec = np.asarray(c["ra"], float), np.asarray(c["dec"], float)
        else:
            ra, dec = c["sky_centroid"].ra.deg, c["sky_centroid"].dec.deg
        ok = np.isfinite(ra) & np.isfinite(dec)
        rows.append(Table({"ra": ra[ok], "dec": dec[ok], "origin": np.full(ok.sum(), k)}))
    t = vstack(rows)
    out = t[first_of_groups(np.asarray(t["ra"]), np.asarray(t["dec"]), dedup_arcsec)]
    out["pos_id"] = [f"g{k:05d}" for k in range(len(out))]
    out.meta = {"provenance": "observed", "source": "union of pipeline catalogs"}
    return out


@contextmanager
def _open(path: str) -> Iterator[tuple[np.ndarray, np.ndarray, WCS, float]]:
    """SCI, ERR (memory-mapped), WCS and pixel scale (arcsec) of a level-3 mosaic."""
    with fits.open(path, memmap=True) as h:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(h["SCI"].header)
        scale = float(np.mean(proj_plane_pixel_scales(wcs))) * 3600.0
        yield h["SCI"].data, h["ERR"].data, wcs, scale


def _stamp(img: np.ndarray, x: float, y: float, half: int) -> tuple[np.ndarray, int, int] | None:
    ny, nx = img.shape
    if not (np.isfinite(x) and np.isfinite(y)):
        return None
    ix, iy = int(round(x)), int(round(y))
    if ix - half < 0 or iy - half < 0 or ix + half >= nx or iy + half >= ny:
        return None
    return (
        np.asarray(img[iy - half : iy + half + 1, ix - half : ix + half + 1], float),
        ix - half,
        iy - half,
    )


def n_exposures(path: str, ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """Exposures contributing at each position: bits set over every plane of the mosaic's ``CON``
    extension (32 exposures per plane); 0 off the image."""
    with fits.open(path, memmap=True) as h:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(h["SCI"].header)
        con = h["CON"].data
        if con.ndim == 2:
            con = con[None]
        x, y = wcs.world_to_pixel_values(ra, dec)
        out = np.zeros(len(ra), int)
        ok = np.isfinite(x) & np.isfinite(y)
        ix, iy = np.round(x[ok]).astype(int), np.round(y[ok]).astype(int)
        inside = (ix >= 0) & (iy >= 0) & (ix < con.shape[2]) & (iy < con.shape[1])
        k = np.flatnonzero(ok)[inside]
        vals = np.asarray(con[:, iy[inside], ix[inside]]).astype(np.uint32)  # planes x n
        bits = np.unpackbits(vals.T.copy().view(np.uint8), axis=1).sum(axis=1)
        out[k] = bits
    return out


def _count_job(job: tuple) -> tuple[str, str, np.ndarray]:
    band, tag, path, ra, dec = job
    return band, tag, n_exposures(path, ra, dec)


def gaussian(
    shape: tuple[int, int], x: float, y: float, sigma_px: float, flux: float
) -> np.ndarray:
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    g = np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / sigma_px**2)
    return flux * g / g.sum()


def measure_band(
    img1: tuple,
    img2: tuple,
    ra: np.ndarray,
    dec: np.ndarray,
    from_e2: np.ndarray,
    r_arcsec: float,
    recentre_arcsec: float | None,
    inject_flux: np.ndarray | None = None,
    inject_sigma_arcsec: float = 0.05,
) -> dict[str, np.ndarray]:
    """Aperture flux and error in both epochs at each position (recentred in the epoch the
    position came from, then measured at the same sky position in both; ``recentre_arcsec``
    None: measured at the given positions, already refined by ``refine_positions``). ``inject_flux``
    (per position, NaN = none) adds a Gaussian at the recentred position in epoch 2."""
    n = len(ra)
    f = {k: np.full(n, np.nan) for k in ("f1", "e1", "f2", "e2")}
    (s1, er1, w1, sc1), (s2, er2, w2, sc2) = img1, img2
    r1, r2 = r_arcsec / sc1, r_arcsec / sc2
    rca = recentre_arcsec or 0.0
    half1 = int(np.ceil(4.0 * r1)) + int(np.ceil(3 * rca / sc1)) + 2
    half2 = int(np.ceil(4.0 * r2)) + int(np.ceil(3 * rca / sc2)) + 2
    x1, y1 = w1.world_to_pixel_values(ra, dec)
    x2, y2 = w2.world_to_pixel_values(ra, dec)
    for k in range(n):
        a = _stamp(s1, x1[k], y1[k], half1)
        b = _stamp(s2, x2[k], y2[k], half2)
        if a is None or b is None:
            continue
        ea = _stamp(er1, x1[k], y1[k], half1)[0]
        eb = _stamp(er2, x2[k], y2[k], half2)[0]
        (sa, xa0, ya0), (sb, xb0, yb0) = a, b
        if recentre_arcsec is None:  # positions already refined (--recentre-in origin)
            rra, rdec = ra[k], dec[k]
        elif from_e2[k]:
            rc = max(1, int(np.ceil(recentre_arcsec / sc2)))
            cx, cy = recentre(sb, x2[k] - xb0, y2[k] - yb0, rc)
            rra, rdec = w2.pixel_to_world_values(cx + xb0, cy + yb0)
        else:
            rc = max(1, int(np.ceil(recentre_arcsec / sc1)))
            cx, cy = recentre(sa, x1[k] - xa0, y1[k] - ya0, rc)
            rra, rdec = w1.pixel_to_world_values(cx + xa0, cy + ya0)
        pa = np.array(w1.world_to_pixel_values(rra, rdec)) - (xa0, ya0)
        pb = np.array(w2.world_to_pixel_values(rra, rdec)) - (xb0, yb0)
        if inject_flux is not None and np.isfinite(inject_flux[k]):
            sb = sb + gaussian(sb.shape, pb[0], pb[1], inject_sigma_arcsec / sc2, inject_flux[k])
        f["f1"][k], f["e1"][k] = aperture_flux(sa, ea, pa[0], pa[1], r1, 2.5 * r1, 4.0 * r1)
        f["f2"][k], f["e2"][k] = aperture_flux(sb, eb, pb[0], pb[1], r2, 2.5 * r2, 4.0 * r2)
    return f


def catalog_bands(paths: list[str], bands: list[str]) -> list[str]:
    """The band of each catalog: the one band name that is a whole token of its file name
    (split at non-alphanumerics, case-insensitive), e.g. ``f115w`` in a pipeline ``_cat.ecsv``."""
    out = []
    for p in paths:
        tokens = set(re.split(r"[^0-9a-z]+", Path(p).name.lower()))
        found = [b for b in bands if b.lower() in tokens]
        if len(found) != 1:
            raise SystemExit(
                f"cannot tell the band of {p} from its name (matches: {found}); "
                "rename it or use --recentre-in band"
            )
        out.append(found[0])
    return out


def refine_positions(
    path: str, ra: np.ndarray, dec: np.ndarray, recentre_arcsec: float
) -> tuple[np.ndarray, np.ndarray]:
    """Recentre each position once in one image (``recentre``); positions off the image keep
    their input. Used by ``--recentre-in origin``: a source is centroided where its catalog saw
    it, so a band in which it is absent cannot drag the aperture onto a neighbour."""
    ra_o, dec_o = np.array(ra, float), np.array(dec, float)
    with _open(path) as (sci, _, w, sc):
        half = max(1, int(np.ceil(recentre_arcsec / sc)))
        stamp_half = 3 * half + 2
        x, y = w.world_to_pixel_values(ra, dec)
        for k in range(len(ra)):
            st = _stamp(sci, x[k], y[k], stamp_half)
            if st is None:
                continue
            a, x0, y0 = st
            cx, cy = recentre(a, x[k] - x0, y[k] - y0, half)
            ra_o[k], dec_o[k] = w.pixel_to_world_values(cx + x0, cy + y0)
    return ra_o, dec_o


def _refine_job(job: tuple) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path, idx, ra, dec, rc = job
    return (idx, *refine_positions(path, ra, dec, rc))


def band_pass(
    m: dict[str, np.ndarray],
    sys_floor: float,
    min_dmag: float,
    min_sigma: float,
    scale: float | None = None,
    zp: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    """``(dmag, calibrated flux significance, pass sign (+1 brighter, -1 fainter, 0), scale, zp)``.

    ``zp`` is rebuilt from ``compare``'s output (raw minus corrected dmag) so it can be passed
    back. Without ``scale``/``zp`` they come from these positions (``compare``'s clipped zero point,
    the robust std of the significances); pass the real grid's values when measuring injections,
    which would otherwise shift both."""
    f1 = m["f1"] if zp is None else m["f1"] * 10 ** (-0.4 * zp)
    e1 = m["e1"] if zp is None else m["e1"] * 10 ** (-0.4 * zp)
    dm, _, sig = compare(f1, e1, m["f2"], m["e2"], sys_floor, len(f1) + 1 if zp is not None else 10)
    if zp is None:
        with np.errstate(divide="ignore", invalid="ignore"):
            raw = -2.5 * np.log10(m["f2"] / m["f1"])
        d = raw - dm
        zp = float(np.nanmedian(d)) if np.isfinite(d).any() else 0.0
    if scale is None:
        s = robust_std(sig, MIN_SCALE_N)
        scale = max(1.0, s) if np.isfinite(s) else 1.0
    sig_c = sig / scale
    with np.errstate(invalid="ignore"):
        # an undefined dmag (flux <= 0 in an epoch) counts only when the other epoch has flux:
        # both negative means the annulus outshines the aperture (structure, not a source)
        seen = np.fmax(f1, m["f2"]) > 0
        ok = (np.abs(sig_c) >= min_sigma) & ((~np.isfinite(dm) & seen) | (np.abs(dm) >= min_dmag))
    sign = np.where(ok & np.isfinite(sig_c), np.sign(sig_c), 0).astype(int)
    return dm, sig_c, sign, scale, zp


def hits(signs: np.ndarray, min_bands: int, any_sign: bool = False) -> np.ndarray:
    """Rows where at least ``min_bands`` bands pass with the same sign (``signs``: rows × bands),
    or, with ``any_sign``, with either sign (a source whose colour changed, e.g. a supernova that
    fades in the blue and brightens in the red)."""
    same = ((signs > 0).sum(axis=1) >= min_bands) | ((signs < 0).sum(axis=1) >= min_bands)
    return (signs != 0).sum(axis=1) >= min_bands if any_sign else same


def group_hits(ra: np.ndarray, dec: np.ndarray, hit: np.ndarray, link_arcsec: float) -> np.ndarray:
    """Friends-of-friends group index of each hit within ``link_arcsec`` (-1 for non-hits)."""
    grp = np.full(len(ra), -1)
    idx = np.flatnonzero(hit)
    if not len(idx):
        return grp
    c = SkyCoord(ra[idx] * u.deg, dec[idx] * u.deg)
    i1, i2, _, _ = c.search_around_sky(c, link_arcsec * u.arcsec)
    lab = np.arange(len(idx))
    changed = True
    while changed:  # propagate the smallest label along every link
        new = lab.copy()
        np.minimum.at(new, i1, lab[i2])
        changed = bool((new != lab).any())
        lab = new
    _, grp[idx] = np.unique(lab, return_inverse=True)
    return grp


def _band_job(job: tuple) -> tuple[str, dict[str, np.ndarray], float]:
    """One band (pool worker): measure, optionally with injected fluxes; returns wall time."""
    band, p1, p2, ra, dec, from_e2, r, rc, flux, sigma_arcsec = job
    t0 = time.perf_counter()
    with _open(p1) as i1, _open(p2) as i2:
        m = measure_band(i1, i2, ra, dec, from_e2, r, rc, flux, sigma_arcsec)
    return band, m, time.perf_counter() - t0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--band", nargs=3, action="append", required=True, metavar=("BAND", "E1", "E2"))
    ap.add_argument("--catalog", action="append", required=True, help="pipeline _cat.ecsv")
    ap.add_argument(
        "--epoch2-catalogs",
        type=int,
        nargs="*",
        required=True,
        help="indices (0-based) of --catalog entries from epoch 2, recentred there (may be empty)",
    )
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--dedup-arcsec", type=float, default=0.2)
    ap.add_argument("--radius-arcsec", type=float, default=0.15)
    ap.add_argument("--recentre-arcsec", type=float, default=0.15)
    ap.add_argument(
        "--catalog-bands",
        nargs="+",
        help="band of each --catalog, in order (default: read from the file names)",
    )
    ap.add_argument(
        "--recentre-in",
        choices=["origin", "band"],
        default="origin",
        help="origin: centroid once in the catalog's own band and epoch, then measure every band "
        "there; band: centroid in each band separately (the 2026-10-15 Sunrise/SMACS runs)",
    )
    ap.add_argument("--min-dmag", type=float, default=0.3)
    ap.add_argument("--min-sigma", type=float, default=5.0)
    ap.add_argument("--min-bands", type=int, default=2)
    ap.add_argument(
        "--same-sign",
        action="store_true",
        help="require the passing bands to share a sign (the 2026-10-15 Sunrise/SMACS rule); "
        "by default a mixed-sign row is also a hit and is marked 'mixed'",
    )
    ap.add_argument("--sys-floor", type=float, default=0.05)
    ap.add_argument("--group-arcsec", type=float, default=1.5, help="hit grouping link length")
    ap.add_argument("--inject", type=int, default=0)
    ap.add_argument("--inject-frac", type=float, nargs="+", default=[1.0])
    ap.add_argument("--inject-floor", type=float, default=10.0)
    ap.add_argument("--inject-fwhm-arcsec", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--coverage",
        action="store_true",
        help="store per-position exposure counts (CON) and split injection recovery by them",
    )
    ap.add_argument("--procs", type=int, default=min(os.cpu_count() or 1, 8))
    args = ap.parse_args(argv)

    cats = [Table.read(c, format="ascii.ecsv") for c in args.catalog]
    pos = union_positions(cats, args.dedup_arcsec)
    from_e2 = np.isin(np.asarray(pos["origin"]), args.epoch2_catalogs)
    ra, dec = np.asarray(pos["ra"]), np.asarray(pos["dec"])
    args.out.mkdir(parents=True, exist_ok=True)
    r, rc = args.radius_arcsec, args.recentre_arcsec
    bands = [b for b, _, _ in args.band]
    signs, real, meta = [], {}, {"bands": {}}
    pool = ProcessPoolExecutor(max_workers=args.procs) if args.procs > 1 else None
    pmap = pool.map if pool is not None else map  # bands are independent
    try:
        if args.recentre_in == "origin" and rc > 0:
            cat_band = args.catalog_bands or catalog_bands(args.catalog, bands)
            if len(cat_band) != len(args.catalog) or not set(cat_band) <= set(bands):
                raise SystemExit("--catalog-bands needs one measured --band name per --catalog")
            paths = {b: (p1, p2) for b, p1, p2 in args.band}
            origin = np.asarray(pos["origin"])
            rjobs = []
            for k, b in enumerate(cat_band):
                idx = np.flatnonzero(origin == k)
                if len(idx):
                    path = paths[b][1 if k in args.epoch2_catalogs else 0]
                    rjobs.append((path, idx, ra[idx], dec[idx], rc))
            ra, dec = ra.copy(), dec.copy()
            for idx, rra, rdec in pmap(_refine_job, rjobs):
                ra[idx], dec[idx] = rra, rdec
            pos["ra_cat"], pos["dec_cat"] = pos["ra"], pos["dec"]
            pos["ra"], pos["dec"] = ra, dec
            # two catalog entries can centroid onto one peak: merge them again
            keep = first_of_groups(ra, dec, args.dedup_arcsec)
            pos, ra, dec, from_e2 = pos[keep], ra[keep], dec[keep], from_e2[keep]
            meta["n_merged_after_refinement"] = int((~keep).sum())
            rc = None  # every band measures at the refined position
        jobs = [(b, p1, p2, ra, dec, from_e2, r, rc, None, 0.0) for b, p1, p2 in args.band]
        done = {b: (m, dt) for b, m, dt in pmap(_band_job, jobs)}
        for band, p1, p2 in args.band:
            m, dt = done[band]
            dm, sig_c, sign, scale, zp = band_pass(m, args.sys_floor, args.min_dmag, args.min_sigma)
            pos[f"{band}_flux1"], pos[f"{band}_flux2"] = m["f1"], m["f2"]
            pos[f"{band}_dmag"], pos[f"{band}_sigma"], pos[f"{band}_pass"] = dm, sig_c, sign
            signs.append(sign)
            real[band] = (m, scale, zp)
            n_meas = int(np.isfinite(m["f1"] + m["f2"]).sum())
            meta["bands"][band] = {
                "e1": Path(p1).name,
                "e2": Path(p2).name,
                "noise_scale": scale,
                "zero_point_mag": zp,
                "n_measured": n_meas,
                "n_pass": int((sign != 0).sum()),
                "positions_per_s": len(pos) / dt,
            }
            print(
                f"{band}: {n_meas} measured, scale {scale:.2f}, {(sign != 0).sum()} pass, "
                f"{len(pos) / dt:.0f} positions/s",
                flush=True,
            )
        if args.coverage:
            cjobs = [
                (b, tag, p, ra, dec)
                for b, p1, p2 in args.band
                for tag, p in (("e1", p1), ("e2", p2))
            ]
            for b, tag, n in pmap(_count_job, cjobs):
                pos[f"{b}_{tag}_n_exp"] = n
        injections, inj_rows = [], []
        if args.inject:  # only where every band was measured, so recovery is per testable position
            every = np.all([np.isfinite(v[0]["f1"] + v[0]["f2"]) for v in real.values()], axis=0)
            rng = np.random.default_rng(args.seed)
            idx = np.flatnonzero(every)
            inj_idx = np.sort(rng.choice(idx, size=min(args.inject, len(idx)), replace=False))
            sigma = args.inject_fwhm_arcsec / 2.3548
            for frac in args.inject_frac:
                fluxes, jobs = {}, []
                for band, p1, p2 in args.band:
                    m = real[band][0]
                    f1 = np.fmax(m["f1"][inj_idx], 0)
                    fluxes[band] = (f1, frac * f1 >= args.inject_floor * np.nanmedian(m["e2"]))
                    flux = np.fmax(frac * f1, args.inject_floor * np.nanmedian(m["e2"]))
                    jobs.append(
                        (
                            band,
                            p1,
                            p2,
                            ra[inj_idx],
                            dec[inj_idx],
                            from_e2[inj_idx],
                            r,
                            rc,
                            flux,
                            sigma,
                        )
                    )
                si = {}
                for band, mi, _ in pmap(_band_job, jobs):
                    _, scale, zp = real[band]
                    si[band] = band_pass(
                        mi, args.sys_floor, args.min_dmag, args.min_sigma, scale, zp
                    )[2]
                si = np.column_stack([si[b] for b in bands])
                rec = (si > 0).sum(1) >= args.min_bands
                # the floor makes faint positions brighten by more than frac; report both subsets
                nominal = np.all([fluxes[b][1] for b in bands], axis=0)
                rows = Table(
                    {
                        "pos_id": pos["pos_id"][inj_idx],
                        "ra": ra[inj_idx],
                        "dec": dec[inj_idx],
                        "frac": np.full(len(inj_idx), frac),
                        "recovered": rec,
                        "nominal": nominal,
                    }
                )
                for j, job in enumerate(jobs):  # jobs follow --band order, like si's columns
                    rows[f"{job[0]}_inj_flux"] = job[8]
                    rows[f"{job[0]}_pass"] = si[:, j]
                by_cov = None
                if args.coverage:
                    e1n = np.min([pos[f"{b}_e1_n_exp"][inj_idx] for b in bands], axis=0)
                    rows["e1_n_exp_min"] = e1n
                    by_cov = {}
                    for key, sel in (("1", e1n == 1), (">=2", e1n >= 2)):
                        by_cov[key] = {
                            "n": int(sel.sum()),
                            "recovered_frac": float(rec[sel].mean()) if sel.any() else None,
                        }
                inj_rows.append(rows)
                injections.append(
                    {
                        "frac": frac,
                        "nominal_dmag": float(-2.5 * np.log10(1 + frac)),
                        "n": len(inj_idx),
                        "recovered_frac": float(rec.mean()),
                        "n_nominal": int(nominal.sum()),
                        "recovered_frac_nominal": float(rec[nominal].mean())
                        if nominal.any()
                        else None,
                        "per_band_pass_frac": {
                            b: float((si[:, j] > 0).mean()) for j, b in enumerate(bands)
                        },
                        **({"by_e1_n_exp_min": by_cov} if by_cov is not None else {}),
                    }
                )
    finally:
        if pool is not None:
            pool.shutdown()
    s = np.column_stack(signs)
    pos["n_pass_bright"], pos["n_pass_faint"] = (s > 0).sum(1), (s < 0).sum(1)
    pos["hit"] = hits(s, args.min_bands, any_sign=not args.same_sign)
    pos["mixed"] = pos["hit"] & ~hits(s, args.min_bands)
    pos["group"] = group_hits(ra, dec, np.asarray(pos["hit"]), args.group_arcsec)
    pos.meta.update(
        {
            "provenance": "derived",
            "thresholds_are": "ASSUMPTION",
            **meta,
            "min_sigma": args.min_sigma,
            "min_dmag": args.min_dmag,
            "min_bands": args.min_bands,
            "recentre_in": args.recentre_in,
            "same_sign": args.same_sign,
        }
    )
    pos.write(args.out / "grid.ecsv", format="ascii.ecsv", overwrite=True)
    pos[pos["hit"]].write(args.out / "hits.ecsv", format="ascii.ecsv", overwrite=True)
    if inj_rows:
        inj = vstack(inj_rows)
        inj.meta = {"provenance": "derived", "note": "pass: +1 brighter, -1 fainter, 0 none"}
        inj.write(args.out / "injections.ecsv", format="ascii.ecsv", overwrite=True)
    measured = [np.isfinite(pos[f"{b}_flux1"] + pos[f"{b}_flux2"]) for b in bands]
    summary = {
        "n_positions": len(pos),
        "n_merged_after_refinement": meta.get("n_merged_after_refinement", 0),
        "n_measured_any_band": int(np.any(measured, axis=0).sum()),
        "n_measured_all_bands": int(np.all(measured, axis=0).sum()),
        "n_hits": int(pos["hit"].sum()),
        "n_hit_groups": int(pos["group"].max() + 1),
        "n_hits_mixed": int(pos["mixed"].sum()),
        "recentre_in": args.recentre_in,
        "same_sign": args.same_sign,
        "bands": meta["bands"],
        "injection": {
            "floor_x_median_error": args.inject_floor,
            "fwhm_arcsec": args.inject_fwhm_arcsec,
            "runs": injections,
            "note": "recovery by the automatic thresholds only; by-eye vetting is not emulated; "
            "injections sit at the measured (recentred or refined) position, so recentring is not "
            "tested; a brightening injected in every band is recovered alike under the same-sign "
            "and mixed-sign rules",
        },
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
