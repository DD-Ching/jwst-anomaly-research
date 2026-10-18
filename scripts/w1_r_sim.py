"""Shear responsivity R of pipeline-like catalogue shapes from GalSim image simulations.

Plan, checks C1–C3 and the frozen rule: ``results/w1_r_metacal/PLAN.md``. Isolated Sérsic galaxies
are sheared by g = ±0.02 (same noise for both signs, 90°-rotated pairs), smoothed by the jwst
detection kernel, segmented at 3σ and measured with photutils ``SourceCatalog`` on the convolved
data (``semimajor_sigma``, ``ellipticity``, ``orientation``), then PSF-corrected by
``exotic_screens.intrinsic_ellipticity``. R per cell = (⟨ε(+g)⟩ − ⟨ε(−g)⟩) / 2g (``simulated``).
``--weight-field`` weights the cells by the used rows of a real field (nearest cell in
log semimajor σ and log S/N of each cell's kept galaxies, measured in their +g1 copy).

    python scripts/w1_r_sim.py --out results/w1_r_metacal/r_sim.json \
        --weight-field goodsn-o001 --work-dir <scratch>/goodsn [--raw-threshold]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing
import os
import resource
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exotic_screens as es  # noqa: E402
from w1_psf_smooth import KERNEL_SIGMA_PX  # noqa: E402

PIX_ARCSEC = 0.031  # NIRCam SW mosaic pixel, ASSUMPTION (as w1_psf_smooth)
PSF_FWHM_ARCSEC = 0.066  # F200W-like, ASSUMPTION (PLAN.md)
STAMP = 64
G_STEP = 0.02
SIGMA_E = 0.25  # intrinsic ellipticity per component, ASSUMPTION
HLR_PX = (1.5, 2.5, 4.0, 6.0, 10.0)
SNR = (10.0, 20.0, 50.0, 200.0)
SERSIC_N = (1.0, 4.0)


def psf(psf_fwhm=PSF_FWHM_ARCSEC):
    import galsim

    if psf_fwhm <= 0:
        return None
    # Airy (6.5 m, 2 µm) ⊗ Gaussian of FWHM psf_fwhm / 2, ASSUMPTION; the default's FWHM is 0.071″
    # (GalSim calculateFWHM, not stored), slightly broader than the nominal F200W 0.066″
    airy = galsim.Airy(lam=2000.0, diam=6.5, scale_unit=galsim.arcsec)
    return galsim.Convolve(airy, galsim.Gaussian(fwhm=0.5 * psf_fwhm))


def draw(gal, the_psf, noise_sigma: float, rng, kernel_sigma=KERNEL_SIGMA_PX):
    """Image of ``gal`` ⊗ PSF plus white noise from ``rng``, and its kernel-smoothed copy."""
    import galsim
    from scipy import ndimage

    obj = galsim.Convolve(gal, the_psf) if the_psf is not None else gal
    im = obj.drawImage(nx=STAMP, ny=STAMP, scale=PIX_ARCSEC).array.astype(float)
    if noise_sigma > 0:
        im = im + rng.normal(0.0, noise_sigma, im.shape)
    conv = ndimage.gaussian_filter(im, kernel_sigma, mode="constant") if kernel_sigma > 0 else im
    return im, conv


def measure(
    im, conv, noise_sigma: float, kernel_sigma=KERNEL_SIGMA_PX, nsig=3.0, npixels=5, raw_thr=False
):
    """(semimajor σ, minor σ, PA in deg in :func:`to_eps`'s frame, S/N) of the central source, or
    None. The threshold is ``nsig`` times the smoothed noise (the plan), or with ``raw_thr`` times
    the unsmoothed noise (jwst ``source_catalog``'s convention; post hoc sensitivity, a higher
    isophote). Noiseless stamps: 1e-3 of the peak."""
    from photutils.segmentation import SourceCatalog, detect_sources

    if noise_sigma > 0:
        # white noise smoothed by a unit-sum Gaussian of σ s: σ / (2 √π s)
        thr = (
            nsig * noise_sigma / (2 * np.sqrt(np.pi) * kernel_sigma)
            if kernel_sigma > 0 and not raw_thr
            else nsig * noise_sigma
        )
    else:
        thr = 1e-3 * conv.max()
    segm = detect_sources(conv, thr, n_pixels=npixels)
    if segm is None:
        return None
    lab = segm.data[STAMP // 2, STAMP // 2]
    if lab == 0:
        return None
    cat = SourceCatalog(im, segm, convolved_data=conv)
    src = cat.select_label(lab)
    # photutils 3 renamed semimajor_sigma (the jwst catalogue's column) to semimajor_axis
    a, b = float(src.semimajor_axis.value), float(src.semiminor_axis.value)
    area = float(src.area.value)
    snr = float(src.segment_flux) / (noise_sigma * np.sqrt(area)) if noise_sigma > 0 else np.inf
    # photutils orientation θ is from +x counter-clockwise; PA = θ maps GalSim's +g1 (stretch along
    # x) to +Re ε and +g2 to +Im ε in intrinsic_ellipticity's frame (C3)
    return a, b, float(src.orientation.value), snr


def to_eps(a, b, pa, psf_sigma):
    return complex(
        es.intrinsic_ellipticity(np.array([a]), np.array([b]), np.array([pa]), psf_sigma)[0]
    )


def psf_sigma_px(the_psf, kernel_sigma=KERNEL_SIGMA_PX) -> float:
    """Minor-axis σ of a noiseless PSF-only stamp, measured like the galaxies (C2's input)."""
    import galsim

    if the_psf is None and kernel_sigma <= 0:
        return 0.0  # nothing to subtract (C1)
    im, conv = draw(galsim.Gaussian(sigma=1e-4, flux=1.0), the_psf, 0.0, None, kernel_sigma)
    m = measure(im, conv, 0.0, kernel_sigma)
    if m is None:
        raise ValueError("the PSF-only stamp has no central source")
    return m[1]


def cell(job) -> dict:
    """R and measured (semimajor σ, S/N) medians of one (n, hlr, snr) cell."""
    import galsim

    n, hlr, snr, npair, seed, psf_fwhm, kernel_sigma, extra_psig, raw_thr = job
    rng = np.random.default_rng(seed)
    the_psf = psf(psf_fwhm)
    psigs = {"stamp": psf_sigma_px(the_psf, kernel_sigma), **extra_psig}
    base = galsim.Sersic(n=n, half_light_radius=hlr * PIX_ARCSEC, flux=1.0)
    # noise from the target S/N of the noiseless galaxy inside its 3σ-like isophote: start from
    # flux / (σ √A) with A = π (2 hlr)², ASSUMPTION; the measured S/N is what the weights use
    area0 = np.pi * (2.0 * max(hlr, 1.0)) ** 2
    noise = 1.0 / (snr * np.sqrt(area0)) if np.isfinite(snr) else 0.0
    resp, sizes, snrs = {k: [] for k in psigs}, [], []
    for _ in range(npair):
        e = rng.normal(0, SIGMA_E, 2)
        while np.hypot(*e) >= 0.8:
            e = rng.normal(0, SIGMA_E, 2)
        e = complex(*e)
        pair = []  # both 90°-rotated copies (rotation flips ε), so intrinsic shape noise cancels
        for rot in (1.0, -1.0):
            gal0 = base.shear(e1=rot * e.real, e2=rot * e.imag)
            nseed = int(rng.integers(2**31))
            moms = {}
            for comp in (0, 1):
                for sgn in (1.0, -1.0):
                    g = sgn * G_STEP
                    gal = gal0.shear(g1=g, g2=0.0) if comp == 0 else gal0.shear(g1=0.0, g2=g)
                    im, conv = draw(gal, the_psf, noise, np.random.default_rng(nseed), kernel_sigma)
                    m = measure(im, conv, noise, kernel_sigma, raw_thr=raw_thr)
                    if m is not None:
                        moms[comp, sgn] = m
            pair.append(moms)
        if any(len(m) < 4 for m in pair):
            continue  # a pair counts only if all eight sheared copies are detected at the centre
        for moms in pair:  # cell placement uses the same kept pairs (their +g1 copies)
            sizes.append(moms[0, 1.0][0])
            snrs.append(moms[0, 1.0][3])
        for k, ps in psigs.items():
            rs = []
            for moms in pair:
                eps = {key: to_eps(m[0], m[1], m[2], ps) for key, m in moms.items()}
                if not all(np.isfinite(x) for x in eps.values()):
                    break  # unresolved after PSF subtraction, as in the screen: drop the whole pair
                r1 = (eps[0, 1.0] - eps[0, -1.0]).real / (2 * G_STEP)
                r2 = (eps[1, 1.0] - eps[1, -1.0]).imag / (2 * G_STEP)
                rs.append(0.5 * (r1 + r2))
            if len(rs) == 2:
                resp[k].append(0.5 * (rs[0] + rs[1]))
    out = {
        "sersic_n": n,
        "hlr_px": hlr,
        "snr_target": snr,
        "psf_sigma_px": psigs,
        "n_pairs": npair,
        "n_detected": len(sizes),  # galaxies of pairs with all eight sheared copies detected
        "semimajor_px_median": float(np.median(sizes)) if sizes else None,
        "snr_median": float(np.median(snrs)) if snrs else None,
        "R": {},
    }
    for k, v in resp.items():
        v = np.asarray(v)
        out["R"][k] = {
            "R": float(np.mean(v)) if len(v) else None,
            "err": float(np.std(v) / np.sqrt(len(v))) if len(v) > 1 else None,
            "n_pairs_used": len(v),
        }
    return out


def used_rows(field: str, work_dir: Path):
    """(semimajor_px, snr, PSF σ) of a blank field's rows that pass the shear screen's S/N ≥ 10,
    resolved and spike-veto cuts (``shear_sources``; the screen's streak mask and PSF χ are not
    applied, so these approximate its used rows). Gaia star files go to ``work_dir``."""
    import inject_radial as ir
    import inject_shear as ish

    spec = ish.field_spec(field)
    model, shapes, extra, _cat, _pz = ir.load_field(field, work_dir, spec)
    psig = es.psf_sigma_px(shapes)
    e, _ = es.shear_sources(model, shapes, psig, 10.0, responsivity=ish.GOODSN_R, extra=extra)
    ok = np.isfinite(e)
    return np.asarray(shapes["semimajor_px"], float)[ok], np.asarray(shapes["snr"], float)[ok], psig


def weighted_r(cells: list[dict], size, snr, key: str, min_used_frac: float = 0.0) -> dict:
    """Mean cell R (PSF-σ choice ``key``) over real rows, each mapped to the nearest detected cell
    in (log semimajor σ, log S/N) of the cells' measured medians. Rows whose cell has fewer than
    two used pairs, or a used-pair fraction below ``min_used_frac``, are left out and counted as
    ``n_rows_unrepresented`` (R there is unknown, not assumed). Cells without any kept galaxy
    have no position, so rows near them go to a neighbour. ``n_rows_extrapolated`` counts rows
    beyond the placed cells' size or S/N range (they take an edge cell's R)."""
    det = [c for c in cells if c["semimajor_px_median"]]
    cs = np.log([c["semimajor_px_median"] for c in det])
    cn = np.log([max(c["snr_median"], 1e-3) for c in det])
    d = (np.log(size)[:, None] - cs[None]) ** 2 + (np.log(snr)[:, None] - cn[None]) ** 2
    counts = np.bincount(np.argmin(d, axis=1), minlength=len(det))
    outside = (
        (np.log(size) > cs.max())
        | (np.log(size) < cs.min())
        | (np.log(snr) > cn.max())
        | (np.log(snr) < cn.min())
    )
    ok = np.array(
        [
            c["R"][key]["n_pairs_used"] >= 2
            and c["R"][key]["n_pairs_used"] >= min_used_frac * c["n_pairs"]
            for c in det
        ]
    )
    r = np.array([c["R"][key]["R"] if o else 0.0 for c, o in zip(det, ok, strict=True)])
    err = np.array([c["R"][key]["err"] if o else 0.0 for c, o in zip(det, ok, strict=True)])
    w = np.where(ok, counts, 0) / max(int(np.where(ok, counts, 0).sum()), 1)
    return {
        "R": float(np.dot(w, r)),
        "err": float(np.sqrt(np.dot(w**2, err**2))),
        "n_rows_used": int(np.where(ok, counts, 0).sum()),
        "n_rows_unrepresented": int(np.where(ok, 0, counts).sum()),
        "n_rows_extrapolated": int(outside.sum()),
        "rows_per_cell": counts.tolist(),
    }


def rule_outcome(r_sim: float, err: float, r_borrowed: float = 0.414) -> str:
    """The frozen rule (PLAN.md): "conservative" when R_sim >= r_borrowed - 2 err (the borrowed R
    does not overstate the response; limits stand), else "optimistic" (re-run with R_sim)."""
    return "conservative" if r_sim >= r_borrowed - 2.0 * err else "optimistic"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--npair", type=int, default=100, help="intrinsic-ellipticity pairs per cell")
    p.add_argument(
        "--weight-field", default=None, help="blank field whose used rows weight the cells"
    )
    p.add_argument("--workers", type=int, default=os.cpu_count())
    p.add_argument(
        "--work-dir", type=Path, default=None, help="Gaia star files of --weight-field (required)"
    )
    p.add_argument("--seed", type=int, default=20261018)
    p.add_argument(
        "--raw-threshold",
        action="store_true",
        help="3σ of the unsmoothed noise on the smoothed image (jwst convention; post hoc)",
    )
    a = p.parse_args(argv)
    if a.weight_field and a.work_dir is None:
        p.error("--weight-field needs --work-dir")
    extra = {}
    if a.weight_field:
        size, snr, psig_real = used_rows(a.weight_field, a.work_dir)
        extra = {"field": psig_real}  # deviation (post hoc): the PSF σ the screen itself subtracts
    jobs = [
        (
            n,
            h,
            s,
            a.npair,
            a.seed + 1000 * i,
            PSF_FWHM_ARCSEC,
            KERNEL_SIGMA_PX,
            extra,
            a.raw_threshold,
        )
        for i, (n, h, s) in enumerate((n, h, s) for n in SERSIC_N for h in HLR_PX for s in SNR)
    ]
    c1_job = (0.5, 8.0, np.inf, 10, a.seed - 1, 0.0, 0.0, {}, False)
    os.environ["OMP_NUM_THREADS"] = "1"  # spawned workers inherit it: no BLAS oversubscription
    t0 = time.perf_counter()
    ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(a.workers, mp_context=ctx) as ex:
        c1_fut = ex.submit(cell, c1_job)
        cells = list(ex.map(cell, jobs))
        c1 = c1_fut.result()
    wall = time.perf_counter() - t0
    cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu_pct = 100 * (cpu.ru_utime + cpu.ru_stime) / (wall * a.workers)
    print(
        f"{len(jobs) + 1} cells in {wall:.0f} s ({(len(jobs) + 1) / wall:.2f} cells/s), "
        f"CPU {cpu_pct:.0f} % of {a.workers} cores",
        file=sys.stderr,
    )
    c1_r = c1["R"]["stamp"]["R"]
    out = {
        "label": "simulated",
        "plan": "results/w1_r_metacal/PLAN.md",
        "g_step": G_STEP,
        "threshold": "raw (post hoc)" if a.raw_threshold else "smoothed (plan)",
        # deviation: the plan names FWHM 0.066″; the Airy ⊗ Gaussian model measures this
        "psf_fwhm_arcsec_measured": float(
            psf().drawImage(nx=64, ny=64, scale=PIX_ARCSEC / 8, method="no_pixel").calculateFWHM()
        ),
        "check_c1": {"R": c1_r, "passed": c1_r is not None and abs(c1_r - 1) < 0.03},
        "cells": cells,
    }
    if a.weight_field:
        import inject_shear as ish

        out["weight_field"] = {
            "field": a.weight_field,
            "n_rows": len(size),
            "psf_sigma_px": psig_real,
        }
        out["R_sim"] = {}
        for key in ("stamp", "field"):
            for n in SERSIC_N:
                sub = [c for c in cells if c["sersic_n"] == n]
                for tag, frac in (("", 0.0), ("_resolved90", 0.9)):  # _resolved90: post hoc
                    res = weighted_r(sub, size, snr, key, frac)
                    res["rule"] = rule_outcome(res["R"], res["err"], ish.GOODSN_R)
                    out["R_sim"][f"psf_{key}_n{int(n)}{tag}"] = res
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "cells"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
