"""Does detection-kernel, jitter and drizzle rounding bring STPSF's χ to the three PSF stars?

Plan and rule S3: ``results/w1_psf_smooth/PLAN.md``. Same 24 STPSF F200W PSFs as
``w1_stpsf_floor.py`` (``DET_DIST``, 10″); each is smoothed per variant at detector sampling,
then the segment χ is measured on the fine isophote grid and matched to each star's
``isophotal_area``.
Star χ and areas come from ``results/w1_stpsf/stars_check.json`` (derived from the catalogue).

Usage::

    STPSF_PATH=<stpsf-data> python scripts/w1_psf_smooth.py \
        --out results/w1_psf_smooth/smooth_chi.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w1_stpsf_floor as wsf  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
PIX_ARCSEC = 0.031  # NIRCam SW detector and mosaic pixel, ASSUMPTION (as w1_stpsf_floor)
# jwst source_catalog/tweakreg_catalog: morphology on data convolved with a Gaussian of FWHM
# kernel_fwhm (default 2.0 px; ASSUMPTION that the Abell 2744 catalogue used the default).
KERNEL_SIGMA_PX = 2.0 / (2 * np.sqrt(2 * np.log(2)))
BOX_SIGMA_PX = 1 / np.sqrt(12)  # one-pixel drizzle box as its second-moment Gaussian, ASSUMPTION
S3_MAX_DIFF = 0.05  # per component, ASSUMPTION (PLAN.md)
# name → extra Gaussian σ components (px), added in quadrature; STPSF already includes 0.007″ jitter
VARIANTS = {
    "none": (),
    "kernel": (KERNEL_SIGMA_PX,),
    "kernel+box": (KERNEL_SIGMA_PX, BOX_SIGMA_PX),
    "kernel+box+jit5": (KERNEL_SIGMA_PX, BOX_SIGMA_PX, 0.005 / PIX_ARCSEC),
    "kernel+box+jit10": (KERNEL_SIGMA_PX, BOX_SIGMA_PX, 0.010 / PIX_ARCSEC),
}


def smooth(im: np.ndarray, sigmas) -> np.ndarray:
    """Convolve with one Gaussian whose σ² is the sum of ``sigmas``² (Gaussians compose)."""
    from scipy import ndimage

    s = float(np.sqrt(np.sum(np.square(sigmas)))) if len(sigmas) else 0.0
    return ndimage.gaussian_filter(im, s, mode="constant") if s > 0 else im


def _one(job) -> dict:
    import stpsf

    det, x, y = job
    nc = stpsf.NIRCam()
    nc.filter, nc.detector, nc.detector_position = wsf.FILTER, det, (x, y)
    im = nc.calc_psf(fov_arcsec=wsf.FOV_ARCSEC, oversample=2)["DET_DIST"].data.astype(float)
    out = {"detector": det, "x": x, "y": y, "variants": {}}
    for name, sig in VARIANTS.items():
        fine = [wsf.isophotal_chi_area(smooth(im, sig), t) for t in wsf.LEVELS_FINE]
        out["variants"][name] = {
            "fine_area_px": [a for _, a, _ in fine],
            "fine_chi_pix": [[c.real, c.imag] for c, _, _ in fine],
            "fine_edge": [e for _, _, e in fine],
        }
    return out


def compute(workers: int) -> list[dict]:
    from concurrent.futures import ProcessPoolExecutor

    jobs = [(d, x, y) for d in wsf.DETECTORS for x, y in wsf.POSITIONS]
    with ProcessPoolExecutor(workers, initializer=wsf._init_worker) as ex:
        return list(ex.map(_one, jobs))


def s3_star(star_chi, model_chi) -> bool:
    d = complex(*star_chi) - complex(*model_chi)
    return bool(abs(d.real) <= S3_MAX_DIFF and abs(d.imag) <= S3_MAX_DIFF)


def evaluate(rows: list[dict], stars: list[dict]) -> dict:
    """Per variant: each star's area-matched model χ (``wsf.area_matched_chi``), the component
    differences and S3 (all three stars within S3_MAX_DIFF per component)."""
    res = {}
    for name in VARIANTS:
        vr = [r["variants"][name] for r in rows]
        per = []
        for s in stars:
            m = wsf.area_matched_chi(vr, s["isophotal_area_px"])
            per.append(
                {
                    "row": s["row"],
                    "isophotal_area_px": s["isophotal_area_px"],
                    "star_chi_pix": s["star_chi_pix"],
                    "star_abs_chi": s["star_abs_chi"],
                    "model": m,
                    "diff_per_component": [
                        s["star_chi_pix"][0] - m["chi_pix_median"][0],
                        s["star_chi_pix"][1] - m["chi_pix_median"][1],
                    ],
                    "s3_star": s3_star(s["star_chi_pix"], m["chi_pix_median"]),
                }
            )
        res[name] = {
            "sigma_total_px": float(np.sqrt(np.sum(np.square(VARIANTS[name])))),
            "stars": per,
            "s3_pass": all(p["s3_star"] for p in per),
        }
    return res


def main(argv=None) -> int:
    import os

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--stars",
        type=Path,
        default=REPO / "results/w1_stpsf/stars_check.json",
        help="star χ, areas",
    )
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    stars = json.loads(a.stars.read_text())["stars"]
    rows = compute(a.workers)
    out = {
        "provenance": "derived (model PSFs: model_prediction)",
        "stpsf_version": __import__("stpsf").__version__,
        "stpsf_data_version": wsf._data_version(),
        "filter": wsf.FILTER,
        "kernel_sigma_px": KERNEL_SIGMA_PX,
        "box_sigma_px": BOX_SIGMA_PX,
        "s3_max_diff": S3_MAX_DIFF,
        "variants": evaluate(rows, stars),
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    for name, v in out["variants"].items():
        d = [f"{abs(p['model']['abs_chi_median']):.3f}/{p['star_abs_chi']:.3f}" for p in v["stars"]]
        print(name, "S3", v["s3_pass"], "model/star |χ|:", " ".join(d))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
