"""Does the model (STPSF) PSF anisotropy explain the star χ_PSF and set Abell 2744's B-mode floor?

Plan and rules: ``results/w1_stpsf/PLAN.md``. ``psf`` computes STPSF F200W PSFs (``DET_DIST``,
detector sampling, 6″) on a 3-position grid of each NIRCam SW detector and their isophotal (segment)
second-moment distortion χ at several isophote levels, rotated to the sky with the mosaic's PC
matrix. ``stars`` compares the w1_psf stars with the model at their isophotal area. ``screen`` runs
the W1 shear screen with each level's median model χ as a constant PSF distortion (the
``w1_psf_floor`` machinery). All outputs are ``derived`` (model PSFs are ``model_prediction``);
thresholds are ASSUMPTIONs.

    STPSF_PATH=<stpsf-data> python scripts/w1_stpsf_floor.py psf \
        --out results/w1_stpsf/model_chi.json
    python scripts/w1_stpsf_floor.py screen --chi results/w1_stpsf/model_chi.json --work <dir> \
        --out results/w1_stpsf/screens_abell2744.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
REPO = Path(__file__).resolve().parents[1]

FILTER = "F200W"
DETECTORS = tuple(f"NRC{m}{i}" for m in "AB" for i in range(1, 5))
POSITIONS = ((200, 200), (1024, 1024), (1848, 1848))  # detector pixels, ASSUMPTION (3 per detector)
LEVELS = (
    1e-1,
    3e-2,
    1e-2,
    3e-3,
    1e-3,
    3e-4,
)  # isophote / peak, ASSUMPTION (catalogue level unknown)
LEVELS_FINE = tuple(float(v) for v in np.geomspace(0.3, 1e-4, 36))  # area-matching grid (post hoc)
FOV_ARCSEC = 10.0  # 6″ cut the spike arms at the faint levels (code review); edges are flagged
# Abell 2744 F200W mosaic jw02561-o001_t003 SCI header (S3 byte-range read 2026-10-17): PC matrix
MOSAIC_PC = ((-0.75691046852068, 0.6535185863032548), (0.6535185863032548, 0.75691046852068))
S1_MAX_DIFF = 0.05  # per component, ASSUMPTION (PLAN.md)
STAR_CHI = (-0.015, -0.156)  # results/w1_psf/psf_abell2744.json chi_psf_median


def isophotal_chi_area(im: np.ndarray, level: float) -> tuple[complex, int, bool]:
    """Distortion χ = (Qxx - Qyy + 2i Qxy)/(Qxx + Qyy) of the flux-weighted moments of the
    connected pixels above ``level`` × peak that contain the peak (a segment, as the catalogue),
    the segment area (px) and whether it touches the stamp edge (then moments depend on the FOV)."""
    from scipy import ndimage

    lab, _ = ndimage.label(im > level * im.max())
    seg = lab == lab[np.unravel_index(np.argmax(im), im.shape)]
    edge = bool(seg[0].any() or seg[-1].any() or seg[:, 0].any() or seg[:, -1].any())
    y, x = np.nonzero(seg)
    w = im[seg]
    xm, ym = np.average(x, weights=w), np.average(y, weights=w)
    qxx = np.average((x - xm) ** 2, weights=w)
    qyy = np.average((y - ym) ** 2, weights=w)
    qxy = np.average((x - xm) * (y - ym), weights=w)
    return complex(qxx - qyy, 2 * qxy) / (qxx + qyy), int(seg.sum()), edge


def isophotal_chi(im: np.ndarray, level: float) -> complex:
    """χ of :func:`isophotal_chi_area`."""
    return isophotal_chi_area(im, level)[0]


def pixel_chi_to_sky(chi: complex, pc=MOSAIC_PC) -> complex:
    """χ in pixel axes (angle from +x toward +y) → χ with PA east of north, through the mosaic's PC
    matrix (world axis 1 = RA, increasing east). ASSUMPTION: detector sci-frame axes are the
    mosaic's pixel axes (both follow the aperture PA; offsets < 1°)."""
    phi = 0.5 * np.angle(chi)
    de = pc[0][0] * np.cos(phi) + pc[0][1] * np.sin(phi)
    dn = pc[1][0] * np.cos(phi) + pc[1][1] * np.sin(phi)
    return abs(chi) * np.exp(2j * np.arctan2(de, dn))


def s1_consistent(chi_model: complex, chi_star=STAR_CHI) -> bool:
    """S1: model and star χ agree within S1_MAX_DIFF in each component."""
    d = chi_model - complex(*chi_star)
    return bool(abs(d.real) <= S1_MAX_DIFF and abs(d.imag) <= S1_MAX_DIFF)


def _median(c) -> complex:
    c = np.asarray(c, complex)
    return complex(np.median(c.real), np.median(c.imag))


def _angle_deg(c: complex) -> float:
    return float(np.mod(0.5 * np.degrees(np.angle(c)), 180.0))


def _init_worker() -> None:
    import os

    os.environ["OMP_NUM_THREADS"] = "1"


def compute_psfs(workers: int) -> list[dict]:
    from concurrent.futures import ProcessPoolExecutor

    jobs = [(d, x, y) for d in DETECTORS for x, y in POSITIONS]
    with ProcessPoolExecutor(workers, initializer=_init_worker) as ex:
        return list(ex.map(_one_psf, jobs))


def _one_psf(job) -> dict:
    import stpsf

    det, x, y = job
    nc = stpsf.NIRCam()
    nc.filter, nc.detector, nc.detector_position = FILTER, det, (x, y)
    im = nc.calc_psf(fov_arcsec=FOV_ARCSEC, oversample=2)["DET_DIST"].data.astype(float)
    coarse = [isophotal_chi_area(im, t) for t in LEVELS]
    fine = [isophotal_chi_area(im, t) for t in LEVELS_FINE]
    return {
        "detector": det,
        "x": x,
        "y": y,
        "chi_pix": [[c.real, c.imag] for c, _, _ in coarse],
        "edge": [e for _, _, e in coarse],
        "fine_area_px": [a for _, a, _ in fine],
        "fine_chi_pix": [[c.real, c.imag] for c, _, _ in fine],
        "fine_edge": [e for _, _, e in fine],
    }


def area_matched_chi(rows: list[dict], area_px: float) -> dict:
    """Post hoc (PLAN.md addendum): for each grid PSF, the pixel χ at the fine level whose segment
    area is closest to ``area_px`` (the catalogue ``isophotal_area``; model and mosaic pixels are
    both ≈ 0.031″); the median over grid PSFs, the per-PSF |χ| range and how many touch the edge."""
    c, lv, edge = [], [], 0
    for r in rows:
        k = int(np.argmin(np.abs(np.log(np.asarray(r["fine_area_px"], float) / area_px))))
        c.append(complex(*r["fine_chi_pix"][k]))
        lv.append(LEVELS_FINE[k])
        edge += bool(r["fine_edge"][k])
    med = _median(c)
    return {
        "level_median": float(np.median(lv)),
        "chi_pix_median": [med.real, med.imag],
        "abs_chi_median": abs(med),
        "angle_pix_deg": _angle_deg(med),
        "abs_chi_min_max": [float(np.abs(c).min()), float(np.abs(c).max())],
        "n_edge": edge,
    }


def summarize(rows: list[dict]) -> dict:
    """Per isophote level: each grid PSF's χ rotated to the sky, their per-component median, |χ|
    range, edge count and S1 (median over grid PSFs, not over screened sources: deviation)."""
    out = []
    for k, t in enumerate(LEVELS):
        c = np.array([pixel_chi_to_sky(complex(*r["chi_pix"][k])) for r in rows])
        sky = _median(c)
        out.append(
            {
                "level": t,
                "chi_sky_median": [sky.real, sky.imag],
                "abs_chi_median": abs(sky),
                "pa_sky_deg": _angle_deg(sky),
                "abs_chi_min_max": [float(np.abs(c).min()), float(np.abs(c).max())],
                "n_edge": int(sum(r["edge"][k] for r in rows)),
                "s1_consistent": s1_consistent(sky),
            }
        )
    return {"levels": out}


def star_check(chi_json: Path, stars_json: Path, cat_path: Path) -> dict:
    """Each w1_psf star's catalogue χ (pixel frame, from ``orientation``) vs the model χ at the
    level whose segment area matches the star's ``isophotal_area`` (post hoc; PLAN.md addendum)."""
    from astropy.table import Table
    from w1_psf_floor import moments_chi

    rows = json.loads(chi_json.read_text())["psfs"]
    cat = Table.read(cat_path)
    out = []
    for s in json.loads(stars_json.read_text())["stars"]:
        r = cat[s["row"]]
        chi = complex(moments_chi(r["semimajor_sigma"], r["semiminor_sigma"], r["orientation"]))
        m = area_matched_chi(rows, float(r["isophotal_area"]))
        out.append(
            {
                "row": s["row"],
                "isophotal_area_px": float(r["isophotal_area"]),
                "isophotal_flux_jy": float(r["isophotal_flux"]),
                "star_chi_pix": [chi.real, chi.imag],
                "star_abs_chi": abs(chi),
                "star_angle_pix_deg": _angle_deg(chi),
                "model": m,
                "diff_per_component": [
                    chi.real - m["chi_pix_median"][0],
                    chi.imag - m["chi_pix_median"][1],
                ],
                "level_times_flux": m["level_median"] * float(r["isophotal_flux"]),
            }
        )
    return {"provenance": "derived", "label": "post hoc (deviation)", "stars": out}


def screens(chi_json: Path, work: Path, name: str = "abell2744") -> dict:
    from types import SimpleNamespace

    import w1_psf_floor as wpf

    ish, es = wpf.ish, wpf.es
    spec = ish.field_spec(name)
    (work / name).mkdir(parents=True, exist_ok=True)
    model, shapes, extra, _, _ = ish.ir.load_field(name, work / name, spec)
    sargs = SimpleNamespace(**ish.SHEAR_DEFAULTS)
    sargs.max_radius = spec["max_radius"]
    psf = es.psf_sigma_px(shapes)
    res = [wpf.screen_summary(es.shear_screen(model, shapes, sargs, psf, extra), "isotropic")]
    for lv in json.loads(chi_json.read_text())["levels"]:
        sargs.psf_chi = complex(*lv["chi_sky_median"])
        r = es.shear_screen(model, shapes, sargs, psf, extra)
        res.append(wpf.screen_summary(r, f"constant model chi, isophote {lv['level']:g} x peak"))
    return {"field": name, "provenance": "derived", "psf_sigma_px": psf, "screens": res}


def _data_version() -> str:
    import os

    v = Path(os.environ.get("STPSF_PATH", "")) / "version.txt"
    return v.read_text().strip() if v.is_file() else "unknown"


def main(argv=None) -> int:
    import os

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("psf")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    st = sub.add_parser("stars", help="area-matched model χ for the w1_psf stars (post hoc)")
    st.add_argument("--chi", type=Path, required=True)
    st.add_argument("--stars", type=Path, default=REPO / "results/w1_psf/psf_abell2744.json")
    st.add_argument("--catalog", type=Path, required=True, help="the F200W _cat.ecsv")
    st.add_argument("--out", type=Path, required=True)
    s = sub.add_parser("screen")
    s.add_argument("--chi", type=Path, required=True)
    s.add_argument("--work", type=Path, required=True)
    s.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.cmd == "psf":
        rows = compute_psfs(a.workers)
        import stpsf

        out = {
            "provenance": "model_prediction",
            "stpsf_version": stpsf.__version__,
            "stpsf_data_version": _data_version(),
            "filter": FILTER,
            "fov_arcsec": FOV_ARCSEC,
            "positions": POSITIONS,
            "levels_fine": LEVELS_FINE,
            "mosaic_pc": MOSAIC_PC,
            "star_chi_median": STAR_CHI,
            "s1_max_diff": S1_MAX_DIFF,
            **summarize(rows),
            "psfs": rows,
        }
        summary = out["levels"]
    elif a.cmd == "stars":
        out = star_check(a.chi, a.stars, a.catalog)
        summary = [{k: v for k, v in r.items() if k != "model"} | r["model"] for r in out["stars"]]
    else:
        out = screens(a.chi, a.work)
        summary = out["screens"]
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
