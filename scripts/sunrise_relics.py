"""WHL0137 (Sunrise): the JWST image list through the four RELICS HST-era lens models.

The RELICS v1 models (Cerny et al. 2018; HLSP DOI 10.17909/T9SP45): Glafic, Lenstool, WSLAP+ and
Zitrin LTM. Their deflection maps are in arcsec at D_LS/D_S = 1 (D-033 for Lenstool), so
alpha(z) = alpha_map * r(z), r = D_ds/D_s. Only small pixel patches around each catalogued
image (Scofield+2025 Table 3, observed) are read with HTTP byte ranges.

Per model and image (model_prediction from the published maps; positions observed):
* the deflection, and mu at z = 6.2 from finite differences of r(6.2) * alpha, checked against
  the published z = 6.2 magnification map at the same pixel (validates units, axes and scaling);
* per system and trial redshift z, the source-plane rms of the images about their mean.

System 1 (the Sunrise arc, z = 6.2) is the known case. System 2 (redshift unknown) does not fit
the JWST-era Scofield+2025 M5-G1 maps at any z (scripts/sunrise_zfit.py); here the question is
whether the HST-era models reproduce it at some z.

Cosmology: flat LCDM, h = 0.7, Omega_m = 0.3; cluster z = 0.566 (D-033).

Usage: python scripts/sunrise_relics.py [--out results/sunrise/relics_positions.json]
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests
from sunrise_scofield import IMAGES, system_of
from sunrise_zfit import LEVEL_FACTOR, Z_GRID, dratio

BASE = "https://archive.stsci.edu/hlsps/relics/whl0137m08/models/{m}/v1/hlsp_relics_model_model_whl0137m08_{m}_v1_{k}.fits"
MODELS = ("glafic", "lenstool", "wslap", "zitrin-ltm-gauss")
KEYS = {"defx": "x-arcsec-deflect", "defy": "y-arcsec-deflect", "mag": "z06p2-magnif"}
Z_ARC = 6.2
MU_CHECK = 100.0  # ASSUMPTION: |mu| below which finite differences must match the published map
HALF = 3  # patch half-width in pixels (7x7 patch)
DTYPES = {-32: ">f4", -64: ">f8"}


def _get(session: requests.Session, url: str, start: int, stop: int) -> bytes:
    """One byte range, retried with back-off on errors, 429/5xx or short reads."""
    for attempt in range(8):
        try:
            r = session.get(url, headers={"Range": f"bytes={start}-{stop - 1}"}, timeout=60)
            if r.status_code == 206 and len(r.content) == stop - start:
                return r.content
        except requests.RequestException:
            pass
        time.sleep(min(2.0**attempt, 30.0))
    raise RuntimeError(f"range {start}-{stop} of {url} failed")


def read_header(session: requests.Session, url: str):
    """Primary header (any number of 2880-byte blocks) -> (header, data offset in bytes)."""
    from astropy.io import fits

    raw = b""
    while True:
        raw += _get(session, url, len(raw), len(raw) + 2880)
        cards = [raw[k : k + 80] for k in range(len(raw) - 2880, len(raw), 80)]
        if any(c.rstrip() == b"END" for c in cards):
            break
        if len(raw) > 20 * 2880:
            raise RuntimeError(f"{url}: no END card")
    hdr = fits.Header.fromstring(raw.decode("ascii"))
    if hdr["NAXIS"] != 2 or hdr["BITPIX"] not in DTYPES:
        raise RuntimeError(f"{url}: unexpected layout")
    return hdr, len(raw)


def frame_check(wcs) -> float:
    """Pixel scale (arcsec) of a north-up, east-left grid; anything else raises."""
    cd = wcs.pixel_scale_matrix * 3600.0
    if not (
        cd[0, 0] < 0 < cd[1, 1]
        and abs(cd[0, 1]) < 1e-6 * abs(cd[0, 0])
        and abs(cd[1, 0]) < 1e-6 * abs(cd[1, 1])
    ):
        raise RuntimeError(f"grid is not north-up, east-left: {cd}")
    if not np.isclose(-cd[0, 0], cd[1, 1], rtol=1e-6):
        raise RuntimeError(f"non-square pixels: {cd}")
    return float(cd[1, 1])


def read_patch(
    session, url: str, hdr, offset: int, i0: int, j0: int, half: int = HALF
) -> np.ndarray:
    """Pixels [j0-half, j0+half] x [i0-half, i0+half] (0-based) of one map."""
    nx, bpp = hdr["NAXIS1"], abs(hdr["BITPIX"]) // 8
    if not (half <= i0 < nx - half and half <= j0 < hdr["NAXIS2"] - half):
        raise RuntimeError(f"patch at ({i0}, {j0}) leaves the map")
    rows = []
    for j in range(j0 - half, j0 + half + 1):
        start = offset + (j * nx + i0 - half) * bpp
        raw = _get(session, url, start, start + (2 * half + 1) * bpp)
        rows.append(np.frombuffer(raw, DTYPES[hdr["BITPIX"]]))
    return np.array(rows, float)


def bilinear(a: np.ndarray, c: float, r: float) -> float:
    """a at fractional (column c, row r)."""
    c0, r0 = int(np.floor(c)), int(np.floor(r))
    fc, fr = c - c0, r - r0
    return float(
        a[r0, c0] * (1 - fc) * (1 - fr)
        + a[r0, c0 + 1] * fc * (1 - fr)
        + a[r0 + 1, c0] * (1 - fc) * fr
        + a[r0 + 1, c0 + 1] * fc * fr
    )


def local_lens(
    defx: np.ndarray, defy: np.ndarray, step: float, c: float, r: float, scale: float
) -> dict:
    """Convergence, magnification and psi at (c, r) for alpha = scale * map.

    +x runs West along +i and +y North along +j (north-up, east-left grid)."""
    dax_dj, dax_di = np.gradient(scale * defx, step)
    day_dj, day_di = np.gradient(scale * defy, step)
    pxx, pyy = bilinear(dax_di, c, r), bilinear(day_dj, c, r)
    pxy = 0.5 * (bilinear(dax_dj, c, r) + bilinear(day_di, c, r))
    det = (1 - pxx) * (1 - pyy) - pxy**2
    return {"kappa": 0.5 * (pxx + pyy), "mu": 1.0 / det, "psi": [pxx, pxy, pyy]}


def measure_model(model: str, workers: int = 4) -> dict:
    from astropy.wcs import WCS

    session = requests.Session()
    urls = {k: BASE.format(m=model, k=v) for k, v in KEYS.items()}
    heads = {k: read_header(session, u) for k, u in urls.items()}
    wcs = WCS(heads["defx"][0])
    step = frame_check(wcs)
    for k in ("defy", "mag"):  # same grid for all three maps
        h = heads[k][0]
        if any(
            h.get(c) != heads["defx"][0].get(c)
            for c in ("NAXIS1", "NAXIS2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2")
        ):
            raise RuntimeError(f"{model}: {k} grid differs from defx")
    r_arc = dratio(Z_ARC)
    ref = IMAGES["1.1e"]

    def one(item):
        name, (ra, dec) = item
        i, j = (float(v) for v in wcs.world_to_pixel_values(ra, dec))
        i0, j0 = int(round(i)), int(round(j))
        p = {k: read_patch(session, urls[k], *heads[k], i0, j0) for k in KEYS}
        c, r = HALF + i - i0, HALF + j - j0
        lens = local_lens(p["defx"], p["defy"], step, c, r, r_arc)
        unit = local_lens(p["defx"], p["defy"], step, c, r, 1.0)
        return {
            "image": name,
            "x_arcsec": (ra - ref[0]) * np.cos(np.deg2rad(dec)) * -3600.0,  # West
            "y_arcsec": (dec - ref[1]) * 3600.0,
            "alpha_x": bilinear(p["defx"], c, r),  # arcsec at D_LS/D_S = 1
            "alpha_y": bilinear(p["defy"], c, r),
            "kappa_fd_z6p2": lens["kappa"],
            "mu_fd_z6p2": lens["mu"],  # signed (parity)
            "psi_unit": unit["psi"],  # [xx, xy, yy] second derivatives at D_LS/D_S = 1
            "mu_map_z6p2": bilinear(p["mag"], c, r),  # the published maps hold |mu|
        }

    with ThreadPoolExecutor(workers) as ex:
        rows = list(ex.map(one, IMAGES.items()))
    return {"model": model, "pixel_arcsec": step, "url_defx": urls["defx"], "images": rows}


def source_rms(rows: list[dict], z: float) -> float:
    """rms distance (arcsec) of the images' source positions from the mean of their own knot
    (1.1, 1.2, 1.7 or 2; `system_of`), pooled over the knots, for a source at z."""
    r = dratio(z)
    d2 = []
    for knot in sorted({system_of(q["image"]) for q in rows}):
        sel = [q for q in rows if system_of(q["image"]) == knot]
        b = np.array(
            [[q["x_arcsec"] - r * q["alpha_x"], q["y_arcsec"] - r * q["alpha_y"]] for q in sel]
        )
        d2.extend(np.sum((b - b.mean(axis=0)) ** 2, axis=1))
    return float(np.sqrt(np.mean(d2)))


def jacobian(q: dict, z: float) -> np.ndarray:
    """Lens Jacobian A = I - r(z) psi at one image."""
    xx, xy, yy = (dratio(z) * v for v in q["psi_unit"])
    return np.array([[1 - xx, -xy], [-xy, 1 - yy]])


def image_plane_offsets(rows: list[dict], z: float) -> list[float]:
    """Linear image-plane offset |A(z)^-1 (beta_i - <beta>)| of each image of one knot, arcsec.

    Approximate: it diverges near critical curves, so it is used only for system 2's
    moderately magnified images."""
    r = dratio(z)
    b = np.array(
        [[q["x_arcsec"] - r * q["alpha_x"], q["y_arcsec"] - r * q["alpha_y"]] for q in rows]
    )
    d = b - b.mean(axis=0)
    out = []
    for q, dq in zip(rows, d, strict=True):
        out.append(float(np.hypot(*np.linalg.solve(jacobian(q, z), dq))))
    return out


def zscan(rows: list[dict], z_grid=Z_GRID) -> np.ndarray:
    return np.array([source_rms(rows, z) for z in z_grid])


def summarise(meas: dict, z_grid=Z_GRID) -> dict:
    rows = meas["images"]
    sys1 = [q for q in rows if system_of(q["image"]).startswith("1.")]
    sys2 = [q for q in rows if system_of(q["image"]) == "2"]
    # validation: away from critical curves (|mu| < MU_CHECK, ASSUMPTION), |finite-difference mu|
    # matches the published map (which holds |mu|)
    ok = [q for q in rows if abs(q["mu_map_z6p2"]) < MU_CHECK]
    ratios = [abs(q["mu_fd_z6p2"]) / q["mu_map_z6p2"] for q in ok]
    s1, s2 = zscan(sys1, z_grid), zscan(sys2, z_grid)
    k1, k2 = int(np.argmin(s1)), int(np.argmin(s2))
    level = LEVEL_FACTOR * float(s1[k1])
    acc = z_grid[s2 <= level]
    acc1 = z_grid[s1 <= level]
    return {
        "model": meas["model"],
        "mu_check_n": len(ok),
        "mu_ratio_median": float(np.median(ratios)) if ratios else float("nan"),
        "mu_ratio_range": [float(min(ratios)), float(max(ratios))] if ratios else None,
        "sys1_rms_at_6p2": source_rms(sys1, Z_ARC),
        "sys1_best_z": float(z_grid[k1]),
        "sys1_best_rms": float(s1[k1]),
        "sys1_accepted_z": [float(acc1.min()), float(acc1.max())],
        "level": level,
        "sys2_best_z": float(z_grid[k2]),
        "sys2_best_rms": float(s2[k2]),
        "sys2_accepted_z": [float(acc.min()), float(acc.max())] if acc.size else None,
        "sys2_image_offsets_at_best_z": dict(
            zip(
                [q["image"] for q in sys2],
                image_plane_offsets(sys2, float(z_grid[k2])),
                strict=True,
            )
        ),
        "sys2_abs_mu_at_best_z": [
            abs(1.0 / np.linalg.det(jacobian(q, float(z_grid[k2])))) for q in sys2
        ],
        "sys2_rms_at": {str(z): source_rms(sys2, z) for z in (1.0, 2.0, 3.0, 6.2)},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=Path("results/sunrise/relics_positions.json"))
    ap.add_argument("--models", nargs="+", default=list(MODELS))
    args = ap.parse_args(argv)
    with ThreadPoolExecutor(len(args.models)) as ex:
        meas = list(ex.map(measure_model, args.models))
    out = {
        "source": "RELICS v1 lens models, HLSP DOI 10.17909/T9SP45, accessed 2026-10-13",
        "access": "HTTP byte ranges, 7x7-pixel patches; checksums not verifiable (no download)",
        "positions": "Scofield+2025 Table 3 (observed, JWST frame)",
        "provenance": (
            "model_prediction (map values); summaries derived; "
            "LEVEL_FACTOR and MU_CHECK are ASSUMPTIONs"
        ),
        "measurements": meas,
        "summary": [summarise(m) for m in meas],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    for s in out["summary"]:
        print(json.dumps(s))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
