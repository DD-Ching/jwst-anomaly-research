"""Second MUSE redshift attempt for the `arcs` anti row MACS0308 2157.

Plan: results/arcs/muse-2157-v2/PLAN.md (recentred apertures, STAT noise, peak-contrast gate M2,
[O II] injection limit). Reuses the stream reader and line functions of `arcs_muse_2157.py`; the
cube is streamed once (DATA then STAT), only small boxes are kept in memory, nothing is stored.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path

import arcs_muse_2157 as v1
import numpy as np
from astropy.wcs import WCS

CONTRAST_MIN = 1.5  # ASSUMPTION (PLAN.md item 3)
CONTRAST_DZ = 0.01
RECENTRE_ARCSEC = 1.0  # ASSUMPTION (PLAN.md item 1)
OII = (3726.0, 3728.8)
OII_SIGMA_A = 1.1  # ASSUMPTION: observed Gaussian sigma per doublet component
RECOVER_DZ = 0.002
Z_INJ = (0.40, 1.45)
N_INJ = 40


def stream_data_stat(chunks: Iterable[bytes], boxes: dict[str, tuple[int, int, int, int]]) -> dict:
    """Stream empty primary + DATA + STAT (float32 cubes of equal shape); keep boxes of both."""
    reader = v1._Reader(chunks)
    v1._read_header(reader)
    out: dict = {}
    for ext in ("DATA", "STAT"):
        hdr = v1._read_header(reader)
        if hdr.get("NAXIS") != 3 or hdr.get("BITPIX") != -32:
            raise ValueError(f"expected a float32 3-D {ext} extension")
        if hdr.get("EXTNAME", ext) != ext:
            raise ValueError(f"expected EXTNAME {ext}, got {hdr.get('EXTNAME')}")
        nx, ny, nl = hdr["NAXIS1"], hdr["NAXIS2"], hdr["NAXIS3"]
        cubes = {
            k: np.empty((nl, y1 - y0, x1 - x0), np.float32) for k, (x0, x1, y0, y1) in boxes.items()
        }
        plane_bytes = nx * ny * 4
        for k in range(nl):
            plane = np.frombuffer(reader.read(plane_bytes), ">f4").reshape(ny, nx)
            for name, (x0, x1, y0, y1) in boxes.items():
                cubes[name][k] = plane[y0:y1, x0:x1]
        pad = (-nl * plane_bytes) % v1.BLOCK
        if pad:
            reader.read(pad)
        out[ext] = {"boxes": cubes, "header": hdr}
    return out


def recentre(
    img: np.ndarray, x: float, y: float, r_px: float, n_iter: int = 3
) -> tuple[float, float]:
    """Flux-weighted centroid of positive pixels within ``r_px`` of the current centre."""
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    x0, y0 = x, y
    for _ in range(n_iter):
        m = ((xx - x0) ** 2 + (yy - y0) ** 2 <= r_px**2) & np.isfinite(img)
        w = np.where(m, img - np.nanmedian(img), 0.0)
        w = np.clip(w, 0.0, None)
        if w.sum() <= 0:
            break
        x0, y0 = float((w * xx).sum() / w.sum()), float((w * yy).sum() / w.sum())
        if (x0 - x) ** 2 + (y0 - y) ** 2 > r_px**2:  # keep within the search radius of the start
            s = r_px / np.hypot(x0 - x, y0 - y)
            x0, y0 = x + (x0 - x) * s, y + (y0 - y) * s
    return x0, y0


def stat_noise(stat: np.ndarray, xc: float, yc: float, r_px: float) -> np.ndarray:
    ny, nx = stat.shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    m = (xx - xc) ** 2 + (yy - yc) ** 2 <= r_px**2
    return np.sqrt(np.nansum(stat[:, m], axis=1))


def peak_contrast(z_grid: np.ndarray, score: np.ndarray, dz: float = CONTRAST_DZ) -> float:
    k = int(np.argmax(score))
    far = np.abs(z_grid - z_grid[k]) > dz
    other = float(np.max(score[far])) if far.any() else 0.0
    return float(score[k] / other) if other > 0 else float("inf")


def classify_v2(z_best: float, contrast: float) -> str:
    """Frozen reading of PLAN.md item 4 (applied only when gate M2 passes)."""
    if not contrast >= CONTRAST_MIN:
        return "undetermined"
    if abs(z_best - v1.Z_CLUSTER) < 0.02:
        return "member"
    if z_best > 0.40:
        return "background"
    return "foreground"


def gate_m2(z_best: float, contrast: float, z_ref: float) -> bool:
    return bool(abs(z_best - z_ref) < 0.01 and contrast >= CONTRAST_MIN)


def oii_model(lam: np.ndarray, z: float, flux: float, sigma_a: float = OII_SIGMA_A) -> np.ndarray:
    """[O II] doublet (1:1) flux density on ``lam`` for total line ``flux`` (density x A)."""
    out = np.zeros_like(lam, dtype=float)
    for rest in OII:
        c = rest * (1 + z)
        out += (
            0.5 * flux * np.exp(-0.5 * ((lam - c) / sigma_a) ** 2) / (sigma_a * np.sqrt(2 * np.pi))
        )
    return out


def recovered(
    lam: np.ndarray,
    spec: np.ndarray,
    noise: np.ndarray,
    mask: np.ndarray,
    z_grid: np.ndarray,
    z_inj: float,
    flux: float,
) -> tuple[bool, bool]:
    """(detected, z_identified) for one injection (PLAN.md amendment 1).

    detected: [O II] window S/N >= 5 at z_inj (the frozen limit criterion);
    z_identified: the scan's z_best within 0.002 of z_inj (diagnostic only).
    """
    resid = v1.continuum_subtract(lam, spec + oii_model(lam, z_inj, flux))
    snr = v1.features_at(lam, resid, noise, z_inj, mask).get("[OII]3727", float("nan"))
    scan = v1.scan_redshift(lam, resid, noise, z_grid, mask)
    return bool(np.isfinite(snr) and snr >= v1.SNR_FEATURE), bool(
        abs(scan["z_best"] - z_inj) < RECOVER_DZ
    )


def f90(levels: list[float], fractions: list[float], need: float = 0.9) -> float | None:
    for lev, fr in sorted(zip(levels, fractions, strict=True)):
        if fr >= need:
            return lev
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("results/arcs/muse-2157-v2"))
    ap.add_argument("--r-arcsec", type=float, default=0.6)
    ap.add_argument("--half-box", type=int, default=40)
    args = ap.parse_args()
    reader = v1._Reader(v1._http_chunks(v1.CUBE_URL, 1 << 16))
    v1._read_header(reader)
    hdr = v1._read_header(reader)
    w = WCS(hdr)
    boxes, starts = {}, {}
    for key, src in (("target", v1.TARGET), ("control", v1.CONTROL)):
        x, y = v1._pix(w, src["ra"], src["dec"])
        h = args.half_box
        x0, y0 = max(0, int(round(x)) - h), max(0, int(round(y)) - h)
        x1, y1 = min(hdr["NAXIS1"], int(round(x)) + h), min(hdr["NAXIS2"], int(round(y)) + h)
        boxes[key] = (x0, x1, y0, y1)
        starts[key] = (x - x0, y - y0)
    started = datetime.now(UTC).isoformat(timespec="seconds")
    res = stream_data_stat(v1._http_chunks(v1.CUBE_URL), boxes)
    lam = v1.wavelengths(res["DATA"]["header"])
    pix_arcsec = abs(hdr["CD2_2"]) * 3600
    r_px = args.r_arcsec / pix_arcsec
    z_grid = np.arange(0.0, 1.6, 0.0002)
    report: dict = {
        "plan": "results/arcs/muse-2157-v2/PLAN.md",
        "cube": v1.CUBE_URL,
        "streamed_utc": started,
        "aperture_r_arcsec": args.r_arcsec,
        "contrast_min": CONTRAST_MIN,
        "provenance": {
            "observed": "MUSE cube DATA and STAT values (ESO phase-3 product)",
            "derived": "recentred apertures, spectra, noise, z scan, contrast, class, injections",
            "simulated": "injected [O II] doublets",
            "assumption": "aperture 0.6 arcsec, recentre radius 1.0 arcsec, contrast 1.5, "
            "S/N 5, line windows, [O II] sigma 1.1 A, sky mask factor 3",
        },
        "rows": {},
    }
    spectra: dict = {"lam": lam.tolist()}
    work = {}
    for key, src in (("control", v1.CONTROL), ("target", v1.TARGET)):
        cube = res["DATA"]["boxes"][key]
        stat = res["STAT"]["boxes"][key]
        img = np.nanmean(cube, axis=0)
        sx, sy = starts[key]
        xc, yc = recentre(img, sx, sy, RECENTRE_ARCSEC / pix_arcsec)
        spec = v1.aperture_spectrum(cube, xc, yc, r_px)
        n_emp = v1.noise_spectrum(cube, r_px, exclude=[(xc, yc, 3 * r_px)])
        n_stat = stat_noise(stat, xc, yc, r_px)
        noise = np.fmax(n_emp, n_stat)
        mask = v1.sky_line_mask(lam, noise)
        resid = v1.continuum_subtract(lam, spec)
        scan = v1.scan_redshift(lam, resid, noise, z_grid, mask)
        con = peak_contrast(z_grid, scan["score"])
        row = {
            "name": src["name"],
            "start_xy_box": [sx, sy],
            "recentred_xy_box": [xc, yc],
            "recentre_shift_arcsec": float(np.hypot(xc - sx, yc - sy) * pix_arcsec),
            "z_best": scan["z_best"],
            "peak_contrast": con,
            "features_at_best": v1.features_at(lam, resid, noise, scan["z_best"], mask),
            "features_at_z_cluster": v1.features_at(lam, resid, noise, v1.Z_CLUSTER, mask),
            "features_at_z_phot": v1.features_at(lam, resid, noise, v1.Z_PHOT, mask),
            "median_noise_ratio_emp_over_stat": float(np.nanmedian(n_emp / n_stat)),
            "continuum_median_snr_per_pix": float(np.nanmedian(spec / noise)),
            "n_sky_masked": int(mask.sum()),
        }
        if key == "control":
            row["z_ref"] = v1.CONTROL["z_ref"]
            row["gate_m2_passed"] = gate_m2(scan["z_best"], con, v1.CONTROL["z_ref"])
        report["rows"][key] = row
        spectra[key] = {
            "flux": spec.tolist(),
            "noise_emp": n_emp.tolist(),
            "noise_stat": n_stat.tolist(),
            "score": scan["score"].tolist(),
        }
        work[key] = (spec, noise, mask, (xc, yc))
    gate = report["rows"]["control"]["gate_m2_passed"]
    t = report["rows"]["target"]
    t["class"] = classify_v2(t["z_best"], t["peak_contrast"]) if gate else "not_read_gate_failed"
    # [O II] injection limit at the target (PLAN.md item 5; run whatever the gate outcome).
    spec, noise, mask, _ = work["target"]
    rng = np.random.default_rng(2157)
    levels = [50.0, 100.0, 200.0, 400.0, 800.0, 1600.0]  # 1e-20 erg/s/cm2 (cube flux unit)
    z_inj = rng.uniform(*Z_INJ, N_INJ)
    fracs, fz = [], []
    for lev in levels:
        ok = [recovered(lam, spec, noise, mask, z_grid, float(z), lev) for z in z_inj]
        fracs.append(float(np.mean([a for a, _ in ok])))
        fz.append(float(np.mean([b for _, b in ok])))
    report["injection"] = {
        "flux_unit": res["DATA"]["header"].get("BUNIT", "") + " x Angstrom",
        "levels": levels,
        "fraction_recovered": fracs,
        "fraction_z_identified_diagnostic": fz,
        "z_inj": z_inj.tolist(),
        "F90": f90(levels, fracs),
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(report, indent=1, default=float) + "\n")
    with gzip.open(args.out / "spectra.json.gz", "wt") as f:
        json.dump(spectra, f)
    _sheet(args.out, res, work, lam, spectra, r_px, z_grid)
    print(
        json.dumps(
            {
                k: {kk: v.get(kk) for kk in ("z_best", "peak_contrast", "class")}
                for k, v in report["rows"].items()
            }
        )
    )
    print(json.dumps(report["injection"]["fraction_recovered"]), report["injection"]["F90"])


def _sheet(out, res, work, lam, spectra, r_px, z_grid) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(12, 10), gridspec_kw={"width_ratios": [1, 3]})
    for row, key in enumerate(("control", "target")):
        img = np.nanmean(res["DATA"]["boxes"][key], axis=0)
        ax = axes[row, 0]
        lo, hi = np.nanpercentile(img, [5, 99.5])
        ax.imshow(img, origin="lower", vmin=lo, vmax=hi, cmap="gray")
        xc, yc = work[key][3]
        ax.add_patch(plt.Circle((xc, yc), r_px, fill=False, color="r"))
        ax.set_title(f"{key}: white light, recentred aperture")
        ax = axes[row, 1]
        f = np.asarray(spectra[key]["flux"])
        ax.plot(lam, f, lw=0.5, color="k")
        ax.plot(lam, work[key][1], lw=0.5, color="0.6")
        for z, col in ((v1.Z_CLUSTER, "b"), (v1.Z_PHOT, "r")):
            for rest in list(v1.EMISSION.values()) + list(v1.ABSORPTION.values()):
                c = rest * (1 + z)
                if lam[0] < c < lam[-1]:
                    ax.axvline(c, color=col, lw=0.4, alpha=0.5)
        ax.set_ylim(np.nanpercentile(f, 1) - 50, np.nanpercentile(f, 99.5) + 50)
        ax.set_title(f"{key}: spectrum (blue z {v1.Z_CLUSTER}, red z {v1.Z_PHOT}); grey noise")
    ax = axes[2, 1]
    for key, col in (("control", "b"), ("target", "r")):
        ax.plot(z_grid, spectra[key]["score"], color=col, lw=0.6, label=key)
    ax.legend()
    ax.set_xlabel("z")
    ax.set_title("feature score")
    axes[2, 0].axis("off")
    fig.tight_layout()
    fig.savefig(out / "contact_sheet.png", dpi=80)


if __name__ == "__main__":
    main()
