"""Emission-line check of ``arcs`` targets in public MUSE cubes via ESO SODA cutouts (feasibility).

Frozen rules: ``results/arcs/abells1063-cats/vetting/PLAN_spectra.md``. For each (target, cube)
pair a small spatial cutout (``--cut-arcsec`` radius, all wavelengths) is fetched from ESO's SODA
sync service into a temporary directory and deleted. An aperture of radius max(0.5 × FWHM, 0.4″) is
summed (DATA; variance from STAT). A 6-channel running boxcar gives S/N per position; the noise is
rescaled by the MAD spread of the whole spectrum when that exceeds 1. A detection is boxcar S/N ≥ 5
inside a target window on a channel not flagged as a sky line; otherwise the 5σ boxcar limit is
reported.

**Deviation from the plan (after the first run, before any result was used):** the planned sky-
annulus flag marked 99.9 % of channels around 1316 because the annulus sits on cluster light
(continuum, not sky lines), and it missed the 5577 Å sky line next to 1557. Instead: a 101-channel
running median is subtracted from the aperture spectrum (continuum and ICL), and a channel is a
sky-line channel when the aperture variance exceeds 1.5 × its own 101-channel running median. Cubes
with no finite pixels in the aperture (ObsCore footprints include masked edges) are reported as
``no_data``. Expected line fluxes come from the F090W magnitude and assumed rest EWs
(model_prediction).

Output: ``results/arcs/abells1063-cats/vetting/spectra.json`` (derived; cube pixels observed).

    python scripts/arcs_vet_spectra.py
"""

from __future__ import annotations

import json
import sys
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from jwst_anomaly import paths

SODA = "https://dataportal.eso.org/dataPortal/soda/sync"
REST = {"OII": 3728.0, "CIII]": 1908.7, "CIV": 1549.0}  # Å, vacuum-ish centroids
BOX = 6  # channels (1.25 Å each): ASSUMPTION, ~ [O II] doublet width at z 1.35
SNR_DET = 5.0  # ASSUMPTION
F090W_AB = {
    1316: 28.86,
    1557: 28.69,
}  # pipeline isophotal_abmag (observed; catalogue sha256 3eecb26c…)
TARGETS = {
    1316: {
        "radec": (342.1744892170039, -44.53093324709693),
        "z": (1.87, 2.36),
        "lines": ["CIII]", "CIV"],
    },
    1557: {"radec": (342.1902009068879, -44.53319879686645), "z": (1.27, 1.43), "lines": ["OII"]},
}
# (dp_id, FWHM arcsec from ObsCore s_resolution, exposure s, targets)
CUBES = [
    ("ADP.2017-03-28T12:46:01.331", 1.081, 15885.3, (1557,)),  # 095.A-0653 combined
    ("ADP.2017-03-23T15:58:03.937", 1.447, 7877.1, (1316,)),  # 60.A-9345 combined
    ("ADP.2025-11-20T15:47:55.061", 0.681, 2459.5, (1316, 1557)),  # 116.28VK single OB
]


def boxcar(x: np.ndarray, n: int) -> np.ndarray:
    c = np.convolve(np.nan_to_num(x), np.ones(n), mode="same")
    return c


def running_median(x: np.ndarray, n: int = 101) -> np.ndarray:
    from scipy.ndimage import median_filter

    return median_filter(np.nan_to_num(x), size=n, mode="nearest")


def mad_sigma(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def measure(path: Path, ra: float, dec: float, fwhm: float) -> dict:
    with fits.open(path) as h:
        data = np.asarray(h["DATA"].data, float)
        stat = np.asarray(h["STAT"].data, float)
        hdr = h["DATA"].header
    w = WCS(hdr)
    nz = data.shape[0]
    lam = w.spectral.pixel_to_world(np.arange(nz)).to_value("Angstrom")
    cel = w.celestial
    x0, y0 = cel.world_to_pixel_values(ra, dec)
    pix = float(np.sqrt(abs(np.linalg.det(cel.pixel_scale_matrix)))) * 3600.0
    yy, xx = np.indices(data.shape[1:])
    r = np.hypot(xx - x0, yy - y0) * pix
    r_ap = max(0.5 * fwhm, 0.4)
    ap = r <= r_ap
    n_ap = int(ap.sum())
    finite = float(np.isfinite(data[:, ap]).mean()) if n_ap else 0.0
    if finite < 0.5:
        return {"no_data": True, "finite_fraction": round(finite, 3), "n_spaxels": n_ap}
    f = np.nansum(data[:, ap], axis=1)
    v = np.nansum(stat[:, ap], axis=1)
    f = f - running_median(f)  # continuum and cluster light
    fb, vb = boxcar(f, BOX), boxcar(v, BOX)
    snr = fb / np.sqrt(vb)
    resc = max(mad_sigma(snr), 1.0)
    snr /= resc
    sky_line = v > 1.5 * running_median(v)
    sky_snr = np.where(boxcar(sky_line.astype(float), BOX) > 0, np.inf, 0.0)
    unit = 1e-20  # MUSE BUNIT 10**(-20)*erg/s/cm**2/Angstrom; boxcar sum × 1.25 Å per channel
    dlam = float(np.median(np.diff(lam)))
    sigma_flux = np.sqrt(vb) * resc * unit * dlam
    return {
        "no_data": False,
        "finite_fraction": round(finite, 3),
        "lam": lam,
        "snr": snr,
        "sky_snr": sky_snr,
        "sigma_flux": sigma_flux,
        "r_ap_arcsec": round(r_ap, 2),
        "n_spaxels": n_ap,
        "noise_rescale": round(resc, 2),
        "dlam": dlam,
    }


def window_result(m: dict, lo: float, hi: float) -> dict:
    sel = (m["lam"] >= lo) & (m["lam"] <= hi) & np.isfinite(m["snr"])
    if not sel.any():
        return {"window_A": [round(lo), round(hi)], "in_range": False}
    clean = sel & (np.abs(m["sky_snr"]) < SNR_DET)
    k = np.flatnonzero(clean)
    best = int(k[np.argmax(m["snr"][k])]) if len(k) else None
    det = best is not None and m["snr"][best] >= SNR_DET
    return {
        "window_A": [round(lo), round(hi)],
        "in_range": True,
        "frac_sky_flagged": round(1 - clean.sum() / sel.sum(), 3),
        "max_snr_clean": round(float(m["snr"][best]), 2) if best is not None else None,
        "max_snr_lambda_A": round(float(m["lam"][best]), 1) if best is not None else None,
        "detected": bool(det),
        "limit_5sigma_cgs_median": float(f"{SNR_DET * np.nanmedian(m['sigma_flux'][clean]):.3g}")
        if clean.any()
        else None,
    }


def expected_fluxes(ab: float, lam_obs: float, z: float) -> dict:
    fnu = 10 ** (-0.4 * (ab + 48.6))
    flam = fnu * 2.998e18 / lam_obs**2
    return {f"EW_rest_{ew}": float(f"{flam * ew * (1 + z):.3g}") for ew in (30, 100, 300)}


def run(cut_arcsec: float) -> dict:
    out = {
        "label": "feasibility check (derived); expected fluxes model_prediction",
        "plan": "PLAN_spectra.md",
        "assumptions": {"boxcar_channels": BOX, "snr_det": SNR_DET, "cut_arcsec": cut_arcsec},
        "f090w_isophotal_ab": F090W_AB,
        "results": [],
    }
    with tempfile.TemporaryDirectory() as tmp:
        for dp, fwhm, texp, tg in CUBES:
            for t in tg:
                ra, dec = TARGETS[t]["radec"]
                url = (
                    f"{SODA}?ID=ivo://eso.org/ID?{dp}&CIRCLE={ra}%20{dec}%20{cut_arcsec / 3600:.5f}"
                )
                p = Path(tmp) / f"{dp}_{t}.fits"
                with urllib.request.urlopen(url, timeout=300) as r, open(p, "wb") as fh:
                    fh.write(r.read())
                size = p.stat().st_size
                m = measure(p, ra, dec, fwhm)
                p.unlink()
                if m["no_data"]:
                    out["results"].append(
                        {"target": t, "cube": dp, "cutout_bytes": size, "no_data": True}
                        | {k: m[k] for k in ("finite_fraction", "n_spaxels")}
                    )
                    continue
                z0, z1 = TARGETS[t]["z"]
                res = {
                    "target": t,
                    "cube": dp,
                    "fwhm_arcsec": fwhm,
                    "exptime_s": texp,
                    "cutout_bytes": size,
                    "r_ap_arcsec": m["r_ap_arcsec"],
                    "n_spaxels": m["n_spaxels"],
                    "noise_rescale": m["noise_rescale"],
                    "lines": {},
                    "whole_range": window_result(m, 4750, 9350),
                }
                for ln in TARGETS[t]["lines"]:
                    lo, hi = REST[ln] * (1 + z0), REST[ln] * (1 + z1)
                    w = window_result(m, max(lo, 4750), min(hi, 9350))
                    zc = 0.5 * (z0 + z1)
                    w["expected_flux_cgs"] = expected_fluxes(F090W_AB[t], REST[ln] * (1 + zc), zc)
                    res["lines"][ln] = w
                out["results"].append(res)
    d = paths.repo_root() / "results" / "arcs" / "abells1063-cats" / "vetting"
    (d / "spectra.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def main(argv=None) -> int:
    import argparse

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cut-arcsec", type=float, default=2.4, help="cutout radius")
    args = p.parse_args(argv)
    print(json.dumps(run(args.cut_arcsec), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
