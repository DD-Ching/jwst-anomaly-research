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

**Fixes after the final review (PLAN_fixes.md):** a channel is valid only when every aperture
spaxel has finite DATA and STAT > 0; boxcars touching invalid or sky-line channels are not usable;
the MAD rescale and the running medians use usable/valid channels only. ``peak_check`` runs the
pre-registered 6000 Å checks for 1316 (``peak_6000.png``).

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
]


def boxcar(x: np.ndarray, n: int) -> np.ndarray:
    c = np.convolve(np.nan_to_num(x), np.ones(n), mode="same")
    return c


def running_median(x: np.ndarray, n: int = 101) -> np.ndarray:
    """NaN-aware running median (NaN channels are ignored; PLAN_fixes.md fix 2)."""
    x = np.asarray(x, float)
    half = n // 2
    out = np.full(x.size, np.nan)
    for k in range(x.size):
        seg = x[max(0, k - half) : k + half + 1]
        if np.isfinite(seg).any():
            out[k] = np.nanmedian(seg)
    return out


def mad_sigma(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if x.size else float("nan")


def aperture_spectrum(data: np.ndarray, stat: np.ndarray, ap: np.ndarray):
    """Aperture flux and variance per channel; a channel is valid only when every aperture spaxel
    has finite DATA and STAT > 0 (PLAN_fixes.md fix 2). Invalid channels get NaN."""
    d, v = data[:, ap], stat[:, ap]
    valid = np.all(np.isfinite(d) & np.isfinite(v) & (v > 0), axis=1)
    f = np.where(valid, np.sum(np.where(np.isfinite(d), d, 0.0), axis=1), np.nan)
    var = np.where(valid, np.sum(np.where(np.isfinite(v), v, 0.0), axis=1), np.nan)
    return f, var, valid


def spectrum_snr(f: np.ndarray, var: np.ndarray, valid: np.ndarray) -> dict:
    """Continuum-subtracted 6-channel boxcar S/N with the rescale computed on usable channels only.

    ``usable`` = boxcar positions whose channels are all valid and not sky-line channels (variance
    > 1.5 × its running median over valid channels). The MAD rescale uses usable positions only
    (PLAN_fixes.md fix 1); other positions get NaN S/N.
    """
    f = np.where(valid, f - running_median(np.where(valid, f, np.nan)), 0.0)
    v = np.where(valid, var, 0.0)
    sky = valid & (var > 1.5 * running_median(np.where(valid, var, np.nan)))
    bad = (~valid) | sky
    usable = boxcar(bad.astype(float), BOX) == 0
    fb, vb = boxcar(f, BOX), boxcar(v, BOX)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = np.where(usable & (vb > 0), fb / np.sqrt(vb), np.nan)
    resc = max(mad_sigma(snr), 1.0)
    return {
        "snr": snr / resc,
        "usable": usable & np.isfinite(snr),
        "sigma_box": np.where(usable, np.sqrt(np.where(vb > 0, vb, np.nan)) * resc, np.nan),
        "resc": resc,
        "frac_invalid_channels": float(np.mean(~valid)),
        "frac_sky_channels": float(np.mean(sky)),
    }


def load_cube(path: Path):
    with fits.open(path) as h:
        data = np.asarray(h["DATA"].data, float)
        stat = np.asarray(h["STAT"].data, float)
        hdr = h["DATA"].header
    w = WCS(hdr)
    lam = w.spectral.pixel_to_world(np.arange(data.shape[0])).to_value("Angstrom")
    cel = w.celestial
    pix = float(np.sqrt(abs(np.linalg.det(cel.pixel_scale_matrix)))) * 3600.0
    return data, stat, lam, cel, pix


def measure(path: Path, ra: float, dec: float, fwhm: float) -> dict:
    data, stat, lam, cel, pix = load_cube(path)
    x0, y0 = cel.world_to_pixel_values(ra, dec)
    yy, xx = np.indices(data.shape[1:])
    r = np.hypot(xx - x0, yy - y0) * pix
    r_ap = max(0.5 * fwhm, 0.4)
    ap = r <= r_ap
    n_ap = int(ap.sum())
    finite = float(np.isfinite(data[:, ap]).mean()) if n_ap else 0.0
    if finite < 0.5:
        return {"no_data": True, "finite_fraction": round(finite, 3), "n_spaxels": n_ap}
    f, var, valid = aperture_spectrum(data, stat, ap)
    sp = spectrum_snr(f, var, valid)
    unit = 1e-20  # MUSE BUNIT 10**(-20)*erg/s/cm**2/Angstrom; boxcar sum × 1.25 Å per channel
    dlam = float(np.median(np.diff(lam)))
    return {
        "no_data": False,
        "finite_fraction": round(finite, 3),
        "lam": lam,
        "snr": sp["snr"],
        "usable": sp["usable"],
        "sigma_flux": sp["sigma_box"] * unit * dlam,
        "r_ap_arcsec": round(r_ap, 2),
        "n_spaxels": n_ap,
        "noise_rescale": round(sp["resc"], 2),
        "frac_invalid_channels": round(sp["frac_invalid_channels"], 4),
        "frac_sky_channels": round(sp["frac_sky_channels"], 4),
        "dlam": dlam,
    }


def window_result(m: dict, lo: float, hi: float) -> dict:
    sel = (m["lam"] >= lo) & (m["lam"] <= hi)
    if not sel.any():
        return {"window_A": [round(lo), round(hi)], "in_range": False}
    clean = sel & m["usable"]
    k = np.flatnonzero(clean)
    best = int(k[np.argmax(m["snr"][k])]) if len(k) else None
    det = best is not None and m["snr"][best] >= SNR_DET
    return {
        "window_A": [round(lo), round(hi)],
        "in_range": True,
        "frac_not_usable": round(1 - clean.sum() / sel.sum(), 3),
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


def exceedance_runs(snr: np.ndarray, usable: np.ndarray, thr: float) -> int:
    """Number of separate runs of usable boxcar positions with S/N ≥ ``thr`` (chance peaks)."""
    hit = usable & (np.nan_to_num(snr, nan=-np.inf) >= thr)
    return int(np.sum(hit[1:] & ~hit[:-1]) + (1 if hit.size and hit[0] else 0))


def peak_check(path: Path, ra: float, dec: float, fwhm: float, lam0: float, window, out_png):
    """The 6000 Å check (PLAN_fixes.md): validity, narrow band, offset apertures, chance rate."""
    data, stat, lam, cel, pix = load_cube(path)
    x0, y0 = cel.world_to_pixel_values(ra, dec)
    yy, xx = np.indices(data.shape[1:])
    r_ap = max(0.5 * fwhm, 0.4)

    def spec_at(dx_arcsec, dy_arcsec):
        ap = np.hypot(xx - (x0 + dx_arcsec / pix), yy - (y0 + dy_arcsec / pix)) * pix <= r_ap
        f, var, valid = aperture_spectrum(data, stat, ap)
        return spectrum_snr(f, var, valid), int(ap.sum())

    tgt, _ = spec_at(0.0, 0.0)
    k = int(np.argmin(np.abs(lam - lam0)))
    win = (lam >= window[0]) & (lam <= window[1])
    snr_k = float(tgt["snr"][k]) if np.isfinite(tgt["snr"][k]) else None
    # narrow band: the boxcar channels at k minus the side-band mean
    lo, hi = k - BOX // 2, k + BOX - BOX // 2
    side = np.r_[k - 50 : k - 20, k + 21 : k + 51]
    nb = np.nansum(data[lo:hi], axis=0) - BOX * np.nanmean(data[side], axis=0)
    from scipy.ndimage import gaussian_filter

    sm = gaussian_filter(np.nan_to_num(nb), fwhm / 2.355 / pix)
    jy, ix = np.unravel_index(np.argmax(sm), sm.shape)
    d_peak = float(np.hypot(ix - x0, jy - y0) * pix)
    offs = []
    for ang in range(0, 360, 60):
        dx, dy = 1.5 * np.cos(np.radians(ang)), 1.5 * np.sin(np.radians(ang))
        sp, n = spec_at(dx, dy)
        offs.append(
            {
                "pa_deg": ang,
                "n_spaxels": n,
                "snr_at_peak": float(sp["snr"][k]) if np.isfinite(sp["snr"][k]) else None,
                "usable_at_peak": bool(sp["usable"][k]),
                "chance_runs_in_window": (
                    exceedance_runs(sp["snr"][win], sp["usable"][win], snr_k)
                    if snr_k is not None
                    else None
                ),
                "usable_positions_in_window": int(np.sum(sp["usable"][win])),
            }
        )
    runs = [o["chance_runs_in_window"] for o in offs if o["chance_runs_in_window"] is not None]
    expected = float(np.mean(runs)) if runs else None
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
    ax[0].imshow(sm, origin="lower", cmap="viridis")
    ax[0].add_patch(plt.Circle((x0, y0), r_ap / pix, fill=False, color="w"))
    ax[0].plot(ix, jy, "r+")
    ax[0].set_title(f"narrow band {lam[lo]:.0f}-{lam[hi - 1]:.0f} Å (smoothed)", fontsize=8)
    sel = (lam > lam0 - 150) & (lam < lam0 + 150)
    ax[1].plot(lam[sel], tgt["snr"][sel], "k-", lw=0.8, label="target")
    for ang in range(0, 360, 60):
        sp, _ = spec_at(1.5 * np.cos(np.radians(ang)), 1.5 * np.sin(np.radians(ang)))
        ax[1].plot(lam[sel], sp["snr"][sel], lw=0.4, alpha=0.6)
    ax[1].axvline(lam0, color="r", lw=0.5)
    ax[1].axhline(SNR_DET, color="grey", lw=0.5, ls="--")
    ax[1].set_xlabel("Å")
    ax[1].set_ylabel("boxcar S/N")
    ax[1].set_title("target (black) and six 1.5″ offset apertures", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=90)
    plt.close(fig)
    checks = {
        "valid_and_unmasked": bool(tgt["usable"][k]),
        "snr_ge_threshold": bool(snr_k is not None and snr_k >= SNR_DET),
        "narrowband_peak_within_r_ap": bool(d_peak <= r_ap),
        "no_offset_aperture_ge_3": all(
            (o["snr_at_peak"] is None) or (o["snr_at_peak"] < 3) for o in offs
        ),
        "expected_chance_peaks_lt_0p05": bool(expected is not None and expected < 0.05),
    }
    return {
        "lambda_A": round(float(lam[k]), 2),
        "snr_after_fixes": None if snr_k is None else round(snr_k, 2),
        "noise_rescale_after_fixes": round(tgt["resc"], 2),
        "narrowband_channels_A": [round(float(lam[lo]), 1), round(float(lam[hi - 1]), 1)],
        "narrowband_peak_offset_arcsec": round(d_peak, 2),
        "r_ap_arcsec": round(r_ap, 2),
        "offset_apertures": offs,
        "expected_chance_peaks_per_window_search": expected,
        "checks": checks,
        "class": "line to be checked"
        if all(checks.values())
        else "not significant under the calibration",
    }


def run(cut_arcsec: float, keep_dir: Path | None = None) -> dict:
    out = {
        "label": "feasibility check (derived); expected fluxes model_prediction",
        "plan": "PLAN_spectra.md; fixes PLAN_fixes.md",
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
                p = (keep_dir or Path(tmp)) / f"{dp}_{t}.fits"
                with urllib.request.urlopen(url, timeout=300) as r, open(p, "wb") as fh:
                    fh.write(r.read())
                size = p.stat().st_size
                m = measure(p, ra, dec, fwhm)
                if t == 1316:
                    peak = peak_check(
                        p,
                        ra,
                        dec,
                        fwhm,
                        6000.3,
                        (
                            REST["CIII]"] * (1 + TARGETS[t]["z"][0]),
                            REST["CIII]"] * (1 + TARGETS[t]["z"][1]),
                        ),
                        paths.repo_root() / "results/arcs/abells1063-cats/vetting/peak_6000.png",
                    )
                    out["peak_6000_1316"] = peak
                if keep_dir is None:
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
                    "frac_invalid_channels": m["frac_invalid_channels"],
                    "frac_sky_channels": m["frac_sky_channels"],
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
    p.add_argument("--keep-dir", type=Path, default=None, help="keep the cutouts here")
    args = p.parse_args(argv)
    print(json.dumps(run(args.cut_arcsec, args.keep_dir), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
