"""Frozen line test for one source in the 4043 F444W grism c1d spectra (see results/grism_4043)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter

LINES_UM = {"HeI_10833": 1.0833, "Pagamma": 1.0941, "Pabeta": 1.2822}  # rest, vacuum
Z_RANGE = (2.68, 3.10)  # DJA eazy 95 % interval for f200w_2915
CONT_PIX = 51  # ASSUMPTION, running-median continuum width (pixels)
BOX_PIX = 5  # ASSUMPTION, boxcar width (pixels)
SNR_MIN = 5.0  # ASSUMPTION
DLAM_MAX_UM = 0.003  # ASSUMPTION, 3 nm
GATE_RECOVERY = 0.90  # ASSUMPTION
INJ_SNR = 7.0  # ASSUMPTION
INJ_SIGMA_UM = 0.0015  # ASSUMPTION, injected Gaussian sigma
SEED = 1


def clean(wave, flux, dq, n_input):
    wave, flux = np.asarray(wave, dtype=float), np.asarray(flux, dtype=float)
    ok = np.isfinite(wave) & np.isfinite(flux) & (np.asarray(dq) == 0) & (np.asarray(n_input) > 0)
    order = np.argsort(wave[ok])
    return wave[ok][order], flux[ok][order]


def boxcar_snr(flux: np.ndarray) -> tuple[np.ndarray, float]:
    """Boxcar sums of the continuum-subtracted flux, divided by their robust sigma."""
    resid = flux - median_filter(flux, size=CONT_PIX, mode="nearest")
    box = np.convolve(resid, np.ones(BOX_PIX), mode="same")
    sigma = 1.4826 * np.median(np.abs(box - np.median(box)))
    return box / sigma, float(sigma)


def line_snr(
    wave: np.ndarray, flux: np.ndarray, rest_um: float, z_range=Z_RANGE
) -> tuple[float, float]:
    """Peak S/N and its wavelength inside the line's window; (nan, nan) if not covered."""
    snr, _ = boxcar_snr(flux)
    lo, hi = rest_um * (1 + z_range[0]), rest_um * (1 + z_range[1])
    w = (wave >= lo) & (wave <= hi)
    if w.sum() < BOX_PIX:
        return float("nan"), float("nan")
    i = np.flatnonzero(w)[np.argmax(snr[w])]
    return float(snr[i]), float(wave[i])


def decide_one(per_grism: dict[str, dict[str, tuple[float, float]]], line: str) -> bool:
    """S/N >= SNR_MIN in every grism and peak wavelengths within DLAM_MAX_UM."""
    vals = [per_grism[g][line] for g in per_grism]
    ok = all(np.isfinite(s) and s >= SNR_MIN for s, _ in vals)
    return bool(ok and (max(w for _, w in vals) - min(w for _, w in vals) <= DLAM_MAX_UM))


def decide(per_grism: dict[str, dict[str, tuple[float, float]]]) -> dict:
    """Detection per line (``decide_one``); any detected line gives ``emission_detected``."""
    detected = {line: decide_one(per_grism, line) for line in LINES_UM}
    return {
        "detected": detected,
        "verdict": "emission_detected" if any(detected.values()) else "no_line_detected",
    }


def covers(wave: np.ndarray, center_um: float) -> bool:
    return bool(wave.min() <= center_um <= wave.max())


def inject(wave, flux, center_um, snr, sigma=None):
    """Add a Gaussian line whose 5-pixel boxcar S/N is ``snr`` (spectrum's robust sigma).

    A centre off the spectrum leaves the flux unchanged (the line is not observed there).
    """
    if not covers(wave, center_um):
        return flux
    if sigma is None:
        _, sigma = boxcar_snr(flux)
    prof = np.exp(-0.5 * ((wave - center_um) / INJ_SIGMA_UM) ** 2)
    k = np.argmin(np.abs(wave - center_um))
    boxed = prof[max(k - BOX_PIX // 2, 0) : k + BOX_PIX // 2 + 1].sum()
    return flux + prof * snr * sigma / boxed


def run_injections(spectra: dict, n: int, snr: float, seed: int) -> dict:
    """Recovery over all n draws (the frozen gate) and over draws on both spectra."""
    rng = np.random.default_rng(seed)
    sig = {g: boxcar_snr(fl)[1] for g, (_, fl) in spectra.items()}
    out = {}
    for line, rest in LINES_UM.items():
        hits = hits_cov = n_cov = 0
        for z in rng.uniform(*Z_RANGE, size=n):
            c = rest * (1 + z)
            both = all(covers(wv, c) for wv, _ in spectra.values())
            per = {}
            for g, (wv, fl) in spectra.items():
                f2 = inject(wv, fl, c, snr, sig[g])
                per[g] = {line: line_snr(wv, f2, rest)}
            hit = decide_one(per, line)
            hits += hit
            n_cov += both
            hits_cov += hit and both
        out[line] = {
            "recovery": hits / n,
            "n_both_covered": n_cov,
            "recovery_both_covered": round(hits_cov / n_cov, 3) if n_cov else None,
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--c1d", nargs=2, type=Path, required=True, help="grismr and grismc c1d files")
    ap.add_argument("--source-id", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-inj", type=int, default=200)
    args = ap.parse_args()
    from astropy.io import fits

    spectra, meta = {}, {}
    for path in args.c1d:
        with fits.open(path) as h:
            d = h[1].data
            row = d[d["SOURCE_ID"] == args.source_id][0]
            g = "grismr" if "grismr" in path.name else "grismc"
            spectra[g] = clean(row["WAVELENGTH"], row["FLUX"], row["DQ"], row["N_INPUT"])
            meta[g] = {
                "file": path.name,
                "cal_ver": h[0].header.get("CAL_VER"),
                "n_pix": len(spectra[g][0]),
                "ra": float(row["SOURCE_RA"]),
                "dec": float(row["SOURCE_DEC"]),
                "wave_um": [
                    round(float(spectra[g][0].min()), 4),
                    round(float(spectra[g][0].max()), 4),
                ],
            }
    gate = run_injections(spectra, args.n_inj, INJ_SNR, seed=SEED)
    gate_pass = all(v["recovery"] >= GATE_RECOVERY for v in gate.values())
    per = {
        g: {ln: line_snr(wv, fl, r) for ln, r in LINES_UM.items()}
        for g, (wv, fl) in spectra.items()
    }
    res = decide(per) if gate_pass else {"verdict": "no_verdict_gate_failed"}
    summary = {
        "source_id": args.source_id,
        "spectra": meta,
        "params": {
            "snr_min": SNR_MIN,
            "dlam_max_um": DLAM_MAX_UM,
            "gate_recovery": GATE_RECOVERY,
            "inj_snr": INJ_SNR,
            "z_range": Z_RANGE,
            "cont_pix": CONT_PIX,
            "box_pix": BOX_PIX,
            "inj_sigma_um": INJ_SIGMA_UM,
            "n_inj": args.n_inj,
            "seed": SEED,
        },
        "injection_recovery": gate,
        "gate_pass": gate_pass,
        "observed_peak": {
            g: {ln: [round(s, 2), round(w, 4)] for ln, (s, w) in v.items()} for g, v in per.items()
        },
        **res,
    }
    args.out.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
