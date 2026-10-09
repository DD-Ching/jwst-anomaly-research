"""W5 feasibility in Euclid Q1: can deeper counts open θ_E < 6′, and would shapes do better?

Measures galaxy densities in the three Euclid Q1 deep fields (MER catalogue, IRSA TAP; counts
only, no rows downloaded) and forecasts the gain of the D-063 count screen over Legacy Surveys DR10
with the clustering-inflated scatter measured there (`results/w5_counts/screen_summary.json`), plus
the S/N of a radial-shear test around a point-like negative-mass lens. Writes
`results/w5_counts/euclid_q1_feasibility.json`. Forecasts are `model_prediction`; the inputs marked
ASSUMPTION below are not measured.

    python scripts/w5_euclid_feasibility.py
"""

from __future__ import annotations

import json
import math
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests

TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
TABLE = "euclid_q1_mer_catalogue"
# Points inside each Euclid Deep Field (Q1); radius of the counting disc in degrees.
FIELDS = {"EDF-F": (52.93, -28.09), "EDF-S": (61.24, -48.42), "EDF-N": (269.73, 66.02)}
RADIUS_DEG = 0.1
VIS_MAGS = (23.5, 24.5, 25.0)
# Extended VIS detections, clean flags (ASSUMPTION: point_like_prob < 0.1 removes stars).
SELECTION = "vis_det = 1 AND spurious_flag = 0 AND det_quality_flag = 0 AND point_like_prob < 0.1"
DR10_DENSITY = 37374.6  # deg⁻², r < 23.5 (results/w5_counts/numcounts.ecsv)
CLUSTERING_RATIO = (0.5, 1.0)  # Euclid / DR10 clustering variance at fixed angle (ASSUMPTION range)
SIGMA_GAMMA = 0.3  # shape noise per component, unweighted SExtractor moments (ASSUMPTION)
DETECT_SNR = 6.0  # S/N needed for a shear detection after trials (ASSUMPTION)
# Point mass γ = (θ_E/θ)² averaged over θ_E < θ < 2 θ_E (model_prediction).
MEAN_GAMMA = 2 * math.log(2) / 3


def vis_flux_ujy(mag: float) -> float:
    return 10 ** ((23.9 - mag) / 2.5)


def count(ra: float, dec: float, mag: float) -> int:
    # CONTAINS uses IRSA's spatial index; a plain RA/Dec box ran > 5 min for 0.25 deg² (2026-10-09).
    q = (
        f"SELECT COUNT(*) AS n FROM {TABLE} WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {ra}, {dec}, {RADIUS_DEG})) AND {SELECTION} "
        f"AND flux_detection_total > {vis_flux_ujy(mag):.4f}"
    )
    for attempt in range(4):
        try:
            r = requests.post(TAP, data={"QUERY": q, "FORMAT": "csv"}, timeout=300)
            r.raise_for_status()
            return int(r.text.strip().splitlines()[-1])
        except (requests.RequestException, ValueError):
            if attempt == 3:
                raise
            time.sleep(2 ** (attempt + 1))
    raise AssertionError("unreachable")


def z_std_dr10(summ: dict, scale: float) -> float:
    """Scatter of the DR10 screen's Z at one filter scale, mean of the two regions (derived)."""
    return float(np.mean([summ[f"{r} {scale}"]["z_std"] for r in ("desA", "desB")]))


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    out = repo / "results" / "w5_counts"
    jobs = [(f, m) for f in FIELDS for m in VIS_MAGS]
    with ThreadPoolExecutor(len(jobs)) as ex:
        ns = list(ex.map(lambda j: count(*FIELDS[j[0]], j[1]), jobs))
    area_deg2 = math.pi * RADIUS_DEG**2
    dens = {f"{f} VIS<{m}": n / area_deg2 for (f, m), n in zip(jobs, ns, strict=True)}
    mean = {m: float(np.mean([dens[f"{f} VIS<{m}"] for f in FIELDS])) for m in VIS_MAGS}
    slope = math.log10(mean[24.5] / mean[23.5]) / 1.0  # d log10 N / dm between 23.5 and 24.5

    summ = json.loads((out / "screen_summary.json").read_text())
    counts_gain = {}
    for scale in (2.0, 4.0, 8.0):
        z_std = z_std_dr10(summ, scale)
        var_dr10 = z_std**2  # in units of the DR10 Poisson variance; clustering part = z_std² − 1
        for c in CLUSTERING_RATIO:
            var_eu = DR10_DENSITY / mean[24.5] + c * (z_std**2 - 1)
            counts_gain[f"{scale:g} arcmin, clustering x{c:g}"] = math.sqrt(var_dr10 / var_eu)
    # DR10 S/N at each scale relative to 8′ (ε ≈ 0.5 there): (θ/8) · z_std(8′) / z_std(θ).
    rel_snr_dr10 = {
        f"{s:g}": (s / 8.0) * z_std_dr10(summ, 8.0) / z_std_dr10(summ, s) for s in (2.0, 4.0, 8.0)
    }
    n_arcmin2 = mean[24.5] / 3600.0
    theta_min = DETECT_SNR * SIGMA_GAMMA / (MEAN_GAMMA * math.sqrt(3 * math.pi * n_arcmin2))
    res = {
        "provenance": {
            "densities": "observed (IRSA TAP counts, Euclid Q1 MER)",
            "forecasts": "model_prediction",
            "assumptions": {
                "selection": SELECTION,
                "clustering_ratio": CLUSTERING_RATIO,
                "sigma_gamma": SIGMA_GAMMA,
                "detect_snr": DETECT_SNR,
            },
        },
        "radius_deg": RADIUS_DEG,
        "density_deg2": dens,
        "mean_density_deg2": {f"VIS<{m}": v for m, v in mean.items()},
        "density_ratio_vis24p5_over_dr10": mean[24.5] / DR10_DENSITY,
        "count_slope_dlog10N_dm_23p5_24p5": slope,
        "counts_snr_gain_over_dr10": counts_gain,
        "dr10_snr_relative_to_8arcmin": rel_snr_dr10,
        "shear_theta_e_min_arcmin": theta_min,
    }
    path = out / "euclid_q1_feasibility.json"
    path.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
