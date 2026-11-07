"""MIRI flux-ratio test of Mahler+2022 ICLv2 multiple images (frozen plan:
``results/miri_confirm/PLAN_glens.md``).

D-032's leave-one-out flux-ratio rule applied to the MIRI merged catalogue of a pipeline run:
nearest MIRI source within 0.5'' of each catalogued image, ``<band>_isophotal_flux`` for F770W and
F1000W, implied source magnitude ``m + 2.5 log10 |mu|`` with mu from ICLv2 ``best.par``.

Usage: python scripts/miri_fluxratio.py RUN_DIR OUT.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

MATCH_ARCSEC = 0.5  # ASSUMPTION (plan)
MIN_SNR = 10.0  # D-032
MAX_ABS_MU = 20.0  # D-032
OUTLIER_MAG = 0.75  # D-032 ASSUMPTION
MAX_AREA_PX = 2000  # ASSUMPTION (plan)
BANDS = ("f770w", "f1000w")


def leave_one_out(src_mag: np.ndarray, system: np.ndarray, usable: np.ndarray) -> np.ndarray:
    """Residual of each usable image against the median of the other usable images of its system."""
    res = np.full(len(src_mag), np.nan)
    for s in np.unique(system):
        idx = np.flatnonzero((system == s) & usable)
        if len(idx) < 2:
            continue
        for i in idx:
            others = idx[idx != i]
            res[i] = src_mag[i] - np.median(src_mag[others])
    return res


def _finite(x):
    """JSON-safe value: non-finite floats become None (strict JSON has no NaN)."""
    if isinstance(x, float) and not np.isfinite(x):
        return None
    if isinstance(x, dict):
        return {k: _finite(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_finite(v) for v in x]
    return x


def main(run_dir: Path, out: Path) -> dict:
    model, files, par = lc.load_model("smacs0723-iclv2")
    images, zml = lc.image_list("smacs0723-iclv2", files, par)
    z = lensmodel.image_redshifts(images, zml)
    ok_z = np.isfinite(z)
    mu = np.full(len(images), np.nan)
    pred = model.evaluate(np.asarray(images["ra"])[ok_z], np.asarray(images["dec"])[ok_z], z[ok_z])
    mu[ok_z] = np.asarray(pred["magnification"], float)

    src = Table.read(run_dir / "smacs0723_miri" / "sources.ecsv")
    ci = SkyCoord(images["ra"], images["dec"], unit="deg")
    cs = SkyCoord(src["ra"], src["dec"], unit="deg")
    j, sep, _ = ci.match_to_catalog_sky(cs)
    sep = sep.arcsec
    matched = sep <= MATCH_ARCSEC
    uids, counts = np.unique(j[matched], return_counts=True)
    shared = np.isin(j, uids[counts > 1]) & matched
    system = np.array([lensmodel.system_key(s) for s in images["system"]])

    rows, summary = [], {"run": run_dir.name, "n_images": len(images), "bands": {}}
    for band in BANDS:
        flux = np.asarray(src[f"{band}_isophotal_flux"], float)[j]
        err = np.asarray(src[f"{band}_isophotal_flux_err"], float)[j]
        area = np.asarray(src[f"{band}_isophotal_area"], float)[j]
        mag = np.asarray(src[f"{band}_isophotal_abmag"], float)[j]
        det = np.asarray(src[f"{band}_detected"], bool)[j]
        snr = np.where(err > 0, flux / err, np.nan)
        usable = (
            matched
            & ~shared
            & det
            & (snr >= MIN_SNR)
            & (np.abs(mu) <= MAX_ABS_MU)
            & (area <= MAX_AREA_PX)
            & np.isfinite(mag)
            & np.isfinite(mu)
        )
        src_mag = mag + 2.5 * np.log10(np.abs(mu))
        res = leave_one_out(src_mag, system, usable)
        tested = np.isfinite(res)
        flagged = tested & (np.abs(res) > OUTLIER_MAG)
        summary["bands"][band] = {
            "matched": int(matched.sum()),
            "usable": int(usable.sum()),
            "tested": int(tested.sum()),
            "systems_tested": sorted(set(system[tested])),
            "rms_mag": float(np.sqrt(np.nanmean(res[tested] ** 2))) if tested.any() else None,
            "max_abs_mag": float(np.nanmax(np.abs(res[tested]))) if tested.any() else None,
            "flagged": [str(images["image_id"][i]) for i in np.flatnonzero(flagged)],
        }
        for i in np.flatnonzero(matched):
            rows.append(
                {
                    "band": band,
                    "image": str(images["image_id"][i]),
                    "system": system[i],
                    "z": float(z[i]),
                    "mu": float(mu[i]),
                    "sep_arcsec": round(float(sep[i]), 3),
                    "source_uid": str(src["source_uid"][j[i]]),
                    "shared": bool(shared[i]),
                    "detected": bool(det[i]),
                    "snr": float(snr[i]),
                    "area_px": float(area[i]),
                    "mag": float(mag[i]),
                    "src_mag": float(src_mag[i]),
                    "usable": bool(usable[i]),
                    "loo_residual_mag": None
                    if not np.isfinite(res[i])
                    else round(float(res[i]), 3),
                }
            )
    summary["rows"] = rows
    out.write_text(json.dumps(_finite(summary), indent=1, default=str, allow_nan=False) + "\n")
    return summary


if __name__ == "__main__":
    s = main(Path(sys.argv[1]), Path(sys.argv[2]))
    print(json.dumps(s["bands"], indent=1))
