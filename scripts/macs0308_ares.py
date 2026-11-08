"""MACS J0308.9+2645: image order and time delays of SN Ares in the pre-JWST RELICS maps.

SN Ares (Larison et al. 2026, arXiv:2610.33161v1) is a lensed SN IIb at z_spec 1.276 in image 1.3
of a triply imaged host; three JWST-era models (not public) say it is the first image to arrive,
the others following after 72-79 yr. This script asks whether the public RELICS GLAFIC v2 and
Zitrin LTM-Gauss v1 maps give the same order. Plan:
``results/lens_consistency/macs0308_ares/PLAN.md`` (frozen before the run).

Reuses the Fermat-potential machinery of ``smacs_at2024jfz.py`` (line integrals of the deflection,
straight vs L-shaped path check). Outputs (``--out``): ``summary.json`` (model_prediction per
model, decided by :func:`decide`) and ``images.png`` (predicted images over each model's κ map).

Usage::

    python scripts/macs0308_ares.py --out results/lens_consistency/macs0308_ares
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402
import smacs_at2024jfz as at  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

MODELS = ("macs0308-relics-glafic", "macs0308-relics-zitrin")
# SN Ares, arXiv:2610.33161v1 Sect. II: 03:08:53.418 +26:45:58.63 (J2000)
ARES_RADEC = (15.0 * (3 + 8 / 60 + 53.418 / 3600), 26 + 45 / 60 + 58.63 / 3600)
Z_S = 1.276  # z_spec of SN Ares (same paper)
DAYS_PER_YR = 365.25
#: Published JWST-era predictions (paper text; not public models): later arrivals in years.
PUBLISHED_LATER_YR = {"glafic_jwst": (72.0, 77.0), "allingham": (76.4, 79.1)}
GATE_S_ARCSEC = 0.3  # ASSUMPTION: re-imaged SN within this of the input position
EARLIER_YR = 1.0  # ASSUMPTION: an image this much earlier than Ares is a conflict
LATER_YR = 10.0  # ASSUMPTION: every other image this much later is agreement
GRID_HALF_WIDTH = 90.0  # arcsec; maps are cropped by NaN outside their footprint
GRID_STEP = 0.1


def decide(
    path_ok: bool,
    sn_sep_arcsec: float,
    other_delays_yr: list[float],
    gate_s: float = GATE_S_ARCSEC,
    earlier: float = EARLIER_YR,
    later: float = LATER_YR,
) -> str:
    """Frozen per-model outcome (PLAN.md).

    ``other_delays_yr``: arrival time of every image except the one at the SN position, minus
    the SN image's (positive = later). Returns ``stopped: path``, ``stopped: no SN image``,
    ``conflicts``, ``agrees`` or ``undetermined``.
    """
    if not path_ok:
        return "stopped: path"
    if not (np.isfinite(sn_sep_arcsec) and sn_sep_arcsec <= gate_s):
        return "stopped: no SN image"
    if any(d <= -earlier for d in other_delays_yr):
        return "conflicts"
    if len(other_delays_yr) >= 2 and all(d >= later for d in other_delays_yr):
        return "agrees"
    return "undetermined"


def run_model(name: str) -> dict:
    model, files, _ = lc.load_model(name)
    grid = lensmodel.DeflectionGrid.compute(model, GRID_HALF_WIDTH, GRID_STEP)
    at.TNS_RADEC = ARES_RADEC  # solve() reads the reference position from this module constant
    res = at.solve(model, grid, Z_S)
    imgs = res.get("images", [])
    ref = [r for r in imgs if r["reference_image"]]
    others = [r for r in imgs if not r["reference_image"]]
    delays_yr = [r["delay_days"] / DAYS_PER_YR for r in others]
    sn_sep = res.get("tns_to_reference_arcsec", float("nan"))
    outcome = decide(bool(res.get("path_check_pass", False)), sn_sep, delays_yr)
    x0, y0 = (float(v[0]) for v in model.to_frame([ARES_RADEC[0]], [ARES_RADEC[1]]))
    for r in imgs:  # drop the AT 2024jfz-specific epoch fields
        for k in ("testable_epochs_mjd", "next_window_mjd"):
            r.pop(k, None)
        r["delay_yr"] = round(r["delay_days"] / DAYS_PER_YR, 2)
        r["sep_from_ares_arcsec"] = r.pop("sep_from_tns_arcsec")
    return {
        "model": name,
        "provenance": "model_prediction",
        "files_sha256": {k: v[1] for k, v in lc.MODELS[name]["files"].items()},
        "ares_frame_xy": [round(x0, 3), round(y0, 3)],
        "z_s": Z_S,
        "dls_ds": res["dls_ds"],
        "beta_xy": res["beta_xy"],
        "n_images": res["n_images"],
        "ares_to_reference_arcsec": sn_sep,
        "reference_magnification": ref[0]["magnification"] if ref else None,
        "days_per_arcsec2": res.get("days_per_arcsec2"),
        "path_check_max_rel_diff": res.get("path_check_max_rel_diff"),
        "path_check_pass": res.get("path_check_pass"),
        "other_delays_yr": [round(d, 2) for d in delays_yr],
        "outcome": outcome,
        "images": imgs,
        "_model": model,
    }


def plot(results: list[dict], out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.io import fits

    fig, axes = plt.subplots(1, len(results), figsize=(6 * len(results), 6))
    for ax, res in zip(np.atleast_1d(axes), results, strict=True):
        files = lc.model_files(res["model"], ("kappa_map",))
        kap = fits.getdata(files["kappa_map"]).astype(float)
        ax.imshow(np.log10(np.clip(kap, 1e-2, None)), origin="lower", cmap="Greys")
        model = res["_model"]
        for r in res["images"]:
            i, j = model.wcs.world_to_pixel_values(r["ra"], r["dec"])
            mark = "*" if r["reference_image"] else "o"
            ax.plot(i, j, mark, ms=12, mfc="none", mec="red")
            ax.annotate(
                f"{r['delay_yr']:+.1f} yr, μ {r['magnification']:+.1f}",
                (i, j),
                color="red",
                fontsize=8,
                xytext=(5, 5),
                textcoords="offset points",
            )
        i0, j0 = model.wcs.world_to_pixel_values(*ARES_RADEC)
        ax.plot(i0, j0, "x", ms=10, color="blue")
        ax.set_title(f"{res['model']}: {res['outcome']} (z_s {Z_S})", fontsize=9)
    fig.suptitle(
        "SN Ares (blue x) and predicted images (red; * = image at the SN), κ map (log)", fontsize=9
    )
    fig.tight_layout()
    fig.savefig(out / "images.png", dpi=80)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    results = [run_model(m) for m in MODELS]
    plot(results, args.out)
    summary = {
        "provenance": "model_prediction",
        "source": "RELICS MACS0308 maps (Acebron et al. 2018) at the SN Ares position "
        "(arXiv:2610.33161v1); plan PLAN.md",
        "ares_radec": [round(v, 6) for v in ARES_RADEC],
        "assumptions": {
            "gate_s_arcsec": GATE_S_ARCSEC,
            "earlier_yr": EARLIER_YR,
            "later_yr": LATER_YR,
            "path_tol": at.PATH_TOL,
            "cosmology": "FlatLambdaCDM H0 70, Om0 0.3 (map cosmology)",
        },
        "published_later_arrivals_yr": PUBLISHED_LATER_YR,
        "models": [{k: v for k, v in r.items() if k != "_model"} for r in results],
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for r in summary["models"]:
        print(
            r["model"],
            r["outcome"],
            r["n_images"],
            r["other_delays_yr"],
            r["ares_to_reference_arcsec"],
            r["path_check_max_rel_diff"],
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
