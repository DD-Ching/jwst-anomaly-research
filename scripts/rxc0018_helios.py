"""RXC J0018.5+1626: image order and next-image prediction of SN Helios in the pre-JWST RELICS map.

SN Helios (Fujimoto et al. 2026, arXiv:2609.30440v2) is a lensed Type II SN at z 3.34 with one
reported image; the JWST-era models in its Table 5 (not public) put the next image 2.32-6.0 yr
after discovery, a few arcsec apart. This script asks whether the public RELICS Lenstool v1 map
predicts the same order. Plan: ``results/lens_consistency/rxc0018_helios/PLAN.md`` (frozen).

Reuses ``macs0308_ares.py`` / ``smacs_at2024jfz.py`` (Fermat potential, path check). Outputs
(``--out``): ``summary.json`` (model_prediction, decided by :func:`decide`), ``posthoc.json`` and
``images.png``.

Usage::

    python scripts/rxc0018_helios.py --out results/lens_consistency/rxc0018_helios
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

MODEL = "rxc0018-relics-lenstool"
# SN Helios, arXiv:2609.30440v2 Sect. 2: 00:18:31.50 +16:26:18.0 (J2000)
HELIOS_RADEC = (15.0 * (18 / 60 + 31.50 / 3600), 16 + 26 / 60 + 18.0 / 3600)
Z_S = 3.34  # z_spec of SN Helios (same paper)
DAYS_PER_YR = 365.25


def _hms(ra_s: float, dec_m: float, dec_s: float, ra_m: int = 18, dec_d: int = 16) -> tuple:
    return (15.0 * (ra_m / 60 + ra_s / 3600), dec_d + dec_m / 60 + dec_s / 3600)


#: Published next image per JWST-era model (Table 5): (RA, Dec) deg, delay yr, mu.
PUBLISHED_NEXT = {
    "ALP": (_hms(31.54, 26, 0.3), 2.32, 7.2),
    "zitrin-analytic": (_hms(31.84, 26, 4.3), 2.70, 4.9),
    "glafic": (_hms(31.78, 25, 59.3), 3.09, 10.9),
    "WSLAP+": (_hms(31.96, 25, 57.0), 4.18, 25.1),
    "LENSTOOL": (_hms(31.73, 25, 58.2), 6.0, 3.1),
}
#: Published Helios magnification per model and past-image delays (yr ago; paper text).
PUBLISHED_HELIOS_MU = {
    "ALP": 11.2,
    "zitrin-analytic": 6.1,
    "glafic": 11.5,
    "WSLAP+": 30.6,
    "LENSTOOL": 5.3,
}
PUBLISHED_PAST_YR = {"ALP": 46.0, "LENSTOOL": 42.7}
GATE_S_ARCSEC = 0.3  # ASSUMPTION: re-imaged SN within this of the input position
EARLY_YR = 1.0  # ASSUMPTION: an earlier image within this is still on the plateau at discovery
BRIGHT_RATIO = 0.5  # ASSUMPTION: |mu| >= this x |mu_Helios| would be seen in the discovery mosaic
NEXT_WINDOW_YR = (1.0, 10.0)  # ASSUMPTION (published 2.32-6.0 yr)
NEXT_RADIUS_ARCSEC = 5.0  # ASSUMPTION ("agree to within a few arcseconds")
GRID_HALF_WIDTH = 90.0  # arcsec
GRID_STEP = 0.1


def published_centre() -> tuple[float, float]:
    """Mean (RA, Dec) in deg of the five published next-image positions."""
    pos = np.array([v[0] for v in PUBLISHED_NEXT.values()])
    return float(pos[:, 0].mean()), float(pos[:, 1].mean())


def decide(
    path_ok: bool,
    sn_sep_arcsec: float,
    others: list[dict],
    early: float = EARLY_YR,
    bright: float = BRIGHT_RATIO,
    window: tuple[float, float] = NEXT_WINDOW_YR,
    radius: float = NEXT_RADIUS_ARCSEC,
    gate_s: float = GATE_S_ARCSEC,
) -> str:
    """Frozen outcome (PLAN.md).

    ``others``: every image except the SN's, as dicts with ``delay_yr`` (minus the SN's; positive =
    later), ``mu_ratio`` (mu / mu_Helios) and ``sep_centre_arcsec`` (to :func:`published_centre`).
    Returns ``stopped: path``, ``stopped: no SN image``, ``conflicts``, ``agrees``, ``differs`` or
    ``undetermined``.
    """
    if not path_ok:
        return "stopped: path"
    if not (np.isfinite(sn_sep_arcsec) and sn_sep_arcsec <= gate_s):
        return "stopped: no SN image"
    if any(-early <= r["delay_yr"] < 0 and abs(r["mu_ratio"]) >= bright for r in others):
        return "conflicts"
    later = [r for r in others if r["delay_yr"] > 0]
    if not later:
        return "undetermined" if not others else "differs"
    first = min(later, key=lambda r: r["delay_yr"])
    in_window = window[0] <= first["delay_yr"] <= window[1]
    if in_window and first["sep_centre_arcsec"] <= radius:
        return "agrees"
    return "differs"


def run_model(name: str = MODEL) -> dict:
    from astropy.coordinates import SkyCoord

    model, _, _ = lc.load_model(name)
    grid = lensmodel.DeflectionGrid.compute(model, GRID_HALF_WIDTH, GRID_STEP)
    at.TNS_RADEC = HELIOS_RADEC  # solve() reads the reference position from this module constant
    res = at.solve(model, grid, Z_S)
    imgs = res.get("images", [])
    centre = SkyCoord(*published_centre(), unit="deg")
    for r in imgs:  # drop the AT 2024jfz-specific epoch fields
        for k in ("testable_epochs_mjd", "next_window_mjd"):
            r.pop(k, None)
        r["delay_yr"] = round(r["delay_days"] / DAYS_PER_YR, 3)
        r["sep_from_helios_arcsec"] = r.pop("sep_from_tns_arcsec")
        sep = SkyCoord(r["ra"], r["dec"], unit="deg").separation(centre).arcsec
        r["sep_centre_arcsec"] = round(float(sep), 3)
    ref = [r for r in imgs if r["reference_image"]]
    others = [r for r in imgs if not r["reference_image"]]
    sn_sep = res.get("tns_to_reference_arcsec", float("nan"))
    outcome = decide(bool(res.get("path_check_pass", False)), sn_sep, others)
    return {
        "model": name,
        "provenance": "model_prediction",
        "files_sha256": {k: v[1] for k, v in lc.MODELS[name]["files"].items()},
        "z_s": Z_S,
        "dls_ds": res["dls_ds"],
        "beta_xy": res["beta_xy"],
        "n_images": res["n_images"],
        "helios_to_reference_arcsec": sn_sep,
        "reference_magnification": ref[0]["magnification"] if ref else None,
        "days_per_arcsec2": res.get("days_per_arcsec2"),
        "path_check_max_rel_diff": res.get("path_check_max_rel_diff"),
        "path_check_pass": res.get("path_check_pass"),
        "outcome": outcome,
        "images": imgs,
        "_model": model,
    }


def plot(res: dict, out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.io import fits

    fig, ax = plt.subplots(figsize=(7, 7))
    kap = fits.getdata(lc.model_files(res["model"], ("kappa_map",))["kappa_map"]).astype(float)
    ax.imshow(np.log10(np.clip(kap, 1e-2, None)), origin="lower", cmap="Greys")
    model = res["_model"]
    for r in res["images"]:
        i, j = model.wcs.world_to_pixel_values(r["ra"], r["dec"])
        ax.plot(i, j, "*" if r["reference_image"] else "o", ms=12, mfc="none", mec="red")
        ax.annotate(
            f"{r['delay_yr']:+.2f} yr, μ {r['magnification']:+.1f}",
            (i, j),
            color="red",
            fontsize=8,
            xytext=(5, 5),
            textcoords="offset points",
        )
    ax.plot(*model.wcs.world_to_pixel_values(*HELIOS_RADEC), "x", ms=10, color="blue")
    for (ra, dec), _, _ in PUBLISHED_NEXT.values():
        ax.plot(*model.wcs.world_to_pixel_values(ra, dec), "+", ms=10, color="green")
    ax.set_title(f"{res['model']}: {res['outcome']} (z_s {Z_S})", fontsize=9)
    fig.suptitle(
        "SN Helios (blue x), RELICS images (red; * = at the SN), published next images (green +)",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(out / "images.png", dpi=80)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    res = run_model()
    plot(res, args.out)
    summary = {
        "provenance": "model_prediction",
        "source": "RELICS RXC J0018 Lenstool v1 map at the SN Helios position "
        "(arXiv:2609.30440v2); plan PLAN.md",
        "helios_radec": [round(v, 6) for v in HELIOS_RADEC],
        "published_next_centre_radec": [round(v, 6) for v in published_centre()],
        "assumptions": {
            "gate_s_arcsec": GATE_S_ARCSEC,
            "early_yr": EARLY_YR,
            "bright_ratio": BRIGHT_RATIO,
            "next_window_yr": NEXT_WINDOW_YR,
            "next_radius_arcsec": NEXT_RADIUS_ARCSEC,
            "path_tol": at.PATH_TOL,
            "cosmology": "FlatLambdaCDM H0 70, Om0 0.3 (map cosmology)",
        },
        "published_next": {
            k: {"radec": [round(c, 6) for c in v[0]], "delay_yr": v[1], "mu": v[2]}
            for k, v in PUBLISHED_NEXT.items()
        },
        "published_helios_mu": PUBLISHED_HELIOS_MU,
        "published_past_yr_ago": PUBLISHED_PAST_YR,
        "model": {k: v for k, v in res.items() if k != "_model"},
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    past = [r for r in res["images"] if not r["reference_image"] and r["delay_yr"] < 0]
    posthoc = {
        "provenance": "derived",
        "label": "post hoc (planned in PLAN.md 'Not done here'; no verdict depends on it)",
        "past_images_yr_ago": [round(-r["delay_yr"], 2) for r in past],
        "past_images_mu": [r["magnification"] for r in past],
        "published_past_yr_ago": PUBLISHED_PAST_YR,
    }
    (args.out / "posthoc.json").write_text(json.dumps(posthoc, indent=2) + "\n")
    print(
        res["outcome"],
        res["n_images"],
        res["helios_to_reference_arcsec"],
        res["path_check_max_rel_diff"],
    )
    for r in res["images"]:
        print(r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
