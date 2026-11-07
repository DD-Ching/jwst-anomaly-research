"""SMACS 0723 ``f200w_1096``: Mahler+2022 ICLv2 shear, magnification and orientation test.

Frozen plan: ``results/smacs_1096_shear/PLAN.md``. Evaluates the pinned ICLv2 model at the
catalogued position of record ``jw02736-o001_t001_nircam_f200w_1096`` for its NED redshift
(1.3618, observed) and z_s = 1, 2, 4, and compares the catalog major-axis PA with the predicted
stretch axis. Shape-noise probability from ``lens_consistency.anti_probability`` (D-024 addendum).

    python scripts/smacs_1096_shear.py --cat-dir <dir with the 2736 o001 CATs> \\
        [--chain-samples 200] [--out results/smacs_1096_shear]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import cutouts, lensmodel  # noqa: E402

RA, DEC = 110.846894, -73.442973  # record position (results/smacs_jwst3/candidates.csv)
Z_NED = 1.3618  # NED redshift of SMACS J0723:[NDA2023] 00908 (observed; candidate record)
Z_TEST = (1.0, 2.0, 4.0)
BANDS = ("f150w", "f200w", "f444w")
CAT = "jw02736-o001_t001_nircam_clear-{b}_cat.ecsv"
S3 = "s3://stpubdata/jwst/public/jw02736/L3/t/o001/"
SHEET_ARCSEC = 6.0
MATCH_ARCSEC = 0.3
# Frozen thresholds (ASSUMPTIONs; PLAN.md).
P_INCONSISTENT = 0.01
G_NO_POWER = 0.05
SIGMA_EPS = 0.26
MIN_ELLIPTICITY = 0.5


def decide(g: float, p: float) -> str:
    """Frozen rule: ``no_power`` if g < G_NO_POWER; else ``inconsistent`` if p < P_INCONSISTENT;
    else ``consistent``. ``g`` is the reduced shear at Z_NED, ``p`` = P(offset >= observed)."""
    if not np.isfinite(g) or g < G_NO_POWER:
        return "no_power"
    return "inconsistent" if p < P_INCONSISTENT else "consistent"


def offset_deg(a: float, b: float) -> float:
    return float(lensmodel.axis_offset_deg(a, b))


def sheet(path: Path, work: Path, shapes: dict, pa_pred: float) -> None:
    """Contact sheet: one north-up-agnostic panel per band (pixel grid, WCS-correct PA lines).

    Red: observed major axis (catalog PA, E of N). Cyan: predicted stretch axis at Z_NED."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    target = Table({"source_uid": ["f200w_1096"], "ra": [RA], "dec": [DEC]})
    fig, axes = plt.subplots(1, len(BANDS), figsize=(4 * len(BANDS), 4.3))
    for ax, b in zip(axes, BANDS, strict=True):
        uri = S3 + CAT.format(b=b).replace("_cat.ecsv", "_i2d.fits")
        tab = cutouts.make_cutouts(uri, target, SHEET_ARCSEC, work, storage_options={"anon": True})
        with fits.open(tab["path"][0]) as h:
            data, wcs = h["SCI"].data.astype(float), WCS(h["SCI"].header)
        lo, hi = np.nanpercentile(data, [5, 99.5])
        ax.imshow(np.clip(data, lo, hi), origin="lower", cmap="gray_r")
        c = SkyCoord(RA, DEC, unit="deg")
        x0, y0 = wcs.world_to_pixel(c)
        for pa, col, lab in (
            ((shapes.get(b, {}) or {}).get("pa_obs_deg"), "r", "observed"),
            (
                pa_pred,
                "c",
                "ICLv2 stretch",
            ),
        ):
            if pa is None:
                continue
            end = c.directional_offset_by(np.deg2rad(pa) * u.rad, 1.5 * u.arcsec)
            x1, y1 = wcs.world_to_pixel(end)
            ax.plot([2 * x0 - x1, x1], [2 * y0 - y1, y1], col, lw=1.2, label=lab)
        n = c.directional_offset_by(0 * u.rad, 1.0 * u.arcsec)
        xn, yn = wcs.world_to_pixel(n)
        ax.annotate("N", (xn, yn), color="y")
        ax.set_title(f'{b}  (6" cutout)')
        ax.legend(loc="lower left", fontsize=7)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("SMACS 0723 f200w_1096: catalog PA (red) vs ICLv2 predicted stretch (cyan)")
    fig.tight_layout()
    fig.savefig(path, dpi=90)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cat-dir", type=Path, required=True)
    ap.add_argument("--chain-samples", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, default=Path("results/smacs_1096_shear"))
    ap.add_argument("--work", type=Path, help="scratch dir for cutouts; enables the contact sheet")
    args = ap.parse_args()

    name = "smacs0723-iclv2"
    model, _, par = lc.load_model(name)
    lc.apply_frame_offset(name, model)
    zs = (Z_NED, *Z_TEST)
    pred = model.evaluate(np.full(len(zs), RA), np.full(len(zs), DEC), np.array(zs))
    r = float(np.hypot(pred["x"][0], pred["y"][0]))

    target = SkyCoord(RA, DEC, unit="deg")
    shapes = {}
    for b in BANDS:
        s = lc.load_shapes(args.cat_dir / CAT.format(b=b))
        sep = SkyCoord(s["ra"], s["dec"], unit="deg").separation(target).arcsec
        i = int(np.argmin(sep))
        if sep[i] > MATCH_ARCSEC:
            shapes[b] = None
            continue
        shapes[b] = {
            "label": int(s["label"][i]),
            "sep_arcsec": round(float(sep[i]), 3),
            "pa_obs_deg": round(float(s["pa_obs"][i]), 2),
            "ellipticity": round(float(s["ellipticity"][i]), 3),
            "semimajor_px": round(float(s["semimajor_px"][i]), 2),
            "snr": round(float(s["snr"][i]), 1),
        }

    g_ned = float(pred["reduced_shear"][0])
    pa_pred = float(pred["tangential_pa"][0])
    per_band = {}
    for b, sh in shapes.items():
        if sh is None:
            continue
        d = offset_deg(sh["pa_obs_deg"], pa_pred)
        p = lc.anti_probability(g_ned, SIGMA_EPS, MIN_ELLIPTICITY, anti_deg=d)
        per_band[b] = {"offset_deg": round(d, 2), "p_offset_ge": round(p, 4)}
        per_band[b]["p_offset_ge_shear0"] = round(
            lc.anti_probability(0.0, SIGMA_EPS, MIN_ELLIPTICITY, anti_deg=d), 4
        )

    # Model uncertainty: predicted PA and g at Z_NED over random bayes.dat rows, if pinned.
    chain = {"note": "no pinned bayes.dat for this model; best model only"}
    if "bayes" in lc.MODELS[name]:
        bayes, _, _, reference = lc.load_chain(name, par)
        rng = np.random.default_rng(args.seed)
        picks = rng.choice(len(bayes), min(args.chain_samples, len(bayes)), replace=False)
        pa_c, g_c, mu_c = [], [], []
        for row in picks:
            pr = lensmodel.posterior_par(par, bayes, int(row), reference)
            m = lensmodel.LensModel.from_par(pr)
            lc.apply_frame_offset(name, m)
            t = m.evaluate([RA], [DEC], [Z_NED])
            pa_c.append(float(t["tangential_pa"][0]))
            g_c.append(float(t["reduced_shear"][0]))
            mu_c.append(float(t["magnification"][0]))
        dpa = [offset_deg(p, pa_pred) for p in pa_c]
        chain = {
            "n_rows": len(picks),
            "seed": args.seed,
            "reduced_shear_p16_50_84": lc._pct(g_c),
            "magnification_p16_50_84": lc._pct(mu_c),
            "pa_offset_from_best_p50_p95_deg": [
                round(float(np.percentile(dpa, 50)), 2),
                round(float(np.percentile(dpa, 95)), 2),
            ],
        }

    primary = per_band.get("f200w")
    outcome = decide(g_ned, primary["p_offset_ge"]) if primary else "undetermined"
    summary = {
        "plan": "results/smacs_1096_shear/PLAN.md",
        "model": name,
        "record": "jw02736-o001_t001_nircam_f200w_1096",
        "ra": RA,
        "dec": DEC,
        "radius_arcsec": round(r, 2),
        "prediction": [
            {
                "z_s": float(pred["z_s"][k]),
                "kappa": round(float(pred["kappa"][k]), 4),
                "gamma": round(float(pred["gamma"][k]), 4),
                "reduced_shear": round(float(pred["reduced_shear"][k]), 4),
                "magnification": round(float(pred["magnification"][k]), 3),
                "tangential_pa_deg": round(float(pred["tangential_pa"][k]), 2),
                "provenance": "model_prediction",
            }
            for k in range(len(zs))
        ],
        "shapes": shapes,
        "shapes_provenance": "observed (jwst pipeline catalog moments)",
        "offsets": per_band,
        "chain": chain,
        "params": {
            "P_INCONSISTENT": P_INCONSISTENT,
            "G_NO_POWER": G_NO_POWER,
            "SIGMA_EPS": SIGMA_EPS,
            "MIN_ELLIPTICITY": MIN_ELLIPTICITY,
            "MATCH_ARCSEC": MATCH_ARCSEC,
            "provenance": "assumption",
        },
        "outcome_f200w": outcome,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if args.work:
        sheet(args.out / "contact_sheet.png", args.work, shapes, pa_pred)


if __name__ == "__main__":
    main()
