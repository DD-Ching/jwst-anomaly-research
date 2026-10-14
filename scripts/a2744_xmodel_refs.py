"""Abell 2744: sibling reference fluxes for the extra images both lens models predict.

The D-041 cross-model run (``a2744_uncover.py --crossmodel``) left five extra images that
Bergamini+23 and UNCOVER v2.0 both predict but that neither forced run could test
(``no_reference``): every catalogued sibling is a fold-pair image with |mu| > 50 (the D-034 cap) or
fails a reference gate. 4.2c was tested separately (``a2744_uncover.py --check-42c``). This script
asks for 3.2, 26.2, 33.1 and 62.1 (and 4.2 for comparison), per model, without the |mu| cap:

- each catalogued sibling's peak flux (0.4" search, 0.2" aperture, 1" median-filtered stamp as the
  forced runs) divided by its |mu| and times |mu_target|, as an S/N against the target's own
  aperture error (``pred_snr_from``), and whether the sibling passes the other D-034 gates;
- ``pred_snr_pair``: the same from the summed sibling flux over the summed |mu| (a fold pair's
  total magnification is less sensitive to where the critical curve falls than either image's);
- the best aperture S/N within the system's forced search radius (``best_snr``), its flux over the
  pair prediction and ``lens_consistency.forced_class`` of the two;
- ``need_mu_factor``: by how much the pair's summed |mu| would have to be overestimated before the
  pair predicts the target at the D-029 ``absent`` floor (10 sigma);
- the nearest catalogued image of the same galaxy's other knots (``cat_sibling_knot``): an extra
  image that lands on a catalogued c image of another knot is that galaxy's third image, seen.

Targets are the ``both_loose`` rows of ``crossmodel.ecsv``. Inputs are the tracked
``images_predicted_*.ecsv`` and ``crossmodel.ecsv`` of the cross-model run (positions and mu are
``model_prediction``) and the Bergamini+23 ``arcs.dat`` positions (``observed``). Photometry is
``derived`` from S3 byte-range stamps of the UNCOVER/GLASS mosaics (no lens maps are needed).
Thresholds are ASSUMPTIONs.

Usage: python scripts/a2744_xmodel_refs.py [--out results/lens_consistency/abell2744-uncover2]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a2744_uncover as au  # noqa: E402
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

SYSTEMS = ("3.2", "4.2", "26.2", "33.1", "62.1")
MODELS = ("bergamini23", "uncover2")
BANDS = ("f150w", "f277w", "f444w")
ABSENT_SNR = 10.0  # D-029 floor for an absent image (ASSUMPTION)
KNOT_MATCH_ARCSEC = 1.5  # an extra this close to another knot's catalogued image is that galaxy's
DEFAULT_IN = Path("results/lens_consistency/abell2744-uncover2")


def galaxy(system: str) -> str:
    """Galaxy id of a knot system: '3.2' -> '3' (Lenstool knots of one source share it)."""
    return system.split(".")[0]


def predictions(sib: dict[str, tuple[float, float]], mu_t: float, e_t: float) -> dict:
    """Predicted target S/N from each sibling ``{id: (flux, |mu|)}`` and from the pair sum.

    Siblings with a non-finite or non-positive flux or |mu| are skipped (a negative peak flux is
    no reference). ``pair`` is NaN when no sibling qualifies or ``e_t`` is unusable."""
    ok = {
        k: (f, m)
        for k, (f, m) in sib.items()
        if np.isfinite(f) and f > 0 and np.isfinite(m) and m > 0
    }
    if not (np.isfinite(e_t) and e_t > 0 and np.isfinite(mu_t)) or not ok:
        return {"from": {}, "pair": np.nan, "pair_flux": np.nan}
    mu_t = abs(mu_t)
    per = {k: f / m * mu_t / e_t for k, (f, m) in ok.items()}
    pair_flux = sum(f for f, _ in ok.values()) / sum(m for _, m in ok.values()) * mu_t
    return {"from": per, "pair": pair_flux / e_t, "pair_flux": pair_flux}


def _r(v, nd: int = 2):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return round(v, nd) if np.isfinite(v) else None


def sibling_fluxes(cat: Table, stamp) -> dict:
    """Peak flux, S/N and D-034 gate result (without the |mu| cap) of each catalogued image."""
    out = {}
    for r in cat:
        img, err, xx, yy = stamp(float(r["ra"]), float(r["dec"]))
        f, e = lc.peak_flux(img, err, xx, yy)
        f0, e0 = lc.aperture_snr(img, err, xx, yy, 0.0, 0.0)
        px, py = lc.peak_position(img, xx, yy)
        fb, _ = lc.aperture_snr(img, err, xx, yy, px, py, r_ap=0.4)
        finite = all(np.isfinite(v) for v in (f, e, f0, e0)) and e > 0 and e0 > 0
        snr, snr0 = (f / e, f0 / e0) if finite else (np.nan, np.nan)
        conc = f / fb if np.isfinite(fb) and fb > 0 else np.nan
        out[str(r["image_id"])] = {
            "flux": f,
            "snr": snr,
            "conc": conc,
            "gates_ok": bool(
                finite
                and snr > 5
                and snr0 > 5
                and f <= lc.MAX_PEAK_OVER_POSITION * max(f0, 0)
                and conc >= lc.MIN_REF_CONCENTRATION
            ),
        }
    return out


def sibling_mu(pred: Table, system: str) -> dict[str, float]:
    """|mu| each model gives the catalogued images of ``system`` (``observed`` rows)."""
    sel = (np.asarray(pred["system"]).astype(str) == system) & (pred["image_class"] == "observed")
    return {str(r["matched_image"]): abs(float(r["magnification"])) for r in pred[sel]}


def targets(pred: Table, xm: Table, model: str, system: str) -> Table:
    """The system's ``no_reference`` extra images that the other model also predicts (the
    ``both_loose`` rows of the cross-model table, matched by position to 0.01")."""
    t = pred[np.asarray(pred["system"]).astype(str) == system]
    fc = np.array([au._col(t, "forced_class", k, "") for k in range(len(t))]).astype(str)
    x = xm[
        (np.asarray(xm["model"]).astype(str) == model)
        & (np.asarray(xm["system"]).astype(str) == system)
        & np.asarray(xm["both_loose"], bool)
    ]
    keep = [
        fc[k] == "no_reference"
        and len(x) > 0
        and float(np.min(au._sep_arcsec(x["ra"], x["dec"], t["ra"][k], t["dec"][k]))) < 0.01
        for k in range(len(t))
    ]
    return t[np.asarray(keep, bool)]


def knot_match(arcs: Table, system: str, ra: float, dec: float) -> tuple[str, float]:
    """Nearest catalogued image of the same galaxy's *other* knots (id, arcsec)."""
    sysid = np.asarray(arcs["system"]).astype(str)
    sel = np.flatnonzero(np.array([galaxy(s) == galaxy(system) for s in sysid]) & (sysid != system))
    if not len(sel):
        return "", np.nan
    d = au._sep_arcsec(np.asarray(arcs["ra"])[sel], np.asarray(arcs["dec"])[sel], ra, dec)
    k = int(np.argmin(d))
    return str(arcs["image_id"][sel[k]]), float(d[k])


def measure(
    pred: Table, xm: Table, model: str, arcs: Table, system: str, stamps: dict
) -> list[dict]:
    """One row per target extra image of ``system`` under one model."""
    cat = arcs[np.asarray(arcs["system"]).astype(str) == system]
    mus = sibling_mu(pred, system)
    rows = []
    for t in targets(pred, xm, model, system):
        ra, dec, mu_t = float(t["ra"]), float(t["dec"]), float(t["magnification"])
        search = float(t["search_arcsec"]) if np.isfinite(t["search_arcsec"]) else 1.0
        kid, ksep = knot_match(arcs, system, ra, dec)
        row = {
            "radec": [round(ra, 6), round(dec, 6)],
            "mu": _r(mu_t),
            "search_arcsec": _r(search),
            "sibling_mu": {k: _r(v, 1) for k, v in mus.items()},
            "cat_sibling_knot": kid,
            "cat_sibling_knot_sep_arcsec": _r(ksep),
        }
        for band, stamp in stamps.items():
            sib = sibling_fluxes(cat, stamp)
            img, err, xx, yy = stamp(ra, dec)
            _, e_t = lc.aperture_snr(img, err, xx, yy, 0.0, 0.0)
            p = predictions({k: (v["flux"], mus.get(k, np.nan)) for k, v in sib.items()}, mu_t, e_t)
            snr, flux, dx, dy = lc.best_within(img, err, xx, yy, search)
            ratio = flux / p["pair_flux"] if p["pair_flux"] > 0 else np.nan
            gated = [p["from"][k] for k in p["from"] if sib[k]["gates_ok"]]
            row[band] = {
                "siblings": {
                    k: {"flux": _r(v["flux"], 4), "snr": _r(v["snr"], 1), "gates_ok": v["gates_ok"]}
                    for k, v in sib.items()
                },
                "pred_snr_from": {k: _r(v, 1) for k, v in p["from"].items()},
                "pred_snr_pair": _r(p["pair"], 1),
                "pred_snr_gated_max": _r(max(gated), 1) if gated else None,
                "need_mu_factor": _r(ABSENT_SNR / p["pair"], 1) if p["pair"] > 0 else None,
                "best_snr": _r(snr, 1),
                "best_offset": [_r(dx, 1), _r(dy, 1)],
                "flux_over_pair_pred": _r(ratio),
                "forced_class": lc.forced_class(p["pair"], snr, ratio)
                if np.isfinite(ratio)
                else "no_reference",
            }
        rows.append(row)
    return rows


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--inputs", type=Path, default=DEFAULT_IN)
    ap.add_argument("--out", type=Path, default=DEFAULT_IN)
    args = ap.parse_args(argv)
    preds = {m: Table.read(args.inputs / f"images_predicted_{m}.ecsv") for m in MODELS}
    for t in preds.values():
        au.check_forced_table(t)
    xm = Table.read(args.inputs / "crossmodel.ecsv")
    _, files, _ = lc.load_model("abell2744-bergamini23")
    arcs = lensmodel.load_lenstool_images(files["arcs.dat"])
    half = lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.2
    opened = {b: lc.image_stamper(au.S3_I2D.format(b), half) for b in BANDS}
    stamps = {b: lc.halo_stamper(o[0], au.HALO_FILTER_ARCSEC) for b, o in opened.items()}
    try:
        runs = {s: {m: measure(preds[m], xm, m, arcs, s, stamps) for m in MODELS} for s in SYSTEMS}
    finally:
        for _, close in opened.values():
            close()
    summary = {
        "provenance": "positions and mu: model_prediction; catalogued positions: observed; "
        "fluxes: derived",
        "images": {b: au.S3_I2D.format(b) for b in BANDS},
        "assumptions": {
            "halo_filter_arcsec": au.HALO_FILTER_ARCSEC,
            "absent_snr": ABSENT_SNR,
            "knot_match_arcsec": KNOT_MATCH_ARCSEC,
            "search": "each row's forced search radius (lens_consistency.system_search_radius)",
            "mu_cap": "none (D-034's |mu| <= 50 cap dropped on purpose)",
        },
        "runs": runs,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "xmodel_refs.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return summary


if __name__ == "__main__":
    main()
