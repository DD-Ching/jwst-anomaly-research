"""Deblend-fragment check for newly admitted MIRI F770W singles (TASKS Then 7, D-014 addendum 4).

Plan: ``results/miri_confirm/PLAN_fragments.md`` (rule frozen before measurement). Reads the MAST L3
F770W ``_cat.ecsv`` and ``_segm.fits`` of SMACS 0723 (program 2736) and the uid lists in
``summary.json``; writes ``summary_fragments.json``.

Provenance: catalog and segmentation map are ``observed`` (pipeline products); adjacency, the
population test and verdicts are ``derived``; the 5x flux factor, 4-connectivity, S/N floor and the
p threshold are ASSUMPTIONs from the plan.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.stats import mannwhitneyu

PREFIX = "jw02736-o002_t001_miri_f770w"
NEW_TOP = (221, 239, 176, 118, 240)  # RESULT.md § Re-rank, new top-20 entries in rank order
FLUX_FACTOR = 5.0  # ASSUMPTION (plan)
MIN_SNR = 5.0  # ASSUMPTION (plan)
P_MAX = 0.01  # ASSUMPTION (plan)
MIN_GROUP = 10  # ASSUMPTION (plan): smallest group per flux tertile


def adjacent_labels(segm: np.ndarray) -> dict[int, set[int]]:
    """Map each segment label to its 4-connected neighbouring labels (0 = background)."""
    adj: dict[int, set[int]] = {}
    for a, b in ((segm[:, :-1], segm[:, 1:]), (segm[:-1, :], segm[1:, :])):
        m = (a != b) & (a > 0) & (b > 0)
        for x, y in {(int(p), int(q)) for p, q in zip(a[m], b[m], strict=True)}:
            adj.setdefault(x, set()).add(y)
            adj.setdefault(y, set()).add(x)
    return adj


def bright_adjacent(labels, flux, adj, factor: float = FLUX_FACTOR) -> np.ndarray:
    """True where a 4-adjacent segment has isophotal flux >= ``factor`` x the source's own."""
    fmap = dict(zip((int(v) for v in labels), (float(v) for v in flux), strict=True))
    out = np.zeros(len(labels), dtype=bool)
    for i, (lab, f) in enumerate(zip(labels, flux, strict=True)):
        nb = [fmap.get(n, np.nan) for n in adj.get(int(lab), ())]
        out[i] = bool(
            np.isfinite(f) and f > 0 and any(np.isfinite(v) and v >= factor * f for v in nb)
        )
    return out


def mechanism_test(ci, flux, flag, p_max: float = P_MAX, min_group: int = MIN_GROUP) -> dict:
    """Frozen rule: flagged sources have higher CI (MWU p < p_max) overall and per tertile."""

    def one(sel):
        a, b = ci[sel & flag], ci[sel & ~flag]
        if len(a) < min_group or len(b) < min_group:
            return {"n_adj": int(len(a)), "n_other": int(len(b)), "tested": False}
        p = float(mannwhitneyu(a, b, alternative="greater").pvalue)
        da, db = float(np.median(a)), float(np.median(b))
        return {
            "n_adj": int(len(a)),
            "n_other": int(len(b)),
            "tested": True,
            "median_adj": da,
            "median_other": db,
            "p": p,
            "pass": bool(p < p_max and da > db),
        }

    overall = one(np.ones(len(ci), dtype=bool))
    edges = np.quantile(flux, [1 / 3, 2 / 3])
    bins = np.digitize(flux, edges)
    tertiles = [one(bins == k) for k in range(3)]
    tested = [t for t in tertiles if t["tested"]]
    supported = bool(overall.get("pass") and tested and all(t["pass"] for t in tested))
    return {"overall": overall, "tertiles": tertiles, "supported": supported}


def annulus_overlap_table(ci, r_eq, r_in: float, ceiling: float) -> dict:
    """Post hoc: CI-at-ceiling counts for footprints reaching / not reaching the annulus."""
    reach = r_eq > r_in
    hi = ci >= ceiling
    return {
        "bkg_inner_radius_pix": round(r_in, 3),
        "ci_50_30_ceiling": round(ceiling, 3),
        "reach": {"n": int(reach.sum()), "ci_at_ceiling": int((reach & hi).sum())},
        "inside": {"n": int((~reach).sum()), "ci_at_ceiling": int((~reach & hi).sum())},
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("data_dir", type=Path, help="directory with the F770W _cat.ecsv and _segm.fits")
    p.add_argument("--summary", type=Path, default=Path("results/miri_confirm/summary.json"))
    p.add_argument("--out", type=Path, default=Path("results/miri_confirm"))
    a = p.parse_args(argv)
    cat = Table.read(a.data_dir / f"{PREFIX}_cat.ecsv")
    with fits.open(a.data_dir / f"{PREFIX}_segm.fits") as h:
        segm = np.asarray(h["SCI"].data if "SCI" in h else h[0].data, dtype=np.int64)
        cal_ver = h[0].header.get("CAL_VER")
    lab = np.asarray(cat["label"], int)
    flux = np.asarray(np.ma.filled(cat["isophotal_flux"].value, np.nan), float)
    ferr = np.asarray(np.ma.filled(cat["isophotal_flux_err"].value, np.nan), float)
    ci = np.asarray(np.ma.filled(cat["CI_50_30"], np.nan), float)
    flag = bright_adjacent(lab, flux, adjacent_labels(segm))

    with np.errstate(divide="ignore", invalid="ignore"):
        good = np.isfinite(ci) & np.isfinite(flux) & (ferr > 0) & (flux / ferr >= MIN_SNR)
    mech = mechanism_test(ci[good], flux[good], flag[good])

    groups = json.loads(a.summary.read_text())["groups"]

    def labels_of(key):
        uids = groups[key]["matched_uids"]
        return [int(u.rsplit("_", 1)[1]) for u in uids if u.startswith(PREFIX + "_")]

    idx = {int(v): i for i, v in enumerate(lab)}

    def rate(labels):
        f = [bool(flag[idx[v]]) for v in labels]
        return {
            "n": len(f),
            "bright_adjacent": int(sum(f)),
            "fraction": round(sum(f) / len(f), 3) if f else None,
        }

    base = labels_of(
        "passed_gate"
    )  # deviation: NIRCam-matched gate passers (all-113 list not stored)
    admitted = labels_of("single_band_only")
    rates = {
        "baseline_gate_matched_f770w": rate(base),
        "admitted_f770w": rate(admitted),
        "new_top20": rate(list(NEW_TOP)),
    }
    trigger = bool(
        rates["admitted_f770w"]["fraction"] > 0.5
        and rates["baseline_gate_matched_f770w"]["fraction"] < 0.25
    )

    per = []
    for v in NEW_TOP:
        i = idx[v]
        nb = int(cat["nn_label"][i])
        verdict = (
            "explained by deblending (processing)"
            if flag[i] and mech["supported"]
            else "to /vet-candidate"
        )
        per.append(
            {
                "label": v,
                "bright_adjacent": bool(flag[i]),
                "ci_50_30": round(float(ci[i]), 3),
                "isophotal_flux": float(flux[i]),
                "nn_label": nb,
                "nn_flux_ratio": round(float(flux[idx[nb]] / flux[i]), 2) if nb in idx else None,
                "nn_dist_pix": round(float(cat["nn_dist"][i]), 2),
                "verdict": verdict,
            }
        )

    # POST HOC (not in the plan): sources whose isophotal footprint reaches the local-background
    # annulus (equivalent radius > its inner radius) vs CI at or above the flat-profile ceiling.
    radii = np.asarray(cat.meta["aperture_params"]["aperture_radii"], float)
    r_in = float(cat.meta["aperture_params"]["bkg_aperture_inner_radius"])
    ceiling = float((radii[1] / radii[0]) ** 2)
    r_eq = np.sqrt(np.asarray(cat["isophotal_area"], float) / np.pi)
    post = annulus_overlap_table(ci[good], r_eq[good], r_in, ceiling)
    post["aperture_radii_pix"] = [round(float(v), 3) for v in radii]
    for d in per:
        i = idx[d["label"]]
        d["r_eq_pix"] = round(float(r_eq[i]), 2)
        d["reaches_bkg_annulus"] = bool(r_eq[i] > r_in)
        d["ci_at_ceiling"] = bool(ci[i] >= ceiling)

    out = {
        "products": PREFIX,
        "cal_ver": cal_ver,
        "n_catalog": len(cat),
        "n_tested": int(good.sum()),
        "rule": {
            "flux_factor": FLUX_FACTOR,
            "connectivity": 4,
            "min_snr": MIN_SNR,
            "p_max": P_MAX,
            "min_group": MIN_GROUP,
        },
        "mechanism": mech,
        "rates": rates,
        "flag_rule_trigger": trigger,
        "new_top20": per,
        "post_hoc_annulus": post,
    }
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "summary_fragments.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
