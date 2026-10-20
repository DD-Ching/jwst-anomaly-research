"""Star-ellipticity (PSF anisotropy) radial-pattern test for the Euclid Q1 shear null (D-067).

D-067's galaxy shapes are deconvolved for PSF *size* only, so PSF anisotropy leaks into them with
factor β = T_psf / T_int. This measures, on the D-067 trial-centre grid of each Deep Field, the
annulus mean radial component of the *stars'* ellipticities (the PSF pattern itself) and compares
β times its largest value with the injected radial shear averaged over the same annulus.
Frozen plan and decision rule: ``results/w5_shear/psf_radial/PLAN.md``.

Star rows: the same Q1 MER catalogue from ESA's Euclid archive (EAS TAP), one cone per Deep Field,
kept only on the field's D-067 tiles (``tileid`` = ``segmentation_map_id`` // 10⁶). IRSA answered
about one tile query per minute on 2026-10-20; EAS gave all three fields in ~8 s. Rows with a NULL
``point_like_prob`` pass EAS's ``> 0.9`` filter and are dropped (IRSA excludes them). Any IRSA star
tiles already cached are compared row by row (``irsa_identity``). Cached under the data root.
Outputs are ``derived``; β, R and the star-count floor are ASSUMPTIONs named below.

    python scripts/w5_euclid_psf_radial.py
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table, vstack
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exotic_screens as es  # noqa: E402
import w5_euclid_shear as w5  # noqa: E402

from jwst_anomaly import paths  # noqa: E402

BETAS = (1.0, 2.3)  # fiducial leakage T_psf/T_int, and its bound at the 1.2 x PSF cut (ASSUMPTION)
POST_HOC_MAX_ABS_E = (
    0.2  # post hoc star clip: rounder than this (binaries, blends removed; ASSUMPTION)
)
MIN_STARS = 5  # centres with fewer annulus stars are skipped (ASSUMPTION)
# Injected radial shear averaged over the annulus 1.5-3 θ_E (area-weighted <(θ_E/θ)²>,
# 2 ln 2 / 6.75), times the responsivity: the signal a flag needs (model_prediction, scale-free).
SIGNAL = w5.RESPONSIVITY * 2.0 * math.log(2.0) / (3.0**2 - 1.5**2)


EAS_TAP = "https://eas.esac.esa.int/tap-server/tap/sync"
# Cone radius per Deep Field: max tile-centre distance from w5.PILOTS + 0.36° (tile half-diagonal).
CONE_DEG = {"EDF-F": 3.0, "EDF-S": 4.7, "EDF-N": 3.8}


def eas_query(name: str) -> str:
    ra, dec = w5.PILOTS[name]
    sel = w5._selection(True)
    return (
        "SELECT segmentation_map_id, right_ascension AS ra, declination AS dec, ellipticity, "
        "position_angle, semimajor_axis, flux_detection_total, point_like_prob "
        "FROM catalogue.mer_catalogue WHERE "
        "1=CONTAINS(POINT('ICRS', right_ascension, declination), "
        f"CIRCLE('ICRS', {ra}, {dec}, {CONE_DEG[name]})) AND {sel}"
    )


def eas_stars(name: str, tiles, cache: Path) -> Table:
    """Stars of one Deep Field from EAS, on the field's D-067 tiles only (cached CSV)."""
    q = eas_query(name)
    t = eas_csv(q, cache / f"eas_stars_{name}_{w5._key(q)}.csv")
    t["tileid"] = np.asarray(t["segmentation_map_id"], np.int64) // 10**6
    keep = ~np.ma.getmaskarray(t["point_like_prob"]) & np.isin(t["tileid"], [x[0] for x in tiles])
    return t[keep]


COVERAGE_CELL_ARCMIN = 0.5  # galaxy count-map cell (server-side GROUP BY on EAS)


def coverage_query(name: str) -> str:
    """Counts of D-067 galaxies (VIS < 24.5, before the resolved cut) in 0.5' cells of one field."""
    ra, dec = w5.PILOTS[name]
    k = 60.0 / COVERAGE_CELL_ARCMIN
    kx = k * math.cos(math.radians(dec))
    sel = w5._selection(False)
    return (
        "SELECT iy, ix, COUNT(*) AS n FROM (SELECT "
        f"FLOOR(declination*{k:g}) AS iy, FLOOR(right_ascension*{kx:.4f}) AS ix "
        "FROM catalogue.mer_catalogue WHERE "
        "1=CONTAINS(POINT('ICRS', right_ascension, declination), "
        f"CIRCLE('ICRS', {ra}, {dec}, {CONE_DEG[name]})) AND {sel}) AS s GROUP BY iy, ix"
    )


def eas_csv(q: str, path: Path) -> Table:
    """EAS TAP sync query as CSV, cached at ``path`` (atomic write)."""
    if not path.exists():
        r = requests.get(
            EAS_TAP,
            params={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q},
            timeout=600,
        )
        r.raise_for_status()
        if r.text.lstrip().startswith("<"):
            raise RuntimeError(f"EAS TAP error: {r.text[:400]}")
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".part")
        tmp.write_text(r.text)
        tmp.replace(path)
    return Table.read(path, format="ascii.csv")


def galaxy_cells(name: str, cache: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Galaxy count-map cell centres (RA, Dec) and counts for the coverage cut."""
    q = coverage_query(name)
    t = eas_csv(q, cache / f"eas_cover_{name}_{w5._key(q)}.csv")
    k = 60.0 / COVERAGE_CELL_ARCMIN
    kx = k * math.cos(math.radians(w5.PILOTS[name][1]))
    ra = (np.asarray(t["ix"], float) + 0.5) / float(f"{kx:.4f}")
    dec = (np.asarray(t["iy"], float) + 0.5) / k
    return ra, dec, np.asarray(t["n"], float)


def annulus_counts(x, y, n, cx, cy, r_in: float, r_out: float) -> np.ndarray:
    """Sum of cell counts ``n`` at (x, y) with r_in <= distance < r_out from each centre."""
    tree = cKDTree(np.c_[x, y])
    out = np.zeros(len(cx))
    for i, (px, py) in enumerate(zip(cx, cy, strict=True)):
        rows = np.asarray(tree.query_ball_point([px, py], r_out), int)
        keep = np.hypot(x[rows] - px, y[rows] - py) >= r_in
        out[i] = n[rows][keep].sum()
    return out


def irsa_identity(stars: Table, cache: Path) -> dict:
    """Row-by-row check against IRSA star tiles cached by ``w5_euclid_shear.py`` (if any)."""
    n_irsa = n_match = n_eas = 0
    max_de = 0.0
    for f in sorted((cache / "tiles").glob("star_*.ecsv")):
        a = Table.read(f)
        b = stars[stars["tileid"] == int(f.name.split("_")[1])]
        n_irsa, n_eas = n_irsa + len(a), n_eas + len(b)
        if not (len(a) and len(b)):
            continue
        d, i = cKDTree(np.c_[b["ra"], b["dec"]]).query(np.c_[a["ra"], a["dec"]])
        ok = d < 1e-6  # deg
        n_match += int(ok.sum())
        de = np.abs(np.asarray(a["ellipticity"])[ok] - np.asarray(b["ellipticity"])[i[ok]])
        max_de = max(max_de, float(de.max(initial=0.0)))
    return {"n_irsa": n_irsa, "n_eas_same_tiles": n_eas, "n_matched": n_match, "max_abs_de": max_de}


def star_ellipticity(stars) -> np.ndarray:
    """Observed star ellipticity (PA east of north; no deconvolution), as ``w5.star_psf``."""
    a = np.asarray(stars["semimajor_axis"], float)
    ell = np.asarray(stars["ellipticity"], float)
    return es.intrinsic_ellipticity(
        a, a * (1.0 - ell), w5.pa_east_of_north(stars["position_angle"]), 0.0
    )


def annulus_radial(x, y, e, cx, cy, r_in: float, r_out: float, min_n: int = MIN_STARS) -> dict:
    """Per centre: mean radial component of ``e`` over sources in r_in <= r < r_out, its error
    (rms / sqrt(n)) and n. Radial is positive, the sign of ``w5.radial_shear``."""
    x, y, e = np.asarray(x, float), np.asarray(y, float), np.asarray(e, complex)
    tree = cKDTree(np.c_[x, y])
    mean, err, n = (np.full(len(cx), np.nan) for _ in range(3))
    for i, (px, py) in enumerate(zip(cx, cy, strict=True)):
        rows = np.asarray(tree.query_ball_point([px, py], r_out), int)
        de, dn = -(x[rows] - px), y[rows] - py
        keep = np.hypot(de, dn) >= r_in
        n[i] = keep.sum()
        if n[i] < min_n:
            continue
        er = np.real(e[rows][keep] * np.exp(-2j * np.arctan2(de[keep], dn[keep])))
        mean[i] = er.mean()
        err[i] = er.std(ddof=1) / math.sqrt(n[i])
    return {"mean": mean, "err": err, "n": n}


def decide(rho: dict[str, float]) -> str:
    """Frozen PLAN rule on ρ(β) = β max|ē_r| / SIGNAL, maximised over fields and θ_E."""
    if rho["2.3"] < 0.2:
        return "PSF gradients cannot fake or hide a flag"
    if rho["1"] < 1.0:
        return "PSF term not negligible at the edge of the resolved cut"
    return "PSF could fake or hide a flag"


def field_stats(name: str, stars: Table, tiles, cells) -> dict:
    ra, dec = np.asarray(stars["ra"], float), np.asarray(stars["dec"], float)
    e = star_ellipticity(stars)
    ok = np.isfinite(e)
    ra0 = float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians([t[1] for t in tiles]))))) % 360)
    dec0 = float(np.mean([t[2] for t in tiles]))
    x, y = w5.tangent_plane(ra[ok], dec[ok], ra0, dec0)
    e = e[ok]
    out = {
        "n_stars": int(ok.sum()),
        "n_tiles_with_stars": int(len(np.unique(stars["tileid"]))),
        "rms_abs_e": float(np.sqrt(np.mean(np.abs(e) ** 2))),
    }
    out["theta_e"] = {}
    for te in w5.SURVEY_THETA_E_ARCSEC:
        step = w5.GRID_STEP_OVER_THETA_E * te
        gx, gy = np.meshgrid(
            np.arange(x.min(), x.max() + step, step), np.arange(y.min(), y.max() + step, step)
        )
        # D-067 coverage cut: annulus galaxy count >= COVERAGE_MIN x the field median
        cra, cdec, cn = cells
        kx, ky = w5.tangent_plane(cra, cdec, ra0, dec0)
        ng = annulus_counts(kx, ky, cn, gx.ravel(), gy.ravel(), 1.5 * te, 3.0 * te)
        covered = ng >= w5.COVERAGE_MIN * np.median(ng[ng > 0])
        cx, cy = gx.ravel()[covered], gy.ravel()[covered]
        r = annulus_radial(x, y, e, cx, cy, 1.5 * te, 3.0 * te)
        v = np.isfinite(r["mean"])
        m, s = r["mean"][v], r["err"][v]
        i = int(np.argmax(np.abs(m)))
        z = m / s
        # Post hoc (not in the frozen plan): significance counts against a Gaussian, the typical
        # 95 % upper bound |ē_r| + 2σ, and the same max statistic on stars with |e| < 0.2.
        clip = np.abs(e) < POST_HOC_MAX_ABS_E
        rc = annulus_radial(x[clip], y[clip], e[clip], cx, cy, 1.5 * te, 3.0 * te)
        vc = np.isfinite(rc["mean"])
        mc, sc = rc["mean"][vc], rc["err"][vc]
        ic = int(np.argmax(np.abs(mc)))
        upper = np.abs(m) + 2.0 * s
        post_hoc = {
            "n_abs_z_gt_3": int(np.sum(np.abs(z) > 3)),
            "gauss_expect_abs_z_gt_3": float(0.0027 * len(z)),
            "max_abs_z": float(np.max(np.abs(z))),
            "median_upper95": float(np.median(upper)),
            "rho_median_upper95": {f"{b:g}": float(b * np.median(upper) / SIGNAL) for b in BETAS},
            "clipped": {
                "n_stars": int(clip.sum()),
                "max_abs_er": float(abs(mc[ic])),
                "max_abs_er_sigma": float(abs(mc[ic]) / sc[ic]),
                "max_abs_z": float(np.max(np.abs(mc / sc))),
                "chi2_over_dof": float(np.mean((mc / sc) ** 2)),
                "median_upper95": float(np.median(np.abs(mc) + 2.0 * sc)),
                "rho": {f"{b:g}": float(b * abs(mc[ic]) / SIGNAL) for b in BETAS},
            },
        }
        out["theta_e"][f"{te:g}"] = {
            "n_grid": int(gx.size),
            "n_covered": int(covered.sum()),
            "n_centres": int(v.sum()),
            "median_n_stars": float(np.median(r["n"][v])),
            "max_abs_er": float(abs(m[i])),
            "max_abs_er_sigma": float(abs(m[i]) / s[i]),
            "mean_er": float(np.mean(m)),
            "median_err": float(np.median(s)),
            # independent-centre approximation: annuli overlap, so this over-counts dof (label)
            "chi2": float(np.sum((m / s) ** 2)),
            "dof": int(v.sum()),
            "rho": {f"{b:g}": float(b * abs(m[i]) / SIGNAL) for b in BETAS},
            "post_hoc": post_hoc,
        }
        print(name, te, json.dumps(out["theta_e"][f"{te:g}"]), flush=True)
    return out


def main() -> int:
    cache = paths.data_root() / "euclid_q1_shear"
    tiles = w5.q1_tiles(cache)
    stars = {n: eas_stars(n, v, cache) for n, v in tiles.items()}
    res = {n: field_stats(n, stars[n], v, galaxy_cells(n, cache)) for n, v in tiles.items()}
    ident = irsa_identity(vstack(list(stars.values())), cache)
    rho = {
        f"{b:g}": max(t["rho"][f"{b:g}"] for f in res.values() for t in f["theta_e"].values())
        for b in BETAS
    }
    summary = {
        "provenance": {
            "stars": "observed (Euclid Q1 MER catalogue, ESA EAS TAP; SExtractor moments)",
            "service": f"{EAS_TAP} table catalogue.mer_catalogue",
            "statistics": "derived",
            "signal": "model_prediction (annulus-mean injected radial shear, R = 0.5)",
            "run_utc": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime()),
            "assumptions": {
                "betas": BETAS,
                "min_stars": MIN_STARS,
                "signal": SIGNAL,
                "star_selection": w5.STAR_SEL,
                "star_mags": w5.STAR_MAGS,
            },
        },
        "irsa_identity": ident,
        "fields": res,
        "rho_max": rho,
        "decision": decide(rho),
    }
    out = paths.repo_root() / "results" / "w5_shear" / "psf_radial" / "summary.json"
    out.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({"rho_max": rho, "decision": summary["decision"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
