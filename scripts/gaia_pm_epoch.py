"""Does Gaia DR3 proper motion (epoch J2016.0) cost JWST star labels? (plan:
``results/gaia_pm_epoch/PLAN.md``; D-012 matches Gaia positions without propagation.)

One Gaia DR3 TAP cone around the JWST pipeline catalogue's sources; every Gaia source with a proper
motion is moved to the JWST epoch (linear, small-angle; parallax ignored) and the nearest catalogue
source to the J2016 and the propagated position is recorded.

    python scripts/gaia_pm_epoch.py CAT.ecsv --epoch 2022.43 --out results/gaia_pm_epoch
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import crossmatch

GAIA_EPOCH = 2016.0
# ASSUMPTION (D-012): radii of the classify rule (gaia_position) and of the query.
GAIA_RADIUS_ARCSEC = 0.3
QUERY_RADIUS_ARCSEC = 0.5
# ASSUMPTION: report every Gaia source moving more than this by the JWST epoch.
REPORT_SHIFT_ARCSEC = 0.15


def propagate(
    ra_deg: np.ndarray, dec_deg: np.ndarray, pmra: np.ndarray, pmdec: np.ndarray, dt_yr: float
) -> tuple[np.ndarray, np.ndarray]:
    """Linear proper-motion propagation; ``pmra`` is μ_α* (includes cos δ), mas/yr."""
    dra = np.asarray(pmra, float) * dt_yr / 3.6e6 / np.cos(np.deg2rad(dec_deg))
    ddec = np.asarray(pmdec, float) * dt_yr / 3.6e6
    return np.asarray(ra_deg, float) + dra, np.asarray(dec_deg, float) + ddec


def label_at(sep_arcsec: float) -> str:
    """Outcome of the Gaia-only D-012 rule for a Gaia source at this separation from its nearest
    JWST source (astrometric stars count within the query radius, position-only within 0.3″)."""
    if sep_arcsec <= GAIA_RADIUS_ARCSEC:
        return "position"
    if sep_arcsec <= QUERY_RADIUS_ARCSEC:
        return "astrometry_only"
    return "missed"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("catalog", type=Path)
    ap.add_argument("--epoch", type=float, required=True, help="JWST epoch (decimal year)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    cat = Table.read(args.catalog)
    jw = cat["sky_centroid"]
    ra_c, dec_c, r = crossmatch.field_circle(jw, QUERY_RADIUS_ARCSEC)
    import pyvo

    q = (
        "SELECT source_id, ra, dec, pmra, pmdec, pmra_error, pmdec_error, pmra_pmdec_corr, "
        "parallax, parallax_error, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE "
        + crossmatch._adql_cone(ra_c, dec_c, r)
    )
    gaia = pyvo.dal.TAPService(crossmatch.GAIA_TAP_URL).run_sync(q, maxrec=20000).to_table()
    dt = args.epoch - GAIA_EPOCH
    has_pm = ~np.ma.getmaskarray(gaia["pmra"]) & ~np.ma.getmaskarray(gaia["pmdec"])
    pmra = np.ma.filled(gaia["pmra"], 0.0).astype(float)
    pmdec = np.ma.filled(gaia["pmdec"], 0.0).astype(float)
    _, pm_snr = crossmatch.gaia_astrometric_snr(
        *(
            np.ma.filled(gaia[c], np.nan).astype(float)
            for c in (
                "parallax",
                "parallax_error",
                "pmra",
                "pmra_error",
                "pmdec",
                "pmdec_error",
                "pmra_pmdec_corr",
            )
        )
    )
    shift = np.hypot(pmra, pmdec) * dt / 1000.0
    ra1, dec1 = propagate(gaia["ra"], gaia["dec"], pmra, pmdec, dt)
    g0 = SkyCoord(gaia["ra"], gaia["dec"], unit="deg")
    g1 = SkyCoord(ra1, dec1, unit="deg")
    i0, s0, _ = g0.match_to_catalog_sky(jw)
    i1, s1, _ = g1.match_to_catalog_sky(jw)
    s0 = s0.arcsec
    s1 = s1.arcsec

    rows = []
    for k in np.flatnonzero(has_pm & (shift > REPORT_SHIFT_ARCSEC)):
        rows.append(
            {
                "source_id": int(gaia["source_id"][k]),
                "G": float(gaia["phot_g_mean_mag"][k]),
                "pm_masyr": float(np.hypot(pmra[k], pmdec[k])),
                "pm_snr": float(pm_snr[k]),
                "astrometric_star": bool(pm_snr[k] >= crossmatch.GAIA_MIN_PM_SNR),
                "shift_arcsec": float(shift[k]),
                "jwst_id_j2016": int(cat["label"][i0[k]]),
                "sep_j2016": float(s0[k]),
                "jwst_id_prop": int(cat["label"][i1[k]]),
                "sep_prop": float(s1[k]),
                "rule_j2016": label_at(s0[k]),
                "rule_prop": label_at(s1[k]),
            }
        )
    slow = has_pm & (shift < 0.05) & (s0 < QUERY_RADIUS_ARCSEC)
    summary = {
        "catalog": args.catalog.name,
        "catalog_sha256": hashlib.sha256(args.catalog.read_bytes()).hexdigest(),
        "jwst_epoch": args.epoch,
        "dt_yr": dt,
        "field_circle_deg": [ra_c, dec_c, r],
        "n_gaia": int(len(gaia)),
        "n_with_pm": int(has_pm.sum()),
        "n_shift_gt_0.15": int((has_pm & (shift > 0.15)).sum()),
        "n_shift_gt_0.3": int((has_pm & (shift > GAIA_RADIUS_ARCSEC)).sum()),
        "n_shift_gt_0.5": int((has_pm & (shift > QUERY_RADIUS_ARCSEC)).sum()),
        "slow_matched_n": int(slow.sum()),
        "slow_matched_median_sep_arcsec": float(np.median(s0[slow])) if slow.any() else None,
        "changed_outcome": [r for r in rows if r["rule_j2016"] != r["rule_prop"]],
        "provenance": "derived (Gaia DR3 observed astrometry, linear propagation)",
    }
    args.out.mkdir(parents=True, exist_ok=True)
    if rows:
        t = Table(rows=rows)
        t.meta.update(
            {"provenance": "derived", "source": "Gaia DR3 TAP cone; JWST pipeline catalogue"}
        )
        t.write(args.out / "fast_movers.ecsv", overwrite=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "changed_outcome"}, indent=1))
    for r_ in rows:
        print(r_)


if __name__ == "__main__":
    main()
