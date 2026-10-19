"""AT 2024jfz (SMACS 0723): host redshift from Caminha+22 MUSE and the DJA v7.4 eazy photo-z.

Plan: ``results/transients/smacs_at2024jfz/host_z/PLAN.md`` (frozen before any catalogue row near
the position was read). Candidates are DJA sources within ``R_MAX`` of the TNS position; the host
is the one with the smallest d / r_half; a trusted MUSE z_spec within ``MUSE_MATCH`` of the host
overrides its photo-z; :func:`classify` maps the host's z interval to the counter-image windows of
#360. Catalogue values are ``observed`` (z_spec) or ``derived`` (photo-z); radii ASSUMPTIONs.

    python scripts/smacs_at2024jfz_hostz.py --out results/transients/smacs_at2024jfz/host_z
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import offplane_potentials as op  # noqa: E402
import w1_stage_inputs as w1  # noqa: E402
from smacs_at2024jfz import TNS_RADEC  # noqa: E402

from jwst_anomaly import paths  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

PHOT = (
    "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/smacs0723-grizli-v7.4-fix_phot.fits",
    "361a4eb944fefa3afbbbfe08192d2d5c5764ae7137e41d77c2ce63c5b4a47166",
)
ZOUT = "smacs0723-grizli-v7.4-fix.eazypy.zout.fits"
R_MAX = 1.0  # arcsec, ASSUMPTION: host search radius
MUSE_MATCH = 0.5  # arcsec, ASSUMPTION: MUSE-to-host match
PIX = 0.04  # arcsec per DJA mosaic pixel (flux_radius unit)
R_FLOOR = 0.1  # arcsec, ASSUMPTION: floor on r_half
DZ_SPEC = 0.01  # ASSUMPTION: z_spec interval half-width
Z_MEMBER = 0.42  # upper edge for "member or foreground" (cluster z 0.388)
TESTED = ((1.35, 1.50), (1.75, 2.50))  # #360: counter-images absent there
UNDETERMINED = (1.55, 1.70)  # #360: undetermined


def classify(lo: float | None, hi: float | None) -> str:
    """Frozen outcome for a host z interval [lo, hi]; ``None`` means no host within R_MAX."""
    if lo is None or hi is None or not (np.isfinite(lo) and np.isfinite(hi)):
        return "untested"
    if hi <= Z_MEMBER:
        return "member_or_fg"
    if any(a <= lo and hi <= b for a, b in TESTED):
        return "tested_window"
    if lo <= UNDETERMINED[1] and hi >= UNDETERMINED[0]:
        return "undetermined_window"
    return "untested"


def pick_host(sep_arcsec: np.ndarray, r_half_arcsec: np.ndarray) -> int | None:
    """Index of the candidate with the smallest d / max(r_half, R_FLOOR) within R_MAX."""
    ok = sep_arcsec <= R_MAX
    if not ok.any():
        return None
    norm = np.where(ok, sep_arcsec / np.maximum(np.nan_to_num(r_half_arcsec), R_FLOOR), np.inf)
    return int(np.argmin(norm))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    target = paths.data_root() / "cache/dja/photoz-v7.4" / ZOUT
    log = w1.stage(target, w1.PHOTOZ[ZOUT])
    if log["status"].startswith("FAILED"):
        raise SystemExit(log["status"])
    zout = Table.read(target)
    phot = Table.read(fetch_catalog(*PHOT), hdu=1)
    if not np.array_equal(np.asarray(phot["id"]), np.asarray(zout["id"])):
        raise SystemExit("DJA phot and zout rows differ")
    spec = op.MUSE["smacs0723-iclv2"]
    muse = op.read_muse(fetch_catalog(spec["url"], spec["sha256"]), spec)

    tns = SkyCoord(*TNS_RADEC, unit="deg")
    sky = SkyCoord(zout["ra"], zout["dec"], unit="deg")
    sep = tns.separation(sky).arcsec
    rh = np.asarray(phot["flux_radius"], float) * PIX
    near = np.flatnonzero(sep <= R_MAX)
    rows = [
        {
            "id": int(zout["id"][i]),
            "sep_arcsec": round(float(sep[i]), 3),
            "r_half_arcsec": round(float(rh[i]), 3),
            "d_over_rhalf": round(float(sep[i] / max(rh[i], R_FLOOR)), 2),
            "z_phot": round(float(zout["z_phot"][i]), 3),
            "z160": round(float(zout["z160"][i]), 3),
            "z840": round(float(zout["z840"][i]), 3),
            "z_spec_dja": round(float(zout["z_spec"][i]), 4),
            "nusefilt": int(zout["nusefilt"][i]),
        }
        for i in near[np.argsort(sep[near])]
    ]
    k = pick_host(sep, rh)
    host, lo, hi, zsrc = None, None, None, None
    muse_rows = []
    msky = SkyCoord(muse["RA"], muse["DEC"], unit="deg")
    msep_tns = tns.separation(msky).arcsec
    for j in np.flatnonzero(msep_tns <= 3.0):  # context only: MUSE sources within 3"
        muse_rows.append(
            {
                "source": str(muse["SOURCE"][j]),
                "sep_tns_arcsec": round(float(msep_tns[j]), 3),
                "z_spec": float(muse["Z_SPEC"][j]),
                "trusted": bool(muse["Z_OK"][j]),
            }
        )
    if k is not None:
        host = rows[[r["id"] for r in rows].index(int(zout["id"][k]))]
        msep = sky[k].separation(msky).arcsec
        m = np.flatnonzero((msep <= MUSE_MATCH) & np.asarray(muse["Z_OK"]) & (muse["Z_SPEC"] > 0))
        if m.size:
            j = m[np.argmin(msep[m])]
            z = float(muse["Z_SPEC"][j])
            lo, hi, zsrc = z - DZ_SPEC, z + DZ_SPEC, f'MUSE {muse["SOURCE"][j]} ({msep[j]:.2f}")'
        else:
            lo, hi, zsrc = host["z160"], host["z840"], "DJA eazy z160-z840"
    res = {
        "provenance": "derived (z_spec observed; photo-z derived; radii ASSUMPTION)",
        "tns_radec": list(TNS_RADEC),
        "rules": {"R_MAX": R_MAX, "MUSE_MATCH": MUSE_MATCH, "R_FLOOR": R_FLOOR, "DZ_SPEC": DZ_SPEC},
        "sources": {"dja_phot_sha256": PHOT[1], "zout": log, "muse": spec["ref"]},
        "candidates": rows,
        "muse_within_3arcsec": muse_rows,
        "host": host,
        "host_z_interval": None if lo is None else [round(lo, 3), round(hi, 3)],
        "host_z_source": zsrc,
        "outcome": classify(lo, hi),
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "result.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("host", "host_z_interval", "host_z_source", "outcome")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
