"""MACS0647 g85: host redshift from the DJA v7.0 photo-z and DJA NIRSpec v4.6 catalogues.

Frozen plan: ``results/transients/macs0647_1433/host/PLAN.md``. Rules H1-H3 are the functions below.
Both inputs are streamed into the scratchpad and deleted after the run (cloud disk rule):

    curl -sS <PHOTOZ url> | tar -xzO macs0647-grizli-v7.0-fix.eazypy.zout.fits > zout.fits
    curl -sS <SPECZ url> | gunzip -c > specz.csv
    python scripts/macs0647_g85_host.py --zout zout.fits --specz specz.csv

"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

OUT = Path("results/transients/macs0647_1433/host")
G85 = (101.99665, 70.250721)  # ICRS, pipeline frame (docs/candidates/macs0647-1433-g85.md)
M_OBS = 26.28  # AB, 0.15" aperture, not aperture-corrected (vetting record)
PHOTOZ = {
    "url": "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/macs0647-grizli-v7.0-fix.photoz.tar.gz",
    "member": "macs0647-grizli-v7.0-fix.eazypy.zout.fits",
    "member_sha256": "ceb199c61f1ca67df6afc34df6ed7ad61267f59ccaf27e9a2fd1e64c2662c04e",
}
SPECZ = {
    "url": "https://s3.amazonaws.com/msaexp-nirspec/extractions/dja_msaexp_emission_lines_v4.6.csv.gz",
    "sha256_gz": "ee64bb460f417d3fcf96fc53818c06ca0e27d784833063ae8eb48feb18c6ddf3",
}
# Frozen rules (PLAN.md; ASSUMPTIONs)
HOST_RADIUS_ARCSEC = 0.5
NEIGHBOUR_RADIUS_ARCSEC = 1.5
SPEC_RADIUS_ARCSEC = 1.0
SPEC_GRADE = 3
Z_CLUSTER = 0.591
MEMBER_WINDOW = 0.1  # in units of (1 + z_cl)
UNCONSTRAINED_WIDTH = 1.5  # z840 - z160
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)


def sep_arcsec(ra, dec, ra0: float = G85[0], dec0: float = G85[1]) -> np.ndarray:
    """Small-angle separation in arcsec (valid well below a degree)."""
    ra = np.asarray(ra, float)
    dec = np.asarray(dec, float)
    dra = ((ra - ra0 + 180.0) % 360.0 - 180.0) * np.cos(np.radians(dec0))
    return np.hypot(dra, dec - dec0) * 3600.0


def select_host(sep: np.ndarray) -> int | None:
    """H1: index of the nearest row within HOST_RADIUS_ARCSEC, or None."""
    sep = np.asarray(sep, float)
    if sep.size == 0 or not np.any(np.isfinite(sep)):
        return None
    i = int(np.nanargmin(sep))
    return i if sep[i] <= HOST_RADIUS_ARCSEC else None


def _member_class(z: float) -> str:
    if abs(z - Z_CLUSTER) <= MEMBER_WINDOW * (1 + Z_CLUSTER):
        return "member"
    return "foreground" if z < Z_CLUSTER else "background"


def classify_photoz(z160: float, z500: float, z840: float) -> str:
    """H2: class from the photo-z percentiles."""
    if not all(np.isfinite([z160, z500, z840])) or z840 - z160 > UNCONSTRAINED_WIDTH:
        return "unconstrained"
    return _member_class(z500)


def spec_class(rows: list[dict]) -> tuple[str, float | None]:
    """H3: class and median z of the grade-3 NIRSpec rows (already within SPEC_RADIUS_ARCSEC)."""
    z = [float(r["z_best"]) for r in rows if _grade(r) == SPEC_GRADE and _finite(r.get("z_best"))]
    if not z:
        return "none", None
    zm = float(np.median(z))
    return _member_class(zm), zm


def _grade(r: dict) -> int | None:
    try:
        return int(float(r.get("grade", "")))
    except ValueError:
        return None


def _finite(v) -> bool:
    try:
        return bool(np.isfinite(float(v)))
    except (TypeError, ValueError):
        return False


def abs_mag(m: float, z: float) -> float:
    """Derived: M = m - DM(z) + 2.5 log10(1 + z) (no K-correction, mu = 1; ASSUMPTIONs)."""
    return float(m - COSMO.distmod(z).to_value(u.mag) + 2.5 * np.log10(1 + z))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _specz_near(path: Path) -> list[dict]:
    out = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if not (_finite(r.get("ra")) and _finite(r.get("dec"))):
                continue
            s = float(sep_arcsec(float(r["ra"]), float(r["dec"])))
            if s <= SPEC_RADIUS_ARCSEC:
                keep = {
                    k: r.get(k) for k in ("file", "srcid", "ra", "dec", "grade", "z_best") if k in r
                }
                keep["sep_arcsec"] = round(s, 3)
                out.append(keep)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zout", type=Path, required=True)
    ap.add_argument("--specz", type=Path, required=True, help="the gunzipped DJA v4.6 csv")
    args = ap.parse_args(argv)

    sha = _sha256(args.zout)
    if sha != PHOTOZ["member_sha256"]:
        raise SystemExit(f"zout sha256 {sha} != pinned {PHOTOZ['member_sha256']}")
    t = Table.read(args.zout)
    sep = sep_arcsec(t["ra"], t["dec"])
    order = np.argsort(sep)
    cols = [
        c
        for c in ("id", "ra", "dec", "z_phot", "z160", "z500", "z840", "z_spec")
        if c in t.colnames
    ]

    def row(i: int) -> dict:
        d = {c: (float(t[c][i]) if c != "id" else int(t[c][i])) for c in cols}
        d["sep_arcsec"] = round(float(sep[i]), 3)
        return d

    neighbours = [row(int(i)) for i in order if sep[i] <= NEIGHBOUR_RADIUS_ARCSEC]
    ih = select_host(sep)
    host = row(ih) if ih is not None else None
    photoz_class = (
        classify_photoz(host["z160"], host["z500"], host["z840"]) if host is not None else "no_host"
    )
    spec_rows = _specz_near(args.specz)
    sclass, zspec = spec_class(spec_rows)
    outcome = sclass if sclass != "none" else photoz_class
    m_abs = None
    if host is not None and np.isfinite(host["z500"]) and host["z500"] > 0:
        m_abs = round(abs_mag(M_OBS, host["z500"]), 2)
    summary = {
        "g85": {"ra": G85[0], "dec": G85[1], "m_obs_ab": M_OBS},
        "inputs": {"photoz": PHOTOZ, "specz": SPECZ, "zout_rows": len(t)},
        "rules": {
            "host_radius_arcsec": HOST_RADIUS_ARCSEC,
            "neighbour_radius_arcsec": NEIGHBOUR_RADIUS_ARCSEC,
            "spec_radius_arcsec": SPEC_RADIUS_ARCSEC,
            "spec_grade": SPEC_GRADE,
            "z_cluster": Z_CLUSTER,
            "member_window_1pz": MEMBER_WINDOW,
            "unconstrained_width": UNCONSTRAINED_WIDTH,
        },
        "host": host,
        "neighbours": neighbours,
        "photoz_class": photoz_class,
        "specz_rows": spec_rows,
        "specz_class": sclass,
        "z_spec": zspec,
        "outcome": outcome,
        "abs_mag_at_z500_mu1_noK": m_abs,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {k: summary[k] for k in ("host", "photoz_class", "specz_class", "outcome")}, indent=1
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
