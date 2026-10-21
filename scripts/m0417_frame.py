"""MACS0417: offset of the Jauzac et al. 2019 image list (RELICS HST frame, an ASSUMPTION tested
by gate F2) to the JWST 1208 F200W mosaic frame, then ``m0417_extra.run_field`` on the shifted list.

Frozen plan: ``results/lens_consistency/macs0417_frame/PLAN.md``. Inputs are the RELICS macs0417m11
HST catalogue and the JWST pipeline ``_cat.ecsv`` (sha256 in SOURCES.md and summary.json).

    python scripts/m0417_frame.py --relics <relics_cat.txt> [--cache data/cache/m0417x]

The JWST catalogue is ``<cache>/jw01208-o006_t003_nircam_clear-f200w_cat.ecsv``, the file
``m0417_extra.run_field`` reads.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import astropy.units as u
import lens_consistency as lc
import m0417_extra as mx
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

OUT = Path("results/lens_consistency/macs0417_frame")
# Frozen rules (PLAN.md; ASSUMPTIONs)
MAX_FWHM = 0.5  # arcsec, RELICS fwhm
MIN_NF5SIG = 3
PASS1, PASS2 = 1.0, 0.3  # arcsec match radii
F1_MIN_N, F1_MAX_ERR, F1_MAX_HALVES = 30, 0.02, 0.03
F2_MAX_MEDIAN = 0.18  # arcsec, 0.6 x the unshifted 0.300"


def read_relics(path: Path) -> Table:
    names = ["id", "ra", "dec", "x", "y", "fwhm", "area", "stel", "ell", "theta", "nf5sig", "nfobs"]
    rows = np.genfromtxt(path, comments="#", usecols=range(len(names)))
    t = Table(rows, names=names)
    return t[(t["fwhm"] <= MAX_FWHM) & (t["nf5sig"] >= MIN_NF5SIG)]


def offsets(
    relics: Table, jwst: SkyCoord, shift=(0.0, 0.0), radius=PASS1, unique: bool = False
) -> Table:
    """(ΔRA cos δ, ΔDec) in arcsec, JWST minus (RELICS + shift), nearest match within ``radius``."""
    r = SkyCoord(relics["ra"] * u.deg, relics["dec"] * u.deg).spherical_offsets_by(
        shift[0] * u.arcsec, shift[1] * u.arcsec
    )
    idx, sep, _ = r.match_to_catalog_sky(jwst)
    ok = sep.arcsec <= radius
    if unique:  # post hoc check: drop JWST sources matched by more than one RELICS source
        used, counts = np.unique(idx[ok], return_counts=True)
        ok &= np.isin(idx, used[counts == 1])
    dra, ddec = r[ok].spherical_offsets_to(jwst[idx[ok]])
    return Table(
        {
            "id": relics["id"][ok].astype(int),
            "dra": dra.to_value(u.arcsec) + shift[0],
            "ddec": ddec.to_value(u.arcsec) + shift[1],
        }
    )


def summarize(t: Table) -> dict:
    n = len(t)
    out = {"n": n}
    if n == 0:  # JSON-safe; gate_f1 fails on n first
        return {**out, "dra": None, "dra_err": None, "ddec": None, "ddec_err": None}
    for k in ("dra", "ddec"):
        v = np.asarray(t[k], float)
        med = float(np.median(v))
        mad = float(np.median(np.abs(v - med)))
        out[k] = med
        out[f"{k}_err"] = 1.2533 * 1.4826 * mad / np.sqrt(n)
    return out


def measure(relics: Table, jwst: SkyCoord) -> dict:
    p1 = summarize(offsets(relics, jwst))
    if p1["n"] == 0:
        return {"pass1": p1, "offset": p1, "halves_max_diff": None, "F1": False}
    t2 = offsets(relics, jwst, (p1["dra"], p1["ddec"]), PASS2)
    p2 = summarize(t2)
    even = summarize(t2[t2["id"] % 2 == 0])
    odd = summarize(t2[t2["id"] % 2 == 1])
    if p2["n"] == 0 or even["n"] == 0 or odd["n"] == 0:
        halves = None
    else:
        halves = max(abs(even["dra"] - odd["dra"]), abs(even["ddec"] - odd["ddec"]))
    f1 = gate_f1(p2, halves)
    unique = summarize(offsets(relics, jwst, (p1["dra"], p1["ddec"]), PASS2, unique=True))
    return {
        "pass1": p1,
        "offset": p2,
        "even": even,
        "odd": odd,
        "halves_max_diff": halves,
        "F1": f1,
        "posthoc_unique_matches": unique,
    }


def gate_f1(offset: dict, halves_max_diff: float) -> bool:
    """Frozen gate F1: N >= 30, error <= 0.02" per axis, halves within 0.03" per axis."""
    return bool(
        halves_max_diff is not None
        and offset["n"] >= F1_MIN_N
        and max(offset["dra_err"], offset["ddec_err"]) <= F1_MAX_ERR
        and halves_max_diff <= F1_MAX_HALVES
    )


def gate_f2(median_abs_offset: float | None) -> bool:
    """Frozen gate F2: median |peak - position| of the shifted list <= 0.18"."""
    return median_abs_offset is not None and bool(median_abs_offset <= F2_MAX_MEDIAN)


def shifted_loader(dra: float, ddec: float, base=None):
    """``base`` (default ``m0417_extra.load_images``) moved by (dRA cos dec, dDec) arcsec."""
    base = base or mx.load_images

    def load(field: str) -> Table:
        t = base(field)
        ra, dec = np.asarray(t["ra"], float), np.asarray(t["dec"], float)
        c = SkyCoord(ra * u.deg, dec * u.deg).spherical_offsets_by(dra * u.arcsec, ddec * u.arcsec)
        t["ra"], t["dec"] = c.ra.deg, c.dec.deg
        return t

    return load


def peak_offsets() -> dict:
    """``m0417_extra_offsets`` diagnostic on the (shifted) MACS0417 list."""
    _, _, obs, tile = mx.FIELDS["macs0417"]
    t = mx.load_images("macs0417")
    stamp, close = lc.image_stamper(
        mx.S3.format(obs, tile, "i2d.fits"), lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.1
    )
    d = []
    try:
        for r in t:
            s, _, xx, yy = stamp(float(r["ra"]), float(r["dec"]))
            d.append(float(np.hypot(*lc.peak_position(s, xx, yy))))
    finally:
        close()
    finite = [v for v in d if np.isfinite(v)]
    return {
        "n_images": len(d),
        "median_abs_offset_arcsec": float(np.median(finite)) if finite else None,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--relics", type=Path, required=True)
    ap.add_argument("--cache", type=Path, default=Path("data/cache/m0417x"))
    args = ap.parse_args(argv)
    relics = read_relics(args.relics)
    _, _, obs, tile = mx.FIELDS["macs0417"]
    # the same catalogue run_field reads
    jwst_cat = args.cache / Path(mx.S3.format(obs, tile, "cat.ecsv")).name
    jw = lc.load_shapes(jwst_cat)
    m = measure(
        relics, SkyCoord(np.asarray(jw["ra"], float) * u.deg, np.asarray(jw["dec"], float) * u.deg)
    )
    m["n_relics_selected"] = len(relics)
    summary = {
        "inputs": {
            "relics_sha256": hashlib.sha256(args.relics.read_bytes()).hexdigest(),
            "jwst_cat": jwst_cat.name,
            "jwst_cat_sha256": hashlib.sha256(jwst_cat.read_bytes()).hexdigest(),
        },
        "measurement": m,
    }
    if m["F1"]:
        saved = mx.load_images, mx.OUT
        mx.load_images = shifted_loader(m["offset"]["dra"], m["offset"]["ddec"])
        try:
            po = peak_offsets()
            po["F2"] = gate_f2(po["median_abs_offset_arcsec"])
            summary["shifted_peak_offsets"] = po
            if po["F2"]:
                mx.OUT = OUT
                summary["macs0417"] = mx.run_field("macs0417", args.cache)
        finally:
            mx.load_images, mx.OUT = saved
    summary["rules"] = {
        "MAX_FWHM": MAX_FWHM,
        "MIN_NF5SIG": MIN_NF5SIG,
        "PASS1": PASS1,
        "PASS2": PASS2,
        "F1": [F1_MIN_N, F1_MAX_ERR, F1_MAX_HALVES],
        "F2_MAX_MEDIAN": F2_MAX_MEDIAN,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=1)[:3000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
