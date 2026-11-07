"""MIRI single-band detections vs NIRCam counterparts in DJA v7.4 (TASKS Then 7, D-014).

Plan: ``results/miri_confirm/PLAN.md`` (rule frozen before measurement). Reads one pipeline run's
``smacs0723_miri`` ``sources`` and ``quality`` tables and the DJA v7.4 SMACS ``*_phot.fits``.
A source is NIRCam-confirmed when a DJA source with S/N >= 5 in at least one NIRCam band (0.5"
aperture) lies within ``MATCH_ARCSEC``. The chance rate comes from the same match at shifted
positions. Writes ``summary.json``.

Provenance: catalogs are ``observed``; matches, rates and the verdict are ``derived``; the radius,
S/N floor, offsets and the 20 % gate are ASSUMPTIONs from the plan.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

SAMPLE = "smacs0723_miri"
NIRCAM_BANDS = ("f090w", "f150w", "f200w", "f277w", "f356w", "f444w")
APERTURE = 1  # DJA ASEC_1 = 0.5" diameter, as D-013
MATCH_ARCSEC = 0.3  # ASSUMPTION (plan): about one F770W FWHM
MIN_SNR = 5.0  # ASSUMPTION: NIRCam detection floor, as quality.min_detection_snr
SHIFTS_ARCSEC = (5.0, 6.0, 7.0, 8.0, 9.0, 10.0, -5.0, -6.5, -8.0, -10.0)  # ASSUMPTION (plan)
MAX_CHANCE_FRACTION = 0.2  # ASSUMPTION (plan): chance rate / real rate gate


def nircam_detected(dja: Table, bands=NIRCAM_BANDS, aperture=APERTURE, min_snr=MIN_SNR):
    """Boolean mask: DJA rows with S/N >= ``min_snr`` in at least one of ``bands``."""
    ok = np.zeros(len(dja), dtype=bool)
    for b in bands:
        f = np.asarray(np.ma.filled(dja[f"{b}_flux_aper_{aperture}"], np.nan), float)
        e = np.asarray(np.ma.filled(dja[f"{b}_fluxerr_aper_{aperture}"], np.nan), float)
        with np.errstate(divide="ignore", invalid="ignore"):
            ok |= np.isfinite(f) & np.isfinite(e) & (e > 0) & (f / e >= min_snr)
    return ok


def match_fraction(ra, dec, ref: SkyCoord, radius_arcsec: float, dra_arcsec: float = 0.0):
    """Per-source match flag and separation (arcsec) to ``ref`` after an RA-direction shift."""
    if len(ra) == 0 or len(ref) == 0:
        return np.zeros(len(ra), dtype=bool), np.full(len(ra), np.nan)
    dec = np.asarray(dec, float)
    ra = np.asarray(ra, float) + dra_arcsec / 3600.0 / np.cos(np.radians(dec))
    _, sep, _ = SkyCoord(ra * u.deg, dec * u.deg).match_to_catalog_sky(ref)
    sep = sep.to_value(u.arcsec)
    return sep <= radius_arcsec, sep


def verdict(real_rate: float, chance_rate: float, gate: float = MAX_CHANCE_FRACTION) -> str:
    """Frozen plan rule: adopt NIRCam confirmation only if chance < gate x real."""
    if not np.isfinite(real_rate) or real_rate <= 0:
        return "not adopted: no real matches"
    return "adopt" if chance_rate < gate * real_rate else "not adopted: chance rate too high"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("run_dir", type=Path)
    p.add_argument("dja", type=Path)
    p.add_argument("--out", type=Path, default=Path("results/miri_confirm"))
    a = p.parse_args(argv)
    src = Table.read(a.run_dir / SAMPLE / "sources.ecsv")
    q = Table.read(a.run_dir / SAMPLE / "quality.ecsv")
    if list(src["source_uid"]) != list(q["source_uid"]):
        raise ValueError("sources and quality rows differ")
    reason = np.asarray(np.ma.filled(q["quality_reason"].astype(str), ""), str)
    dja = Table.read(a.dja)
    det = nircam_detected(dja)
    ref = SkyCoord(np.asarray(dja["ra"])[det] * u.deg, np.asarray(dja["dec"])[det] * u.deg)
    anchor = np.asarray(src["anchor_band"]).astype(str)
    groups = {
        "single_band_only": reason == "single_band",
        "single_band_any": np.char.find(reason, "single_band") >= 0,
        "passed_gate": reason == "",
    }
    out: dict = {
        "run": a.run_dir.name,
        "dja": a.dja.name,
        "n_sources": len(src),
        "n_dja_nircam_detected": int(det.sum()),
        "rule": {
            "match_arcsec": MATCH_ARCSEC,
            "min_snr": MIN_SNR,
            "shifts_arcsec": list(SHIFTS_ARCSEC),
            "max_chance_fraction": MAX_CHANCE_FRACTION,
        },
        "groups": {},
    }
    for name, m in groups.items():
        ra, dec = np.asarray(src["ra"])[m], np.asarray(src["dec"])[m]
        hit, sep = match_fraction(ra, dec, ref, MATCH_ARCSEC)
        chance = [match_fraction(ra, dec, ref, MATCH_ARCSEC, s)[0].mean() for s in SHIFTS_ARCSEC]
        real, ch = float(hit.mean()) if len(hit) else float("nan"), float(np.mean(chance))
        by_band = {
            b: {
                "n": int((anchor[m] == b).sum()),
                "matched": int(hit[anchor[m] == b].sum()),
                "rate": round(float(hit[anchor[m] == b].mean()), 3),
            }
            for b in sorted(set(anchor[m]))
        }
        out["groups"][name] = {
            "n": int(m.sum()),
            "matched": int(hit.sum()),
            "real_rate": real,
            "chance_rate": ch,
            "chance_rate_std": float(np.std(chance)),
            "median_sep_matched_arcsec": float(np.median(sep[hit])) if hit.any() else None,
            "by_anchor_band": by_band,
            "verdict": verdict(real, ch),
            "matched_uids": [str(x) for x in np.asarray(src["source_uid"])[m][hit]],
        }
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    for name, g in out["groups"].items():
        print(name, {k: g[k] for k in ("n", "matched", "real_rate", "chance_rate", "verdict")})
        print("   ", g["by_anchor_band"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
