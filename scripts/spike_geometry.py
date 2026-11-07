"""Spike-geometry check for single-band elongated detections (frozen plan:
``results/spike_geometry/PLAN.md``; TASKS Then 6, D-014 "Revisit if").

A source is ``spike_geom`` when it is detected in one band only (``n_bands == 1``), its pipeline
ellipticity in that band is >= ``E_MIN``, a source at least ``DMAG`` brighter in the same band lies
within ``R_MAX``, and the major axis (``<band>_sky_orientation``, deg E of N) points at that
neighbour within ``DTHETA``. The null rotates every major axis by 90 deg (positions and neighbours
kept), which estimates chance alignment.

Inputs: one pipeline run directory per sample (``sources.ecsv``; ``quality.ecsv`` when present for
the gate status). Known negatives: Mahler+2022 ICLv2 SMACS multiple images
(``lens_consistency.image_list``) matched to the SMACS merged catalog within ``MATCH_ARCSEC``.

Provenance: catalog columns are ``observed`` (MAST pipeline); flags and rates are ``derived``; all
thresholds are ASSUMPTIONs frozen in the plan.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

E_MIN = 0.8  # ASSUMPTION (plan): "highly elongated"
DMAG = 3.0  # ASSUMPTION (plan): neighbour at least this much brighter, same band
R_MAX = 10.0  # ASSUMPTION (plan): arcsec, search radius for the bright neighbour
DTHETA = 10.0  # ASSUMPTION (plan): deg, major axis vs bearing to the neighbour
MATCH_ARCSEC = 0.3  # ASSUMPTION: catalogued image to merged-catalog source


def axis_diff(a, b):
    """Smallest angle between two axes (deg, mod 180), in [0, 90]."""
    d = np.abs((np.asarray(a, float) - np.asarray(b, float)) % 180.0)
    return np.minimum(d, 180.0 - d)


def band_mag(t: Table, band: str) -> np.ndarray:
    """Total AB mag in ``band``, isophotal where the total is not finite; NaN if undetected."""
    tot = np.asarray(t[f"{band}_aper_total_abmag"], float)
    iso = np.asarray(t[f"{band}_isophotal_abmag"], float)
    m = np.where(np.isfinite(tot), tot, iso)
    return np.where(np.asarray(t[f"{band}_detected"], bool), m, np.nan)


def detected_band(t: Table, bands: list[str]) -> np.ndarray:
    """The single detected band of each row ('' when n_bands != 1)."""
    out = np.full(len(t), "", dtype=object)
    for b in bands:
        det = np.asarray(t[f"{b}_detected"], bool)
        out[det & (np.asarray(t["n_bands"]) == 1)] = b
    return out


def spike_flags(t: Table, rows: np.ndarray, bands_of: np.ndarray, rotate: float = 0.0,
                require_single: bool = True) -> list[dict]:
    """Evaluate the rule on ``rows``; ``bands_of[i]`` is the band whose shape is tested."""
    coords = SkyCoord(t["ra"], t["dec"], unit="deg")
    mags = {}
    out = []
    for i, b in zip(rows, bands_of):
        rec = {"row": int(i), "source_uid": str(t["source_uid"][i]), "band": b, "flag": False}
        if not b or (require_single and int(t["n_bands"][i]) != 1):
            out.append(rec)
            continue
        if b not in mags:
            mags[b] = band_mag(t, b)
        m = mags[b]
        ell = float(t[f"{b}_ellipticity"][i])
        pa = float(t[f"{b}_sky_orientation"][i]) + rotate
        rec.update(ellipticity=ell, mag=float(m[i]))
        sep = coords[i].separation(coords).arcsec
        cand = np.flatnonzero((sep > 0) & (sep <= R_MAX) & (m <= m[i] - DMAG))
        best = None
        for j in cand:
            bearing = coords[i].position_angle(coords[j]).deg
            d = float(axis_diff(pa, bearing))
            if best is None or d < best[1]:
                best = (int(j), d, float(sep[j]), float(m[j]))
        if best is not None:
            rec.update(neighbour=str(t["source_uid"][best[0]]), dtheta=best[1],
                       neighbour_sep=best[2], neighbour_mag=best[3])
        rec["flag"] = bool(
            np.isfinite(ell) and ell >= E_MIN and best is not None and best[1] <= DTHETA
            and np.isfinite(m[i])
        )
        out.append(rec)
    return out


def bands_in(t: Table) -> list[str]:
    return sorted(c[: -len("_detected")] for c in t.colnames if c.endswith("_detected"))


def gate_pass(run_sample: Path) -> set[str] | None:
    q = run_sample / "quality.ecsv"
    if not q.exists():
        return None
    qt = Table.read(q)
    col = "quality_ok" if "quality_ok" in qt.colnames else None
    if col is None:
        return None
    return {str(u) for u, p in zip(qt["source_uid"], qt[col]) if bool(p)}


def screen(run_sample: Path) -> dict:
    t = Table.read(run_sample / "sources.ecsv")
    bands = bands_in(t)
    single = np.flatnonzero(np.asarray(t["n_bands"]) == 1)
    bof = detected_band(t, bands)[single]
    real = spike_flags(t, single, bof)
    rot = spike_flags(t, single, bof, rotate=90.0)
    passed = gate_pass(run_sample)
    flagged = [r for r in real if r["flag"]]
    for r in flagged:
        r["gate_pass"] = None if passed is None else r["source_uid"] in passed
    n_rot = sum(r["flag"] for r in rot)
    return {
        "sample": run_sample.name,
        "run": run_sample.parent.name,
        "n_sources": len(t),
        "n_single_band": int(single.size),
        "n_flagged": len(flagged),
        "n_flagged_rotated": int(n_rot),
        "useful": bool(len(flagged) > 3 * n_rot),
        "n_flagged_gate_pass": None if passed is None else sum(bool(r["gate_pass"]) for r in flagged),
        "flagged": flagged,
    }


def find_uid(t: Table, label_suffix: str) -> int | None:
    hit = [i for i, u in enumerate(t["source_uid"]) if str(u).endswith(label_suffix)]
    return hit[0] if hit else None


def arcs_negative(smacs_sample: Path) -> dict:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import lens_consistency as lc

    _, files, par = lc.load_model("smacs0723-iclv2")
    images = lc.image_list("smacs0723-iclv2", files, par)[0]
    t = Table.read(smacs_sample / "sources.ecsv")
    bands = bands_in(t)
    cs = SkyCoord(t["ra"], t["dec"], unit="deg")
    ci = SkyCoord(images["ra"], images["dec"], unit="deg")
    idx, sep, _ = ci.match_to_catalog_sky(cs)
    ok = sep.arcsec <= MATCH_ARCSEC
    rows = idx[ok]
    ref = str(t["ref_band"][0]).lower() if "ref_band" in t.colnames else "f200w"
    # Shape band: the single detected band if n_bands == 1, else the reference band (or the anchor).
    sb = detected_band(t, bands)
    shape_band = []
    for i in rows:
        if sb[i]:
            shape_band.append(sb[i])
        elif bool(t[f"{ref}_detected"][i]):
            shape_band.append(ref)
        else:
            shape_band.append(str(t["anchor_band"][i]).lower())
    full = spike_flags(t, rows, np.array(shape_band, dtype=object))
    geom = spike_flags(t, rows, np.array(shape_band, dtype=object), require_single=False)
    ids = np.asarray(images["image_id"])[ok]
    return {
        "n_images": len(images),
        "n_matched": int(ok.sum()),
        "n_single_band": int(sum(int(t["n_bands"][i]) == 1 for i in rows)),
        "n_flagged_full_rule": int(sum(r["flag"] for r in full)),
        # Post hoc (not in the plan): the geometry conditions 2-4 alone, on every matched image.
        "n_flagged_geometry_only": int(sum(r["flag"] for r in geom)),
        "flagged_full_rule": [str(ids[k]) for k, r in enumerate(full) if r["flag"]],
        "flagged_geometry_only": [
            {"image": str(ids[k]), **{c: r.get(c) for c in ("source_uid", "ellipticity", "dtheta",
                                                         "neighbour_sep", "neighbour_mag")}}
            for k, r in enumerate(geom) if r["flag"]
        ],
    }


def main(argv=None) -> dict:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--ceers", type=Path, required=True, help="run dir of ceers_t021")
    p.add_argument("--smacs", type=Path, required=True, help="run dir of smacs0723_nircam")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--positive", default="_f200w_2842", help="known spike streak (uid suffix)")
    a = p.parse_args(argv)
    ceers = screen(a.ceers)
    smacs = screen(a.smacs)
    tc = Table.read(a.ceers / "sources.ecsv")
    ip = find_uid(tc, a.positive)
    pos = None
    if ip is not None:
        uid = str(tc["source_uid"][ip])
        hit = [r for r in ceers["flagged"] if r["source_uid"] == uid]
        rec = spike_flags(tc, np.array([ip]), detected_band(tc, bands_in(tc))[[ip]])[0]
        pos = {"source_uid": uid, "flagged": bool(hit), "record": rec,
               "ra": float(tc["ra"][ip]), "dec": float(tc["dec"][ip])}
    summary = {
        "plan": "results/spike_geometry/PLAN.md",
        "thresholds": {"E_MIN": E_MIN, "DMAG": DMAG, "R_MAX": R_MAX, "DTHETA": DTHETA,
                       "MATCH_ARCSEC": MATCH_ARCSEC},
        "positive": pos,
        "arcs_negative": arcs_negative(a.smacs),
        "ceers": ceers,
        "smacs": smacs,
    }
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    return summary


if __name__ == "__main__":
    s = main()
    print(json.dumps({k: v for k, v in s.items() if k not in ("ceers", "smacs")}, indent=1,
                     default=float)[:3000])
    for k in ("ceers", "smacs"):
        print(k, {x: s[k][x] for x in s[k] if x != "flagged"})
