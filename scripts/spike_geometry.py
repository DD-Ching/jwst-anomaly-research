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


def spike_flags(
    t: Table,
    rows: np.ndarray,
    bands_of: np.ndarray,
    rotate: float = 0.0,
    require_single: bool = True,
) -> list[dict]:
    """Evaluate the rule on ``rows``; ``bands_of[i]`` is the band whose shape is tested."""
    coords = SkyCoord(t["ra"], t["dec"], unit="deg")
    mags = {}
    out = []
    for i, b in zip(rows, bands_of, strict=True):
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
            rec.update(
                neighbour=str(t["source_uid"][best[0]]),
                dtheta=best[1],
                neighbour_sep=best[2],
                neighbour_mag=best[3],
            )
        rec["flag"] = bool(
            np.isfinite(ell)
            and ell >= E_MIN
            and best is not None
            and best[1] <= DTHETA
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
    return {str(u) for u, p in zip(qt["source_uid"], qt[col], strict=True) if bool(p)}


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
        "n_flagged_gate_pass": None
        if passed is None
        else sum(bool(r["gate_pass"]) for r in flagged),
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
            {
                "image": str(ids[k]),
                **{
                    c: r.get(c)
                    for c in (
                        "source_uid",
                        "ellipticity",
                        "dtheta",
                        "neighbour_sep",
                        "neighbour_mag",
                    )
                },
            }
            for k, r in enumerate(geom)
            if r["flag"]
        ],
    }


SEG_LEN = 64  # ASSUMPTION (addendum B, arm 2): row/column-segment length in block pixels
ENRICH = 3.0  # ASSUMPTION (addendum B): single-band in-mask fraction vs >= 3-band fraction


def streak_verdict(positive_in_mask: bool, n_images_in_mask: int, fields: list[dict]) -> dict:
    """Addendum B adoption rule: all three conditions must hold.

    ``fields``: per field ``{"single_in": k, "single_n": n, "multi_in": k3, "multi_n": n3}``
    (multi = detected in >= 3 bands). Condition 3 is judged only in fields with ``single_in`` >= 1
    and requires single_in/single_n >= ENRICH * multi_in/multi_n there (at least one such field).
    """
    judged = [f for f in fields if f["single_in"] >= 1]
    c3 = bool(judged) and all(
        f["single_in"] / f["single_n"] >= ENRICH * f["multi_in"] / max(f["multi_n"], 1)
        for f in judged
    )
    cond = {"positive": bool(positive_in_mask), "arcs": n_images_in_mask <= 1, "enrichment": c3}
    return {**cond, "adopt": all(cond.values())}


def streak_run(run_samples: dict[str, Path], positive_uid: str, out: Path, processes: int) -> dict:
    """Addendum B: W1 trail masks (arm 1) and trail + segment stripes (arm 2) per field."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import w1_streak_mask as wm

    from jwst_anomaly.pipeline import l3_image_uri

    cloud = {
        "s3_bucket": "stpubdata",
        "l3_key_pattern": "jwst/public/jw{proposal:05d}/L3/t/o{observation:03d}/{obs_id}_{suffix}",
    }
    obs = {
        "ceers_t021": "jw01345-o001_t021_nircam_clear-f200w",
        "smacs0723_nircam": "jw02736-o001_t001_nircam_clear-f200w",
    }
    arms: dict[str, dict] = {"trail": {}, "trail_seg64": {}}
    pos_in = {"trail": False, "trail_seg64": False}
    img_in = {"trail": 0, "trail_seg64": 0}
    for name, rs in run_samples.items():
        t = Table.read(rs / "sources.ecsv")
        uri = l3_image_uri(cloud, obs[name])
        img, wcs = wm.stream_binned(uri)
        z, wz = wm.normalise(img), wcs[:: wm.REBIN, :: wm.REBIN]
        trail, cat = wm.find_trails(z, processes)
        seg = trail | wm.to_mask_grid(wm.segment_stripes(img, SEG_LEN), trail.shape)
        nb = np.asarray(t["n_bands"])
        for arm, mask in (("trail", trail), ("trail_seg64", seg)):
            inm = wm.veto_rows(t["ra"], t["dec"], mask, wz)
            rec = {
                "i2d": uri,
                "mask_pixels": int(mask.sum()),
                "n_trails_status2": 0
                if cat is None
                else int((cat["status"] == wm.KEEP_STATUS).sum()),
                "single_in": int((inm & (nb == 1)).sum()),
                "single_n": int((nb == 1).sum()),
                "multi_in": int((inm & (nb >= 3)).sum()),
                "multi_n": int((nb >= 3).sum()),
            }
            if name == "ceers_t021":
                ip = find_uid(t, positive_uid)
                pos_in[arm] = bool(ip is not None and inm[ip])
            if name == "smacs0723_nircam":
                _, files, par = __import__("lens_consistency").load_model("smacs0723-iclv2")
                images = __import__("lens_consistency").image_list("smacs0723-iclv2", files, par)[0]
                img_in[arm] = int(wm.veto_rows(images["ra"], images["dec"], mask, wz).sum())
                rec["n_iclv2_images_in_mask"] = img_in[arm]
            arms[arm][name] = rec
            wm._sheet(out / f"{name}_{arm}_mask.png", img, mask, f"{name} {arm}")
    res = {
        "plan": "results/spike_geometry/PLAN.md addendum B",
        "veto_arcsec": wm.VETO_ARCSEC,
        "seg_len": SEG_LEN,
        "enrich": ENRICH,
        "arms": {},
    }
    for arm in arms:
        res["arms"][arm] = {
            "fields": arms[arm],
            "positive_in_mask": pos_in[arm],
            "verdict": streak_verdict(pos_in[arm], img_in[arm], list(arms[arm].values())),
        }
    (out / "summary_streak.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def main(argv=None) -> dict:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--ceers", type=Path, required=True, help="run dir of ceers_t021")
    p.add_argument("--smacs", type=Path, required=True, help="run dir of smacs0723_nircam")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--positive", default="_f200w_2842", help="known spike streak (uid suffix)")
    p.add_argument("--streak", action="store_true", help="run addendum B (streak masks) only")
    p.add_argument("--processes", type=int, default=4)
    a = p.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    if a.streak:
        return streak_run(
            {"ceers_t021": a.ceers, "smacs0723_nircam": a.smacs}, a.positive, a.out, a.processes
        )
    ceers = screen(a.ceers)
    smacs = screen(a.smacs)
    tc = Table.read(a.ceers / "sources.ecsv")
    ip = find_uid(tc, a.positive)
    pos = None
    if ip is not None:
        uid = str(tc["source_uid"][ip])
        hit = [r for r in ceers["flagged"] if r["source_uid"] == uid]
        rec = spike_flags(tc, np.array([ip]), detected_band(tc, bands_in(tc))[[ip]])[0]
        pos = {
            "source_uid": uid,
            "flagged": bool(hit),
            "record": rec,
            "ra": float(tc["ra"][ip]),
            "dec": float(tc["dec"][ip]),
        }
    summary = {
        "plan": "results/spike_geometry/PLAN.md",
        "thresholds": {
            "E_MIN": E_MIN,
            "DMAG": DMAG,
            "R_MAX": R_MAX,
            "DTHETA": DTHETA,
            "MATCH_ARCSEC": MATCH_ARCSEC,
        },
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
    if "arms" in s:
        print(json.dumps(s, indent=1))
        sys.exit(0)
    print(
        json.dumps(
            {k: v for k, v in s.items() if k not in ("ceers", "smacs")}, indent=1, default=float
        )[:3000]
    )
    for k in ("ceers", "smacs"):
        print(k, {x: s[k][x] for x in s[k] if x != "flagged"})
