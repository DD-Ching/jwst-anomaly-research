"""W1/W2 in published lens catalogues: lenses with no visible deflector (D-056).

Merges lenscat, Euclid Q1 and SuGOHI (``jwst_anomaly.lenscats``), asks Legacy Surveys DR10
(Tractor, NOIRLab Astro Data Lab TAP, batched box queries) whether a galaxy sits at each
galaxy-scale lens position, calibrates the lens magnitude expected for an Einstein radius on the
visible lenses (SIS + Faber-Jackson, ASSUMPTION), and lists the systems without a visible
deflector together with the cheapest ordinary flags. Every threshold is an ASSUMPTION
(``Params``); outputs are ``derived``. Nothing here is evidence of exotic physics.

    python scripts/w12_lenscats.py screen --out outputs/w12_lenscats
    python scripts/w12_lenscats.py vet --out outputs/w12_lenscats      # clusters, cutouts, sheet
"""

from __future__ import annotations

import argparse
import io
import json
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table, vstack

from jwst_anomaly import lenscats, paths, schema
from jwst_anomaly.photometry import fetch_catalog

TAP = "https://datalab.noirlab.edu/tap/sync"
TRACTOR_COLS = "ls_id,ra,dec,type,mag_g,mag_r,mag_z,shape_r,galdepth_z,nobs_z,maskbits,ref_cat"
CUTOUT = "https://www.legacysurvey.org/viewer/cutout.jpg"


@dataclass(frozen=True)
class Params:
    merge_radius: float = 3.0  # arcsec, friends-of-friends (SuGOHI merges at 3'')
    theta_e_max: float = 3.0  # arcsec, galaxy scale
    lens_radius: float = 1.5  # arcsec: an extended Tractor source this close is the deflector
    box: float = 5.0  # arcsec half-width of each Tractor box query (coverage and depth)
    batch: int = 300  # boxes per TAP query
    theta_e_default: float = 1.0  # arcsec when no catalogue gives theta_E
    theta_e_floor: float = 0.5  # conservative small Einstein radius for the limit
    z_l_typical: float = 0.5  # lens redshift of the "typical" variant when none is given
    z_s_default: float = 3.0  # conservative (faint-lens) source redshift when none is given
    calib_z_s: float = 2.0  # source redshift for calibration lenses without one
    faint_sigma: float = 2.0  # required lens = calibration + faint_sigma * rms
    margin: float = 0.5  # mag: detectable if required mag < 5-sigma depth - margin
    calib_sep: float = 0.5  # arcsec: calibration lenses have their galaxy this close
    image_radius: float = 3.0  # arcsec: point images of a quasar lens around the position
    sep_min: float = 2.0  # arcsec: image pairs closer than this hide the lens in LS seeing
    image_exclusion: float = 0.5  # arcsec: a deflector candidate is not an image
    n_rounded_sample: int = 32  # random rounded-position systems put on a contact sheet


def fetch_all() -> dict[str, Path]:
    return {k: fetch_catalog(v["url"], v["sha256"]) for k, v in lenscats.CATALOGUES.items()}


def load_systems(p: Params) -> tuple[Table, dict]:
    f = fetch_all()
    tabs = [
        lenscats.read_lenscat(f["lenscat"]),
        lenscats.read_euclid_q1(f["euclid_q1"], f["euclid_q1_mass"], f["euclid_q1_sersic_mag"]),
        lenscats.read_sugohi(f["sugohi"]),
    ]
    sizes = {t.meta["source"]: len(t) for t in tabs}
    return lenscats.merge(tabs, p.merge_radius), sizes


def box_terms(ra: float, dec: float, half_arcsec: float) -> list[str]:
    """ADQL conditions for one box; split at RA 0/360, full RA range near the poles."""
    h = half_arcsec / 3600
    dpart = f"dec BETWEEN {max(dec - h, -90):.7f} AND {min(dec + h, 90):.7f}"
    if abs(dec) + h >= 89.9:
        return [f"({dpart})"]
    hr = h / np.cos(np.radians(abs(dec) + h))
    lo, hi = ra - hr, ra + hr
    if lo < 0:
        ranges = [(0.0, hi), (lo + 360, 360.0)]
    elif hi >= 360:
        ranges = [(lo, 360.0), (0.0, hi - 360)]
    else:
        ranges = [(lo, hi)]
    return [f"(ra BETWEEN {a:.7f} AND {b:.7f} AND {dpart})" for a, b in ranges]


def tap_boxes(ra, dec, half_arcsec: float) -> Table:
    """One ADQL query for many boxes (the Data Lab TAP rejects uploads and q3c in ADQL)."""
    terms = [t for r, d in zip(ra, dec, strict=True) for t in box_terms(r, d, half_arcsec)]
    q = f"SELECT {TRACTOR_COLS} FROM ls_dr10.tractor WHERE " + " OR ".join(terms)
    err = ""
    for attempt in range(4):
        try:
            r = requests.post(
                TAP,
                data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q},
                timeout=900,
            )
            if r.status_code == 200 and r.text.startswith("ls_id"):
                return Table.read(io.BytesIO(r.content), format="ascii.csv", guess=False)
            err = r.text[:300]
        except requests.RequestException as e:
            err = str(e)
        time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"Data Lab TAP failed: {err}")


def _keys(ra, dec) -> np.ndarray:
    return np.array([f"{r:.6f},{d:.6f}" for r, d in zip(ra, dec, strict=True)])


def query_tractor(systems: Table, p: Params, cache: Path) -> Table:
    """Tractor sources in a box around every system, batched; positions already queried
    (``queried.ecsv``) are skipped and every cached batch is read back (a superset is harmless:
    classification is by position)."""
    import hashlib

    cache.mkdir(parents=True, exist_ok=True)
    qf = cache / "queried.ecsv"
    done = set(Table.read(qf)["key"]) if qf.exists() else set()
    keys = _keys(systems["ra"], systems["dec"])
    todo = np.flatnonzero(~np.isin(keys, list(done)))
    for i in range(0, len(todo), p.batch):
        s = systems[todo[i : i + p.batch]]
        bk = keys[todo[i : i + p.batch]]
        h = hashlib.sha1((";".join(bk) + f"|{p.box}").encode()).hexdigest()[:16]
        t = tap_boxes(np.asarray(s["ra"]), np.asarray(s["dec"]), p.box)
        t.write(cache / f"tractor_{h}.ecsv", overwrite=True)
        done.update(bk)
        Table({"key": sorted(done)}).write(qf, overwrite=True)
        print(f"tractor batch {i}/{len(todo)}: {len(t)} rows", flush=True)
    parts = []
    for f in sorted(cache.glob("tractor_*.ecsv")):
        t = Table.read(f)
        if len(t):
            t["type"] = t["type"].astype(str)
            t["ref_cat"] = t["ref_cat"].astype(str)
            parts.append(t)
    out = vstack(parts, metadata_conflicts="silent") if parts else Table()
    out = Table(out, masked=False) if len(out) else out
    for c in ("mag_g", "mag_r", "mag_z", "galdepth_z", "shape_r"):
        if c in out.colnames:
            out[c] = np.asarray(
                out[c].filled(np.nan) if hasattr(out[c], "filled") else out[c], float
            )
    ids, first = np.unique(np.asarray(out["ls_id"]), return_index=True)
    out = out[np.sort(first)]
    out.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source="Legacy Surveys DR10 Tractor (Data Lab TAP)",
    )
    return out


def calibrate(systems: Table, defl: Table, p: Params) -> lenscats.FJCalibration:
    """FJ calibration (LS z) on galaxy-selected lenses whose catalogued lens position has an
    extended Tractor source within ``calib_sep`` (the lens galaxy, not an arc or a quasar image)
    and a theta_E and z_lens."""
    vis = np.asarray(defl["visible_deflector"], bool) & (
        np.asarray(defl["sep_ext"], float) <= p.calib_sep
    )
    mag = np.asarray(defl["mag_ext"], float)
    th = np.asarray(systems["theta_e"], float)
    zl = np.asarray(systems["z_lens"], float)
    zs = np.asarray(systems["z_source"], float)
    zs = np.where(np.isfinite(zs), zs, p.calib_z_s)
    use = vis & np.isfinite(mag) & np.isfinite(th) & (th > 0) & (th <= p.theta_e_max)
    use &= np.isfinite(zl) & (zl > 0.05) & (zl < zs - 0.1)
    use &= np.asarray(systems["selection"]) == "galaxy"
    return lenscats.fit_fj(th[use], zl[use], zs[use], mag[use], band="ls_z")


def required_mags(t: Table, theta_e, cal: lenscats.FJCalibration, p: Params) -> dict:
    """Typical, conservative (max over z_l, z_s default 3, + faint_sigma rms) and floor
    (conservative with theta_E floor where theta_E is unknown) required lens magnitudes."""
    th = np.asarray(theta_e, float)
    zs_c = np.where(np.isfinite(np.asarray(t["z_source"], float)), t["z_source"], p.z_s_default)
    th_used = np.where(np.isfinite(th), th, p.theta_e_default)
    cons, zcons = lenscats.required_lens_mag(th_used, zs_c, cal, n_sigma_faint=p.faint_sigma)
    floor, _ = lenscats.required_lens_mag(
        np.where(np.isfinite(th), th, p.theta_e_floor), zs_c, cal, n_sigma_faint=p.faint_sigma
    )
    zl = np.asarray(t["z_lens"], float)
    zl_t = np.where(np.isfinite(zl), zl, p.z_l_typical)
    zs_t = np.where(np.isfinite(np.asarray(t["z_source"], float)), t["z_source"], p.calib_z_s)
    zl_t = np.where(zl_t < zs_t - 0.1, zl_t, np.nan)
    typ = cal.mag(lenscats.sis_sigma(th_used, zl_t, zs_t), zl_t)
    # faintest magnitude a deflector may have and still count as the lens: typical + 2 rms,
    # else the conservative value, else -inf (nothing qualifies; never +inf)
    lim = typ + p.faint_sigma * cal.rms
    lim = np.where(np.isfinite(lim), lim, cons)
    lim = np.where(np.isfinite(lim), lim, -np.inf)
    return {
        "theta_e_used": th_used,
        "req_mag_z_typical": typ,
        "req_mag_z": cons,
        "req_z_l": zcons,
        "req_mag_z_floor": floor,
        "defl_mag_max": lim,
    }


BRICKS_QUERY = (
    "SELECT brickname, ra1, ra2, dec1, dec2, nexp_r, nexp_z, galdepth_z FROM ls_dr10.bricks_s "
    "WHERE survey_primary = 1 AND nexp_z > 0"
)


def load_bricks(out: Path) -> tuple[Table, dict]:
    """DR10 (DECam) brick summary: footprint and depth independent of detected sources
    (``ls_dr10.tractor`` is DECam-only, so the DR9 north bricks are not used). Cached once."""
    import hashlib

    f = out / "bricks_dr10_south.csv"
    if not f.exists():
        r = requests.post(
            TAP,
            data={
                "REQUEST": "doQuery",
                "LANG": "ADQL",
                "FORMAT": "csv",
                "MAXREC": "1000000",
                "QUERY": BRICKS_QUERY,
            },
            timeout=900,
        )
        if r.status_code != 200 or not r.text.startswith("brickname"):
            raise RuntimeError(f"brick query failed: {r.text[:300]}")
        f.write_bytes(r.content)
    data = f.read_bytes()
    meta = {
        "query": BRICKS_QUERY,
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
    }
    t = Table.read(f, format="ascii.csv", guess=False)
    t.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source="Legacy Surveys DR10 brick summary ls_dr10.bricks_s (Data Lab TAP)",
    )
    meta["rows"] = len(t)
    return t, meta


def load_sources(out: Path) -> Table:
    parts = [Table.read(f) for f in sorted((out / "tractor_cache").glob("tractor_*.ecsv"))]
    tr = vstack([t for t in parts if len(t)], metadata_conflicts="silent")
    tr = Table(tr, masked=False)
    tr["type"] = np.array([str(x) for x in tr["type"]])
    for c in ("mag_z", "galdepth_z"):
        tr[c] = np.asarray(np.ma.filled(np.ma.asarray(tr[c], float), np.nan))
    tr["maskbits"] = np.asarray(np.ma.filled(np.ma.asarray(tr["maskbits"], int), 0))
    _, first = np.unique(np.asarray(tr["ls_id"]), return_index=True)
    return tr[np.sort(first)]


def deflector_test(t: Table, src: Table, p: Params, images: Table) -> Table:
    """The per-class test a deflector without light could fail (D-056 reviews).

    - galaxy-selected: "insensitive" (the finder needed a galaxy at the position);
    - submm: "insensitive" (single-dish centroids are uncertain by more than theta_E);
    - quasar: :func:`lenscats.quasar_pair_test` on ``images``;
    - radio (interferometric positions, ASSUMPTION error < theta_E): an extended source within
      ``lens_radius`` bright enough for the lens -> "deflector"; a fainter one -> "faint
      galaxy"; else, with optical point sources, the quasar test; else "none".
    "faint galaxy" is undecided: a galaxy sits where the lens should be, and whether it is
    luminous enough depends on the Faber-Jackson scatter beyond the 2-rms margin, so the system
    can neither show nor exclude a dark deflector. Returns ``test_status``, ``n_bright_gal``.
    """
    sel = np.asarray(t["selection"])
    lim = np.asarray(t["defl_mag_max"], float)
    pair = lenscats.quasar_pair_test(
        t, src, lim, images, p.image_radius, p.sep_min, p.image_exclusion
    )
    ngal, _ = lenscats.bright_galaxy_near(t, src, lim, p.lens_radius)
    nany, _ = lenscats.bright_galaxy_near(t, src, np.full(len(t), np.inf), p.lens_radius)
    status = np.array(["insensitive"] * len(t), dtype=object)
    q = sel == "quasar"
    status[q] = np.asarray(pair["status"], dtype=object)[q]
    rad = np.where(
        ngal > 0,
        "deflector",
        np.where(
            nany > 0,
            "faint galaxy",
            np.where(np.asarray(images["n_images"]) > 0, np.asarray(pair["status"]), "none"),
        ),
    )
    r = sel == "radio"
    status[r] = rad[r]
    return Table({"test_status": status, "n_bright_gal": np.asarray(ngal, int)})


def flags_undecided(t: Table, src: Table, p: Params) -> np.ndarray:
    """Position and imaging problems that make the deflector test uninformative, evaluated for
    every covered system whatever the test said."""
    near = lenscats._neighbours(t, src, p.image_radius)
    masked = np.array(
        [bool(np.any((np.asarray(src["maskbits"])[idx] & MASK_BAD) > 0)) for idx, _ in near]
    )
    reasons = []
    for row, m in zip(t, masked, strict=True):
        r = []
        if row["name_offset"] > 5:
            r.append("position error (name vs RA/Dec)")
        if row["pos_quantum"] > p.lens_radius:
            r.append("position rounded")
        if row["catalogues"] == "euclid_q1":
            r.append("Euclid host position")
        if m:
            r.append("LS maskbits")
        reasons.append("; ".join(r))
    return np.array(reasons, dtype=object)


def cmd_screen(args, p: Params) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    systems, sizes = load_systems(p)
    counts = {"entries": sizes, "systems": len(systems)}
    galaxy = ~np.asarray(systems["any_group_scale"], bool) & ~(
        np.asarray(systems["theta_e"], float) > p.theta_e_max
    )
    counts["galaxy_scale"] = int(galaxy.sum())
    counts["galaxy_scale_no_lens_z_no_lens_mag"] = int(
        lenscats.select_no_lens_info(systems, p.theta_e_max).sum()
    )
    gs = systems[galaxy]
    gs["selection"] = [selection_class(r) for r in gs]
    bricks, bmeta = load_bricks(out)
    counts["bricks"] = bmeta
    query_tractor(gs, p, out / "tractor_cache")
    src = load_sources(out)
    defl = lenscats.classify_deflectors(gs, src, p.lens_radius, p.box)
    foot = lenscats.brick_coverage(gs, bricks)
    defl["covered"] = foot["covered"]  # footprint, not "a source was detected nearby"
    defl["depth_z"] = foot["depth_z"]
    cov = np.asarray(defl["covered"], bool)
    cal = calibrate(gs, defl, p)

    res = gs.copy()
    for c in defl.colnames[1:]:
        res[c] = defl[c]
    res["name_offset"] = [
        lenscats.name_position_offset(n, r, d)
        for n, r, d in zip(res["name"], res["ra"], res["dec"], strict=True)
    ]
    res["pos_quantum"] = [
        lenscats.position_quantum_arcsec(r, d) for r, d in zip(res["ra"], res["dec"], strict=True)
    ]
    images = lenscats.pair_images(res, src, p.image_radius)
    for c in images.colnames:
        res[f"pair_{c}"] = images[c]
    # theta_E: the catalogue's; for quasar and radio systems without one, half the image
    # separation; galaxy and sub-mm systems keep the catalogue value (or the default)
    th = np.asarray(res["theta_e"], float)
    qr = np.isin(res["selection"], ("quasar", "radio")) & ~np.isfinite(th)
    th = np.where(qr, np.asarray(images["theta_e"], float), th)
    for k, v in required_mags(res, th, cal, p).items():
        res[k] = v
    test = deflector_test(res, src, p, images)
    for c in test.colnames:
        res[c] = test[c]
    res["visible_bright"] = np.asarray(res["n_bright_gal"]) > 0
    depth = np.asarray(res["depth_z"], float)
    res["detectable_typical"] = res["req_mag_z_typical"] < depth - p.margin
    res["detectable_floor"] = res["req_mag_z_floor"] < depth - p.margin
    res["undecided"] = flags_undecided(res, src, p)
    res.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="; ".join(
            (
                str(res.meta.get("source", "merged lens catalogues")),
                str(src.meta.get("source", "LS DR10 Tractor")),
                "Legacy Surveys DR10 brick summary",
            )
        ),
        params=asdict(p),
        calibration=asdict(cal),
    )
    res.write(out / "systems.ecsv", overwrite=True)

    sel = np.asarray(res["selection"])
    vis_any = np.asarray(res["visible_deflector"], bool)
    vis = np.asarray(res["visible_bright"], bool)
    counts.update(
        ls_covered=int(cov.sum()),
        ls_extended_within_lens_radius=int((cov & vis_any).sum()),
        ls_bright_enough_within_lens_radius=int((cov & vis).sum()),
        selection_galaxy_scale={c: int((sel == c).sum()) for c in SELECTIONS},
        selection_covered={c: int((cov & (sel == c)).sum()) for c in SELECTIONS},
        median_depth_z=float(np.nanmedian(depth[cov])),
        depth_z_p05_p50_p95=[float(x) for x in np.nanpercentile(depth[cov], [5, 50, 95])],
        median_req_mag_z_typical=float(np.nanmedian(np.asarray(res["req_mag_z_typical"])[cov])),
        median_req_mag_z_floor=float(np.nanmedian(np.asarray(res["req_mag_z_floor"])[cov])),
        covered_detectable_fraction_typical=float(np.mean(res["detectable_typical"][cov])),
        covered_detectable_fraction_floor=float(np.mean(res["detectable_floor"][cov])),
        undecided_covered=int((cov & (np.asarray(res["undecided"]) != "")).sum()),
        test_status_by_selection={
            c: dict(
                zip(
                    *[
                        x.tolist()
                        for x in np.unique(
                            np.asarray(res["test_status"])[cov & (sel == c)], return_counts=True
                        )
                    ],
                    strict=True,
                )
            )
            for c in SELECTIONS
        },
    )
    summary = {
        "counts": counts,
        "calibration": asdict(cal),
        "params": asdict(p),
        "wall_s": round(time.time() - t0, 1),
        "provenance": schema.Provenance.DERIVED.value,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["counts"], indent=1))


SELECTIONS = ("galaxy", "quasar", "radio", "submm")

# Cluster and object catalogues for vetting (CDS XMatch, one request per catalogue).
VET_XMATCH = {
    "redmapper": ("vizier:J/ApJS/224/1/cat_dr8", 120.0),  # Rykoff+2016 SDSS DR8 clusters
    "whl": ("vizier:J/ApJS/199/34/table1", 120.0),  # Wen, Han & Liu 2012 SDSS-III clusters
    "simbad": ("simbad", 5.0),
}
# Legacy Surveys maskbits that degrade detection at a position (BRIGHT, MEDIUM, GALAXY, CLUSTER).
MASK_BAD = (1 << 1) | (1 << 11) | (1 << 12) | (1 << 13)


_SIMBAD_GALAXY_TYPES = {"Galaxy", "LensingG", "BrightestCG", "RadioG", "EmissionG", "PairG"}
_SIMBAD_GALAXY_TYPES |= {"GtowardsCl", "GtowardsGroup", "GinCl", "GinGroup", "GinPair", "LSB_G"}
_SIMBAD_GALAXY_TYPES |= {"StarburstG", "HIIG", "Compact_Gr_G", "EllipticalG", "SpiralG"}


def is_simbad_galaxy(main_id: str, main_type: str) -> bool:
    """A SIMBAD galaxy, or a lens-galaxy component named "... G" / "G1" (CASTLES-style)."""
    return main_type in _SIMBAD_GALAXY_TYPES or bool(re.search(r"\sG\d?$", main_id))


SIMBAD_CLUSTER_TYPES = {"ClG", "GroupG", "CGG", "PCG", "protoClG", "BrightestCG", "SuperClG"}
# Single-dish sub-mm lens surveys: the catalogued position is a centroid uncertain by arcsec.
_SUBMM = re.compile(
    r"\bSPT|ACT-S|H-?ATLAS|HeLMS|HELMS|\bHERS\d|HerBS|NKC2016|ACS2016|PLCK|Negrello|"
    r"Amvrosiadis|Nayyeri"
)
# Radio-interferometric lens surveys (VLA / VLBI positions, ASSUMPTION: error < theta_E).
_RADIO = re.compile(r"\bCLASS\b|JVAS|\bMG\d|\bB\d{4}[+-]\d|SML2019|mJIVE|MJV\d")


# Systems the catalogue tests leave open, settled by reading the discovery paper (D-056).
LITERATURE = {
    "J1329+4325": "[SML2019] MJV16999: mJIVE-20 VLBI milli-lens candidate rejected as a core-jet "
    "source by Spingola et al. 2019 (arXiv:1811.09152, sect. 4.1.12)",
    "221216-010345": "HSC J2212-0103: likely lensed quasar; He et al. 2025 (arXiv:2509.03858, "
    "lens-light table) fit a lens galaxy in HSC with i = 22.40, r = 23.30",
}


def selection_class(r) -> str:
    """How a system was found: "quasar" (lensed-quasar searches), "radio" (interferometric
    radio searches), "submm" (single-dish sub-mm flux selection) or "galaxy" (finders run
    around galaxies, or unknown). Only quasar and radio searches can select a lens without a
    visible deflector *and* locate it to better than theta_E."""
    text = " ".join(str(r[c]) for c in ("name", "entries", "refs", "system_type"))
    if _SUBMM.search(text):
        return "submm"
    if _RADIO.search(text):
        return "radio"
    if "lensedquasars" in text or "GQ" in str(r["system_type"]):
        return "quasar"
    return "galaxy"


def xmatch(t: Table, cat: str, radius: float) -> Table:
    from astropy import units as u
    from astroquery.xmatch import XMatchClass

    up = Table({"sid": np.asarray(t["system_id"], str), "ra0": t["ra"], "dec0": t["dec"]})
    xm = XMatchClass()
    xm.TIMEOUT = 600
    return xm.query(
        cat1=up, cat2=cat, max_distance=radius * u.arcsec, colRA1="ra0", colDec1="dec0", cache=False
    )


def add_xmatch(s: Table, p: Params) -> None:
    """SIMBAD / cluster columns; a failed service leaves "not run" (never crashes)."""
    for c in ("redmapper", "whl", "simbad", "simbad_galaxy_3as", "simbad_cluster"):
        s[c] = np.full(len(s), "not run" if len(s) else "", dtype=object)
    if not len(s):
        return
    for key, (cat, rad) in VET_XMATCH.items():
        try:
            m = xmatch(s, cat, rad)
        except Exception as e:  # noqa: BLE001 - record and continue
            print(f"xmatch {key} failed: {e}", flush=True)
            continue
        m.sort("angDist")
        hit: dict[str, str] = {}
        gal: dict[str, str] = {}
        clu: dict[str, str] = {}
        for row in m:
            sid = str(row["sid"])
            if key == "simbad":
                label = f"{row['main_id']} [{row['main_type']}]"
                if row["angDist"] <= p.theta_e_max and is_simbad_galaxy(
                    str(row["main_id"]), str(row["main_type"])
                ):
                    gal.setdefault(sid, label)
                if str(row["main_type"]) in SIMBAD_CLUSTER_TYPES:
                    clu.setdefault(sid, label)
            else:
                skip = {"angDist", "sid", "ra0", "dec0", "_RAJ2000", "_DEJ2000", "recno"}
                idcol = next(c for c in m.colnames if c not in skip)
                label = f"{row[idcol]} @{row['angDist']:.0f}''"
            hit.setdefault(sid, label)
        s[key] = [hit.get(str(x), "") for x in s["system_id"]]
        if key == "simbad":
            s["simbad_galaxy_3as"] = [gal.get(str(x), "") for x in s["system_id"]]
            s["simbad_cluster"] = [clu.get(str(x), "") for x in s["system_id"]]


def fetch_cutout(ra: float, dec: float, path: Path, size: int = 80) -> bool:
    if path.exists():
        return True
    try:
        r = requests.get(
            CUTOUT,
            params={"ra": ra, "dec": dec, "layer": "ls-dr10", "pixscale": 0.262, "size": size},
            timeout=60,
        )
    except requests.RequestException:
        return False
    if r.status_code != 200 or not r.content[:3] == b"\xff\xd8\xff":
        return False
    path.write_bytes(r.content)
    time.sleep(0.5)  # politeness towards the viewer
    return True


def contact_sheet(t: Table, cut_dir: Path, out_png: Path, title: str, ncol: int = 8) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    n = len(t)
    nrow = max(1, int(np.ceil(n / ncol)))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.0 * ncol, 2.25 * nrow), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, row in zip(axes.ravel(), t, strict=False):
        f = cut_dir / f"{row['system_id']}.jpg"
        if f.exists():
            img = mpimg.imread(f)
            ax.imshow(img, origin="upper")
            c = img.shape[0] / 2
            ax.add_patch(plt.Circle((c, c), 1.5 / 0.262, fill=False, color="w", lw=0.6, ls=":"))
        ax.set_title(f"{row['system_id']} {str(row['name'])[:16]}\n{row['label'][:28]}", fontsize=6)
    fig.suptitle(f"{title}: LS DR10 grz, 21'' (dotted 1.5'')", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def _count(values) -> dict:
    vals, cnt = np.unique(np.asarray(list(values), dtype=str), return_counts=True)
    return {str(v): int(c) for v, c in zip(vals, cnt, strict=True)}


def poisson95(k: int) -> float:
    from scipy.stats import chi2

    return float(0.5 * chi2.ppf(0.95, 2 * k + 2))


def cmd_vet(args, p: Params) -> None:
    """Vet every sensitive system whose test found no deflector; limits per selection class.

    The limits assume the pair / radio test is complete for the decided systems: a dark lens
    in them would give "none". Not simulated: blending of a lens with the images, Tractor
    typing of a compact lens as PSF, seeing and depth variations at the position (a realistic
    injection needs re-running Tractor on images, which this offline step cannot do)."""
    out = Path(args.out)
    res = Table.read(out / "systems.ecsv")
    cov = np.asarray(res["covered"], bool)
    sens = np.isin(res["selection"], ("quasar", "radio"))
    decided = cov & sens & (np.asarray(res["undecided"]) == "")
    decided &= np.isin(res["test_status"], ("deflector", "none"))  # "faint galaxy" undecided
    cand = res[decided & (res["test_status"] == "none")]
    add_xmatch(cand, p)

    def ran(v) -> bool:
        return str(v) not in ("", "not run")

    tests = (  # ordinary explanations of a deflector-less configuration, cheapest first
        (
            "SIMBAD galaxy within 3'' (literature lens galaxy)",
            lambda r: ran(r["simbad_galaxy_3as"]),
        ),
        ("lens redshift published (deflector observed)", lambda r: bool(r["lens_z_known"])),
        (
            "cluster: redMaPPer/WHL within 2' or SIMBAD cluster within 5''",
            lambda r: ran(r["redmapper"]) or ran(r["whl"]) or ran(r["simbad_cluster"]),
        ),
        ("literature: not a lens / explained", lambda r: str(r["name"]) in LITERATURE),
    )
    flags = [[n for n, f in tests if f(r)] for r in cand]
    cand["flags"] = ["; ".join(f) for f in flags]
    cand["label"] = [f[0] if f else "unexplained by catalogue tests" for f in flags]
    cand.meta.update(provenance=schema.Provenance.DERIVED.value)
    cand.write(out / "candidates_vetted.ecsv", overwrite=True)

    # cutouts: every candidate, a random sample of rounded positions, and sensitive controls
    cut = out / "cutouts"
    cut.mkdir(exist_ok=True)
    rng = np.random.default_rng(56)
    rounded = np.flatnonzero(cov & (np.asarray(res["pos_quantum"]) > p.lens_radius))
    rsample = res[np.sort(rng.choice(rounded, min(len(rounded), p.n_rounded_sample), False))]
    rsample["label"] = [f"rounded {q:.1f}''" for q in rsample["pos_quantum"]]
    ctrl_idx = np.flatnonzero(decided & (res["test_status"] == "deflector"))
    ctrl = res[np.sort(rng.choice(ctrl_idx, min(len(ctrl_idx), args.n_control), False))]
    ctrl["label"] = ["deflector found"] * len(ctrl)
    sheets = {"candidates": cand, "rounded_sample": rsample, "controls": ctrl}
    for name, t in sheets.items():
        for r in t:
            fetch_cutout(float(r["ra"]), float(r["dec"]), cut / f"{r['system_id']}.jpg")
        for i in range(0, len(t), 64):
            contact_sheet(t[i : i + 64], cut, out / f"sheet_{name}_{i // 64:02d}.png", name)

    # a candidate counts against a variant and class if no test explains it and an ordinary lens
    # would have been detectable in that variant
    open_ = np.asarray(cand["label"]) == "unexplained by catalogue tests"
    csel = np.asarray(cand["selection"])
    summary = json.loads((out / "summary.json").read_text())
    summary["vetting"] = {
        "sensitive_covered": int((cov & sens).sum()),
        "sensitive_decided": int(decided.sum()),
        "no_deflector": len(cand),
        "labels": dict(
            zip(
                *[x.tolist() for x in np.unique(np.asarray(cand["label"]), return_counts=True)],
                strict=True,
            )
        )
        if len(cand)
        else {},
        "sensitive_undecided_by_reason": _count(
            x for u in np.asarray(res["undecided"])[cov & sens] for x in str(u).split("; ") if x
        ),
        "sensitive_status": _count(np.asarray(res["test_status"])[cov & sens]),
        "rounded_covered": int(len(rounded)),
        "rounded_sample_inspected": int(len(rsample)),
        "controls_inspected": int(len(ctrl)),
    }
    lim = {"open_candidates": int(open_.sum()), "efficiency_assumed": 1.0}
    for variant, col in (("typical", "detectable_typical"), ("conservative", "detectable_floor")):
        for c in ("quasar", "radio", "all"):
            m = decided & np.asarray(res[col], bool)
            mc = open_ & np.asarray(cand[col], bool)
            if c != "all":
                m &= np.asarray(res["selection"]) == c
                mc &= csel == c
            n, k = int(m.sum()), int(mc.sum())
            lim[f"{variant}_{c}_N"] = n
            lim[f"{variant}_{c}_k"] = k
            lim[f"{variant}_{c}_f95"] = poisson95(k) / n if n else None
    summary["limits"] = lim
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: summary[k] for k in ("vetting", "limits")}, indent=1))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["screen", "vet"])
    ap.add_argument("--n-control", type=int, default=32, help="random deflector-found cutouts")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    args.out = args.out or paths.outputs_dir() / "w12_lenscats"
    p = Params()
    {"screen": cmd_screen, "vet": cmd_vet}[args.cmd](args, p)


if __name__ == "__main__":
    main()
