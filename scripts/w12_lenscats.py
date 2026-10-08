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
    offset_radius: float = 5.0  # arcsec: catalogue position error plus theta_E (vetting)


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


def tap_boxes(ra, dec, half_arcsec: float) -> Table:
    """One ADQL query for many boxes (the Data Lab TAP rejects uploads and q3c in ADQL)."""
    h = half_arcsec / 3600
    terms = []
    for r, d in zip(ra, dec, strict=True):
        hr = h / max(np.cos(np.radians(d)), 1e-3)
        terms.append(
            f"(ra BETWEEN {r - hr:.7f} AND {r + hr:.7f}"
            f" AND dec BETWEEN {d - h:.7f} AND {d + h:.7f})"
        )
    q = f"SELECT {TRACTOR_COLS} FROM ls_dr10.tractor WHERE " + " OR ".join(terms)
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
    """FJ calibration on visible lenses with theta_E and z_lens (LS z-band magnitude)."""
    j = {s: i for i, s in enumerate(defl["system_id"])}
    idx = [j[s] for s in systems["system_id"]]
    vis = np.asarray(defl["visible_deflector"], bool)[idx]
    mag = np.asarray(defl["mag_ext"], float)[idx]
    th = np.asarray(systems["theta_e"], float)
    zl = np.asarray(systems["z_lens"], float)
    zs = np.asarray(systems["z_source"], float)
    zs = np.where(np.isfinite(zs), zs, p.calib_z_s)
    use = vis & np.isfinite(mag) & np.isfinite(th) & (th > 0) & (th <= p.theta_e_max)
    use &= np.isfinite(zl) & (zl > 0.05) & (zl < zs - 0.1)
    return lenscats.fit_fj(th[use], zl[use], zs[use], mag[use], band="ls_z")


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
    no_info = lenscats.select_no_lens_info(systems, p.theta_e_max)
    counts["galaxy_scale_no_lens_z_no_lens_mag"] = int(no_info.sum())
    gs = systems[galaxy]
    src = query_tractor(gs, p, out / "tractor_cache")
    defl = lenscats.classify_deflectors(gs, src, p.lens_radius, p.box)
    cov = np.asarray(defl["covered"], bool)
    vis = np.asarray(defl["visible_deflector"], bool)
    counts["ls_covered"] = int(cov.sum())
    counts["ls_visible_deflector"] = int((cov & vis).sum())
    counts["ls_no_extended_within_lens_radius"] = int((cov & ~vis).sum())
    ni = no_info[galaxy]
    counts["no_lens_info_covered"] = int((ni & cov).sum())
    counts["no_lens_info_no_extended"] = int((ni & cov & ~vis).sum())

    cal = calibrate(gs, defl, p)
    th = np.asarray(gs["theta_e"], float)
    zs = np.asarray(gs["z_source"], float)
    zs = np.where(np.isfinite(zs), zs, p.z_s_default)
    th_used = np.where(np.isfinite(th), th, p.theta_e_default)
    req, zreq = lenscats.required_lens_mag(th_used, zs, cal, n_sigma_faint=p.faint_sigma)
    req_floor, _ = lenscats.required_lens_mag(
        np.where(np.isfinite(th), th, p.theta_e_floor), zs, cal, n_sigma_faint=p.faint_sigma
    )
    # typical ordinary lens: catalogue z_l (else z_l_typical), z_s (else calib_z_s), no extra
    # faintness; the conservative columns above maximise over z_l and add faint_sigma * rms
    zl = np.asarray(gs["z_lens"], float)
    zl_t = np.where(np.isfinite(zl), zl, p.z_l_typical)
    zs_t = np.where(np.isfinite(np.asarray(gs["z_source"], float)), gs["z_source"], p.calib_z_s)
    zl_t = np.where(zl_t < zs_t - 0.1, zl_t, np.nan)
    req_typ = cal.mag(lenscats.sis_sigma(th_used, zl_t, zs_t), zl_t)
    depth = np.asarray(defl["depth_z"], float)
    res = gs.copy()
    for c in defl.colnames[1:]:
        res[c] = defl[c]
    res["theta_e_used"] = th_used
    res["z_s_used"] = zs
    res["req_mag_z"] = req
    res["req_z_l"] = zreq
    res["req_mag_z_floor"] = req_floor
    res["detectable"] = req < depth - p.margin
    res["detectable_floor"] = req_floor < depth - p.margin
    res["req_mag_z_typical"] = req_typ
    res["detectable_typical"] = req_typ < depth - p.margin
    res["name_offset"] = [
        lenscats.name_position_offset(n, r, d)
        for n, r, d in zip(res["name"], res["ra"], res["dec"], strict=True)
    ]
    res["pos_quantum"] = [
        lenscats.position_quantum_arcsec(r, d) for r, d in zip(res["ra"], res["dec"], strict=True)
    ]
    res["selection"] = [selection_class(r) for r in res]
    res.meta.update(
        provenance=schema.Provenance.DERIVED.value, params=asdict(p), calibration=asdict(cal)
    )
    res.write(out / "systems.ecsv", overwrite=True)
    surv = res[cov & ~vis]
    surv.write(out / "no_visible_deflector.ecsv", overwrite=True)
    counts["no_extended_and_detectable"] = int(np.sum(surv["detectable"]))
    counts["no_extended_and_detectable_floor"] = int(np.sum(surv["detectable_floor"]))
    counts["covered_detectable_fraction"] = float(np.mean(res["detectable"][cov]))
    counts["covered_detectable_fraction_floor"] = float(np.mean(res["detectable_floor"][cov]))
    counts["no_extended_and_detectable_typical"] = int(np.sum(surv["detectable_typical"]))
    counts["covered_detectable_fraction_typical"] = float(np.mean(res["detectable_typical"][cov]))
    counts["median_depth_z"] = float(np.nanmedian(depth[cov]))
    summary = {
        "counts": counts,
        "calibration": asdict(cal),
        "params": asdict(p),
        "wall_s": round(time.time() - t0, 1),
        "provenance": schema.Provenance.DERIVED.value,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


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
# Sub-mm / radio lens surveys: the catalogued position is the source centroid, not the lens.
_SUBMM_RADIO = re.compile(
    r"\bSPT|ACT-S|H-?ATLAS|HeLMS|HELMS|\bHERS\d|HerBS|NKC2016|ACS2016|PLCK|\bCLASS\b|JVAS|"
    r"\bMG\d|Negrello|Amvrosiadis|Nayyeri"
)


# Systems the catalogue tests leave open, settled by reading the discovery paper (D-056).
LITERATURE = {
    "J1329+4325": "[SML2019] MJV16999: mJIVE-20 VLBI milli-lens candidate rejected as a core-jet "
    "source by Spingola et al. 2019 (arXiv:1811.09152, sect. 4.1.12)",
}


def selection_class(r) -> str:
    """How a system was found: "quasar" (lensed-quasar searches), "submm/radio" (flux- or
    radio-selected; lens light plays no part) or "galaxy" (finders run around galaxies, or not
    known). Only the first two could have selected a lens without a visible deflector."""
    text = " ".join(str(r[c]) for c in ("name", "entries", "refs", "system_type"))
    if _SUBMM_RADIO.search(text):
        return "submm/radio"
    if "lensedquasars" in text or "GQ" in str(r["system_type"]):
        return "quasar"
    return "galaxy"


def is_submm_radio(r) -> bool:
    text = " ".join(str(r[c]) for c in ("name", "entries", "refs", "simbad") if c in r.colnames)
    return bool(_SUBMM_RADIO.search(text))


def xmatch(t: Table, cat: str, radius: float) -> Table:
    from astropy import units as u
    from astroquery.xmatch import XMatchClass

    up = Table({"sid": np.asarray(t["system_id"], str), "ra0": t["ra"], "dec0": t["dec"]})
    xm = XMatchClass()
    xm.TIMEOUT = 600
    return xm.query(
        cat1=up, cat2=cat, max_distance=radius * u.arcsec, colRA1="ra0", colDec1="dec0", cache=False
    )


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


def contact_sheet(t: Table, cut_dir: Path, out_png: Path, ncol: int = 8) -> None:
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
            r = row["theta_e_used"] / 0.262
            ax.add_patch(plt.Circle((c, c), 1.5 / 0.262, fill=False, color="w", lw=0.6, ls=":"))
            ax.add_patch(plt.Circle((c, c), r, fill=False, color="c", lw=0.5))
        ax.set_title(
            f"{row['system_id']} {str(row['name'])[:16]}\n{row['verdict'][:26]}", fontsize=6
        )
    fig.suptitle(
        "W1/W2 lens-catalogue survivors: LS DR10 grz, 21'' (dotted 1.5'', cyan theta_E)", fontsize=8
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def cmd_vet(args, p: Params) -> None:
    """Cheapest-first ordinary tests on every system without a visible deflector."""
    out = Path(args.out)
    s = Table.read(out / "no_visible_deflector.ecsv")
    parts = [Table.read(f) for f in sorted((out / "tractor_cache").glob("tractor_*.ecsv"))]
    tr = vstack([t for t in parts if len(t)], metadata_conflicts="silent")
    tr["type"] = tr["type"].astype(str)
    for c in ("mag_z", "galdepth_z"):
        tr[c] = np.asarray(np.ma.filled(np.ma.asarray(tr[c], float), np.nan))
    from scipy.spatial import cKDTree

    xyz_t = lenscats._unit(np.asarray(tr["ra"], float), np.asarray(tr["dec"], float))
    tree = cKDTree(xyz_t)
    xyz_s = lenscats._unit(np.asarray(s["ra"], float), np.asarray(s["dec"], float))
    ext = np.isin(np.asarray(tr["type"]), lenscats.EXTENDED_TYPES)
    magz = np.asarray(tr["mag_z"], float)
    rad = np.pi / 180 / 3600
    # bright enough to be the deflector: typical required magnitude + faint_sigma * rms
    lim = np.asarray(s["req_mag_z_typical"], float) + p.faint_sigma * s.meta["calibration"]["rms"]
    lim = np.where(np.isfinite(lim), lim, np.inf)
    n_gal, n_pt = [], []
    for x, m in zip(xyz_s, lim, strict=True):
        near = np.asarray(tree.query_ball_point(x, p.offset_radius * rad), int)
        sep = np.linalg.norm(xyz_t[near] - x, axis=1) / rad if len(near) else np.zeros(0)
        bright = magz[near] <= m
        n_gal.append(int(np.sum(ext[near] & bright)))
        n_pt.append(int(np.sum(~ext[near] & (sep <= p.theta_e_max))))
    s["n_gal_offset"] = n_gal  # extended, bright enough, within offset_radius
    s["n_point_3as"] = n_pt  # point sources (any magnitude) within theta_e_max
    # maskbits at the nearest source
    d, k = cKDTree(xyz_t).query(
        lenscats._unit(np.asarray(s["ra"], float), np.asarray(s["dec"], float))
    )
    mb = np.asarray(np.ma.filled(np.ma.asarray(tr["maskbits"], int), 0))[k]
    s["masked"] = (mb & MASK_BAD) > 0
    s["nearest_sep"] = d * 3600 * 180 / np.pi
    s["nearest_type"] = np.asarray(tr["type"])[k]
    s["nearest_ref_cat"] = np.asarray(np.ma.filled(np.ma.asarray(tr["ref_cat"]).astype(str), ""))[k]
    s["simbad_galaxy_3as"] = np.full(len(s), "", dtype=object)
    for key, (cat, rad) in VET_XMATCH.items():
        try:
            m = xmatch(s, cat, rad)
        except Exception as e:  # noqa: BLE001 - record and continue; the column says "not run"
            print(f"xmatch {key} failed: {e}", flush=True)
            s[key] = np.full(len(s), "not run", dtype=object)
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
        s[key] = np.array([hit.get(str(x), "") for x in s["system_id"]], dtype=object)
        if key == "simbad":
            s["simbad_galaxy_3as"] = [gal.get(str(x), "") for x in s["system_id"]]
            s["simbad_cluster"] = [clu.get(str(x), "") for x in s["system_id"]]
        print(f"xmatch {key}: {len(m)} rows, {len(hit)} systems", flush=True)

    def ran(v) -> bool:
        return str(v) not in ("", "not run")

    tests = (  # cheapest first; the verdict is the first that applies, ``flags`` lists all
        ("catalogue position error (name vs RA/Dec)", lambda r: r["name_offset"] > 5),
        ("catalogue position rounded to > 1.5''", lambda r: r["pos_quantum"] > p.lens_radius),
        (
            "Euclid position may be off by up to 10'' (README)",
            lambda r: r["catalogues"] == "euclid_q1",
        ),
        ("submm/radio centroid position (arcsec uncertainty)", is_submm_radio),
        ("galaxy bright enough for the lens within 5''", lambda r: r["n_gal_offset"] > 0),
        ("LS maskbits: bright star / large galaxy / cluster", lambda r: bool(r["masked"])),
        ("point sources within 3'': lens blended or PSF-typed", lambda r: r["n_point_3as"] > 0),
        (
            "SIMBAD galaxy within 3'' (literature lens galaxy)",
            lambda r: ran(r["simbad_galaxy_3as"]),
        ),
        (
            "lens redshift published: deflector observed, not at this position",
            lambda r: bool(r["lens_z_known"]),
        ),
        (
            "cluster: redMaPPer/WHL within 2' or SIMBAD cluster within 5''",
            lambda r: ran(r["redmapper"]) or ran(r["whl"]) or ran(r["simbad_cluster"]),
        ),
        ("required lens below LS depth (conservative)", lambda r: not r["detectable"]),
        ("literature: not a lens / explained", lambda r: str(r["name"]) in LITERATURE),
    )
    flags = [[name for name, f in tests if f(r)] for r in s]
    s["flags"] = np.array(["; ".join(f) for f in flags], dtype=object)
    s["verdict"] = np.array(
        [f[0] if f else "unexplained by catalogue tests" for f in flags], dtype=object
    )
    s.meta.update(provenance=schema.Provenance.DERIVED.value)
    s.write(out / "vetted.ecsv", overwrite=True)
    vals, cnt = np.unique(s["verdict"], return_counts=True)
    print(dict(zip(vals.tolist(), cnt.tolist(), strict=True)))
    # cutouts: systems left with no flag but the conservative depth test, where a typical ordinary
    # lens would be detected (eyes needed), then a random control sample of the flagged rest
    eyes = ("unexplained by catalogue tests", "required lens below LS depth (conservative)")
    det = np.asarray(s["detectable_typical"], bool) & np.isin(s["verdict"], eyes)
    rest = np.flatnonzero(~det)
    rng = np.random.default_rng(56)
    ctrl = rng.choice(rest, size=min(len(rest), args.n_control), replace=False)
    look = vstack([s[det], s[np.sort(ctrl)]])[: args.max_cutouts]
    print(f"cutouts: {int(det.sum())} to inspect + {len(ctrl)} controls", flush=True)
    cut = out / "cutouts"
    cut.mkdir(exist_ok=True)
    for r in look:
        fetch_cutout(float(r["ra"]), float(r["dec"]), cut / f"{r['system_id']}.jpg")
    for i in range(0, len(look), 64):
        contact_sheet(look[i : i + 64], cut, out / f"contact_sheet_{i // 64:02d}.png")
    summary = json.loads((out / "summary.json").read_text())
    summary["vetting"] = dict(zip(vals.tolist(), [int(c) for c in cnt], strict=True))
    summary["limits"] = limits(out, s, p)
    (out / "summary.json").write_text(json.dumps(summary, indent=2))


UNDECIDED = (  # the LS test says nothing about the deflector for these
    "catalogue position error (name vs RA/Dec)",
    "catalogue position rounded to > 1.5''",
    "Euclid position may be off by up to 10'' (README)",
    "submm/radio centroid position (arcsec uncertainty)",
    "LS maskbits: bright star / large galaxy / cluster",
    "point sources within 3'': lens blended or PSF-typed",
)
POISSON95_ZERO = 2.996  # 95 % upper limit on a Poisson mean for zero events


def limits(out: Path, vetted: Table, p: Params) -> dict:
    """Upper limits on the dark-deflector fraction among tested galaxy-scale lenses (derived).

    Tested: LS-covered systems where the test is decisive, i.e. a deflector was found at the
    position, within the offset radius, in the literature or in a cluster catalogue, or none was
    found although the position is good, unmasked and not blended (the ``UNDECIDED`` verdicts are
    excluded). A dark lens is counted only if ``detectable`` (an ordinary lens would have been
    seen); with zero unexplained systems the 95 % limit is 2.996 / N."""
    allsys = Table.read(out / "systems.ecsv")
    cov = np.asarray(allsys["covered"], bool)
    und = dict(zip(vetted["system_id"], np.isin(vetted["verdict"], UNDECIDED), strict=True))
    decisive = cov & ~np.array([und.get(x, False) for x in allsys["system_id"]])
    unexplained = set(vetted["system_id"][vetted["verdict"] == "unexplained by catalogue tests"])
    k = sum(1 for x in allsys["system_id"] if x in unexplained)
    res = {"unexplained": k, "covered": int(cov.sum()), "decisive": int(decisive.sum())}
    sel = np.asarray(allsys["selection"])
    for variant, col in (("typical", "detectable_typical"), ("conservative", "detectable")):
        det = decisive & np.asarray(allsys[col], bool)
        for name, m in (
            ("all", np.ones(len(sel), bool)),
            ("lens_light_independent", sel != "galaxy"),
        ):
            n = int((det & m).sum())
            res[f"{variant}_{name}_N"] = n
            res[f"{variant}_{name}_f95"] = POISSON95_ZERO / n if n and k == 0 else None
    res["selection_counts_covered"] = {
        c: int((cov & (sel == c)).sum()) for c in ("galaxy", "quasar", "submm/radio")
    }
    res["selection_counts_decisive"] = {
        c: int((decisive & (sel == c)).sum()) for c in ("galaxy", "quasar", "submm/radio")
    }
    return res


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["screen", "vet"])
    ap.add_argument("--n-control", type=int, default=32, help="random other cutouts")
    ap.add_argument("--max-cutouts", type=int, default=600)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    args.out = args.out or paths.outputs_dir() / "w12_lenscats"
    p = Params()
    {"screen": cmd_screen, "vet": cmd_vet}[args.cmd](args, p)


if __name__ == "__main__":
    main()
