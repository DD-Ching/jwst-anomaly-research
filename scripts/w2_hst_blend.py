"""W2 at HST resolution for the D-056 "blended" and "too close" lensed quasars (line C, 2026-10-21).

Stages (run in order; each writes small tracked tables to ``results/w2_hst_blend/``):

    python scripts/w2_hst_blend.py classes  --work <dir>   # rebuild D-056 class membership
    python scripts/w2_hst_blend.py coverage --work <dir>   # HST imaging per system (MAST)
    python scripts/w2_hst_blend.py select                  # frozen selection rule -> selection.ecsv

``classes`` rebuilds only what the selection needs: the galaxy-scale quasar and radio systems in
the DR10 footprint and their pair-test status (``lenscats.quasar_pair_test`` via
``w12_lenscats.deflector_test``). "blended" and "too close" depend only on each system's own
Tractor sources, not on the Faber-Jackson calibration, so Tractor is queried for these systems
only and the calibration is not rebuilt (the deflector / faint galaxy / none labels of this
rebuild are therefore not the D-056 ones and are not used). Outputs are ``derived``; catalogue
values are ``observed``. Nothing here is evidence of exotic physics.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w12_lenscats as lc  # noqa: E402

from jwst_anomaly import lenscats, paths, schema  # noqa: E402

RES = paths.repo_root() / "results" / "w2_hst_blend"
UNDECIDABLE = ("blended", "too close")


def log_download(log: Path | None, url: str, path: Path, note: str = "") -> None:
    """One line per download: UTC time, URL, bytes, sha256 (the line's reproducibility record)."""
    if log is None:
        return
    data = Path(path).read_bytes()
    sha = hashlib.sha256(data).hexdigest()
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(log, "a") as fh:
        fh.write(f"{stamp}\t{url}\t{len(data)}\t{sha}\t{Path(path).name}\t{note}\n")


def cmd_classes(args) -> None:
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    log = Path(args.log) if args.log else None
    p = lc.Params()
    f = lc.fetch_all()
    for k, v in lenscats.CATALOGUES.items():
        log_download(log, v["url"], f[k], f"catalogue {k} (sha256 pinned by lenscats)")
    systems, sizes = lc.load_systems(p)
    galaxy = ~np.asarray(systems["any_group_scale"], bool) & ~(
        np.asarray(systems["theta_e"], float) > p.theta_e_max
    )
    gs = systems[galaxy]
    gs["selection"] = [lc.selection_class(r) for r in gs]
    sens = gs[np.isin(np.asarray(gs["selection"]), ("quasar", "radio"))]
    had_bricks = (work / "bricks_dr10_south.csv").exists()
    bricks, bmeta = lc.load_bricks(work)
    if not had_bricks:
        log_download(log, lc.TAP, work / "bricks_dr10_south.csv", "ls_dr10.bricks_s TAP")
    foot = lenscats.brick_coverage(sens, bricks)
    cov = sens[np.asarray(foot["covered"], bool)]
    cov["depth_z"] = np.asarray(foot["depth_z"], float)[np.asarray(foot["covered"], bool)]
    cache = work / "tractor_cache"
    before = set(cache.glob("tractor_*.ecsv")) if cache.exists() else set()
    lc.query_tractor(cov, p, cache)
    for fn in sorted(set(cache.glob("tractor_*.ecsv")) - before):
        log_download(log, lc.TAP, fn, "ls_dr10.tractor box batch (file = re-serialised response)")
    src = lc.load_sources(work)
    images = lenscats.pair_images(cov, src, p.image_radius)
    # blended / too close do not depend on the deflector magnitude limit; -inf keeps every
    # candidate "faint", so deflector/none labels here are NOT the D-056 ones (not used)
    cov["defl_mag_max"] = np.full(len(cov), -np.inf)
    test = lc.deflector_test(cov, src, p, images)
    out = Table(
        {
            "system_id": np.asarray(cov["system_id"], str),
            "name": np.asarray(cov["name"], str),
            "ra": np.asarray(cov["ra"], float),
            "dec": np.asarray(cov["dec"], float),
            "selection": np.asarray(cov["selection"], str),
            "status": np.asarray(test["test_status"], str),
            "ls_n_psf": np.asarray(images["n_images"], int),
            "ls_sep": np.asarray(images["sep"], float),
            "theta_e": np.asarray(cov["theta_e"], float),
            "z_lens": np.asarray(cov["z_lens"], float),
            "lens_z_known": np.asarray(cov["lens_z_known"], bool),
            "z_source": np.asarray(cov["z_source"], float),
            "catalogues": np.asarray(cov["catalogues"], str),
            "refs": np.asarray(cov["refs"], str),
        }
    )
    status = np.asarray(out["status"])
    sel = np.asarray(out["selection"])
    counts = {
        "galaxy_scale": int(galaxy.sum()),
        "quasar_radio_galaxy_scale": len(sens),
        "quasar_radio_covered": len(cov),
        "status_by_selection": {
            c: {s: int(((sel == c) & (status == s)).sum()) for s in np.unique(status)}
            for c in ("quasar", "radio")
        },
        "quasar_blended_or_too_close": int(
            ((sel == "quasar") & np.isin(status, UNDECIDABLE)).sum()
        ),
        "d056_expected": {
            "quasar": {"blended": 242, "too close": 65},
            "radio": {"blended": 1, "too close": 1},
        },
        "bricks": bmeta,
        "note": "deflector/faint galaxy/none labels use defl_mag_max = -inf (no FJ calibration) "
        "and are not the D-056 labels; only blended/too close are used",
    }
    out = out[np.isin(status, UNDECIDABLE) & (sel == "quasar")]
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="w12_lenscats.deflector_test on lenscat 1.1.3 + Euclid Q1 + SuGOHI (sha256-pinned), "
        "LS DR10 bricks_s and Tractor boxes (Data Lab TAP); scripts/w2_hst_blend.py classes",
    )
    RES.mkdir(parents=True, exist_ok=True)
    out.write(RES / "classes.ecsv", overwrite=True)
    (RES / "classes_summary.json").write_text(json.dumps(counts, indent=2) + "\n")
    print(json.dumps(counts, indent=1))


MAST_TAP = "https://mast.stsci.edu/vo-tap/api/v0.1/caom/sync"
COVERAGE_FILTERS = ("F814W", "F160W")  # the chain's band (D-064 tools) and a redder second band
XMATCH_RADIUS_DEG = 0.03  # observation centre within 108": every ACS/WFC or WFC3/UVIS pointing that
# can contain the target (fields 202" / 162"); containment is then tested on s_region
COVERAGE_INSTRUMENTS = ("ACS/WFC", "WFC3/UVIS", "WFC3/IR")


def point_in_region(ra: float, dec: float, s_region: str) -> bool:
    """Whether (ra, dec) lies inside an STC-S ``POLYGON`` (any of several) of ``s_region``, tested
    in the tangent plane at the point (fields are < 4', so the projection error is negligible)."""
    toks = str(s_region).replace("ICRS", " ").replace("J2000", " ").split()
    polys, cur = [], None
    for t in toks:
        u = t.upper()
        if u in ("POLYGON", "CIRCLE", "BOX", "POSITION", "UNION", "INTERSECTION"):
            if cur:
                polys.append(cur)
            cur = [] if u == "POLYGON" else None
            continue
        if cur is not None:
            try:
                cur.append(float(t))
            except ValueError:
                pass
    if cur:
        polys.append(cur)
    for v in polys:
        if len(v) < 6 or len(v) % 2:
            continue
        x = (((np.asarray(v[0::2]) - ra + 180) % 360) - 180) * np.cos(np.radians(dec))
        y = np.asarray(v[1::2]) - dec
        inside = False
        j = len(x) - 1
        for i in range(len(x)):  # ray casting from the origin along +x
            if (y[i] > 0) != (y[j] > 0):
                xc = x[i] + (0 - y[i]) * (x[j] - x[i]) / (y[j] - y[i])
                if xc > 0:
                    inside = not inside
            j = i
        if inside:
            return True
    return False


def cmd_coverage(args) -> None:
    """HST imaging at every class member from one MAST CAOM TAP query per filter (all calibrated
    level-3 HST images in that filter; the TAP rejects OR-ed geometry, and the position services
    take ~6 s per object), matched locally: observation centre within ``XMATCH_RADIUS_DEG``, then
    the position inside ``s_region``. ``observed``."""
    import io

    import requests
    from astropy.io.votable import parse_single_table

    log = Path(args.log) if args.log else None
    cls = Table.read(RES / "classes.ecsv")
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    rows = []
    cra, cdec = np.radians(np.asarray(cls["ra"], float)), np.radians(np.asarray(cls["dec"], float))
    cxyz = np.c_[np.cos(cdec) * np.cos(cra), np.cos(cdec) * np.sin(cra), np.sin(cdec)]
    for filt in COVERAGE_FILTERS:
        q = (
            "SELECT obs_id, target_name, s_ra, s_dec, instrument_name, filters, t_exptime, "
            "proposal_id, s_region FROM dbo.ObsPointing WHERE obs_collection='HST' AND "
            f"dataproduct_type='image' AND calib_level=3 AND filters='{filt}'"
        )
        fn = work / f"caom_hst_{filt}.vot"
        if not fn.exists():
            err = ""
            for attempt in range(4):
                try:
                    r = requests.post(
                        MAST_TAP,
                        data={"REQUEST": "doQuery", "LANG": "ADQL", "MAXREC": "500000", "QUERY": q},
                        timeout=900,
                    )
                    if r.status_code == 200 and b'QUERY_STATUS" value="ERROR' not in r.content:
                        fn.write_bytes(r.content)
                        log_download(log, MAST_TAP, fn, f"CAOM TAP: {q}")
                        break
                    err = f"{r.status_code} {r.text[:300]}"
                except requests.RequestException as e:  # retried, never recorded as "no data"
                    err = str(e)[:300]
                time.sleep(10 * (attempt + 1))
            else:
                raise RuntimeError(f"MAST TAP failed: {err}")
        t = parse_single_table(str(fn)).to_table()
        # this service flags OVERFLOW even on complete results: compare with a COUNT instead
        cq = q.replace(q[q.index("SELECT ") + 7 : q.index(" FROM")], "COUNT(*)")
        rc = requests.post(
            MAST_TAP, data={"REQUEST": "doQuery", "LANG": "ADQL", "QUERY": cq}, timeout=300
        )
        n_all = int(parse_single_table(io.BytesIO(rc.content)).to_table()[0][0])
        if n_all != len(t):
            raise RuntimeError(f"{fn.name}: {len(t)} rows of {n_all} (truncated)")
        ok = np.isin(np.asarray(t["instrument_name"]).astype(str), COVERAGE_INSTRUMENTS)
        t = t[ok & np.isfinite(np.asarray(t["s_ra"], float))]
        ora, odec = (
            np.radians(np.asarray(t["s_ra"], float)),
            np.radians(np.asarray(t["s_dec"], float)),
        )
        oxyz = np.c_[np.cos(odec) * np.cos(ora), np.cos(odec) * np.sin(ora), np.sin(odec)]
        from scipy.spatial import cKDTree

        tree = cKDTree(oxyz)
        rad = np.radians(XMATCH_RADIUS_DEG)
        for i, idx in enumerate(tree.query_ball_point(cxyz, rad)):
            for j in idx:
                x = t[j]
                d = float(np.degrees(np.arccos(np.clip(cxyz[i] @ oxyz[j], -1, 1))) * 3600)
                rows.append(
                    {
                        "system_id": str(cls["system_id"][i]),
                        "obs_id": str(x["obs_id"]),
                        "instrument": str(x["instrument_name"]),
                        "filters": str(x["filters"]),
                        "t_exptime": float(x["t_exptime"]),
                        "proposal_id": str(x["proposal_id"]),
                        "target_name": str(x["target_name"]),
                        "dist_arcsec": d,
                        "contains": point_in_region(
                            float(cls["ra"][i]), float(cls["dec"][i]), str(x["s_region"])
                        ),
                    }
                )
    t = Table(rows=rows) if rows else Table(names=("system_id", "obs_id"))
    t.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=f"MAST CAOM TAP dbo.ObsPointing, HST calib_level 3 images in {COVERAGE_FILTERS}, "
        f"instruments {COVERAGE_INSTRUMENTS}; centre within {XMATCH_RADIUS_DEG} deg; "
        "contains = position inside s_region",
    )
    t.write(RES / "hst_coverage.ecsv", overwrite=True)
    c = t[np.asarray(t["contains"], bool)] if len(t) else t
    print(len(t), "observations;", len(set(c["system_id"])), "systems contained")


def cmd_select(args) -> None:
    import w2_hst_blend_chain as ch

    p = ch.Params()
    sel = ch.select_systems(
        Table.read(RES / "classes.ecsv"), Table.read(RES / "hst_coverage.ecsv"), p
    )
    sel.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="w2_hst_blend_chain.select_systems on classes.ecsv and hst_coverage.ecsv",
    )
    sel.write(RES / "selection.ecsv", overwrite=True)
    sel.pprint(max_lines=40, max_width=200)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("classes", "coverage"):
        a = sub.add_parser(name)
        a.add_argument("--work", required=True)
        a.add_argument("--log", default="")
    sub.add_parser("select")
    args = ap.parse_args(argv)
    {"classes": cmd_classes, "coverage": cmd_coverage, "select": cmd_select}[args.cmd](args)


if __name__ == "__main__":
    main()
