"""Stage the inputs ``inject_radial.load_field`` reads for the W1 fields into the data root.

For each field (``inject_shear.field_spec``) this fetches, into the paths ``load_field`` expects:

- the MAST F200W pipeline catalogue ``spec["catalog"]`` from the public S3 bucket
  (``stpubdata``), pinned by the sha256 in ``data/manifests/*.ecsv`` or, for catalogues not in a
  manifest, the value recorded in SOURCES.md / ``results/arcs/*/summary.json``;
- the photo-z file ``spec["photoz"][1]``: a DJA eazy zout (a plain file, or one member streamed
  out of a ``.photoz.tar.gz`` that is never stored), a CANUCS DR1 photometry catalogue, or
  (MACS0717) a zout that ``astrodeep_zout.py`` writes from a live VizieR query, checked against
  a sha256 of its column values (``astrodeep_zout.content_sha256``), not of the FITS bytes.

Downloads reuse ``photometry.fetch_catalog`` and ``orphan_pairs.fetch_tar_member`` (verified,
``.part`` files removed on failure); the verified file is then moved to the expected name. Files
without a recorded sha256 are fetched unverified and listed at the end with the sha256 they got.
Lens models load from the repository or their pinned URLs, except maps pinned as zip members
(RXJ2129), which are staged with ``arcs_rxj2129.stage_members``. The Gaia star files are made by
``load_field`` itself.

    python scripts/w1_stage_inputs.py --fields abell2744 macs0416 [--check OUT_DIR]

``--check`` then runs ``load_field`` on each field (Gaia star files go to ``OUT_DIR/<field>``)
and prints the number of shape rows.
"""

from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))

import arcs_rxj2129  # noqa: E402
import astrodeep_zout  # noqa: E402
import inject_shear as ish  # noqa: E402
import orphan_pairs as op  # noqa: E402

from jwst_anomaly import acquire, paths  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

ir = ish.ir
S3 = "https://stpubdata.s3.amazonaws.com/jwst/public/{prog}/L3/t/{obs}/{name}"
DJA = op._DJA
_CANUCS_URL = (
    "https://archive.stsci.edu/hlsps/canucs/dr1/{0}/clu/"
    "hlsp_canucs_jwst-hst_multi_{0}-clu_multi_v1_photometry-cat.fits.gz"
)
MAX_BYTES = 100_000_000  # every file and tarball here is < 80 MB (SOURCES.md)

# F200W catalogues not in data/manifests: sha256 as recorded in SOURCES.md (MACS0417, MACS1423,
# MACS0647) and results/arcs/rxj2129/summary.json (RXJ2129).
CATALOG_SHA256 = {
    "jw01208-o006_t003_nircam_clear-f200w_cat.ecsv": (
        "e55dffb990a8e8124e37658933a47952a77b357b36dc5c8f523a5b487d423c24"
    ),
    "jw01208-o010_t005_nircam_clear-f200w_cat.ecsv": (
        "1a968574d0f206c72a74554c1b3cbbf71e9a5219a70a553a1458e80dd4d5bfa3"
    ),
    "jw01433-o010_t001_nircam_clear-f200w_cat.ecsv": (
        "473e32d23c3fbc857e89d6b48196130d77670be2b9df2b4dfae5df82e71e7a7e"
    ),
    "jw02767-o002_t001_nircam_clear-f200w_cat.ecsv": (
        "2ed1a2c7666734483974e568e1c3b332eea223ea487b4cd279dd89b4e9261f53"
    ),
    # CANUCS NCF F150W (results/w1_shear_ncf), as fetched 2026-10-18
    "jw01208-o023_t002_nircam_clear-f150w_cat.ecsv": (
        "6ac987242db0435afb0aac549b825d12a7f67868e93cf23e63aef0dada8ec38c"
    ),
    "jw01208-o029_t004_nircam_clear-f150w_cat.ecsv": (
        "e331b754f054cf3c31eb84096ff7dd3131b2aeb6f595f711953244ddad2c672d"
    ),
    "jw01208-o020_t001_nircam_clear-f150w_cat.ecsv": (
        "6be0bacd8cd437bbf557076a099888770c301ea2b81b835ddde84f12b7b9328c"
    ),
    "jw01208-o026_t003_nircam_clear-f150w_cat.ecsv": (
        "b3d71250aec8bbe965c1c4908958e5d6366af7c8e25abf587006de65cee8383e"
    ),
    "jw01208-o032_t005_nircam_clear-f150w_cat.ecsv": (
        "fb2f118625ae5a4a8c1107e546c5c19dfb7f2c4d1125d4c24d003ff4677898d2"
    ),
}


def _tar(stem: str, sha: str | None, member_sha: str) -> dict:
    member = f"{stem}.eazypy.zout.fits"
    url = f"{DJA}{stem}.photoz.tar.gz"
    return {"url": url, "sha256": sha, "member": member, "member_sha256": member_sha}


# Photo-z files by the basename load_field expects; url / sha256 as recorded in SOURCES.md,
# docs/fields/*.md and orphan_pairs.FIELDS (``member``: streamed out of a DJA tarball).
PHOTOZ = {
    "smacs0723-grizli-v7.4-fix.eazypy.zout.fits": _tar(
        "smacs0723-grizli-v7.4-fix",
        "13b5b6f40048b1f474827d93ca52dc94f0200261187c59ca397d70382c8c2821",
        "82db8910962f3050ddb9f8209aafa487fb30f3cf46551d57ef3824d7d2c355ea",
    ),
    "elgordo-grizli-v7.0-fix.eazypy.zout.fits": _tar(
        "elgordo-grizli-v7.0-fix",
        "9a4a39cee4b154c34040e05591a151597e85d05d0f580cbabdd15519488baade",
        "efd051e4e688ed46f7279e08cab19ec4179e114705d38f879ffb8feff230edd0",
    ),
    "abell2744clu-grizli-v7.2-fix.eazypy.zout.fits": op.FIELDS["abell2744"]["photoz"],
    "abells1063-grizli-v7.5-fix.eazypy.zout.fits": _tar(
        "abells1063-grizli-v7.5-fix",
        "088a1954e357d96e9e3b0ab4af1f474eb9d35676774c1ad990f40140e81b1319",
        "03bc4c289982dc50a4540305b8858b279aca7c7cdb64c809d69a855e631b6428",
    ),
    "sunrise-grizli-v7.5-fix.eazypy.zout.fits": _tar(
        "sunrise-grizli-v7.5-fix",
        "a99b604f7d7418271d94f4665e2e8ab2839f67edf53f28d93e0b9a02cfffc6f1",
        "672d4bdf31926209633c1a7feeb436bab3e5b08e847a079521820f3409f548e0",
    ),
    "rxj2129-grizli-v7.0-fix.eazypy.zout.fits": _tar(
        "rxj2129-grizli-v7.0-fix",
        None,  # tarball sha256 not recorded (SOURCES.md pins only the member)
        "651c0b24f09bafa27c36a2285f081d42884be3c6dcbbb1708014339afa726217",
    ),
    "macs0647-grizli-v7.0-fix.eazypy.zout.fits": _tar(
        "macs0647-grizli-v7.0-fix",
        "d56246b571ccfc78a54338a7ab16009fa8b73b2cdbcfefba4a2617c0058cdf0c",
        "ceb199c61f1ca67df6afc34df6ed7ad61267f59ccaf27e9a2fd1e64c2662c04e",
    ),
}
# No DJA photo-z for MACS0717 (2026-10-17): ASTRODEEP-FF from VizieR, written by astrodeep_zout.py
PHOTOZ["m0717clu_zout.fits"] = {
    "astrodeep": "m0717",
    "sha256": "ef1a24d3440e70dd537db605b94e8a5f66080a3541878384aa27d5c3a4de7ca0",
}
for _f, _sha in (
    ("macs0417", "ca6b48aa6404444d2e0b8965a35376ff032917071ea64cbf6306d06ee57d90dd"),
    ("macs1423", "c257d8dff95319f4eea981cc838269c04aff0a79bac603c7ed64100deac622b0"),
):
    PHOTOZ[Path(ir._CANUCS.format(_f)).name] = {"url": _CANUCS_URL.format(_f), "sha256": _sha}
for _f in ("macs0416", "macs1149", "abell370"):
    _url, _sha = op.FIELDS[_f]["catalogue"]
    PHOTOZ[_url.rsplit("/", 1)[-1]] = {"url": _url, "sha256": _sha}
for _spec in op.DEEP_FIELDS.values():  # CANUCS NCF photometry (results/w1_shear_ncf)
    _url, _sha = _spec["catalogue"]
    PHOTOZ.setdefault(_url.rsplit("/", 1)[-1], {"url": _url, "sha256": _sha})


class _Unpinned(str):
    """Stands in for a sha256 that was never recorded: every digest compares equal, so the fetch
    helpers download without verifying. The sha256 actually obtained is reported afterwards."""

    def lower(self):
        return self

    def __eq__(self, other):
        return True

    def __ne__(self, other):
        return False

    __hash__ = str.__hash__


def manifest_sha256() -> dict[str, str]:
    """productFilename -> sha256 over every download manifest in ``data/manifests``."""
    out = {}
    for path in sorted(paths.manifests_dir().glob("*.ecsv")):
        try:
            t = acquire.read_manifest(path)
        except ValueError:  # product listings and other tables: no sha256 column
            continue
        out.update(zip(t["productFilename"], t["sha256"], strict=True))
    return out


def catalog_source(rel: str, pins: dict[str, str]) -> dict:
    name = Path(rel).name
    obs = name.split("_")[0]  # jw01208-o006
    prog, o = obs.split("-")
    return {"url": S3.format(prog=prog, obs=o, name=name), "sha256": pins.get(name)}


def stage(target: Path, src: dict) -> dict:
    """Fetch ``src`` to ``target`` (verified where a sha256 is recorded); returns a log row."""
    pinned = src["member_sha256"] if "member" in src else src["sha256"]
    row = {"path": str(target), "url": src["url"], "pinned": pinned is not None, "bytes": 0}
    if "member" in src:
        row["archive_pinned"] = src["sha256"] is not None
    if target.exists() and (pinned is None or acquire.sha256_file(target) == pinned):
        row.update(status="present", sha256=acquire.sha256_file(target))
        return row
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        if "member" in src:
            spec = {
                **src,
                "sha256": src["sha256"] or _Unpinned(),
                "member_sha256": pinned or _Unpinned(),
                "max_bytes": MAX_BYTES,
            }
            got = op.fetch_tar_member(spec, cache_dir=target.parent)
        else:
            sha = pinned or _Unpinned()
            got = fetch_catalog(src["url"], sha, cache_dir=target.parent, max_bytes=MAX_BYTES)
        os.replace(got, target)
        row["bytes"] = target.stat().st_size  # the kept file (a tar member: the member only)
    except Exception as exc:  # noqa: BLE001 - report and continue with the other files
        row.update(status=f"FAILED: {type(exc).__name__}: {exc}")
        return row
    row.update(status="fetched", sha256=acquire.sha256_file(target))
    return row


def stage_zip_members(files: dict[str, tuple[str, str]]) -> dict:
    """Lens-model maps pinned as ``<zip url>#<member>`` (RXJ2129): range-read out of the zip
    by ``arcs_rxj2129.stage_members`` into the cache ``load_model`` reads. ``bytes`` is the
    staged (uncompressed) size, an upper bound on the member bytes read."""
    cache = paths.cache_dir() / "external"
    targets = [cache / f"{sha[:12]}_{u.rsplit('/', 1)[-1]}" for u, sha in files.values()]
    missing = [t for t in targets if not t.exists()]
    row = {"path": ", ".join(t.name for t in targets), "url": next(iter(files.values()))[0]}
    row.update(pinned=True, bytes=0)
    try:
        arcs_rxj2129.stage_members(files)
    except Exception as exc:  # noqa: BLE001 - report and continue with the other files
        return row | {"status": f"FAILED: {type(exc).__name__}: {exc}"}
    row["bytes"] = sum(t.stat().st_size for t in missing)
    return row | {"status": "fetched" if missing else "present"}


def _content_sha256(path: Path) -> str | None:
    try:
        return astrodeep_zout.content_sha256(Table.read(path))
    except Exception:  # noqa: BLE001 - an unreadable file is re-fetched
        return None


def stage_astrodeep(target: Path, src: dict) -> dict:
    """Write the ASTRODEEP-FF zout (``astrodeep_zout.write``, from VizieR) and check its values.

    The pin is ``astrodeep_zout.content_sha256`` (column values), so a new astropy version that
    serialises identical content differently still passes.
    """
    url = f"vizier:{astrodeep_zout.CATALOG}/{src['astrodeep']}"
    row = {"path": str(target), "url": url, "pinned": True, "bytes": 0}
    if target.exists() and _content_sha256(target) == src["sha256"]:
        return {**row, "status": "present", "sha256": src["sha256"]}
    try:
        got = astrodeep_zout.write(src["astrodeep"], target, expect=src["sha256"])
    except Exception as exc:  # noqa: BLE001 - report and continue with the other files
        return {**row, "status": f"FAILED: {type(exc).__name__}: {exc}"}
    return {**row, "status": "fetched", "sha256": got, "bytes": target.stat().st_size}


def run_job(job: tuple[str, Path | None, dict]) -> dict:
    _, target, src = job
    if "astrodeep" in src:
        return stage_astrodeep(target, src)
    return stage_zip_members(src["zip_members"]) if target is None else stage(target, src)


def plan(fields: list[str]) -> list[tuple[str, Path | None, dict]]:
    root = paths.data_root()
    pins = {**manifest_sha256(), **CATALOG_SHA256}
    jobs = []
    for name in fields:
        spec = ish.field_spec(name)
        jobs.append((name, root / spec["catalog"], catalog_source(spec["catalog"], pins)))
        if spec.get("photoz"):
            rel = spec["photoz"][1]
            src = PHOTOZ.get(Path(rel).name)
            if src is None:
                raise KeyError(f"{name}: no recorded source for {rel}")
            jobs.append((name, root / rel, src))
        files = ir.lc.MODELS[spec["model"]].get("files", {}) if spec["model"] else {}
        if any("#" in str(v[0]) for v in files.values()):
            jobs.append((name, None, {"zip_members": files}))
    return jobs


def check(fields: list[str], out: Path) -> dict[str, str]:
    res = {}
    for name in fields:
        d = out / name
        d.mkdir(parents=True, exist_ok=True)
        try:
            _, shapes, extra, _, photoz = ir.load_field(name, d, ish.field_spec(name))
        except Exception as exc:  # noqa: BLE001
            res[name] = f"FAILED: {type(exc).__name__}: {exc}"
        else:
            stars = "-" if extra is None else len(extra)
            res[name] = f"OK shapes={len(shapes)} photoz={bool(photoz)} gaia_stars={stars}"
        print(f"load_field {name}: {res[name]}", flush=True)
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fields", nargs="+", default=ish.FIELD_NAMES, choices=ish.FIELD_NAMES)
    ap.add_argument("--workers", type=int, default=8, help="concurrent downloads")
    ap.add_argument("--check", type=Path, help="run load_field per field, outputs here")
    args = ap.parse_args(argv)

    jobs = plan(args.fields)
    with ThreadPoolExecutor(args.workers) as pool:
        rows = list(pool.map(run_job, jobs))
    failed = set()
    for (name, _, _), row in zip(jobs, rows, strict=True):
        print(f"{name:11s} {row['status']:8s} {row['bytes']:>11,d} B  {row['path']}", flush=True)
        if row["status"].startswith("FAILED"):
            failed.add(name)
    unpinned = [r for r in rows if not r["pinned"] or r.get("archive_pinned") is False]
    for r in unpinned:
        what = "file" if not r["pinned"] else "archive (member verified)"
        print(f"no recorded sha256 for the {what}: {r['url']} (got {r.get('sha256', '-')})")
    print(f"total streamed: {sum(r['bytes'] for r in rows):,d} B")
    if args.check:
        res = check([f for f in args.fields if f not in failed], args.check)
        failed |= {f for f, s in res.items() if s.startswith("FAILED")}
    if failed:
        print(f"failed fields: {sorted(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
