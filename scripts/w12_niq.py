"""W1/W2 in rejected lensed-quasar pairs: the D-056 dark-deflector test on pairs that lens searches
rejected because no lens galaxy was seen (D-064).

Inputs (VizieR, pinned by sha256 in ``data/manifests/w12_niq_inputs.json``):
- Lemon et al. 2023 (MNRAS 520, 3305) Gaia-selected candidates: classes "UQP" (unclassified
  quasar pair: same redshift, no lens seen; "akin to NIQs") and "QSO pair" are the *rejected*
  sample; "lens", "quad" and "lens (?)" are the *control* sample of real lenses.
- SQLS (Inada et al. 2008, 2010, 2012) candidate tables: comments "no lens(ing) object" are
  rejected; "SDSS lens" / "known lens" are controls.

Each system goes through ``w12_lenscats``: DR10 brick coverage and depth, Tractor boxes, the
quasar pair test (two PSF images >= 2" apart, a deflector between them) and the required lens
magnitude from the D-056 Faber-Jackson calibration (fixed here, ``D056_FJ``). The control sample
measures the test's efficiency for ordinary lenses, which D-056 had to assume. Outputs are
``derived``; a "none" is an ordinary binary quasar or a lens below the depth until vetted, never
evidence of a dark lens.

    python scripts/w12_niq.py screen --out outputs/w12_niq
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w12_lenscats as w12  # noqa: E402

from jwst_anomaly import __version__, lenscats, paths, schema  # noqa: E402

MAX_ROWS = 5000
VIZIER = (
    "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source={}&-out.max="
    + str(MAX_ROWS)
    + "&-out.all"
)
# (VizieR table, release label); SQLS tables: DR3 (2008), DR5 (2010), DR7 (2012)
INPUTS = {
    "J/MNRAS/520/3305/table1": "Lemon2023",
    "J/AJ/135/496/table2": "SQLS-DR3",
    "J/AJ/135/496/table3": "SQLS-DR3",
    "J/AJ/140/403/table2": "SQLS-DR5",
    "J/AJ/140/403/table3": "SQLS-DR5",
    "J/AJ/143/119/table3": "SQLS-DR7",
    "J/AJ/143/119/table4": "SQLS-DR7",
}
HENNAWI = "J/AJ/131/1/binqso"  # binary-quasar catalogue for vetting (Hennawi et al. 2006)
LEMON_REJECTED = {"UQP", "UQP (?)", "QSO pair"}
LEMON_CONTROL = {"lens", "quad", "lens (?)"}
# D-056 calibration (LS z; 605 galaxy-selected lenses; docs/exotic_limits.md), reused unchanged
D056_FJ = lenscats.FJCalibration(a=20.41, k=0.67, slope=-10.0, rms=0.88, n=605, band="ls_z")
# lensed images share a colour (ASSUMPTION: |delta(g - z)| <= 0.5 allows microlensing, dust and
# variability; the tolerance of quasar_pair_test's further-image rule)
COLOUR_TOL = 0.5
# two redshifts differing by more than 0.01 (1 + z) (~3,000 km/s, beyond broad-line redshift
# errors; ASSUMPTION) belong to two quasars, not two images
DZ_TOL = 0.01
# catalogue positions are one image (SQLS: the SDSS quasar), so the second image lies inside the
# 3" image search only for pairs up to image_radius (galaxy scale, theta_E <= 1.5"); wider pairs
# are dropped (ASSUMPTION)
SEP_MAX = 3.0
# the LS image pair must be the catalogued pair (ASSUMPTION: separations agree within 0.5")
SEP_MATCH = 0.5
DEDUP_ARCSEC = 3.0  # same system in two tables (ASSUMPTION; SQLS positions are one image)


MANIFEST = "w12_niq_inputs.json"


def pinned() -> dict[str, str]:
    """Pinned hashes from the tracked manifest: VizieR tables by ID, plus ``"bricks"``."""
    man = paths.manifests_dir() / MANIFEST
    if not man.exists():
        return {}
    m = json.loads(man.read_text())
    pins = {r["source"]: r["sha256"] for r in m["inputs"]}
    pins["bricks"] = m["bricks"]["sha256"]
    if "tractor" in m:
        pins["tractor"] = m["tractor"]["sha256"]
    return pins


def check_pin(key: str, sha: str, pins: dict[str, str] | None) -> None:
    """Refuse a mismatch or a missing pin; ``pins`` is None under ``--repin``."""
    if pins is None:
        return
    if key not in pins:
        raise RuntimeError(f"{key}: no pinned sha256 in {MANIFEST}; inspect it, then --repin")
    if pins[key] != sha:
        raise RuntimeError(
            f"{key}: sha256 {sha} != pinned {pins[key]} (input changed or a stale cache); "
            "inspect, then rerun with --repin"
        )


def data_sha256(data: bytes) -> str:
    """sha256 of an ASU-TSV without its '#' header lines, which carry the request time."""
    body = b"\n".join(ln for ln in data.splitlines() if not ln.startswith(b"#"))
    return hashlib.sha256(body).hexdigest()


def check_vizier(data: bytes, source: str, path: Path) -> None:
    """Refuse an ASU-TSV that is an error page, empty or truncated at -out.max (the notices sit
    in '#' lines, which the hash ignores)."""
    lines = data.decode("utf-8", "replace").splitlines()
    notes = [
        ln for ln in lines if ln.startswith("#") and re.search(r"error|overflow|truncat", ln, re.I)
    ]
    n_data = len([ln for ln in lines if ln and not ln.startswith("#")]) - 3
    if notes or n_data < 1 or n_data >= MAX_ROWS:
        raise RuntimeError(f"{source}: bad VizieR response in {path} ({n_data} rows; {notes[:2]})")


def fetch(source: str, out: Path, pins: dict[str, str] | None) -> tuple[Path, dict]:
    """Download (or reuse) one VizieR table; refuse it if it differs from the pin."""
    f = out / (source.replace("/", "_") + ".tsv")
    if not f.exists():
        r = requests.get(VIZIER.format(source), timeout=300)
        r.raise_for_status()
        f.write_bytes(r.content)
    data = f.read_bytes()
    check_vizier(data, source, f)
    sha = data_sha256(data)
    check_pin(source, sha, pins)
    retrieved = dt.datetime.fromtimestamp(f.stat().st_mtime, dt.UTC).isoformat(timespec="seconds")
    return f, {
        "source": source,
        "url": VIZIER.format(source),
        "file": f.name,
        "bytes": len(data),
        "sha256": sha,
        "sha256_of": "data lines (ASU-TSV without '#' header lines)",
        "retrieved_utc": retrieved,
        "pipeline_version": __version__,
    }


def read_tsv(path: Path) -> Table:
    """VizieR ASU-TSV: header, units and dashes rows, then data (all columns as strings)."""
    lines = [
        ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")
    ]
    hdr = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[3:]]
    rows = [r + [""] * (len(hdr) - len(r)) for r in rows]
    return Table(rows=[[c.strip() for c in r] for r in rows], names=hdr, dtype=[str] * len(hdr))


def _float(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def sdss_name_radec(name: str) -> tuple[float, float]:
    """RA, Dec (deg) from an SDSS name ``Jhhmmss.ss+ddmmss.s``."""
    m = re.match(r"J?(\d{2})(\d{2})(\d{2}\.\d+)([+-])(\d{2})(\d{2})(\d{2}(?:\.\d+)?)", name)
    if not m:
        return np.nan, np.nan
    hh, mm, ss, sg, dd, dm, ds = m.groups()
    ra = 15 * (int(hh) + int(mm) / 60 + float(ss) / 3600)
    dec = int(dd) + int(dm) / 60 + float(ds) / 3600
    return ra, (-dec if sg == "-" else dec)


def comment_z_pair(comment: str) -> tuple[float, float]:
    """The two redshifts an SQLS comment quotes, e.g. "QSO pair (z=1.686, 1.600)"; NaN if none."""
    m = re.search(r"z\s*=\s*([0-9.]+)\s*,\s*([0-9.]+)", comment)
    return (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)


def different_redshift(z1, z2) -> np.ndarray:
    """Images of one source share its redshift: |dz| / (1 + z) > DZ_TOL means two sources."""
    z1, z2 = np.asarray(z1, float), np.asarray(z2, float)
    with np.errstate(invalid="ignore"):
        return np.abs(z1 - z2) / (1 + z1) > DZ_TOL


def sqls_group(comment: str) -> str:
    """rejected: no lens object, or a quasar pair / binary (as Lemon's "QSO pair" class);
    control: a catalogued lens."""
    c = comment.lower()
    if "sdss lens" in c or "known lens" in c:
        return "control"  # checked first: a comment naming a lens never makes a rejection
    if "no lens" in c or "qso pair" in c or "binary" in c:
        return "rejected"
    return ""


def build_sample(tables: dict[str, Table]) -> Table:
    """Rejected and control pairs from all inputs, merged within DEDUP_ARCSEC, then restricted to
    known catalogued separations <= SEP_MAX. ``meta["dropped"]``: rows dropped per reason."""
    rows, dropped = [], {}
    for source, t in tables.items():
        label = INPUTS[source]
        for r in t:
            if label == "Lemon2023":
                cls = r["Class"]
                group = (
                    "rejected"
                    if cls in LEMON_REJECTED
                    else "control"
                    if cls in LEMON_CONTROL
                    else ""
                )
                ra, dec, name, z, sep = (
                    _float(r["RAJ2000"]),
                    _float(r["DEJ2000"]),
                    r["Name"],
                    _float(r["z"]),
                    _float(r["Sep"]),
                )
                # z2 is a second *quasar* redshift only when n_z2 says "zqso=" (otherwise it is
                # z_lens or zgal); flagged (uncertain) values are not used
                z2 = (
                    _float(r["z2"])
                    if r.get("n_z2", "").startswith("zqso")
                    and not r.get("f_z2", "").strip()
                    and not r.get("f_z", "").strip()
                    else np.nan
                )
                comment = cls
            else:
                name = r["SDSS"]
                comment = r["Com"]
                group = sqls_group(comment)
                ra, dec = (
                    (_float(r["_RA"]), _float(r["_DE"]))
                    if "_RA" in t.colnames
                    else (np.nan, np.nan)
                )
                if not np.isfinite(ra):
                    ra, dec = sdss_name_radec(name)
                z, sep = _float(r["z"]), _float(r["theta"])
                z2 = comment_z_pair(comment)[1]
            if group and not np.isfinite(ra):
                dropped[label] = dropped.get(label, 0) + 1
            elif group:
                rows.append((name, ra, dec, z, z2, sep, group, label, comment))
    s = Table(
        rows=rows,
        names=("name", "ra", "dec", "z_source", "z2", "sep_cat", "group", "catalogue", "comment"),
    )
    s["comment"] = s["comment"].astype(object)  # notes are appended below (no truncation)
    # merge first (a lens row that is dropped below still promotes its group to control), then
    # drop component rows of one cluster-scale lens, unknown separations and wide pairs
    n_rows = len(s)
    s = dedup(s)
    n_dup = n_rows - len(s)
    comp = np.array(["component" in str(c).lower() for c in s["comment"]], bool)
    nosep = ~np.isfinite(np.asarray(s["sep_cat"], float))
    wide = np.asarray(s["sep_cat"], float) > SEP_MAX
    s = s[~comp & ~nosep & ~wide]
    s.meta["dropped"] = {
        "no coordinates": dropped,
        "duplicates": int(n_dup),
        "cluster-lens component rows": int(comp.sum()),
        "no catalogued separation": int((~comp & nosep).sum()),
        f"sep_cat > {SEP_MAX} arcsec": int((~comp & ~nosep & wide).sum()),
    }
    s["selection"] = "quasar"
    s["z_lens"] = np.nan
    s["theta_e"] = s["sep_cat"] / 2  # SIS model_prediction from the catalogued separation
    s["name_offset"] = 0.0
    s["pos_quantum"] = [
        lenscats.position_quantum_arcsec(a, d) for a, d in zip(s["ra"], s["dec"], strict=True)
    ]
    s["catalogues"] = s["catalogue"]
    s.meta.update(
        provenance=schema.Provenance.OBSERVED.value, source="; ".join(f"VizieR {k}" for k in tables)
    )
    return s


def dedup(s: Table) -> Table:
    """Merge entries within DEDUP_ARCSEC transitively (connected components). The kept row is
    Lemon's, else the newest SQLS release's; a component holding any catalogued lens is a control
    (noted in ``comment``), so a known lens never enters the rejected sample."""
    if len(s) < 2:
        return s
    order = {"Lemon2023": 0, "SQLS-DR7": 1, "SQLS-DR5": 2, "SQLS-DR3": 3}
    s = s[np.argsort([order[c] for c in s["catalogue"]], kind="stable")]
    c = SkyCoord(s["ra"], s["dec"], unit="deg")
    i, j, _, _ = c.search_around_sky(c, DEDUP_ARCSEC * u.arcsec)
    graph = csr_matrix((np.ones(len(i)), (i, j)), shape=(len(s), len(s)))
    _, label = connected_components(graph, directed=False)
    keep = []
    for lab in np.unique(label):
        members = np.flatnonzero(label == lab)  # sorted: members[0] has the preferred catalogue
        k = members[0]
        if s["group"][k] == "rejected" and np.any(s["group"][members] == "control"):
            s["group"][k] = "control"
            s["comment"][k] += " | listed as lens elsewhere"
        for m in members[1:]:  # keep the merged rows' vetting information
            s["comment"][k] += f" | {s['catalogue'][m]}: {s['comment'][m]}"
            if not np.isfinite(s["z2"][k]) and np.isfinite(s["z2"][m]):
                s["z2"][k] = s["z2"][m]
            if not np.isfinite(s["sep_cat"][k]) and np.isfinite(s["sep_cat"][m]):
                s["sep_cat"][k] = s["sep_cat"][m]
        keep.append(k)
    return s[np.sort(keep)]


def pair_colour_difference(images: Table, src: Table) -> np.ndarray:
    """g - z of image 1 minus image 2 (LS Tractor); NaN when either is missing. ``derived``."""
    if not len(src) or "mag_g" not in src.colnames:
        return np.full(len(images), np.nan)
    gz = np.asarray(src["mag_g"], float) - np.asarray(src["mag_z"], float)
    out = np.full(len(images), np.nan)
    for k, (i1, i2) in enumerate(zip(images["img1"], images["img2"], strict=True)):
        if i1 >= 0 and i2 >= 0:
            out[k] = gz[i1] - gz[i2]
    return out


def sheet(t: Table, out: Path, path: Path, title: str) -> None:
    """Contact sheet (w12_lenscats layout), stored as a JPEG to stay well under 1 MB."""
    from PIL import Image

    cut = out / "cutouts"
    cut.mkdir(exist_ok=True)
    t = t.copy()
    t["system_id"] = t["name"]
    t["label"] = [
        f"{c[:18]} {s} dgz={g:.2f}"
        for c, s, g in zip(t["comment"], t["test_status"], t["dgz"], strict=True)
    ]
    for r in t:
        w12.fetch_cutout(float(r["ra"]), float(r["dec"]), cut / f"{r['name']}.jpg")
    png = out / (path.stem + ".png")
    w12.contact_sheet(t, cut, png, title)
    Image.open(png).convert("RGB").save(path, quality=80)


def _coords(ra_col, dec_col) -> tuple[np.ndarray, np.ndarray]:
    """Degrees from decimal or sexagesimal ("h:m:s" / "d:m:s", as VizieR writes RA1/DE1) strings."""
    ra = np.array([_float(x) for x in ra_col])
    de = np.array([_float(x) for x in dec_col])
    rs = np.array([str(x).strip() for x in ra_col])
    ds = np.array([str(x).strip() for x in dec_col])
    todo = ~(np.isfinite(ra) & np.isfinite(de)) & (rs != "") & (ds != "")
    if todo.any():
        try:
            c = SkyCoord(rs[todo], ds[todo], unit=("hourangle", "deg"))
            ra[todo], de[todo] = c.ra.deg, c.dec.deg
        except ValueError:
            pass  # unparseable: stays NaN (binary_match raises if nothing parses)
    return ra, de


def binary_match(s: Table, binq: Table, radius: float = 3.0) -> np.ndarray:
    """Hennawi et al. 2006 binary-quasar entry within ``radius`` of the position (either quasar).

    Raises if the catalogue has rows but no parseable coordinate pair, so a format change cannot
    silently give "no binary"."""
    hit = np.zeros(len(s), bool)
    if not len(binq) or not len(s):
        return hit
    c = SkyCoord(s["ra"], s["dec"], unit="deg")
    pairs = [(rc, rc.replace("RA", "DE")) for rc in binq.colnames if rc.startswith(("_RA", "RA"))]
    parsed = 0
    for rc, dc in pairs:
        if dc not in binq.colnames:
            continue
        ra, de = _coords(binq[rc], binq[dc])
        ok = np.isfinite(ra) & np.isfinite(de)
        if not ok.any():
            continue
        parsed += 1
        _, d, _ = c.match_to_catalog_sky(SkyCoord(ra[ok], de[ok], unit="deg"))
        hit |= d.arcsec < radius
    if not parsed:
        raise ValueError(f"no coordinates parsed from {binq.colnames}")
    return hit


def cmd_screen(args) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    p = w12.Params()
    tables, manifest = {}, []
    pins = None if args.repin else pinned()
    for src in INPUTS:
        f, m = fetch(src, out, pins)
        tables[src] = read_tsv(f)
        m["rows"] = len(tables[src])
        manifest.append(m)
    fb, mb = fetch(HENNAWI, out, pins)
    binq = read_tsv(fb)
    mb["rows"] = len(binq)
    manifest.append(mb)
    s = build_sample(tables)
    bricks, bmeta = w12.load_bricks(out)
    check_pin("bricks", bmeta["sha256"], pins)
    cov = lenscats.brick_coverage(s, bricks)
    for c in cov.colnames:
        s[c] = cov[c]
    covered = np.asarray(s["covered"], bool)
    sc = s[covered]
    src = w12.query_tractor(sc, p, out / "tractor_cache")
    tractor_sha = hashlib.sha256(
        "\n".join(
            f"{i},{r:.7f},{d:.7f},{t},{z:.4f}"
            for i, r, d, t, z in sorted(
                zip(src["ls_id"], src["ra"], src["dec"], src["type"], src["mag_z"], strict=True)
            )
        ).encode()
    ).hexdigest()
    check_pin("tractor", tractor_sha, pins)
    req = w12.required_mags(sc, sc["theta_e"], D056_FJ, p)
    for k, v in req.items():
        sc[k] = v
    images = lenscats.pair_images(sc, src, p.image_radius)
    sc["sep_ls"] = images["sep"]
    sc["dgz"] = pair_colour_difference(images, src)
    sc["colour_match"] = np.abs(sc["dgz"]) <= COLOUR_TOL
    sc["colour_mismatch"] = np.abs(sc["dgz"]) > COLOUR_TOL  # NaN (no colour) is neither
    sc["pair_match"] = np.abs(sc["sep_ls"] - sc["sep_cat"]) <= SEP_MATCH
    sc["test_status"] = w12.deflector_test(sc, src, p, images)["test_status"]
    sc["undecided_flags"] = w12.flags_undecided(sc, src, p)
    sc["detectable_typical"] = sc["req_mag_z_typical"] < sc["depth_z"] - p.margin
    sc["detectable_conservative"] = sc["req_mag_z"] < sc["depth_z"] - p.margin
    sc["hennawi_binary"] = binary_match(sc, binq)
    sc["different_z"] = different_redshift(sc["z_source"], sc["z2"])
    decided = np.isin(sc["test_status"], ["deflector", "none"]) & (sc["undecided_flags"] == "")
    decided &= np.asarray(sc["pair_match"], bool)  # otherwise the test ran on another pair
    sc["decided"] = decided
    summary = {
        "inputs": manifest,
        "bricks": bmeta,
        "fj_calibration": D056_FJ.__dict__,
        "params": p.__dict__,
        "n_systems": len(s),
        "dropped": s.meta["dropped"],
        "n_covered": int(covered.sum()),
    }
    for g in ("rejected", "control"):
        m = sc["group"] == g
        st = {
            k: int(v)
            for k, v in zip(*np.unique(sc["test_status"][m], return_counts=True), strict=True)
        }
        d = m & decided
        nd = int((d & (sc["test_status"] == "deflector")).sum())
        nn = int((d & (sc["test_status"] == "none")).sum())
        summary[g] = {
            "n_all": int((s["group"] == g).sum()),
            "n_covered": int(m.sum()),
            "status": st,
            "decided": int(d.sum()),
            "deflector": nd,
            "none": nn,
            "none_colour_match": int(
                (d & (sc["test_status"] == "none") & sc["colour_match"]).sum()
            ),
            "none_colour_mismatch": int(
                (d & (sc["test_status"] == "none") & sc["colour_mismatch"]).sum()
            ),
            "pair_mismatch": int(
                (m & np.isin(sc["test_status"], ["deflector", "none"]) & ~sc["pair_match"]).sum()
            ),
            "none_typical_detectable": int(
                (d & (sc["test_status"] == "none") & sc["detectable_typical"]).sum()
            ),
        }
    ctl = summary["control"]
    summary["control_efficiency"] = ctl["deflector"] / ctl["decided"] if ctl["decided"] else None
    sc.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="scripts/w12_niq.py screen (D-056 test on rejected pairs)",
    )
    res = paths.repo_root() / "results" / "w12_niq"
    res.mkdir(parents=True, exist_ok=True)
    keep = [
        "name",
        "ra",
        "dec",
        "z_source",
        "sep_cat",
        "sep_ls",
        "group",
        "catalogue",
        "comment",
        "depth_z",
        "req_mag_z_typical",
        "req_mag_z",
        "test_status",
        "undecided_flags",
        "decided",
        "dgz",
        "colour_match",
        "colour_mismatch",
        "pair_match",
        "detectable_typical",
        "detectable_conservative",
        "hennawi_binary",
        "z2",
        "different_z",
    ]
    sc[keep].write(res / "systems.ecsv", overwrite=True)
    (res / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    if args.repin:  # the manifest is the pin: rewritten only on purpose
        man = paths.manifests_dir() / MANIFEST
        tmeta = {
            "sha256": tractor_sha,
            "rows": len(src),
            "sha256_of": "sorted ls_id,ra,dec,type,mag_z",
            "source": "Legacy Surveys DR10 Tractor boxes (Data Lab TAP)",
            "box_arcsec": p.box,
        }
        man.write_text(
            json.dumps({"inputs": manifest, "bricks": bmeta, "tractor": tmeta}, indent=2) + "\n"
        )
    print(
        json.dumps(
            {
                k: v
                for k, v in summary.items()
                if k in ("rejected", "control", "control_efficiency", "n_systems", "n_covered")
            },
            indent=2,
        )
    )
    if args.sheet:
        for g in ("rejected", "control"):
            sel = sc[(sc["group"] == g) & sc["decided"]]
            if len(sel):
                sheet(sel, out, res / f"contact_sheet_{g}.jpg", f"Decided {g} pairs (LS DR10 grz)")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("screen")
    s.add_argument("--out", default=str(paths.outputs_dir() / "w12_niq"))
    s.add_argument("--sheet", action="store_true", help="cutouts + contact sheets of decided pairs")
    s.add_argument(
        "--repin", action="store_true", help="accept new input hashes (after inspection)"
    )
    args = ap.parse_args(argv)
    {"screen": cmd_screen}[args.cmd](args)


if __name__ == "__main__":
    main()
