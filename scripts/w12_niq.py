"""W1/W2 in rejected lensed-quasar pairs: the D-056 dark-deflector test on pairs that lens searches
rejected because no lens galaxy was seen (D-064).

Inputs (VizieR, pinned by sha256 in ``data/manifests/w12_niq_inputs.json``):
- Lemon et al. 2023 (MNRAS 520, 3305) Gaia-selected candidates: classes "UQP" (unclassified
  quasar pair: same redshift, no lens seen; "akin to NIQs") and "QSO pair" are the *rejected*
  sample; "lens" and "quad" are the *control* sample of real lenses; classes with "?" are
  undecided; other classes (QSO + star, projected, ...) veto a rejection of the same system.
- SQLS (Inada et al. 2008, 2010, 2012) candidate tables: comments "no lens(ing) object", "QSO
  pair" and "binary" are rejected; "SDSS lens" / "known lens" are controls; QSO+star, different
  SED and similar companion rows are dropped without vetoing (they describe another companion).

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
from dataclasses import asdict, dataclass
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
LEMON_COLUMNS = {"Name", "RAJ2000", "DEJ2000", "z", "f_z", "Sep", "Class", "z2", "n_z2", "f_z2"}
SQLS_COLUMNS = {"SDSS", "z", "theta", "Com"}
HENNAWI = "J/AJ/131/1/binqso"  # binary-quasar catalogue for vetting (Hennawi et al. 2006)
LEMON_REJECTED = {"UQP", "QSO pair"}  # "UQP (?)" is undecided, like "lens (?)"
LEMON_CONTROL = {"lens", "quad"}  # "lens (?)" is undecided: neither control nor veto
# lenses the quasar pair test cannot decide (extended images): they promote a merged rejection to
# control, then leave the sample
LEMON_OTHER_LENS = {"lensed gal."}
# SQLS comments that classify a candidate as something other than a quasar pair or a lens
# ("Sngle QSO" is a typo in J/AJ/140/403/table2)
SQLS_NONPAIR = re.compile(
    r"qso\s*\+\s*(star|galaxy|unknown)|different sed|not qso|s(i)?ngle qso", re.I
)
# D-056 calibration (LS z; 605 galaxy-selected lenses; docs/exotic_limits.md), reused unchanged
D056_FJ = lenscats.FJCalibration(a=20.41, k=0.67, slope=-10.0, rms=0.88, n=605, band="ls_z")


@dataclass(frozen=True)
class NiqParams:
    """Sample and vetting thresholds (ASSUMPTIONs; recorded in summary.json)."""

    # lensed images share a colour: |delta(g - z)| <= colour_tol allows microlensing, dust and
    # variability (the tolerance of quasar_pair_test's further-image rule)
    colour_tol: float = 0.5
    # two redshifts differing by more than dz_tol (1 + z) (~3,000 km/s, beyond broad-line
    # redshift errors) belong to two quasars, not two images
    dz_tol: float = 0.01
    # catalogue positions are one image (SQLS: the SDSS quasar), so the second image lies inside
    # the 3" image search only for pairs up to image_radius; wider pairs are dropped
    sep_max: float = w12.Params().image_radius
    # the LS image pair must be the catalogued pair: separations agree within sep_match
    sep_match: float = 0.5
    # entries closer than this in two tables are one system (merged transitively)
    dedup_arcsec: float = w12.Params().merge_radius
    # Hennawi et al. 2006 binary match radius around the catalogued position
    binary_arcsec: float = 3.0
    # an SQLS companion row belongs to the primary row above it when their distance equals theta
    # within this tolerance
    companion_tol: float = 1.0
    # a Hennawi binary is the catalogued pair when its theta agrees with the catalogued separation
    # within this tolerance
    binary_theta_tol: float = 1.0


P = NiqParams()
MANIFEST = "w12_niq_inputs.json"
# all Tractor columns the D-056 chain and the vetting read (w12_lenscats.TRACTOR_COLS minus none)
TRACTOR_PIN_COLUMNS = tuple(w12.TRACTOR_COLS.split(","))


def pinned() -> dict[str, str]:
    """Pinned hashes from the tracked manifest: VizieR tables by ID, plus ``"bricks"``."""
    man = paths.manifests_dir() / MANIFEST
    if not man.exists():
        return {}
    m = json.loads(man.read_text())
    pins = {r["source"]: r["sha256"] for r in m["inputs"]}
    for key in ("bricks", "tractor"):
        if key in m:
            pins[key] = m[key]["sha256"]
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
    # VizieR status lines only ('#INFO QUERY_STATUS=...', '#++++' notices), never descriptions
    notes = [
        ln
        for ln in lines
        if re.match(r"#INFO\s+QUERY_STATUS\s*=\s*(ERROR|OVERFLOW)", ln, re.I)
        or (ln.startswith("#++") and re.search(r"error|truncat", ln, re.I))
    ]
    n_data = len([ln for ln in lines if ln and not ln.startswith("#")]) - 3
    if notes or n_data < 1 or n_data >= MAX_ROWS:
        raise RuntimeError(f"{source}: bad VizieR response in {path} ({n_data} rows; {notes[:2]})")


def fetch(source: str, out: Path, pins: dict[str, str] | None) -> tuple[Path, dict]:
    """Download (or reuse) one VizieR table; refuse it if it differs from the pin."""
    f = out / (source.replace("/", "_") + ".tsv")
    if pins is None and f.exists():
        f.unlink()  # --repin pins a fresh download, never a stale cache
    if not f.exists():
        r = requests.get(VIZIER.format(source), timeout=300)
        r.raise_for_status()
        check_vizier(r.content, source, f)  # before caching: a bad response is never stored
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
    """VizieR ASU-TSV, one resource: header, units and dashes rows, then data (strings). Raises on
    a second resource or a row wider than the header (trailing empty fields are dropped)."""
    lines = [
        ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")
    ]
    hdr = lines[0].split("\t")
    rows = []
    for ln in lines[3:]:
        r = ln.split("\t")
        while len(r) > len(hdr) and not r[-1].strip():
            r.pop()
        if len(r) > len(hdr) or r[: len(hdr)] == hdr or set(ln.replace("\t", "")) == {"-"}:
            raise ValueError(f"{path}: unexpected row (second resource or wide row): {ln[:80]}")
        rows.append(r + [""] * (len(hdr) - len(r)))
    return Table(rows=[[c.strip() for c in r] for r in rows], names=hdr, dtype=[str] * len(hdr))


def _float(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def sdss_name_radec(name: str) -> tuple[float, float]:
    """RA, Dec (deg) from an SDSS name ``Jhhmmss.ss+ddmmss.s``."""
    m = re.match(r"J?(\d{2})(\d{2})(\d{2}(?:\.\d+)?)([+-])(\d{2})(\d{2})(\d{2}(?:\.\d+)?)?", name)
    if not m:
        return np.nan, np.nan
    hh, mm, ss, sg, dd, dm, ds = m.groups()
    ra = 15 * (int(hh) + int(mm) / 60 + float(ss) / 3600)
    dec = int(dd) + int(dm) / 60 + (float(ds) if ds else 0.0) / 3600
    return ra, (-dec if sg == "-" else dec)


def comment_z_pair(comment: str) -> tuple[float, float]:
    """The two redshifts an SQLS comment quotes, e.g. "QSO pair (z=1.686, 1.600)"; NaN if none."""
    m = re.search(r"z\s*=\s*(\d+(?:\.\d+)?)\s*,\s*(\d+(?:\.\d+)?)", comment)
    return (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)


def different_redshift(z1, z2) -> np.ndarray:
    """Images of one source share its redshift: |dz| / (1 + z) > P.dz_tol means two sources."""
    z1, z2 = np.asarray(z1, float), np.asarray(z2, float)
    with np.errstate(invalid="ignore"):
        return np.abs(z1 - z2) / (1 + z1) > P.dz_tol


def sqls_radec(t: Table, r) -> tuple[float, float]:
    ra, dec = (_float(r["_RA"]), _float(r["_DE"])) if "_RA" in t.colnames else (np.nan, np.nan)
    return (ra, dec) if np.isfinite(ra) else sdss_name_radec(r["SDSS"])


def _sqls_z(r) -> float:
    """An SQLS row's redshift, NaN when blank or flagged (f_z, e.g. "(" = uncertain)."""
    flag = r["f_z"].strip() if "f_z" in r.colnames else ""
    return np.nan if flag else _float(r["z"])


def primary_rows(t: Table) -> np.ndarray:
    """Index of the nearest earlier row without theta (the primary quasar) for each row; -1 if
    none. Earlier companions of the same primary are skipped."""
    has_theta = np.array([np.isfinite(_float(x)) for x in t["theta"]])
    idx = np.where(~has_theta, np.arange(len(t)), -1)
    last = np.maximum.accumulate(idx) if len(t) else idx
    return np.concatenate([[-1], last[:-1]]) if len(t) else last


def sqls_redshifts(t: Table, r, prim: np.ndarray | None = None) -> tuple[float, float]:
    """(quasar z, second z) for an SQLS candidate row.

    Pair-format tables (DR3 table3, DR5 table3, DR7 table4) list the primary quasar (its z, no
    theta) on the row before the companion (theta, Com and the companion's own z, if measured).
    The primary's z is the source redshift and the companion's z the second one; two redshifts
    quoted in the comment win over both."""
    za, zb = comment_z_pair(r["Com"])
    if np.isfinite(za):
        return za, zb
    theta = _float(r["theta"])
    k = (primary_rows(t) if prim is None else prim)[r.index]
    if k >= 0 and not t[k]["Com"].strip() and np.isfinite(theta):
        (ra1, de1), (ra2, de2) = sqls_radec(t, t[k]), sqls_radec(t, r)
        d = SkyCoord(ra1, de1, unit="deg").separation(SkyCoord(ra2, de2, unit="deg")).arcsec
        if abs(d - theta) < P.companion_tol:
            return _sqls_z(t[k]), _sqls_z(r)
    return _sqls_z(r), np.nan


def lemon_second_qso_z(r) -> float:
    """Lemon et al. 2023 ``z2``: a second quasar redshift when ``n_z2`` is blank (or "zqso="),
    otherwise a lens or galaxy redshift ("z_lens=", "zgal="); flagged values are not used."""
    n = r["n_z2"].strip()  # columns are required (build_sample checks LEMON_COLUMNS)
    if (n == "" or n.startswith("zqso")) and not (r["f_z"].strip() or r["f_z2"].strip()):
        return _float(r["z2"])
    return np.nan


def sqls_group(comment: str) -> str:
    """rejected: no lens object, or a quasar pair / binary (as Lemon's "QSO pair" class);
    control: a catalogued lens."""
    c = comment.lower()
    if re.search(r"(?<!not a )(?<!not )(?<!no )\b(sdss|known) lens(ed)?\b", c):
        return "control"  # first: a comment naming a lens never makes a rejection
    if SQLS_NONPAIR.search(c):
        return "nonpair"  # second: a star/galaxy/SED classification vetoes a rejection
    if re.search(r"\bno lens|\bqso pair\b|(?<!not )(?<!not a )\bbinary\b", c):
        return "rejected"
    return ""


def build_sample(tables: dict[str, Table]) -> Table:
    """Rejected and control pairs from all inputs, merged within P.dedup_arcsec, then restricted to
    known catalogued separations <= P.sep_max. ``meta["dropped"]``: rows dropped per reason."""
    rows, dropped = [], {}
    for source, t in tables.items():
        label = INPUTS[source]
        prim = primary_rows(t) if label != "Lemon2023" else None
        need = {"Lemon2023": LEMON_COLUMNS}.get(label, SQLS_COLUMNS) - set(t.colnames)
        if need:
            raise ValueError(f"{source}: missing columns {sorted(need)}")
        for r in t:
            if label == "Lemon2023":
                cls = r["Class"]
                group = (
                    "rejected"
                    if cls in LEMON_REJECTED
                    else "control"
                    if cls in LEMON_CONTROL | LEMON_OTHER_LENS
                    else ""
                    if not cls or "?" in cls  # undecided classes neither count nor veto
                    else "nonpair"
                )
                ra, dec, name, z, sep = (
                    _float(r["RAJ2000"]),
                    _float(r["DEJ2000"]),
                    r["Name"],
                    np.nan if r["f_z"].strip() else _float(r["z"]),  # flagged z not used
                    _float(r["Sep"]),
                )
                z2 = lemon_second_qso_z(r)
                comment = cls
            else:
                name = r["SDSS"]
                comment = r["Com"]
                group = sqls_group(comment)
                ra, dec = sqls_radec(t, r)
                sep = _float(r["theta"])
                z, z2 = sqls_redshifts(t, r, prim)
            if group and not np.isfinite(ra):
                dropped[label] = dropped.get(label, 0) + 1
            elif group:
                # rows of one cluster-scale lens, and lensed galaxies: lenses outside the test
                comp = "component" in comment.lower() or comment in LEMON_OTHER_LENS
                rows.append((name, ra, dec, z, z2, sep, group, label, comment, comp))
    s = Table(
        rows=rows,
        names=(
            "name",
            "ra",
            "dec",
            "z_source",
            "z2",
            "sep_cat",
            "group",
            "catalogue",
            "comment",
            "component",
        ),
    )
    s["comment"] = s["comment"].astype(object)  # notes are appended below (no truncation)
    s["different_z"] = different_redshift(s["z_source"], s["z2"])  # per row, OR-ed in dedup
    # merge first (a lens row that is dropped below still promotes its group to control), then
    # drop component rows of one cluster-scale lens, unknown separations and wide pairs
    n_rows = len(s)
    s = dedup(s)
    n_dup = n_rows - len(s)
    nonpair = np.asarray(s["group"]) == "nonpair"
    s = s[~nonpair]
    comp = np.asarray(s["component"], bool)  # the kept row's own flag
    nosep = ~np.isfinite(np.asarray(s["sep_cat"], float))
    wide = np.asarray(s["sep_cat"], float) > P.sep_max
    s = s[~comp & ~nosep & ~wide]
    s.meta["dropped"] = {
        "no coordinates": dropped,
        "duplicates": int(n_dup),
        "non-pair classifications (incl. vetoed rejections)": int(nonpair.sum()),
        "lens rows outside the pair test (cluster components, lensed galaxies)": int(comp.sum()),
        "no catalogued separation": int((~comp & nosep).sum()),
        f"sep_cat > {P.sep_max} arcsec": int((~comp & ~nosep & wide).sum()),
    }
    s["selection"] = "quasar"
    s["z_lens"] = np.nan
    s["theta_e"] = s["sep_cat"] / 2  # SIS model_prediction from the catalogued separation
    s["name_offset"] = [  # designation vs position: w12.flags_undecided flags > 5"
        lenscats.name_position_offset(n, a, d)
        for n, a, d in zip(s["name"], s["ra"], s["dec"], strict=True)
    ]
    s["pos_quantum"] = [
        lenscats.position_quantum_arcsec(a, d) for a, d in zip(s["ra"], s["dec"], strict=True)
    ]
    s["catalogues"] = s["catalogue"]
    s.meta.update(  # merged, reclassified and with a model theta_e: derived from observed rows
        provenance=schema.Provenance.DERIVED.value,
        source="build_sample of " + "; ".join(f"VizieR {k}" for k in tables),
    )
    return s


def dedup(s: Table) -> Table:
    """Merge entries within P.dedup_arcsec transitively (connected components).

    The merged system is a control if any member is a lens (a known lens never enters the
    rejected sample); else rejected if any member is, unless Lemon classified the system as a
    non-pair (star, galaxy, projected), which vetoes it; else a non-pair (dropped later). The
    kept row is the first member of the winning group in catalogue order (Lemon, then the newest
    SQLS release); its position, separation and redshifts are never mixed with other members'.
    Other members add their comments, and ``different_z`` is OR-ed over the system."""
    if len(s) < 2:
        return s
    order = {"Lemon2023": 0, "SQLS-DR7": 1, "SQLS-DR5": 2, "SQLS-DR3": 3}
    s = s[np.argsort([order[c] for c in s["catalogue"]], kind="stable")]
    c = SkyCoord(s["ra"], s["dec"], unit="deg")
    i, j, _, _ = c.search_around_sky(c, P.dedup_arcsec * u.arcsec)
    graph = csr_matrix((np.ones(len(i)), (i, j)), shape=(len(s), len(s)))
    _, label = connected_components(graph, directed=False)
    grp = np.array(s["group"], copy=True)  # original classes (s["group"] is relabelled below)
    keep = []
    for lab in np.unique(label):
        members = np.flatnonzero(label == lab)  # sorted by catalogue preference
        lens = (grp[members] == "control") | np.asarray(s["component"][members], bool)
        rej = grp[members] == "rejected"
        lemon_nonpair = (grp[members] == "nonpair") & (s["catalogue"][members] == "Lemon2023")
        sep = np.asarray(s["sep_cat"], float)[members]
        testable = ~np.asarray(s["component"][members], bool) & (sep <= P.sep_max)
        if lens.any():
            k = members[np.argmax(lens & testable) if (lens & testable).any() else np.argmax(lens)]
            new, note = "control", " | also listed as a rejection" if rej.any() else ""
        elif rej.any() and not lemon_nonpair.any():
            usable = rej & testable  # the member describing the testable pair
            k, new, note = members[np.argmax(usable if usable.any() else rej)], "rejected", ""
        else:
            k, new, note = members[0], "nonpair", ""
        s["group"][k] = new
        s["comment"][k] += note
        for m in members[members != k]:
            s["comment"][k] += f" | {s['catalogue'][m]}: {s['comment'][m]}"
        same = members[grp[members] == grp[k]]  # rows of the kept row's own class
        # two redshifts from rows describing the same pair (same class), never from non-pair
        # companions or other objects in the group
        s["different_z"][k] = bool(np.any(np.asarray(s["different_z"][same], bool)))
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
    t["system_id"] = [  # cutout cache key: name and position (a moved row gets a new image)
        f"{n}_{a:.5f}_{d:+.5f}" for n, a, d in zip(t["name"], t["ra"], t["dec"], strict=True)
    ]
    t["label"] = [
        f"{c[:18]} {s} dgz={g:.2f}"
        for c, s, g in zip(t["comment"], t["test_status"], t["dgz"], strict=True)
    ]
    for r in t:
        w12.fetch_cutout(float(r["ra"]), float(r["dec"]), cut / f"{r['system_id']}.jpg")
    png = out / (path.stem + ".png")
    w12.contact_sheet(t, cut, png, title)
    Image.open(png).convert("RGB").save(path, quality=80)


def _coords(ra_col, dec_col) -> tuple[np.ndarray, np.ndarray]:
    """Degrees from decimal or sexagesimal ("h:m:s" / "d:m:s", as VizieR writes RA1/DE1) strings.
    Raises on any non-empty entry that parses as neither (no silent "no match")."""
    ra = np.array([_float(x) for x in ra_col])
    de = np.array([_float(x) for x in dec_col])
    rs = np.array([str(x).strip() for x in ra_col])
    ds = np.array([str(x).strip() for x in dec_col])
    filled = (rs != "") & (ds != "")
    if np.any(np.isfinite(ra) & ~np.isfinite(de) & filled):
        raise ValueError("decimal RA with a sexagesimal Dec: ambiguous units")
    # sexagesimal RA, with a Dec in d:m:s or whole degrees ("+05" floats but is degrees anyway)
    todo = np.flatnonzero(~np.isfinite(ra) & filled)
    if len(todo):
        try:
            c = SkyCoord(rs[todo], ds[todo], unit=("hourangle", "deg"))
        except ValueError as e:
            raise ValueError(f"unparseable coordinates among {len(todo)} entries: {e}") from e
        ra[todo], de[todo] = c.ra.deg, c.dec.deg
    return ra, de


def binary_match(s: Table, binq: Table, radius: float | None = None) -> np.ndarray:
    """Hennawi et al. 2006 binary-quasar entry within ``radius`` of the position (either quasar).

    Raises if the catalogue has rows but no parseable coordinate pair, so a format change cannot
    silently give "no binary"."""
    radius = P.binary_arcsec if radius is None else radius
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
        idx, d, _ = c.match_to_catalog_sky(SkyCoord(ra[ok], de[ok], unit="deg"))
        near = d.arcsec < radius
        if "theta" in binq.colnames and "sep_cat" in s.colnames:
            # the binary must be the catalogued pair, not a wide binary next to it
            th = np.array([_float(x) for x in binq["theta"]])[ok][idx]
            sep = np.asarray(s["sep_cat"], float)
            agree = np.abs(th - sep) < P.binary_theta_tol
            near &= agree | ~(np.isfinite(th) & np.isfinite(sep))
        hit |= near
    if not parsed:
        raise ValueError(f"no coordinates parsed from {binq.colnames}")
    return hit


def tractor_times(out: Path) -> str:
    """Oldest and newest retrieval times (UTC, file mtimes) of the cached Tractor batches."""
    t = sorted(f.stat().st_mtime for f in (out / "tractor_cache").glob("tractor_*.ecsv"))
    if not t:
        return ""
    iso = [
        dt.datetime.fromtimestamp(x, dt.UTC).isoformat(timespec="seconds") for x in (t[0], t[-1])
    ]
    return iso[0] if iso[0] == iso[1] else f"{iso[0]} .. {iso[1]}"


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
    # row order and float format of a TAP CSV are not guaranteed: pin the sorted, formatted rows
    bmeta["sha256_raw"] = bmeta["sha256"]
    bmeta["sha256"] = hashlib.sha256(
        "\n".join(
            sorted(
                f"{r['brickname']},{float(r['ra1']):.6f},{float(r['ra2']):.6f},"
                f"{float(r['dec1']):.6f},{float(r['dec2']):.6f},{int(r['nexp_r'])},"
                f"{int(r['nexp_z'])},{float(r['galdepth_z']):.4f}"
                for r in bricks
            )
        ).encode()
    ).hexdigest()
    bmeta["sha256_of"] = "sorted rows: brickname,ra1,ra2,dec1,dec2,nexp_r,nexp_z,galdepth_z"
    bmeta.update(
        uri=f"{w12.TAP} (ADQL: {w12.BRICKS_QUERY})",
        file="bricks_dr10_south.csv",
        retrieved_utc=dt.datetime.fromtimestamp(
            (out / "bricks_dr10_south.csv").stat().st_mtime, dt.UTC
        ).isoformat(timespec="seconds"),
        pipeline_version=__version__,
    )
    check_pin("bricks", bmeta["sha256"], pins)
    cov = lenscats.brick_coverage(s, bricks)
    for c in cov.colnames:
        s[c] = cov[c]
    covered = np.asarray(s["covered"], bool)
    sc = s[covered]
    src = w12.query_tractor(sc, p, out / "tractor_cache")
    # pin only the rows within a box of the current sample (the cache may hold earlier batches)
    near = np.zeros(len(src), bool)
    if len(src) and len(sc) and "ra" in src.colnames:
        cs = SkyCoord(np.asarray(src["ra"], float), np.asarray(src["dec"], float), unit="deg")
        cp = SkyCoord(np.asarray(sc["ra"], float), np.asarray(sc["dec"], float), unit="deg")
        icp, ics, _, _ = cs.search_around_sky(cp, p.box * np.sqrt(2) * u.arcsec)
        dra, ddec = cp[icp].spherical_offsets_to(cs[ics])  # the query boxes, not a circle
        inbox = (np.abs(dra.arcsec) <= p.box) & (np.abs(ddec.arcsec) <= p.box)
        near[np.unique(ics[inbox])] = True
    src = src[near] if len(src) else src
    tractor_sha = hashlib.sha256(
        "\n".join(
            # every Tractor column the test reads (pair test, colour, required depth, mask flags)
            ",".join(f"{v:.7f}" if isinstance(v, float) else str(v) for v in row)
            for row in sorted(
                zip(*(np.asarray(src[c]).tolist() for c in TRACTOR_PIN_COLUMNS), strict=True)
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
    sc["colour_match"] = np.abs(sc["dgz"]) <= P.colour_tol
    sc["colour_mismatch"] = np.abs(sc["dgz"]) > P.colour_tol  # NaN (no colour) is neither
    sc["pair_match"] = np.abs(sc["sep_ls"] - sc["sep_cat"]) <= P.sep_match
    sc["test_status"] = w12.deflector_test(sc, src, p, images)["test_status"]
    sc["undecided_flags"] = w12.flags_undecided(sc, src, p)
    sc["detectable_typical"] = sc["req_mag_z_typical"] < sc["depth_z"] - p.margin
    sc["detectable_conservative"] = sc["req_mag_z"] < sc["depth_z"] - p.margin
    sc["hennawi_binary"] = binary_match(sc, binq)
    decided = np.isin(sc["test_status"], ["deflector", "none"]) & (sc["undecided_flags"] == "")
    decided &= np.asarray(sc["pair_match"], bool)  # otherwise the test ran on another pair
    sc["decided"] = decided
    summary = {
        "inputs": manifest,
        "bricks": bmeta,
        "fj_calibration": D056_FJ.__dict__,
        "params": p.__dict__,
        "niq_params": asdict(P),
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
        none = d & (sc["test_status"] == "none")
        cm = none & sc["colour_match"]
        summary[g] = {
            "n_all": int((s["group"] == g).sum()),
            "n_covered": int(m.sum()),
            "status": st,
            "decided": int(d.sum()),
            "deflector": int((d & (sc["test_status"] == "deflector")).sum()),
            "none": int(none.sum()),
            "none_colour_match": int(cm.sum()),
            "none_colour_mismatch": int((none & sc["colour_mismatch"]).sum()),
            "none_no_colour": int((none & ~sc["colour_match"] & ~sc["colour_mismatch"]).sum()),
            "pair_mismatch": int(
                (m & np.isin(sc["test_status"], ["deflector", "none"]) & ~sc["pair_match"]).sum()
            ),
            "none_typical_detectable": int((none & sc["detectable_typical"]).sum()),
            "none_colour_match_binary": int((cm & sc["hennawi_binary"]).sum()),
            "none_colour_match_different_z": int((cm & sc["different_z"]).sum()),
            "none_untestable": int((cm & ~sc["hennawi_binary"] & ~sc["different_z"]).sum()),
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
            "sha256_of": "sorted rows of " + ",".join(TRACTOR_PIN_COLUMNS),
            "source": "Legacy Surveys DR10 Tractor boxes (Data Lab TAP)",
            "uri": f"{w12.TAP} (ls_dr10.tractor, {w12.TRACTOR_COLS})",
            "file": "tractor_cache/tractor_*.ecsv (rows near the sample)",
            "retrieved_utc": tractor_times(out),
            "pipeline_version": __version__,
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
