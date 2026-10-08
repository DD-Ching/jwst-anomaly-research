"""W1/W2 in rejected lensed-quasar pairs: the D-056 dark-deflector test on pairs that lens searches
rejected because no lens galaxy was seen (D-TBD).

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
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w12_lenscats as w12  # noqa: E402

from jwst_anomaly import lenscats, paths, schema  # noqa: E402

VIZIER = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source={}&-out.max=5000&-out.all"
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
DEDUP_ARCSEC = 3.0  # same system in two tables (ASSUMPTION; SQLS positions are one image)


def fetch(source: str, out: Path) -> tuple[Path, dict]:
    f = out / (source.replace("/", "_") + ".tsv")
    if not f.exists():
        r = requests.get(VIZIER.format(source), timeout=300)
        r.raise_for_status()
        f.write_bytes(r.content)
    data = f.read_bytes()
    return f, {
        "source": source,
        "url": VIZIER.format(source),
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
    }


def read_tsv(path: Path) -> Table:
    """VizieR ASU-TSV: header, units and dashes rows, then data (all columns as strings)."""
    lines = [ln for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]
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


def sqls_group(comment: str) -> str:
    c = comment.lower()
    if "no lens" in c:
        return "rejected"
    if "sdss lens" in c or "known lens" in c:
        return "control"
    return ""


def build_sample(tables: dict[str, Table]) -> Table:
    rows = []
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
            if group and np.isfinite(ra):
                rows.append((name, ra, dec, z, sep, group, label, comment))
    s = Table(
        rows=rows,
        names=("name", "ra", "dec", "z_source", "sep_cat", "group", "catalogue", "comment"),
    )
    s["comment"] = s["comment"].astype(object)  # notes are appended below (no truncation)
    # deduplicate across tables/releases within DEDUP_ARCSEC: keep Lemon, then the newest SQLS
    order = {"Lemon2023": 0, "SQLS-DR7": 1, "SQLS-DR5": 2, "SQLS-DR3": 3}
    s = s[np.argsort([order[c] for c in s["catalogue"]], kind="stable")]
    c = SkyCoord(s["ra"], s["dec"], unit="deg")
    keep = np.ones(len(s), bool)
    for i in range(len(s)):
        if keep[i]:
            dup = (c[i].separation(c).arcsec < DEDUP_ARCSEC) & (np.arange(len(s)) > i)
            # a pair listed as rejected in one table and as a lens in another stays a control
            if np.any(dup & (s["group"] == "control")) and s["group"][i] == "rejected":
                s["group"][i] = "control"
                s["comment"][i] += " | listed as lens elsewhere"
            keep &= ~dup
    s = s[keep]
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


def pair_colour_difference(images: Table, src: Table) -> np.ndarray:
    """g - z of image 1 minus image 2 (LS Tractor); NaN when either is missing. ``derived``."""
    gz = lenscats._gz(src) if len(src) else np.zeros(0)
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
    for i, (r, d) in enumerate(zip(ra_col, dec_col, strict=True)):
        if not (np.isfinite(ra[i]) and np.isfinite(de[i])) and str(r).strip() and str(d).strip():
            try:
                c = SkyCoord(str(r).strip(), str(d).strip(), unit=("hourangle", "deg"))
                ra[i], de[i] = c.ra.deg, c.dec.deg
            except ValueError:
                pass
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
    for src in INPUTS:
        f, m = fetch(src, out)
        tables[src] = read_tsv(f)
        m["rows"] = len(tables[src])
        manifest.append(m)
    fb, mb = fetch(HENNAWI, out)
    binq = read_tsv(fb)
    mb["rows"] = len(binq)
    manifest.append(mb)
    s = build_sample(tables)
    bricks, bmeta = w12.load_bricks(out)
    cov = lenscats.brick_coverage(s, bricks)
    for c in cov.colnames:
        s[c] = cov[c]
    covered = np.asarray(s["covered"], bool)
    sc = s[covered]
    src = w12.query_tractor(sc, p, out / "tractor_cache")
    req = w12.required_mags(sc, sc["theta_e"], D056_FJ, p)
    for k, v in req.items():
        sc[k] = v
    images = lenscats.pair_images(sc, src, p.image_radius)
    sc["sep_ls"] = images["sep"]
    sc["dgz"] = pair_colour_difference(images, src)
    sc["colour_match"] = np.abs(sc["dgz"]) <= COLOUR_TOL  # no colour counts as no match
    sc["test_status"] = w12.deflector_test(sc, src, p, images)["test_status"]
    sc["undecided_flags"] = w12.flags_undecided(sc, src, p)
    sc["detectable_typical"] = sc["req_mag_z_typical"] < sc["depth_z"] - p.margin
    sc["detectable_conservative"] = sc["req_mag_z"] < sc["depth_z"] - p.margin
    sc["hennawi_binary"] = binary_match(sc, binq)
    decided = np.isin(sc["test_status"], ["deflector", "none"]) & (sc["undecided_flags"] == "")
    sc["decided"] = decided
    summary = {
        "inputs": manifest,
        "bricks": bmeta,
        "fj_calibration": D056_FJ.__dict__,
        "params": p.__dict__,
        "n_systems": len(s),
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
        "detectable_typical",
        "detectable_conservative",
        "hennawi_binary",
    ]
    sc[keep].write(res / "systems.ecsv", overwrite=True)
    (res / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    man = paths.manifests_dir() / "w12_niq_inputs.json"
    man.write_text(json.dumps({"inputs": manifest, "bricks": bmeta}, indent=2) + "\n")
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
    s.add_argument("--sheet", action="store_true", help="cutouts + contact sheet of 'none' pairs")
    args = ap.parse_args(argv)
    {"screen": cmd_screen}[args.cmd](args)


if __name__ == "__main__":
    main()
