"""W2 screen of the 73 lensed quasars with published time delays (Gaia GraL X, Table 2).

Question: does any system have both a reliable measured delay and high-resolution imaging deep
enough to show an ordinary lens galaxy, yet no detected deflector? Published products only.

    python scripts/w2_gralx_delay_screen.py parse    # Table 2 from the pinned arXiv source
    python scripts/w2_gralx_delay_screen.py coords   # SIMBAD (Harvard mirror) or name-derived
    python scripts/w2_gralx_delay_screen.py dedup    # against D-056 / D-064 / later W2 records
    python scripts/w2_gralx_delay_screen.py assemble # deflectors.json (curated) -> summary.json

`deflectors.json` is curated by hand from the original papers (one source per system); this script
only joins and counts it. Labels: delays and classes are published (observed); matches derived.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import tarfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

from jwst_anomaly import paths, photometry

SRC = "https://arxiv.org/src/2603.12810v1"  # Gaia GraL X (Ducourant et al.), A&A, arXiv v1
SRC_SHA = "270f955601361c71a4278df6b61c1eb3e0eb8efcd585b717fa857973c3cb187e"
TEX = "aa58049-25corr.tex"
SIMBAD = "https://simbad.cfa.harvard.edu/simbad/sim-id"  # CDS host resets; Harvard mirror
PAIRS = ("AB", "AC", "AD", "BC", "BD", "CD")
# Names whose plain form SIMBAD resolves wrongly or not at all (checked by hand, 2026-11-03).
SIMBAD_ALIAS = {
    "B1608+656": "CLASS B1608+656",  # plain "B1608+656" resolves to 3C 330 G2
    "QJ0158-4325": "QSO J0158-4325",
    "Q0142-100": "Q0142-100",
}
B1950 = {
    "Q0142-100", "B0218+357", "HE0047-1756", "HE0230-2130", "HE0435-1223", "HS0818+1227",
    "APM08279+5255", "SBS0909+532", "Q0957+561", "B1030+074", "HE1104-1805", "PG1115+080",
    "H1413+117", "B1422+231", "SBS1520+530", "B1600+434", "B1608+656", "PKS1830-211",
    "Q2138-431", "HE2149-2745", "HS2209+1914",
}  # fmt: skip
PREFIX = re.compile(r"^(DES|PS|WISE|WGD|WG|GRAL|GraL|SDSS|2MASS|2M|VDES|RXJ|RX|FBQ|WFI|QJ|Q)")
JNAME = re.compile(r"J?(\d{2})(\d{2})(\d{2}(?:\.\d+)?)?([+-])(\d{2})(\d{2})?(\d{2}(?:\.\d+)?)?")


def out_dir() -> Path:
    d = paths.repo_root() / "results" / "w2_gralx_delay_screen"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write(name: str, obj) -> Path:
    p = out_dir() / name
    p.write_text(json.dumps(obj, indent=1, ensure_ascii=False) + "\n")
    return p


# ----------------------------------------------------------------------------- Table 2


def parse_rows(tex: str) -> list[dict]:
    """Rows of the time-delay longtable, printed and %-commented alike (the arXiv PDF prints a
    portion; the source carries all 73)."""
    start = tex.index(r"\label{timedelay_table}")
    end = tex.index(r"\end{longtable}", start)
    rows = []
    for line in tex[start:end].splitlines():
        s = line.strip().lstrip("%").strip()
        if "&" not in s or s.startswith("Name") or "multicolumn" in s:
            continue
        parts = [p.strip() for p in s.split("&")]
        if len(parts) < 8:
            continue
        delays = {}
        for k, p in zip(PAIRS, parts[1:7], strict=True):
            m = re.search(r"([-+]?\d+\.?\d*)\s*\\pm\s*(\d+\.?\d*)", p)
            v = re.search(r"([-+]?\d+\.?\d*)", p)
            if m:
                delays[k] = {"dt": float(m.group(1)), "err": float(m.group(2))}
            elif v:
                delays[k] = {"dt": float(v.group(1)), "err": None}
        author = re.sub(r"\\\\.*$", "", parts[7]).replace("\\_", "_").strip()
        rows.append({"name": parts[0], "delays": delays, "ref": author})
    return rows


def parse() -> dict:
    tar = photometry.fetch_catalog(SRC, SRC_SHA)
    with tarfile.open(tar, "r:*") as t:
        tex = t.extractfile(TEX).read().decode("utf-8", "ignore")
    rows = parse_rows(tex)
    res = {
        "label": "published (Gaia GraL X Table 2, transcribed from the arXiv v1 LaTeX source)",
        "source": f"{SRC} (sha256 {SRC_SHA}), {TEX}",
        "n_systems": len(rows),
        "n_delays_with_error": sum(
            1 for r in rows for d in r["delays"].values() if d["err"] is not None
        ),
        "n_delays_without_error": sum(
            1 for r in rows for d in r["delays"].values() if d["err"] is None
        ),
        "rows": rows,
    }
    _write("table2.json", res)
    return res


# ----------------------------------------------------------------------------- coordinates


def _simbad(ident: str):
    url = f"{SIMBAD}?" + urllib.parse.urlencode({"Ident": ident, "output.format": "ASCII"})
    with urllib.request.urlopen(url, timeout=30) as r:  # noqa: S310
        s = r.read().decode("utf-8", "ignore")
    m = re.search(r"Coordinates\(ICRS[^:]*:\s*([0-9 .]+?)\s+([+-][0-9 .]+?)\s+\(", s)
    if not m:
        return None
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(m.group(1).strip(), m.group(2).strip(), unit=(u.hourangle, u.deg))
    mid = re.search(r"Object (.+?)  ---", s)
    return float(c.ra.deg), float(c.dec.deg), mid.group(1).strip() if mid else ""


def name_coords(name: str):
    """Coarse J2000 position from a J-designation (HHMM[SS]±DD[MM[SS]]); None for B1950 names."""
    key = name.replace(" ", "").replace("*", "")
    if key in B1950:
        return None
    m = JNAME.match(PREFIX.sub("", key))
    if not m:
        return None
    hh, mm, ss, sg, dd, dm, ds = m.groups()
    ra = 15.0 * (int(hh) + int(mm) / 60 + (float(ss) if ss else 0.0) / 3600)
    dec = int(dd) + (int(dm) / 60 if dm else 0.0) + (float(ds) if ds else 0.0) / 3600
    return ra, (-dec if sg == "-" else dec)


def coords() -> dict:
    rows = json.loads((out_dir() / "table2.json").read_text())["rows"]
    out = {}
    for r in rows:
        n = r["name"].replace("*", "").strip()
        nc = name_coords(n)
        if nc is not None and n not in SIMBAD_ALIAS:
            out[r["name"]] = {"ra": nc[0], "dec": nc[1], "source": "name (coarse)"}
            continue
        hit = None
        for ident in (SIMBAD_ALIAS.get(n, n), "QSO B" + re.sub(r"^[A-Za-z]+", "", n)):
            try:
                hit = _simbad(ident)
            except OSError:
                hit = None
            if hit:
                out[r["name"]] = {"ra": hit[0], "dec": hit[1], "source": f"SIMBAD {ident}"}
                out[r["name"]]["simbad_main_id"] = hit[2]
                break
            time.sleep(0.3)
        if not hit:
            out[r["name"]] = None
    _write("coords.json", {"label": "derived (SIMBAD Harvard mirror or name)", "coords": out})
    return out


# ----------------------------------------------------------------------------- dedup


def _desig(ra: float, dec: float) -> tuple[int, int]:
    h = ra / 15.0
    a = abs(dec)
    return int(h) * 60 + int((h - int(h)) * 60), int(np.sign(dec) or 1) * (
        int(a) * 60 + int((a - int(a)) * 60)
    )


def prior_systems() -> list[tuple]:
    from astropy.table import Table

    root = paths.repo_root() / "results"
    prior = []
    for r in Table.read(root / "w12_niq" / "systems.ecsv"):
        prior.append(("w12_niq (D-064)", str(r["name"]), float(r["ra"]), float(r["dec"]),
                      f"{r['group']}/{r['catalogue']}/{r['comment']}"))  # fmt: skip
    for s in json.loads((root / "w2_desi_pairs" / "pairs.json").read_text())["systems"]:
        if s.get("ra") is not None:
            prior.append(("w2_desi_pairs", s["name"], float(s["ra"]), float(s["dec"]),
                          f"{s.get('role')}/{s.get('comment')}"))  # fmt: skip
    prior.append(("w2_j2244_hst", "J2244-0550", 341.0473, -5.8472, "rejected UQP"))
    prior.append(("w2_j0740_gaia_delay", "J0740+2926", 115.0560, 29.4464, "binary / NIQ"))
    return prior


def dedup() -> dict:
    cs = json.loads((out_dir() / "coords.json").read_text())["coords"]
    prior = prior_systems()
    keys = [_desig(p[2], p[3]) for p in prior]
    res = {}
    for n, c in cs.items():
        if not c:
            res[n] = None
            continue
        a = _desig(c["ra"], c["dec"])
        hits = sorted(
            {
                (prior[i][0], prior[i][1], prior[i][4])
                for i, b in enumerate(keys)
                if abs(a[0] - b[0]) <= 1 and abs(a[1] - b[1]) <= 1
            }
        )
        res[n] = [list(h) for h in hits]
    _write(
        "dedup.json",
        {
            "label": "derived: HHMM±DDMM designation match (±1 in the last digit) against the "
            "repo's earlier W1/W2 records; D-056's decided systems are not stored as a list, "
            "only those named in DECISIONS.md",
            "matches": res,
            "n_with_prior": sum(1 for v in res.values() if v),
        },
    )
    return res


# ----------------------------------------------------------------------------- assemble

CATEGORIES = ("galaxy_seen", "quantitative_limit", "shallow_or_blended", "catalogue_blank")


def assemble() -> dict:
    t2 = json.loads((out_dir() / "table2.json").read_text())
    dfl = json.loads((out_dir() / "deflectors.json").read_text())["systems"]
    dd = json.loads((out_dir() / "dedup.json").read_text())["matches"]
    names = [r["name"] for r in t2["rows"]]
    missing = [n for n in names if n not in dfl]
    counts = {k: 0 for k in CATEGORIES}
    rows = []
    for r in t2["rows"]:
        d = dfl.get(r["name"], {})
        cat = d.get("category", "catalogue_blank")
        counts[cat] = counts.get(cat, 0) + 1
        qualifies = (
            cat == "quantitative_limit"
            and d.get("delay_reliable") is True
            and d.get("imaging_high_res") is True
        )
        rows.append({"name": r["name"], "category": cat, "qualifies": qualifies,
                     "prior_repo": dd.get(r["name"]) or []})  # fmt: skip
    res = {
        "label": "derived from published classes (deflectors.json) and Table 2",
        "n_systems": len(names),
        "not_curated": missing,
        "counts": counts,
        "qualifying": [x["name"] for x in rows if x["qualifies"]],
        "seen_with_z_lens": sum(
            1
            for d in dfl.values()
            if d["category"] == "galaxy_seen" and d.get("z_lens") is not None
        ),
        "seen_ground_only_without_z_lens": sorted(
            n
            for n, d in dfl.items()
            if d["category"] == "galaxy_seen"
            and not d.get("imaging_high_res")
            and d.get("z_lens") is None
        ),
        "seen_automated_fit_only": sorted(
            n
            for n, d in dfl.items()
            if "automated" in (d.get("note", "") + d.get("instrument", ""))
        ),
        "table2_delay_not_a_measurement": sorted(
            n
            for n, d in dfl.items()
            if d.get("delay_check", {}).get("reliable_lens_delay") is False
        ),
        "delay_originals_read": sorted(n for n, d in dfl.items() if "delay_check" in d),
        "rows": rows,
    }
    _write("summary.json", res)
    return res


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("cmd", choices=["parse", "coords", "dedup", "assemble"])
    a = p.parse_args(argv)
    fn = {"parse": parse, "coords": coords, "dedup": dedup, "assemble": assemble}[a.cmd]
    r = fn()
    print(
        json.dumps(r if a.cmd != "parse" else {k: r[k] for k in r if k != "rows"}, indent=1)[:3000]
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
