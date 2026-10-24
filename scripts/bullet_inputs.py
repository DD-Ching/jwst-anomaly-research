"""Bullet Cluster: input check for the reference-catalogue extra-image route (no residual read).

The three inputs the route needs (CHANGELOG 2026-10-23 MACS0647 handoff): a public JWST-era image
list, maps of a model fitted to it, and a reference catalogue that deblends arcs from cluster light.
All three come from the CANUCS Bullet Cluster DR1 on CANFAR (GO 4598; Rihtarsic et al. 2026,
arXiv:2601.22245): the Lenstool gold multiple-image list (the model's constraints) is parsed by
``parse_lenstool_images`` and written to ``bullet_gold_images.ecsv``; the full graded list is
parsed by ``parse_multim_full`` and checked against it; the photometry catalogue's FITS header
(HTTP range, first ``HEAD_BYTES`` only) is checked for the columns of the frozen CANUCS DR1 mapping
(``../macs0417_canucsref/PLAN.md``); the lens-model directories are listed (names and sizes only, no
map is read). Writes ``summary.json``.

    python scripts/bullet_inputs.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import elgordo_djaref as eg
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.table import Table

OUT = Path("results/lens_consistency/bullet_inputs")
IMAGES = "bullet_gold_images.ecsv"
VAULT = "https://www.canfar.net/storage/vault"
DR1 = "canucs/BULLET_DR1"
GOLD = f"{VAULT}/file/{DR1}/lensing/bullet_gold_v4.dat"
FULL = f"{VAULT}/file/{DR1}/lensing/bullet_multim_full.dat"
README = f"{VAULT}/file/{DR1}/lensing/README.dat"
PHOT = (
    f"{VAULT}/file/{DR1}/catalogs/"
    "hlsp_canucs_jwst-hst_multi_bullet-clu_multi_v1_photometry-cat.fits"
)
MAPS_DIR = "lensing/best_fit_model/best_fit_maps"
SAMPLES_DIR = "lensing/sample_model/sample_maps"
DEFLECTION_MAPS = ("bulletclu-dpl1-best-50mas.fits.gz", "bulletclu-dpl2-best-50mas.fits.gz")
MAX_SEP_ARCSEC = 0.01  # ASSUMPTION: the two lists are one table written twice
HEAD_BYTES = 230_400  # 80 FITS blocks: the primary and the table header (323 columns)
# The frozen CANUCS DR1 mapping's inputs (../macs0417_canucsref/PLAN.md "Column mapping")
MAPPING_COLS = (
    "RA",
    "DEC",
    "X",
    "Y",
    "A",
    "B",
    "PA",
    "FLUX_KRON_TOTAL_F200W",
    "FLUXERR_KRON_TOTAL_F200W",
    "FLAG_POINTSRC",
)
GRADES = ("gold", "silver", "bronze")
_ID = re.compile(r"^([A-Za-z]+\d+)([a-z]*)\.(\d+)$")


def split_id(image_id: str) -> tuple[str, str]:
    """``K1b.2`` -> (system ``K1b``, galaxy ``K1``); anything else raises."""
    m = _ID.match(image_id)
    if not m:
        raise ValueError(f"image id not parsed: {image_id!r}")
    return m.group(1) + m.group(2), m.group(1)


def _rows(text: str) -> list[list[str]]:
    rows = [ln.split() for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    ids = [p[0] for p in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate image id")
    return rows


def parse_lenstool_images(text: str) -> Table:
    """Lenstool multfile rows ``id ra dec a b theta z mag`` (ICRS degrees); z must be > 0."""
    ids, ra, dec, z = [], [], [], []
    for p in _rows(text):
        if len(p) != 8:
            raise ValueError(f"multfile row not parsed: {' '.join(p)[:80]!r}")
        split_id(p[0])
        ids.append(p[0])
        ra.append(float(p[1]))
        dec.append(float(p[2]))
        z.append(float(p[6]))
    t = Table(
        [ids, ra, dec, z, np.ones(len(ids), bool)],
        names=("image_id", "ra", "dec", "z", "z_spec"),
    )
    if not len(t) or not (t["z"] > 0).all():
        raise ValueError("empty list or a non-positive redshift")
    t["ra"].unit = t["dec"].unit = "deg"
    return t


def parse_multim_full(text: str) -> Table:
    """``bullet_multim_full.dat``: header ``#ID RA(ICRS) DEC(ICRS) specz quality``; every row
    needs 5 fields and a grade in ``GRADES``, and a gold row a redshift > 0."""
    rows = _rows(text)
    for p in rows:
        if len(p) != 5 or p[4] not in GRADES or (p[4] == "gold" and not float(p[3]) > 0):
            raise ValueError(f"row not parsed: {' '.join(p)[:80]!r}")
        split_id(p[0])
    return Table(
        [[p[0] for p in rows]]
        + [[float(p[k]) for p in rows] for k in (1, 2, 3)]
        + [[p[4] for p in rows]],
        names=("image_id", "ra", "dec", "z", "grade"),
    )


def gold_vs_full(gold: Table, full: Table) -> dict:
    """Same ids as the full list's gold rows; position and z differences over the shared ids.
    ``passes``: identical id sets, every separation below ``MAX_SEP_ARCSEC`` and equal z."""
    g = {r["image_id"]: r for r in gold}
    f = {r["image_id"]: r for r in full if r["grade"] == "gold"}
    shared = sorted(set(g) & set(f))
    if not shared:
        raise ValueError("no shared gold id")
    a = SkyCoord([g[k]["ra"] for k in shared], [g[k]["dec"] for k in shared], unit="deg")
    b = SkyCoord([f[k]["ra"] for k in shared], [f[k]["dec"] for k in shared], unit="deg")
    dz = np.abs([g[k]["z"] - f[k]["z"] for k in shared])
    sep = float(a.separation(b).arcsec.max())
    return {
        "n_gold_list": len(g),
        "n_full_gold": len(f),
        "n_shared": len(shared),
        "only_in_gold_list": sorted(set(g) - set(f)),
        "only_in_full_gold": sorted(set(f) - set(g)),
        "max_sep_arcsec": sep,
        "max_abs_dz": float(dz.max()),
        "passes": set(g) == set(f) and sep < MAX_SEP_ARCSEC and float(dz.max()) == 0.0,
    }


def counts(t: Table) -> dict:
    sys_gal = [split_id(i) for i in t["image_id"]]
    return {
        "images": len(t),
        "systems": len({s for s, _ in sys_gal}),
        "galaxies": len({g for _, g in sys_gal}),
        "distinct_z": len(set(np.round(np.asarray(t["z"], float), 4))),
        "z_min": float(np.min(t["z"])),
        "z_max": float(np.max(t["z"])),
    }


def header_check(blob: bytes) -> dict:
    """The first BINTABLE header in the first bytes of the catalogue: rows, columns, the
    mapping's columns present, and their flux unit."""
    hdr = next((h for h in eg.header_cards(blob) if h.get("XTENSION") == "BINTABLE"), None)
    if hdr is None:
        raise ValueError("no complete BINTABLE header in the bytes read")
    n = int(hdr["TFIELDS"])
    cols = [hdr[f"TTYPE{i}"].strip() for i in range(1, n + 1)]
    units = {hdr[f"TTYPE{i}"].strip(): hdr.get(f"TUNIT{i}", "").strip() for i in range(1, n + 1)}
    missing = [c for c in MAPPING_COLS if c not in cols]
    return {
        "rows": int(hdr["NAXIS2"]),
        "n_columns": n,
        "mapping_columns_missing": missing,
        "flux_unit_f200w": units.get("FLUX_KRON_TOTAL_F200W"),
        "passes": not missing and units.get("FLUX_KRON_TOTAL_F200W") == "nJy",
    }


def _get(session: requests.Session, url: str) -> requests.Response:
    r = session.get(url, timeout=120)
    r.raise_for_status()
    return r


def fetch_head(session: requests.Session, url: str = PHOT) -> tuple[bytes, int]:
    """First ``HEAD_BYTES`` and the total size; a server that ignores the range, or a short
    read, is refused."""
    with session.get(
        url, headers={"Range": f"bytes=0-{HEAD_BYTES - 1}"}, timeout=120, stream=True
    ) as r:
        r.raise_for_status()
        if r.status_code != 206:
            raise RuntimeError(f"{url}: range request not honoured (HTTP {r.status_code})")
        blob = r.raw.read(HEAD_BYTES, decode_content=True)
        if len(blob) != HEAD_BYTES:
            raise RuntimeError(f"{url}: {len(blob)} of {HEAD_BYTES} bytes read")
        total = r.headers.get("Content-Range", "").rsplit("/", 1)[-1]
    return blob, int(total) if total.isdigit() else -1


_ROW = re.compile(r'rows\.push\(\["","([^"]+)","([^"]+)"')


def list_dir(session: requests.Session, path: str) -> dict[str, str]:
    """Names and stated sizes of a CANFAR vault directory listing page; an empty page raises."""
    rows = dict(_ROW.findall(_get(session, f"{VAULT}/list/{DR1}/{path}").text))
    if not rows:
        raise ValueError(f"{path}: no listing row parsed")
    return rows


def maps_check(maps: dict[str, str]) -> dict:
    """Both best-fit deflection maps are listed."""
    missing = [m for m in DEFLECTION_MAPS if m not in maps]
    return {"deflection_maps_missing": missing, "passes": not missing}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.parse_args(argv)
    session = requests.Session()
    texts = {
        k: _get(session, u).content for k, u in (("gold", GOLD), ("full", FULL), ("readme", README))
    }
    gold = parse_lenstool_images(texts["gold"].decode())
    full = parse_multim_full(texts["full"].decode())
    blob, total = fetch_head(session)
    maps = list_dir(session, MAPS_DIR)
    samples = list_dir(session, SAMPLES_DIR)
    summary = {
        "sources": {
            "release": "CANUCS Bullet Cluster DR1 (GO 4598), https://niriss.github.io/"
            "data_release_bullet_cluster.html; Rihtarsic et al. 2026, arXiv:2601.22245",
            "files_sha256": {
                u: hashlib.sha256(texts[k]).hexdigest()
                for k, u in (("gold", GOLD), ("full", FULL), ("readme", README))
            },
        },
        "image_list_gold": counts(gold),
        "image_list_full": {g: int((full["grade"] == g).sum()) for g in GRADES},
        "gold_vs_full": gold_vs_full(gold, full),
        "reference_catalogue": {"url": PHOT, "bytes": total, **header_check(blob)},
        "model_best_fit_maps": maps,
        "model_maps_check": maps_check(maps),
        "model_sample_maps": {
            "listing_page_entries": len(samples),  # one listing page; may be paginated
            "examples": dict(sorted(samples.items())[:3]),
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    gold.meta = {
        "provenance": "observed",
        "source": f"{GOLD} (Lenstool multfile, the model's gold constraints; z_spec = true: "
        "the paper states spectroscopic system redshifts, ASSUMPTION per row)",
    }
    gold.write(OUT / IMAGES, format="ascii.ecsv", overwrite=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in summary if k != "model_sample_maps"}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
