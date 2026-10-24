"""MACS J1931.8-2635: input check for the fitted-model route (no model prediction read).

The three inputs the route needs (CHANGELOG 2026-10-24 Bullet handoff): a public JWST-era image
list, maps of a model fitted to it, and a reference catalogue that deblends arcs from cluster light.
The model is the VENUS AstroLensPy v2.1 model (Allingham et al. 2026, arXiv:2602.14074) on Zenodo
(10.5281/zenodo.18646525). Its image list is parsed by ``parse_arcs`` and written to
``m1931_images.ecsv``; the spectroscopic families are matched to the record's MUSE redshift list;
the FITS headers of the two deflection maps are read from the first ``HEAD_BYTES`` of each zip
(HTTP range, inflated; no map value is read); the DJA v7 imaging index is searched for the field.
Writes ``summary.json``.

    python scripts/macs1931_inputs.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import struct
import zlib
from pathlib import Path

import elgordo_djaref as eg
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.table import Table

OUT = Path("results/lens_consistency/macs1931_inputs")
IMAGES = "m1931_images.ecsv"
RECORD = "https://zenodo.org/api/records/18646525"
FILES = "https://zenodo.org/records/18646525/files"
README = "README_Allingham_m1931_AstroLensPy_Allingham26.txt"
ARCS = "m1931_AstroLensPy_v2.1_arcs_WCS.txt"
MUSE = "m1931_redshifts_MUSE.txt"
CHIRES = "chires.txt"
DEFLECTION_MAPS = ("Alpha_x.fits.zip", "Alpha_y.fits.zip")
DJA_INDEX = "https://dawn-cph.github.io/dja/imaging/v7/"
DJA_NAMES = ("macs1931", "m1931", "j1931")  # field names DJA uses: cluster name or j<RA><Dec>
# The arcs file header: "Systems 1-7, 100 and 101 are spec-z identified"
SPEC_SYSTEMS = ("1", "2", "3", "4", "5", "6", "7", "100", "101")
MUSE_MATCH_ARCSEC = 0.5  # ASSUMPTION: an image and its MUSE detection within 0.5 arcsec
MUSE_DZ = 0.01  # ASSUMPTION: the list z and the MUSE z agree to 0.01
MUSE_MIN_CONF = 2  # the README: flag > 2 trusted, flag 2 "with caution"
HEAD_BYTES = 200_000  # compressed bytes; inflate to well over the two FITS headers
README_UNIT = "arcsec"  # the README: "deflection field map (x and y axis, in arcsec)"
_ID = re.compile(r"^(c?)(\d+)\.(\d+)$")


def parse_arcs(text: str) -> Table:
    """AstroLensPy arcs rows ``id ra dec a b theta z mag``. A ``c`` id prefix marks a candidate
    the model did not use (``used`` false; the README: 51 images used); ``z`` 0 means free."""
    ids, system, ra, dec, z, used = [], [], [], [], [], []
    for ln in text.splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        p = ln.split()
        m = _ID.match(p[0]) if p else None
        if len(p) != 8 or not m:
            raise ValueError(f"arcs row not parsed: {ln[:80]!r}")
        ids.append(p[0])
        system.append(m.group(2))
        ra.append(float(p[1]))
        dec.append(float(p[2]))
        z.append(float(p[6]))
        used.append(not m.group(1))
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate image id")
    t = Table(
        [ids, system, ra, dec, z, used], names=("image_id", "system", "ra", "dec", "z", "used")
    )
    t["z_spec"] = [s in SPEC_SYSTEMS for s in t["system"]]
    if (t["z"] < 0).any() or not (t["z"][t["z_spec"]] > 0).all():
        raise ValueError("negative z, or a spectroscopic family without z")
    t["ra"].unit = t["dec"].unit = "deg"
    return t


def list_counts(t: Table) -> dict:
    u = t[t["used"]]
    spec = u[u["z_spec"]]
    one_z = [
        s
        for s in sorted(set(spec["system"]), key=int)
        if len(set(spec["z"][spec["system"] == s])) == 1
    ]
    photo = u[~u["z_spec"]]
    mixed = [
        s
        for s in sorted(set(photo["system"]), key=int)
        if len(set(photo["z"][photo["system"] == s])) > 1
    ]
    return {
        "rows": len(t),
        "used_images": len(u),
        "used_systems": len(set(u["system"])),
        "spec_images": len(spec),
        "spec_systems": len(set(spec["system"])),
        "spec_systems_single_z": len(one_z),
        "candidates_not_used": [str(i) for i in t["image_id"][~t["used"]]],
        "photo_systems": len(set(photo["system"])),
        "photo_systems_mixed_z": mixed,
        "photo_images_z0": int((photo["z"] == 0).sum()),
    }


def parse_muse(text: str) -> Table:
    """MUSE list rows ``Id RA Dec a b theta z mag conf xcorr``."""
    rows = [ln.split() for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not rows or any(len(p) != 10 for p in rows):
        raise ValueError("MUSE row not parsed")
    return Table(
        [[p[0] for p in rows]]
        + [[float(p[k]) for p in rows] for k in (1, 2, 6)]
        + [[int(p[8]) for p in rows]],
        names=("muse_id", "ra", "dec", "z", "conf"),
    )


def muse_check(t: Table, muse: Table) -> dict:
    """Each used spectroscopic image: nearest MUSE entry (conf >= ``MUSE_MIN_CONF``) within
    ``MUSE_MATCH_ARCSEC`` with |dz| <= ``MUSE_DZ``. ``passes``: every spectroscopic system has
    at least one matched image."""
    s = t[t["used"] & t["z_spec"]]
    m = muse[muse["conf"] >= MUSE_MIN_CONF]
    a = SkyCoord(s["ra"], s["dec"], unit="deg")
    b = SkyCoord(m["ra"], m["dec"], unit="deg")
    idx, sep, _ = a.match_to_catalog_sky(b)
    dz = np.abs(np.asarray(m["z"])[idx] - np.asarray(s["z"]))
    ok = (sep.arcsec <= MUSE_MATCH_ARCSEC) & (dz <= MUSE_DZ)
    systems = sorted(set(s["system"]), key=int)
    matched = sorted({str(x) for x in s["system"][ok]}, key=int)
    return {
        "muse_entries": len(muse),
        "muse_entries_conf_ok": len(m),
        "spec_images": len(s),
        "images_matched": int(ok.sum()),
        "systems_matched": matched,
        "systems_unmatched": [x for x in systems if x not in matched],
        "max_abs_dz_matched": float(dz[ok].max()) if ok.any() else None,
        "passes": matched == systems,
    }


def chires_rms(text: str, used_ids: list[str]) -> dict:
    """RMS over the used images of each per-image residual column of ``chires.txt``."""
    lines = [ln.split() for ln in text.splitlines() if ln.strip()]
    head, rows = lines[0], {p[0]: [float(v) for v in p[1:]] for p in lines[1:] if len(p) == 4}
    cols = [c for c in head[1:] if not c.startswith("[")]
    have = [i for i in used_ids if i in rows]
    if not have:
        raise ValueError("no used image found in chires.txt")
    v = np.array([rows[i] for i in have])
    return {
        "images": len(have),
        "missing": [i for i in used_ids if i not in rows],
        **{c: round(float(np.sqrt(np.mean(v[:, k] ** 2))), 4) for k, c in enumerate(cols)},
    }


def zip_first_member_header(blob: bytes) -> dict:
    """FITS headers of the first member of a zip, from its leading bytes (deflate or stored)."""
    sig, flags, meth, nl, el = struct.unpack("<I2xHH16xHH", blob[:30])
    if sig != 0x04034B50:
        raise ValueError("not a zip local file header")
    if flags & 1 or meth not in (0, 8):
        raise ValueError(f"zip member encrypted or compressed with method {meth}")
    data = blob[30 + nl + el :]
    raw = zlib.decompressobj(-15).decompress(data) if meth == 8 else data
    heads = eg.header_cards(raw)
    img = next((h for h in heads if h.get("NAXIS", 0) == 2), None)
    if img is None:
        raise ValueError("no complete image header in the bytes read")
    keys = ("NAXIS1", "NAXIS2", "CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2")
    keys += ("CDELT1", "CDELT2", "RADESYS", "BUNIT", "MODEL", "Z_CL", "EXTNAME")
    out = {"member": blob[30 : 30 + nl].decode(), **{k: img.get(k) for k in keys}}
    # COMMENT text wraps at 72 characters, often mid-word: pad each card, then collapse spaces
    out["comment"] = " ".join("".join(str(c).ljust(72) for c in img.get("COMMENT", [])).split())
    out["pixel_arcsec"] = round(abs(float(img["CDELT2"])) * 3600, 6)
    return out


def unit_check(headers: list[dict]) -> dict:
    """The deflection unit in the map headers against the README; ``agree`` false is a conflict
    a later plan must resolve before reading any value."""
    units = sorted({str(h["BUNIT"]).strip() for h in headers})
    return {"bunit": units, "readme": README_UNIT, "agree": units == [README_UNIT]}


def dja_check(index_html: str) -> dict:
    # full field names, which may hold '-' or '_' (cluster-j002748p2616, spitzer_idf)
    fields = sorted(set(re.findall(r"([a-z0-9][a-z0-9_-]*?)-grizli-v7", index_html.lower())))
    if not fields:
        raise ValueError("no field parsed from the DJA index")
    hits = [f for f in fields if any(n in f for n in DJA_NAMES)]
    return {"index": DJA_INDEX, "fields_listed": len(fields), "matches": hits, "passes": bool(hits)}


def _get(session: requests.Session, url: str, **kw) -> requests.Response:
    r = session.get(url, timeout=120, **kw)
    r.raise_for_status()
    return r


def fetch_head(session: requests.Session, url: str) -> tuple[bytes, int]:
    """First ``HEAD_BYTES``; a server that ignores the range, or a short read, is refused."""
    with session.get(
        url, headers={"Range": f"bytes=0-{HEAD_BYTES - 1}"}, timeout=120, stream=True
    ) as r:
        r.raise_for_status()
        if r.status_code != 206:
            raise RuntimeError(f"{url}: range request not honoured (HTTP {r.status_code})")
        blob = r.raw.read(HEAD_BYTES, decode_content=True)
        total = r.headers.get("Content-Range", "").rsplit("/", 1)[-1]
    if len(blob) != HEAD_BYTES:
        raise RuntimeError(f"{url}: {len(blob)} of {HEAD_BYTES} bytes read")
    return blob, int(total) if total.isdigit() else -1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.parse_args(argv)
    session = requests.Session()
    record = _get(session, RECORD).json()
    listed = {f["key"]: {"bytes": f["size"], "checksum": f["checksum"]} for f in record["files"]}
    texts = {
        k: _get(session, f"{FILES}/{k}?download=1").content for k in (README, ARCS, MUSE, CHIRES)
    }
    for k, blob in texts.items():
        if "md5:" + hashlib.md5(blob).hexdigest() != listed[k]["checksum"]:
            raise RuntimeError(f"{k}: md5 differs from the record")
    arcs = parse_arcs(texts[ARCS].decode())
    muse = parse_muse(texts[MUSE].decode())
    headers = []
    for name in DEFLECTION_MAPS:
        blob, total = fetch_head(session, f"{FILES}/{name}?download=1")
        headers.append({"file": name, "zip_bytes": total, **zip_first_member_header(blob)})
    dja = dja_check(_get(session, DJA_INDEX).text)
    used_ids = [str(i) for i in arcs["image_id"][arcs["used"]]]
    summary = {
        "sources": {
            "record": f"{RECORD} (doi {record.get('doi')}, published "
            f"{record['metadata'].get('publication_date')})",
            "paper": "Allingham et al. 2026, arXiv:2602.14074 (AstroLensPy, VENUS)",
            "files_sha256": {k: hashlib.sha256(v).hexdigest() for k, v in texts.items()},
        },
        "record_files": len(listed),
        "deflection_maps_listed": {m: listed.get(m, {}).get("bytes") for m in DEFLECTION_MAPS},
        "image_list": list_counts(arcs),
        "muse_check": muse_check(arcs, muse),
        "chires_rms_arcsec": chires_rms(texts[CHIRES].decode(), used_ids),
        "deflection_headers": headers,
        "deflection_unit": unit_check(headers),
        "reference_catalogue_dja_v7": dja,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    arcs.meta = {
        "provenance": "observed (positions); z: spec-z for z_spec rows, else photo-z or "
        "model-optimised per the arcs file header; 0 = free",
        "source": f"{FILES}/{ARCS}",
    }
    arcs.write(OUT / IMAGES, format="ascii.ecsv", overwrite=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
