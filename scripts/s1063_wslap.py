"""Abell S1063: input check, unit check U and spec-z gate K of the WSLAP+ JWST 2025 deflection maps.

Frozen plan: ``results/lens_consistency/s1063_wslap/PLAN.md``. The WSLAP+ model (Zenodo
10.5281/zenodo.17726321, CC-BY-SA-4.0) is one 1 GB zip; only its two alpha members are fetched, each
by one HTTP range request over its deflated bytes, inflated while streaming and CRC-32-checked
against the zip directory. Checks: C (coverage), U (arcsec vs pixel reading at z_spec), K
(``m0417_images.gate_k``, mean source), K-WLS (reported), R (a DJA v7 imaging field). Image list:
the Chan+2020 systems with a MUSE z_spec. No extra-image, flux or flag reading.

    python scripts/s1063_wslap.py --maps DIR [--grid-step 0.25]

The maps are deleted after the run unless ``--keep``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import struct
import warnings
import zlib
from pathlib import Path

import a2744_uncover as au
import bullet_gatek as bg
import lens_consistency as lc
import m0417_images as mi
import macs1931_inputs as m19
import numpy as np
import requests
import spt0546_inputs as sp
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/s1063_wslap")
IMAGES = Path("results/lens_consistency/abells1063_images/chan20_images.ecsv")
ZIP_URL = (
    "https://zenodo.org/api/records/17726321/files/"
    "LensModel_AS1063_WSLAPplus_JWST_Case3_2025.zip/content"
)
ZIP_DIR = "LensModel_AS1063_WSLAPplus_JWST_Case3_2025/"
MAPS = ("AlphaX_DdsDs1_arcsec.fits", "AlphaY_DdsDs1_arcsec.fits")
Z_LENS = 0.348  # lensmodel HFF table (abells1063)
COSMOLOGY = (70.0, 0.3)  # flat; ASSUMPTION (PLAN.md), only D_ds / D_s enters
PIXEL_ARCSEC = 0.06  # header |CDELT| 1.66667e-5 deg
READINGS = {"A_arcsec": 1.0, "P_pixels": PIXEL_ARCSEC}  # map value -> arcsec
U_MAX_MEDIAN_ARCSEC = 1.0  # ASSUMPTION (PLAN.md rule 2)
GRID_STEP = 0.25
DJA_NAMES = ("abells1063",)
# D-035 addendum b (results/lens_consistency/abells1063_d035/RESULT.md), for comparison only
CATS_V41_RMS_ARCSEC = 18.5
CHUNK = 1 << 20


def spec_images(path: Path = IMAGES) -> Table:
    """The z_spec systems of the Chan+2020 list (``dup_position`` rows dropped), each at its z."""
    t = Table.read(path)
    t = t[np.asarray(t["z_spec"], bool) & ~np.asarray(t["dup_position"], bool)]
    t["system"] = [str(i).split(".")[0] for i in t["image_id"]]
    t["a"] = np.ones(len(t))
    return t


# --- zip members by HTTP range --------------------------------------------------------------


def _range(session: requests.Session, start: int, end: int) -> bytes:
    r = session.get(ZIP_URL, headers={"Range": f"bytes={start}-{end}"}, timeout=120)
    r.raise_for_status()
    if r.status_code != 206:
        raise RuntimeError(f"range request answered {r.status_code}, not 206")
    return r.content


def parse_central_directory(blob: bytes) -> dict[str, dict]:
    """Zip central-directory records (ZIP64 extra fields resolved): name -> crc, sizes, offset."""
    out, k = {}, 0
    while k + 46 <= len(blob) and blob[k : k + 4] == b"PK\x01\x02":
        (method, _, _, crc, csize, usize, nlen, xlen, clen) = struct.unpack(
            "<HHHIIIHHH", blob[k + 10 : k + 34]
        )
        offset = struct.unpack("<I", blob[k + 42 : k + 46])[0]
        name = blob[k + 46 : k + 46 + nlen].decode()
        extra = blob[k + 46 + nlen : k + 46 + nlen + xlen]
        j = 0
        while j + 4 <= len(extra):
            hid, hlen = struct.unpack("<HH", extra[j : j + 4])
            if hid == 1:  # ZIP64: present fields in order usize, csize, offset
                vals = list(struct.unpack(f"<{hlen // 8}Q", extra[j + 4 : j + 4 + hlen]))
                if usize == 0xFFFFFFFF:
                    usize = vals.pop(0)
                if csize == 0xFFFFFFFF:
                    csize = vals.pop(0)
                if offset == 0xFFFFFFFF:
                    offset = vals.pop(0)
            j += 4 + hlen
        out[name] = {
            "method": method,
            "crc": crc,
            "csize": csize,
            "usize": usize,
            "offset": offset,
        }
        k += 46 + nlen + xlen + clen
    if not out:
        raise ValueError("no central-directory record parsed")
    return out


def zip_directory(session: requests.Session) -> tuple[int, dict[str, dict]]:
    size = int(session.head(ZIP_URL, allow_redirects=True, timeout=60).headers["Content-Length"])
    tail = _range(session, max(0, size - 65536), size - 1)
    k = tail.rfind(b"PK\x06\x06")  # ZIP64 end of central directory
    if k >= 0:
        cd_size, cd_off = struct.unpack("<QQ", tail[k + 40 : k + 56])
    else:
        k = tail.rfind(b"PK\x05\x06")
        cd_size, cd_off = struct.unpack("<II", tail[k + 12 : k + 20])
    return size, parse_central_directory(_range(session, cd_off, cd_off + cd_size - 1))


def fetch_member(session: requests.Session, entry: dict, dest: Path) -> str:
    """Inflate one deflated member to ``dest`` while streaming; CRC-32 and size checked.
    Returns the sha256 of the inflated file."""
    if entry["method"] != 8:
        raise ValueError("member is not deflated")
    head = _range(session, entry["offset"], entry["offset"] + 29)
    nlen, xlen = struct.unpack("<HH", head[26:30])
    start = entry["offset"] + 30 + nlen + xlen
    r = session.get(
        ZIP_URL,
        headers={"Range": f"bytes={start}-{start + entry['csize'] - 1}"},
        stream=True,
        timeout=300,
    )
    r.raise_for_status()
    dec, crc, n, sha = zlib.decompressobj(-15), 0, 0, hashlib.sha256()
    with dest.open("wb") as f:
        for chunk in r.iter_content(CHUNK):
            out = dec.decompress(chunk)
            crc, n = zlib.crc32(out, crc), n + len(out)
            sha.update(out)
            f.write(out)
        out = dec.flush()
        crc, n = zlib.crc32(out, crc), n + len(out)
        sha.update(out)
        f.write(out)
    if n != entry["usize"] or crc != entry["crc"]:
        dest.unlink(missing_ok=True)
        raise RuntimeError(f"{dest.name}: size or CRC-32 differs from the zip directory")
    return sha.hexdigest()


# --- checks ---------------------------------------------------------------------------------


def unit_values(model, images: Table) -> dict[str, float | None]:
    """Per system at its z_spec: rms of the ``backtrace_images`` ``dtheta_arcsec``."""
    d = np.asarray(lensmodel.backtrace_images(model, images, {})["dtheta_arcsec"], float)
    sysv = np.asarray(images["system"], str)
    out = {}
    for s in dict.fromkeys(sysv):
        r = mi.rms_of(d[sysv == s])
        out[s] = r if np.isfinite(r) else None
    return out


def unit_check(values: dict[str, dict[str, float | None]]) -> dict:
    """PLAN.md rule 2: accepted when the median per-system rms is <= U_MAX_MEDIAN_ARCSEC (every
    system finite); passes when exactly one reading is accepted."""
    med = {}
    for r, d in values.items():
        v = np.array([np.nan if x is None else x for x in d.values()], float)
        med[r] = float(np.median(v)) if v.size and np.all(np.isfinite(v)) else float("nan")
    accepted = [r for r, v in med.items() if np.isfinite(v) and v <= U_MAX_MEDIAN_ARCSEC]
    return {
        "median_arcsec": med,
        "accepted": accepted,
        "passes": len(accepted) == 1,
        "reading": accepted[0] if len(accepted) == 1 else None,
    }


def dja_check(index_html: str) -> dict:
    """PLAN.md rule 5 (the ``spt0546_inputs.dja_check`` parser with this field's name)."""
    names = sorted(set(re.findall(r"([a-z0-9][a-z0-9_-]*?)-grizli-v7", index_html.lower())))
    if not names:
        raise ValueError("no field parsed from the DJA index")
    hits = [f for f in names if any(n in f for n in DJA_NAMES)]
    return {
        "index": m19.DJA_INDEX,
        "fields_listed": len(names),
        "names_searched": list(DJA_NAMES),
        "matches": hits,
        "passes": bool(hits),
    }


def build_model(ax, ay, header, scale: float, sha: str, name: str) -> lensmodel.MapLensModel:
    """Map values x ``scale`` (arcsec) are the deflection at D_ds / D_s = 1 (record)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header)
    model = lensmodel.MapLensModel(
        ax * scale,
        ay * scale,
        wcs,
        Z_LENS,
        FlatLambdaCDM(H0=COSMOLOGY[0], Om0=COSMOLOGY[1]),
        source=name,
        sha256=sha,
    )
    if not np.isclose(model.pixel_arcsec, PIXEL_ARCSEC, rtol=1e-4):
        raise ValueError(f'{name}: map pixel {model.pixel_arcsec}" != {PIXEL_ARCSEC}"')
    return model


def _read(path: Path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with fits.open(path) as hdul:
            return np.asarray(hdul[0].data, float), hdul[0].header.copy()


def k_block(model, images: Table, grid_step: float) -> dict:
    half = lc.grid_half_width(model, images)
    grid = lensmodel.DeflectionGrid.compute(model, half, grid_step)
    bt = lensmodel.backtrace_images(model, images, {})
    ip = lensmodel.imageplane_residuals(model, grid, bt)
    ipw = au.imageplane(model, grid, bt)
    return {
        "grid_half_width_arcsec": half,
        "ip": ip,
        "ip_wls": ipw,
        "gate_K": mi.gate_k(ip),
        "K_wls": mi.gate_k(ipw),
        "per_system": bg.per_system(ip),
        "per_system_wls": bg.per_system(ipw),
    }


def run(maps: Path, out: Path = OUT, grid_step: float = GRID_STEP) -> dict:
    session = requests.Session()
    size, directory = zip_directory(session)
    maps.mkdir(parents=True, exist_ok=True)
    shas = {m: fetch_member(session, directory[ZIP_DIR + m], maps / m) for m in MAPS}
    images = spec_images()
    (ax, hx), (ay, _) = (_read(maps / m) for m in MAPS)
    model_sha = hashlib.sha256((shas[MAPS[0]] + shas[MAPS[1]]).encode()).hexdigest()
    models = {
        r: build_model(ax, ay, hx, s, model_sha, f"abells1063-wslap ({r})")
        for r, s in READINGS.items()
    }
    summary = {
        "plan": str(out / "PLAN.md"),
        "zip_bytes": size,
        "members": {
            m: {
                "sha256": shas[m],
                "bytes": directory[ZIP_DIR + m]["usize"],
                "deflated_bytes": directory[ZIP_DIR + m]["csize"],
                "crc32": f"{directory[ZIP_DIR + m]['crc']:08x}",
            }
            for m in MAPS
        },
        "cosmology": {"H0": COSMOLOGY[0], "Om0": COSMOLOGY[1], "z_lens": Z_LENS},
        "image_list": {
            "source": str(IMAGES),
            "images": len(images),
            "systems": len(set(images["system"])),
        },
        "map_header": {
            "shape": list(models["A_arcsec"].shape),
            "pixel_arcsec": round(models["A_arcsec"].pixel_arcsec, 6),
            "bunit": hx.get("BUNIT"),
        },
        "check_C": sp.coverage(models["A_arcsec"], images),
    }
    if summary["check_C"]["passes"]:
        vals = {r: unit_values(m, images) for r, m in models.items()}
        u = unit_check(vals)
        summary["check_U"] = {**u, "per_system": vals}
        if u["passes"]:
            k = k_block(models[u["reading"]], images, grid_step)
            lc._write(k["ip"], out / "imageplane_zspec.ecsv")
            lc._write(k["ip_wls"], out / "imageplane_zspec_wls.ecsv")
            summary["gate_K"] = k["gate_K"]
            summary["reported"] = {
                "grid_step_arcsec": grid_step,
                "grid_half_width_arcsec": k["grid_half_width_arcsec"],
                "K_wls_not_a_gate": k["K_wls"],
                "per_system": k["per_system"],
                "per_system_wls": k["per_system_wls"],
                "cats_v41_rms_arcsec_d035b": CATS_V41_RMS_ARCSEC,
            }
    summary["check_R"] = dja_check(m19._get(session, m19.DJA_INDEX).text)
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--grid-step", type=float, default=GRID_STEP)
    ap.add_argument("--keep", action="store_true", help="keep the fetched maps")
    a = ap.parse_args(argv)
    try:
        s = run(a.maps, a.out, a.grid_step)
    finally:
        if not a.keep:
            for f in MAPS:
                (a.maps / f).unlink(missing_ok=True)
    print(json.dumps({k: v for k, v in s.items() if k != "reported"}, indent=1)[:4000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
