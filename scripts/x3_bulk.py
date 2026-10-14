"""X3: stream the Gaia DR3 bulk epoch-photometry files that cover GAPS, keeping chosen sources.

Unit file: docs/state/units/X3.md (G1). The bulk files hold one row per source with array-valued
per-transit columns. Each kept row is rewritten as the per-transit CSV that DataLink returns
(``DATAMODEL_STANDARD``), so ``gaia_mulens.parse_epoch_csv`` and ``x3_g1.parse_bands`` run unchanged
on bulk and DataLink data. DataLink's ``rejected_by_photometry`` column is not in the bulk files; it
is rebuilt as ``g_transit_n_obs == 0`` (ASSUMPTION, checked on one transit; the ``verify``
subcommand compares the parsed light curves of every cached dev source in a file).

Files are streamed (gunzip on the fly, MD5 checked against the release's ``_MD5SUM.txt`` values
below) and never stored (CLAUDE.md, Cloud disk decision).

Subcommands:
  ``verify``  stream one file; compare bulk vs cached DataLink light curves of dev sources in it.
  ``heldout-ids``  held-out catalogue rows (G <= 19, MOD(source_id, 10) != 0): ids, position, RUWE.
  ``screen``  run the frozen G1 screen (x3_g1.screen_full) on held-out sources, file by file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import sys
import time
import urllib.request
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import gaia_mulens as gm

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


g1 = _load("x3_g1")
x3 = g1.x3

BASE_URL = "https://cdn.gea.esac.esa.int/Gaia/gdr3/Photometry/epoch_photometry/"
# The five files whose level-8 nested HEALPix ranges cover the 5.6 deg GAPS cone
# (name, MD5 from _MD5SUM.txt, size in bytes); listed 2026-10-14.
FILES = (
    ("EpochPhotometry_040908-041937.csv.gz", "fef1d3a4143c77aa049c6af1c0b2ee89", 938261666),
    ("EpochPhotometry_041938-042766.csv.gz", "1c146b3b1085641e47ef1ec525912a2b", 1407649739),
    ("EpochPhotometry_042767-043485.csv.gz", "57a46b66cefda42e29db3ca5177155c6", 2217846486),
    ("EpochPhotometry_043486-044091.csv.gz", "f46c522910e7e2f7757be0bb4ef89b5a", 1148105703),
    ("EpochPhotometry_044092-044563.csv.gz", "4a063775af8c29cb6c3cb1aad3ebdde5", 967524071),
)

# Per-transit columns written for the parsers (DataLink names).
OUT_COLS = (
    "source_id",
    "transit_id",
    "g_transit_time",
    "g_transit_flux",
    "g_transit_flux_error",
    "g_transit_flux_over_error",
    "g_transit_mag",
    "bp_obs_time",
    "bp_flux",
    "bp_flux_error",
    "bp_mag",
    "rp_obs_time",
    "rp_flux",
    "rp_flux_error",
    "rp_mag",
    "variability_flag_g_reject",
    "variability_flag_bp_reject",
    "variability_flag_rp_reject",
    "rejected_by_photometry",
)


def file_range(name: str) -> tuple[int, int]:
    lo, hi = name.split("_")[1].split(".")[0].split("-")
    return int(lo), int(hi)


def level8(source_id: int) -> int:
    return int(source_id) >> 43


def _arr(text: str) -> list[str]:
    """Array cell ("[1.0, NaN, ...]" or "") -> list of tokens."""
    s = text.strip()
    if not s:
        return []
    s = s.strip("[]")
    return [x.strip().strip('"') for x in s.split(",")] if s else []


def _val(tok: str | None) -> str:
    if tok is None or tok in ("NaN", "nan", "null", ""):
        return ""
    if tok in ("True", "False"):
        return tok.lower()
    return tok


def to_datalink_csv(row: dict) -> str:
    """One bulk row -> the per-transit CSV text DataLink returns (the columns in OUT_COLS)."""
    tid = _arr(row["transit_id"])
    n = len(tid)
    cols = {c: _arr(row.get(c, "")) for c in OUT_COLS if c not in ("source_id", "transit_id")}
    nobs = _arr(row["g_transit_n_obs"])
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(OUT_COLS)
    for i in range(n):
        rec = [row["source_id"], tid[i]]
        for c in OUT_COLS[2:]:
            if c == "rejected_by_photometry":
                rec.append("true" if i < len(nobs) and nobs[i] in ("0", "") else "false")
                continue
            v = cols[c]
            rec.append(_val(v[i] if i < len(v) else None))
        w.writerow(rec)
    return buf.getvalue()


def stream(name: str, keep, *, check_md5: bool = True):
    """Yield (source_id, DataLink-style CSV text) for rows whose source_id is in ``keep``.

    ``keep`` is a set of ints. After the last row the stream's MD5 is compared with FILES; a
    mismatch raises. Returns per-file statistics through StopIteration's value.
    """
    md5_ref = {n: m for n, m, _ in FILES}[name]
    csv.field_size_limit(1 << 30)
    t1 = time.time()
    resp = urllib.request.urlopen(BASE_URL + name, timeout=600)  # noqa: S310
    dec = zlib.decompressobj(wbits=31)
    md5 = hashlib.md5()  # noqa: S324 (archive checksum, not security)
    header = None
    buf = b""
    stats = {"bytes": 0, "rows": 0, "kept": 0}
    while True:
        chunk = resp.read(1 << 22)
        if not chunk:
            break
        md5.update(chunk)
        stats["bytes"] += len(chunk)
        buf += dec.decompress(chunk)
        lines = buf.split(b"\n")
        buf = lines.pop()
        for line in lines:
            if not line or line.startswith(b"#"):
                continue
            if line.startswith(b"solution_id"):
                header = line.decode().split(",")
                continue
            stats["rows"] += 1
            sid = int(line.split(b",", 2)[1])
            if sid not in keep:
                continue
            stats["kept"] += 1
            vals = next(csv.reader([line.decode()]))
            yield sid, to_datalink_csv(dict(zip(header, vals, strict=True)))
    resp.close()
    stats["seconds"] = round(time.time() - t1, 1)
    stats["md5_ok"] = md5.hexdigest() == md5_ref
    if check_md5 and not stats["md5_ok"]:
        raise RuntimeError(f"{name}: MD5 mismatch")
    return stats


# ----------------------------------------------------------------------------- verify on dev
def _lc_arrays(text: str, sid: str):
    tab = gm.parse_epoch_csv(text)
    bands = g1.parse_bands(text)
    return (
        np.asarray(tab["time"], float),
        np.asarray(tab["mag"], float),
        np.asarray(tab["mag_err"], float),
        bands,
    )


def compare(text_bulk: str, text_dl: str, sid: str) -> dict:
    a, b = _lc_arrays(text_bulk, sid), _lc_arrays(text_dl, sid)
    rec = {"source_id": sid, "n_g_bulk": a[0].size, "n_g_datalink": b[0].size}
    same = all(np.array_equal(x, y, equal_nan=True) for x, y in zip(a[:3], b[:3], strict=True))
    for band in ("BP", "RP"):
        same &= all(
            np.array_equal(x, y, equal_nan=True)
            for x, y in zip(a[3][band], b[3][band], strict=True)
        )
    rec["identical"] = bool(same)
    return rec


def verify(name: str) -> dict:
    lo, hi = file_range(name)
    survey = x3.GapsDev()
    cached = {
        int(p.stem)
        for p in survey._cache.glob("*.csv")
        if p.stem.isdigit() and lo <= level8(int(p.stem)) <= hi
    }
    rows, gen = [], stream(name, cached)
    while True:
        try:
            sid, text = next(gen)
        except StopIteration as stop:
            stats = stop.value
            break
        rows.append(compare(text, survey._lc_path(str(sid)).read_text(), str(sid)))
    rec = {
        "file": name,
        "n_cached_dev_in_file": len(cached),
        "n_found_in_bulk": len(rows),
        "n_identical": sum(r["identical"] for r in rows),
        "differing": [r for r in rows if not r["identical"]][:20],
        "stream": stats,
        "provenance": "observed (Gaia DR3 bulk vs DataLink epoch photometry, dev sources only)",
    }
    (x3.results_dir() / "g1_bulk_verify.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


# ----------------------------------------------------------------------------- held-out ids
HELDOUT_QUERY = (
    "SELECT source_id, ra, dec, phot_g_mean_mag, ruwe, ipd_frac_multi_peak "
    "FROM gaiadr3.gaia_source "
    "WHERE source_id BETWEEN {lo} AND {hi} AND in_andromeda_survey = 'True' "
    "AND MOD(source_id, 10) <> 0 AND phot_g_mean_mag <= {gmax}"
)
TAP_CAP = 50000


def _ids_range(level: int, pix: int, gmax: float) -> list[str]:
    shift = 35 + 2 * (12 - level)
    q = HELDOUT_QUERY.format(lo=pix << shift, hi=((pix + 1) << shift) - 1, gmax=gmax)
    lines = gm._tap_csv(q).strip().splitlines()
    if len(lines) - 1 >= TAP_CAP:  # split into the four children
        return [ln for c in range(4) for ln in _ids_range(level + 1, 4 * pix + c, gmax)]
    return lines


def heldout_ids(gmax: float = 19.0) -> Table:
    """Catalogue attributes of held-out sources (no light curves are read here)."""
    path = x3.out_dir() / f"heldout_ids_g{gmax:g}.csv"
    if not path.exists():
        rows, header = [], None
        for p in x3.gaps_pixels():
            lines = _ids_range(x3.HPX_LEVEL, p, gmax)
            if lines and lines[0].startswith("source_id"):
                header = lines[0]
            rows.extend(ln for ln in lines if not ln.startswith("source_id"))
        gm._write_atomic(path, ("\n".join([header, *rows]) + "\n").encode())
    t = Table.read(path, format="ascii.csv")
    t.meta.update(provenance="observed", source="Gaia archive TAP: " + HELDOUT_QUERY)
    return t


# ----------------------------------------------------------------------------- held-out screen
_CTX: dict = {}


def _init(ctx):
    _CTX.update(ctx)


def _screen_one(args):
    sid, text = args
    pos, veto, scales = _CTX["pos"][sid], _CTX["veto"][sid], _CTX["scales"]
    row = {"kind": "heldout", "event_id": str(sid)}
    try:
        tab = gm.parse_epoch_csv(text)
        row["n_g"] = len(tab)
        if len(tab) < x3.MIN_POINTS:
            row.update({"skipped_min_points": True, "final": False, "error": ""})
            return row, None
        lc = x3.wg.gaia_light_curve(tab, str(sid), pos)
        row.update(g1.screen_full(lc, g1.parse_bands(text), veto, scales))
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001
        row.update({"final": False, "error": repr(exc)[:200]})
    # keep the light curve of anything that reached the fitter flag, for vetting
    return row, (text if row.get("fit_flag") else None)


def screen_heldout(procs: int, files: list[str]) -> None:
    ids = heldout_ids(19.0)
    veto = x3.v1_veto(ids)
    pos = {
        int(s): (float(r), float(d))
        for s, r, d in zip(ids["source_id"], ids["ra"], ids["dec"], strict=True)
    }
    keep = set(pos)
    ctx = {"pos": pos, "veto": {int(k): v for k, v in veto.items()}, "scales": g1.calib_scales()}
    keep_dir = x3.paths.data_root() / "raw" / "gaia_dr3_gaps_heldout_flags"
    keep_dir.mkdir(parents=True, exist_ok=True)
    for name in files:
        out = x3.out_dir() / f"heldout_{name.split('_')[1].split('.')[0]}.ecsv"
        if out.exists():
            print(f"{name}: done already")
            continue
        rows, gen = [], stream(name, keep)
        stats = {}

        def feed(gen=gen):
            nonlocal stats
            while True:
                try:
                    yield next(gen)
                except StopIteration as stop:
                    stats = stop.value or {}
                    return

        with Pool(procs, initializer=_init, initargs=(ctx,)) as pool:
            for row, text in pool.imap_unordered(_screen_one, feed(), 64):
                rows.append(row)
                if text is not None:
                    (keep_dir / f"{row['event_id']}.csv").write_text(text)
        g1.write_rows(rows, out, meta={"file": name, **stats, "provenance": "observed (held-out)"})
        print(f"{name}: {stats} screened {len(rows)}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("cmd", choices=("verify", "heldout-ids", "screen"))
    ap.add_argument("--file", default=FILES[0][0])
    ap.add_argument("--procs", type=int, default=3)
    a = ap.parse_args(argv)
    if a.cmd == "verify":
        print(json.dumps(verify(a.file), indent=1))
    elif a.cmd == "heldout-ids":
        t = heldout_ids(19.0)
        print(f"{len(t)} held-out sources; {sum(x3.v1_veto(t).values())} v1-vetoed")
    else:
        screen_heldout(a.procs, [n for n, _, _ in FILES])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
