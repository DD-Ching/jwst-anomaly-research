"""S1 burst twins, stage 1: stream Fermi GBM bcat files and reduce them to compact three-band light curves.

Writes one tracked table per trigger year, ``results/s1_twins/lc_<YYYY>.ecsv.gz`` (< 1 MB each), plus the
catalogue subset ``results/s1_twins/catalogue.ecsv.gz``. Raw FITS files are held in memory only and never stored.
Resumable: a year whose chunk exists is skipped.

  python scripts/s1_ingest.py [--years 2019 2020] [--io 12] [--cpu 4]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import hashlib
import io
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "s1_twins"
TAP = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
FTP = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/bursts"
CAT_COLS = (
    "trigger_name, ra, dec, error_radius, trigger_time, t90, t90_error, t90_start, fluence, fluence_error, "
    "flux_64, flux_256, flnc_best_fitting_model, flnc_comp_epeak, flnc_band_epeak, flnc_comp_ergflnc, "
    "bcatalog, last_modified"
)


def write_ecsv_gz(t: Table, path: Path) -> None:
    """astropy does not compress on write; gzip the ECSV text explicitly (mtime 0 for reproducible bytes)."""
    buf = io.StringIO()
    t.write(buf, format="ascii.ecsv")
    with open(path, "wb") as fh, gzip.GzipFile(fileobj=fh, mode="wb", mtime=0) as gz:
        gz.write(buf.getvalue().encode())


def _init_worker():
    os.environ["OMP_NUM_THREADS"] = "1"


def get(session: requests.Session, url: str, tries: int = 6) -> requests.Response:
    """GET with exponential back-off on 429/5xx and connection errors."""
    delay = 2.0
    for i in range(tries):
        try:
            r = session.get(url, timeout=120)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                r.raise_for_status()
            if r.status_code not in (429, 500, 502, 503, 504):
                r.raise_for_status()
        except (requests.ConnectionError, requests.Timeout):
            if i == tries - 1:
                raise
        time.sleep(delay)
        delay *= 2
    raise RuntimeError(f"giving up on {url}")


def fetch_catalogue() -> Table:
    from jwst_anomaly.schema import Provenance

    r = requests.get(
        TAP,
        params={"REQUEST": "doQuery", "LANG": "ADQL", "QUERY": f"SELECT {CAT_COLS} FROM fermigbrst"},
        timeout=300,
    )
    r.raise_for_status()
    t = Table.read(io.BytesIO(r.content), format="votable")
    for c in t.colnames:
        if t[c].dtype.kind == "O":
            t[c] = t[c].astype(str)
    t.sort("trigger_name")
    t.meta = {
        "provenance": str(Provenance.OBSERVED),
        "source": f"HEASARC TAP {TAP} fermigbrst ({len(t)} rows), fetched "
        f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}; response sha256 {hashlib.sha256(r.content).hexdigest()}",
    }
    return t


def download(session, name: str):
    """List the burst's 'current' directory and download its bcat file into memory."""
    year = "20" + name[2:4]
    base = f"{FTP}/{year}/{name}/current/"
    try:
        listing = get(session, base).text
    except Exception as exc:  # noqa: BLE001 - recorded per burst
        return name, None, None, f"listing: {exc}"[:120]
    files = sorted(set(re.findall(r'href="(glg_bcat_all_' + name + r'_v\d+\.fit)"', listing)))
    if not files:
        return name, None, None, "no bcat file"
    fn = files[-1]
    try:
        content = get(session, base + fn).content
    except Exception as exc:  # noqa: BLE001
        return name, fn, None, f"download: {exc}"[:120]
    return name, fn, content, ""


def reduce_one(args):
    """Worker: bytes -> encoded row dict (runs in a process pool)."""
    from astropy.io import fits

    from jwst_anomaly import burst_twins as bt

    name, fn, content, t90_start, t90 = args
    sha = hashlib.sha256(content).hexdigest()
    try:
        with fits.open(io.BytesIO(content), memmap=False) as h:
            red = bt.reduce_bcat(h, t90_start, t90)
        units, fs, es = bt.encode(red["flux"], red["err"])
        return {
            "trigger_name": name,
            "bcat_file": fn,
            "sha256": sha,
            "size": len(content),
            "dets": red["dets"],
            "t_lo": round(red["t_lo"], 4),
            "dt": red["dt"],
            "nbin": red["flux"].shape[1],
            "u1": units[0],
            "u2": units[1],
            "u3": units[2],
            "f1": fs[0],
            "f2": fs[1],
            "f3": fs[2],
            "e1": es[0],
            "e2": es[1],
            "e3": es[2],
            "status": "ok",
        }
    except Exception as exc:  # noqa: BLE001
        return {"trigger_name": name, "bcat_file": fn, "sha256": sha, "size": len(content), "status": str(exc)[:120]}


def empty_row(name, fn, status):
    return {"trigger_name": name, "bcat_file": fn or "", "sha256": "", "size": 0, "status": status}


def run_year(year: str, cat: Table, n_io: int, pool: cf.ProcessPoolExecutor) -> None:
    from jwst_anomaly.schema import Provenance

    sel = cat[[str(n)[2:4] == year[2:] for n in cat["trigger_name"]]]
    t0, cpu0 = time.time(), os.times()
    rows, nbytes = [], 0
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(pool_connections=n_io, pool_maxsize=n_io)
    session.mount("https://", adapter)
    info = {str(r["trigger_name"]): (float(r["t90_start"]), float(r["t90"])) for r in sel}
    futs = []
    with cf.ThreadPoolExecutor(n_io) as tp:
        for name, fn, content, err in tp.map(lambda n: download(session, n), list(info)):
            if content is None:
                rows.append(empty_row(name, fn, err))
                continue
            nbytes += len(content)
            futs.append(pool.submit(reduce_one, (name, fn, content, *info[name])))
        for f in cf.as_completed(futs):
            rows.append(f.result())
    dt = time.time() - t0
    cpu1 = os.times()
    cpu = (cpu1.children_user + cpu1.children_system - cpu0.children_user - cpu0.children_system) + (
        cpu1.user + cpu1.system - cpu0.user - cpu0.system
    )
    keys = [
        "trigger_name",
        "bcat_file",
        "sha256",
        "size",
        "status",
        "dets",
        "t_lo",
        "dt",
        "nbin",
        "u1",
        "u2",
        "u3",
        "f1",
        "f2",
        "f3",
        "e1",
        "e2",
        "e3",
    ]
    defaults = {"dets": "", "t_lo": np.nan, "dt": np.nan, "nbin": 0, "u1": np.nan, "u2": np.nan, "u3": np.nan}
    cols = {k: [r.get(k, defaults.get(k, "")) for r in rows] for k in keys}
    t = Table(cols)
    t.sort("trigger_name")
    t.meta = {
        "provenance": str(Provenance.DERIVED),
        "source": f"Fermi GBM bcat files ({FTP}/<year>/<bn>/current/glg_bcat_all_*.fit; sha256 per row); "
        "INCLUDED NaI detectors inverse-variance combined; bands keV 10-50/50-300/300-1000 (CTIME channels "
        "1-2/3-4/5-6); flux f<b> and error e<b> are integer tenths of unit u<b> (ph cm^-2 s^-1); bin i covers "
        "[t_lo + i dt, t_lo + (i+1) dt) s from trigger; scripts/s1_ingest.py",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"lc_{year}.ecsv.gz"
    write_ecsv_gz(t, path)
    ok = sum(1 for r in rows if r.get("status") == "ok")
    ncores = os.cpu_count() or 1
    print(
        f"[{year}] {len(rows)} bursts ({ok} ok) in {dt:.0f} s: {len(rows) / dt:.2f} items/s, "
        f"{nbytes / dt / 1e6:.1f} MB/s, CPU {100 * cpu / dt / ncores:.0f} % of {ncores} cores, "
        f"chunk {path.stat().st_size / 1e3:.0f} kB",
        flush=True,
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--years", nargs="*")
    ap.add_argument("--io", type=int, default=12)
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args(argv)
    sys.path.insert(0, str(ROOT / "src"))
    OUT.mkdir(parents=True, exist_ok=True)
    cat_path = OUT / "catalogue.ecsv.gz"
    if cat_path.exists():
        cat = Table.read(cat_path, format="ascii.ecsv")
    else:
        cat = fetch_catalogue()
        write_ecsv_gz(cat, cat_path)
        print(f"catalogue: {len(cat)} bursts -> {cat_path}", flush=True)
    years = a.years or sorted({"20" + str(n)[2:4] for n in cat["trigger_name"]})
    with cf.ProcessPoolExecutor(a.cpu, initializer=_init_worker) as pool:
        for y in years:
            if (OUT / f"lc_{y}.ecsv.gz").exists():
                print(f"[{y}] exists, skipped", flush=True)
                continue
            run_year(y, cat, a.io, pool)


if __name__ == "__main__":
    main()
