"""Neutrino Frontier benchmark 2: time-integrated SkyLLH fits on all IceTracks-DR2 seasons.

Benchmark 1 (``nf_txs_benchmark.py``, D-080) used one season (IC86_IV). This one uses every season
(IC40, IC59, IC79, IC86_I-XI), so it also checks the IC40/IC59/IC79 IRFs and SkyLLH's multi-dataset
likelihood. Targets: arXiv:2605.19040 Table 8 ("DR2" columns, SkyLLH on IceTracks-DR2) for the four
most significant catalogue sources, at the Table 4 positions:
NGC 1068 80.1 / 3.2, PKS 1424+240 75.6 / 3.6, TXS 0506+056 8.8 / 2.0, GB6 J1542+6129 49.6 / 3.8
(n_s / gamma). Agreement: |dn_s| <= max(1, 0.05 n_s) and |dgamma| <= 0.1 (ASSUMPTIONs; Table 8
rounds gamma to 0.1, so 0.1 is the rounding floor).

1. ``fetch`` lists the Dataverse dataset version (1.0: the IRFs the paper used; 3.1: binning fix),
   streams every events / uptime / IRF file (4 smearing matrices, 598 MB each in 1.0), checks
   Dataverse's md5 and pins sha256 in ``data/manifests/nf_skyllh_dr2_all.ecsv``; a later run checks
   the pin.
2. ``run`` fits n_s and gamma for each source (SkyLLH time-integrated point-source analysis).
3. ``--cleanup`` deletes the raw files afterwards (cloud-disk decision).

Writes ``results/nf/ti_skyllh_benchmark.json``.

  python scripts/nf_ti_benchmark.py [--irfs v1 v2] [--sources "NGC 1068"] [--cleanup]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from nf_txs_benchmark import HEADERS, SAMPLES, URL, _sha256_stream  # noqa: E402

from jwst_anomaly.acquire import sha256_file  # noqa: E402

OUT = ROOT / "results" / "nf" / "ti_skyllh_benchmark.json"
MANIFEST = ROOT / "data" / "manifests" / "nf_skyllh_dr2_all.ecsv"
DOI = "doi:10.7910/DVN/MMIIZA"
LIST_URL = (
    "https://dataverse.harvard.edu/api/datasets/:persistentId/versions/{}/files?persistentId=" + DOI
)
#: Dataverse dataset version per IRF version (see data/manifests/nf_skyllh_dr2.ecsv)
DV_VERSION = {"v1": "1.0", "v2": "3.1"}
#: Table 4 positions (deg) and Table 8 DR2 fits of arXiv:2605.19040
SOURCES = {
    "NGC 1068": {"ra": 40.67, "dec": -0.01, "ns": 80.1, "gamma": 3.2, "p_pre": 1.3e-7},
    "PKS 1424+240": {"ra": 216.76, "dec": 23.8, "ns": 75.6, "gamma": 3.6, "p_pre": 8.0e-5},
    "TXS 0506+056": {"ra": 77.35, "dec": 5.7, "ns": 8.8, "gamma": 2.0, "p_pre": 5.0e-3},
    "GB6 J1542+6129": {"ra": 235.75, "dec": 61.5, "ns": 49.6, "gamma": 3.8, "p_pre": 2.0e-3},
}
GAMMA_TOL = 0.1


def agrees(fit: dict, target: dict) -> bool:
    """ASSUMPTION: |dn_s| <= max(1, 0.05 n_s) and |dgamma| <= 0.1."""
    ok_ns = abs(fit["ns"] - target["ns"]) <= max(1.0, 0.05 * target["ns"])
    return ok_ns and abs(fit["gamma"] - target["gamma"]) <= GAMMA_TOL


def dataverse_files(dv_version: str) -> list[dict]:
    """Events, uptime and IRF CSVs of one Dataverse dataset version (plots and readme skipped)."""
    import requests

    r = requests.get(LIST_URL.format(dv_version), timeout=120, headers=HEADERS)
    r.raise_for_status()
    rows = []
    for f in r.json()["data"]:
        df, folder = f["dataFile"], f.get("directoryLabel", "")
        name = df.get("originalFileName", df["filename"])
        if folder in ("events", "uptime", "irfs") and name.endswith(".csv"):
            size = df.get("originalFileSize", df["filesize"])
            rows.append(
                {"id": df["id"], "file": f"{folder}/{name}", "md5": df["md5"], "size": size}
            )
    return rows


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(base: Path, version: str) -> None:
    """Stream the files of ``version``; keep a file whose md5 matches; check / extend the pin."""
    from astropy.table import Table

    pins = {}
    if MANIFEST.exists():
        man = Table.read(MANIFEST, format="ascii.ecsv")
        pins = {(r["irf_version"], r["file"]): r["sha256"] for r in man}
    new = []
    for row in dataverse_files(DV_VERSION[version]):
        dest = base / SAMPLES[version][1] / row["file"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        t = time.time()
        sha = None
        if not (dest.exists() and _md5(dest) == row["md5"]):
            sha = _sha256_stream(URL.format(row["id"]), dest)
            if _md5(dest) != row["md5"]:
                dest.unlink()
                raise SystemExit(f"{row['file']} ({version}): md5 differs from Dataverse's")
            print(f"{version} {row['file']}: {row['size'] / 1e6:.1f} MB, {time.time() - t:.0f} s")
        sha = sha or sha256_file(dest)
        pinned = pins.get((version, row["file"]))
        if pinned is not None and pinned != sha:
            raise SystemExit(f"{row['file']} ({version}): sha256 {sha} != pinned {pinned}")
        if pinned is None:
            new.append(
                {
                    "irf_version": version,
                    "dataverse_file_id": int(row["id"]),
                    "file": row["file"],
                    "md5": row["md5"],
                    "sha256": sha,
                    "size": int(row["size"]),
                    "retrieved_utc": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime()),
                    "dataverse_version": DV_VERSION[version],
                }
            )
    if new:
        write_manifest(new)


def write_manifest(rows: list[dict]) -> None:
    from astropy.table import Table, vstack

    t = Table(rows=rows)
    if MANIFEST.exists():
        t = vstack([Table.read(MANIFEST, format="ascii.ecsv"), t])
    t.meta = {
        "provenance": "observed",
        "source": f"IceCube IceTracks-DR2, Harvard Dataverse {DOI} (CC0); "
        "v1 = dataset version 1.0, v2 = 3.1; md5 as listed by Dataverse, "
        "sha256 of the downloaded original file (the pin)",
        "note": "scripts/nf_ti_benchmark.py; raw files are deleted after use (cloud-disk "
        "decision); 4 smearing matrices of 0.6-0.8 GB each: the stated reason is benchmark 2 "
        "(TASKS Neutrino Frontier item 2), which needs every season's IRFs",
    }
    t.sort(["irf_version", "file"])
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    t.write(MANIFEST, format="ascii.ecsv", overwrite=True)


def fit_sources(base: Path, version: str, names: list[str]) -> dict:
    from skyllh.analyses.i3.publicdata_ps.time_integrated_ps import create_analysis
    from skyllh.core.config import Config
    from skyllh.core.random import RandomStateService
    from skyllh.core.source_model import PointLikeSource
    from skyllh.datasets import create_datasets

    cfg = Config()
    ds = create_datasets(SAMPLES[version][0], cfg=cfg, base_path=str(base.resolve()))
    out = {}
    for name in names:
        s = SOURCES[name]
        t = time.time()
        src = PointLikeSource(ra=np.deg2rad(s["ra"]), dec=np.deg2rad(s["dec"]))
        ana = create_analysis(cfg=cfg, datasets=ds, source=src)
        ts, par, _ = ana.unblind(minimizer_rss=RandomStateService(1))
        fit = {"ns": round(float(par["ns"]), 3), "gamma": round(float(par["gamma"]), 3)}
        out[name] = {
            "datasets": [d.name for d in ds],
            "ts": round(float(ts), 3),
            **fit,
            "agrees_with_paper_dr2": agrees(fit, s),
            # approximation, not a measurement (the paper's p-values come from trials)
            "p_wilks_half_chi2_2dof": float(0.5 * stats.chi2.sf(ts, 2)),
            "seconds": round(time.time() - t, 1),
        }
        print(version, name, out[name], flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--irfs", nargs="+", choices=sorted(SAMPLES), default=["v1"])
    ap.add_argument("--sources", nargs="+", choices=sorted(SOURCES), default=list(SOURCES))
    ap.add_argument("--cleanup", action="store_true", help="delete the raw files afterwards")
    a = ap.parse_args(argv)
    from jwst_anomaly.paths import data_root

    base = data_root() / "nf_skyllh"
    t0 = time.time()
    fits = {}
    for v in a.irfs:  # one IRF version on disk at a time when cleaning up
        try:
            fetch(base, v)
            fits[v] = fit_sources(base, v, a.sources)
        finally:
            if a.cleanup:
                shutil.rmtree(base / SAMPLES[v][1], ignore_errors=True)
    from importlib.metadata import version

    out = {
        "test": (
            "Neutrino Frontier benchmark 2: time-integrated SkyLLH fits, all IceTracks-DR2 seasons"
        ),
        "provenance": "derived (observed IceTracks-DR2 events; IRFs are model_prediction)",
        "paper": {
            "ref": "arXiv:2605.19040 Table 4 (positions) and Table 8 (DR2 fits)",
            "sources": SOURCES,
        },
        "assumptions": {"ns_tolerance": "max(1, 0.05 n_s)", "gamma_tolerance": GAMMA_TOL},
        "skyllh_version": version("skyllh"),
        "fits": fits,
        "seconds": round(time.time() - t0, 1),
    }
    if OUT.exists():  # a partial run keeps the other fits
        old = json.loads(OUT.read_text()).get("fits", {})
        out["fits"] = {v: {**old.get(v, {}), **fits.get(v, {})} for v in {*old, *fits}}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
