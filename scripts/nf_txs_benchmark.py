"""Neutrino Frontier benchmark 1: TXS 0506+056 2014-15 box flare with SkyLLH on IceTracks-DR2.

Known case before any SkyLLH use on new questions (docs/neutrino_frontier/data_audit.md,
"Benchmarks").
IceTracks-DR2 (arXiv:2605.19040) section 5 fixes the box to the Table 6 best fit (T0 = MJD 57020,
dT = 185 d, IC86_IV only) and fits n_s and gamma with SkyLLH: n_s = 12.7, gamma = 2.3 (internal
IceCube tools: 12.56 / 2.26). Table 6 rounds T0 and dT to whole days, and a box flare fit puts its
edges on event times, so the script fits two windows:

- ``rounded``: [T0 - dT/2, T0 + dT/2] as printed;
- ``event_edges``: each rounded edge snapped to the nearest on-source event (within ``SNAP_DEG`` of
  the source) no more than ``SNAP_DAYS`` away (ASSUMPTIONs; rounding moves an edge by <= 0.75 d).

1. ``fetch`` streams the IC86_IV events and uptime plus the IC86 IRFs (smearing matrix 817 MB, v2;
   598 MB, v1) from Harvard Dataverse into the data root and checks each against
   ``data/manifests/nf_skyllh_dr2.ecsv`` (sha256 pin; Dataverse's md5 is listed too).
2. ``run`` fits both windows under each IRF version and, with ``--trials N``, estimates the
   fixed-window background TS distribution (SkyLLH's scrambled background, empirical p).
3. ``--cleanup`` deletes the raw files afterwards (cloud-disk decision).

Writes ``results/nf/txs_skyllh_benchmark.json``.

  python scripts/nf_txs_benchmark.py [--irfs v2 v1] [--trials 50000] [--cpu 4] [--cleanup]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from pathlib import Path

if __name__ == "__main__":  # SkyLLH trial pool: one thread per process (set before numpy loads)
    os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from jwst_anomaly.acquire import sha256_file  # noqa: E402
from jwst_anomaly.burst_twins import angsep_deg  # noqa: E402
from jwst_anomaly.neutrino_frontier import DR2_COLUMNS  # noqa: E402

OUT = ROOT / "results" / "nf" / "txs_skyllh_benchmark.json"
MANIFEST = ROOT / "data" / "manifests" / "nf_skyllh_dr2.ecsv"
URL = "https://dataverse.harvard.edu/api/access/datafile/{}?format=original"
#: SkyLLH sample name and sub-directory per IRF version (Dataverse 1.0 vs >= 2.0 IRF binning)
SAMPLES = {
    "v1": ("IceTracks-DR2-v1", "icecube_14year_ps"),
    "v2": ("IceTracks-DR2", "icecube_14year_ps_v2"),
}
SEASON = "IC86_IV"
#: Table 6 / section 5 of arXiv:2605.19040 (the paper's catalogue position, not the VLBI one
#: nf_ghost_dr2 uses; rounded window) and the targets
TXS = {"ra": 77.35, "dec": 5.7, "t0": 57020.0, "dt": 185.0}
PAPER = {
    "skyllh": {"ns": 12.7, "gamma": 2.3},
    "internal": {"ns": 12.56, "gamma": 2.26, "p_pre": 4.18e-3},
}
#: ASSUMPTIONs: edge snapping (on-source radius, maximal shift) and the agreement tolerance
SNAP_DEG = 1.0
SNAP_DAYS = 1.0
TOL = {"ns": 1.0, "gamma": 0.1}


def rounded_window(t0: float, dt: float) -> tuple[float, float]:
    return t0 - dt / 2, t0 + dt / 2


def snap_edges(
    mjd: np.ndarray,
    sep_deg: np.ndarray,
    window: tuple[float, float],
    snap_deg: float = SNAP_DEG,
    snap_days: float = SNAP_DAYS,
) -> tuple[tuple[float, float], tuple[bool, bool]]:
    """Move each window edge to the nearest on-source event within ``snap_days``; an edge with no
    such event stays where it is. Returns the edges and which of them were snapped (snapped
    edges land exactly on event times; ``fit_windows`` pads them by 1e-3 d so the edge events
    are inside the box)."""
    on = mjd[sep_deg < snap_deg]
    edges, snapped = [], []
    for edge in window:
        d = np.abs(on - edge)
        hit = bool(d.size and d.min() <= snap_days)
        edges.append(float(on[np.argmin(d)]) if hit else float(edge))
        snapped.append(hit)
    if edges[0] >= edges[1]:
        raise ValueError(f"snapped window {edges} is empty")
    return (edges[0], edges[1]), (snapped[0], snapped[1])


def read_events(path: Path) -> np.ndarray:
    """Original DR2 events CSV as an array in ``DR2_COLUMNS`` order (header checked)."""
    with open(path) as f:
        head = f.readline().lstrip("#").split()
    names = [h.split("[")[0].replace("log10(E/GeV)", "log10e").lower() for h in head]
    if names != list(DR2_COLUMNS):
        raise SystemExit(f"{path}: unexpected header {head}")
    return np.loadtxt(path)


def agrees(fit: dict, target: dict, tol: dict = TOL) -> bool:
    return all(abs(fit[k] - target[k]) <= tol[k] for k in tol)


def _sha256_stream(url: str, dest: Path, tries: int = 5) -> str:
    """Stream ``url`` to ``dest`` with back-off on 429/5xx and on dropped or stalled transfers
    (scripts/CLAUDE.md); returns the sha256. A failed attempt leaves no partial file."""
    import requests

    tmp = dest.with_suffix(dest.suffix + ".part")
    retry = (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError)
    for k in range(tries):
        try:
            with requests.get(url, stream=True, timeout=900) as r:
                if r.status_code != 429 and r.status_code < 500:
                    r.raise_for_status()
                    h = hashlib.sha256()
                    with open(tmp, "wb") as f:
                        for chunk in r.iter_content(1 << 22):
                            h.update(chunk)
                            f.write(chunk)
                    tmp.replace(dest)
                    return h.hexdigest()
        except retry:
            if k == tries - 1:
                raise
        finally:
            tmp.unlink(missing_ok=True)  # any failure (incl. a full disk) leaves no partial file
        time.sleep(2 ** (k + 1))
    raise SystemExit(f"{url}: still failing after {tries} tries")


def fetch(base: Path, versions: list[str]) -> None:
    """Download the manifest's files for ``versions``; a file whose sha256 matches is kept."""
    from astropy.table import Table

    man = Table.read(MANIFEST, format="ascii.ecsv")
    for row in man:
        if row["irf_version"] not in versions:
            continue
        dest = base / SAMPLES[row["irf_version"]][1] / row["file"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and sha256_file(dest) == row["sha256"]:
            continue
        t = time.time()
        got = _sha256_stream(URL.format(int(row["dataverse_file_id"])), dest)
        if got != row["sha256"]:
            dest.unlink()
            what = f"{row['file']} ({row['irf_version']})"
            raise SystemExit(f"{what}: sha256 {got} != pinned {row['sha256']}")
        mb = row["size"] / 1e6
        print(
            f"{row['irf_version']} {row['file']}: {mb:.1f} MB, {time.time() - t:.0f} s", flush=True
        )


def fit_windows(base: Path, version: str, trials: int, cpu: int) -> dict:
    from skyllh.analyses.i3.publicdata_ps.time_dependent_ps import (
        change_signal_time_pdf_of_llhratio_function,
        create_analysis,
    )
    from skyllh.core.config import Config
    from skyllh.core.random import RandomStateService
    from skyllh.core.source_model import PointLikeSource
    from skyllh.datasets import create_datasets

    cfg = Config()
    ds = create_datasets(SAMPLES[version][0], cfg=cfg, base_path=str(base.resolve()), names=SEASON)
    ds = ds if isinstance(ds, list) else [ds]
    ev = read_events(base / SAMPLES[version][1] / "events" / f"{SEASON}_exp.csv")
    sep = angsep_deg(ev[:, 6], ev[:, 7], TXS["ra"], TXS["dec"])
    rounded = rounded_window(TXS["t0"], TXS["dt"])
    (lo, hi), snapped = snap_edges(ev[:, 3], sep, rounded)
    on = sep < SNAP_DEG
    pad = 1e-3
    windows = {"rounded": rounded, "event_edges": (lo - pad * snapped[0], hi + pad * snapped[1])}
    src = PointLikeSource(ra=np.deg2rad(TXS["ra"]), dec=np.deg2rad(TXS["dec"]))
    ana = create_analysis(
        cfg=cfg, datasets=ds, source=src, box={"start": rounded[0], "stop": rounded[1]}
    )
    out = {}
    for name, (start, stop) in windows.items():
        change_signal_time_pdf_of_llhratio_function(ana, box={"start": start, "stop": stop})
        ts, par, _ = ana.unblind(minimizer_rss=RandomStateService(1))
        fit = {"ns": round(float(par["ns"]), 3), "gamma": round(float(par["gamma"]), 3)}
        w = {
            "start_mjd": round(start, 4),
            "stop_mjd": round(stop, 4),
            "t0_mjd": round((start + stop) / 2, 3),
            "dt_days": round(stop - start, 3),
            "edges_snapped": list(snapped) if name == "event_edges" else [False, False],
            "n_on_source_in_box": int(((ev[:, 3] >= start) & (ev[:, 3] <= stop) & on).sum()),
            "ts": round(float(ts), 3),
            **fit,
            "agrees_with_paper_skyllh": agrees(fit, PAPER["skyllh"]),
            # approximation, not a measurement: half the background trials have TS = 0
            "p_wilks_half_chi2_2dof": float(0.5 * stats.chi2.sf(ts, 2)),
        }
        if trials:
            t = time.time()
            bg = ana.do_trials(rss=RandomStateService(20261010), n=trials, ncpu=cpu, mean_n_sig=0)
            bts = np.asarray(bg["ts"], dtype=float)
            k = int((bts >= ts).sum())
            w["bkg_trials"] = {
                "n": trials,
                "n_ge_obs": k,
                "p_empirical": (k + 1) / (trials + 1),  # an upper bound when n_ge_obs = 0
                "frac_ts_zero": round(float((bts <= 0).mean()), 4),
                "ts_median": round(float(np.median(bts)), 3),
                "seconds": round(time.time() - t, 1),
            }
        out[name] = w
        print(version, name, w, flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--irfs", nargs="+", choices=sorted(SAMPLES), default=["v2", "v1"])
    ap.add_argument("--trials", type=int, default=0, help="background trials per window (v2 only)")
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1, help="trial processes")
    ap.add_argument("--cleanup", action="store_true", help="delete the raw files afterwards")
    a = ap.parse_args(argv)
    from jwst_anomaly.paths import data_root

    base = data_root() / "nf_skyllh"
    t0 = time.time()
    fits = {}
    for v in a.irfs:  # one IRF version on disk at a time when cleaning up
        sub = base / SAMPLES[v][1]
        try:
            fetch(base, [v])
            fits[v] = fit_windows(base, v, a.trials if v == "v2" else 0, a.cpu)
        finally:
            if a.cleanup:
                shutil.rmtree(sub, ignore_errors=True)
    from importlib.metadata import version

    out = {
        "test": (
            "Neutrino Frontier benchmark 1: TXS 0506+056 2014-15 box flare, SkyLLH on IceTracks-DR2"
        ),
        "provenance": "derived (observed IceTracks-DR2 IC86_IV; IRFs are model_prediction)",
        "source": {"ra_deg": TXS["ra"], "dec_deg": TXS["dec"], "season": SEASON},
        "paper": {
            "ref": "arXiv:2605.19040 Table 6 and section 5",
            **PAPER,
            "t0": 57020.0,
            "dt": 185.0,
        },
        "assumptions": {"snap_deg": SNAP_DEG, "snap_days": SNAP_DAYS, "agreement_tolerance": TOL},
        "skyllh_version": version("skyllh"),
        "fits": fits,
        "seconds": round(time.time() - t0, 1),
    }
    if OUT.exists():  # a partial run (one IRF version) keeps the other version's fits
        out["fits"] = {**json.loads(OUT.read_text()).get("fits", {}), **fits}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
