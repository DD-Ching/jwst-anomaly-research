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
   ``data/manifests/nf_skyllh_dr2.ecsv`` (Dataverse md5 of the original CSV).
2. ``run`` fits both windows under each IRF version and, with ``--trials N``, estimates the
   fixed-window background TS distribution (SkyLLH's scrambled background, empirical p).
3. ``--cleanup`` deletes the raw files afterwards (cloud-disk decision).

Writes ``results/nf/txs_skyllh_benchmark.json``.

  python scripts/nf_txs_benchmark.py [--irfs v2 v1] [--trials 2000] [--cpu 4] [--cleanup]
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

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

OUT = ROOT / "results" / "nf" / "txs_skyllh_benchmark.json"
MANIFEST = ROOT / "data" / "manifests" / "nf_skyllh_dr2.ecsv"
URL = "https://dataverse.harvard.edu/api/access/datafile/{}?format=original"
#: SkyLLH sample name and sub-directory per IRF version (Dataverse 1.0 vs >= 2.0 IRF binning)
SAMPLES = {
    "v1": ("IceTracks-DR2-v1", "icecube_14year_ps"),
    "v2": ("IceTracks-DR2", "icecube_14year_ps_v2"),
}
SEASON = "IC86_IV"
#: Table 6 / section 5 of arXiv:2605.19040 (catalogue position, rounded window) and the targets
TXS = {"ra": 77.35, "dec": 5.7, "t0": 57020.0, "dt": 185.0}
PAPER = {
    "skyllh": {"ns": 12.7, "gamma": 2.3},
    "internal": {"ns": 12.56, "gamma": 2.26, "p_pre": 4.18e-3},
}
#: ASSUMPTIONs: edge snapping (on-source radius, maximal shift) and the agreement tolerance
SNAP_DEG = 1.0
SNAP_DAYS = 1.0
TOL = {"ns": 1.0, "gamma": 0.1}


def angsep_deg(ra1, dec1, ra2, dec2):
    """Great-circle separation in degrees (inputs in degrees)."""
    r1, d1, r2, d2 = map(np.deg2rad, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.rad2deg(np.arccos(np.clip(c, -1.0, 1.0)))


def rounded_window(t0: float, dt: float) -> tuple[float, float]:
    return t0 - dt / 2, t0 + dt / 2


def snap_edges(
    mjd: np.ndarray,
    sep_deg: np.ndarray,
    window: tuple[float, float],
    snap_deg: float = SNAP_DEG,
    snap_days: float = SNAP_DAYS,
) -> tuple[float, float]:
    """Move each window edge to the nearest on-source event within ``snap_days``; an edge with no
    such event stays where it is. Edges land exactly on event times (``run`` pads them by 1e-3 d
    so the edge events are inside the box)."""
    on = mjd[sep_deg < snap_deg]
    out = []
    for edge in window:
        d = np.abs(on - edge)
        out.append(float(on[np.argmin(d)]) if d.size and d.min() <= snap_days else float(edge))
    if out[0] >= out[1]:
        raise ValueError(f"snapped window {out} is empty")
    return out[0], out[1]


def agrees(fit: dict, target: dict, tol: dict = TOL) -> bool:
    return all(abs(fit[k] - target[k]) <= tol[k] for k in tol)


def _md5_stream(url: str, dest: Path, tries: int = 5) -> str:
    """Stream ``url`` to ``dest`` with back-off on 429/5xx (scripts/CLAUDE.md); returns the md5."""
    import requests

    for k in range(tries):
        try:
            with requests.get(url, stream=True, timeout=900) as r:
                if r.status_code != 429 and r.status_code < 500:
                    r.raise_for_status()
                    h = hashlib.md5()
                    tmp = dest.with_suffix(dest.suffix + ".part")
                    with open(tmp, "wb") as f:
                        for chunk in r.iter_content(1 << 22):
                            h.update(chunk)
                            f.write(chunk)
                    tmp.replace(dest)
                    return h.hexdigest()
        except requests.ConnectionError:
            if k == tries - 1:
                raise
        time.sleep(2 ** (k + 1))
    raise SystemExit(f"{url}: still failing after {tries} tries")


def fetch(base: Path, versions: list[str]) -> None:
    """Download the manifest's files for ``versions`` (skipping files whose md5 already matches)."""
    from astropy.table import Table

    man = Table.read(MANIFEST, format="ascii.ecsv")
    for row in man:
        if row["irf_version"] not in versions:
            continue
        dest = base / SAMPLES[row["irf_version"]][1] / row["file"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and hashlib.md5(dest.read_bytes()).hexdigest() == row["md5"]:
            continue
        t = time.time()
        got = _md5_stream(URL.format(int(row["dataverse_file_id"])), dest)
        if got != row["md5"]:
            dest.unlink()
            what = f"{row['file']} ({row['irf_version']})"
            raise SystemExit(f"{what}: md5 {got} != pinned {row['md5']}")
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
    ev = np.loadtxt(base / SAMPLES[version][1] / "events" / f"{SEASON}_exp.csv")
    sep = angsep_deg(ev[:, 6], ev[:, 7], TXS["ra"], TXS["dec"])
    rounded = rounded_window(TXS["t0"], TXS["dt"])
    snapped = snap_edges(ev[:, 3], sep, rounded)
    on = sep < SNAP_DEG
    windows = {"rounded": rounded, "event_edges": (snapped[0] - 1e-3, snapped[1] + 1e-3)}
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
    ap.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--cleanup", action="store_true", help="delete the raw files afterwards")
    a = ap.parse_args(argv)
    from jwst_anomaly.paths import data_root

    base = data_root() / "nf_skyllh"
    t0 = time.time()
    fetch(base, a.irfs)
    fits = {v: fit_windows(base, v, a.trials if v == "v2" else 0, a.cpu) for v in a.irfs}
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
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    if a.cleanup:
        shutil.rmtree(base)
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
