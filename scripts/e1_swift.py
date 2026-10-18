"""E1 + Swift BAT: Swift BAT GRBs as a fifth node set in the E1 lag x separation pair counts
(hypothesis; owner idea 4, D-074). Frozen plan: ``results/e1_swift/PLAN.md``.

Tested family: BAT-BAT, BAT-ICECAT, BAT-CHIME (same, wide) and BAT-GW (all) over the five D-074 lag
bins. BAT-GBM is the positive control (the same GRBs), not in the family.

Inputs (~5 MB) go to ``$JWST_ANOMALY_DATA/e1_swift/`` and are pinned by sha256 in
``data/manifests/e1_swift.ecsv`` (own pins: the live GBM catalogue has moved since D-074's pin).

  python scripts/e1_swift.py [--jit 10000] [--jitday 10000] [--inject 40] [--cpu 4] [--refresh]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import e1_events as e1  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402
from jwst_anomaly.paths import data_root  # noqa: E402
from jwst_anomaly.schema import Provenance  # noqa: E402

OUT = ROOT / "results" / "e1_swift"
MANIFEST = ROOT / "data" / "manifests" / "e1_swift.ecsv"
SOURCES = {
    "swiftgrbba.vot": "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync?REQUEST=doQuery&LANG=ADQL&QUERY="
    "SELECT%20name,time,trigger_obs,ra,dec,error_radius%20FROM%20swiftgrbba",
    **e1.SOURCES,
}
#: Tested channels (frozen plan), then the positive-control channel.
TESTED = (("BAT", "BAT"), ("BAT", "ICECAT"), ("BAT", "CHIME"), ("BAT", "GW"))
CONTROL = ("BAT", "GBM")
CHANNELS = TESTED + (CONTROL,)
#: ASSUMPTION: Swift orbital period (s), ~95 min; it decays over the mission, so k x period keeps
#: orbit phase only approximately.
SWIFT_ORBIT = 5700.0
P = en.Params(quantum_s={"GBM": en.FERMI_ORBIT, "BAT": SWIFT_ORBIT})
P_DAY = en.Params(
    quantum_s={
        "GBM": en.FERMI_ORBIT,
        "BAT": SWIFT_ORBIT,
        "ICECAT": en.DAY,
        "CHIME": en.DAY,
        "GW": en.DAY,
    }
)
#: ASSUMPTIONs of the frozen control rule.
CONTROL_MIN_PAIRS = 50
CONTROL_MAX_P = 1.8e-5
#: ASSUMPTION: family-wise 3 sigma (one-sided), as D-074.
DETECT_P = e1.DETECT_P
#: 2-D Gaussian: r90 = 2.146 sigma.
R90_TO_SIGMA = 2.146


def bat_events(t: Table) -> Table:
    """swiftgrbba rows -> events: on-board Swift triggers with finite time, position and BAT error;
    1-sigma = error_radius (90 %, arcmin) / 2.146, in degrees."""

    def col(c):
        v = t[c]
        return np.ma.filled(np.ma.asarray(v, dtype=float), np.nan)

    obs = np.char.strip(np.asarray(t["trigger_obs"], dtype=str))
    mjd, ra, dec, er = col("time"), col("ra"), col("dec"), col("error_radius")
    keep = (
        (obs == "Swift") & np.isfinite(mjd) & np.isfinite(ra) & np.isfinite(dec) & np.isfinite(er)
    )
    out = en._events(
        np.asarray(t["name"], dtype=str)[keep],
        "BAT",
        mjd[keep],
        ra[keep],
        dec[keep],
        er[keep] / 60.0 / R90_TO_SIGMA,
        "HEASARC swiftgrbba (trigger_obs == Swift)",
    )
    return out


def fetch(refresh: bool) -> dict[str, Path]:
    d = data_root() / "e1_swift"
    d.mkdir(parents=True, exist_ok=True)
    pinned = {}
    if MANIFEST.exists() and not refresh:
        m = Table.read(MANIFEST, format="ascii.ecsv")
        pinned = {r["file"]: r["sha256"] for r in m}
    rows, paths = [], {}
    for fn, url in SOURCES.items():
        path = d / fn
        if not path.exists() or refresh:
            for k in range(5):  # HEASARC resets connections now and then
                try:
                    r = requests.get(url, timeout=600)
                    r.raise_for_status()
                    break
                except requests.RequestException:
                    if k == 4:
                        raise
                    time.sleep(2 ** (k + 1))
            path.write_bytes(r.content)
        content = path.read_bytes()
        h = e1.sha256(content)
        if fn in pinned and pinned[fn] != h:
            raise SystemExit(f"{fn}: sha256 {h} != pinned {pinned[fn]} (rerun with --refresh)")
        rows.append((url, fn, h, len(content)))
        paths[fn] = path
    if refresh or not MANIFEST.exists():
        t = Table(rows=rows, names=("uri", "file", "sha256", "size"))
        t.meta = {
            "provenance": str(Provenance.OBSERVED),
            "source": "HEASARC TAP swiftgrbba (Swift GRB compilation by burst advocate) plus the "
            "four D-074 E1 catalogues (scripts/e1_events.py SOURCES), re-pinned for this run",
            "retrieved": time.strftime("%Y-%m-%d", time.gmtime()),
            "pipeline_version": "0.0.1",
            "reason": "~5 MB total; event tables only",
        }
        t.write(MANIFEST, format="ascii.ecsv", overwrite=True)
    return paths


def load(paths: dict[str, Path]) -> dict[str, Table]:
    ev = e1.load(paths)
    ev["BAT"] = bat_events(Table.read(paths["swiftgrbba.vot"], format="votable"))
    return ev


# ------------------------------------------------------------------------------------- counting
_W: dict = {}


def count_all(samples: dict[str, en.Sample], p: en.Params) -> np.ndarray:
    """Counts, shape (n_channels, n_windows, 3)."""
    return np.stack(
        [en.count_channel(samples[a], samples[b], p, same=(a == b)) for a, b in CHANNELS]
    )


def _init(samples):
    os.environ["OMP_NUM_THREADS"] = "1"
    # only the catalogues some channel uses are scrambled and counted
    _W["s"] = {k: samples[k] for k in {c for ch in CHANNELS for c in ch}}


def _null_chunk(args):
    kind, seeds = args
    s = _W["s"]
    q = P_DAY if kind == "jitday" else P
    out = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        out.append(count_all({k: en.scramble_jit(v, rng, q) for k, v in s.items()}, P))
    return np.stack(out)


def run_null(kind: str, n: int, samples, cpu: int, base_seed: int) -> np.ndarray:
    seeds = np.arange(n) + base_seed
    chunks = [(kind, seeds[i : i + 250]) for i in range(0, n, 250)]
    t0 = time.time()
    if cpu <= 1:
        _init(samples)
        res = [_null_chunk(c) for c in chunks]
    else:
        with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(samples,)) as ex:
            res = list(ex.map(_null_chunk, chunks))
    print(f"null {kind}: {n} scrambles in {time.time() - t0:.0f} s", flush=True)
    return np.concatenate(res)


def test_mask(p: en.Params) -> np.ndarray:
    """Tested family: lag bins only; same/wide for localized channels, all for BAT-GW; the
    control channel is never tested."""
    nwin, nlag = len(en.windows_for(p)), len(p.lag_edges) - 1
    m = np.zeros((len(CHANNELS), nwin, 3), bool)
    for c, ch in enumerate(CHANNELS):
        if ch == CONTROL:
            continue
        if "GW" in ch:
            m[c, :nlag, 2] = True
        else:
            m[c, :nlag, :2] = True
    return m


# ------------------------------------------------------------------------------------- frozen rules


def control_passes(obs_pairs: int, null_sums: np.ndarray) -> tuple[bool, float]:
    """Frozen control rule: >= CONTROL_MIN_PAIRS same-direction BAT-GBM pairs at 0-100 s and the
    summed count's analytic p (jit null ensemble of the same sum) <= CONTROL_MAX_P."""
    pa = float(en.analytic_p(np.array([obs_pairs]), np.asarray(null_sums, float)[:, None])[0])
    return (obs_pairs >= CONTROL_MIN_PAIRS) and (pa <= CONTROL_MAX_P), pa


def outcome(cell_p: np.ndarray, n_cells: int) -> str:
    """Frozen outcome: 'excess' if any tested cell's analytic p x n_cells <= DETECT_P."""
    return "excess" if float(np.min(cell_p)) * n_cells <= DETECT_P else "null"


# ------------------------------------------------------------------------------------- main


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--jit", type=int, default=10000)
    ap.add_argument("--jitday", type=int, default=10000)
    ap.add_argument("--inject", type=int, default=40, help="trials per (channel, bin, n)")
    ap.add_argument("--cpu", type=int, default=1)
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args(argv)
    p = P
    OUT.mkdir(parents=True, exist_ok=True)
    paths = fetch(a.refresh)
    ev = load(paths)
    samples = {k: en.Sample.from_table(v) for k, v in ev.items()}
    for k, v in ev.items():
        print(k, len(v), f"MJD {v['mjd'].min():.1f}-{v['mjd'].max():.1f}", flush=True)

    obs = count_all(samples, p)
    jit = run_null("jit", a.jit, samples, a.cpu, 2_000_000)
    jday = run_null("jitday", a.jitday, samples, a.cpu, 3_000_000)
    mask = test_mask(p)
    n_cells = int(mask.sum())

    # positive control first: no result if it fails
    cc = CHANNELS.index(CONTROL)
    c_obs = int(obs[cc, 0, 0] + obs[cc, 1, 0])
    c_null = jit[:, cc, 0, 0] + jit[:, cc, 1, 0]
    c_ok, c_p = control_passes(c_obs, c_null)

    wins = [w for w, _ in en.windows_for(p)]
    labels = en.lag_labels(p) + wins[len(p.lag_edges) - 1 :]
    p_jit, z_jit, pa_jit = en.empirical_p(obs, jit), en.z_score(obs, jit), en.analytic_p(obs, jit)
    p_jday, z_jday = en.empirical_p(obs, jday), en.z_score(obs, jday)
    ul = en.upper_limit_95(obs, jit)
    rows = []
    for c, (ca, cb) in enumerate(CHANNELS):
        for k, lab in enumerate(labels):
            for j, cl in enumerate(e1.CLASSES):
                if ("GW" in (ca, cb)) != (cl == "all") and cl != "all":
                    continue
                rows.append(
                    (
                        f"{ca}-{cb}", lab, cl, bool(mask[c, k, j]), int(obs[c, k, j]),
                        round(float(jit[:, c, k, j].mean()), 2),
                        round(float(jit[:, c, k, j].std()), 2),
                        round(float(z_jit[c, k, j]), 2), float(p_jit[c, k, j]),
                        float(pa_jit[c, k, j]),
                        round(float(jday[:, c, k, j].mean()), 2),
                        round(float(z_jday[c, k, j]), 2), float(p_jday[c, k, j]),
                        round(float(ul[c, k, j]), 2),
                    )
                )  # fmt: skip
    counts = Table(
        rows=rows,
        names=(
            "channel", "lag", "cls", "tested", "obs", "jit_mean", "jit_sd", "z_jit", "p_jit",
            "p_jit_analytic", "jday_mean", "z_jday", "p_jday", "ul95_pairs",
        ),
    )  # fmt: skip
    counts.meta = {
        "provenance": str(Provenance.DERIVED),
        "source": "scripts/e1_swift.py on data/manifests/e1_swift.ecsv; "
        f"jit null {a.jit}, jitday null {a.jitday}",
    }
    counts.write(OUT / "counts.ecsv", format="ascii.ecsv", overwrite=True)

    pa_t = pa_jit[mask]
    cells = [(f"{x}-{y}", lab, cl) for x, y in CHANNELS for lab in labels for cl in e1.CLASSES]
    tested_cells = [c for c, m in zip(cells, mask.ravel(), strict=True) if m]
    jmin = int(np.argmin(pa_t))
    res = outcome(pa_t, n_cells)

    dup = {}
    for ca, cb in (("BAT", "BAT"), ("BAT", "ICECAT"), ("BAT", "CHIME"), CONTROL):
        t = e1.close_pairs(
            ev[ca], ev[cb], p, ca == cb, 100.0 if (ca, cb) == CONTROL else p.lag_edges[-1]
        )
        t.write(OUT / f"same_dir_pairs_{ca}-{cb}.ecsv", format="ascii.ecsv", overwrite=True)
        dup[f"{ca}-{cb}"] = len(t)

    # injections and limits through e1_events' code with this run's channel list (in-process only)
    e1.CHANNELS = TESTED
    inj, p_star = e1.injections(
        samples, p, jit[:, : len(TESTED)], mask[: len(TESTED)], a.inject, np.random.default_rng(11)
    )
    inj.write(OUT / "injections.ecsv", format="ascii.ecsv", overwrite=True)
    lim = e1.limits(counts, inj, ev, p)
    lim.write(OUT / "limits.ecsv", format="ascii.ecsv", overwrite=True)

    summary = {
        "provenance": "derived",
        "plan": "results/e1_swift/PLAN.md",
        "events": {k: len(v) for k, v in ev.items()},
        "bat_mjd_range": [float(ev["BAT"]["mjd"].min()), float(ev["BAT"]["mjd"].max())],
        "n_jit": a.jit,
        "n_jitday": a.jitday,
        "tested_cells": n_cells,
        "positive_control_bat_gbm_same_0_100s": {
            "obs": c_obs,
            "jit_mean": float(c_null.mean()),
            "jit_sd": float(c_null.std()),
            "analytic_p": c_p,
            "passes": bool(c_ok),
            "rule": f">= {CONTROL_MIN_PAIRS} pairs and analytic p <= {CONTROL_MAX_P}",
        },
        "outcome": res if c_ok else "no result (control failed)",
        "min_cell": {
            "cell": " ".join(tested_cells[jmin]),
            "analytic_p": float(pa_t[jmin]),
            "bonferroni_p": float(min(1.0, pa_t[jmin] * n_cells)),
        },
        "global_jit_pooled": dict(
            zip(("min_cell_p", "trials_corrected_p"), en.global_p(obs, jit, mask), strict=True)
        ),
        "global_jitday_pooled": dict(
            zip(("min_cell_p", "trials_corrected_p"), en.global_p(obs, jday, mask), strict=True)
        ),
        "injection_detection_p_star": p_star,
        "same_direction_pairs": dup,
        "params": {k: (list(v) if isinstance(v, tuple) else v) for k, v in p.__dict__.items()},
        "swift_orbit_s": SWIFT_ORBIT,
    }
    with open(OUT / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
        fh.write("\n")
    print(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "events",
                    "positive_control_bat_gbm_same_0_100s",
                    "outcome",
                    "min_cell",
                    "global_jit_pooled",
                )
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
