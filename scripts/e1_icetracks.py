"""E1 + IceTracks-DR2: northern IceCube tracks as a seventh node set in the E1 lag x separation pair
counts (hypothesis; owner idea 4, D-074). Frozen plan: ``results/e1_icetracks/PLAN.md``.

Tested family: DR2-GBM, DR2-BAT, DR2-CHIME (same, wide) and DR2-GW (all) over the five D-074 lag
bins. DR2-ICECAT same-direction pairs at 0-10 s (the same neutrinos) are the positive control.

Inputs: the IceTracks-DR2 files pinned in ``data/manifests/nf_icetracks_dr2.ecsv`` (D-079, read by
``nf_ghost_dr2.fetch``) and the five catalogues pinned in ``data/manifests/e1_swift.ecsv``.

  python scripts/e1_icetracks.py [--jit 10000] [--jitday 10000] [--inject 40] [--cpu 4]
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
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import e1_events as e1  # noqa: E402
import e1_swift as sw  # noqa: E402
import nf_ghost_dr2 as dr2  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402
from jwst_anomaly import neutrino_frontier as nf  # noqa: E402
from jwst_anomaly.schema import Provenance  # noqa: E402

OUT = ROOT / "results" / "e1_icetracks"
#: ASSUMPTIONs of the node set (frozen plan, as D-079).
DEC_MIN = -5.0
LOG10E_MIN = 4.0
TESTED = (("DR2", "GBM"), ("DR2", "BAT"), ("DR2", "CHIME"), ("DR2", "GW"))
CONTROL = ("DR2", "ICECAT")
CHANNELS = TESTED + (CONTROL,)
P = sw.P
P_DAY = sw.P_DAY
#: ASSUMPTIONs of the frozen control rule.
CONTROL_MIN_PAIRS = 10
CONTROL_MAX_P = 1.8e-5
DETECT_P = e1.DETECT_P


def dr2_events(ev: np.ndarray) -> Table:
    """Reduced DR2 rows (run, event, subevent, mjd, log10E, AngErr, ra, dec, ...) -> events: the
    first subevent of each (run, event) readout, kept if Dec > DEC_MIN and log10 E >= LOG10E_MIN;
    1-sigma = AngErr (deg)."""
    e = ev[np.lexsort((ev[:, 2], ev[:, 1], ev[:, 0]))]
    _, first = np.unique(np.round(e[:, :2]).astype(np.int64), axis=0, return_index=True)
    e = e[np.sort(first)]
    e = e[(e[:, 7] > DEC_MIN) & (e[:, 4] >= LOG10E_MIN)]
    e = e[np.argsort(e[:, 3], kind="stable")]
    names = [f"r{int(r)}e{int(x)}" for r, x in np.round(e[:, :2])]
    return en._events(
        names, "DR2", e[:, 3], e[:, 6], e[:, 7], e[:, 5],
        f"IceTracks-DR2 v3.1 (Dec > {DEC_MIN}, log10E >= {LOG10E_MIN})",
    )  # fmt: skip


# ------------------------------------------------------------------------------------- counting
_W: dict = {}


def count_all(samples: dict[str, en.Sample], p: en.Params) -> np.ndarray:
    """Counts, shape (n_channels, n_windows, 3)."""
    return np.stack([en.count_channel(samples[a], samples[b], p, same=False) for a, b in CHANNELS])


def scramble(kind: str, s: en.Sample, rng, up: np.ndarray) -> en.Sample:
    """jit / jitday: DR2 by the uptime-aware +-3 d jitter (D-079); others as #351."""
    if s.cat == "DR2":
        return nf.jitter_uptime(s, up[:, 0], up[:, 1], P.jitter_s / en.DAY, rng)
    return en.scramble_jit(s, rng, P_DAY if kind == "jitday" else P)


def _init(samples, up):
    os.environ["OMP_NUM_THREADS"] = "1"
    # fixed order, so a seed gives the same draws in every interpreter (str hashing is salted)
    _W["s"] = {k: samples[k] for k in sorted({c for ch in CHANNELS for c in ch})}
    _W["up"] = up


def _null_chunk(args):
    kind, seeds = args
    out = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        out.append(count_all({k: scramble(kind, v, rng, _W["up"]) for k, v in _W["s"].items()}, P))
    return np.stack(out)


def run_null(kind: str, n: int, samples, up, cpu: int, base_seed: int) -> np.ndarray:
    seeds = np.arange(n) + base_seed
    chunks = [(kind, seeds[i : i + 250]) for i in range(0, n, 250)]
    t0 = time.time()
    if cpu <= 1:
        _init(samples, up)
        res = [_null_chunk(c) for c in chunks]
    else:
        with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(samples, up)) as ex:
            res = list(ex.map(_null_chunk, chunks))
    print(f"null {kind}: {n} scrambles in {time.time() - t0:.0f} s", flush=True)
    return np.concatenate(res)


def test_mask(p: en.Params) -> np.ndarray:
    """Tested family: lag bins only; same/wide for localized channels, all for DR2-GW; the control
    channel is never tested."""
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


def control_passes(obs_pairs: int, null_counts: np.ndarray) -> tuple[bool, float]:
    """Frozen control rule: >= CONTROL_MIN_PAIRS same-direction DR2-ICECAT pairs at 0-10 s and
    their analytic p (jit null ensemble of the same cell) <= CONTROL_MAX_P."""
    pa = float(en.analytic_p(np.array([obs_pairs]), np.asarray(null_counts, float)[:, None])[0])
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
    a = ap.parse_args(argv)
    p = P
    OUT.mkdir(parents=True, exist_ok=True)
    raw, up = dr2.fetch()
    ev = sw.load(sw.fetch(False))
    ev["DR2"] = dr2_events(raw)
    samples = {k: en.Sample.from_table(v) for k, v in ev.items()}
    for k, v in ev.items():
        print(k, len(v), f"MJD {v['mjd'].min():.1f}-{v['mjd'].max():.1f}", flush=True)

    obs = count_all(samples, p)
    jit = run_null("jit", a.jit, samples, up, a.cpu, 6_000_000)
    jday = run_null("jitday", a.jitday, samples, up, a.cpu, 7_000_000)
    mask = test_mask(p)
    n_cells = int(mask.sum())

    # positive control first: no result if it fails
    cc = CHANNELS.index(CONTROL)
    c_obs = int(obs[cc, 0, 0])
    c_ok, c_p = control_passes(c_obs, jit[:, cc, 0, 0])

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
        "source": "scripts/e1_icetracks.py on data/manifests/nf_icetracks_dr2.ecsv and "
        f"data/manifests/e1_swift.ecsv; jit null {a.jit}, jitday null {a.jitday}",
    }
    counts.write(OUT / "counts.ecsv", format="ascii.ecsv", overwrite=True)

    pa_t = pa_jit[mask]
    cells = [(f"{x}-{y}", lab, cl) for x, y in CHANNELS for lab in labels for cl in e1.CLASSES]
    tested_cells = [c for c, m in zip(cells, mask.ravel(), strict=True) if m]
    jmin = int(np.argmin(pa_t))
    res = outcome(pa_t, n_cells)

    dup = {}
    for ca, cb in (("DR2", "GBM"), ("DR2", "BAT"), ("DR2", "CHIME"), CONTROL):
        t = e1.close_pairs(ev[ca], ev[cb], p, False, 10.0 if (ca, cb) == CONTROL else 3600.0)
        t.write(OUT / f"same_dir_pairs_{ca}-{cb}.ecsv", format="ascii.ecsv", overwrite=True)
        dup[f"{ca}-{cb}"] = len(t)

    # injections and limits through e1_events' code with this run's channel list (in-process only)
    e1.CHANNELS = TESTED
    inj, p_star = e1.injections(
        samples, p, jit[:, : len(TESTED)], mask[: len(TESTED)], a.inject, np.random.default_rng(13)
    )
    inj.write(OUT / "injections.ecsv", format="ascii.ecsv", overwrite=True)
    lim = e1.limits(counts, inj, ev, p)
    lim.write(OUT / "limits.ecsv", format="ascii.ecsv", overwrite=True)

    summary = {
        "provenance": "derived",
        "plan": "results/e1_icetracks/PLAN.md",
        "events": {k: len(v) for k, v in ev.items()},
        "dr2_mjd_range": [float(ev["DR2"]["mjd"].min()), float(ev["DR2"]["mjd"].max())],
        "n_jit": a.jit,
        "n_jitday": a.jitday,
        "tested_cells": n_cells,
        "positive_control_dr2_icecat_same_0_10s": {
            "obs": c_obs,
            "jit_mean": float(jit[:, cc, 0, 0].mean()),
            "jit_sd": float(jit[:, cc, 0, 0].std()),
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
        "node_set": {"dec_min": DEC_MIN, "log10e_min": LOG10E_MIN, "one_per_readout": True},
        "params": {k: (list(v) if isinstance(v, tuple) else v) for k, v in p.__dict__.items()},
        "jitday_quantum_s": P_DAY.quantum_s,
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
                    "positive_control_dr2_icecat_same_0_10s",
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
