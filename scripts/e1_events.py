"""E1 causal event network: pair-count excess versus time lag and separation between Fermi GBM
bursts, IceCube ICECAT-1 tracks, GWTC events and CHIME/FRB Catalog 2 bursts (hypothesis; owner
idea 4).

Inputs (~5 MB) are downloaded to ``$JWST_ANOMALY_DATA/e1_events/`` and pinned by sha256 in
``data/manifests/e1_events.ecsv`` (a mismatch is refused unless ``--refresh``). Outputs (derived,
small) go to ``results/e1_events/``.

  python scripts/e1_events.py [--perm 10000] [--jit 10000] [--inject 40] [--cpu 4] [--refresh]
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import io
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

from jwst_anomaly import event_network as en  # noqa: E402
from jwst_anomaly.paths import data_root  # noqa: E402
from jwst_anomaly.schema import Provenance  # noqa: E402

OUT = ROOT / "results" / "e1_events"
MANIFEST = ROOT / "data" / "manifests" / "e1_events.ecsv"
SOURCES = {
    "fermigbrst.vot": "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync?REQUEST=doQuery&LANG=ADQL&QUERY="
    "SELECT%20trigger_name,ra,dec,error_radius,trigger_time,t90,fluence,last_modified%20FROM%20fermigbrst",
    "icecat1_gold_bronze_tracks.csv": "https://dataverse.harvard.edu/api/access/datafile/7502710?format=original",
    "gwtc_events.csv": "https://gwosc.org/api/v2/catalogs/GWTC/events?include-default-parameters=true&format=csv",
    "chimefrbcat2.csv": "https://www.canfar.net/storage/vault/file/AstroDataCitationDOI/CISTI.CANFAR/25.0066/"
    "data/table/chimefrbcat2.csv",
}
#: Channels (A, B). Same-catalogue channels count each unordered pair once.
CHANNELS = (
    ("GBM", "GBM"),
    ("GBM", "ICECAT"),
    ("GBM", "GW"),
    ("ICECAT", "ICECAT"),
    ("ICECAT", "GW"),
    ("CHIME", "CHIME"),
    ("CHIME", "GBM"),
    ("CHIME", "ICECAT"),
    ("CHIME", "GW"),
)
CLASSES = ("same", "wide", "all")
#: GW170817 / GRB 170817A positive control (Abbott+2017, ApJL 848, L13: GRB onset 1.74 +- 0.05 s
#: after merger).
PC_GW, PC_GBM = "GW170817", "bn170817529"
#: NGC 4993 (host of GW170817; Coulter+2017, Science 358, 1556), used only to check the GBM row's
#: position.
NGC4993 = (197.450374, -23.381495)
#: Vetting null "jitday": like "jit" but ground-instrument and GW shifts are whole solar days
#: (+- slop), so time-of-day structure (daytime RFI, daily maintenance, duty cycles) is kept.
P_DAY = en.Params(
    quantum_s={"GBM": en.FERMI_ORBIT, "ICECAT": en.DAY, "CHIME": en.DAY, "GW": en.DAY}
)
#: ASSUMPTION: global 3-sigma (one-sided) defines "detected" for injections.
DETECT_P = 1.35e-3


def sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def fetch(refresh: bool) -> dict[str, Path]:
    d = data_root() / "e1_events"
    d.mkdir(parents=True, exist_ok=True)
    pinned = {}
    if MANIFEST.exists() and not refresh:
        m = Table.read(MANIFEST, format="ascii.ecsv")
        pinned = {r["file"]: r["sha256"] for r in m}
    rows, paths = [], {}
    for fn, url in SOURCES.items():
        path = d / fn
        if not path.exists() or refresh:
            r = requests.get(url, timeout=600)
            r.raise_for_status()
            path.write_bytes(r.content)
        content = path.read_bytes()
        h = sha256(content)
        if fn in pinned and pinned[fn] != h:
            raise SystemExit(
                f"{fn}: sha256 {h} != pinned {pinned[fn]} (rerun with --refresh to re-pin)"
            )
        rows.append((url, fn, h, len(content)))
        paths[fn] = path
    if refresh or not MANIFEST.exists():
        t = Table(rows=rows, names=("uri", "file", "sha256", "size"))
        t.meta = {
            "provenance": str(Provenance.OBSERVED),
            "source": "E1 event catalogues: HEASARC TAP fermigbrst; ICECAT-1 v4 "
            "(doi:10.7910/DVN/SCRUCD, IceCube_Gold_Bronze_Tracks original CSV); GWOSC GWTC "
            "cumulative event list; CHIME/FRB Catalog 2 "
            "(doi:10.11570/25.0066)",
            "retrieved": time.strftime("%Y-%m-%d", time.gmtime()),
            "pipeline_version": "0.0.1",
            "reason": "~4.7 MB total; event tables only (no sky maps, no exposure file)",
        }
        t.write(MANIFEST, format="ascii.ecsv", overwrite=True)
    return paths


def load(paths: dict[str, Path], chime_all: bool = False) -> dict[str, Table]:
    import pandas as pd

    gbm = Table.read(paths["fermigbrst.vot"], format="votable")
    for c in gbm.colnames:
        if gbm[c].dtype.kind == "O":
            gbm[c] = gbm[c].astype(str)
    return {
        "GBM": en.gbm_events(gbm),
        "ICECAT": en.icecat_events(pd.read_csv(paths["icecat1_gold_bronze_tracks.csv"])),
        "GW": en.gw_events(pd.read_csv(paths["gwtc_events.csv"])),
        "CHIME": en.chime_events(
            pd.read_csv(paths["chimefrbcat2.csv"]), one_per_source=not chime_all
        ),
    }


# ------------------------------------------------------------------------------------- counting

_W: dict = {}


def _init(samples, p):
    os.environ["OMP_NUM_THREADS"] = "1"
    _W["s"], _W["p"] = samples, p


def count_all(samples: dict[str, en.Sample], p: en.Params) -> np.ndarray:
    """Counts, shape (n_channels, n_windows, 3)."""
    return np.stack(
        [en.count_channel(samples[a], samples[b], p, same=(a == b)) for a, b in CHANNELS]
    )


def _null_chunk(args):
    kind, seeds = args
    s, p = _W["s"], _W["p"]
    out = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        if kind == "perm":
            sc = {k: en.scramble_perm(v, rng) for k, v in s.items()}
        elif kind == "jitday":
            sc = {k: en.scramble_jit(v, rng, P_DAY) for k, v in s.items()}
        else:
            sc = {k: en.scramble_jit(v, rng, p) for k, v in s.items()}
        out.append(count_all(sc, p))
    return np.stack(out)


def run_null(kind: str, n: int, samples, p, cpu: int, base_seed: int) -> np.ndarray:
    seeds = np.arange(n) + base_seed
    chunks = [(kind, seeds[i : i + 250]) for i in range(0, n, 250)]
    t0 = time.time()
    if cpu <= 1:  # shared machines: run in-process (OMP_NUM_THREADS=1 is the caller's job)
        _init(samples, p)
        res = [_null_chunk(c) for c in chunks]
    else:
        with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(samples, p)) as ex:
            res = list(ex.map(_null_chunk, chunks))
    print(f"null {kind}: {n} scrambles in {time.time() - t0:.0f} s", flush=True)
    return np.concatenate(res)


def test_mask(p: en.Params) -> np.ndarray:
    """Cells in the pre-registered trials family: lag bins only; same/wide for localized channels,
    all for GW channels."""
    nwin = len(en.windows_for(p))
    nlag = len(p.lag_edges) - 1
    m = np.zeros((len(CHANNELS), nwin, 3), bool)
    for c, (a, b) in enumerate(CHANNELS):
        if "GW" in (a, b):
            m[c, :nlag, 2] = True
        else:
            m[c, :nlag, :2] = True
    return m


# ------------------------------------------------------------------------------------- controls


def positive_control(ev: dict[str, Table]) -> dict:
    g = ev["GBM"][ev["GBM"]["name"] == PC_GBM][0]
    w = ev["GW"][ev["GW"]["name"] == PC_GW][0]
    return {
        "gw": PC_GW,
        "gbm": PC_GBM,
        "dt_s": float((g["mjd"] - w["mjd"]) * en.DAY),
        "published_dt_s": 1.74,
        "gbm_row_sep_from_ngc4993_deg": float(en.sep_deg(g["ra"], g["dec"], *NGC4993)),
        "gbm_row_sigma_deg": float(g["sigma"]),
        "note": "GWTC CSV has no sky position, so GW channels are lag-only (class 'all'); "
        "the GBM catalogue "
        "row for 170817A carries the optical-counterpart position (error_radius = 0).",
    }


def close_pairs(
    a: Table, b: Table, p: en.Params, same: bool, max_lag_s: float, cls_keep=(0,)
) -> Table:
    """List pairs within max_lag in the given separation classes (duplicate / re-trigger check)."""
    sa, sb = en.Sample.from_table(a), en.Sample.from_table(b)
    i, j = en.pairs_within(sa.mjd, sb.mjd, max_lag_s, same)
    sep = en.sep_deg(sa.ra[i], sa.dec[i], sb.ra[j], sb.dec[j])
    cls = en.sep_class(sep, sa.sigma[i], sb.sigma[j], p)
    m = np.isin(cls, cls_keep)
    t = Table(
        {
            "a": a["name"][i[m]],
            "b": b["name"][j[m]],
            "lag_s": np.round(np.abs(sb.mjd[j[m]] - sa.mjd[i[m]]) * en.DAY, 3),
            "sep_deg": np.round(sep[m], 3),
            "sigma_comb_deg": np.round(np.hypot(sa.sigma[i[m]], sb.sigma[j[m]]), 3),
        }
    )
    t.sort("lag_s")
    t.meta = {
        "provenance": str(Provenance.DERIVED),
        "source": f"{a.meta['source']} x {b.meta['source']}",
    }
    return t


# ------------------------------------------------------------------------------------- injection


def injections(samples, p, null_jit, mask, n_trials, rng) -> Table:
    """Inject n wide-separation pairs at a known lag bin into channel (A, B); measure the recovered
    count excess and whether the global trials-corrected p (jit null) falls below DETECT_P."""
    obs0 = count_all(samples, p)
    nlag = len(p.lag_edges) - 1
    rows = []
    # precompute per-scramble min-p distribution once (detection threshold on the min per-cell p)
    flat = null_jit[:, mask]
    srt = np.sort(flat, axis=0)
    nn = flat.shape[0]
    p_null = np.empty_like(flat, dtype=float)
    for c in range(flat.shape[1]):
        p_null[:, c] = (nn - np.searchsorted(srt[:, c], flat[:, c], side="left")) / nn
    min_null = np.sort(p_null.min(axis=1))
    for c, (a, b) in enumerate(CHANNELS):
        cls = 2 if "GW" in (a, b) else 1
        for k in range(nlag):
            lo, hi = p.lag_edges[k], p.lag_edges[k + 1]
            base = float(null_jit[:, c, k, cls].mean())
            sd = float(null_jit[:, c, k, cls].std())
            grid = sorted(
                {max(1, int(round(x))) for x in (np.array([1, 2, 3, 5, 8]) * max(sd, 1.0))}
            )
            for n in grid:
                det, rec = 0, []
                for _ in range(n_trials):
                    sb = en.inject_pairs(samples[a], samples[b], n, lo, hi, p, a == b, rng)
                    # only channel c is recounted (the moved B events also enter other channels;
                    # those changes are second order and ignored)
                    o = obs0.copy()
                    o[c] = en.count_channel(sb if a == b else samples[a], sb, p, same=(a == b))
                    rec.append(o[c, k, cls] - obs0[c, k, cls])
                    pc = (1.0 + (flat >= o[mask][None]).sum(axis=0)) / (nn + 1.0)
                    g = (1.0 + np.searchsorted(min_null, pc.min(), side="right")) / (nn + 1.0)
                    det += g < DETECT_P
                rows.append(
                    (
                        f"{a}-{b}",
                        CLASSES[cls],
                        k,
                        n,
                        base,
                        sd,
                        float(np.mean(rec)) / n,
                        det / n_trials,
                    )
                )
    t = Table(
        rows=rows,
        names=("channel", "cls", "lag_bin", "n_inj", "null_mean", "null_sd", "eff", "det_frac"),
    )
    t.meta = {
        "provenance": str(Provenance.SIMULATED),
        "source": "synthetic wide-separation lagged pairs injected into the real catalogues "
        "(event_network.inject_pairs); detection = global trials-corrected p (jit null) < "
        f"{DETECT_P}",
    }
    return t


# ------------------------------------------------------------------------------------- main


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm", type=int, default=10000)
    ap.add_argument("--jit", type=int, default=10000)
    ap.add_argument("--jitday", type=int, default=10000)
    ap.add_argument("--inject", type=int, default=40, help="trials per (channel, bin, n)")
    ap.add_argument("--cpu", type=int, default=1)
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args(argv)
    p = en.Params()
    OUT.mkdir(parents=True, exist_ok=True)
    paths = fetch(a.refresh)
    ev = load(paths)
    samples = {k: en.Sample.from_table(v) for k, v in ev.items()}
    for k, v in ev.items():
        print(k, len(v), f"MJD {v['mjd'].min():.1f}-{v['mjd'].max():.1f}", flush=True)

    obs = count_all(samples, p)
    perm = run_null("perm", a.perm, samples, p, a.cpu, 1_000_000)
    jit = run_null("jit", a.jit, samples, p, a.cpu, 2_000_000)
    jday = run_null("jitday", a.jitday, samples, p, a.cpu, 3_000_000)
    mask = test_mask(p)

    wins = [w for w, _ in en.windows_for(p)]
    labels = en.lag_labels(p) + wins[len(p.lag_edges) - 1 :]
    p_perm, p_jit = en.empirical_p(obs, perm), en.empirical_p(obs, jit)
    z_perm, z_jit = en.z_score(obs, perm), en.z_score(obs, jit)
    ul = en.upper_limit_95(obs, jit)
    p_jday, z_jday = en.empirical_p(obs, jday), en.z_score(obs, jday)
    rows = []
    for c, (ca, cb) in enumerate(CHANNELS):
        for k, lab in enumerate(labels):
            for j, cl in enumerate(CLASSES):
                if ("GW" in (ca, cb)) != (cl == "all") and cl != "all":
                    continue
                rows.append(
                    (
                        f"{ca}-{cb}",
                        lab,
                        cl,
                        bool(mask[c, k, j]),
                        int(obs[c, k, j]),
                        round(float(perm[:, c, k, j].mean()), 2),
                        round(float(perm[:, c, k, j].std()), 2),
                        round(float(z_perm[c, k, j]), 2),
                        float(p_perm[c, k, j]),
                        round(float(jit[:, c, k, j].mean()), 2),
                        round(float(jit[:, c, k, j].std()), 2),
                        round(float(z_jit[c, k, j]), 2),
                        float(p_jit[c, k, j]),
                        round(float(jday[:, c, k, j].mean()), 2),
                        round(float(z_jday[c, k, j]), 2),
                        float(p_jday[c, k, j]),
                        round(float(ul[c, k, j]), 2),
                    )
                )
    counts = Table(
        rows=rows,
        names=(
            "channel", "lag", "cls", "tested", "obs", "perm_mean", "perm_sd", "z_perm", "p_perm",
            "jit_mean", "jit_sd", "z_jit", "p_jit", "jday_mean", "z_jday", "p_jday", "ul95_pairs",
        ),
    )  # fmt: skip
    counts.meta = {
        "provenance": str(Provenance.DERIVED),
        "source": "scripts/e1_events.py on data/manifests/e1_events.ecsv; "
        f"perm null {a.perm}, jit null {a.jit}, jitday null {a.jitday}",
    }
    counts.write(OUT / "counts.ecsv", format="ascii.ecsv", overwrite=True)

    gp_perm = en.global_p(obs, perm, mask & (np.arange(3) < 2)[None, None, :])
    gp_jit = en.global_p(obs, jit, mask)
    gp_jday = en.global_p(obs, jday, mask)
    core = mask & np.array(["CHIME" not in ch for ch in CHANNELS])[:, None, None]
    gp_core = en.global_p(obs, jit, core)

    # positive controls: GW170817 (jit null) and CHIME with every repeater burst as its own node
    pc = positive_control(ev)
    ev_all = load(paths, chime_all=True)
    s_all = en.Sample.from_table(ev_all["CHIME"])
    o_rep = en.count_channel(s_all, s_all, p, same=True)
    rng = np.random.default_rng(7)
    rep_null = np.stack(
        [en.count_channel(*(2 * [en.scramble_perm(s_all, rng)]), p, same=True) for _ in range(200)]
    )
    rep = {
        "n_events": len(s_all.mjd),
        "same_obs": o_rep[: len(p.lag_edges) - 1, 0].tolist(),
        "same_perm_mean": rep_null[:, : len(p.lag_edges) - 1, 0].mean(axis=0).round(2).tolist(),
        "same_perm_sd": rep_null[:, : len(p.lag_edges) - 1, 0].std(axis=0).round(2).tolist(),
    }

    dup = {}
    for ca, cb in (
        ("GBM", "GBM"),
        ("ICECAT", "ICECAT"),
        ("CHIME", "CHIME"),
        ("GBM", "ICECAT"),
        ("CHIME", "GBM"),
    ):
        t = close_pairs(ev[ca], ev[cb], p, ca == cb, p.lag_edges[-1])
        t.write(OUT / f"same_dir_pairs_{ca}-{cb}.ecsv", format="ascii.ecsv", overwrite=True)
        dup[f"{ca}-{cb}"] = len(t)

    inj = injections(samples, p, jit, mask, a.inject, np.random.default_rng(11))
    inj.write(OUT / "injections.ecsv", format="ascii.ecsv", overwrite=True)

    summary = {
        "provenance": "derived",
        "events": {k: len(v) for k, v in ev.items()},
        "n_perm": a.perm,
        "n_jit": a.jit,
        "tested_cells": int(mask.sum()),
        "global_perm": {"min_cell_p": gp_perm[0], "trials_corrected_p": gp_perm[1]},
        "global_jit": {"min_cell_p": gp_jit[0], "trials_corrected_p": gp_jit[1]},
        "global_jit_core_5_channels": {"min_cell_p": gp_core[0], "trials_corrected_p": gp_core[1]},
        "global_jitday": {"min_cell_p": gp_jday[0], "trials_corrected_p": gp_jday[1]},
        "n_jitday": a.jitday,
        "positive_control_gw170817": pc,
        "positive_control_chime_repeaters": rep,
        "same_direction_pairs_within_7d": dup,
        "params": {k: (list(v) if isinstance(v, tuple) else v) for k, v in p.__dict__.items()},
    }
    with open(OUT / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
        fh.write("\n")
    buf = io.StringIO()
    counts[counts["tested"]].write(buf, format="ascii.fixed_width_two_line")
    print(buf.getvalue())
    print(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "global_perm",
                    "global_jit",
                    "global_jit_core_5_channels",
                    "global_jitday",
                )
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
