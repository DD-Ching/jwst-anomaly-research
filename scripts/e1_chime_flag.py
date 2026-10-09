"""E1: ordinary-explanation tests on the open CHIME-CHIME 1 h-1 d wide flag (no uptime file needed).

(a) per-day concentration of the excess and a 30-day block jackknife;
(b) epochs and seasons (calendar year, Catalog 1 period via ``catalog1_flag``, quarter of year);
(c) property independence of the excess pairs (DM, fluence, S/N, Dec, sidereal-phase difference)
    against pairs from scrambled catalogues;
(d) a rate-modulated null: times redrawn from a smooth per-day rate estimated from the catalogue
    itself (running mean of daily counts), calibrated by injection.

Writes ``results/e1_events/chime_flag_tests.json``. Single process.
``--cell 100s-1h`` runs only (d) on the 100 s - 1 h wide cell, with a second variant that keeps
the catalogue's time-of-day distribution (``keep_tod``), and writes
``results/e1_events/chime_rate_null_100s_1h.json``.

  python scripts/e1_chime_flag.py [--n 1000] [--cell 1h-1d|100s-1h]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import e1_events as E  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402

OUT = ROOT / "results" / "e1_events" / "chime_flag_tests.json"
OUT_100S = ROOT / "results" / "e1_events" / "chime_rate_null_100s_1h.json"
K = 3  # 1 h - 1 d lag bin
CELLS = {"1h-1d": 3, "100s-1h": 2}
WIDE = 1


def chime_frame(paths) -> pd.DataFrame:
    """The CHIME rows behind en.chime_events (one node per source), in the same (time) order."""
    df = pd.read_csv(paths["chimefrbcat2.csv"])
    df = df[df["sub_num"] == 0].sort_values("mjd_400")
    key = df["repeater_name"].fillna(df["tns_name"]).astype(str)
    df = df[~key.duplicated()]
    return df.sort_values("mjd_400", kind="stable").reset_index(drop=True)


def cell_pairs(s: en.Sample, p: en.Params):
    """Index pairs (i, j) in the 1 h - 1 d wide cell."""
    lo, hi = p.lag_edges[K], p.lag_edges[K + 1]
    i, j = en.pairs_within(s.mjd, s.mjd, hi, True)
    lag = np.abs(s.mjd[j] - s.mjd[i]) * en.DAY
    sep = en.sep_deg(s.ra[i], s.dec[i], s.ra[j], s.dec[j])
    cls = en.sep_class(sep, s.sigma[i], s.sigma[j], p)
    m = (lag >= lo) & (lag < hi) & (cls == WIDE)
    return i[m], j[m]


def count(s: en.Sample, p: en.Params, k: int = K) -> int:
    return int(en.count_channel(s, s, p, True, windows=[("c", [p.lag_edges[k : k + 2]])])[0, 1])


def jit_null(s, p, n, rng):
    return np.array([count(en.scramble_jit(s, rng, p), p) for _ in range(n)])


def z_of(o, null):
    return float((o - null.mean()) / null.std()) if null.std() > 0 else 0.0


def subset(s: en.Sample, m: np.ndarray) -> en.Sample:
    return en.Sample(s.cat, s.mjd[m], s.ra[m], s.dec[m], s.sigma[m], s.year[m])


def rate_null(s: en.Sample, rng, window_days: int, keep_tod: bool = False) -> en.Sample:
    """Permute events, then redraw every time from a smooth per-day rate: the running mean (window
    days) of the catalogue's own daily counts; uniform time within the drawn day, or with
    ``keep_tod`` the permuted event's own UTC time of day (keeps the daily duty cycle). Dec and
    hour angle kept (with_times)."""
    day = np.floor(s.mjd).astype(int)
    d0 = day.min()
    cnt = np.bincount(day - d0).astype(float)
    kern = np.ones(window_days) / window_days
    rate = np.convolve(cnt, kern, mode="same")
    prob = rate / rate.sum()
    perm = rng.permutation(len(s.mjd))
    days = rng.choice(len(prob), size=len(s.mjd), p=prob) + d0
    new = np.empty_like(s.mjd)
    tod = np.mod(s.mjd[perm], 1.0) if keep_tod else rng.uniform(0, 1, len(s.mjd))
    new[perm] = days + tod
    return s.with_times(new)


def rate_null_test(s, p, k, n, rng, windows=(7, 3), keep_tod=False, n_inj=300) -> dict:
    """Observed vs rate-modulated null in lag bin k (wide), and the same after injecting n_inj
    wide pairs in that bin (the share of injected excess the null keeps is its calibration)."""
    rn = {}
    for w in windows:
        res = {}
        for ninj in (0, n_inj):
            sb = (
                s
                if ninj == 0
                else en.inject_pairs(s, s, ninj, *p.lag_edges[k : k + 2], p, True, rng)
            )
            ob = count(sb, p, k)
            nb = np.array([count(rate_null(sb, rng, w, keep_tod), p, k) for _ in range(n)])
            pa = float(en.analytic_p(np.array([ob]), nb[:, None])[0])
            res[str(ninj)] = {
                "obs": ob,
                "null_mean": round(float(nb.mean()), 1),
                "null_sd": round(float(nb.std()), 1),
                "z": round(z_of(ob, nb), 2),
                "excess_kept": round(float(ob - nb.mean()), 1),
                "analytic_p": pa,
                "bonferroni_45": min(1.0, pa * 45),
            }
        rn[f"running_mean_{w}d"] = res
    return rn


def main_100s(s: en.Sample, p: en.Params, n: int, t0: float) -> int:
    k = CELLS["100s-1h"]
    rng = np.random.default_rng(5_100_000)
    o = count(s, p, k)
    null = np.array([count(en.scramble_jit(s, rng, p), p, k) for _ in range(n)])
    out: dict = {
        "cell": "CHIME-CHIME 100s-1h wide",
        "n_scrambles": n,
        "full": {
            "obs": o,
            "jit_mean": round(float(null.mean()), 1),
            "z_jit": round(z_of(o, null), 2),
        },
        "d_rate_modulated_null": rate_null_test(s, p, k, n, rng),
        "d_rate_modulated_null_keep_tod": rate_null_test(
            s, p, k, n, rng, windows=(7,), keep_tod=True
        ),
        "runtime_s": None,
    }
    out["runtime_s"] = round(time.time() - t0)
    with open(OUT_100S, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps(out, indent=1))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--cell", choices=tuple(CELLS), default="1h-1d")
    a = ap.parse_args(argv)
    t0 = time.time()
    p = en.Params()
    paths = E.fetch(False)
    df = chime_frame(paths)
    s = en.Sample.from_table(en.chime_events(pd.read_csv(paths["chimefrbcat2.csv"])))
    assert np.allclose(np.sort(df["mjd_400"].to_numpy()), s.mjd)
    if a.cell == "100s-1h":
        return main_100s(s, p, a.n, t0)
    rng = np.random.default_rng(5_000_000)
    out: dict = {"cell": "CHIME-CHIME 1h-1d wide", "n_scrambles": a.n}

    o = count(s, p)
    null = jit_null(s, p, a.n, rng)
    out["full"] = {
        "obs": o,
        "jit_mean": round(float(null.mean()), 1),
        "z_jit": round(z_of(o, null), 2),
    }

    # (a) per-day concentration: pairs with both members on the top x % days (by burst count)
    day = np.floor(s.mjd).astype(int)
    udays, dcount = np.unique(day, return_counts=True)
    order = np.argsort(-dcount, kind="stable")
    i, j = cell_pairs(s, p)
    conc = {}
    scr = [en.scramble_jit(s, rng, p) for _ in range(min(a.n, 300))]
    for frac in (0.01, 0.05, 0.10):
        top = set(udays[order[: max(1, int(round(frac * len(udays))))]].tolist())

        def n_top(x, top=top):
            ii, jj = cell_pairs(x, p)
            d1 = np.floor(x.mjd[ii]).astype(int)
            d2 = np.floor(x.mjd[jj]).astype(int)
            return int(np.count_nonzero(np.isin(d1, list(top)) | np.isin(d2, list(top))))

        ot = n_top(s)
        nt = np.array([n_top(x) for x in scr])
        excess_total = o - null.mean()
        conc[f"top_{int(frac * 100)}pct_days"] = {
            "n_days": len(top),
            "obs_pairs_touching": ot,
            "null_pairs_touching": round(float(nt.mean()), 1),
            "share_of_excess": round(float((ot - nt.mean()) / excess_total), 3),
        }
        # drop the bursts on those days entirely and re-test the cell
        keep = ~np.isin(day, list(top))
        sk = subset(s, keep)
        ok = count(sk, p)
        nk = jit_null(sk, p, 300, rng)
        conc[f"top_{int(frac * 100)}pct_days"].update(
            {"drop_days_n_bursts": int((~keep).sum()), "drop_days_obs": ok,
             "drop_days_jit_mean": round(float(nk.mean()), 1),
             "drop_days_z": round(z_of(ok, nk), 2)}
        )  # fmt: skip
    out["a_day_concentration"] = conc
    # jackknife: drop 30-day blocks
    blocks = (day - day.min()) // 30
    jk = []
    for b in np.unique(blocks):
        m = blocks != b
        sb = subset(s, m)
        ob = count(sb, p)
        nb = jit_null(sb, p, 150, rng)
        jk.append((int(b), int((~m).sum()), round(z_of(ob, nb), 2)))
    zs = np.array([x[2] for x in jk])
    out["a_jackknife_30d"] = {
        "n_blocks": len(jk),
        "z_min": float(zs.min()),
        "z_median": float(np.median(zs)),
        "z_max": float(zs.max()),
        "largest_drop": sorted(jk, key=lambda x: x[2])[:5],
        "note": "block = 30-day span index from the first burst; (block, bursts dropped, z)",
    }

    # (b) epochs and seasons
    from astropy.time import Time

    dt = Time(s.mjd, format="mjd").datetime64
    month = (dt.astype("datetime64[M]").astype(int) % 12) + 1
    splits = {f"year_{y}": s.year == y for y in np.unique(s.year)}
    c1 = df.set_index("tns_name").loc[[*pd.Index(df["tns_name"])]]["catalog1_flag"].to_numpy()
    splits["catalog1_period"] = c1 == 1
    splits["after_catalog1"] = c1 != 1
    for q, ms in enumerate(((12, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11))):
        splits[f"season_{['DJF', 'MAM', 'JJA', 'SON'][q]}"] = np.isin(month, ms)
    ep = {}
    for name, m in splits.items():
        sb = subset(s, m)
        ob = count(sb, p)
        nb = jit_null(sb, p, 200, rng)
        ep[name] = {"n": int(m.sum()), "obs": ob, "jit_mean": round(float(nb.mean()), 1),
                    "z": round(z_of(ob, nb), 2)}  # fmt: skip
    out["b_epochs_seasons"] = ep

    # (c) property independence: observed cell pairs vs scrambled cell pairs
    props = {
        "dm": df["dm_fitb"].to_numpy(float),
        "log_fluence": np.log10(np.clip(df["fluence"].to_numpy(float), 1e-3, None)),
        "snr": df["bonsai_snr"].to_numpy(float),
        "dec": s.dec,
    }
    from scipy import stats

    def diffs(x, ii, jj):
        res = {k: np.abs(v[ii] - v[jj]) for k, v in props.items()}
        lst = np.mod(en.gmst_deg(x.mjd[ii]) - en.gmst_deg(x.mjd[jj]) + 180, 360) - 180
        res["sidereal_phase_deg"] = np.abs(lst)
        return res

    obs_d = diffs(s, i, j)
    # scrambles move times only, so properties stay with their events
    sc_d = {k: [] for k in obs_d}
    for x in scr[:100]:
        ii, jj = cell_pairs(x, p)
        for k, v in diffs(x, ii, jj).items():
            sc_d[k].append(v)
    pr = {}
    for k in obs_d:
        ref = np.concatenate(sc_d[k])
        ov, rv = obs_d[k][np.isfinite(obs_d[k])], ref[np.isfinite(ref)]
        pr[k] = {
            "obs_median": round(float(np.median(ov)), 3),
            "scrambled_median": round(float(np.median(rv)), 3),
            "ks_p": float(stats.ks_2samp(ov, rv).pvalue),
        }
    out["c_property_differences"] = pr

    # (d) rate-modulated null, calibrated by injection
    out["d_rate_modulated_null"] = rate_null_test(s, p, K, a.n, rng)
    out["runtime_s"] = round(time.time() - t0)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
