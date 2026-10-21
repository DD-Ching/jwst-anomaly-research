"""S1 burst twins on CHIME/FRB Catalog 2: driver of the <= 30-event validation (not a search).

Frozen plan: ``results/s1_chime/PLAN.md``; rules and adapter: ``scripts/s1_chime_lib.py``.

  python -I scripts/s1_chime.py select --cat <chimefrbcat2.csv> --listing <listing.tsv>
  python -I scripts/s1_chime.py fetch  --raw <empty dir> --log <downloads.log>  (stream, delete)
  python -I scripts/s1_chime.py run    [--plots <dir>]

``select`` writes ``sample.ecsv``; ``fetch`` writes ``profiles.json`` (reduced sub-band profiles,
derived) and ``data/manifests/s1_chime_dynspec.ecsv``; ``run`` writes ``pairs.ecsv``, ``null.json``,
``injections.ecsv.gz``, ``summary.json`` and ``contact_sheet.png``. A match is an anomaly, never
evidence of new physics.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import hashlib
import io
import json
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import s1_chime_lib as lib  # noqa: E402

OUT = ROOT / "results" / "s1_chime"
MANIFEST = ROOT / "data" / "manifests" / "s1_chime_dynspec.ecsv"
BASE = "https://cadc-west-01.canfar.net/vault/files/AstroDataCitationDOI/CISTI.CANFAR/25.0066/data/dynamic_spectra/data/"
P = lib.Params()
SEED = 20261021


def write_json(obj, path: Path) -> None:
    path.write_text(json.dumps(obj, indent=1, default=float) + "\n")


def write_ecsv_gz(t: Table, path: Path) -> None:
    buf = io.StringIO()
    t.write(buf, format="ascii.ecsv")
    with open(path, "wb") as fh, gzip.GzipFile(fileobj=fh, mode="wb", mtime=0) as gz:
        gz.write(buf.getvalue().encode())


# ------------------------------------------------------------------------------------- select


def cmd_select(a) -> None:
    import pandas as pd

    cat = pd.read_csv(a.cat, low_memory=False)
    lst = pd.read_csv(a.listing, sep="\t", header=None, names=["file", "size", "md5"])
    s = lib.select_sample(cat, set(lst.file))
    s["file"] = s.tns_name + "_stokesi_dynamic_spectrum.h5"
    s = s.merge(lst[["file", "size"]], on="file", how="left")
    t = Table.from_pandas(
        s[["tns_name", "group", "source", "repeater_name", "n_comp", "snr", "file", "size"]].fillna(
            {"repeater_name": ""}
        )
    )
    t.meta["provenance"] = "observed"
    t.meta["source"] = (
        "CHIME/FRB Catalog 2 table (doi:10.11570/25.0066, sha256 5108ada7...) and the VOSpace "
        "listing of data/dynamic_spectra/data"
    )
    t.meta["rule"] = lib.select_sample.__doc__.split("\n\n")[1].strip()
    OUT.mkdir(parents=True, exist_ok=True)
    t.write(OUT / "sample.ecsv", format="ascii.ecsv", overwrite=True)
    print(t)
    print("total bytes", int(np.sum(t["size"])))


# ------------------------------------------------------------------------------------- fetch


def _fetch_one(args):
    name, raw, log = args
    dest = Path(raw) / name
    url = BASE + name
    if dest.exists():  # left by an interrupted run: reuse if its bytes match the logged sha256
        h = hashlib.sha256(dest.read_bytes()).hexdigest()
        for line in Path(log).read_text().splitlines():
            w = line.split("\t")
            if len(w) >= 4 and w[1] == url and w[3] == h:
                return name, {"size": int(w[2]), "sha256": h, "when": w[0], "dest": str(dest)}, None
    for attempt in range(3):
        r = subprocess.run(
            ["curl", "-sSL", "--fail", "--max-time", "600", "-o", str(dest), url],
            capture_output=True,
            text=True,
        )
        if r.returncode == 0:
            break
        time.sleep(10 * 3**attempt)
    else:
        return name, None, f"curl failed: {r.stderr.strip()[:200]}"
    h = hashlib.sha256(dest.read_bytes()).hexdigest()
    size = dest.stat().st_size
    when = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(log, "a") as fh:
        fh.write(
            f"{when}\t{url}\t{size}\t{h}\t"
            "S1 CHIME frozen sample (results/s1_chime/PLAN.md); reduced then deleted\n"
        )
    return name, {"size": size, "sha256": h, "when": when, "dest": str(dest)}, None


def _pack(red: dict) -> dict:
    keep = {
        k: red[k]
        for k in (
            "tns_name",
            "crop",
            "region",
            "n_time",
            "dm_incoherent",
            "res_time",
            "beam_number",
            "cal_source",
            "has_calibration",
        )
    }
    keep["band_ok"] = red["band_ok"].tolist()
    keep["n_valid_chan"] = red["n_valid_chan"].tolist()
    keep["cal_shape"] = [
        None if not np.isfinite(x) else round(float(x), 5) for x in red["cal_shape"]
    ]
    keep["err"] = [float(f"{x:.5g}") for x in red["err"][:, 0]]
    keep["flux"] = [[float(f"{x:.4g}") for x in row] for row in red["flux"]]
    keep["model"] = [[float(f"{x:.4g}") for x in row] for row in red["model"]]
    return keep


def unpack(d: dict) -> dict:
    flux = np.array(d["flux"], float)
    err = np.repeat(np.array(d["err"], float)[:, None], flux.shape[1], axis=1)
    return {
        "tns_name": d["tns_name"],
        "flux": flux,
        "model": np.array(d["model"], float),
        "err": err,
        "band_ok": np.array(d["band_ok"], bool),
        "crop": tuple(d["crop"]),
        "region": tuple(d["region"]),
        "n_time": int(d["n_time"]),
        "cal_shape": np.array([np.nan if x is None else x for x in d["cal_shape"]], float),
    }


def cmd_fetch(a) -> None:
    sample = Table.read(OUT / "sample.ecsv", format="ascii.ecsv")
    raw = Path(a.raw)
    raw.mkdir(parents=True, exist_ok=True)
    names = list(sample["file"])
    out, rows, fails = {}, [], []
    t0 = time.time()
    with cf.ThreadPoolExecutor(4) as ex:
        for name, info, err in ex.map(_fetch_one, [(n, raw, a.log) for n in names]):
            if info is None:
                fails.append((name, err))
                print("FAIL", name, err, flush=True)
                continue
            red = lib.reduce_file(info["dest"], P)
            if not a.keep_raw:
                Path(info["dest"]).unlink()
            out[red["tns_name"]] = _pack(red)
            rows.append((BASE + name, name, info["size"], info["sha256"], info["when"]))
            print(
                f"{name} {info['size'] / 1e6:.1f} MB, nt={red['n_time']}, crop={red['crop']}",
                flush=True,
            )
    print(
        f"{len(rows)} files, {sum(r[2] for r in rows) / 1e6:.1f} MB in {time.time() - t0:.0f} s; "
        f"{len(fails)} failed"
    )
    write_json(
        {
            "provenance": "derived",
            "source": "CHIME/FRB Catalog 2 dynamic spectra; scripts/s1_chime_lib.reduce_file",
            "units": "Jy per sub-band (mean over valid channels); err = off-pulse 1.4826 MAD",
            "bands_mhz": P.bands_mhz,
            "res_time_s": lib.RES_TIME_S,
            "failed": fails,
            "bursts": out,
        },
        OUT / "profiles.json",
    )
    m = Table(rows=rows, names=("data_uri", "file_name", "size", "sha256", "retrieved_utc"))
    m["pipeline_version"] = "CHIME/FRB Catalog 2 release (doi:10.11570/25.0066)"
    m.meta["provenance"] = "observed"
    m.meta["source"] = (
        "CHIME/FRB Catalog 2 Stokes-I dynamic spectra (CANFAR VOSpace, route of issue #153)"
    )
    m.meta["reason"] = (
        "S1 CHIME validation (<= 30 events, results/s1_chime/PLAN.md); streamed, reduced to "
        "results/s1_chime/profiles.json, deleted"
    )
    m.sort("file_name")
    m.write(MANIFEST, format="ascii.ecsv", overwrite=True)


# ------------------------------------------------------------------------------------- run


def load_profiles():
    prof = json.loads((OUT / "profiles.json").read_text())
    sample = Table.read(OUT / "sample.ecsv", format="ascii.ecsv")
    reds, meta = [], []
    for row in sample:
        d = prof["bursts"].get(str(row["tns_name"]))
        if d is None:
            continue
        reds.append(unpack(d))
        meta.append(
            {
                "tns_name": str(row["tns_name"]),
                "source": str(row["source"]),
                "group": str(row["group"]),
            }
        )
    return reds, meta


def cmd_run(a) -> None:
    reds, meta = load_profiles()
    n = len(reds)
    npul = [lib.n_pulses(r, P) for r in reds]
    elig = [lib.is_eligible(r, P) for r in reds]
    ei = [i for i in range(n) if elig[i]]
    print(f"{n} reduced, {len(ei)} eligible (>= {P.min_pulses} pulses at 5 sigma)")
    pairs = [(i, j) for k, i in enumerate(ei) for j in ei[k + 1 :]]
    chance = [(i, j) for i, j in pairs if meta[i]["source"] != meta[j]["source"]]

    # --- surrogate null for rho* (D-071 rule, on the different-source pairs) ---
    rng = np.random.default_rng(SEED)
    maxes = (
        np.array([lib.surrogate_rho_max(reds, chance, rng) for _ in range(P.n_surrogates)])
        if chance
        else np.array([])
    )
    rho_star = float(np.quantile(maxes, P.rho_star_q)) if maxes.size else float("nan")
    write_json(
        {
            "provenance": "derived",
            "rule": "rho* = 0.95 quantile of the per-catalogue max rho over 100 pulse-shuffled "
            "surrogate catalogues (different-source eligible pairs)",
            "n_surrogates": int(maxes.size),
            "per_catalogue_max": [round(float(x), 5) for x in maxes],
            "rho_star": rho_star,
        },
        OUT / "null.json",
    )
    print(f"rho* = {rho_star:.4f}")

    # --- real pairs: chance (ii) and same-source (iii) ---
    rows = []
    for i, j in pairs:
        rho, lag = lib.pair_rho(reds[i], reds[j])
        c = lib.chi2_union(reds[i], reds[j], lag, P)
        c2, d2, p2 = lib.chi2_d071(reds[i], reds[j], lag)
        flag = bool(rho > rho_star)
        rows.append(
            {
                "a": meta[i]["tns_name"],
                "b": meta[j]["tns_name"],
                "kind": "same_source"
                if meta[i]["source"] == meta[j]["source"]
                else "different_source",
                "source_a": meta[i]["source"],
                "source_b": meta[j]["source"],
                "rho": rho,
                "lag": lag,
                "flag": flag,
                "testable": c["testable"],
                "chi2": c["chi2"],
                "dof": c["dof"],
                "p": c["p"],
                "n_signal_bands": c["n_signal_bands"],
                "gains_ok": c["gains_ok"],
                "chi2_pass": c["passed"],
                "passed": bool(flag and c["passed"]),
                "p_d071": p2,
                "gain_snr": ",".join(
                    "nan" if not np.isfinite(x) else f"{x:.1f}" for x in c["gain_snr"]
                ),
            }
        )
    pt = Table(rows=rows)
    pt.meta["provenance"] = "derived"
    pt.meta["source"] = "results/s1_chime/profiles.json; chain of scripts/s1_chime_lib.py"
    pt.write(OUT / "pairs.ecsv", format="ascii.ecsv", overwrite=True)

    # --- injections (i): every ordered (A eligible, host B != A) x set x ratio, one draw each ---
    cal = np.array([r["cal_shape"] for r in reds])
    irng = np.random.default_rng(SEED + 1)
    # deviation (post hoc, before any outcome): a host with no usable band (no calibration group)
    # cannot carry a copy in Jy, so it is not used as a host
    hosts = [i for i in range(n) if reds[i]["band_ok"].any()]
    inj = []
    t0 = time.time()
    for kind in lib.GAIN_SETS:
        for ratio in P.ratios:
            for ia in ei:
                for ib in hosts:
                    if ib == ia:
                        continue
                    g = lib.draw_gains(kind, irng, cal, P)
                    sh = float(irng.uniform(-P.shift_max_bins, P.shift_max_bins))
                    cp = lib.make_copy(reds[ia], reds[ib], ratio, g, sh, kind == "P_path", P)
                    r = lib.run_chain(reds[ia], cp, rho_star, P)
                    inj.append(
                        {
                            "set": kind,
                            "ratio": ratio,
                            "a": meta[ia]["tns_name"],
                            "host": meta[ib]["tns_name"],
                            "shift": sh,
                            "true_lag": cp["true_lag"],
                            "lag": r["lag"],
                            "rho": r["rho"],
                            "elig_b": r["elig_b"],
                            "flag": r["flag"],
                            "chi2_testable": r["chi2"]["testable"],
                            "chi2_p": r["chi2"]["p"],
                            "gains_ok": r["chi2"]["gains_ok"],
                            "passed": r["passed"],
                            "g": ",".join(f"{x:.3f}" for x in g),
                        }
                    )
    print(f"{len(inj)} injections in {time.time() - t0:.0f} s")
    it = Table(rows=inj)
    it.meta["provenance"] = "simulated"
    it.meta["source"] = (
        "copies of fitburst models in other sample bursts' real residual noise "
        "(scripts/s1_chime_lib.make_copy)"
    )
    write_ecsv_gz(it, OUT / "injections.ecsv.gz")

    # --- summary and verdict ---
    def frac(m):
        return float(np.mean(m)) if len(m) else float("nan")

    eff = {}
    for kind in lib.GAIN_SETS:
        for ratio in P.ratios:
            m = (it["set"] == kind) & (it["ratio"] == ratio)
            eff[f"{kind}_r{ratio}"] = {
                "n": int(m.sum()),
                "eligible_copy": frac(it["elig_b"][m]),
                "flagged": frac(it["flag"][m]),
                "passed": frac(it["passed"][m]),
                "chi2_pass_given_flag": frac(it["passed"][m & it["flag"]]),
                "rho_median": float(np.median(it["rho"][m])),
                "lag_within_2": frac(np.abs(it["lag"][m] - it["true_lag"][m]) <= 2),
            }
    ch = pt[pt["kind"] == "different_source"]
    rp = pt[pt["kind"] == "same_source"]
    q = float(np.quantile(ch["rho"], P.t3_quantile)) if len(ch) else float("nan")
    m1 = (it["set"] == "G_cal") & (it["ratio"] == 1.0)
    vin = {
        "n_eligible": len(ei),
        "eff_gcal_r1": eff["G_cal_r1.0"]["passed"],
        "n_chance_pairs": len(ch),
        "n_chance_pass": int(np.sum(ch["passed"])) if len(ch) else 0,
        "rho_inj_gcal_r1_median": float(np.median(it["rho"][m1])),
        "rho_chance_q": q,
        "n_rep_pairs": len(rp),
        "n_rep_pass": int(np.sum(rp["passed"])) if len(rp) else 0,
    }
    v = lib.verdict(vin, P)

    def side(tab):
        if not len(tab):
            return {"n": 0}
        return {
            "n": len(tab),
            "rho_median": float(np.median(tab["rho"])),
            "rho_max": float(np.max(tab["rho"])),
            "n_flag": int(np.sum(tab["flag"])),
            "n_chi2_pass_any_rho": int(np.sum(tab["chi2_pass"])),
            "n_passed_chain": int(np.sum(tab["passed"])),
            "n_untestable": int(np.sum(~tab["testable"])),
        }

    cal_ok = cal[np.isfinite(cal).all(1)]
    summary = {
        "provenance": "derived (real pairs) and simulated (injections)",
        "plan": "results/s1_chime/PLAN.md",
        "n_reduced": n,
        "n_eligible": len(ei),
        "per_burst": [
            dict(meta[i], n_pulses=npul[i], eligible=elig[i], band_ok=reds[i]["band_ok"].tolist())
            for i in range(n)
        ],
        "rho_star": rho_star,
        "chance_pairs": side(ch),
        "same_source_pairs": side(rp),
        "same_source_by_repeater": {
            s: side(rp[rp["source_a"] == s]) for s in sorted(set(rp["source_a"]))
        }
        if len(rp)
        else {},
        "injections": eff,
        "cal_shape_rms_per_band": float(np.sqrt(np.mean(cal_ok**2))) if len(cal_ok) else None,
        "cal_shape_pair_diff_rms_per_band": float(
            np.sqrt(np.mean([(x - y) ** 2 for k, x in enumerate(cal_ok) for y in cal_ok[k + 1 :]]))
        )
        if len(cal_ok) > 1
        else None,
        "verdict_inputs": vin,
        "verdict": v,
    }
    write_json(summary, OUT / "summary.json")
    print(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "n_eligible",
                    "rho_star",
                    "chance_pairs",
                    "same_source_pairs",
                    "verdict_inputs",
                    "verdict",
                )
            },
            indent=1,
            default=float,
        )
    )
    if a.plots is not None:
        contact_sheet(reds, meta, elig, pt, OUT / "contact_sheet.png")


def contact_sheet(reds, meta, elig, pt, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    idx = {m["tns_name"]: k for k, m in enumerate(meta)}
    order = np.argsort(-np.asarray(pt["rho"]))
    show = []
    for kind in ("different_source", "same_source"):
        sel = [k for k in order if pt["kind"][k] == kind][:6]
        show += sel
    nb = len(reds)
    ncol = 6
    nrow_b = int(np.ceil(nb / ncol))
    nrow_p = int(np.ceil(len(show) / ncol))
    fig, axes = plt.subplots(nrow_b + nrow_p, ncol, figsize=(3 * ncol, 2.1 * (nrow_b + nrow_p)))
    axes = np.atleast_2d(axes)
    cols = ["#1b5e9b", "#3a9b3a", "#c9861a", "#b0302a"]
    for k in range(nrow_b * ncol):
        ax = axes[k // ncol, k % ncol]
        if k >= nb:
            ax.axis("off")
            continue
        f, e = lib.crop(reds[k])
        for b in range(f.shape[0]):
            if reds[k]["band_ok"][b]:
                ax.plot(f[b] / max(e[b, 0], 1e-30) + 6 * b, color=cols[b], lw=0.8)
        ax.set_title(
            f"{meta[k]['tns_name']} {'E' if elig[k] else '-'} "
            f"{meta[k]['source'][-8:] if meta[k]['group'] == 'repeater' else ''}",
            fontsize=7,
        )
        ax.tick_params(labelsize=6)
    for q, k in enumerate(show):
        ax = axes[nrow_b + q // ncol, q % ncol]
        ra, rb = reds[idx[str(pt["a"][k])]], reds[idx[str(pt["b"][k])]]
        lag = int(pt["lag"][k])
        tb = np.arange(rb["n_time"])
        ta = np.arange(ra["n_time"]) + lag
        for b in range(4):
            if ra["band_ok"][b] and rb["band_ok"][b]:
                sa, sb = ra["flux"][b] / ra["err"][b, 0], rb["flux"][b] / rb["err"][b, 0]
                ax.plot(tb, sb + 8 * b, color=cols[b], lw=0.8)
                ax.plot(
                    ta,
                    sa / max(sa.max(), 1e-9) * max(sb.max(), 1e-9) + 8 * b,
                    color="k",
                    lw=0.5,
                    alpha=0.6,
                )
        lo = min(ra["crop"][0] + lag, rb["crop"][0]) - 5
        hi = max(ra["crop"][1] + lag, rb["crop"][1]) + 5
        ax.set_xlim(lo, hi)
        ax.set_title(
            f"{pt['kind'][k][:4]} {pt['a'][k][3:]}/{pt['b'][k][3:]} "
            f"rho={pt['rho'][k]:.2f} p={pt['p'][k]:.1e}",
            fontsize=6.5,
        )
        ax.tick_params(labelsize=6)
    for q in range(len(show), nrow_p * ncol):
        axes[nrow_b + q // ncol, q % ncol].axis("off")
    fig.suptitle(
        "S1 CHIME sample: crop sub-band S/N (400-500 blue ... 700-800 red); "
        "bottom: top-rho pairs (B colour, A black, scaled)",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=55)
    plt.close(fig)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--cat", required=True)
    s.add_argument("--listing", required=True)
    f = sub.add_parser("fetch")
    f.add_argument("--raw", required=True)
    f.add_argument("--log", required=True)
    f.add_argument("--keep-raw", action="store_true", help="keep raw files (delete by hand)")
    r = sub.add_parser("run")
    r.add_argument("--plots", default=None)
    a = ap.parse_args(argv)
    {"select": cmd_select, "fetch": cmd_fetch, "run": cmd_run}[a.cmd](a)


if __name__ == "__main__":
    main()
