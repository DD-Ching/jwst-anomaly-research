"""D1 distance self-consistency (owner idea 1; hypothesis round 1, docs/hypotheses/round-1/).

Does any single sightline carry two independent distance measures that no ordinary model reconciles?

Lenses (H0LiCOW public posteriors; TDCOSMO 2025 SDSS1206): per lens with a joint (D_d, D_dt)
posterior,
R = D_dt / ((1 + z_d) D_d) = D_s / D_ds (H0-free). Statistics per lens:
  A. prior-predictive pull vs flat LCDM, Omega_m ~ U(0.1, 0.5);
  B. the same in flat wCDM, w ~ U(-2, -0.5);
  C. leave-one-out (LOO) pull vs the other lenses (Omega_m and a common ln-scale free: absorbs a
     sample-wide lambda_MST);
  D. for all 6 lenses, LOO pull of D_dt vs the H0 + Omega_m fit to the other five.
Null: every permutation of the redshift pairs across lenses. Injection: scale one lens's D_dt by f.
FRBs (second probe): DM_obs vs the Macquart-relation predictive distribution per localized host.

    python scripts/d1_distance.py lenses --h0licow DIR --tdcosmo DIR
    python scripts/d1_distance.py tdcosmo --tdcosmo DIR   # statistic D, TDCOSMO 2025 power-law
    python scripts/d1_distance.py frb --frb DIR
Inputs are shallow clones of the three public repos; every file used is pinned by sha256.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import distance_consistency as dc
from jwst_anomaly import paths
from jwst_anomaly.schema import Provenance

OUT = paths.repo_root() / "results" / "d1_distance"
MANIFEST = paths.repo_root() / "data" / "manifests" / "d1_distance.ecsv"
H0LICOW = "https://github.com/shsuyu/H0LiCOW-public"
H0LICOW_COMMIT = "57cf97357e16a94ba0065472bd4185cfbc604aa3"  # 2025-05-14
TDCOSMO = "https://github.com/TDCOSMO/TDCOSMO2025_public"
TDCOSMO_COMMIT = "d7f38db341f68be1df0d9ac1fc528c45113f94cf"  # 2026-01-21
FRBREPO = "https://github.com/FRBs/FRB"
FRB_COMMIT = "996fcda9b0b22431e3171208b4c8e1cf2798e823"  # 2026-05-06
PYGEDM_SDIST = "https://files.pythonhosted.org/packages/source/p/pygedm/pygedm-3.3.0.tar.gz"

# (repo, commit, path) -> sha256
FILES = {
    "HE0435": (
        H0LICOW,
        H0LICOW_COMMIT,
        "h0licow_distance_chains/HE0435_Ddt_AO+HST.dat",
        "f710c5001b2cb0b31909924ac37a5f29298904a79b2c00472df41eb6a3f440da",
    ),
    "J1206": (
        H0LICOW,
        H0LICOW_COMMIT,
        "h0licow_distance_chains/J1206_final.csv",
        "bda60de628658ec43e95a812b8c43b9003acc51ede2b8fb745ce1b3c0acb8dab",
    ),
    "PG1115": (
        H0LICOW,
        H0LICOW_COMMIT,
        "h0licow_distance_chains/PG1115_AO+HST_Dd_Ddt.dat",
        "4f88ce31b7e0c33aee8cc618e8dabd3c2c6fc1e7f458d221a1e11cb5b848091d",
    ),
    "RXJ1131": (
        H0LICOW,
        H0LICOW_COMMIT,
        "h0licow_distance_chains/RXJ1131_AO+HST_Dd_Ddt.dat",
        "9183ff074f0d514dd12c3a3dd94a0a5dc60f96447a21084ff10d67147822f96d",
    ),
    "WFI2033": (
        H0LICOW,
        H0LICOW_COMMIT,
        "h0licow_distance_chains/wfi2033_dt_bic.dat",
        "fb0dc9af101f99ab57bcd1eb3fa38c08937dffa5569479079b6487324e211ad9",
    ),
    "B1608": (
        H0LICOW,
        H0LICOW_COMMIT,
        "MontePython_cosmo_sampling/data/timedelay_6lenses/B1608_Dd_Ddt_params.dat",
        "494da93a660138f6fa083b36915cc669d8c76975dc2597d802fc4fa0cb360028",
    ),
    "H0LiCOW_likelihood": (
        H0LICOW,
        H0LICOW_COMMIT,
        "MontePython_cosmo_sampling/likelihoods/timedelay_6lenses/__init__.py",
        "778c1953a2d626b24b29025f7ec2c6711a498c203ffc8bb2eb0ca9b679e012f6",
    ),
    "TDCOSMO_J1206_Dd": (
        TDCOSMO,
        TDCOSMO_COMMIT,
        "TDCOSMO_sample/TDCOSMO_data/SDSS1206+4332/final_D_d.npy",
        "fe30f78ea15d31f04857fc21ac639a46f7f9362f75b1a3887b1261fdac776666",
    ),
    "TDCOSMO_J1206_Ddt": (
        TDCOSMO,
        TDCOSMO_COMMIT,
        "TDCOSMO_sample/TDCOSMO_data/SDSS1206+4332/final_D_dt.npy",
        "d39ec666058ce285d80ae8e49fe0cde4259ebf0ae152553f74ac11536f8c7592",
    ),
    "TDCOSMO_yaml": (
        TDCOSMO,
        TDCOSMO_COMMIT,
        "TDCOSMO_sample/tdcosmo_sample.yaml",
        "494ca352dccae3a6fe39b6e579c403f0147e07ade603b0bb525d4d268e5ecf60",
    ),
    "FRB_hosts": (
        FRBREPO,
        FRB_COMMIT,
        "frb/data/Galaxies/public_hosts.csv",
        "9349eadc88d1742f5e2a4088b5d07dbc2d71c77c4af22bd49f75a43ef5d1b2d3",
    ),
}
# Redshifts: H0LiCOW MontePython likelihood (timedelay_6lenses/__init__.py);
# identical in TDCOSMO 2025 tdcosmo_sample.yaml.
REDSHIFTS = {
    "B1608": (0.6304, 1.394),
    "RXJ1131": (0.295, 0.654),
    "PG1115": (0.311, 1.722),
    "J1206": (0.745, 1.789),
    "HE0435": (0.4546, 1.693),
    "WFI2033": (0.6575, 1.662),
}
JOINT = ["B1608", "RXJ1131", "PG1115", "J1206"]
ALL = JOINT + ["HE0435", "WFI2033"]
INJ_FACTORS = (
    0.15,
    0.2,
    0.25,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    1.1,
    1.25,
    1.5,
    2.0,
    2.5,
    3.0,
    4.0,
    5.0,
    6.0,
    8.0,
)  # 1.0 = baseline with the injection seed


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pinned(key: str, roots: dict[str, Path], rows: list) -> Path:
    repo, commit, rel, digest = FILES[key]
    path = roots[repo] / rel
    got = sha256(path)
    if got != digest:
        raise SystemExit(f"sha256 mismatch for {path}: {got} != {digest}")
    rows.append((f"{repo}/blob/{commit}/{rel}", digest, path.stat().st_size, commit))
    return path


def load_lenses(roots: dict[str, Path], rng: np.random.Generator, n_keep: int, rows: list) -> dict:
    """{lens: dict(ddt, dd or None, w or None)}: posterior samples (derived)."""
    import pandas as pd

    out = {}
    for name in ("RXJ1131", "PG1115"):  # columns: Dd, Ddt
        arr = pd.read_csv(
            pinned(name, roots, rows), sep=r"\s+", comment="#", header=None
        ).to_numpy()
        out[name] = {"dd": arr[:, 0], "ddt": arr[:, 1], "w": None}
    j = pd.read_csv(pinned("J1206", roots, rows))
    out["J1206"] = {"dd": j["dd"].to_numpy(), "ddt": j["ddt"].to_numpy(), "w": None}
    # TDCOSMO 2025 SDSS1206 files: check they are the same posterior (recorded in the summary).
    td_dd = np.load(pinned("TDCOSMO_J1206_Dd", roots, rows))
    td_ddt = np.load(pinned("TDCOSMO_J1206_Ddt", roots, rows))
    pinned("TDCOSMO_yaml", roots, rows)
    pinned("H0LiCOW_likelihood", roots, rows)  # source of the lens redshifts (REDSHIFTS)
    same = (
        len(td_dd) == len(j)
        and np.allclose(np.sort(td_dd), np.sort(j["dd"]))
        and np.allclose(np.sort(td_ddt), np.sort(j["ddt"]))
    )
    out["J1206"]["tdcosmo_identical"] = bool(same)
    he = pd.read_csv(pinned("HE0435", roots, rows), header=None).to_numpy()[:, 0]
    out["HE0435"] = {"dd": None, "ddt": he, "w": None}
    wf = pd.read_csv(pinned("WFI2033", roots, rows))
    out["WFI2033"] = {"dd": None, "ddt": wf["Dt"].to_numpy(), "w": wf["weight"].to_numpy()}
    # B1608: analytic shifted log-normals (Suyu+2010 D_dt; Jee+2019 D_d), independent.
    b = np.loadtxt(pinned("B1608", roots, rows), comments="#")
    n = 2_000_000
    ddt = b[2] + np.exp(rng.normal(b[3], b[4], n))
    dd = b[5] + np.exp(rng.normal(np.log(b[6]), b[7], n))
    out["B1608"] = {"dd": dd, "ddt": ddt, "w": None, "analytic": True}
    # subsample long chains (keep weights aligned)
    for d in out.values():
        m = len(d["ddt"])
        if m > n_keep:
            idx = rng.choice(m, n_keep, replace=False)
            for k in ("dd", "ddt", "w"):
                if d.get(k) is not None:
                    d[k] = d[k][idx]
    return out


def ratio_stats(
    lens: dict, rng, params: dc.Params, redshifts: dict, scale: dict | None = None
) -> list[dict]:
    """Statistics A, B, C for the joint lenses. scale = {lens: factor on D_dt} for injections."""
    scale = scale or {}
    ln_obs = []
    for name in JOINT:
        zd, _ = redshifts[name]
        r = dc.ratio_from_samples(lens[name]["dd"], lens[name]["ddt"] * scale.get(name, 1.0), zd)
        ln_obs.append(np.log(r))
    zd = np.array([redshifts[n][0] for n in JOINT])
    zs = np.array([redshifts[n][1] for n in JOINT])
    om, grid = dc.grid_ln_model(zd, zs, "ratio", params)
    lns = np.linspace(-params.lns_half_width, params.lns_half_width, params.n_lns)
    loo = dc.leave_one_out(ln_obs, grid, om, lns, rng, params)
    loo_fixed = dc.leave_one_out(ln_obs, grid, om, None, rng, params)
    rows = []
    for k, name in enumerate(JOINT):
        a = dc.pull(ln_obs[k], dc.prior_predictive_ratio(zd[k], zs[k], rng, params), rng)
        b = dc.pull(ln_obs[k], dc.prior_predictive_ratio(zd[k], zs[k], rng, params, wcdm=True), rng)
        q = np.percentile(np.exp(ln_obs[k]), [16, 50, 84])
        rows.append(
            {
                "lens": name,
                "z_d": zd[k],
                "z_s": zs[k],
                "R_p16": q[0],
                "R_p50": q[1],
                "R_p84": q[2],
                "R_lcdm_om0.3": float(dc.ratio_model(zd[k], zs[k], 0.3)[0]),
                "z_A_lcdm": a["z"],
                "zemp_A_lcdm": a["z_emp"],
                "dlnR_A": a["delta_ln"],
                "sig_lnR_A": a["sigma_ln"],
                "z_B_wcdm": b["z"],
                "zemp_B_wcdm": b["z_emp"],
                "z_C_loo": loo[k]["z"],
                "zemp_C_loo": loo[k]["z_emp"],
                "dlnR_C": loo[k]["delta_ln"],
                "sig_lnR_C": loo[k]["sigma_ln"],
                "z_C0_loo_fixed": loo_fixed[k]["z"],
            }
        )
    return rows


def ddt_loo(
    lens: dict,
    rng,
    params: dc.Params,
    redshifts: dict,
    scale: dict | None = None,
    names: list[str] = ALL,
) -> list[dict]:
    """Statistic D: LOO pull of ln D_dt vs H0 + Omega_m fitted to the other lenses."""
    scale = scale or {}
    ln_obs = [np.log(lens[n]["ddt"] * scale.get(n, 1.0)) for n in names]
    w = [lens[n]["w"] for n in names]
    zd = np.array([redshifts[n][0] for n in names])
    zs = np.array([redshifts[n][1] for n in names])
    om, grid = dc.grid_ln_model(zd, zs, "ddt", params)
    ln_s = np.linspace(
        np.log(dc.C_KMS / params.h0_range[1]), np.log(dc.C_KMS / params.h0_range[0]), params.n_lns
    )
    loo = dc.leave_one_out(ln_obs, grid, om, ln_s, rng, params, weights=w)
    return [
        {
            "lens": n,
            "z_D_loo": r["z"],
            "zemp_D_loo": r["z_emp"],
            "dlnDdt_D": r["delta_ln"],
            "sig_lnDdt_D": r["sigma_ln"],
            "H0_loo_others": float(dc.C_KMS / np.exp(r["lns_post_mean"])),
        }
        for n, r in zip(names, loo, strict=True)
    ]


def _crossing(sub: Table, col: str, thr: float, below: bool) -> float | None:
    """Factor nearest 1 at which the injected |z| first reaches thr (None if outside the grid)."""
    f = np.asarray(sub["factor"], float)
    z = np.abs(np.asarray(sub[col], float))
    idx = np.where(f <= 1.0)[0][::-1] if below else np.where(f >= 1.0)[0]
    prev = None
    for i in idx:
        if z[i] >= thr:
            if prev is None:
                return float(f[i])
            lf = np.interp(thr, [z[prev], z[i]], [np.log(f[prev]), np.log(f[i])])
            return round(float(np.exp(lf)), 3)
        prev = i
    return None


def run_lenses(args) -> None:
    params = dc.Params(n_pred=args.n_pred, n_obs=args.n_obs)
    rng = np.random.default_rng(args.seed)
    roots = {H0LICOW: Path(args.h0licow), TDCOSMO: Path(args.tdcosmo)}
    manifest_rows: list = []
    lens = load_lenses(roots, rng, params.n_obs, manifest_rows)
    OUT.mkdir(parents=True, exist_ok=True)

    ratio = ratio_stats(lens, rng, params, REDSHIFTS)
    ddt = ddt_loo(lens, rng, params, REDSHIFTS)
    by = {r["lens"]: r for r in ddt}
    stat_cols = ["z_A_lcdm", "z_B_wcdm", "z_C_loo", "z_D_loo"]
    n_trials = len(JOINT) * 3 + len(
        ALL
    )  # ASSUMPTION: Sidak over every (lens, statistic) pull reported
    thr = dc.local_sigma_threshold(n_trials, params.flag_sigma_global)

    rows = []
    for name in ALL:
        r = next(
            (x for x in ratio if x["lens"] == name),
            {"lens": name, "z_d": REDSHIFTS[name][0], "z_s": REDSHIFTS[name][1]},
        )
        r = {**r, **by[name]}
        zmax = max(
            abs(r.get(c, 0.0)) for c in stat_cols if c in r and np.isfinite(r.get(c, np.nan))
        )
        r["max_abs_z"] = zmax
        r["global_sigma"] = dc.global_sigma(zmax, n_trials)
        r["flag"] = bool(zmax >= thr)
        r["has_Dd"] = name in JOINT
        rows.append(r)
    t = Table(
        rows=[{k: r.get(k, np.nan) for k in rows[0].keys() | set().union(*rows)} for r in rows]
    )
    order = [
        "lens",
        "z_d",
        "z_s",
        "has_Dd",
        "R_p16",
        "R_p50",
        "R_p84",
        "R_lcdm_om0.3",
        "z_A_lcdm",
        "zemp_A_lcdm",
        "z_B_wcdm",
        "zemp_B_wcdm",
        "z_C_loo",
        "zemp_C_loo",
        "z_C0_loo_fixed",
        "dlnR_A",
        "sig_lnR_A",
        "dlnR_C",
        "sig_lnR_C",
        "z_D_loo",
        "zemp_D_loo",
        "dlnDdt_D",
        "sig_lnDdt_D",
        "H0_loo_others",
        "max_abs_z",
        "global_sigma",
        "flag",
    ]
    t = t[order]
    for c in t.colnames:
        if t[c].dtype.kind == "f":
            t[c] = np.round(t[c], 4)
    t.meta["provenance"] = str(Provenance.MODEL_PREDICTION)
    t.meta["source"] = (
        f"H0LiCOW-public@{H0LICOW_COMMIT[:7]} distance posteriors (derived) + TDCOSMO2025_public@"
        f"{TDCOSMO_COMMIT[:7]}; predictions astropy.cosmology flat (w)CDM; scripts/d1_distance.py"
    )
    t.meta["n_trials"] = n_trials
    t.meta["local_sigma_threshold"] = thr
    t.meta["seed"] = args.seed
    t.write(OUT / "lenses.ecsv", overwrite=True)

    # ---- null: every permutation of the redshift pairs across lenses ----
    nparams = dc.Params(n_pred=args.n_null_pred, n_obs=args.n_obs)
    null = {"ratio_max_abs_z_C": [], "ratio_max_abs_z_A": [], "ddt_max_abs_z_D": []}
    perms4 = [p for p in itertools.permutations(range(len(JOINT))) if p != tuple(range(len(JOINT)))]
    for p in perms4:
        z = {JOINT[i]: REDSHIFTS[JOINT[p[i]]] for i in range(len(JOINT))}
        rr = ratio_stats(lens, rng, nparams, {**REDSHIFTS, **z})
        null["ratio_max_abs_z_C"].append(max(abs(x["z_C_loo"]) for x in rr))
        null["ratio_max_abs_z_A"].append(max(abs(x["z_A_lcdm"]) for x in rr))
    perms6 = [p for p in itertools.permutations(range(len(ALL))) if p != tuple(range(len(ALL)))]
    for p in [
        perms6[i] for i in rng.choice(len(perms6), min(args.n_null6, len(perms6)), replace=False)
    ]:
        z = {ALL[i]: REDSHIFTS[ALL[p[i]]] for i in range(len(ALL))}
        dd = ddt_loo(lens, rng, nparams, z)
        null["ddt_max_abs_z_D"].append(max(abs(x["z_D_loo"]) for x in dd))
    obs = {
        "ratio_max_abs_z_C": max(abs(r["z_C_loo"]) for r in ratio),
        "ratio_max_abs_z_A": max(abs(r["z_A_lcdm"]) for r in ratio),
        "ddt_max_abs_z_D": max(abs(r["z_D_loo"]) for r in ddt),
    }
    null_summary = {
        k: {
            "observed": obs[k],
            "n_perm": len(v),
            "null_median": float(np.median(v)),
            "frac_null_ge_observed": float(np.mean(np.array(v) >= obs[k])),
            "null_values": [round(x, 3) for x in v],
        }
        for k, v in null.items()
    }

    # ---- injection: scale one lens's D_dt by f, recover the shift ----
    inj_rows = []
    base_r = {r["lens"]: r for r in ratio}
    base_d = {r["lens"]: r for r in ddt}
    for name in ALL:
        for f in INJ_FACTORS:
            rec = {"lens": name, "factor": f, "ln_factor": float(np.log(f))}
            if name in JOINT:
                rr = {x["lens"]: x for x in ratio_stats(lens, rng, nparams, REDSHIFTS, {name: f})}[
                    name
                ]
                rec.update(
                    {
                        "dlnR_A_recovered": rr["dlnR_A"] - base_r[name]["dlnR_A"],
                        "z_A": rr["z_A_lcdm"],
                        "dlnR_C_recovered": rr["dlnR_C"] - base_r[name]["dlnR_C"],
                        "z_C": rr["z_C_loo"],
                    }
                )
            dd = {x["lens"]: x for x in ddt_loo(lens, rng, nparams, REDSHIFTS, {name: f})}[name]
            rec.update(
                {
                    "dlnDdt_D_recovered": dd["dlnDdt_D"] - base_d[name]["dlnDdt_D"],
                    "z_D": dd["z_D_loo"],
                }
            )
            inj_rows.append(rec)
    it = Table(
        rows=[
            {k: r.get(k, np.nan) for k in inj_rows[0].keys() | set().union(*inj_rows)}
            for r in inj_rows
        ]
    )
    it = it[
        [
            "lens",
            "factor",
            "ln_factor",
            "dlnR_A_recovered",
            "z_A",
            "dlnR_C_recovered",
            "z_C",
            "dlnDdt_D_recovered",
            "z_D",
        ]
    ]
    for c in it.colnames:
        if it[c].dtype.kind == "f":
            it[c] = np.round(it[c], 4)
    it.meta["provenance"] = str(Provenance.SIMULATED)
    it.meta["source"] = (
        "lenses.ecsv inputs, one lens's D_dt samples x `factor` (scripts/d1_distance.py)"
    )

    # injected flag through the whole chain (max |pull| over the statistics vs the threshold)
    zcols = ["z_A", "z_C", "z_D"]
    it["flag"] = [bool(np.nanmax(np.abs([row[c] for c in zcols])) >= thr) for row in it]
    it.write(OUT / "injections.ecsv", overwrite=True)

    # minimum detectable factor per lens and statistic: where the injected |z| crosses thr
    # (interpolated in ln f between the bracketing injections; the baseline pull is included)
    sens = {}
    for name in ALL:
        sub = it[it["lens"] == name]
        sub.sort("factor")
        s = {}
        for c in zcols:
            if np.all(np.isnan(sub[c])):
                continue
            s[c] = {
                "down": _crossing(sub, c, thr, below=True),
                "up": _crossing(sub, c, thr, below=False),
            }
        sens[name] = s

    summary = {
        "provenance": str(Provenance.MODEL_PREDICTION),
        "n_trials": n_trials,
        "local_sigma_threshold": thr,
        "flagged": [r["lens"] for r in rows if r["flag"]],
        "max_abs_z": {r["lens"]: round(float(r["max_abs_z"]), 3) for r in rows},
        "null": null_summary,
        "min_factor_at_threshold_from_injections": sens,
        "j1206_tdcosmo2025_identical_to_h0licow": lens["J1206"]["tdcosmo_identical"],
        "b1608_dd_ddt_correlation": "none (analytic independent fits)",
        "params": {k: v for k, v in vars(params).items()},
    }
    (OUT / "lens_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    write_manifest(manifest_rows, "lenses")
    plot(t, OUT / "lenses.png")
    print(t)
    print(
        json.dumps({k: v for k, v in summary.items() if k != "params"}, indent=1, default=str)[
            :3000
        ]
    )


def plot(t: Table, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    j = t[t["has_Dd"]]
    y = np.arange(len(j))
    ax[0].errorbar(
        j["R_p50"],
        y,
        xerr=[j["R_p50"] - j["R_p16"], j["R_p84"] - j["R_p50"]],
        fmt="o",
        color="k",
        label="observed R (16-84 %)",
    )
    ax[0].plot(j["R_lcdm_om0.3"], y, "s", mfc="none", color="C1", label="flat LCDM, Om = 0.3")
    ax[0].set_yticks(y, j["lens"])
    ax[0].set_xlabel("R = D_dt / ((1+z_d) D_d) = D_s / D_ds")
    ax[0].legend(fontsize=7)
    yy = np.arange(len(t))
    for k, (c, lab) in enumerate(
        (
            ("z_A_lcdm", "A prior LCDM"),
            ("z_B_wcdm", "B prior wCDM"),
            ("z_C_loo", "C LOO R"),
            ("z_D_loo", "D LOO D_dt"),
        )
    ):
        ax[1].plot(t[c], yy + 0.15 * (k - 1.5), "o", label=lab)
    thr = t.meta["local_sigma_threshold"]
    for s in (-thr, thr):
        ax[1].axvline(s, color="r", ls="--", lw=0.8)
    ax[1].axvline(0, color="0.6", lw=0.6)
    ax[1].set_yticks(yy, t["lens"])
    ax[1].set_xlabel("pull (sigma); dashed = 5 sigma trials-corrected")
    ax[1].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=110)


def write_manifest(rows: list, tag: str) -> None:
    import datetime as _dt

    now = _dt.datetime.now(_dt.UTC).strftime("%Y-%m-%d")
    new = Table(
        rows=[(u, s, n, now, c, tag) for u, s, n, c in rows],
        names=("url", "sha256", "size_bytes", "retrieved_utc", "commit", "used_by"),
    )
    if MANIFEST.exists():
        old = Table.read(MANIFEST)
        keep = [i for i, url in enumerate(old["url"]) if url not in set(new["url"])]
        from astropy.table import vstack

        new = vstack([old[keep], new])
    new.meta["provenance"] = str(Provenance.OBSERVED)
    new.meta["source"] = (
        "scripts/d1_distance.py FILES (clones of H0LiCOW-public, TDCOSMO2025_public, FRBs/FRB)"
    )
    new.write(MANIFEST, overwrite=True)


def ymw16_dm_ism(ra_deg, dec_deg, dist_pc: float = dc.FRB_DEFAULT.ymw16_dist_pc) -> np.ndarray:
    """YMW16 (Yao, Manchester & Wang 2017) Milky-Way DM to `dist_pc` along each (ra, dec), pc cm^-3.

    Uses the compiled `ymw16` extension of pygedm directly (`dmdtau`, ndir = 2: distance to DM,
    mode 1: Galactic). The `pygedm` package itself is not imported: its YT2020 halo module calls
    `scipy.integrate.simps`, which recent SciPy removed.
    """
    import importlib.util

    from astropy import units as u
    from astropy.coordinates import SkyCoord

    spec = importlib.util.find_spec("pygedm")
    if spec is None or not spec.submodule_search_locations or not importlib.util.find_spec("ymw16"):
        raise SystemExit(
            "YMW16 needs pygedm's ymw16 extension (results/d1_distance/README.md has a recipe)"
        )
    import ymw16  # compiled extension shipped by pygedm

    datapath = str(spec.submodule_search_locations[0])
    g = SkyCoord(np.atleast_1d(ra_deg) * u.deg, np.atleast_1d(dec_deg) * u.deg).galactic
    return np.array(
        [
            ymw16.dmdtau(float(gl), float(gb), float(dist_pc), 0.0, 2, 1, 0, datapath, "")["DM"]
            for gl, gb in zip(g.l.deg, g.b.deg, strict=True)
        ]
    )


def _js(d: dict) -> str:
    """Canonical form for comparing params dicts (tuples written to ECSV come back as lists)."""
    return json.dumps(d, sort_keys=True, default=str)


def compare_ism(t: Table, ne: Table, model: str, thr: float) -> dict:
    """Join one ISM model's FRB table to the NE2001 one; write `frb_ism_compare.ecsv`.

    A burst counts as flagged if either Galactic model flags it (model error is not averaged)."""
    from astropy.table import join

    keep = ["frb", "z", "dm_obs", "dm_ism_ne2001", "z_low", "z_high", "flag_low", "flag_high"]
    j = join(
        ne[keep],
        t["frb", f"dm_ism_{model}", "z_low", "z_high", "flag_low", "flag_high"],
        keys="frb",
        metadata_conflicts="silent",
        table_names=["ne2001", model],
    )
    r = np.asarray(j[f"dm_ism_{model}"] / j["dm_ism_ne2001"], dtype=float)
    j["ism_ratio"] = np.round(r, 3)
    j["flag_either"] = (
        j["flag_low_ne2001"]
        | j["flag_high_ne2001"]
        | j[f"flag_low_{model}"]
        | j[f"flag_high_{model}"]
    )
    j.sort(f"z_low_{model}", reverse=True)
    j.meta["provenance"] = str(Provenance.MODEL_PREDICTION)
    j.meta["source"] = f"frb.ecsv (NE2001) joined to frb_{model}.ecsv; scripts/d1_distance.py frb"
    j.meta["local_one_sided_sigma_threshold"] = thr
    j.write(OUT / "frb_ism_compare.ecsv", overwrite=True)
    dz = np.maximum(
        np.abs(j[f"z_low_{model}"] - j["z_low_ne2001"]),
        np.abs(j[f"z_high_{model}"] - j["z_high_ne2001"]),
    )
    return {
        "n_joined": len(j),
        "flag_either": [str(x) for x in j["frb"][j["flag_either"]]],
        "ism_ratio_median": round(float(np.median(r)), 4),
        "ism_ratio_16_84": [round(float(x), 4) for x in np.percentile(r, [16, 84])],
        "ism_ratio_range": [round(float(r.min()), 4), round(float(r.max()), 4)],
        "n_ism_ratio_beyond_frac_err": int(np.sum(np.abs(r - 1) > dc.FRB_DEFAULT.ism_frac_err)),
        "max_abs_pull_change": [str(j["frb"][np.argmax(dz)]), float(np.max(dz))],
    }


def run_frb(args) -> None:
    import csv

    root = Path(args.frb)
    rows_m: list = []
    hosts = pinned("FRB_hosts", {FRBREPO: root}, rows_m)
    p = dc.FRBParams()
    inputs, missing = [], []
    with open(hosts, newline="") as f:
        for h in csv.DictReader(f):
            name = h["FRB"].strip()
            try:
                z = float(h["z"])
            except ValueError:
                missing.append((name, "no z"))
                continue
            js = [
                root / "frb/data/FRBs" / f"FRB{name}.json",
                root / "frb/data/FRBs" / f"FRB{name.rstrip('A')}.json",
            ]
            js = [x for x in js if x.exists()]
            if not js:
                missing.append((name, "no json"))
                continue
            d = json.loads(js[0].read_text())
            rows_m.append(
                (
                    f"{FRBREPO}/blob/{FRB_COMMIT}/frb/data/FRBs/{js[0].name}",
                    sha256(js[0]),
                    js[0].stat().st_size,
                    FRB_COMMIT,
                )
            )
            need = ("DM", "DMISM") if args.ism == "ne2001" else ("DM", "ra", "dec")
            if any(k not in d for k in need) or z <= 0:
                missing.append((name, f"no {'/'.join(need)} or z<=0"))
                continue
            ism = d["DMISM"]["value"] if args.ism == "ne2001" else (d["ra"], d["dec"])
            inputs.append((name, z, d["DM"]["value"], ism, h["Projects"]))
    if args.ism == "ymw16":  # one vectorised coordinate transform for all bursts
        ra, dec = np.array([x[3] for x in inputs], dtype=float).T
        dm_ism = ymw16_dm_ism(ra, dec, p.ymw16_dist_pc)
        inputs = [(*x[:3], float(v), x[4]) for x, v in zip(inputs, dm_ism, strict=True)]
    n = len(inputs)
    n_trials = 2 * n  # one low and one high one-sided test per burst
    thr = dc.local_sigma_threshold(n_trials, p.flag_sigma_global, one_sided=True)

    def flags(r):
        return bool(r["z_low"] >= thr), bool(r["z_high"] >= thr)

    ism_col = f"dm_ism_{args.ism}"
    sfx = "" if args.ism == "ne2001" else f"_{args.ism}"
    out, inj = [], []
    for name, z, dm, ism, proj in inputs:
        pred = dc.FRBPredictive(ism, z, p)
        r = dc.frb_predictive_tails(dm, ism, z, p, mean_dm=pred.mean_dm)
        lo, hi = pred.detect_limits(thr)
        fl, fh = flags(r)
        out.append(
            {
                "frb": name,
                "z": z,
                "dm_obs": dm,
                ism_col: ism,
                "projects": proj,
                **r,
                "dm_low_detect": lo,
                "dm_high_detect": hi,
                "flag_low": fl,
                "flag_high": fh,
            }
        )
        # injection through the whole chain: an observed DM just beyond each detection limit
        for side, dm_inj in (("low", 0.9 * lo), ("high", 1.1 * hi)):
            if not np.isfinite(dm_inj) or dm_inj <= 0:
                inj.append({"frb": name, "side": side, "dm_injected": dm_inj, "physical": False})
                continue
            ri = dc.frb_predictive_tails(dm_inj, ism, z, p, mean_dm=pred.mean_dm)
            il, ih = flags(ri)
            inj.append(
                {
                    "frb": name,
                    "side": side,
                    "dm_injected": dm_inj,
                    "physical": True,
                    "z_low": ri["z_low"],
                    "z_high": ri["z_high"],
                    "flag_low": il,
                    "flag_high": ih,
                }
            )
    t = Table(rows=out)
    for c in t.colnames:
        if t[c].dtype.kind == "f" and not c.startswith("p_"):
            t[c] = np.round(t[c], 4)
    t.sort("z_low", reverse=True)
    t.meta["provenance"] = str(Provenance.MODEL_PREDICTION)
    t.meta["source"] = (
        f"FRBs/FRB@{FRB_COMMIT[:7]} public_hosts.csv + frb/data/FRBs/*.json (DM; "
        + (
            "DMISM as stored, NE2001 via frb/mw.py ismDM"
            if args.ism == "ne2001"
            else f"DM_ISM from YMW16 (pygedm ymw16 extension) to {p.ymw16_dist_pc:.0f} pc at ra/dec"
        )
        + "); Macquart+2020 DM_cosmic PDF, log-normal host; grid "
        "convolution (distance_consistency.FRBPredictive); scripts/d1_distance.py frb"
    )
    t.meta["n_trials"] = n_trials
    t.meta["local_one_sided_sigma_threshold"] = thr
    t.meta["skipped"] = [f"{a}: {b}" for a, b in missing]
    t.meta["params"] = {k: v for k, v in vars(p).items()}
    t.write(OUT / f"frb{sfx}.ecsv", overwrite=True)

    it = Table(
        rows=[{k: r.get(k, np.nan) for k in inj[-1].keys() | set().union(*inj)} for r in inj]
    )
    it = it["frb", "side", "dm_injected", "physical", "z_low", "z_high", "flag_low", "flag_high"]
    it.meta["provenance"] = str(Provenance.SIMULATED)
    it.meta["source"] = (
        f"frb{sfx}.ecsv inputs with DM_obs replaced by 0.9 x dm_low_detect / 1.1 x dm_high_detect"
    )
    it.meta["local_one_sided_sigma_threshold"] = thr
    it.write(OUT / f"frb_injections{sfx}.ecsv", overwrite=True)
    if args.ism == "ne2001":
        write_manifest(rows_m, "frb")
    else:  # same FRB inputs; record only the YMW16 code
        sd = Path(args.pygedm_sdist) if args.pygedm_sdist else None
        if sd and sd.exists():
            write_manifest([(PYGEDM_SDIST, sha256(sd), sd.stat().st_size, "3.3.0")], "frb_ymw16")

    phys = it[it["physical"]]
    lo_i, hi_i = phys[phys["side"] == "low"], phys[phys["side"] == "high"]
    low_ok = np.array([x > 0 for x in t["dm_low_detect"]])
    summary = {
        "provenance": str(Provenance.MODEL_PREDICTION),
        "n_frb": n,
        "n_trials": n_trials,
        "local_one_sided_sigma_threshold": thr,
        "flag_low": [str(x) for x in t["frb"][t["flag_low"]]],
        "flag_high": [str(x) for x in t["frb"][t["flag_high"]]],
        "max_z_low": [str(t["frb"][0]), float(t["z_low"][0])],
        "max_z_high": [str(t["frb"][np.argmax(t["z_high"])]), float(np.max(t["z_high"]))],
        "low_flag_reachable_dm_gt_0": int(low_ok.sum()),
        "median_dm_low_detect_over_dm_ism": float(
            np.median(np.asarray(t["dm_low_detect"])[low_ok] / np.asarray(t[ism_col])[low_ok])
        ),
        "median_dm_high_detect": float(np.nanmedian(t["dm_high_detect"])),
        "injections_low_recovered": f"{int(lo_i['flag_low'].sum())}/{len(lo_i)}",
        "injections_high_recovered": f"{int(hi_i['flag_high'].sum())}/{len(hi_i)}",
        "injections_wrong_side_flags": int(lo_i["flag_high"].sum() + hi_i["flag_low"].sum()),
    }
    summary["ism_model"] = args.ism
    if args.ism != "ne2001":
        ne = Table.read(OUT / "frb.ecsv") if (OUT / "frb.ecsv").exists() else None
        ne_params = dict(ne.meta.get("params", {})) if ne is not None else {}
        cur = dict(t.meta["params"])
        if ne is None:
            reason = "frb.ecsv missing: run the ne2001 path first"
        elif sorted(ne["frb"]) != sorted(t["frb"]) or _js(ne_params) != _js(cur):
            reason = "frb.ecsv has other bursts or params: re-run the ne2001 path first"
        else:
            reason = None
        if reason:
            print("vs_ne2001 skipped:", reason)
            summary["vs_ne2001"] = f"skipped: {reason}"
        else:
            summary["vs_ne2001"] = compare_ism(t, ne, args.ism, thr)
    (OUT / f"frb_summary{sfx}.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(t["frb", "z", "dm_obs", ism_col, "dm_pred_median", "z_low", "z_high"][:8])
    print(json.dumps(summary, indent=1))


# ---- TDCOSMO 2025 power-law D_dt chains (D-073 addendum 2: statistic D on 8 lenses) ----
_TD = "TDCOSMO_sample/TDCOSMO_data/"
J1206_KHIST = (
    "kappahist_J1206_measured_5innermask_nobeta_zgap-1.0_-1.0_fiducial_120_gal_120_zoverr_45_gal_45_"
    "zoverr_24_med_increments4_4_4_4.cat"
)
WGD2038_KHIST = (
    "kappahist_2038_measured_3innermask_nobeta_removehandpicked_zgap-1.0_-1.0_fiducial_120_gal_120_"
    "oneoverr_45_gal_45_zoverr_22.5_med_increments2_2_2_2_emptymsk.cat"
)
# lens -> (D_dt^model file, kappa_ext file, sha256 of each); the pairing is the one in TDCOSMO's
# tdcosmo_sample.ipynb (the power-law models behind the TDCOSMO 2025 likelihoods).
TDCOSMO_PL = {
    "B1608": (
        "B1608+656/B1608_Dtmod_n5e5.dat",
        "049903c193a9cd8d588f99d0dcc722073cbe3645553435d67d56153d58223a19",
        "B1608+656/B1608_kext.txt",
        "eab36e676bc1ef629d293b7ef24f25376c2f4f9cf60a0cd709c0731b7431833d",
    ),
    "RXJ1131": (
        "RXJ1131-1231/rxj_powerlaw_Ddt.dat",
        "1c715fe33eb4ca8dd6eace556198e96d794d5f53a01a173c95f489023a47cfff",
        "RXJ1131-1231/kappa_powerlaw_rxj.dat",
        "9e5d313f6634072bf724e7ec72ad7a799d647f35186db98d0c1ac0d8969682c5",
    ),
    "PG1115": (
        "PG1115+080/pg_powerlaw_Ddt.dat",
        "3b51740e01d62ea63d5c2345a75476569a1f7c4b875790d614767689bcd0f7b8",
        "PG1115+080/kappa_powerlaw_pg.dat",
        "7a6bb6ef03d5b9bf1a3c3cee83c4baf5f84ec22e2ef4e527e2264f6422398fd7",
    ),
    "J1206": (
        "SDSS1206+4332/angular_diameter_pre_LOS_power_law.txt",
        "32bd7ec755a7e4c128d14cc56d5f7738881df5911e2f67ead2e5d66802c315cd",
        f"SDSS1206+4332/{J1206_KHIST}",
        "bc498dbebbc99f03520963e01921ae49f829585fe44367732a178fdc7d046760",
    ),
    "HE0435": (
        "HE0435-1223/he_powerlaw_Ddt.dat",
        "891cec93d84277a5abd3fef84bbc46e00921dad61520231a97eb115f74245e1c",
        "HE0435-1223/kappa_powerlaw_he.dat",
        "d4d09087f7b213532368ad3606fbf6f423886f071afb9331ad4e2437b6eaed11",
    ),
    "WFI2033": (
        "WFI2033-4723/wfi2033_pl_dt_nokext.dat",
        "27de2eb954b29b631e81a85eea1b5ac11f3c34baddbe976a6e3e25538e8443cc",
        "WFI2033-4723/wfi2033_kext_bic_pdf.dat",
        "faf3eb5165f8a72453d7a01ec552491675b0c3cdd27b1fd234a338efd7c6c5b8",
    ),
    "DES0408": (
        "DES0408-5354/power_law_dist_post_no_kext.txt",
        "1d7167d7cc695eb406062db0a741cdb47879ce6a8e8b53df932067de326e9e5d",
        "DES0408-5354/kext_sampled.txt",
        "3b08ec4e566daac686665639d05a5a818e5f84b9686df6f7236c131594aef048",
    ),
    "WGD2038": (
        "WGD2038-4008/desj2038_pl_nokext_nokin_dt_weight.csv",
        "ac816f55de9b4b1d41635c80a4d1a6db3acf1f0bc2bc6c74772a912796931d53",
        f"WGD2038-4008/{WGD2038_KHIST}",
        "b7288742f9c596570d947e223dc64fb4d76e34dc3966f6a851b054c524b79ca7",
    ),
}
# Redshifts from TDCOSMO_sample/tdcosmo_sample.yaml (the 6 H0LiCOW lenses agree with REDSHIFTS).
TD_REDSHIFTS = {**REDSHIFTS, "DES0408": (0.597, 2.375), "WGD2038": (0.2283, 0.777)}
TD_ALL = ALL + ["DES0408", "WGD2038"]
KEXT_VARIANTS = ("kext", "nokext")  # per-lens TDCOSMO kappa_ext PDFs, or kappa_ext = 0 everywhere


def _pinned_td(root: Path, rel: str, digest: str, rows: list) -> Path:
    path = root / _TD / rel
    got = sha256(path)
    if got != digest:
        raise SystemExit(f"sha256 mismatch for {path}: {got} != {digest}")
    rows.append(
        (f"{TDCOSMO}/blob/{TDCOSMO_COMMIT}/{_TD}{rel}", digest, path.stat().st_size, TDCOSMO_COMMIT)
    )
    return path


def load_tdcosmo_pl(root: Path, rng: np.random.Generator, n_keep: int, rows: list) -> dict:
    """{lens: {"ddt_model", "kappa"}}: n_keep equal-weight draws each (derived posteriors).

    Readers follow tdcosmo_sample.ipynb: single-column chains with a header line; WFI2033 and
    WGD2038 carry weights; B1608 is (weight, Dt^mod, gamma', kappa); SDSS1206 is a pickle of
    (D_d D_s / D_ds, D_d, kappa_pert), D_dt = (1 + z_d) x first / (1 + kappa_pert); the J1206 and
    WGD2038 kappa_ext files are histograms on [-0.1, 1] and [-0.2, 1]."""
    import pandas as pd

    def col(path, skip=1, sep=","):
        return pd.read_csv(path, sep=sep, header=None, skiprows=skip, comment="#").to_numpy()

    out = {}
    for name, (fd, dd_sha, fk, k_sha) in TDCOSMO_PL.items():
        pd_ = _pinned_td(root, fd, dd_sha, rows)
        pk = _pinned_td(root, fk, k_sha, rows)
        w = None
        if name == "B1608":
            a = col(pd_, skip=1, sep=r"\s+")
            ddt, w = a[:, 1], a[:, 0]
        elif name == "J1206":
            first, _, kpert = dc.load_array_pickle(pd_)
            ddt = (1 + TD_REDSHIFTS[name][0]) * first / (1 + kpert)
            ddt = ddt[ddt > 0]
        elif name in ("WFI2033", "WGD2038"):
            a = col(pd_)
            ddt, w = a[:, 0], a[:, 1]
        elif name == "DES0408":
            ddt = col(pd_, skip=0)[:, 0]
        else:
            ddt = col(pd_)[:, 0]
        if name in ("J1206", "WGD2038"):
            pdf = np.loadtxt(pk, usecols=[0], comments="#")
            lo = -0.1 if name == "J1206" else -0.2
            kappa = dc.draw_from_hist(pdf, np.linspace(lo, 1.0, len(pdf) + 1), n_keep, rng)
        elif name == "WFI2033":
            a = col(pk)
            kappa = dc.resample(a[:, 0], n_keep, rng, a[:, 1])
        elif name in ("B1608", "DES0408"):
            kappa = dc.resample(col(pk, skip=0)[:, 0], n_keep, rng)
        else:
            kappa = dc.resample(col(pk)[:, 0], n_keep, rng)
        out[name] = {
            "ddt_model": dc.resample(ddt, n_keep, rng, w),
            "kappa": kappa,
            "n_chain": int(len(ddt)),
        }
    return out


def _td_lens(td: dict, variant: str) -> dict:
    """lens dict for ddt_loo: D_dt with the variant's kappa_ext treatment, equal weights."""
    return {
        n: {
            "ddt": d["ddt_model"]
            if variant == "nokext"
            else dc.ddt_with_kext(d["ddt_model"], d["kappa"]),
            "w": None,
        }
        for n, d in td.items()
    }


def run_tdcosmo(args) -> None:
    params = dc.Params(n_pred=args.n_pred, n_obs=args.n_obs)
    nparams = dc.Params(n_pred=args.n_null_pred, n_obs=args.n_obs)
    rng = np.random.default_rng(args.seed)
    rows_m: list = []
    td = load_tdcosmo_pl(Path(args.tdcosmo), rng, params.n_obs, rows_m)
    OUT.mkdir(parents=True, exist_ok=True)
    n_trials = len(TD_ALL) * len(KEXT_VARIANTS)  # ASSUMPTION: Sidak over (lens, variant) pulls
    thr = dc.local_sigma_threshold(n_trials, params.flag_sigma_global)

    rows, base, null, sens = [], {}, {}, {}
    inj_rows = []
    for v in KEXT_VARIANTS:
        lens = _td_lens(td, v)
        res = ddt_loo(lens, rng, params, TD_REDSHIFTS, names=TD_ALL)
        base[v] = {r["lens"]: r for r in res}
        for r in res:
            q = np.percentile(lens[r["lens"]]["ddt"], [16, 50, 84])
            k = np.percentile(td[r["lens"]]["kappa"], [16, 50, 84]) if v == "kext" else (0, 0, 0)
            rows.append(
                {
                    "lens": r["lens"],
                    "variant": v,
                    "z_d": TD_REDSHIFTS[r["lens"]][0],
                    "z_s": TD_REDSHIFTS[r["lens"]][1],
                    "Ddt_p16": q[0],
                    "Ddt_p50": q[1],
                    "Ddt_p84": q[2],
                    "kext_p16": k[0],
                    "kext_p50": k[1],
                    "kext_p84": k[2],
                    "n_chain": td[r["lens"]]["n_chain"],
                    **{c: r[c] for c in r if c != "lens"},
                    "global_sigma": dc.global_sigma(r["z_D_loo"], n_trials),
                    "flag": bool(abs(r["z_D_loo"]) >= thr),
                }
            )
        # null: random permutations of the redshift pairs across the 8 lenses
        perms = [p for p in itertools.permutations(range(len(TD_ALL))) if list(p) != sorted(p)]
        vals = []
        for i in rng.choice(len(perms), args.n_null, replace=False):
            p = perms[i]
            z = {TD_ALL[j]: TD_REDSHIFTS[TD_ALL[p[j]]] for j in range(len(TD_ALL))}
            vals.append(
                max(abs(x["z_D_loo"]) for x in ddt_loo(lens, rng, nparams, z, names=TD_ALL))
            )
        obs = max(abs(r["z_D_loo"]) for r in res)
        null[v] = {
            "observed_max_abs_z": round(obs, 3),
            "n_perm": len(vals),
            "null_median": round(float(np.median(vals)), 3),
            "frac_null_ge_observed": float(np.mean(np.array(vals) >= obs)),
        }
        # injection: one lens's D_dt x f through the same statistic
        for name in TD_ALL:
            for f in INJ_FACTORS:
                dd = {
                    x["lens"]: x
                    for x in ddt_loo(lens, rng, nparams, TD_REDSHIFTS, {name: f}, names=TD_ALL)
                }[name]
                inj_rows.append(
                    {
                        "lens": name,
                        "variant": v,
                        "factor": f,
                        "dlnDdt_D_recovered": dd["dlnDdt_D"] - base[v][name]["dlnDdt_D"],
                        "z_D": dd["z_D_loo"],
                        "flag": bool(abs(dd["z_D_loo"]) >= thr),
                    }
                )
    t = Table(rows=rows)
    it = Table(rows=inj_rows)
    for tab in (t, it):
        for c in tab.colnames:
            if tab[c].dtype.kind == "f":
                tab[c] = np.round(tab[c], 4)
    src = (
        f"TDCOSMO2025_public@{TDCOSMO_COMMIT[:7]} power-law D_dt chains, kappa_ext PDFs "
        "(derived); "
        "flat LCDM predictions (astropy.cosmology); scripts/d1_distance.py tdcosmo"
    )
    t.meta.update(
        {
            "provenance": str(Provenance.MODEL_PREDICTION),
            "source": src,
            "n_trials": n_trials,
            "local_sigma_threshold": thr,
            "seed": args.seed,
        }
    )
    it.meta.update(
        {
            "provenance": str(Provenance.SIMULATED),
            "source": "tdcosmo_lenses.ecsv inputs, one lens's D_dt x `factor`",
        }
    )
    t.write(OUT / "tdcosmo_lenses.ecsv", overwrite=True)
    it.write(OUT / "tdcosmo_injections.ecsv", overwrite=True)
    for v in KEXT_VARIANTS:
        sv = it[it["variant"] == v]
        sens[v] = {}
        for name in TD_ALL:
            sub = sv[sv["lens"] == name]
            sub.sort("factor")
            sens[v][name] = {
                "down": _crossing(sub, "z_D", thr, below=True),
                "up": _crossing(sub, "z_D", thr, below=False),
            }
    summary = {
        "provenance": str(Provenance.MODEL_PREDICTION),
        "n_trials": n_trials,
        "local_sigma_threshold": round(thr, 3),
        "flagged": [f"{r['lens']}/{r['variant']}" for r in rows if r["flag"]],
        "max_abs_z": {f"{r['lens']}/{r['variant']}": round(float(r["z_D_loo"]), 3) for r in rows},
        "H0_loo_others": {
            f"{r['lens']}/{r['variant']}": round(float(r["H0_loo_others"]), 1) for r in rows
        },
        "null": null,
        "min_factor_at_threshold_from_injections": sens,
        "params": {k: v for k, v in vars(params).items()},
    }
    (OUT / "tdcosmo_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    write_manifest(rows_m, "tdcosmo")
    print(t["lens", "variant", "Ddt_p50", "kext_p50", "z_D_loo", "H0_loo_others", "flag"])
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k not in ("params", "max_abs_z")},
            indent=1,
            default=str,
        )
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = paths.data_root() / "d1"
    a = sub.add_parser("lenses")
    a.add_argument("--h0licow", default=str(d / "H0LiCOW-public"))
    a.add_argument("--tdcosmo", default=str(d / "TDCOSMO2025_public"))
    a.add_argument("--n-pred", type=int, default=200_000)
    a.add_argument("--n-obs", type=int, default=200_000)
    a.add_argument("--n-null-pred", type=int, default=50_000)
    a.add_argument("--n-null6", type=int, default=100)
    a.add_argument("--seed", type=int, default=20261009)
    c = sub.add_parser("tdcosmo")
    c.add_argument("--tdcosmo", default=str(d / "TDCOSMO2025_public"))
    c.add_argument("--n-pred", type=int, default=200_000)
    c.add_argument("--n-obs", type=int, default=200_000)
    c.add_argument("--n-null-pred", type=int, default=50_000)
    c.add_argument("--n-null", type=int, default=100)
    c.add_argument("--seed", type=int, default=20261009)
    b = sub.add_parser("frb")
    b.add_argument("--frb", default=str(d / "FRB"))
    b.add_argument("--seed", type=int, default=20261009)
    b.add_argument(
        "--ism",
        choices=("ne2001", "ymw16"),
        default="ne2001",
        help="Galactic DM model: ne2001 (DMISM stored in the FRB repo) or ymw16 (needs pygedm)",
    )
    b.add_argument(
        "--pygedm-sdist", help="pygedm 3.3.0 sdist used for --ism ymw16 (manifest sha256)"
    )
    args = ap.parse_args()
    {"lenses": run_lenses, "tdcosmo": run_tdcosmo, "frb": run_frb}[args.cmd](args)


if __name__ == "__main__":
    main()
