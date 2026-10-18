"""Achromaticity test of a two-epoch fade (plan: ``results/transients/confirm_r2/B_choice.md``).

Gravitational magnification of a point or uniform-colour source by a compact deflector (point
mass, Ellis wormhole, negative-mass point lens; ``exotic_sim``'s Kitamura-Nakajima-Asada family)
multiplies every band by the same factor, so the epoch-to-epoch flux ratio is band-independent. An
intrinsic transient (e.g. a supernova) changes colour. This script

- ``calibrate``: per band, the matched-control median ratio m_b and robust scatter s_b between the
  two mosaics (controls with S/N >= ``MIN_SNR`` at both epochs) from a ``hit_lightcurve.py`` table;
- ``simulate``: the false-rejection rate (equal ratios) and the power against supernova templates
  (sncosmo, run in an isolated environment) at the real epochs, bands and errors;
- ``test``: the frozen rule on the target.

Every output is ``derived`` (``simulated`` for the simulation); thresholds are ASSUMPTIONs in the
plan.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_dist

MIN_SNR = 5.0  # control S/N at both epochs (plan)
P_REJECT = 0.01  # achromatic rejected below this (plan)
MIN_CONTROLS = 20


def ratio_err(f1, e1, f2, e2):
    """Ratio f2/f1 and its first-order error."""
    r = f2 / f1
    return r, np.abs(r) * np.hypot(e1 / f1, e2 / f2)


def achromatic_chi2(rho, sig):
    """chi^2, dof and p of a constant across bands (weighted mean)."""
    rho, sig = np.asarray(rho, float), np.asarray(sig, float)
    w = 1.0 / sig**2
    mean = (w * rho).sum() / w.sum()
    c = float((w * (rho - mean) ** 2).sum())
    dof = len(rho) - 1
    return c, dof, float(chi2_dist.sf(c, dof)), float(mean)


def calibrate(table, target: str, ep1: dict, ep2: dict) -> dict:
    """``ep1``/``ep2``: {band: obs_id}. Returns per-band control and target ratios."""
    names = np.asarray(table["name"]).astype(str)
    obs = np.asarray(table["obs_id"]).astype(str)
    f, e = np.asarray(table["flux_njy"], float), np.asarray(table["err_njy"], float)
    out = {}
    for b in ep1:

        def get(name, oid):
            k = np.flatnonzero((names == name) & (obs == oid))
            return (f[k[0]], e[k[0]]) if len(k) else (np.nan, np.nan)

        ctrl = []
        for n in sorted(set(names) - {target}):
            f1, e1 = get(n, ep1[b])
            f2, e2 = get(n, ep2[b])
            if np.isfinite([f1, f2, e1, e2]).all() and f1 / e1 >= MIN_SNR and f2 / e2 >= MIN_SNR:
                ctrl.append(f2 / f1)
        ctrl = np.asarray(ctrl)
        med = float(np.median(ctrl)) if len(ctrl) else np.nan
        mad = float(1.4826 * np.median(np.abs(ctrl - med))) if len(ctrl) else np.nan
        f1, e1 = get(target, ep1[b])
        f2, e2 = get(target, ep2[b])
        out[b] = {
            "n_controls": int(len(ctrl)),
            "control_median": med,
            "control_scatter": mad,
            "target": [f1, e1, f2, e2],
        }
    return out


def rule(cal: dict) -> dict:
    """The frozen decision on the target."""
    if any(c["n_controls"] < MIN_CONTROLS for c in cal.values()):
        return {"decision": "not decided (fewer than 20 usable controls in a band)"}
    rho, sig = [], []
    for c in cal.values():
        r, er = ratio_err(*c["target"])
        rho.append(r / c["control_median"])
        sig.append(np.hypot(er, c["control_scatter"]) / c["control_median"])
    chi2, dof, p, mean = achromatic_chi2(rho, sig)
    dec = (
        "chromatic beyond the calibrated noise: achromatic magnification disfavoured"
        if p < P_REJECT
        else "achromatic within errors: not distinguished"
    )
    return {
        "rho": dict(zip(cal, rho, strict=True)),
        "sigma": dict(zip(cal, sig, strict=True)),
        "chi2": chi2,
        "dof": dof,
        "p": p,
        "weighted_mean_ratio": mean,
        "decision": dec,
    }


def simulate(cal: dict, epochs: tuple[float, float], n: int, seed: int, data_dir: str) -> dict:
    """False-rejection rate (equal ratios) and power against SN templates (``simulated``)."""
    import sncosmo

    sncosmo.conf.data_dir = data_dir
    rng = np.random.default_rng(seed)
    bands = list(cal)
    f1 = np.array([cal[b]["target"][0] for b in bands])
    e1 = np.array([cal[b]["target"][1] for b in bands])
    e2 = np.array([cal[b]["target"][3] for b in bands])
    s = np.array([cal[b]["control_scatter"] / cal[b]["control_median"] for b in bands])

    def p_of(t1, t2):
        a1 = t1 + rng.normal(0, e1)
        a2 = t2 + rng.normal(0, e2)
        r, er = ratio_err(a1, e1, a2, e2)
        r = r * (1 + rng.normal(0, s))  # calibration scatter
        return achromatic_chi2(r, np.hypot(er, s))[2]

    rates = {}
    ps = [p_of(f1, f1 * rng.uniform(0.3, 1.5)) for _ in range(n)]
    rates["equal_ratio_false_rejection"] = float(np.mean(np.array(ps) < P_REJECT))
    dt = epochs[1] - epochs[0]
    for src in ("hsiao", "nugent-sn2p", "nugent-sn2l", "nugent-sn1bc"):
        model = sncosmo.Model(source=src)
        ps = []
        while len(ps) < n:
            z = rng.uniform(0.3, 1.5)
            model.set(z=z, t0=0.0, amplitude=1.0)
            ph = rng.uniform(-10, 40)
            t1 = ph * (1 + z)
            lo, hi = model.mintime(), model.maxtime()
            if not (lo <= t1 and t1 + dt <= hi):
                continue
            try:
                m1 = np.array([model.bandflux(b.lower(), t1) for b in bands])
                m2 = np.array([model.bandflux(b.lower(), t1 + dt) for b in bands])
            except ValueError:  # band outside the template's rest range at this z
                continue
            if not (np.all(m1 > 0) and np.all(np.isfinite(m2))):
                continue
            k = f1[bands.index("F150W")] / m1[bands.index("F150W")]
            ps.append(p_of(m1 * k, m2 * k))
        rates[f"power_{src}"] = float(np.mean(np.array(ps) < P_REJECT))
    rates["n_per_class"] = n
    return rates


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("action", choices=["calibrate", "simulate", "test"])
    ap.add_argument("--lc", type=Path, required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--ep1", required=True, help="JSON {band: obs_id}")
    ap.add_argument("--ep2", required=True, help="JSON {band: obs_id}")
    ap.add_argument("--epochs", type=float, nargs=2, default=None, help="MJD of the two epochs")
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--sncosmo-data", default=None)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    from astropy.table import Table

    cal = calibrate(Table.read(a.lc), a.target, json.loads(a.ep1), json.loads(a.ep2))
    if a.action == "calibrate":
        res = {k: {kk: vv for kk, vv in v.items() if kk != "target"} for k, v in cal.items()}
    elif a.action == "simulate":
        res = simulate(cal, tuple(a.epochs), a.n, a.seed, a.sncosmo_data)
    else:
        res = {"calibration": cal, **rule(cal)}
    res = {"label": "simulated" if a.action == "simulate" else "derived", **res}
    a.out.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
