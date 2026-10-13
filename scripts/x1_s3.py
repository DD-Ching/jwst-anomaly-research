"""X1 S3: the 2018-07-22 GRAVITY flare (CDS J/A+A/618/L10 avg.dat) vs far-side and ordinary spots.

Exploratory (docs/state/units/X1.md "Independence of data": no GRAVITY flare is held out).
Frozen before this comparison (X1.md, S2): a far-side spot at a = 2.1 M has a loop radius
<= 8.5 uas, the ordinary spot >= 14.5 uas (same static-light nuisance).

Data (processed numerical values): 10 epochs, t [min], x, y [mas] offsets from Sgr A* with errors,
flux relative to S2. Model per epoch, with the observed flux F(t) (intrinsic evolution included):
  c(t) = c0 + R_theta H [ (F(t) - S) c_spot(phi0 + 2 pi t / P) + S d ] / F(t)
where c_spot is the delayed k = 0 total-light centroid (scripts/x1_tracer.py), H = diag(1, +-1) the
sense of rotation, R_theta the sky position angle, S = fs * mean(F) static flux (fs 0-0.5) at
offset d (0-20 uas). Free: c0 (analytic, weighted), theta, sense, phi0, fs, d. Only the astrometry
enters chi^2; the flux contrast is not tested here because one flare over ~half an orbit cannot
separate orbital modulation from intrinsic evolution (data class insufficient for C).

  python scripts/x1_s3.py   -> results/x1/s3_avg_fit.json

Every number is derived (our computation) under unestablished premises (X1.md); a good or bad fit
is not evidence for or against a wormhole.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import urllib.request
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh_farside as wf  # noqa: E402
import x1_s2  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "x1"
URL = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/618/L10/avg.dat"
SHA256 = "7d1729d345ea397d3c09223fcbc6777b9c675f614e66490c8c10c649f5deb05f"
PERIODS = (57.0, 60.0, 63.0)  # min, P_pol = 60 +- 3 (premise P_pol = P_orb)
# ASSUMPTION: nuisance grids (S2 ranges) and fit grids
FS = (0.0, 0.25, 0.5)
D_UAS = (0.0, 10.0, 20.0)
PSI = np.linspace(0.0, 2 * math.pi, 8, endpoint=False)
N_GRID = 72  # theta and phi0 grid points


def data_path() -> Path:
    root = Path(os.environ.get("JWST_ANOMALY_DATA", ROOT / "data"))
    p = root / "raw" / "x1" / "avg.dat"
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(URL, timeout=60) as r:
            p.write_bytes(r.read())
    digest = hashlib.sha256(p.read_bytes()).hexdigest()
    if digest != SHA256:
        raise RuntimeError(f"{p}: sha256 {digest} != pinned {SHA256}")
    return p


def load() -> dict:
    a = np.loadtxt(data_path())
    return {
        "t": a[:, 0],
        "xy": a[:, [1, 3]] * 1e3,  # uas
        "err": a[:, [2, 4]] * 1e3,
        "flux": a[:, 5],
    }


def best_chi2(track: np.ndarray, d: dict, period: float, d_uas=D_UAS) -> dict:
    """Minimum chi^2 over theta, sense, phi0 (grids), fs, d, psi (grids); centre analytic."""
    ph = track[:, 0]
    cs = track[:, 1:3] * x1_s2.THETA_G
    w = 1.0 / d["err"] ** 2
    f = d["flux"]
    best = {"chi2": math.inf}
    phis = np.linspace(0.0, 2 * math.pi, N_GRID, endpoint=False)
    thetas = np.linspace(0.0, 2 * math.pi, N_GRID, endpoint=False)
    # spot centroid at each epoch for every phi0: (n_phi0, n_epoch, 2), periodic interpolation
    arg = (phis[:, None] + 2 * math.pi * d["t"][None, :] / period) % (2 * math.pi)
    pp = np.concatenate([ph, ph + 2 * math.pi])
    spot = np.stack(
        [np.interp(arg, pp, np.tile(cs[:, k], 2)) for k in (0, 1)], axis=-1
    )  # (P, E, 2)
    rot = np.stack(
        [np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]]) for t in thetas]
    )  # (T, 2, 2)
    for fs in FS:
        s = fs * f.mean()
        if np.any(f - s <= 0):
            continue
        for dd in d_uas:
            for psi in PSI if dd > 0 else (0.0,):
                off = dd * np.array([math.cos(psi), math.sin(psi)])
                for sense in (1.0, -1.0):
                    m = ((f - s)[None, :, None] * spot + s * off) / f[None, :, None]
                    m = m * np.array([1.0, sense])
                    mod = np.einsum("tij,pej->tpei", rot, m)  # (T, P, E, 2)
                    res = d["xy"][None, None] - mod
                    c0 = (res * w).sum(2, keepdims=True) / w.sum(0)
                    chi2 = (((res - c0) ** 2) * w).sum((2, 3))
                    k = np.unravel_index(np.argmin(chi2), chi2.shape)
                    if chi2[k] < best["chi2"]:
                        best = {
                            "chi2": float(chi2[k]),
                            "theta_deg": math.degrees(thetas[k[0]]),
                            "phi0_deg": math.degrees(phis[k[1]]),
                            "sense": sense,
                            "static_fraction": fs,
                            "static_offset_uas": dd,
                            "static_psi_deg": math.degrees(psi),
                            "centre_uas": c0[k[0], k[1], 0].tolist(),
                        }
    return best


def main() -> int:
    d = load()
    w = 1.0 / d["err"] ** 2
    c0 = (d["xy"] * w).sum(0) / w.sum(0)
    chi2_static = float((((d["xy"] - c0) ** 2) * w).sum())
    n = d["xy"].size
    # free: centre (2), theta, phi0, sense, fs, d, psi -> count 7 continuous-ish; ASSUMPTION
    dof = n - 7
    rows = []
    for p in PERIODS:
        rk = wf.kepler_radius(p, x1_s2.M_SUN)
        models = [("ordinary", 0.0, e, e ** (2 / 3) * rk, "near") for e in x1_s2.ETAS]
        models += [("far", a, 1.0, rk, "far") for a in x1_s2.A_FAR if rk > a]
        for kind, a, eta, re, side in models:
            tr, nw = x1_s2.track(a, re, side, eta)
            b = best_chi2(tr, d, p)
            rows.append(
                {"model": kind, "a_M": a, "eta": eta, "P_min": p, "R_M": re, "warnings": nw}
                | b
                | {"p_value": float(stats.chi2.sf(b["chi2"], dof))}
            )
    by_a = {}
    for r in rows:
        if r["model"] == "far":
            o = min(
                q["chi2"] for q in rows if q["model"] == "ordinary" and q["P_min"] == r["P_min"]
            )
            by_a.setdefault(r["a_M"], []).append(r["chi2"] - o)
    delta = {str(a): min(v) for a, v in by_a.items()}
    # outside the pre-written nuisance range: off-centre static light is the known escape (#172)
    rk = wf.kepler_radius(60.0, x1_s2.M_SUN)
    sens = []
    for dmax in (40.0, 60.0, 100.0):
        grid = (0.0, 0.5 * dmax, dmax)
        far = best_chi2(x1_s2.track(2.1, rk, "far", 1.0)[0], d, 60.0, grid)
        ordi = best_chi2(x1_s2.track(0.0, rk, "near", 1.0)[0], d, 60.0, grid)
        sens.append(
            {
                "static_offset_max_uas": dmax,
                "far_a2.1_chi2": far["chi2"],
                "ordinary_chi2": ordi["chi2"],
            }
            | {
                "far_best_offset_uas": far["static_offset_uas"],
                "far_best_fs": far["static_fraction"],
            }
        )
    out = {
        "labels": "derived (our computation); data: processed numerical values (CDS"
        " J/A+A/618/L10 avg.dat); exploratory, no held-out flare",
        "data": {"url": URL, "sha256": SHA256, "n_epochs": int(len(d["t"]))},
        "chi2_static_point": chi2_static,
        "p_static_point": float(stats.chi2.sf(chi2_static, n - 2)),
        "dof_models": dof,
        "min_delta_chi2_far_minus_best_ordinary_same_P": delta,
        "note_error_scaling": "delta chi2 scales as 1/k^2 if all errors are underestimated by k",
        "sensitivity_static_offset_P60": sens,
        "fits": rows,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s3_avg_fit.json").write_text(json.dumps(out, indent=1) + "\n")
    for r in rows:
        print(
            f"{r['model']:8s} a={r['a_M']:.1f} eta={r['eta']:.1f} P={r['P_min']:.0f} "
            f"chi2={r['chi2']:.1f} p={r['p_value']:.2g} fs={r['static_fraction']}"
        )
    print("static point chi2", chi2_static, "p", out["p_static_point"])
    print("min delta chi2 far - ordinary", delta)
    print(json.dumps(sens, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
