"""X4 S1b step D: precessing, higher-mode ordinary templates (IMRPhenomXPHM), representative cases.

Unit file: docs/state/units/X4.md (S1b; frozen before computing). `simulated` / `derived`; no data.
Needs lalsuite.

Target: IMRPhenomXPHM (q 1.5, aligned spins 0.2, inclination 60 deg, reference phase 0), detector
response h_plus (polarisation angle 0), multiplied by the negative-lens factor F_neg.

Ordinary templates: IMRPhenomXPHM with chirp mass, q <= 15, both spin vectors (Cartesian,
|chi| <= 0.99), inclination, reference phase and polarisation angle psi
(h = cos 2psi h_plus + sin 2psi h_cross); arrival time and amplitude maximised; NO free overall
phase (the real part of the correlation is maximised, because with higher modes a constant phase
rotation is not a physical change).
  UA  aligned spins only (higher modes, no precession);
  UP  full precessing spins;
  PP  UP x positive point lens: waveform fitted with the lens fixed at step B's best lens, then the
      lens polished with the waveform fixed; two alternations.

  python scripts/x4_s1b_xphm.py            # -> results/x4/s1b_xphm.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


c = _load("x4_s1b_checks")
b = c.b
INC_TARGET = np.pi / 3
D_CASES = ("R2", "R3", "R4")  # R1 is undetermined after step C (P fit-error check failed)
PAD = 32
CHI_MAX, Q_MAX, MC_FRAC = 0.99, 15.0, 0.60


def xphm(mc, q, s1, s2, inc, phi, psi):
    import lal
    import lalsimulation as ls

    m1, m2 = b.masses(mc, q)
    hp, hc = ls.SimInspiralChooseFDWaveform(
        m1 * lal.MSUN_SI,
        m2 * lal.MSUN_SI,
        *s1,
        *s2,
        1e8 * lal.PC_SI,
        inc,
        phi,
        0,
        0,
        0,
        b.DF,
        b.F_LOW,
        b.F_HIGH,
        b.F_LOW,
        lal.CreateDict(),
        ls.GetApproximantFromString("IMRPhenomXPHM"),
    )
    n = min(hp.data.length, b._FREQS.size)
    h = np.zeros(b._FREQS.size, complex)
    h[:n] = np.cos(2 * psi) * np.asarray(hp.data.data[:n]) + np.sin(2 * psi) * np.asarray(
        hc.data.data[:n]
    )
    h[~b._BAND] = 0.0
    return h


def overlap_real(a, t):
    """max over arrival time of Re<a, t> / norms (amplitude free, no free overall phase)."""
    w = b.weight()
    prod = t * np.conj(a) * w
    n = PAD * prod.size
    corr = np.fft.ifft(prod, n) * n
    na, nt = np.sum(np.abs(a) ** 2 * w), np.sum(np.abs(t) ** 2 * w)
    if na <= 0 or nt <= 0:
        return 0.0
    return float(np.max(corr.real) / np.sqrt(na * nt))


def unpack(x, aligned):
    if aligned:
        lmc, q, c1z, c2z, inc, phi, psi = x
        s1, s2 = (0.0, 0.0, c1z), (0.0, 0.0, c2z)
    else:
        lmc, q, a1x, a1y, a1z, a2x, a2y, a2z, inc, phi, psi = x
        s1, s2 = (a1x, a1y, a1z), (a2x, a2y, a2z)
    return np.exp(lmc), q, s1, s2, inc, phi, psi


def obj(x, tgt, mc0, aligned, lensf=None):
    mc, q, s1, s2, inc, phi, psi = unpack(x, aligned)
    if not (
        abs(mc / mc0 - 1) <= MC_FRAC
        and 1.0 <= q <= Q_MAX
        and 0.0 <= inc <= np.pi
        and np.linalg.norm(s1) <= CHI_MAX
        and np.linalg.norm(s2) <= CHI_MAX
    ):
        return 1.0
    try:
        h = xphm(mc, q, s1, s2, inc, phi, psi)
    except Exception:  # noqa: BLE001
        return 1.0
    if lensf is not None:
        h = h * lensf
    return 1.0 - overlap_real(h, tgt)


def nm(fun, x0, steps, maxfev=1500, restarts=3):
    x0 = np.asarray(x0, float)
    best_f, best_x, hist, nfev = 1.0, x0, [], 0
    x = x0
    for k in range(restarts + 1):
        simplex = np.vstack([x] + [x + np.eye(x.size)[i] * steps[i] for i in range(x.size)])
        r = minimize(
            fun,
            x,
            method="Nelder-Mead",
            options={
                "maxfev": maxfev if k == 0 else 400,
                "xatol": 1e-6,
                "fatol": 1e-9,
                "initial_simplex": simplex,
            },
        )
        nfev += r.nfev
        improve = best_f - float(r.fun)
        if r.fun < best_f:
            best_f, best_x = float(r.fun), r.x
        hist.append(best_f)
        x = best_x
        if k > 0 and improve < c.TOL:
            return best_f, best_x, hist, nfev, True
    return best_f, best_x, hist, nfev, False


A_STEPS = np.array([0.03, 0.5, 0.2, 0.2, 0.3, 0.5, 0.3])
P_STEPS = np.array([0.03, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.5, 0.3])


def target(case):
    mtot, ml, y = c.CASES[case]
    mc0 = b.chirp(mtot, b.Q_INJ)
    h = xphm(mc0, b.Q_INJ, (0, 0, b.CHI_INJ), (0, 0, b.CHI_INJ), INC_TARGET, 0.0, 0.0)
    return mc0, h, h * b.lens(ml, y, -1)


def starts(case, aligned, rng):
    mc0 = b.chirp(c.CASES[case][0], b.Q_INJ)
    base = [np.log(mc0), b.Q_INJ, b.CHI_INJ, b.CHI_INJ, INC_TARGET, 0.0, 0.0]
    out = [base]
    fits = json.loads((b.OUT / "s1b_fits.json").read_text())
    ub = max(
        (r for r in fits["U"] if r["case"] == case and r["q_max"] == Q_MAX),
        key=lambda r: r["FF_pad64"],
    )
    out.append([*ub["x"], INC_TARGET, 0.0, 0.0])
    for _ in range(3):
        out.append(
            [
                np.log(mc0 * rng.uniform(0.7, 1.3)),
                rng.uniform(1, 12),
                *rng.uniform(-0.7, 0.7, 2),
                rng.uniform(0.2, 2.9),
                rng.uniform(0, 2 * np.pi),
                rng.uniform(0, np.pi),
            ]
        )
    if aligned:
        return out
    full = []
    for s in out:
        lmc, q, c1, c2, inc, phi, psi = s
        inpl = rng.uniform(-0.3, 0.3, 4) if s is not out[0] else np.zeros(4)
        full.append([lmc, q, inpl[0], inpl[1], c1, inpl[2], inpl[3], c2, inc, phi, psi])
    return full


def job(args):
    case, model, x0, lensp = args
    t1 = time.time()
    mc0, h_unl, tgt = target(case)
    if model == "sanity":  # unlensed XPHM target fitted by UP: must reach FF ~ 1
        f, x, hist, nfev, conv = nm(lambda z: obj(z, h_unl, mc0, False), x0, P_STEPS)
    elif model in ("UA", "UP"):
        aligned = model == "UA"
        f, x, hist, nfev, conv = nm(
            lambda z: obj(z, tgt, mc0, aligned), x0, A_STEPS if aligned else P_STEPS
        )
    else:  # PP: alternate waveform (lens fixed) and lens (waveform fixed)
        lml, y = np.log(lensp[0]), lensp[1]
        x, f, hist, nfev, conv = np.asarray(x0, float), 1.0, [], 0, False
        for _ in range(2):
            lf = b.lens_exact(np.exp(lml), y, +1)
            f, x, h1, n1, conv = nm(
                lambda z, lf=lf: obj(z, tgt, mc0, False, lf), x, P_STEPS, maxfev=1200, restarts=2
            )
            nfev += n1

            def lobj(v, x=x):
                if not (0.0 <= v[0] <= np.log(1e4) and 0.01 <= v[1] <= 5.0):
                    return 1.0
                return obj(x, tgt, mc0, False, b.lens_exact(np.exp(v[0]), v[1], +1))

            r = minimize(
                lobj,
                [lml, y],
                method="Nelder-Mead",
                options={
                    "maxfev": 40,
                    "xatol": 1e-4,
                    "fatol": 1e-7,
                    "initial_simplex": [[lml, y], [lml + 0.1, y], [lml, y + 0.05]],
                },
            )
            nfev += r.nfev
            lml, y = float(r.x[0]), float(r.x[1])
            f = min(f, float(r.fun))
            hist.append(f)
        lensp = (float(np.exp(lml)), y)
    mc, q, s1, s2, inc, phi, psi = unpack(x, model == "UA")
    return {
        "case": case,
        "model": model,
        "FF": 1.0 - f,
        "x": [float(v) for v in x],
        "lens": lensp,
        "q": float(q),
        "chi1": [float(v) for v in s1],
        "chi2": [float(v) for v in s2],
        "inclination": float(inc),
        "nfev": int(nfev),
        "history": hist,
        "converged_restart": bool(conv),
        "q_on_bound": bool(q >= 0.99 * Q_MAX),
        "seconds": round(time.time() - t1, 1),
    }


def main(argv=None) -> int:
    rng = np.random.default_rng(31)
    fits = json.loads((b.OUT / "s1b_fits.json").read_text())
    jobs = []
    # cases that passed the numerical stage, plus the absorbed control (X4.md S1b, frozen order)
    for case in D_CASES:
        for s in starts(case, True, rng):
            jobs.append((case, "UA", s, None))
        for s in starts(case, False, rng):
            jobs.append((case, "UP", s, None))
        pb = max((r for r in fits["P"] if r["case"] == case), key=lambda r: r["FF_pad64"])
        for s in starts(case, False, rng)[:2]:
            jobs.append((case, "PP", s, (pb["M_Lz"], pb["y"])))
    # tool sanity: an unlensed XPHM target must be recovered (FF ~ 1) from a perturbed start
    for case in ("R2", "R3"):
        mc0 = b.chirp(c.CASES[case][0], b.Q_INJ)
        jobs.append(
            (
                case,
                "sanity",
                [np.log(mc0) + 0.02, 2.0, 0.1, -0.1, 0.1, 0.0, 0.1, 0.3, 1.2, 0.5, 0.3],
                None,
            )
        )
    t0 = time.time()
    with Pool(4) as pool:
        rows = []
        for r in pool.imap_unordered(job, jobs):
            rows.append(r)
            print(
                r["case"],
                r["model"],
                round(r["FF"], 5),
                "q",
                round(r["q"], 2),
                "nfev",
                r["nfev"],
                "conv",
                r["converged_restart"],
                f"({r['seconds']} s)",
                flush=True,
            )
    best = {}
    for r in rows:
        k = f"{r['case']}:{r['model']}"
        if k not in best or r["FF"] > best[k]["FF"]:
            best[k] = r
    rec = {
        "label": "simulated; IMRPhenomXPHM target and templates (S1b step D)",
        "wall_s": round(time.time() - t0, 1),
        "rows": rows,
        "best": best,
    }
    (b.OUT / "s1b_xphm.json").write_text(json.dumps(rec, indent=1, default=float) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
