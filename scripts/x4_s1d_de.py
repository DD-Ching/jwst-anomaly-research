"""X4 S1d: R1 with a reference global search (SciPy differential evolution), XPHM family only.

Unit file: docs/state/units/X4.md (S1d; reuse check and plan frozen before computing).
`simulated` / `derived`; no data. Needs lalsuite. One fitter serves validation and case fits.

Score: the S1c family-D score (noise-weighted real-part overlap, time maximised) with the
polarisation angle and amplitude maximised analytically: per lag, max_{a,b} of the overlap of
a h+ + b hx is sqrt(v^T G^-1 v) / |t|, v = (Re<t,h+>, Re<t,hx>), G the real 2x2 Gram matrix
(the polarisation maximisation used in precessing searches; equal to a free psi in [0, pi)).

  python scripts/x4_s1d_de.py validate   # -> results/x4/s1d_validate.json
  python scripts/x4_s1d_de.py case       # -> results/x4/s1d_case.json
  python scripts/x4_s1d_de.py summary    # -> results/x4/s1d_summary.json
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import scipy.fft as sf
from scipy.optimize import differential_evolution, minimize

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


m = _load("x4_s1c_r1")
b, d = m.b, m.d
MC0 = m.MC0
TOL = 1e-5
BOUNDS = [
    (np.log(0.4 * MC0), np.log(1.6 * MC0)),
    (1.0, 15.0),
    (0.0, 0.99),
    (-1.0, 1.0),
    (0.0, 2 * np.pi),
    (0.0, 0.99),
    (-1.0, 1.0),
    (0.0, 2 * np.pi),
    (0.0, np.pi),
    (0.0, 2 * np.pi),
]
STEPS = np.array([0.03, 0.5, 0.1, 0.2, 0.5, 0.1, 0.2, 0.5, 0.3, 0.5])
DE_U = {"popsize": 12, "maxiter": 200}
DE_P = {"popsize": 12, "maxiter": 100}
OUT = b.OUT


# ------------------------------------------------------------------ waveforms and score
def to_cart(x):
    """Box parameters -> (mc, q, s1, s2, inc, phi); spins as magnitude, cos tilt, azimuth."""
    out = []
    for a, ct, ph in (x[2:5], x[5:8]):
        st = np.sqrt(max(0.0, 1.0 - ct * ct))
        out.append((a * st * np.cos(ph), a * st * np.sin(ph), a * ct))
    return float(np.exp(x[0])), float(x[1]), out[0], out[1], float(x[8]), float(x[9])


def from_cart(mc, q, s1, s2, inc, phi):
    x = [np.log(mc), q]
    for s in (s1, s2):
        a = float(np.linalg.norm(s))
        x += [
            min(a, 0.99),
            (s[2] / a) if a > 0 else 0.0,
            float(np.arctan2(s[1], s[0]) % (2 * np.pi)),
        ]
    return np.array(x + [inc, phi % (2 * np.pi)], float)


def hpc(mc, q, s1, s2, inc, phi):
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
    a = np.zeros(b._FREQS.size, complex)
    c = np.zeros(b._FREQS.size, complex)
    a[:n] = hp.data.data[:n]
    c[:n] = hc.data.data[:n]
    a[~b._BAND] = 0.0
    c[~b._BAND] = 0.0
    return a, c


def stat(t, hp, hc, pad=32):
    """Overlap maximised over time, polarisation angle and amplitude (no free overall phase)."""
    w = b.weight()
    n = sf.next_fast_len(pad * t.size)
    vp = (sf.ifft(t * np.conj(hp) * w, n) * n).real
    vc = (sf.ifft(t * np.conj(hc) * w, n) * n).real
    g11 = np.sum(np.abs(hp) ** 2 * w)
    g22 = np.sum(np.abs(hc) ** 2 * w)
    g12 = np.sum((hp * np.conj(hc)).real * w)
    det = g11 * g22 - g12 * g12
    nt = np.sum(np.abs(t) ** 2 * w)
    if det <= 1e-12 * g11 * g22 or nt <= 0:  # h+ and hx parallel: one polarisation only
        if g11 <= 0:
            return 0.0
        return float(max(np.max(np.abs(vp)) / np.sqrt(g11 * nt), 0.0))
    q = (g22 * vp * vp - 2 * g12 * vp * vc + g11 * vc * vc) / det
    return float(np.sqrt(max(np.max(q), 0.0) / nt))


_S: dict = {}


def _init(t, lensf):
    _S["t"], _S["lf"] = t, lensf


def objective(x, t, lensf=None):
    x = np.asarray(x, float)
    for v, (lo, hi) in zip(x, BOUNDS, strict=True):
        if not lo <= v <= hi:
            return 1.0
    try:
        hp, hc = hpc(*to_cart(x))
    except Exception:  # noqa: BLE001
        return 1.0
    if lensf is not None:
        hp, hc = hp * lensf, hc * lensf
    return 1.0 - stat(t, hp, hc)


def _obj_worker(x):
    return objective(x, _S["t"], _S["lf"])


def run_de(t, lensf, seed, cfg, x0=None):
    t1 = time.time()
    with Pool(4, initializer=_init, initargs=(t, lensf)) as pool:
        r = differential_evolution(
            _obj_worker,
            BOUNDS,
            popsize=cfg["popsize"],
            maxiter=cfg["maxiter"],
            mutation=(0.5, 1.0),
            recombination=0.7,
            init="sobol",
            seed=seed,
            tol=1e-4,
            workers=pool.map,
            updating="deferred",
            polish=False,
            x0=x0,
        )
    return r, {
        "de_fun": float(r.fun),
        "de_message": str(r.message),
        "de_success": bool(r.success),
        "de_nit": int(r.nit),
        "de_nfev": int(r.nfev),
        "de_spread": float(np.std(r.population_energies)),
        "seconds": round(time.time() - t1, 1),
    }


def polish(t, lensf, x):
    r = m.c.nm_restart(
        lambda z: objective(z, t, lensf), x, STEPS, maxfev=400, max_restarts=4, restart_fev=150
    )
    return r


def bounds_hit(x, lensp=None):
    hits = []
    if x[1] >= 0.99 * 15:
        hits.append("q_max")
    if x[1] <= 1.01:
        hits.append("q_min")
    if max(x[2], x[5]) >= 0.98:
        hits.append("chi")
    if not (0.02 < x[8] < np.pi - 0.02):
        hits.append("inc")
    if abs(np.exp(x[0]) / MC0 - 1) >= 0.99 * 0.6:
        hits.append("mc")
    if lensp is not None:
        if lensp[0] >= 0.99e4 or lensp[0] <= 1.01:
            hits.append("M_Lz")
        if lensp[1] >= 0.99 * 5 or lensp[1] <= 0.0101:
            hits.append("y")
    return hits


# ------------------------------------------------------------------ the one fitter
def fit_u(t, seed, x0=None):
    r, info = run_de(t, None, seed, DE_U, x0)
    p = polish(t, None, r.x)
    x = np.asarray(p["x"]) if p["f"] < r.fun else r.x
    f = min(p["f"], float(r.fun))
    info.update(
        {
            "polish_gain": float(r.fun) - f,
            "polish_converged": p["converged_restart"],
            "polish_nfev": p["nfev"],
        }
    )
    return {"f": f, "x": [float(v) for v in x], "lens": None, "log": [info]}


def fit_p(t, seed, u, x0_saved=None):
    g = m.grid()
    x, f, lensp, node_prev, log = np.asarray(u["x"]), 1.0, None, None, []
    if x0_saved is not None:  # (x, lens) of a saved solution: score it first
        fs = objective(x0_saved[0], t, b.lens_exact(*x0_saved[1], +1))
        log.append({"step": "saved", "f": fs})
        if fs < f:
            x, f, lensp = np.asarray(x0_saved[0]), fs, tuple(x0_saved[1])
    for cyc in range(3):
        hp, hc = hpc(*to_cart(x))
        scan = np.array([stat(t, hp * g["F"][k], hc * g["F"][k]) for k in range(len(g["nodes"]))])
        k = int(np.argmax(scan))
        if 1.0 - scan[k] < f:
            f, lensp = 1.0 - float(scan[k]), g["nodes"][k]
        if k == node_prev:
            log.append({"step": f"cycle{cyc}", "node_repeat": True, "f": f})
            break
        r, info = run_de(t, g["F"][k], seed + 100 * (cyc + 1), DE_P, x0=x)
        if r.fun < f:
            x, f, lensp = r.x, float(r.fun), g["nodes"][k]
        info.update({"step": f"cycle{cyc}", "node": list(map(float, g["nodes"][k])), "f": f})
        log.append(info)
        node_prev = k
    lml, y = float(np.log(lensp[0])), float(lensp[1])
    for rnd in range(3):
        gain = 0.0

        def lobj(v, x=x):
            if not (0.0 <= v[0] <= np.log(1e4) and 0.01 <= v[1] <= 5.0):
                return 1.0
            return objective(x, t, b.lens_exact(float(np.exp(v[0])), float(v[1]), +1))

        rl = minimize(
            lobj,
            [lml, y],
            method="Nelder-Mead",
            options={
                "maxfev": 40,
                "xatol": 1e-4,
                "fatol": 1e-7,
                "initial_simplex": [[lml, y], [lml + 0.1, y], [lml, y * 1.15]],
            },
        )
        if rl.fun < f:
            gain += f - float(rl.fun)
            f, lml, y = float(rl.fun), float(rl.x[0]), float(rl.x[1])
        p = polish(t, b.lens_exact(float(np.exp(lml)), y, +1), x)
        if p["f"] < f:
            gain += f - p["f"]
            f, x = p["f"], np.asarray(p["x"])
        log.append(
            {
                "step": f"polish{rnd}",
                "f": f,
                "gain": gain,
                "intrinsic_converged": p["converged_restart"],
                "lens_nfev": int(rl.nfev),
            }
        )
        if gain < TOL:
            break
    return {"f": f, "x": [float(v) for v in x], "lens": [float(np.exp(lml)), y], "log": log}


def finish(t, r, model):
    x = np.asarray(r["x"])
    hp, hc = hpc(*to_cart(x))
    if r["lens"] is not None:
        lf = b.lens_exact(*r["lens"], +1)
        hp, hc = hp * lf, hc * lf
    r["FF"] = 1.0 - r["f"]
    r["FF_pad64"] = stat(t, hp, hc, pad=64)
    r["bounds"] = bounds_hit(x, r["lens"])
    r["model"] = model
    mc, q, s1, s2, inc, phi = to_cart(x)
    r["physical"] = {"mc": mc, "q": q, "chi1": list(s1), "chi2": list(s2), "inc": inc, "phi": phi}
    return r


# ------------------------------------------------------------------ validation
V7 = (np.log(MC0), 2.0, (0.5, 0.2, 0.1), (-0.3, 0.1, 0.4), 2.0, 1.3, 1.0)


def injections():
    out = {}
    for name in ("V5", "V6"):
        _, _, x, lensp = m.valid_truth(name)
        mc, q, s1, s2, inc, phi, psi = d.unpack(x, False)
        out[name] = ((mc, q, s1, s2, inc, phi, psi), lensp)
    out["V7"] = ((float(np.exp(V7[0])), *V7[1:]), None)
    best = json.loads((OUT / "s1c_summary.json").read_text())["D"]["P"]
    mc, q, s1, s2, inc, phi, psi = d.unpack(best["x"], False)
    out["V9"] = ((mc, q, s1, s2, inc, phi, psi), tuple(best["lens"]))
    return out


def injection_target(spec):
    (mc, q, s1, s2, inc, phi, psi), lensp = spec
    hp, hc = hpc(mc, q, tuple(s1), tuple(s2), inc, phi)
    t = np.cos(2 * psi) * hp + np.sin(2 * psi) * hc
    lf = None if lensp is None else b.lens_exact(*lensp, +1)
    if lf is not None:
        t, hp, hc = t * lf, hp * lf, hc * lf
    return t, 1.0 - stat(t, hp, hc)


def save(name, rows, label):
    (OUT / name).write_text(
        json.dumps({"label": label, "rows": rows}, indent=1, default=float) + "\n"
    )


def main_validate():
    rows = []
    for name, spec in injections().items():
        t, truth_mm = injection_target(spec)
        for seed in (11, 12):
            t0 = time.time()
            u = finish(t, fit_u(t, seed), "U")
            res = u
            if spec[1] is not None:
                res = finish(t, fit_p(t, seed, u), "P")
                res["u_stage_FF"] = u["FF"]
            res.update({"name": name, "seed": seed, "truth_mismatch": truth_mm})
            res["wall_s"] = round(time.time() - t0, 1)
            rows.append(res)
            print(
                name,
                seed,
                res["model"],
                round(res["FF"], 6),
                res["bounds"],
                res["wall_s"],
                flush=True,
            )
            save("s1d_validate.json", rows, "simulated; S1d tool validation (SciPy DE)")


def saved_d():
    """Every saved R1 family-D solution (S1c), re-scored with the S1d score on the same target."""
    _, t = m.target("D")
    out = []
    for r in json.loads((OUT / "s1c_case_DA.json").read_text())["rows"]:
        if r["family"] != "D":
            continue
        mc, q, s1, s2, inc, phi, _psi = d.unpack(r["x"], False)
        x = from_cart(mc, q, s1, s2, inc, phi)
        lf = None if r["lens"] is None else b.lens_exact(*r["lens"], +1)
        f = objective(x, t, lf)
        out.append(
            {
                "model": r["model"],
                "start": r["start_label"],
                "x": [float(v) for v in x],
                "lens": r["lens"],
                "FF_s1c_pad64": r["FF_pad64"],
                "FF_s1d": 1.0 - f,
            }
        )
    return out


def _attr_job(args):
    model, start = args
    _, t = m.target("D")
    r = m.fit(t, "A", model, start)
    r["FF_pad64"] = m.final_score(t, r["x"], "A", r["lens"])
    r["model"] = model
    return r


def main_case():
    _, t = m.target("D")
    saved = saved_d()
    rows = {"saved_rescored": saved, "fits": [], "attribution_xas": []}
    label = "simulated; S1d R1 case fits (SciPy DE, XPHM) and attribution"

    def dump():
        (OUT / "s1d_case.json").write_text(
            json.dumps({"label": label, **rows}, indent=1, default=float) + "\n"
        )

    dump()
    bu = max((z for z in saved if z["model"] == "U"), key=lambda z: z["FF_s1d"])
    bp = max((z for z in saved if z["model"] == "P"), key=lambda z: z["FF_s1d"])
    for seed in (21, 22, 23):
        t0 = time.time()
        u = finish(t, fit_u(t, seed, x0=np.array(bu["x"])), "U")
        u.update({"seed": seed, "wall_s": round(time.time() - t0, 1)})
        rows["fits"].append(u)
        print("U", seed, round(u["FF_pad64"], 6), u["bounds"], u["wall_s"], flush=True)
        dump()
        t0 = time.time()
        p = finish(t, fit_p(t, seed, u, x0_saved=(np.array(bp["x"]), bp["lens"])), "P")
        p.update({"seed": seed, "wall_s": round(time.time() - t0, 1)})
        rows["fits"].append(p)
        print("P", seed, round(p["FF_pad64"], 6), p["bounds"], p["lens"], p["wall_s"], flush=True)
        dump()
    jobs = [(mdl, s0) for mdl in ("U", "P") for s0 in m.generic_starts("A")]
    with Pool(4) as pool:
        for r in pool.imap_unordered(_attr_job, jobs):
            rows["attribution_xas"].append(r)
            print("XAS-on-XPHM", r["model"], round(r["FF_pad64"], 6), r["bounds"], flush=True)
            dump()


def summary():
    v = json.loads((OUT / "s1d_validate.json").read_text())["rows"]
    case_file = OUT / "s1d_case.json"  # absent when the tool fails validation (the round stops)
    c = (
        json.loads(case_file.read_text())
        if case_file.exists()
        else {"fits": [], "saved_rescored": []}
    )
    h, t = m.target("D")
    w = b.weight()
    ratio = float(np.sqrt(np.sum(np.abs(t) ** 2 * w) / np.sum(np.abs(h) ** 2 * w)))
    out = {"snr_ratio_lensed_over_unlensed": ratio, "validation": {}, "case": {}}
    for name in sorted({r["name"] for r in v}):
        rr = [r for r in v if r["name"] == name]
        out["validation"][name] = {
            "model": rr[0]["model"],
            "eps_per_seed": {str(r["seed"]): 1.0 - r["FF"] for r in rr},
            "truth_mismatch": rr[0]["truth_mismatch"],
        }
    out["tool_pass"] = all(
        e <= 5e-4 for z in out["validation"].values() for e in z["eps_per_seed"].values()
    )
    if not c["fits"]:
        out["case"] = "not run: the tool failed its frozen validation rule"
        (OUT / "s1d_summary.json").write_text(json.dumps(out, indent=1) + "\n")
        return out
    for model in ("U", "P"):
        fits = [r for r in c["fits"] if r["model"] == model]
        saved = [z["FF_s1d"] for z in c["saved_rescored"] if z["model"] == model]
        ffs = sorted((r["FF_pad64"] for r in fits), reverse=True)
        best = max(fits, key=lambda r: r["FF_pad64"]) if fits else None
        ff = max(ffs + saved) if (ffs or saved) else None
        rho = b.rho_req(1.0 - ff) if ff is not None else None
        out["case"][model] = {
            "FF_seeds_pad64": ffs,
            "best_saved_rescored": max(saved) if saved else None,
            "FF_best": ff,
            "rho_req": rho,
            "intrinsic_snr_needed": None if rho is None else rho / ratio,
            "top2_agree_2e-4": len(ffs) >= 2 and ffs[0] - ffs[1] < 2e-4,
            "best_physical": None if best is None else best["physical"],
            "best_lens": None if best is None else best["lens"],
            "best_bounds": None if best is None else best["bounds"],
        }
        a = [r for r in c.get("attribution_xas", []) if r["model"] == model]
        if a:
            out["case"][model]["xas_templates_on_xphm_target"] = max(r["FF_pad64"] for r in a)
    (OUT / "s1d_summary.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "validate":
        main_validate()
    elif argv and argv[0] == "case":
        main_case()
    elif argv and argv[0] == "summary":
        print(json.dumps(summary(), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
