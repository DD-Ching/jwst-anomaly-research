"""X4 S1c: R1 only -- one validated ordinary fitter, a frozen lens-grid diagnostic, and XPHM.

Unit file: docs/state/units/X4.md (S1c; frozen before computing). `simulated` / `derived`; no data.
Needs lalsuite. One fitter (`fit`) serves the tool validation and the case fit alike.

Families: A = IMRPhenomXAS (aligned spins, dominant mode; |correlation| maximised, phase free);
D = IMRPhenomXPHM (precessing, higher modes; real part maximised, no free overall phase).
Models: U unlensed; P x positive point lens (M_Lz 1-1e4 Msun, y 0.01-5).

P fit: cycles of {scan of the frozen lens grid with the intrinsics fixed -> best node; intrinsic
Nelder-Mead (restarted to convergence) with that lens fixed} until the grid-best node repeats and a
cycle gains < TOL, then a continuous lens polish (2-D Nelder-Mead, exact lens factors) alternated
with the intrinsic fit.

Speed (checked, see X4.md): FFT length next_fast_len(pad * N) (never shorter; FF change <= 1.6e-7);
a grid scan shares the waveform, target product and node norms and loops the FFT per node (a
batched FFT was measured slower here).

  python scripts/x4_s1c_r1.py validate A,D        # -> results/x4/s1c_validate_AD.json
  python scripts/x4_s1c_r1.py case D,A D         # -> results/x4/s1c_case_DA.json (as run)
  (the third argument lists the families that also start from the saved R1 solutions)
  python scripts/x4_s1c_r1.py summary            # -> results/x4/s1c_summary.json
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


d = _load("x4_s1b_xphm")
c = d.c
b = d.b

TOL = 1e-5
Q_MAX, CHI_MAX, MC_FRAC = 15.0, 0.99, 0.60
FIT_PAD, SCORE_PAD = 32, 64
GRID_M = np.logspace(0.0, 4.0, 41)
GRID_Y = np.logspace(np.log10(0.01), np.log10(5.0), 31)
R1 = (30.0, 300.0, 0.3)
MC0 = b.chirp(R1[0], b.Q_INJ)
NM_BUDGET = {"A": (400, 4, 150), "D": (1500, 6, 400)}  # first run, max restarts, restart fev
STEPS = {"A": np.array([0.03, 0.5, 0.2, 0.2]), "D": d.P_STEPS}
_G: dict = {}


# ------------------------------------------------------------------ waveforms and overlaps
def in_domain(x, fam):
    mc = np.exp(x[0])
    if not (abs(mc / MC0 - 1) <= MC_FRAC and 1.0 <= x[1] <= Q_MAX):
        return False
    if fam == "A":
        return abs(x[2]) <= CHI_MAX and abs(x[3]) <= CHI_MAX
    return (
        np.linalg.norm(x[2:5]) <= CHI_MAX
        and np.linalg.norm(x[5:8]) <= CHI_MAX
        and 0.0 <= x[8] <= np.pi
    )


def wave(x, fam):
    if fam == "A":
        return b.waveform(float(np.exp(x[0])), float(x[1]), float(x[2]), float(x[3]))
    mc, q, s1, s2, inc, phi, psi = d.unpack(x, False)
    return d.xphm(mc, q, s1, s2, inc, phi, psi)


def _nfft(pad):
    return sf.next_fast_len(pad * b._FREQS.size)


def ov(t, h, fam, pad=FIT_PAD):
    w = b.weight()
    prod = t * np.conj(h) * w
    n = _nfft(pad)
    corr = sf.ifft(prod, n) * n
    nt, nh = np.sum(np.abs(t) ** 2 * w), np.sum(np.abs(h) ** 2 * w)
    if nt <= 0 or nh <= 0:
        return 0.0
    peak = np.max(np.abs(corr)) if fam == "A" else np.max(corr.real)
    return float(peak / np.sqrt(nt * nh))


def grid():
    if "F" not in _G:
        nodes = [(m, y) for m in GRID_M for y in GRID_Y]
        _G["nodes"] = nodes
        _G["F"] = np.stack([b.lens(float(m), float(y), +1) for m, y in nodes])
        _G["F2"] = np.abs(_G["F"]) ** 2
    return _G


def grid_scan(t, h, fam, pad=FIT_PAD):
    """Overlap of t with h x F_k for every grid node k (intrinsics fixed). Returns the map."""
    g = grid()
    w = b.weight()
    base = t * np.conj(h) * w
    nt = np.sum(np.abs(t) ** 2 * w)
    nh = g["F2"] @ (np.abs(h) ** 2 * w)  # all node norms in one product
    n = _nfft(pad)
    out = np.zeros(len(g["nodes"]))
    for k in range(out.size):
        corr = sf.ifft(base * np.conj(g["F"][k]), n) * n
        peak = np.max(np.abs(corr)) if fam == "A" else np.max(corr.real)
        out[k] = peak / np.sqrt(nt * nh[k]) if nh[k] > 0 else 0.0
    return out


# ------------------------------------------------------------------ the one fitter
def nm_conv(fun, x0, fam):
    first, max_r, rfev = NM_BUDGET[fam]
    return c.nm_restart(fun, x0, STEPS[fam], maxfev=first, max_restarts=max_r, restart_fev=rfev)


def obj(x, t, fam, lensf=None):
    x = np.asarray(x, float)
    if not in_domain(x, fam):
        return 1.0
    try:
        h = wave(x, fam)
    except Exception:  # noqa: BLE001
        return 1.0
    if lensf is not None:
        h = h * lensf
    return 1.0 - ov(t, h, fam)


def lens_polish(x, t, fam, lml, y):
    def lobj(v):
        if not (0.0 <= v[0] <= np.log(1e4) and 0.01 <= v[1] <= 5.0):
            return 1.0
        return obj(x, t, fam, b.lens_exact(float(np.exp(v[0])), float(v[1]), +1))

    r = minimize(
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
    return float(r.fun), float(r.x[0]), float(r.x[1]), int(r.nfev), str(r.message)


def fit(t, fam, model, start, lens0=None):
    """The single fitter (validation and case). Returns f = 1 - FF and diagnostics."""
    t1 = time.time()
    x = np.asarray(start, float)
    log = []
    if model == "U":
        r = nm_conv(lambda z: obj(z, t, fam), x, fam)
        f, x = r["f"], np.asarray(r["x"])
        log.append({"step": "intrinsic", **{k: r[k] for k in ("nfev", "converged_restart")}})
        lensp, conv = None, r["converged_restart"]
    else:
        g = grid()
        node_prev, f, conv, lensp = None, 1.0, False, lens0
        if lens0 is not None:  # a saved solution's lens: fit the intrinsics there first
            r = nm_conv(lambda z: obj(z, t, fam, b.lens_exact(*lens0, +1)), x, fam)
            x, f = np.asarray(r["x"]), r["f"]
        for cyc in range(5):
            m = grid_scan(t, wave(x, fam), fam)
            k = int(np.argmax(m))
            lf = g["F"][k]
            r = nm_conv(lambda z, lf=lf: obj(z, t, fam, lf), x, fam)
            gain = f - r["f"]
            if r["f"] < f:
                f, x = r["f"], np.asarray(r["x"])
                lensp = g["nodes"][k]
            log.append(
                {
                    "step": f"cycle{cyc}",
                    "node": [float(v) for v in g["nodes"][k]],
                    "grid_best": float(m[k]),
                    "f": f,
                    "nfev": r["nfev"],
                    "converged_restart": r["converged_restart"],
                }
            )
            if k == node_prev and gain < TOL:
                conv = True
                break
            node_prev = k
        _G["last_map"] = m
        if lensp is None:  # no cycle improved on f = 1 (not reached in the recorded runs)
            lensp = g["nodes"][k]
        lml, y = float(np.log(lensp[0])), float(lensp[1])
        pconv = False
        for rnd in range(3):
            fl, lml2, y2, nfl, msg = lens_polish(x, t, fam, lml, y)
            gain = 0.0
            if fl < f:
                gain, f, lml, y = f - fl, fl, lml2, y2
            lf = b.lens_exact(float(np.exp(lml)), y, +1)
            r = nm_conv(lambda z, lf=lf: obj(z, t, fam, lf), x, fam)
            if r["f"] < f:
                gain += f - r["f"]
                f, x = r["f"], np.asarray(r["x"])
            log.append({"step": f"polish{rnd}", "lens": [float(np.exp(lml)), y], "f": f})
            if gain < TOL:
                pconv = True
                break
        lensp = (float(np.exp(lml)), y)
        conv = conv and pconv
    return {
        "f": float(f),
        "FF": 1.0 - float(f),
        "x": [float(v) for v in x],
        "lens": None if lensp is None else [float(v) for v in lensp],
        "converged": bool(conv),
        "bounds": bounds(x, fam, lensp),
        "log": log,
        "seconds": round(time.time() - t1, 1),
    }


def bounds(x, fam, lensp):
    hits = []
    if x[1] >= 0.99 * Q_MAX:
        hits.append("q_max")
    if abs(np.exp(x[0]) / MC0 - 1) >= 0.99 * MC_FRAC:
        hits.append("mc")
    spins = (x[2:4],) if fam == "A" else (x[2:5], x[5:8])
    if fam == "A" and max(abs(x[2]), abs(x[3])) >= 0.99 * CHI_MAX:
        hits.append("chi")
    if fam == "D" and max(np.linalg.norm(s) for s in spins) >= 0.99 * CHI_MAX:
        hits.append("chi")
    if lensp is not None:
        if lensp[0] >= 0.99e4 or lensp[0] <= 1.01:
            hits.append("M_Lz")
        if lensp[1] >= 0.99 * 5.0 or lensp[1] <= 0.0101:
            hits.append("y")
    return hits


def final_score(t, x, fam, lensp):
    h = wave(np.asarray(x), fam)
    if lensp is not None:
        h = h * b.lens_exact(*lensp, +1)
    return ov(t, h, fam, pad=SCORE_PAD)


# ------------------------------------------------------------------ targets, starts
def target(fam):
    if fam == "A":
        h = b.waveform(MC0, b.Q_INJ, b.CHI_INJ, b.CHI_INJ)
    else:
        h = d.xphm(MC0, b.Q_INJ, (0, 0, b.CHI_INJ), (0, 0, b.CHI_INJ), d.INC_TARGET, 0.0, 0.0)
    return h, h * b.lens(R1[1], R1[2], -1)


def generic_starts(fam):
    rng = np.random.default_rng(2026)
    out = []
    for _ in range(8 if fam == "A" else 6):
        s = [np.log(MC0) + rng.uniform(-0.2, 0.2), rng.uniform(1, 14)]
        if fam == "A":
            s += list(rng.uniform(-0.9, 0.9, 2))
        else:
            ip = rng.uniform(-0.4, 0.4, 4)
            cz = rng.uniform(-0.9, 0.9, 2)
            s += [ip[0], ip[1], cz[0], ip[2], ip[3], cz[1]]
            s += [rng.uniform(0.2, 2.9), rng.uniform(0, 2 * np.pi), rng.uniform(0, np.pi)]
        out.append(s)
    return out


def saved_starts():
    """Every saved R1 best solution (S1b steps A and B), family A, best (pad-64 FF) first."""
    out = []
    for r in json.loads((b.OUT / "s1b_rescore.json").read_text())["rows"]:
        if r["case"] == "R1":
            out.append((r["model"], r["params"], r["lens"], r["FF_pad64"]))
    for m in ("U", "P"):
        for r in json.loads((b.OUT / "s1b_fits.json").read_text())[m]:
            if r["case"] == "R1" and r.get("q_max", Q_MAX) == Q_MAX:
                xx = r["x"]
                ln = [float(np.exp(xx[4])), xx[5]] if m == "P" else None
                out.append((m, xx[:4], ln, r["FF_pad64"]))
    out.sort(key=lambda z: -z[3])
    seen, uniq = set(), []
    for m, p, ln, _ff in out:
        key = (m, tuple(round(v, 2) for v in p))
        if key not in seen:
            seen.add(key)
            uniq.append((m, p, ln))
    return uniq


def a_to_d(p):
    return [p[0], p[1], 0.0, 0.0, p[2], 0.0, 0.0, p[3], d.INC_TARGET, 0.0, 0.0]


# ------------------------------------------------------------------ validation
VALID = {
    "V1": ("A", "P", [2.547, 15.0, 0.632, 0.99], (2441.0, 0.096)),
    "V2": ("A", "P", [None, 4.0, 0.5, -0.3], (1000.0, 0.2)),
    "V3": ("A", "P", [None, 2.0, 0.2, 0.2], (300.0, 0.5)),
    "V4": ("A", "U", "s1b_U_best", None),
    "V5": ("D", "P", [None, 3.0, 0.3, 0.0, 0.2, 0.0, 0.3, 0.1, 1.0, 0.5, 0.3], (1000.0, 0.2)),
    "V6": ("D", "U", [None, 3.0, 0.3, 0.0, 0.2, 0.0, 0.3, 0.1, 1.0, 0.5, 0.3], None),
}


def valid_truth(name):
    fam, model, x, lensp = VALID[name]
    if x == "s1b_U_best":
        u = [r for r in json.loads((b.OUT / "s1b_fits.json").read_text())["U"] if r["case"] == "R1"]
        u = [r for r in u if r.get("q_max", Q_MAX) == Q_MAX]
        x = max(u, key=lambda r: r["FF_pad64"])["x"]
    x = [np.log(MC0) if v is None else v for v in x]
    return fam, model, np.asarray(x, float), lensp


def _valid_job(args):
    name, i = args
    fam, model, x, lensp = valid_truth(name)
    t = wave(x, fam)
    if lensp is not None:
        t = t * b.lens_exact(*lensp, +1)
    truth_check = obj(x, t, fam, None if lensp is None else b.lens_exact(*lensp, +1))
    r = fit(t, fam, model, generic_starts(fam)[i])
    r.update({"name": name, "start": i, "truth_mismatch": truth_check})
    return r


def _case_job(args):
    fam, model, start, lens0, label = args
    _, t = target(fam)
    r = fit(t, fam, model, start, lens0)
    r["FF_pad64"] = final_score(t, r["x"], fam, r["lens"])
    r.update({"family": fam, "model": model, "start_label": label})
    if model == "P":
        r["final_grid_map"] = [float(v) for v in _G.get("last_map", [])]
    return r


def run(jobs, fn, out_name, label):
    t0 = time.time()
    rows = []
    with Pool(4) as pool:
        for r in pool.imap_unordered(fn, jobs):
            rows.append(r)
            partial = {"label": label + " (partial)", "rows": rows}
            (b.OUT / out_name).write_text(json.dumps(partial, default=float) + "\n")
            print(
                r.get("name", r.get("family")),
                r.get("model", ""),
                r.get("start", r.get("start_label")),
                round(r["FF"], 5),
                "conv",
                r["converged"],
                r["bounds"],
                f"({r['seconds']} s)",
                flush=True,
            )
    best = {}
    for r in rows:  # keep the final grid map only for the best P fit per family (file size)
        if "final_grid_map" in r:
            k = r["family"]
            if k not in best or r["FF"] > best[k]["FF"]:
                best[k] = r
    for r in rows:
        if "final_grid_map" in r and best.get(r["family"]) is not r:
            del r["final_grid_map"]
        elif "final_grid_map" in r:
            r["final_grid_map"] = [round(v, 5) for v in r["final_grid_map"]]
    rec = {"label": label, "wall_s": round(time.time() - t0, 1), "rows": rows}
    (b.OUT / out_name).write_text(json.dumps(rec, indent=1, default=float) + "\n")


VALID_FOR = {
    ("A", "P"): ["V1", "V2", "V3"],
    ("A", "U"): ["V4"],
    ("D", "P"): ["V5"],
    ("D", "U"): ["V6"],
}


def summary() -> dict:
    cs = json.loads((b.OUT / "s1c_case_DA.json").read_text())["rows"]
    vs = json.loads((b.OUT / "s1c_validate_AD.json").read_text())["rows"]
    out = {}
    for fam in ("A", "D"):
        h, t = target(fam)
        w = b.weight()
        ratio = float(np.sqrt(np.sum(np.abs(t) ** 2 * w) / np.sum(np.abs(h) ** 2 * w)))
        out[fam] = {"snr_ratio_lensed_over_unlensed": ratio}
        for model in ("U", "P"):
            rr = [r for r in cs if r["family"] == fam and r["model"] == model]
            best = max(rr, key=lambda r: r["FF_pad64"])
            ff = best["FF_pad64"]
            rho = b.rho_req(1 - ff)
            eps = {
                n: 1 - max(r["FF"] for r in vs if r["name"] == n) for n in VALID_FOR[(fam, model)]
            }
            agree = sum(1 for r in rr if best["FF_pad64"] - r["FF_pad64"] < 2e-4)
            out[fam][model] = {
                "FF_pad64": ff,
                "rho_req": rho,
                "intrinsic_snr_needed": rho / ratio,
                "starts_within_2e-4": f"{agree}/{len(rr)}",
                "best_start": best["start_label"],
                "converged": best["converged"],
                "bounds": best["bounds"],
                "x": best["x"],
                "lens": best["lens"],
                "validation_eps": eps,
                "eps_limit": 0.1 * (1 - ff),
                "validation_pass": all(e <= 0.1 * (1 - ff) for e in eps.values()),
            }
    (b.OUT / "s1c_summary.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    what = argv[0] if argv else "validate"
    if what == "summary":
        print(json.dumps(summary(), indent=1))
        return 0
    fams = argv[1].split(",") if len(argv) > 1 else ["A", "D"]
    # families that also start from the saved R1 solutions (default: all; budget cut, see X4.md)
    with_saved = argv[2].split(",") if len(argv) > 2 else fams
    if what == "validate":
        jobs = [
            (n, i)
            for n, v in VALID.items()
            if v[0] in fams
            for i in range(len(generic_starts(v[0])))
        ]
        run(jobs, _valid_job, f"s1c_validate_{''.join(fams)}.json", "simulated; S1c validation")
    else:
        jobs = []
        for fam in fams:
            for model in ("U", "P"):
                for i, s in enumerate(generic_starts(fam)):
                    jobs.append((fam, model, s, None, f"generic{i}"))
            saved = saved_starts() if fam in with_saved else []
            if fam == "D":  # the two best saved solutions per model (cost of an XPHM fit)
                saved = [z for m in ("U", "P") for z in [z for z in saved if z[0] == m][:2]]
            for m, p, ln in saved:
                s = p if fam == "A" else a_to_d(p)
                jobs.append((fam, m, s, ln if m == "P" else None, f"saved_{m}"))
        run(jobs, _case_job, f"s1c_case_{''.join(fams)}.json", "simulated; S1c R1 case fits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
