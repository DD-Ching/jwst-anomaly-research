"""X4 S1: can ordinary binaries absorb the GW signature of a negative-mass point lens?

Unit file: docs/state/units/X4.md (S1; configuration frozen there before computing). Everything
here is `simulated` or `derived`; no strain data is used.

For each of 36 injections (negative point lens x IMRPhenomXAS binary) the best fitting factor (FF)
is found for two ordinary models, with amplitude, arrival time and phase maximised analytically:
  U  unlensed binary, free (chirp mass, mass ratio, chi1z, chi2z);
  P  positive point-mass lens x binary, the same intrinsics plus (M_Lz, y): a 10 x 10 lens grid with
     intrinsic Nelder-Mead at every node, then a 5 x 5 grid around the best node (main run); the
     frozen 6-parameter Nelder-Mead from the best three nodes runs in the ``refine`` stage.
Lens factors use the S0 closed form (mpmath) on every frequency, with a geometric-optics fallback
only where mpmath's 1F1 does not converge (very large w).

Screening threshold for the next stage (NOT a detection threshold or a calibrated false-alarm rate):
rho_req = 3 / sqrt(1 - FF^2) for the observed (lensed) SNR; pass if rho_req <= 30 against both U
and P and rho_req <= 100 x (lensed / unlensed SNR ratio).

Needs lalsuite (lalsimulation) in the running interpreter; it is not a project dependency.

  python scripts/x4_s1_degeneracy.py validate     # tool checks (FF of matching models)
  python scripts/x4_s1_degeneracy.py run          # 36 cases -> results/x4/s1_cases.json
  python scripts/x4_s1_degeneracy.py refine       # passing cases: frozen 6-parameter P refinement
  python scripts/x4_s1_degeneracy.py refine_wide  # the same with mass ratio up to 15 (U and P)
  python scripts/x4_s1_degeneracy.py ubound       # U alone, q <= 15 and chirp mass +-60 %
  python scripts/x4_s1_degeneracy.py summary      # best ordinary fit per case over all stages
  python scripts/x4_s1_degeneracy.py robust       # re-fit passing cases: other starts, wider bounds
  python scripts/x4_s1_degeneracy.py crosscheck   # 3 key cases fitted with IMRPhenomD templates
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "x4"
_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


gl = _load("s0_gw_neglens")

# ---- frozen configuration (docs/state/units/X4.md, S1) ----
F_LOW, F_HIGH, DF = 20.0, 1024.0, 0.25
SOURCES = (30.0, 60.0, 200.0)  # detector-frame total mass (Msun)
Q_INJ, CHI_INJ = 1.5, 0.2
LENS_M = (100.0, 300.0, 1000.0, 3000.0)
LENS_Y = (0.3, 1.0, 1.5)
MC_FRAC = 0.30  # chirp-mass bound +-30 % (robustness: 60 %)
Q_MAX, CHI_MAX = 8.0, 0.99
BOUNDS = {"q_max": Q_MAX}  # widened to 15 by the ubound check only
P_GRID_M = np.logspace(0.0, 4.0, 10)
P_GRID_Y = np.linspace(0.01, 5.0, 10)
APPROX = "IMRPhenomXAS"  # injections, and fits unless overridden
FIT_APPROX = {"name": APPROX}
# overlap time-lag resolution: pad 4 in the main run (0.24 ms lags), pad 16 in the refine stage
OVERLAP_PAD = {"pad": 4}
CROSSCHECK_CASES = ((30.0, 300.0, 1.0), (60.0, 1000.0, 1.5), (200.0, 300.0, 1.0))

_FREQS = np.arange(0.0, F_HIGH + DF / 2, DF)
_BAND = _FREQS >= F_LOW
_CACHE: dict = {}


def psd() -> np.ndarray:
    if "psd" not in _CACHE:
        import lal
        import lalsimulation as ls

        s = lal.CreateREAL8FrequencySeries("psd", 0, 0.0, DF, lal.HertzUnit, _FREQS.size)
        ls.SimNoisePSDaLIGODesignSensitivityP1200087(s, F_LOW)
        _CACHE["psd"] = np.array(s.data.data)
    return _CACHE["psd"]


def weight() -> np.ndarray:
    if "w" not in _CACHE:
        p = psd()
        ok = _BAND & (p > 0)  # the design curve is 0 at and below f_low
        w = np.zeros_like(_FREQS)
        w[ok] = 1.0 / p[ok]
        _CACHE["w"] = w
    return _CACHE["w"]


def masses(mc: float, q: float) -> tuple[float, float]:
    m = mc * (1.0 + q) ** 1.2 / q**0.6
    return m * q / (1.0 + q), m / (1.0 + q)


def chirp(mtot: float, q: float) -> float:
    m1, m2 = mtot * q / (1 + q), mtot / (1 + q)
    return (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2


def waveform(mc: float, q: float, c1: float, c2: float, approx: str | None = None) -> np.ndarray:
    import lal
    import lalsimulation as ls

    m1, m2 = masses(mc, q)
    hp, _ = ls.SimInspiralChooseFDWaveform(
        m1 * lal.MSUN_SI,
        m2 * lal.MSUN_SI,
        0,
        0,
        c1,
        0,
        0,
        c2,
        1e8 * lal.PC_SI,
        0,
        0,
        0,
        0,
        0,
        DF,
        F_LOW,
        F_HIGH,
        F_LOW,
        lal.CreateDict(),
        ls.GetApproximantFromString(approx or FIT_APPROX["name"]),
    )
    h = np.zeros(_FREQS.size, complex)
    n = min(hp.data.length, _FREQS.size)
    h[:n] = hp.data.data[:n]
    h[~_BAND] = 0.0
    return h


def amp_geometric(w: float, y: float, s: int) -> complex:
    """Geometric-optics F with the closed form's phase convention (time origin = unlensed arrival).

    Used only where mpmath's 1F1 fails to converge (w of several hundred and up); it agrees with the
    exact F to <= 1e-3 at w >= 150 (tests/test_x4_s1.py)."""
    if s > 0:
        r = np.sqrt(y * y + 4.0)
        xs = ((y + r) / 2.0, (y - r) / 2.0)
    else:
        if y < 2.0:
            return 0.0 + 0.0j
        r = np.sqrt(y * y - 4.0)
        xs = ((y + r) / 2.0, (y - r) / 2.0)
    out = 0.0 + 0.0j
    for x in xs:
        mu = 1.0 / (1.0 - x**-4) if s > 0 else 1.0 / ((1.0 + x**-2) * (1.0 - x**-2))
        t = (x - y) ** 2 / 2.0 - s * np.log(abs(x))
        out += np.sqrt(abs(mu)) * np.exp(1j * w * t - 1j * np.pi * (0.0 if mu > 0 else 0.5))
    return out


def amp(w: float, y: float, s: int) -> complex:
    # Far in the geometric regime (w > 300 and w y^2 > 6000) mpmath's 1F1 is slow or fails, and
    # geometric optics agrees with the exact form to <= 1e-4 there (checked 2026-10-15); use it.
    # Only the refine stage reaches this region; main-run arrays were computed before this rule.
    if w > 300.0 and w * y * y > 6000.0:
        _CACHE["n_geometric"] = _CACHE.get("n_geometric", 0) + 1
        return amp_geometric(w, y, s)
    try:
        return gl.amp_factor(w, y, s)
    except Exception:  # noqa: BLE001  (mpmath NoConvergence at very large w)
        _CACHE["n_geometric"] = _CACHE.get("n_geometric", 0) + 1
        return amp_geometric(w, y, s)


def lens_cache_dir() -> Path:
    """On-disk cache of exact lens factors, shared by all stages and workers (data root, untracked).

    The directory name encodes the frozen frequency grid; arrays are only reused on that grid.
    """
    sys.path.insert(0, str(ROOT / "src"))
    from jwst_anomaly import paths

    d = paths.data_root() / "derived" / "x4" / f"lens_f{F_LOW:g}-{F_HIGH:g}_df{DF:g}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def lens(m_lz: float, y: float, s: int) -> np.ndarray:
    key = (round(m_lz, 6), round(y, 6), s)
    if key in _CACHE:
        return _CACHE[key]
    path = lens_cache_dir() / f"s{s:+d}_m{key[0]:.6f}_y{key[1]:.6f}.npy"
    if path.exists():
        _CACHE[key] = np.load(path)
        return _CACHE[key]
    f = _FREQS
    out = np.ones(f.size, complex)
    idx = np.where(_BAND)[0]
    out[idx] = [amp(gl.W_PER_MSUN_HZ * m_lz * f[i], y, s) for i in idx]
    tmp = path.with_suffix(f".{os.getpid()}.tmp.npy")
    np.save(tmp, out)
    os.replace(tmp, path)  # atomic: concurrent workers never read a partial file
    _CACHE[key] = out
    return out


def overlap(a: np.ndarray, b: np.ndarray, pad: int | None = None) -> float:
    """max over time and phase of the normalised inner product (amplitude free)."""
    w = weight()
    prod = a * np.conj(b) * w
    n = (pad or OVERLAP_PAD["pad"]) * prod.size
    corr = np.fft.ifft(prod, n) * n
    na = np.sum(np.abs(a) ** 2 * w)
    nb = np.sum(np.abs(b) ** 2 * w)
    if na <= 0 or nb <= 0:
        return 0.0
    return float(np.max(np.abs(corr)) / np.sqrt(na * nb))


# ----------------------------------------------------------------------------- fits
def _intrinsic_obj(x, target, lensf, mc0, frac):
    lmc, q, c1, c2 = x
    mc = np.exp(lmc)
    if not (
        abs(mc / mc0 - 1) <= frac
        and 1.0 <= q <= BOUNDS["q_max"]
        and abs(c1) <= CHI_MAX
        and abs(c2) <= CHI_MAX
    ):
        return 1.0
    try:
        h = waveform(mc, q, c1, c2)
    except Exception:  # noqa: BLE001  (waveform outside its domain)
        return 1.0
    if lensf is not None:
        h = h * lensf
    return 1.0 - overlap(target, h)


def fit_intrinsic(target, lensf, mc0, starts, frac=MC_FRAC, maxfev=300):
    best = (1.0, None)
    for s in starts:
        r = minimize(
            _intrinsic_obj,
            s,
            args=(target, lensf, mc0, frac),
            method="Nelder-Mead",
            options={
                "maxfev": maxfev,
                "xatol": 1e-4,
                "fatol": 1e-7,
                "initial_simplex": _simplex(s),
            },
        )
        if r.fun < best[0]:
            best = (float(r.fun), [float(v) for v in r.x])
    return best


def _simplex(s):
    s = np.asarray(s, float)
    steps = np.array([0.03, 0.3, 0.2, 0.2])
    return np.vstack([s] + [s + np.eye(4)[i] * steps[i] for i in range(4)])


def coarse_start(target, mc0, frac):
    best = (1.0, [np.log(mc0), 1.0, 0.0, 0.0])
    for mc in mc0 * np.linspace(1 - frac, 1 + frac, 13):
        for q in (1.0, 1.5, 2.5, 4.0, 7.0):
            mm = _intrinsic_obj([np.log(mc), q, 0.0, 0.0], target, None, mc0, frac)
            if mm < best[0]:
                best = (mm, [np.log(mc), q, 0.0, 0.0])
    return best[1]


def fit_unlensed(target, mc0, rng, frac=MC_FRAC, n_random=2):
    starts = [[np.log(mc0), Q_INJ, CHI_INJ, CHI_INJ], coarse_start(target, mc0, frac)]
    for _ in range(n_random):
        starts.append(
            [
                np.log(mc0 * rng.uniform(1 - frac, 1 + frac)),
                rng.uniform(1, 5),
                *rng.uniform(-0.8, 0.8, 2),
            ]
        )
    mm, x = fit_intrinsic(target, None, mc0, starts, frac)
    return {"mismatch": mm, "params": x}


def fit_positive(target, mc0, start, frac=MC_FRAC, grid_m=P_GRID_M, grid_y=P_GRID_Y):
    nodes = []
    for m in grid_m:
        for y in grid_y:
            mm, x = fit_intrinsic(target, lens(m, y, +1), mc0, [start], frac, maxfev=150)
            nodes.append((mm, m, y, x))
    nodes.sort(key=lambda t: t[0])
    best = nodes[0]
    # 5 x 5 refinement around the best node (half the grid spacing)
    lm, y0 = np.log10(best[1]), best[2]
    dm = np.diff(np.log10(grid_m))[0] / 2
    dy = np.diff(grid_y)[0] / 2
    for m in 10 ** np.linspace(lm - dm, lm + dm, 5):
        for y in np.linspace(max(0.01, y0 - dy), y0 + dy, 5):
            mm, x = fit_intrinsic(target, lens(m, y, +1), mc0, [best[3], start], frac, maxfev=200)
            if mm < best[0]:
                best = (mm, float(m), float(y), x)
    return {
        "mismatch": best[0],
        "M_Lz": best[1],
        "y": best[2],
        "params": best[3],
        "grid_best": nodes[:3],
    }


def lens_exact(m_lz: float, y: float, s: int) -> np.ndarray:
    """Exact lens factor at an arbitrary (M_Lz, y), not cached on disk (continuous refinement)."""
    out = np.ones(_FREQS.size, complex)
    idx = np.where(_BAND)[0]
    out[idx] = [amp(gl.W_PER_MSUN_HZ * m_lz * _FREQS[i], y, s) for i in idx]
    return out


def _p6_obj(x, target, mc0, frac, memo):
    lmc, q, c1, c2, lml, y = x
    if not (0.0 <= lml <= np.log(1e4) and 0.01 <= y <= 5.0):
        return 1.0
    key = (round(lml, 7), round(y, 7))
    if key not in memo:
        memo.clear()  # one lens array at a time keeps memory flat
        memo[key] = lens_exact(np.exp(lml), y, +1)
    return _intrinsic_obj([lmc, q, c1, c2], target, memo[key], mc0, frac)


def refine_positive(target, mc0, nodes, frac=MC_FRAC, maxfev=80):
    """Frozen S1 refinement: Nelder-Mead over all six parameters from the best three grid nodes."""
    best = (1.0, None)
    for _mm, m, y, x in nodes[:3]:
        s0 = np.array([*x, np.log(m), y], float)
        steps = np.array([0.03, 0.3, 0.2, 0.2, 0.2, 0.1])
        simplex = np.vstack([s0] + [s0 + np.eye(6)[i] * steps[i] for i in range(6)])
        r = minimize(
            _p6_obj,
            s0,
            args=(target, mc0, frac, {}),
            method="Nelder-Mead",
            options={"maxfev": maxfev, "xatol": 1e-4, "fatol": 1e-7, "initial_simplex": simplex},
        )
        if r.fun < best[0]:
            best = (float(r.fun), [float(v) for v in r.x])
    return {
        "mismatch": best[0],
        "M_Lz": float(np.exp(best[1][4])),
        "y": best[1][5],
        "params": best[1][:4],
    }


def refine_case(args):
    """Refine stage for one case: pad-16 overlaps, U polish, frozen 6-parameter P refinement."""
    mtot, m_lz, y, prev_u, prev_p = args
    t1 = time.time()
    OVERLAP_PAD["pad"] = 16
    mc0 = chirp(mtot, Q_INJ)
    h = waveform(mc0, Q_INJ, CHI_INJ, CHI_INJ, approx=APPROX)
    target = h * lens(m_lz, y, -1)
    w = weight()
    ratio = float(np.sqrt(np.sum(np.abs(target) ** 2 * w) / np.sum(np.abs(h) ** 2 * w)))
    mm_u, xu = fit_intrinsic(target, None, mc0, [prev_u], MC_FRAC)
    p_grid = fit_positive(target, mc0, xu)  # grid nodes recomputed at pad 16
    p6 = refine_positive(target, mc0, p_grid["grid_best"])
    p_best = (
        p6
        if p6["mismatch"] < p_grid["mismatch"]
        else {k: v for k, v in p_grid.items() if k != "grid_best"}
    )
    ru, rp = rho_req(mm_u), rho_req(p_best["mismatch"])
    worst = max(ru, rp)
    return {
        "M_total": mtot,
        "M_Lz": m_lz,
        "y": y,
        "snr_ratio_lensed_over_unlensed": ratio,
        "U": {"mismatch": mm_u, "params": xu},
        "P_grid_pad16": {k: v for k, v in p_grid.items() if k != "grid_best"},
        "P_6d": p6,
        "P": p_best,
        "rho_req_U": ru,
        "rho_req_P": rp,
        "rho_req": worst,
        "intrinsic_snr_needed": worst / ratio,
        "passes_screen": bool(worst <= 30.0 and worst <= 100.0 * ratio),
        "seconds": round(time.time() - t1, 1),
        "stage1": {"U": prev_u, "P": prev_p},
    }


def _ubound_job(c):
    """U re-fit with widened bounds (set in the parent before the fork)."""
    mc0 = chirp(c["M_total"], Q_INJ)
    h = waveform(mc0, Q_INJ, CHI_INJ, CHI_INJ, approx=APPROX)
    target = h * lens(c["M_Lz"], c["y"], -1)
    u = fit_unlensed(target, mc0, np.random.default_rng(1500), frac=0.60, n_random=4)
    return {
        "case": [c["M_total"], c["M_Lz"], c["y"]],
        "U_wide": u,
        "rho_req_U": rho_req(u["mismatch"]),
    }


def rho_req(mm: float) -> float:
    o = 1.0 - mm
    return float(3.0 / np.sqrt(max(1.0 - o * o, 1e-24)))


# ----------------------------------------------------------------------------- cases
def cases():
    return [(m, ml, y) for m in SOURCES for ml in LENS_M for y in LENS_Y]


def run_case(args):
    mtot, m_lz, y, seed, frac = args
    t1 = time.time()
    n_geo0 = _CACHE.get("n_geometric", 0)
    rng = np.random.default_rng(seed)
    mc0 = chirp(mtot, Q_INJ)
    h = waveform(mc0, Q_INJ, CHI_INJ, CHI_INJ, approx=APPROX)
    fneg = lens(m_lz, y, -1)
    target = h * fneg
    w = weight()
    ratio = float(np.sqrt(np.sum(np.abs(target) ** 2 * w) / np.sum(np.abs(h) ** 2 * w)))
    u = fit_unlensed(target, mc0, rng, frac)
    p = fit_positive(target, mc0, u["params"], frac)
    ru, rp = rho_req(u["mismatch"]), rho_req(p["mismatch"])
    worst = max(ru, rp)
    return {
        "M_total": mtot,
        "M_Lz": m_lz,
        "y": y,
        "snr_ratio_lensed_over_unlensed": ratio,
        "U": u,
        "P": {k: v for k, v in p.items() if k != "grid_best"},
        "rho_req_U": ru,
        "rho_req_P": rp,
        "rho_req": worst,
        "intrinsic_snr_needed": worst / ratio if ratio > 0 else None,
        "passes_screen": bool(worst <= 30.0 and worst <= 100.0 * ratio),
        "seconds": round(time.time() - t1, 1),
        "n_geometric_fallback_in_case": _CACHE.get("n_geometric", 0) - n_geo0,
        "chirp_bound_frac": frac,
    }


def _lens_job(args):
    m, y, s = args
    before = _CACHE.get("n_geometric", 0)
    arr = lens(m, y, s)
    return (round(m, 6), round(y, 6), s), arr, _CACHE.get("n_geometric", 0) - before


def precompute(m_list, y_list, s, procs=4):
    """Exact lens factors computed in parallel, then cached here for the forked case workers."""
    jobs = [(m, y, s) for m in m_list for y in y_list]
    with Pool(procs) as pool:
        for key, arr, n_geo in pool.imap_unordered(_lens_job, jobs):
            _CACHE[key] = arr
            _CACHE["n_geometric"] = _CACHE.get("n_geometric", 0) + n_geo


def main(argv=None) -> int:
    cmd = (argv or sys.argv[1:] or ["run"])[0]
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if cmd == "validate":
        rng = np.random.default_rng(3)
        res = {}
        for mtot in (60.0, 200.0):
            mc0 = chirp(mtot, Q_INJ)
            h = waveform(mc0, Q_INJ, CHI_INJ, CHI_INJ)
            # U on an unlensed injection, started away from the truth (no injection start)
            starts = [coarse_start(h, mc0, MC_FRAC)]
            mm, x = fit_intrinsic(h, None, mc0, starts)
            res[f"U_on_unlensed_M{mtot:g}"] = {"FF": 1 - mm, "params": x, "pass": 1 - mm >= 0.999}
            # P on an off-grid positive-lens injection
            target = h * lens(250.0, 0.7, +1)
            u = fit_unlensed(target, mc0, rng)
            p = fit_positive(target, mc0, u["params"])
            res[f"P_on_positive_M{mtot:g}"] = {
                "FF": 1 - p["mismatch"],
                "M_Lz": p["M_Lz"],
                "y": p["y"],
                "FF_U": 1 - u["mismatch"],
                "pass": 1 - p["mismatch"] >= 0.995,
            }
            print(json.dumps(res, default=float)[-400:], flush=True)
        res["n_geometric_fallback"] = _CACHE.get("n_geometric", 0)
        (OUT / "s1_validation.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    elif cmd in ("refine", "refine_wide"):
        if cmd == "refine_wide":  # U and P mass ratio up to 15 (U fits pile up at q = 8)
            BOUNDS["q_max"] = 15.0
        prev = json.loads((OUT / "s1_cases.json").read_text())["cases"]
        jobs = [
            (c["M_total"], c["M_Lz"], c["y"], c["U"]["params"], c["P"])
            for c in prev
            if c["passes_screen"]
        ]
        with Pool(4) as pool:
            rows = []
            for r in pool.imap_unordered(refine_case, jobs):
                rows.append(r)
                print(
                    f"M={r['M_total']:.0f} ML={r['M_Lz']:.0f} y={r['y']}"
                    f" FF_U={1 - r['U']['mismatch']:.5f}"
                    f" FF_P={1 - r['P']['mismatch']:.5f} rho_req={r['rho_req']:.1f}"
                    f" pass={r['passes_screen']} ({r['seconds']} s)",
                    flush=True,
                )
        rows.sort(key=lambda r: (r["M_total"], r["M_Lz"], r["y"]))
        rec = {
            "label": "simulated; refine stage (frozen 6-parameter P refinement, pad-16 overlaps)",
            "wall_s": round(time.time() - t0, 1),
            "cases": rows,
        }
        rec["q_max"] = BOUNDS["q_max"]
        name = "s1_refine.json" if cmd == "refine" else "s1_refine_wide.json"
        (OUT / name).write_text(json.dumps(rec, indent=1, default=float) + "\n")
    elif cmd == "ubound":
        # Many U fits end at q = 8: re-fit U on passing cases with q <= 15, chirp mass +-60 %.
        BOUNDS["q_max"] = 15.0
        OVERLAP_PAD["pad"] = 16
        prev = json.loads((OUT / "s1_refine.json").read_text())["cases"]

        with Pool(4) as pool:
            rows = pool.map(_ubound_job, prev)
        (OUT / "s1_ubound.json").write_text(
            json.dumps(
                {"label": "simulated; U re-fit with q <= 15, Mc +-60 %", "cases": rows},
                indent=1,
                default=float,
            )
            + "\n"
        )
    elif cmd == "summary":
        # Best ordinary fit per passing case over every stage on disk (no new fitting). Fits made
        # with q <= 8 stay valid under q <= 15, so the q <= 15 column takes the best of all.
        def ff(c, m):
            return 1.0 - c[m]["mismatch"]

        def load(name):
            return json.loads((OUT / name).read_text())["cases"]

        def key(c):
            return (c["M_total"], c["M_Lz"], c["y"])

        q8 = {}
        for name in ("s1_cases.json", "s1_robust.json", "s1_refine.json"):
            for c in load(name):
                b = q8.setdefault(
                    key(c), {"U": 0.0, "P": 0.0, "ratio": c["snr_ratio_lensed_over_unlensed"]}
                )
                b["U"], b["P"] = max(b["U"], ff(c, "U")), max(b["P"], ff(c, "P"))
        q15 = {k: dict(v) for k, v in q8.items()}
        for c in load("s1_refine_wide.json"):
            b = q15[key(c)]
            b["U"], b["P"] = max(b["U"], ff(c, "U")), max(b["P"], ff(c, "P"))
        for c in load("s1_ubound.json"):
            b = q15[tuple(c["case"])]
            b["U"] = max(b["U"], 1.0 - c["U_wide"]["mismatch"])
        passing = {key(c) for c in load("s1_refine.json")}
        out = {
            "label": "derived: best ordinary FF over all S1 stages (q <= 8 and q <= 15)",
            "cases": [],
        }
        for k in sorted(passing):
            row = {"case": list(k)}
            for tag, d in (("q8", q8[k]), ("q15", q15[k])):
                worst = max(rho_req(1 - d["U"]), rho_req(1 - d["P"]))
                row[tag] = {
                    "FF_U": d["U"],
                    "FF_P": d["P"],
                    "rho_req": worst,
                    "intrinsic_snr_needed": worst / d["ratio"],
                    "passes_screen": bool(worst <= 30.0 and worst <= 100.0 * d["ratio"]),
                }
            out["cases"].append(row)
        (OUT / "s1_summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
        for r in out["cases"]:
            print(
                r["case"],
                {t: round(r[t]["rho_req"], 1) for t in ("q8", "q15")},
                {t: r[t]["passes_screen"] for t in ("q8", "q15")},
            )
    elif cmd == "crosscheck":
        FIT_APPROX["name"] = "IMRPhenomD"
        jobs = [(m, ml, y, 700 + i, MC_FRAC) for i, (m, ml, y) in enumerate(CROSSCHECK_CASES)]
        with Pool(3) as pool:
            rows = sorted(pool.map(run_case, jobs), key=lambda r: (r["M_total"], r["M_Lz"], r["y"]))
        rec = {
            "label": "simulated; injection IMRPhenomXAS, fits IMRPhenomD (waveform cross-check)",
            "cases": rows,
        }
        (OUT / "s1_crosscheck.json").write_text(json.dumps(rec, indent=1, default=float) + "\n")
    elif cmd in ("run", "robust"):
        # exact lens factors shared by forked workers
        precompute(P_GRID_M, P_GRID_Y, +1)
        precompute(LENS_M, LENS_Y, -1)
        print(f"lens factors ready in {time.time() - t0:.0f} s", flush=True)
        if cmd == "run":
            jobs = [(m, ml, y, 100 + i, MC_FRAC) for i, (m, ml, y) in enumerate(cases())]
            name = "s1_cases.json"
        else:
            prev = json.loads((OUT / "s1_cases.json").read_text())["cases"]
            sel = [c for c in prev if c["passes_screen"]]
            jobs = [
                (c["M_total"], c["M_Lz"], c["y"], seed, 0.60)
                for c in sel
                for seed in (901, 902, 903)
            ]
            name = "s1_robust.json"
        with Pool(4) as pool:
            rows = []
            for r in pool.imap_unordered(run_case, jobs):
                rows.append(r)
                print(
                    f"M={r['M_total']:.0f} ML={r['M_Lz']:.0f} y={r['y']}"
                    f" ratio={r['snr_ratio_lensed_over_unlensed']:.3f}"
                    f" FF_U={1 - r['U']['mismatch']:.5f} FF_P={1 - r['P']['mismatch']:.5f}"
                    f" rho_req={r['rho_req']:.1f} pass={r['passes_screen']} ({r['seconds']} s)",
                    flush=True,
                )
        rows.sort(key=lambda r: (r["M_total"], r["M_Lz"], r["y"]))
        rec = {
            "label": "simulated injections; derived fitting factors; no data",
            "waveform": APPROX,
            "psd": "lalsimulation.SimNoisePSDaLIGODesignSensitivityP1200087 (LIGO-P1200087)",
            "band_hz": [F_LOW, F_HIGH, DF],
            "screen": (
                "rho_req = 3/sqrt(1-FF^2) <= 30 vs U and P, and <= 100 x SNR ratio (screening only)"
            ),
            "wall_s": round(time.time() - t0, 1),
            "n_geometric_fallback_parent": _CACHE.get("n_geometric", 0),
            "cases": rows,
        }
        (OUT / name).write_text(json.dumps(rec, indent=1, default=float) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
