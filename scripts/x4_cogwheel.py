"""X4 cogwheel round: cogwheel's published refinement as the U fitter (thin interface only).

Unit file: docs/state/units/X4.md (cogwheel round; plan frozen before checks and runs). `simulated`/
`derived`; no data. Needs cogwheel 1.6.0 (PyPI `cogwheel-pe`), lalsuite and mpmath in a separate
environment. Inputs: results/x4/handoff/targets.npz and the X4 score code; truths only for checks
and for scoring, never as a start.

  python scripts/x4_cogwheel.py checks      # -> results/x4/handoff/cogwheel_checks.json (F1-F3)
  python scripts/x4_cogwheel.py one C1 71   # one run (30-min cap applied by `runs`)
  python scripts/x4_cogwheel.py runs        # C1 and C2 -> results/x4/handoff/cogwheel_runs.json
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_DIR))
import x4_handoff as H  # noqa: E402

k, e, b = H._modules()
OUT_DIR = H.OUT_DIR
DET = "H"
TGPS = 1e9  # fixed and arbitrary: sets the detector orientation only
F_REF = 20.0
Q_MAX = 15.0
CAP_S = 30 * 60
SEEDS = {"C1": 71, "C2": 72}
EPS_PASS = 5e-4
REL_PASS = 0.1
FF_CASE = 0.915149
TOL_F1, TOL_F2, TOL_F3 = 1e-8, 1e-6, 1e-3
APPROX = "IMRPhenomXPHM"


def event(name):
    from cogwheel.data import EventData

    z = np.load(OUT_DIR / "targets.npz")
    return EventData(
        eventname=f"x4_{name}",
        frequencies=z["freqs_hz"],
        strain=z[f"{name}_target"][np.newaxis],
        wht_filter=np.sqrt(z["weight_1_over_psd"])[np.newaxis],
        detector_names=DET,
        tgps=TGPS,
        tcoarse=0.0,
    )


def to_par(x):
    """Our box point -> cogwheel intrinsic parameters (in-plane spins in its N frame)."""
    from cogwheel.waveform import inplane_spins_xy_to_xy_n

    mc, q, s1, s2, inc, phi = e.to_cart(np.asarray(x, float))
    m1, m2 = b.masses(mc, q)
    d = {"m1": m1, "m2": m2, "iota": inc, "phi_ref": phi, "f_ref": F_REF, "l1": 0.0, "l2": 0.0}
    for i, s in ((1, s1), (2, s2)):
        d[f"s{i}x"], d[f"s{i}y"], d[f"s{i}z"] = map(float, s)
    inplane_spins_xy_to_xy_n(d)
    for key in ("s1x", "s1y", "s2x", "s2y"):
        d.pop(key)
    return d


def physical(par):
    """cogwheel parameters -> (mc, q, s1, s2, inc, phi) in our frame, spin magnitudes unclipped.
    In-plane spins rotate by -phi_ref (cogwheel's `inplane_spins_xy_n_to_xy`, written out here
    so the PyCBC environment needs no cogwheel)."""
    m1, m2, phi = par["m1"], par["m2"], par["phi_ref"]
    assert m1 >= m2  # cogwheel samples q = m2 / m1 <= 1
    c, s_ = np.cos(phi), np.sin(phi)
    spins = []
    for i in (1, 2):
        xn, yn = par[f"s{i}x_n"], par[f"s{i}y_n"]
        spins.append((float(c * xn + s_ * yn), float(-s_ * xn + c * yn), float(par[f"s{i}z"])))
    mc = (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2
    return float(mc), float(m1 / m2), spins[0], spins[1], float(par["iota"]), float(phi)


def to_box(par):
    """Box point and whether it lies in our box (|chi| <= 0.99, q <= 15, Mc range)."""
    mc, q, s1, s2, inc, phi = physical(par)
    x = e.from_cart(mc, q, s1, s2, inc, phi)  # clips |chi| to 0.99
    inside = (
        np.linalg.norm(s1) <= 0.99
        and np.linalg.norm(s2) <= 0.99
        and q <= Q_MAX
        and e.BOUNDS[0][0] <= x[0] <= e.BOUNDS[0][1]
    )
    return x, bool(inside), [float(np.linalg.norm(s1)), float(np.linalg.norm(s2))]


def box_diff(x1, x2):
    dx = np.abs(np.asarray(x1) - np.asarray(x2))
    for i in (4, 7, 9):
        dx[i] = min(dx[i], 2 * np.pi - dx[i])
    return float(dx.max())


def check_points():
    pts = H.points()
    lo = np.array([bd[0] for bd in e.BOUNDS])
    hi = np.array([bd[1] for bd in e.BOUNDS])
    mid = 0.5 * (lo + hi)
    p_q = mid.copy()
    p_q[[1, 2, 5]] = 15.0, 0.99, 0.99
    p_t = mid.copy()
    p_t[[2, 3, 5, 6]] = 0.9, 0.999, 0.9, -0.999
    p_i = mid.copy()
    p_i[[2, 5, 8]] = 0.8, 0.8, 0.01
    p_j = p_i.copy()
    p_j[8] = np.pi - 0.01
    out = {}
    for n, p in pts.items():
        out[f"{n}_truth"] = p["x_truth"]
        out[f"{n}_found"] = p["x_found"]
    out.update({"q15_chi099": p_q, "tilts_pm0999": p_t, "inc_001": p_i, "inc_pi_m001": p_j})
    return out


def wf_generator(harmonic_modes=None):
    from cogwheel.waveform import WaveformGenerator

    return WaveformGenerator.from_event_data(event("C1"), APPROX, harmonic_modes=harmonic_modes)


def wf_mismatch(gen, x):
    """1 - normalised overlap (no maximisation) of cogwheel's h+ and hx against ours, 100 Mpc."""
    w = b.weight()
    hp, hc = e.hpc(*e.to_cart(x))
    par = to_par(x) | {"d_luminosity": 100.0}
    f = b._FREQS[b._BAND]
    cp, cc = gen.get_hplus_hcross(f, par)
    out = {}
    for lab, ours, cog in (("plus", hp[b._BAND], cp), ("cross", hc[b._BAND], cc)):
        ww = w[b._BAND]
        o = np.sum((ours * np.conj(cog)).real * ww)
        na, nb = np.sqrt(np.sum(np.abs(ours) ** 2 * ww)), np.sqrt(np.sum(np.abs(cog) ** 2 * ww))
        out[lab] = {"one_minus_overlap": float(1 - o / (na * nb)), "norm_ratio": float(nb / na)}
    return out


def no_multiband_mismatch(x):
    """Diagnosis: cogwheel vs lal XPHM without multibanding, and ours vs that (100 Mpc)."""
    import lal
    import lalsimulation as ls
    from cogwheel.waveform import compute_hplus_hcross

    mc, q, s1, s2, inc, phi = e.to_cart(np.asarray(x, float))
    m1, m2 = b.masses(mc, q)
    dic = lal.CreateDict()
    ls.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(dic, 0)
    ls.SimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(dic, 0)
    hp, hc = ls.SimInspiralChooseFDWaveform(
        m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s1, *s2, 1e8 * lal.PC_SI, inc, phi, 0, 0, 0,
        b.DF, b.F_LOW, b.F_HIGH, b.F_LOW, dic, ls.GetApproximantFromString(APPROX),
    )  # fmt: skip
    n, band, w = b._FREQS.size, b._BAND, b.weight()[b._BAND]
    ap, ac = hp.data.data[:n][band], hc.data.data[:n][band]
    cp, cc = compute_hplus_hcross(b._FREQS[band], to_par(x) | {"d_luminosity": 100.0}, APPROX)
    op, _ = e.hpc(mc, q, s1, s2, inc, phi)

    def mm(u, v):
        o = np.sum((u * np.conj(v)).real * w)
        return float(1 - o / np.sqrt(np.sum(np.abs(u) ** 2 * w) * np.sum(np.abs(v) ** 2 * w)))

    return {
        "cog_vs_noMB_plus": mm(ap, cp),
        "cog_vs_noMB_cross": mm(ac, cc),
        "ours_vs_noMB_plus": mm(op[band], ap),
    }


def injection_par(name, x_truth, psi_t, ra=1.0, dec=0.3):
    """Full cogwheel par_dic reproducing our target at the detector: t = cos2psi h+ + sin2psi hx."""
    from cogwheel.gw_utils import fplus_fcross, time_delay_from_geocenter

    f0 = fplus_fcross(DET, ra, dec, 0.0, TGPS)[:, 0]
    f4 = fplus_fcross(DET, ra, dec, np.pi / 4, TGPS)[:, 0]
    m = np.column_stack([f0, f4])  # F(psi) = cos2psi F(0) + sin2psi F(pi/4)
    v = np.linalg.solve(m, [np.cos(2 * psi_t), np.sin(2 * psi_t)])
    psi = 0.5 * np.arctan2(v[1], v[0]) % np.pi
    dl = 100.0 / np.linalg.norm(v)
    tg = -float(time_delay_from_geocenter(DET, ra, dec, TGPS)[0])
    return to_par(x_truth) | {
        "ra": ra,
        "dec": dec,
        "psi": float(psi),
        "t_geocenter": tg,
        "d_luminosity": float(dl),
    }


def rho2(name):
    ev = event(name)
    return float(4 * ev.df * np.sum(np.abs(ev.strain[0] * ev.wht_filter[0]) ** 2))


def checks():
    from cogwheel.likelihood import RelativeBinningLikelihood
    from cogwheel.waveform import APPROXIMANTS

    pts = check_points()
    f1 = {lab: box_diff(x, to_box(to_par(x))[0]) for lab, x in pts.items()}
    gen = wf_generator()
    f2 = {lab: wf_mismatch(gen, x) for lab, x in pts.items()}
    worst2 = max(v[p]["one_minus_overlap"] for v in f2.values() for p in ("plus", "cross"))
    f3 = {}
    for n, p in H.points().items():
        psi_t = p["spec"][0][6]
        par = injection_par(n, p["x_truth"], psi_t)
        ev = event(n)
        lk = RelativeBinningLikelihood(ev, wf_generator(), par_dic_0=par, pn_phase_tol=0.05)
        r2 = rho2(n)
        exact, rb = float(lk.lnlike_fft(par)), float(lk.lnlike(par))
        drift = float(lk.asd_drift[0])
        lk.asd_drift = [1.0]  # diagnosis (post hoc): cogwheel's ASD-drift estimate switched off
        exact1, rb1 = float(lk.lnlike_fft(par)), float(lk.lnlike(par))
        f3[n] = {
            "rho2": r2,
            "half_rho2_FF2": 0.5 * r2,
            "lnL_exact": exact,
            "lnL_rb": rb,
            "rel_exact": abs(exact - 0.5 * r2) / (0.5 * r2),
            "rel_rb": abs(rb - 0.5 * r2) / (0.5 * r2),
            "asd_drift_default": drift,
            "rel_exact_drift1": abs(exact1 - 0.5 * r2) / (0.5 * r2),
            "rel_rb_drift1": abs(rb1 - 0.5 * r2) / (0.5 * r2),
            "par": par,
        }
    f2b = {lab: no_multiband_mismatch(x) for lab, x in pts.items()}
    out = {
        "label": "derived; cogwheel 1.6.0 fixed-point checks against the X4 score (no search)",
        "versions": _versions(),
        "cogwheel_xphm_modes": APPROXIMANTS[APPROX].harmonic_modes,
        "points_box": {lab: list(map(float, x)) for lab, x in pts.items()},
        "F1_round_trip_max_abs": f1,
        "F1_pass": max(f1.values()) <= TOL_F1,
        "F2_waveform": f2,
        "F2_worst": worst2,
        "F2_pass": worst2 <= TOL_F2,
        "F2_diagnosis_post_hoc": {
            "what": "1 - overlap of cogwheel's h+ / hx against lal XPHM with multibanding off "
            "(PhenomXHM and PhenomXPHM ThresholdMband = 0), and of our default (multibanded) "
            "h+ against that",
            "points": f2b,
        },
        "F3_score": f3,
        "F3_pass": all(v["rel_exact"] <= TOL_F3 for v in f3.values()),
        "rules": {"F1": TOL_F1, "F2": TOL_F2, "F3": TOL_F3},
    }
    return out


def _versions():
    import cogwheel
    import lal
    import lalsimulation as ls
    import scipy

    return {
        "cogwheel": cogwheel.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "lal": lal.__version__,
        "lalsimulation": ls.__version__,
    }


def gate(path=None):
    """Names of the fixed-point checks (F1-F3) that did not pass; a search may start only if empty.
    Added after this round's runs were started by hand although F2 and F3 had failed."""
    path = OUT_DIR / "cogwheel_checks.json" if path is None else Path(path)
    if not path.exists():
        return ["no checks file"]
    out = json.loads(path.read_text())
    return [f for f in ("F1", "F2", "F3") if out.get(f"{f}_pass") is not True]


def run_one(name, seed):
    """Published workflow: Posterior.from_event (RWF aligned-spin start) then one DE refinement."""
    import cogwheel.posterior as cp
    from cogwheel.posterior import Posterior

    rec = {}
    real_de = cp.optimize.differential_evolution

    def de_spy(*a, **kw):  # records SciPy's result; passes everything through unchanged
        r = real_de(*a, **kw)
        rec.update({"nfev": int(r.nfev), "nit": int(r.nit), "message": str(r.message)})
        rec["success"] = bool(r.success)
        return r

    cp.optimize.differential_evolution = de_spy
    t0 = time.time()
    ev = event(name)
    post = Posterior.from_event(
        ev,
        mchirp_guess=e.MC0,
        approximant=APPROX,
        prior_class="IASPrior",
        # deviation: the frozen q_min = 1/15 failed (C2's aligned-spin start had q 17.9, outside
        # the prior cube); cogwheel's default q_min 0.05 (q <= 20) contains our box instead
        ref_wf_finder_kwargs={"mchirp_range": (0.4 * e.MC0, 1.6 * e.MC0), "f_ref": F_REF},
    )
    t_rwf = time.time() - t0
    par0 = dict(post.likelihood.par_dic_0)
    lnl0 = float(post.likelihood.lnlike(par0))
    n_calls = [0]
    real_lnlike = post.likelihood.lnlike

    def counted(par_dic):
        n_calls[0] += 1
        return real_lnlike(par_dic)

    post.likelihood.lnlike = counted
    post.refine_reference_waveform(rng=seed)
    post.likelihood.lnlike = real_lnlike
    par = dict(post.likelihood.par_dic_0)
    out = {
        "name": name,
        "seed": seed,
        "start": "cogwheel ReferenceWaveformFinder (aligned-spin XAS, (2,2)); "
        "no truth, no saved solution",
        "rwf_s": round(t_rwf, 1),
        "lnL_rb_start": lnl0,
        "lnL_rb_final": float(post.likelihood.lnlike(par)),
        "lnL_exact_final": float(post.likelihood.lnlike_fft(par)),
        "asd_drift": [float(a) for a in post.likelihood.asd_drift],
        "lnlike_calls": n_calls[0],
        "de": rec,
        "wall_s": round(time.time() - t0, 1),
        "par_start": par0,
        "par_final": par,
    }
    return score(name, out)


def score(name, out):
    t, _ = k.target(name)
    for lab in ("start", "final"):
        x, inside, mags = to_box(out[f"par_{lab}"])
        mc, q, s1, s2, inc, phi = physical(out[f"par_{lab}"])
        hp, hc = e.hpc(mc, q, s1, s2, inc, phi)  # unclipped spins: the point cogwheel found
        ff = e.stat(t, hp, hc, pad=64)
        out[lab] = {
            "x_box": list(map(float, x)),
            "inside_box": inside,
            "chi_mags": mags,
            "q": q,
            "FF_pad64": ff,
            "eps": 1 - ff,
        }
        if not inside:  # projection onto |chi| <= 0.99 (reported, not used for the pass)
            out[lab]["FF_pad64_projected"] = e.stat(t, *e.hpc(*e.to_cart(x)), pad=64)
    r2 = rho2(name)
    lnl1 = out["lnL_exact_final"] * out["asd_drift"][0] ** 2  # lnL at asd_drift = 1
    out["lnL_exact_final_drift1"] = lnl1
    out["FF_implied_by_exact_lnL"] = float(np.sqrt(max(0.0, 2 * lnl1 / r2)))
    eps = out["final"]["eps"]
    out["pass"] = bool(
        out["final"]["inside_box"] and eps <= EPS_PASS and eps <= REL_PASS * (1 - FF_CASE)
    )
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    cmd = argv[0] if argv else ""
    if cmd == "checks":
        out = checks()
        (OUT_DIR / "cogwheel_checks.json").write_text(
            json.dumps(out, indent=1, default=float) + "\n"
        )
        print(json.dumps({k_: out[k_] for k_ in out if k_.endswith(("pass", "worst"))}))
    elif cmd in ("one", "runs") and gate():
        print(f"refused: fixed-point checks not passed: {gate()}", file=sys.stderr)
        return 2
    elif cmd == "one":
        r = run_one(argv[1], int(argv[2]))
        (OUT_DIR / f".cogwheel_{argv[1]}.json").write_text(json.dumps(r, default=float) + "\n")
        print(r["name"], r["final"]["eps"], r["de"], r["wall_s"], flush=True)
    elif cmd == "pycbc":  # PyCBC environment: the #447 check at cogwheel's start and final points
        p = OUT_DIR / "cogwheel_runs.json"
        runs = json.loads(p.read_text())
        for r in runs["rows"]:
            if "final" not in r:
                continue
            t, _ = k.target(r["name"])
            for lab in ("start", "final"):
                mc, q, s1, s2, inc, phi = physical(r[f"par_{lab}"])
                ph = {"mc": mc, "q": q, "chi1": s1, "chi2": s2, "inc": inc, "phi": phi}
                ff = H.pycbc_ff(t, ph, pad=64)
                r[lab]["FF_pycbc"] = ff
                r[lab]["pycbc_minus_ours"] = ff - r[lab]["FF_pad64"]
        p.write_text(json.dumps(runs, indent=1, default=float) + "\n")
    elif cmd == "runs":
        procs = {
            n: subprocess.Popen([sys.executable, __file__, "one", n, str(sd)])
            for n, sd in SEEDS.items()
        }
        rows, t0 = [], time.time()
        for n, p in procs.items():
            try:
                code = p.wait(timeout=max(1.0, CAP_S - (time.time() - t0)))
            except subprocess.TimeoutExpired:
                p.kill()
                code = "killed at the 30-min cap"
            f = OUT_DIR / f".cogwheel_{n}.json"
            if code == 0 and f.exists():
                rows.append(json.loads(f.read_text()))
                f.unlink()
            else:
                rows.append({"name": n, "seed": SEEDS[n], "status": f"not finished ({code})"})
        out = {
            "label": "simulated; cogwheel 1.6.0 refine_reference_waveform on C1, C2 (development "
            "comparison on already-seen cases, not a validation)",
            "versions": _versions(),
            "rows": rows,
        }
        (OUT_DIR / "cogwheel_runs.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    else:
        print(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
