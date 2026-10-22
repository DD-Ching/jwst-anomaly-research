"""X4 cogwheel exact-score round: fixed-point checks of a thin exact-likelihood interface.

No search. Unit file: docs/state/units/X4.md (cogwheel exact-score round; plan frozen before this
script ran). `simulated`/`derived`; no data. Both sides use lal IMRPhenomXPHM with multibanding off.
Needs cogwheel 1.6.0, lalsuite 7.26.16 and mpmath (cogwheel environment). Every point is
scored only; nothing is optimised.

  python scripts/x4_cog_exact.py checks   # -> results/x4/cog_exact/checks.json (frozen)
  python scripts/x4_cog_exact.py recheck  # -> results/x4/cog_exact/review_recheck.json (post hoc)
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.fft as sf

_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_DIR))
import x4_cogwheel as C  # noqa: E402

H, k, e, b = C.H, C.k, C.e, C.b
OUT_DIR = _DIR.parent / "results" / "x4" / "cog_exact"
PAD = 64
DT0 = 2.0  # half-width of cogwheel's t_refdet window: covers the 4 s cyclic period
Q_MIN = 1 / 15
CHI_MAX = 0.99
CHI_TOL = 1e-12  # post hoc (review): rounding of |chi| = 0.99 through the spin rotation
BOX_SEEDS = range(1001, 1013)
TOL_F1, TOL_F2, TOL_F3 = 1e-8, 1e-6, 1e-6
RA, DEC = 1.0, 0.3


def nomb_dict():
    import lal
    import lalsimulation as ls

    dic = lal.CreateDict()
    ls.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(dic, 0)
    ls.SimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(dic, 0)
    return dic


def hpc(mc, q, s1, s2, inc, phi):
    """Our h+, hx (100 Mpc) as in x4_s1d_de.hpc, with multibanding off."""
    import lal
    import lalsimulation as ls

    m1, m2 = b.masses(mc, q)
    hp, hc = ls.SimInspiralChooseFDWaveform(
        m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s1, *s2, 1e8 * lal.PC_SI, inc, phi, 0, 0, 0,
        b.DF, b.F_LOW, b.F_HIGH, b.F_LOW, nomb_dict(), ls.GetApproximantFromString(C.APPROX),
    )  # fmt: skip
    n = min(hp.data.length, b._FREQS.size)
    a = np.zeros(b._FREQS.size, complex)
    c = np.zeros(b._FREQS.size, complex)
    a[:n], c[:n] = hp.data.data[:n], hc.data.data[:n]
    a[~b._BAND] = 0.0
    c[~b._BAND] = 0.0
    return a, c


def target(name):
    (mc, q, s1, s2, inc, phi, psi), _ = k.c_spec(name)[0]
    hp, hc = hpc(mc, q, tuple(s1), tuple(s2), inc, phi)
    return np.cos(2 * psi) * hp + np.sin(2 * psi) * hc


def stat_argmax(t, hp, hc, pad=PAD):
    """e.stat (non-degenerate branch) plus its arg max: lag index, n, and (a, c) = G^-1 v there."""
    w = b.weight()
    n = sf.next_fast_len(pad * t.size)
    vp = (sf.ifft(t * np.conj(hp) * w, n) * n).real
    vc = (sf.ifft(t * np.conj(hc) * w, n) * n).real
    g = np.array(
        [
            [np.sum(np.abs(hp) ** 2 * w), np.sum((hp * np.conj(hc)).real * w)],
            [np.sum((hp * np.conj(hc)).real * w), np.sum(np.abs(hc) ** 2 * w)],
        ]
    )
    det = g[0, 0] * g[1, 1] - g[0, 1] ** 2
    assert det > 1e-12 * g[0, 0] * g[1, 1]
    nt = np.sum(np.abs(t) ** 2 * w)
    qv = (g[1, 1] * vp * vp - 2 * g[0, 1] * vp * vc + g[0, 0] * vc * vc) / det
    i = int(np.argmax(qv))
    ac = np.linalg.solve(g, [vp[i], vc[i]])
    return float(np.sqrt(max(qv[i], 0.0) / nt)), i, n, ac, float(ac @ [vp[i], vc[i]])


def points():
    lo = np.array([bd[0] for bd in e.BOUNDS])
    hi = np.array([bd[1] for bd in e.BOUNDS])
    out = {}
    for s in BOX_SEEDS:
        x = lo + np.random.default_rng(s).random(lo.size) * (hi - lo)
        x[1], x[2], x[5] = min(x[1], 14.9), min(x[2], 0.98), min(x[5], 0.98)
        out[f"box_{s}"] = x
    mid = 0.5 * (lo + hi)
    p = {lab: mid.copy() for lab in ("q1499_chi0989", "tilts_pm0999", "inc_001", "inc_pi_m001")}
    p["q1499_chi0989"][[1, 2, 5]] = 14.99, 0.989, 0.989
    p["tilts_pm0999"][[2, 3, 5, 6]] = 0.98, 0.999, 0.98, -0.999
    p["inc_001"][[2, 5, 8]] = 0.8, 0.8, 0.01
    p["inc_pi_m001"][[2, 5, 8]] = 0.8, 0.8, np.pi - 0.01
    out.update(p)
    return out


def event(name, t):
    from cogwheel.data import EventData

    return EventData(
        eventname=f"x4ex_{name}",
        frequencies=b._FREQS.copy(),
        strain=t[np.newaxis],
        wht_filter=np.sqrt(b.weight())[np.newaxis],
        detector_names=C.DET,
        tgps=C.TGPS,
        tcoarse=0.0,
    )


def over_cap(par_dic):
    """True if a spin exceeds |chi| 0.99. Tolerance added after review: without it, 730 of 2000
    points at |chi| = 0.99 exactly were rejected by rounding in the N-frame rotation."""
    return any(
        np.hypot(np.hypot(par_dic[f"s{i}x_n"], par_dic[f"s{i}y_n"]), par_dic[f"s{i}z"])
        > CHI_MAX + CHI_TOL
        for i in (1, 2)
    )


def exact_likelihood_class():
    """The thin interface: cogwheel's exact FFT-grid lnL as `lnlike`, ASD drift 1, our spin cap."""
    from cogwheel.likelihood import CBCLikelihood

    class ExactLikelihood(CBCLikelihood):
        def __init__(self, event_data, waveform_generator, par_dic_0):
            super().__init__(event_data, waveform_generator)
            self.asd_drift = np.ones(len(event_data.detector_names))
            self.par_dic_0 = par_dic_0

        def lnlike(self, par_dic):
            if over_cap(par_dic):
                return -np.inf
            return self.lnlike_fft(par_dic)

    return ExactLikelihood


def prior(f_avg):
    from cogwheel.gw_prior import IASPrior

    return IASPrior(
        tgps=C.TGPS, ref_det_name=C.DET, detector_pair=C.DET, f_avg=f_avg, f_ref=C.F_REF,
        mchirp_range=(0.4 * e.MC0, 1.6 * e.MC0), q_min=Q_MIN, t0_refdet=0.0, dt0=DT0,
    )  # fmt: skip


def extrinsic(ac, i, n):
    """(a, c) at our lag i -> psi, d_luminosity, t_geocenter at the fixed sky (one detector)."""
    from cogwheel.gw_utils import fplus_fcross, time_delay_from_geocenter

    f0 = fplus_fcross(C.DET, RA, DEC, 0.0, C.TGPS)[:, 0]
    f4 = fplus_fcross(C.DET, RA, DEC, np.pi / 4, C.TGPS)[:, 0]
    v = np.linalg.solve(np.column_stack([f0, f4]), ac)  # F(psi) = cos2psi F(0) + sin2psi F(pi/4)
    period = 1 / b.DF
    t = i * period / n
    t = t - period if t > period / 2 else t
    tg = t - float(time_delay_from_geocenter(C.DET, RA, DEC, C.TGPS)[0])
    return {
        "ra": RA,
        "dec": DEC,
        "psi": float(0.5 * np.arctan2(v[1], v[0]) % np.pi),
        "d_luminosity": float(100.0 / np.linalg.norm(v)),
        "t_geocenter": tg,
        "t_refdet": t,
    }


def check_point(lk, pri, gen, t, x):
    hp, hc = hpc(*e.to_cart(x))
    par_i = C.to_par(x)
    out = {"x_box": list(map(float, x))}
    out["F1"] = C.box_diff(x, C.to_box(par_i)[0])
    f = b._FREQS[b._BAND]
    cp, cc = gen.get_hplus_hcross(f, par_i | {"d_luminosity": 100.0})
    w = b.weight()[b._BAND]
    mm = []
    for ours, cog in ((hp[b._BAND], cp), (hc[b._BAND], cc)):
        o = np.sum((ours * np.conj(cog)).real * w)
        mm.append(
            float(1 - o / np.sqrt(np.sum(np.abs(ours) ** 2 * w) * np.sum(np.abs(cog) ** 2 * w)))
        )
    out["F2"] = max(mm)
    ff, i, n, ac, a_v = stat_argmax(t, hp, hc)
    assert abs(ff - e.stat(t, hp, hc, pad=PAD)) <= 1e-12
    ext = extrinsic(ac, i, n)
    par = par_i | {key: ext[key] for key in ("ra", "dec", "psi", "d_luminosity", "t_geocenter")}
    h_f = lk._get_h_f(par)
    dh = float(lk._compute_d_h(h_f)[0].real)
    hh = float(lk._compute_h_h(h_f)[0])
    dd = float(
        4
        * lk.event_data.df
        * np.sum(np.abs(lk.event_data.strain[0] * lk.event_data.wht_filter[0]) ** 2)
    )
    ff_cog = dh / np.sqrt(dd * hh)
    lnl = float(lk.lnlike(par))
    out.update(
        {
            "FF_ours": ff,
            "FF_cog": float(ff_cog),
            "F3": float(abs(ff - ff_cog)),
            "lag_index": i,
            "t_refdet": ext["t_refdet"],
            "psi": ext["psi"],
            "d_luminosity": ext["d_luminosity"],
            "dh_over_4df_minus_ours": float(dh / (4 * b.DF) - a_v),
            "dh_rel_diff": float((dh / (4 * b.DF) - a_v) / a_v),
            "two_lnL_over_rho2_minus_FF2": float(2 * lnl / dd - ff * ff),
        }
    )
    sampled = pri.inverse_transform(**(par | {"f_ref": C.F_REF, "l1": 0.0, "l2": 0.0}))
    vals = np.array([sampled[p_] for p_ in pri.sampled_params])
    in_cube = bool(np.all(vals >= pri.cubemin) and np.all(vals <= pri.cubemin + pri.cubesize))
    back = pri.transform(**sampled)
    out["F1b_in_cube"] = in_cube
    cube = zip(pri.sampled_params, vals, pri.cubemin, pri.cubesize, strict=True)
    out["F1b_out_of_cube_params"] = [p_ for p_, v_, lo_, s_ in cube if not lo_ <= v_ <= lo_ + s_]
    out["d_hat"] = float(sampled["d_hat"])
    out["F1b_round_trip"] = C.box_diff(x, C.to_box(back)[0])
    return out


def checks():
    from cogwheel.posterior import Posterior
    from cogwheel.waveform import WaveformGenerator

    pts = H.points()
    fixed = points()
    out = {
        "label": "derived; X4 cogwheel exact-score round: fixed-point checks, no optimisation; "
        "lal XPHM with multibanding off on both sides; comparison targets rebuilt from the "
        "C1/C2 specs",
        "versions": C._versions(),
        "rules": {"F1": TOL_F1, "F2": TOL_F2, "F3": TOL_F3, "F1b": "in cube and round trip <= F1"},
        "settings": {"pad": PAD, "dt0_s": DT0, "q_min": Q_MIN, "chi_max": CHI_MAX, "ra": RA,
                     "dec": DEC, "asd_drift": 1.0,
                     "likelihood": "CBCLikelihood.lnlike_fft (no relative binning)"},
        "cases": {},
    }  # fmt: skip
    for name in ("C1", "C2"):
        t = target(name)
        ev = event(name, t)
        gen = WaveformGenerator.from_event_data(ev, C.APPROX)
        x_mid = np.array([0.5 * (bd[0] + bd[1]) for bd in e.BOUNDS])
        lk = exact_likelihood_class()(ev, gen, C.to_par(x_mid) | extrinsic([1.0, 0.0], 0, 1))
        pri = prior(f_avg=float(lk.get_average_frequency(lk.par_dic_0)))
        post = Posterior(pri, lk)  # constructs; refine_reference_waveform is not called
        case_pts = {"truth": pts[name]["x_truth"], "found_s1k": pts[name]["x_found"]} | fixed
        rows = {lab: check_point(lk, pri, gen, t, x) for lab, x in case_pts.items()}
        par_t = C.to_par(pts[name]["x_truth"]) | extrinsic([1.0, 0.0], 0, 1)
        t0 = time.perf_counter()
        for _ in range(20):
            lk.lnlike(par_t)
        sec_per_call = (time.perf_counter() - t0) / 20
        # grid loss of our objective (reported, not compared): pad 64 vs pad 512 at the truth
        hp, hc = hpc(*e.to_cart(pts[name]["x_truth"]))
        out["cases"][name] = {
            "posterior_constructed": repr(post),
            "exact_lnlike_s_per_call": sec_per_call,
            "rho2": float(4 * b.DF * np.sum(np.abs(t) ** 2 * b.weight())),
            "grid_loss_truth_pad64_vs_512": e.stat(t, hp, hc, pad=512) - e.stat(t, hp, hc, pad=PAD),
            "points": rows,
        }
    allrows = [r for c in out["cases"].values() for r in c["points"].values()]
    out["worst"] = {
        "F1": max(r["F1"] for r in allrows),
        "F1b_round_trip": max(r["F1b_round_trip"] for r in allrows),
        "F2": max(r["F2"] for r in allrows),
        "F3": max(r["F3"] for r in allrows),
        "dh_rel_diff_abs": max(abs(r["dh_rel_diff"]) for r in allrows),
    }
    out["pass"] = {
        "F1": out["worst"]["F1"] <= TOL_F1,
        "F1b": all(r["F1b_in_cube"] for r in allrows) and out["worst"]["F1b_round_trip"] <= TOL_F1,
        "F2": out["worst"]["F2"] <= TOL_F2,
        "F3": out["worst"]["F3"] <= TOL_F3,
    }
    return out


def recheck():
    """Post hoc, after the independent review (the frozen checks.json is not rewritten): the
    spin-cap edge, exact-boundary points, the prior-transform and unfold path a search would call,
    and our grid loss (pad 64 vs 512) at every frozen point. Scores only; nothing is optimised."""
    from cogwheel.waveform import WaveformGenerator

    lo = np.array([bd[0] for bd in e.BOUNDS])
    hi = np.array([bd[1] for bd in e.BOUNDS])
    rng = np.random.default_rng(2001)
    n_rej = 0
    for _ in range(2000):
        x = lo + rng.random(lo.size) * (hi - lo)
        x[2] = x[5] = CHI_MAX
        n_rej += over_cap(C.to_par(x))
    mid = 0.5 * (lo + hi)
    edge = {lab: mid.copy() for lab in ("q15_chi099", "tilts_pm1_chi099")}
    edge["q15_chi099"][[1, 2, 5]] = 15.0, CHI_MAX, CHI_MAX
    edge["tilts_pm1_chi099"][[2, 3, 5, 6]] = CHI_MAX, 1.0, CHI_MAX, -1.0
    out = {
        "label": "derived; post hoc review re-check of the X4 cogwheel exact-score round (not part "
        "of the frozen checks; checks.json unchanged)",
        "versions": C._versions(),
        "chi_tol": CHI_TOL,
        "cap_edge_rejected_of_2000": int(n_rej),
        "cases": {},
    }
    pts = H.points()
    for name in ("C1", "C2"):
        t = target(name)
        ev = event(name, t)
        gen = WaveformGenerator.from_event_data(ev, C.APPROX)
        x_mid = np.array([0.5 * (bd[0] + bd[1]) for bd in e.BOUNDS])
        lk = exact_likelihood_class()(ev, gen, C.to_par(x_mid) | extrinsic([1.0, 0.0], 0, 1))
        pri = prior(f_avg=float(lk.get_average_frequency(lk.par_dic_0)))
        unfold_max = pri.unfold_apply(lambda *v, lk=lk, pri=pri: lk.lnlike(pri.transform(*v)))
        frozen = {"truth": pts[name]["x_truth"], "found_s1k": pts[name]["x_found"]} | points()
        rows = {}
        for lab, x in (frozen | edge).items():
            hp, hc = hpc(*e.to_cart(x))
            ff, i, n, ac, _ = stat_argmax(t, hp, hc)
            ext = extrinsic(ac, i, n)
            par = C.to_par(x) | {
                key: ext[key] for key in ("ra", "dec", "psi", "d_luminosity", "t_geocenter")
            }
            lnl = float(lk.lnlike(par))
            sampled = pri.inverse_transform(**(par | {"f_ref": C.F_REF, "l1": 0.0, "l2": 0.0}))
            vals = [sampled[p_] for p_ in pri.sampled_params]
            lnl_path = float(lk.lnlike(pri.transform(*vals)))
            folded = pri.fold(*vals)
            r = {
                "post_hoc_edge_point": lab in edge,
                "FF_ours": ff,
                "lnL": lnl,
                "lnL_via_prior_transform_rel": (lnl_path - lnl) / abs(lnl),
                "max_over_unfolds_minus_lnL": float(max(unfold_max(*folded))) - lnl,
                "grid_loss_pad64_vs_512": e.stat(t, hp, hc, pad=512) - ff,
            }
            if lab in edge:
                r |= check_point(lk, pri, gen, t, x)
            rows[lab] = r
        out["cases"][name] = rows
    allr = [r for c in out["cases"].values() for r in c.values()]
    edge_r = [r for r in allr if r["post_hoc_edge_point"]]
    out["worst"] = {
        "lnL_via_prior_transform_rel": max(abs(r["lnL_via_prior_transform_rel"]) for r in allr),
        "max_over_unfolds_minus_lnL": max(r["max_over_unfolds_minus_lnL"] for r in allr),
        "grid_loss_pad64_vs_512": max(r["grid_loss_pad64_vs_512"] for r in allr),
        "edge_F1": max(r["F1"] for r in edge_r),
        "edge_F2": max(r["F2"] for r in edge_r),
        "edge_F3": max(r["F3"] for r in edge_r),
        "edge_lnL_finite": all(np.isfinite(r["lnL"]) for r in edge_r),
        "edge_in_cube": all(r["F1b_in_cube"] for r in edge_r),
    }
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["checks"]:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        out = checks()
        (OUT_DIR / "checks.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
        print(json.dumps({"worst": out["worst"], "pass": out["pass"]}))
    elif argv[:1] == ["recheck"]:
        out = recheck()
        (OUT_DIR / "review_recheck.json").write_text(
            json.dumps(out, indent=1, default=float) + "\n"
        )
        print(json.dumps(out["worst"] | {"cap_edge_rejected": out["cap_edge_rejected_of_2000"]}))
    else:
        print(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
