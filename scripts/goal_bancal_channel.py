"""When can a remote receiver read the sender's bit in a Bancal-type v-causal model?

Goal round 2 (docs/goal_round2/README.md). Model class: Bancal, Pironio, Acin, Liang, Scarani &
Gisin, Nature Phys. 8, 867 (2012), arXiv:1110.3795 (v2), "v-causal models". All offline:

1. Positive control: the paper's 4-qubit state and measurements (App. B) give S = 7.2014 > 7
   for the Bell expression of Lemma 1 (Eq. 4). Standard quantum mechanics (QM) is the
   no-signalling baseline.
2. Signalling magnitude (linear programme). In configuration R (A < D < (B ~ C)), App. C makes
   the BC correlations conditioned on (a, x, d, w) local, and App. D makes the ABD and ACD
   marginals equal to the quantum ones. The programme finds the smallest setting dependence of
   the BCD marginal on x (signalling A -> D') and of the ABC marginal on w (signalling D -> A')
   that any distribution with those properties must have, and whether either alone can vanish.
3. Event timing: the paper's Fig. 3 layout and other layouts meeting conditions i)-iii) of
   App. E, with sender A and D a distance L apart. For each, the earliest time at which the
   B, C, D records (A -> D') or A, B, C records (D -> A') can be together at the receiver's
   site, counted from the sender's measurement.

Units: c = 1, distances in L, times in L / c, influence speed v = r c. Outcome 0 <-> +1.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass

import numpy as np

BITS = (0, 1)
SX = np.array([[0.0, 1.0], [1.0, 0.0]])
SZ = np.array([[1.0, 0.0], [0.0, -1.0]])
HAD = (SX + SZ) / math.sqrt(2.0)
U_ROT = math.cos(4 * math.pi / 5) * SZ - math.sin(4 * math.pi / 5) * SX

# App. B: |Psi> and the observables (U is real, so U^dagger = U^T).
PSI_AMPLITUDES = {
    "0000": 17 / 60,
    "0011": 1 / 3,
    "0101": -1 / math.sqrt(8),
    "0110": 1 / 10,
    "1000": 1 / 4,
    "1011": -1 / 2,
    "1101": -1 / 3,
    "1110": 1 / 2,
}
OBS = {
    "A": (-U_ROT @ SX @ U_ROT.T, U_ROT @ SZ @ U_ROT.T),
    "B": (HAD, -SX @ HAD @ SX),
    "C": (SZ, -SX),
    "D": (-SZ, -SX),
}

# Lemma 1, Eq. (4): coefficient, then (party, setting) pairs of the correlator.
S_TERMS = [
    (-3, (("A", 0),)), (-1, (("B", 0),)), (-1, (("B", 1),)), (-1, (("C", 0),)), (-3, (("D", 0),)),
    (-1, (("A", 1), ("B", 0))), (-1, (("A", 1), ("B", 1))), (1, (("A", 0), ("C", 0))),
    (2, (("A", 1), ("C", 0))), (1, (("A", 0), ("D", 0))), (1, (("B", 0), ("D", 1))),
    (-1, (("B", 1), ("D", 1))), (-1, (("C", 0), ("D", 0))), (-2, (("C", 1), ("D", 1))),
    (1, (("A", 0), ("B", 0), ("D", 0))), (1, (("A", 0), ("B", 0), ("D", 1))),
    (1, (("A", 0), ("B", 1), ("D", 0))), (-1, (("A", 0), ("B", 1), ("D", 1))),
    (-1, (("A", 1), ("B", 0), ("D", 0))), (-1, (("A", 1), ("B", 1), ("D", 0))),
    (1, (("A", 0), ("C", 0), ("D", 0))), (2, (("A", 1), ("C", 0), ("D", 0))),
    (-2, (("A", 0), ("C", 1), ("D", 1))),
]  # fmt: skip
S_LOCAL_BOUND = 7.0


def quantum_state() -> np.ndarray:
    psi = np.zeros(16)
    for k, v in PSI_AMPLITUDES.items():
        psi[int(k, 2)] = v
    return psi


def _projector(obs: np.ndarray, outcome: int) -> np.ndarray:
    return 0.5 * (np.eye(2) + (1 - 2 * outcome) * obs)


def quantum_distribution() -> np.ndarray:
    """P[a, b, c, d, x, y, z, w] of the App. B state and observables (standard QM)."""
    psi = quantum_state()
    p = np.zeros((2,) * 8)
    for a, b, c, d, x, y, z, w in itertools.product(BITS, repeat=8):
        op = np.kron(
            np.kron(_projector(OBS["A"][x], a), _projector(OBS["B"][y], b)),
            np.kron(_projector(OBS["C"][z], c), _projector(OBS["D"][w], d)),
        )
        p[a, b, c, d, x, y, z, w] = psi @ op @ psi
    return p


def bell_value(p: np.ndarray) -> float:
    """S of Eq. (4), with each correlator taken from P averaged over the settings of absent parties.
    For a no-signalling P (QM) this is the unique marginal."""
    idx = {"A": 0, "B": 1, "C": 2, "D": 3}
    total = 0.0
    for coef, parts in S_TERMS:
        val = 0.0
        settings = {idx[k]: s for k, s in parts}
        for out in itertools.product(BITS, repeat=4):
            sign = (-1) ** sum(out[idx[k]] for k, _ in parts)
            # average over the free settings of absent parties
            free = [i for i in range(4) if i not in settings]
            acc = 0.0
            for fs in itertools.product(BITS, repeat=len(free)):
                st = [0, 0, 0, 0]
                for i, s in settings.items():
                    st[i] = s
                for i, s in zip(free, fs, strict=True):
                    st[i] = s
                acc += p[out + tuple(st)]
            val += sign * acc / (2 ** len(free))
        total += coef * val
    return total


def quantum_is_no_signalling(p: np.ndarray, tol: float = 1e-12) -> bool:
    """True if, for every party k, the outcome marginal of the other three does not depend on k's
    setting (axes 0-3 outcomes a, b, c, d; axes 4-7 settings x, y, z, w)."""
    for k in range(4):
        m = p.sum(axis=k)  # setting axis 4 + k becomes 3 + k
        if not np.allclose(m.take(0, axis=3 + k), m.take(1, axis=3 + k), atol=tol, rtol=0):
            return False
    return True


# ----------------------------------------------------------------------------- signalling LP


@dataclass(frozen=True)
class LPResult:
    status: str
    min_bcd_on_x: float | None  # mean total-variation dependence of P(bcd|xyzw) on x
    min_abc_on_w: float | None  # mean total-variation dependence of P(abc|xyzw) on w
    min_of_max: float | None  # min over models of max(the two)
    min_of_sum: float | None
    min_bcd_on_x_if_abc_silent: float | None  # A -> D' signal forced when D -> A' carries none
    min_abc_on_w_if_bcd_silent: float | None
    ns_feasible: bool  # can a fully no-signalling P satisfy the constraints? (Lemma 1: no)
    min_max_per_setting: float | None  # min over models of the largest TV over direction, settings


def _vid(a, b, c, d, x, y, z, w) -> int:
    return ((((((a * 2 + b) * 2 + c) * 2 + d) * 2 + x) * 2 + y) * 2 + z) * 2 + w


def signalling_lp(pq: np.ndarray | None = None) -> LPResult:
    """Linear programme over P(abcd|xyzw) >= 0 (256 variables) with:
    - normalisation for each setting;
    - ABD marginal = quantum ABD for every z, ACD marginal = quantum ACD for every y (App. D);
    - for each (a, x, d, w): the BC correlations conditioned on AD are local (App. C). Given the two
      marginal conditions, this is the 8 CHSH facets in joint form,
      |sum_yz s_yz E_yz(a,x,d,w)| <= 2 P(ad|xw), for sign patterns with an odd number of minuses.
    Objectives: the mean over the remaining settings of the total-variation distance of the BCD
    marginal between x = 0 and x = 1, and of the ABC marginal between w = 0 and w = 1."""
    from scipy.optimize import linprog

    pq = quantum_distribution() if pq is None else pq
    n_p = 256
    # auxiliary |difference| variables: BCD diffs (b,c,d,y,z,w) and ABC diffs (a,b,c,x,y,z)
    n_dx, n_dw = 64, 64
    n = n_p + n_dx + n_dw + 1  # + t for min-max
    eq_rows, eq_b, ub_rows, ub_b = [], [], [], []

    def row():
        return np.zeros(n)

    for x, y, z, w in itertools.product(BITS, repeat=4):
        r = row()
        for a, b, c, d in itertools.product(BITS, repeat=4):
            r[_vid(a, b, c, d, x, y, z, w)] = 1.0
        eq_rows.append(r)
        eq_b.append(1.0)
    for a, b, d, x, y, w, z in itertools.product(BITS, repeat=7):
        r = row()
        for c in BITS:
            r[_vid(a, b, c, d, x, y, z, w)] = 1.0
        eq_rows.append(r)
        eq_b.append(pq[a, b, :, d, x, y, z, w].sum())
    for a, c, d, x, z, w, y in itertools.product(BITS, repeat=7):
        r = row()
        for b in BITS:
            r[_vid(a, b, c, d, x, y, z, w)] = 1.0
        eq_rows.append(r)
        eq_b.append(pq[a, :, c, d, x, y, z, w].sum())
    signs = [s for s in itertools.product((1, -1), repeat=4) if s.count(-1) % 2 == 1]
    for a, d, x, w in itertools.product(BITS, repeat=4):
        p_ad = pq[a, :, :, d, x, 0, 0, w].sum()
        for s in signs:  # the 8 odd-minus patterns already contain each facet's negation
            r = row()
            for k, (y, z) in enumerate(itertools.product(BITS, repeat=2)):
                for b, c in itertools.product(BITS, repeat=2):
                    r[_vid(a, b, c, d, x, y, z, w)] += s[k] * (-1) ** (b + c)
            ub_rows.append(r)
            ub_b.append(2.0 * p_ad)
    # |BCD(x=0) - BCD(x=1)| <= u, |ABC(w=0) - ABC(w=1)| <= u
    k = n_p
    for b, c, d, y, z, w in itertools.product(BITS, repeat=6):
        for sgn in (1, -1):
            r = row()
            for a in BITS:
                r[_vid(a, b, c, d, 0, y, z, w)] += sgn
                r[_vid(a, b, c, d, 1, y, z, w)] -= sgn
            r[k] = -1.0
            ub_rows.append(r)
            ub_b.append(0.0)
        k += 1
    for a, b, c, x, y, z in itertools.product(BITS, repeat=6):
        for sgn in (1, -1):
            r = row()
            for d in BITS:
                r[_vid(a, b, c, d, x, y, z, 0)] += sgn
                r[_vid(a, b, c, d, x, y, z, 1)] -= sgn
            r[k] = -1.0
            ub_rows.append(r)
            ub_b.append(0.0)
        k += 1
    # mean TV over the 8 setting triples = (1/2) * sum |diff| / 8
    tv_x = np.zeros(n)
    tv_x[n_p : n_p + n_dx] = 0.5 / 8
    tv_w = np.zeros(n)
    tv_w[n_p + n_dx : n_p + n_dx + n_dw] = 0.5 / 8
    a_eq, b_eq = np.array(eq_rows), np.array(eq_b)
    a_ub, b_ub = np.array(ub_rows), np.array(ub_b)
    bounds = [(0, None)] * n

    def solve(c, extra_ub=(), extra_b=(), extra_eq=(), extra_eqb=()):
        aub = np.vstack([a_ub, *extra_ub]) if extra_ub else a_ub
        bub = np.concatenate([b_ub, list(extra_b)]) if extra_ub else b_ub
        aeq = np.vstack([a_eq, *extra_eq]) if extra_eq else a_eq
        beq = np.concatenate([b_eq, list(extra_eqb)]) if extra_eq else b_eq
        return linprog(c, A_ub=aub, b_ub=bub, A_eq=aeq, b_eq=beq, bounds=bounds, method="highs")

    r1 = solve(tv_x)
    r2 = solve(tv_w)
    t = np.zeros(n)
    t[-1] = 1.0
    r3 = solve(t, extra_ub=[tv_x - t, tv_w - t], extra_b=[0.0, 0.0])
    r4 = solve(tv_x + tv_w)
    r6 = solve(tv_x, extra_ub=[tv_w], extra_b=[0.0])
    r7 = solve(tv_w, extra_ub=[tv_x], extra_b=[0.0])
    # no-signalling feasibility: all diffs forced to zero
    zero = np.zeros(n)
    rows = []
    for j in range(n_p, n_p + n_dx + n_dw):
        e = np.zeros(n)
        e[j] = 1.0
        rows.append(e)
    r5 = solve(zero, extra_eq=rows, extra_eqb=[0.0] * len(rows))
    # largest TV over (direction, setting triple): diff index j = outcome * 8 + triple
    per_setting = []
    for base, size in ((n_p, n_dx), (n_p + n_dx, n_dw)):
        for triple in range(8):
            e = np.zeros(n)
            for j in range(triple, size, 8):
                e[base + j] = 0.5
            per_setting.append(e - t)
    r8 = solve(t, extra_ub=per_setting, extra_b=[0.0] * len(per_setting))

    solves = {"r1": r1, "r2": r2, "r3": r3, "r4": r4, "r6": r6, "r7": r7, "r8": r8}
    failed = {k: v.message for k, v in solves.items() if v.status != 0}
    if failed or r5.status not in (0, 2):  # 2 = infeasible, the expected Lemma 1 outcome
        raise RuntimeError(f"LP did not finish: {failed}, ns check: {r5.message}")

    def val(res):
        return float(res.fun)

    return LPResult(
        status=" | ".join(f"{k}: {v.message}" for k, v in {**solves, "r5": r5}.items()),
        min_bcd_on_x=val(r1),
        min_abc_on_w=val(r2),
        min_of_max=val(r3),
        min_of_sum=val(r4),
        min_bcd_on_x_if_abc_silent=val(r6),
        min_abc_on_w_if_bcd_silent=val(r7),
        ns_feasible=r5.status == 0,
        min_max_per_setting=val(r8),
    )


# ----------------------------------------------------------------------------- timing


@dataclass(frozen=True)
class Layout:
    """1-D layout in the privileged frame: positions (units of L) and measurement times (L / c)."""

    name: str
    x: dict
    t: dict


def fig3_layout(r: float) -> Layout:
    """Bancal et al. Fig. 3 caption: dB = (1 + 1/r)/4 + 1/(1 + r), dC = 3(1 + 1/r)/4 - 1/(1 + r),
    dD = 1; tA = 0, tB = tC = 2 / (1 + r), tD = 1 / r (c = 1)."""
    if r <= 1:
        raise ValueError("need v > c")
    db = 0.25 * (1 + 1 / r) + 1 / (1 + r)
    dc = 0.75 * (1 + 1 / r) - 1 / (1 + r)
    return Layout(
        "Bancal Fig. 3",
        {"A": 0.0, "B": db, "C": dc, "D": 1.0},
        {"A": 0.0, "B": 2 / (1 + r), "C": 2 / (1 + r), "D": 1 / r},
    )


def earliest_bc_layout(name: str, xb: float, xc: float, r: float) -> Layout:
    """A at 0 (t = 0), D at 1 (t = 1/r); B and C measured together at the earliest time that puts
    both inside the future v-cones of A and D (the ordering A < D < (B, C) of Fig. 3)."""
    t = max(max(x / r, 1 / r + abs(1 - x) / r) for x in (xb, xc)) * (1 + 1e-9)
    x = {"A": 0.0, "B": xb, "C": xc, "D": 1.0}
    return Layout(name, x, {"A": 0, "B": t, "C": t, "D": 1 / r})


def central_layout(r: float, half_gap: float = 1e-3) -> Layout:
    """B and C just either side of the midpoint, measured as soon as A's and D's v-cones allow."""
    return earliest_bc_layout("B, C at the midpoint", 0.5 - half_gap, 0.5 + half_gap, r)


def clustered_layout(r: float, offset: float = 1e-3) -> Layout:
    """B and C next to the receiver D (needed for low latency A -> D')."""
    return earliest_bc_layout("B, C next to D", 1.0 - offset, 1.0 - 2 * offset, r)


def in_future_vcone(e1, e2, lay: Layout, speed: float) -> bool:
    """True if event e2 lies in the closed future speed-cone of e1."""
    return lay.t[e2] - lay.t[e1] >= abs(lay.x[e2] - lay.x[e1]) / speed - 1e-15


def conditions(lay: Layout, r: float) -> dict:
    """App. E conditions i)-ii) and the v-cone ordering A < D < (B, C) of Fig. 3 (B, C in the future
    v-cones of A and D), plus B ~ C (outside each other's v-cones)."""
    t, x = lay.t, lay.x
    return {
        "i_times": max(t["A"], t["D"]) < min(t["B"], t["C"]),
        "A<D": in_future_vcone("A", "D", lay, r),
        "A,D<B": in_future_vcone("A", "B", lay, r) and in_future_vcone("D", "B", lay, r),
        "A,D<C": in_future_vcone("A", "C", lay, r) and in_future_vcone("D", "C", lay, r),
        "ii_B~C": abs(t["B"] - t["C"]) < abs(x["B"] - x["C"]) / r,
    }


def earliest_collection(lay: Layout, parties, site: float) -> float:
    """Earliest time at position `site` when the light-speed records of `parties` are all there."""
    return max(lay.t[p] + abs(site - lay.x[p]) for p in parties)


def channel_timing(lay: Layout) -> dict:
    """A -> D' (BCD records) decoded at D's site; D -> A' (ABC records) decoded at A's site. Latency
    is counted from the sender's measurement (A at t_A; D at t_D). 'superluminal' means the
    collection event lies outside the sender's future light cone (App. E condition iii)."""
    out = {}
    for name, sender, rec, recv_site in (
        ("A_to_Dprime", "A", ("B", "C", "D"), lay.x["D"]),
        ("D_to_Aprime", "D", ("A", "B", "C"), lay.x["A"]),
    ):
        t_ready = earliest_collection(lay, rec, recv_site)
        light = lay.t[sender] + abs(recv_site - lay.x[sender])
        out[name] = {
            "records_complete_at": t_ready,
            "latency_from_sender": t_ready - lay.t[sender],
            "light_time": light - lay.t[sender],
            "superluminal": t_ready < light,
            "effective_speed_over_c": abs(recv_site - lay.x[sender]) / (t_ready - lay.t[sender]),
        }
    out["both_superluminal_iii"] = (
        out["A_to_Dprime"]["superluminal"] and out["D_to_Aprime"]["superluminal"]
    )
    return out


def best_guaranteed_latency(r: float, n: int = 201) -> dict:
    """Grid search over B, C positions in [0, 1] (measured at the earliest time the v-cone
    ordering allows, B ~ C) for the layout that minimises the WORSE of the two channel latencies
    while keeping both superluminal (the theorem does not say which direction signals)."""
    best = None
    xs = np.linspace(0.0, 1.0, n)
    for xb in xs:
        for xc in xs:
            if xc <= xb:
                continue
            lay = earliest_bc_layout("grid", xb, xc, r)
            cond = conditions(lay, r)
            if not all(cond.values()):
                continue
            ch = channel_timing(lay)
            if not ch["both_superluminal_iii"]:
                continue
            worst = max(
                ch["A_to_Dprime"]["latency_from_sender"], ch["D_to_Aprime"]["latency_from_sender"]
            )
            if best is None or worst < best[0]:
                best = (worst, xb, xc)
    out = {"analytic_lower_bound_L_over_c": 0.5 * (1 + 1 / r), "grid_n": n}  # section 3.1
    if best:
        out.update({"worst_latency_L_over_c": best[0], "xB": best[1], "xC": best[2]})
    return out


def trials_for_error(tv: float, eps: float) -> float:
    """Trials sufficient to decide one bit at error <= eps from i.i.d. trials whose laws differ by
    total variation tv, using the best event and a threshold (Hoeffding): n = 2 ln(1/eps) / tv^2."""
    return math.inf if tv <= 0 else 2.0 * math.log(1.0 / eps) / tv**2


LY_S = 365.25 * 86400.0  # light-seconds per light-year
C_M_S = 299_792_458.0
BETA_CMB = 369.82e3 / C_M_S  # Planck 2018 dipole / c (ASSUMPTION: candidate frame, as in round 1)
V_MIN_UNEXCLUDED = 5e6  # c; published bipartite bound in near-CMB frames (Cocciaro et al. 2018)


def r_trial_fraction(r: float, beta: float = BETA_CMB) -> float:
    """Fraction of randomly timed trials that land in R when the frame is unknown: the B ~ C window
    d / v over the simultaneity ambiguity beta d / c, doubled for either sign: f ~ 2 / (r beta)."""
    return min(1.0, 2.0 / (r * beta))


def bc_sync_window_s(separation_ly: float, r: float) -> float:
    """B ~ C needs |t_B - t_C| < |x_B - x_C| / v in the privileged frame: the sync window."""
    return separation_ly * LY_S / r


def frame_sync_ambiguity_s(separation_ly: float, beta: float = BETA_CMB) -> float:
    """Simultaneity offset between the lab and a frame moving at beta along the B-C line: about
    beta * d / c. B ~ C cannot be arranged unless the frame is known to better than the window."""
    return beta * separation_ly * LY_S


def stock_per_second(bits_per_s: float, trials_per_bit: float) -> float:
    """4-party entangled states consumed per second. Each trial uses one state, whose B, C and D
    shares had to travel to their stations at <= c before use (no faster replenishment exists)."""
    return bits_per_s * trials_per_bit


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--distance-ly", type=float, default=20.0)
    ap.add_argument("--r", type=float, nargs="+", default=[2.0, 10.0, 1e4, 6.3e8])
    a = ap.parse_args(argv)
    pq = quantum_distribution()
    lp = signalling_lp(pq)
    out = {
        "provenance": "derived (exact LP and kinematics) under the stated model assumptions",
        "positive_control_S_quantum": bell_value(pq),
        "S_local_bound": S_LOCAL_BOUND,
        "quantum_is_no_signalling": quantum_is_no_signalling(pq),
        "lp": lp.__dict__,
        "timing": {},
    }
    tv_guaranteed = lp.min_of_max  # at least one direction carries at least this (mean TV)
    out["statistics"] = {
        "tv_at_least_one_direction": tv_guaranteed,
        "trials_per_bit_eps_1e-3": trials_for_error(tv_guaranteed, 1e-3),
        "trials_per_bit_eps_1e-6": trials_for_error(tv_guaranteed, 1e-6),
        "states_per_s_for_8_kbit_per_s_voice_eps_1e-3": stock_per_second(
            8e3, trials_for_error(tv_guaranteed, 1e-3)
        ),
        "max_tv_best_setting_pair": lp.min_max_per_setting,
        "unknown_frame_dilution": {
            f"r={r:g}": {
                "R_trial_fraction": r_trial_fraction(r),
                "trials_per_bit_eps_1e-3": trials_for_error(
                    r_trial_fraction(r) * tv_guaranteed, 1e-3
                ),
            }
            for r in (V_MIN_UNEXCLUDED, 6.3e8)
        },
    }
    out["sync"] = {
        f"B-C {d} ly apart, r={r:g}": {
            "window_s": bc_sync_window_s(d, r),
            "frame_ambiguity_s_at_CMB_beta": frame_sync_ambiguity_s(d),
        }
        for d in (0.04, 1e3 / 9.4607e15)
        for r in (1e4, 6.3e8)
    }
    for r in a.r:
        rows = {}
        for lay in (fig3_layout(r), central_layout(r), clustered_layout(r)):
            ch = channel_timing(lay)
            rows[lay.name] = {
                "conditions": conditions(lay, r),
                "channels": ch,
                "A_to_Dprime_latency_years": ch["A_to_Dprime"]["latency_from_sender"]
                * a.distance_ly,
                "D_to_Aprime_latency_years": ch["D_to_Aprime"]["latency_from_sender"]
                * a.distance_ly,
                "positions_ly": {k: v * a.distance_ly for k, v in lay.x.items()},
            }
        g = best_guaranteed_latency(r)
        w = g.get("worst_latency_L_over_c")
        rows["best_guaranteed_grid"] = {
            **g,
            "worst_latency_years": None if w is None else w * a.distance_ly,
        }
        out["timing"][f"r={r:g}"] = rows
    print(json.dumps(out, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
