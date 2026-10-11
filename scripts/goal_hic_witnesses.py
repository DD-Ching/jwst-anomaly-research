"""Which finite-speed hidden-influence witness is the best experimental target?

Goal round 3 (docs/goal_round3/README.md). Compares, offline and under the same assumptions:
- Bancal et al. 2012 (Nature Phys. 8, 867; arXiv:1110.3795v2), Eq. 4: S <= 7, 4-qubit state of
  App. B (scripts/goal_bancal_channel.py);
- Li, Hu, Deng & Scarani 2026 (arXiv:2608.05271v1), Eq. 5: LC4 four-qubit linear cluster,
  S4 <= 6, and Eq. 6: LC5 five-qubit linear cluster, S5 <= 10.

All three are no-BC witnesses of the same set HIC_E = NS & CL(BC|E) (Li et al. Eq. 1; Bancal
App. C-D), with early side E = AD (AD E for LC5) and blind late pair B ~ C.

Computed here [derived]:
1. Positive controls: quantum values 7.2014, 4 + 2 sqrt 2, 6 + 4 sqrt 2 and the white-noise
   thresholds 0.9720, 0.8787, 0.8579 quoted by Li et al.
2. Thresholds of each witness under white noise, local dephasing, local depolarising noise,
   symmetric readout flips, and detection loss with no-click events assigned to a fixed outcome
   (no postselection; the assignment is a local relabelling, so the bound stays valid).
3. The exact linear programme over full distributions: the lowest white-noise visibility at
   which the observed marginals still force signalling (the best threshold any no-BC witness can
   reach for that state and those measurements), and the smallest forced signal (total variation)
   when the marginals are reproduced.
4. Trials for a given significance from the exact per-trial variance of the Experiment-v1
   estimator (Bancal App. E: ABD records from B < C trials, ACD records from C < B trials).
5. The worst-case fraction of randomly timed trials with B ~ C in an unknown frame
   (goal_bancal_channel.r_fraction_table).

Outcome bit 0 <-> eigenvalue +1. Parties are qubits in the order given by each Witness.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

_here = Path(__file__).resolve().parent
if "goal_bancal_channel" in sys.modules:
    bancal = sys.modules["goal_bancal_channel"]
else:
    _spec = importlib.util.spec_from_file_location(
        "goal_bancal_channel", _here / "goal_bancal_channel.py"
    )
    bancal = importlib.util.module_from_spec(_spec)
    sys.modules["goal_bancal_channel"] = bancal
    _spec.loader.exec_module(bancal)

I2 = np.eye(2)
SX, SZ = bancal.SX, bancal.SZ
SY = np.array([[0.0, -1j], [1j, 0.0]])


# ----------------------------------------------------------------------------- witnesses


@dataclass
class Witness:
    name: str
    parties: str  # one letter per qubit, in qubit order
    psi: np.ndarray
    obs: dict  # party -> (observable for setting 0, for setting 1)
    terms: list  # (coefficient, ((party, setting), ...))
    bound: float
    early: str
    source: str
    quoted: dict = field(default_factory=dict)  # values quoted in the source (positive control)

    @property
    def n(self) -> int:
        return len(self.parties)

    def index(self, party: str) -> int:
        return self.parties.index(party)


def _cz(n: int, i: int, j: int) -> np.ndarray:
    d = 2**n
    diag = np.ones(d)
    for k in range(d):
        if (k >> (n - 1 - i)) & 1 and (k >> (n - 1 - j)) & 1:
            diag[k] = -1.0
    return diag


def linear_cluster(n: int) -> np.ndarray:
    """CZ_{12} CZ_{23} ... |+>^n (Li et al. LC4, LC5 definitions)."""
    psi = np.ones(2**n) / math.sqrt(2**n)
    for i in range(n - 1):
        psi = psi * _cz(n, i, i + 1)
    return psi


B_ROT = ((SZ + SX) / math.sqrt(2.0), (SZ - SX) / math.sqrt(2.0))  # Li et al. Eq. S23


def _terms(spec: str) -> list:
    """'+1 A0 B0; -2 C1 D1' -> [(1, (('A', 0), ('B', 0))), (-2, ...)]."""
    out = []
    for part in spec.split(";"):
        tok = part.split()
        out.append((float(tok[0]), tuple((t[0], int(t[1])) for t in tok[1:])))
    return out


def bancal_witness() -> Witness:
    return Witness(
        name="Bancal2012",
        parties="ABCD",
        psi=bancal.quantum_state(),
        obs=bancal.OBS,
        terms=[(float(c), p) for c, p in bancal.S_TERMS],
        bound=bancal.S_LOCAL_BOUND,
        early="AD",
        source="Bancal et al. 2012, Eq. 4 and App. B (arXiv:1110.3795v2)",
        quoted={"quantum": 7.2014, "white_threshold": 0.9720},
    )


def lc4_witness() -> Witness:
    return Witness(
        name="LC4",
        parties="ABCD",
        psi=linear_cluster(4),
        obs={"A": (SX, SZ), "B": B_ROT, "C": (SZ, SX), "D": (SX, SZ)},
        terms=_terms("1 A0 B0; 1 A0 B1; 1 A1 B0 D0; -1 A1 B1 D0; 2 C0 D0; 2 A0 C1 D1"),
        bound=6.0,
        early="AD",
        source="Li et al. 2026, Eq. 5 (arXiv:2608.05271v1)",
        quoted={"quantum": 4 + 2 * math.sqrt(2), "white_threshold": 0.8787},
    )


def lc5_witness() -> Witness:
    return Witness(
        name="LC5",
        parties="ABCDE",
        psi=linear_cluster(5),
        obs={"A": (SX, SZ), "B": B_ROT, "C": (SZ, SX), "D": (SX, SZ), "E": (SZ, SX)},
        terms=_terms(
            "1 A0 B0; 1 A0 B1; -1 B0 D1; -1 B1 D1; 1 B0 E1; 1 B1 E1; 2 A0 C1 D1; 2 A0 C1 E1;"
            " -2 D1 E1; 1 A0 B0 D1 E1; 1 A0 B1 D1 E1; 2 A1 B0 D0 E0; -2 A1 B1 D0 E0; 4 C0 D0 E0"
        ),
        bound=10.0,
        early="ADE",
        source="Li et al. 2026, Eq. 6 (arXiv:2608.05271v1)",
        quoted={"quantum": 6 + 4 * math.sqrt(2), "white_threshold": 0.8579},
    )


WITNESSES = (bancal_witness, lc4_witness, lc5_witness)


# ----------------------------------------------------------------------------- noise


def _kron(ops) -> np.ndarray:
    out = np.array([[1.0]])
    for op in ops:
        out = np.kron(out, op)
    return out


def density(w: Witness) -> np.ndarray:
    return np.outer(w.psi, w.psi.conj())


def white_distribution(p_ideal: np.ndarray, vis: float) -> np.ndarray:
    """P under white noise: affine in vis, uniform outcomes 1 / 2^n for every setting."""
    n = p_ideal.ndim // 2
    return vis * p_ideal + (1 - vis) / 2**n


def white(rho: np.ndarray, vis: float) -> np.ndarray:
    d = rho.shape[0]
    return vis * rho + (1 - vis) * np.eye(d) / d


def local_channel(rho: np.ndarray, kraus: list) -> np.ndarray:
    """Apply the same single-qubit channel (Kraus operators) to every qubit."""
    n = int(round(math.log2(rho.shape[0])))
    for q in range(n):
        new = np.zeros_like(rho, dtype=complex)
        for k in kraus:
            op = _kron([k if i == q else I2 for i in range(n)])
            new += op @ rho @ op.conj().T
        rho = new
    return rho


def dephasing(rho: np.ndarray, p: float) -> np.ndarray:
    """Z flip with probability p on every qubit (pure T2 dephasing)."""
    return local_channel(rho, [math.sqrt(1 - p) * I2, math.sqrt(p) * SZ])


def depolarising(rho: np.ndarray, q: float) -> np.ndarray:
    """rho -> (1 - q) rho + q I/2 on every qubit."""
    k = [math.sqrt(1 - 3 * q / 4) * I2] + [math.sqrt(q / 4) * s for s in (SX, SY, SZ)]
    return local_channel(rho, k)


def amplitude_damping(rho: np.ndarray, g: float) -> np.ndarray:
    """T1 decay with probability g on every qubit (|1> -> |0>)."""
    k0 = np.array([[1.0, 0.0], [0.0, math.sqrt(1 - g)]])
    k1 = np.array([[0.0, math.sqrt(g)], [0.0, 0.0]])
    return local_channel(rho, [k0, k1])


# ----------------------------------------------------------------------------- correlators


def expectation(rho: np.ndarray, w: Witness, factors: dict) -> float:
    """<prod_k F_k> with F_k given per party (identity for absent parties)."""
    op = _kron([factors.get(p, I2) for p in w.parties])
    return float(np.real(np.trace(rho @ op)))


def witness_value(rho: np.ndarray, w: Witness, readout=None) -> float:
    """S with each party's observable O replaced by a O + b I (readout model; see
    `readout_affine`). readout: dict (party, setting) -> (a, b); None = ideal."""
    total = 0.0
    for coef, parts in w.terms:
        factors = {}
        for p, s in parts:
            a, b = (1.0, 0.0) if readout is None else readout[(p, s)]
            factors[p] = a * w.obs[p][s] + b * I2
        total += coef * expectation(rho, w, factors)
    return total


def readout_affine(w: Witness, flip: float = 0.0, eta: float = 1.0, assign=None) -> dict:
    """Per (party, setting): symmetric flip probability `flip`, then detection efficiency `eta`
    with no-click recorded as the outcome assign[(party, setting)] in {+1, -1} (default +1).
    Observed +-1 variable = a O + b with a = eta (1 - 2 flip), b = (1 - eta) * assigned value."""
    out = {}
    for p in w.parties:
        for s in (0, 1):
            v = 1 if assign is None else assign[(p, s)]
            out[(p, s)] = (eta * (1 - 2 * flip), (1 - eta) * v)
    return out


def distribution(rho: np.ndarray, w: Witness, readout=None) -> np.ndarray:
    """P[o_1..o_n, s_1..s_n] for rho, the witness observables and the readout model."""
    n = w.n
    proj = {}
    for p in w.parties:
        for s in (0, 1):
            a, b = (1.0, 0.0) if readout is None else readout[(p, s)]
            o = a * w.obs[p][s] + b * I2  # observed +-1 variable as an operator (POVM mean)
            proj[(p, s)] = (0.5 * (I2 + o), 0.5 * (I2 - o))
    out = np.zeros((2,) * (2 * n))
    for st in itertools.product((0, 1), repeat=n):
        for oc in itertools.product((0, 1), repeat=n):
            op = _kron([proj[(p, st[i])][oc[i]] for i, p in enumerate(w.parties)])
            out[oc + st] = float(np.real(np.trace(rho @ op)))
    return out


# ----------------------------------------------------------------------------- thresholds


def bisect(fun, lo: float, hi: float, tol: float = 1e-7) -> float | None:
    """Root of fun on [lo, hi] (fun(lo) and fun(hi) of opposite sign), else None."""
    flo, fhi = fun(lo), fun(hi)
    if flo * fhi > 0:
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = fun(mid)
        if (fm > 0) == (fhi > 0):
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def noise_thresholds(w: Witness) -> dict:
    """Largest noise at which the witness is still violated, per noise model (others ideal)."""
    rho = density(w)

    def gap(r, readout=None):
        return witness_value(r, w, readout) - w.bound

    out = {
        "white_visibility_min": bisect(lambda v: gap(white(rho, v)), 0.0, 1.0),
        "dephasing_p_max": bisect(lambda p: gap(dephasing(rho, p)), 0.0, 0.5),
        "depolarising_q_max": bisect(lambda q: gap(depolarising(rho, q)), 0.0, 1.0),
        "amplitude_damping_g_max": bisect(lambda g: gap(amplitude_damping(rho, g)), 0.0, 1.0),
        "readout_flip_max": bisect(lambda f: gap(rho, readout_affine(w, flip=f)), 0.0, 0.5),
    }
    out["detection_efficiency_min"] = detection_threshold(w, rho)
    # combined: a realistic state visibility plus readout flips or detection loss
    for v in (0.98, 0.95):
        rv = white(rho, v)
        out[f"at_white_{v}"] = {
            "readout_flip_max": bisect(lambda f, r=rv: gap(r, readout_affine(w, flip=f)), 0.0, 0.5),
            "detection_efficiency_min": detection_threshold(w, rv)["eta_min"],
        }
    return out


def detection_threshold(w: Witness, rho: np.ndarray, grid: int = 2001) -> dict:
    """Lowest per-detector efficiency eta such that some no-click assignment (each party and
    setting records +1 or -1 on no click; no postselection) violates the bound at eta and at every
    efficiency above it. Each assignment's S(eta) is a polynomial (degree <= max body); it is
    evaluated exactly on a grid over [0, 1] for all 4^n assignments at once, and the crossing below
    the violated tail is refined by bisection."""
    keys = [(p, s) for p in w.parties for s in (0, 1)]
    signs = np.array(list(itertools.product((1, -1), repeat=len(keys))), dtype=float)
    col = {k: i for i, k in enumerate(keys)}
    etas = np.linspace(0.0, 1.0, grid)

    corr = {}
    for _, parts in w.terms:
        for r in range(len(parts) + 1):
            for sub in itertools.combinations(parts, r):
                if sub not in corr:
                    corr[sub] = expectation(rho, w, {p: w.obs[p][s] for p, s in sub})

    def values(eta: np.ndarray, sg: np.ndarray) -> np.ndarray:  # (assignments, len(eta))
        total = np.zeros((len(sg), len(eta)))
        for coef, parts in w.terms:
            for r in range(len(parts) + 1):
                for sub in itertools.combinations(parts, r):
                    rest = [x for x in parts if x not in sub]
                    sgn = np.prod(sg[:, [col[x] for x in rest]], axis=1) if rest else 1.0
                    weight = eta**r * (1 - eta) ** len(rest)
                    total += coef * corr[sub] * np.outer(np.ones(len(sg)) * sgn, weight)
        return total

    viol = values(etas, signs) > w.bound
    tail = np.flip(np.logical_and.accumulate(np.flip(viol, axis=1), axis=1), axis=1)
    if not tail[:, -1].any():
        return {"eta_min": None}
    first = np.where(tail.any(axis=1), tail.argmax(axis=1), grid)
    k = int(first.argmin())
    i0 = int(first[k])
    if i0 == 0:
        eta = 0.0
    else:
        row = signs[k : k + 1]

        def fk(e):
            return float(values(np.array([e]), row)[0, 0] - w.bound)

        eta = bisect(fk, float(etas[i0 - 1]), float(etas[i0]), tol=1e-7)
    assign = {f"{p}{s}": int(v) for (p, s), v in zip(keys, signs[k], strict=True)}
    return {"eta_min": eta, "assignment": assign}


# ----------------------------------------------------------------------------- exact LP


def _flat(n: int, oc, st) -> int:
    return int(np.ravel_multi_index(tuple(oc) + tuple(st), (2,) * (2 * n)))


def hic_lp(p_obs: np.ndarray, w: Witness, objective: str = "minmax") -> dict:
    """LP over full P(o|s) >= 0 with: the ABD.. and ACD.. marginals (all parties but C, all but
    B) equal to p_obs for every setting (Bancal App. D; Li et al. Lemma 1); BC local conditioned
    on every early (o_E, s_E) (eight CHSH facets in joint form; App. C, Li Eq. 3).
    objective 'feasible_ns': can such a P be fully no-signalling? 'minmax': the smallest largest
    signal, over early senders X, of the mean total variation of the other parties' marginal
    between X's two settings."""
    n = w.n
    ib, ic = w.index("B"), w.index("C")
    early = [w.index(x) for x in w.early]
    nv = 2 ** (2 * n)
    eq_r, eq_c, eq_v, eq_b = [], [], [], []
    row = 0
    for drop in (ic, ib):
        marg = p_obs.sum(axis=drop)
        for st in itertools.product((0, 1), repeat=n):
            for ok in itertools.product((0, 1), repeat=n - 1):
                for od in (0, 1):
                    oc = list(ok)
                    oc.insert(drop, od)
                    eq_r.append(row)
                    eq_c.append(_flat(n, oc, st))
                    eq_v.append(1.0)
                # marg keeps the n - 1 outcome axes, then all n setting axes
                eq_b.append(float(marg[tuple(ok) + tuple(st)]))
                row += 1
    ub_r, ub_c, ub_v, ub_b = [], [], [], []
    urow = 0
    signs = [s for s in itertools.product((1, -1), repeat=4) if s.count(-1) % 2 == 1]
    p_early = p_obs.sum(axis=(ib, ic))  # outcome axes: others; settings: all n
    others = [k for k in range(n) if k not in (ib, ic)]
    for oe in itertools.product((0, 1), repeat=len(others)):
        for se in itertools.product((0, 1), repeat=len(others)):
            st0 = [0] * n
            for k, s in zip(others, se, strict=True):
                st0[k] = s
            pe = float(p_early[tuple(oe) + tuple(st0)])
            for sg in signs:
                for k_yz, (y, z) in enumerate(itertools.product((0, 1), repeat=2)):
                    st = list(st0)
                    st[ib], st[ic] = y, z
                    for b, c in itertools.product((0, 1), repeat=2):
                        oc = [0] * n
                        for k, o in zip(others, oe, strict=True):
                            oc[k] = o
                        oc[ib], oc[ic] = b, c
                        ub_r.append(urow)
                        ub_c.append(_flat(n, oc, st))
                        ub_v.append(sg[k_yz] * (-1) ** (b + c))
                ub_b.append(2.0 * pe)
                urow += 1
    # signal variables: for each early sender X, |marginal(others | s, x=0) - (x=1)| <= u
    n_sig = 0
    sig_groups = []
    for x in early:
        rest = [k for k in range(n) if k != x]
        cols = []
        for orest in itertools.product((0, 1), repeat=n - 1):
            for srest in itertools.product((0, 1), repeat=n - 1):
                u = nv + n_sig
                n_sig += 1
                cols.append(u)
                for sgn in (1, -1):
                    for ox in (0, 1):
                        for sx, val in ((0, sgn), (1, -sgn)):
                            oc = [0] * n
                            st = [0] * n
                            for k, o, s in zip(rest, orest, srest, strict=True):
                                oc[k], st[k] = o, s
                            oc[x], st[x] = ox, sx
                            ub_r.append(urow)
                            ub_c.append(_flat(n, oc, st))
                            ub_v.append(float(val))
                    ub_r.append(urow)
                    ub_c.append(u)
                    ub_v.append(-1.0)
                    ub_b.append(0.0)
                    urow += 1
        sig_groups.append((w.parties[x], cols, 2 ** (n - 1)))
    t_col = nv + n_sig
    ncol = t_col + 1
    a_eq = sparse.csr_matrix((eq_v, (eq_r, eq_c)), shape=(row, ncol))
    b_eq = np.array(eq_b)

    def tv_row(cols, n_settings):
        r = np.zeros(ncol)
        r[cols] = 0.5 / n_settings
        return r

    base_ub = sparse.csr_matrix((ub_v, (ub_r, ub_c)), shape=(urow, ncol))
    bounds = [(0, None)] * ncol
    if objective == "feasible_ns":
        c = np.zeros(ncol)
        extra = sparse.csr_matrix(
            (np.ones(n_sig), (np.arange(n_sig), nv + np.arange(n_sig))), shape=(n_sig, ncol)
        )
        res = linprog(
            c,
            A_ub=base_ub,
            b_ub=np.array(ub_b),
            A_eq=sparse.vstack([a_eq, extra]),
            b_eq=np.concatenate([b_eq, np.zeros(n_sig)]),
            bounds=bounds,
            method="highs",
        )
        if res.status not in (0, 2):
            raise RuntimeError(f"LP did not finish: {res.message}")
        return {"ns_feasible": res.status == 0}
    rows = [tv_row(cols, ns) for _, cols, ns in sig_groups]
    tmat = np.zeros((len(rows), ncol))
    tmat[:, t_col] = -1.0
    a_ub = sparse.vstack([base_ub, sparse.csr_matrix(np.array(rows) + tmat)])
    c = np.zeros(ncol)
    c[t_col] = 1.0
    res = linprog(
        c,
        A_ub=a_ub,
        b_ub=np.concatenate([ub_b, np.zeros(len(rows))]),
        A_eq=a_eq,
        b_eq=b_eq,
        bounds=bounds,
        method="highs",
    )
    if res.status != 0:
        raise RuntimeError(f"LP did not finish: {res.message}")
    per = {
        name: float(np.array(r) @ res.x) for (name, _, _), r in zip(sig_groups, rows, strict=True)
    }
    return {"min_max_signal": float(res.fun), "signals_at_optimum": per}


def lp_white_threshold(w: Witness, tol: float = 1e-4) -> float:
    """Lowest white-noise visibility at which the marginals cannot come from a no-signalling
    HIC distribution: the best threshold any no-BC witness reaches for this state and setting."""
    rho = density(w)

    p_ideal = distribution(rho, w)

    def forced(v):
        ns = hic_lp(white_distribution(p_ideal, v), w, "feasible_ns")["ns_feasible"]
        return -1.0 if ns else 1.0

    lo, hi = 0.5, 1.0
    if forced(hi) < 0:
        return None
    if forced(lo) > 0:
        raise ValueError("marginals still force signalling at visibility 0.5; widen the search")
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if forced(mid) > 0:
            hi = mid
        else:
            lo = mid
    return hi


# ----------------------------------------------------------------------------- statistics


def v1_estimator_moments(p_obs: np.ndarray, w: Witness) -> dict:
    """Experiment v1 (Bancal App. E) with uniformly random settings: each trial is a B-family
    trial (ABD.. recorded, B before C) or a C-family trial with probability 1/2. Terms with B use
    B-family trials, terms with C use C-family trials, terms with neither use all trials. The
    per-trial estimator sum_j w_j 1[term j observable] prod(outcomes) / pi_j is unbiased; its
    exact mean, variance and range follow from p_obs (i.i.d. trials ASSUMED)."""
    n = w.n
    info = []
    for coef, parts in w.terms:
        names = {p for p, _ in parts}
        if {"B", "C"} <= names:
            raise ValueError("Experiment v1 cannot estimate a term containing both B and C")
        fam = "B" if "B" in names else ("C" if "C" in names else None)
        pi = (0.5 if fam else 1.0) * 0.5 ** len(parts)
        info.append((coef, parts, fam, pi))
    m1 = m2 = 0.0
    lo, hi = math.inf, -math.inf
    for fam_trial in ("B", "C"):
        for st in itertools.product((0, 1), repeat=n):
            for oc in itertools.product((0, 1), repeat=n):
                prob = 0.5 * 0.5**n * p_obs[oc + st]
                val = 0.0
                for coef, parts, fam, pi in info:
                    if fam not in (None, fam_trial):
                        continue
                    if any(st[w.index(p)] != s for p, s in parts):
                        continue
                    sign = (-1) ** sum(oc[w.index(p)] for p, _ in parts)
                    val += coef * sign / pi
                lo, hi = min(lo, val), max(hi, val)
                m1 += prob * val
                m2 += prob * val * val
    return {"mean": m1, "var": m2 - m1 * m1, "range": hi - lo}


def trials_needed(mom: dict, bound: float, z: float = 5.0, p_value: float = 1e-6) -> dict:
    """Trials for a z-sigma excess over the bound (Gaussian, i.i.d.) and for a p-value with a
    bounded-increment (Azuma-Hoeffding) test that also holds with memory between trials, given
    the per-trial range. Both are sufficient for these specific tests, not minimal."""
    gap = mom["mean"] - bound
    if gap <= 0:
        return {"gap": gap, "gaussian_z": math.inf, "azuma": math.inf}
    return {
        "gap": gap,
        "gaussian_z": z * z * mom["var"] / gap**2,
        "azuma": mom["range"] ** 2 * math.log(1 / p_value) / (2 * gap**2),
    }


# ----------------------------------------------------------------------------- main


def analyse(w: Witness, with_lp: bool = True) -> dict:
    rho = density(w)
    p_ideal = distribution(rho, w)
    out = {
        "source": w.source,
        "parties": w.parties,
        "early": w.early,
        "bound": w.bound,
        "quantum_value": (sq := witness_value(rho, w)),
        "quoted": w.quoted,
        "terms": len(w.terms),
        "max_body": max(len(p) for _, p in w.terms),
        "relative_gap": sq / w.bound - 1,
        "thresholds": noise_thresholds(w),
    }
    mom = v1_estimator_moments(p_ideal, w)
    out["quantum_value_from_distribution"] = mom["mean"]  # Experiment-v1 estimator mean
    stats = {"ideal": {**mom, **trials_needed(mom, w.bound)}}
    for v in (0.98, 0.95, 0.92):
        m = v1_estimator_moments(white_distribution(p_ideal, v), w)
        stats[f"white_{v}"] = {**m, **trials_needed(m, w.bound)}
    out["v1_statistics"] = stats
    if with_lp:
        out["lp_white_threshold"] = lp_white_threshold(w)
        sig = {}
        for v in (1.0, 0.99, 0.97, 0.95, 0.9):
            sig[f"white_{v}"] = hic_lp(white_distribution(p_ideal, v), w)
        sig["dephasing_0.02"] = hic_lp(distribution(dephasing(rho, 0.02), w), w)
        sig["readout_flip_0.01"] = hic_lp(distribution(rho, w, readout_affine(w, flip=0.01)), w)
        out["forced_signal"] = sig
    return out


def _finite(x):
    """JSON has no inf / nan: write them as null."""
    if isinstance(x, dict):
        return {k: _finite(v) for k, v in x.items()}
    if isinstance(x, list | tuple):
        return [_finite(v) for v in x]
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-lp5", action="store_true", help="skip the 5-party LPs (slow)")
    a = ap.parse_args(argv)
    out = {
        "provenance": "derived (exact quantum predictions, exact LP) under the stated model "
        "assumptions; no data",
        "witnesses": {},
    }
    for make in WITNESSES:
        w = make()
        out["witnesses"][w.name] = analyse(w, with_lp=not (a.no_lp5 and w.n == 5))
    out["r_fraction_worst_case"] = bancal.r_fraction_table()
    print(json.dumps(_finite(out), indent=1, default=str, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
