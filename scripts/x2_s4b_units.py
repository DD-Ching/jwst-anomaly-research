"""X2 S4b: physical strain and energy flux of the S4 echo readouts (one cell, GW150914 dt = 0.45 s).

Derivation and pre-written checks: docs/state/units/X2_s4_matter.md, section S4b. H2 (Yang & Lan
Eq. 34) is the Moncrief / Regge-Wheeler function: strain h ~ (1/r) int H2 dt far from the source,
GW energy flux ~ H2^2, and the physical EM flux has the same weight in the symmetrised form.

Runs (NED magnetic + phantom, dx = 0.1, observers at r* = +1500 and -1500, five round trips):
  A  input H2 = p   (S4's ringdown packet, unchanged initial data)
  B  input H2 = g   (smooth Gaussian-windowed packet)
  C  input H2 = g'  (the same strain g; checks h = int H2 dt with no detrending)
  D  input H2 = p with the coupling removed (V22 alone; diagnostic, not a self-consistent model)

  python scripts/x2_s4b_units.py [--from-npz]   -> results/x2/s4b_units.json
  (about 15 minutes on 4 processes; waveforms go to $X2_SCRATCH, never to the repository)
"""

from __future__ import annotations

import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import x2_axial_ned as an  # noqa: E402
import x2_echo_delay as xe  # noqa: E402

OUT = ROOT / "results" / "x2"
W0, TAU = 0.5307, 12.2  # S3 / S4 ringdown packet
GC, GS = 40.0, 8.0  # smooth packet centre and width (M)
X_OBS = 1500.0
DX = 0.1
N_TRIPS = 5


def cell_a():
    return 2 + xe.a_for_delay(0.45 / (an.M_DET * an.T_SUN))


def packet(kind, u):
    if kind == "p":
        return np.where(u > 0, np.exp(-u / TAU) * np.sin(W0 * u), 0.0)
    s = u - GC
    env = np.exp(-(s**2) / (2 * GS**2))
    if kind == "g":
        return env * np.sin(W0 * s)
    if kind == "gdot":  # d/du of g
        return env * (W0 * np.cos(W0 * s) - s / GS**2 * np.sin(W0 * s))
    raise ValueError(kind)


def run(job):
    """Inward packet in the GW channel, 25 M inside the GW peak (as S4's echo_run); waveforms."""
    name, model, kind = job
    a = cell_a()
    x_pk = an.gw_peak_x("af", a)  # the GW-channel diagonal is V_AF in both models
    x0 = x_pk - 25.0

    def init(x, dt):
        u0 = x0 - x
        return [(packet(kind, u0), packet(kind, u0 - dt))]

    t_end = N_TRIPS * 1.04 * 4 * x_pk + X_OBS + 400.0
    t, rec = an.evolve(model, a, 6000.0, 7000.0, DX, t_end, init, [X_OBS, -X_OBS])
    return name, {"t": t, "rec": rec, "x_pk": x_pk, "x0": x0, "a": a}


def windows(res, n=2):
    """S4's pulse windows at the + observer: t_ref + k 4L +- 0.45 4L, t_ref = free flight of the
    ringdown packet's energy centroid (the same windows are used for every run)."""
    x_pk, x0 = res["x_pk"], res["x0"]
    uu = np.arange(0, 40 * TAU, 0.01)
    pe = packet("p", uu) ** 2
    t_ref = X_OBS - x0 + float(np.trapezoid(uu * pe, uu) / np.trapezoid(pe, uu))
    dte = 4 * x_pk
    return [(t_ref + k * dte - 0.45 * dte, t_ref + k * dte + 0.45 * dte) for k in range(1, n + 1)]


def energy(y, t):
    return float(np.trapezoid(y * y, t))


def _budget(rec, t, e0, weight):
    dt = t[1] - t[0]
    f = (lambda y: y) if weight == "psi2" else (lambda y: np.gradient(y, dt))
    names = ["plus_GW", "plus_EM", "minus_GW", "minus_EM"]
    vals = {}
    for o, side in ((0, "plus"), (1, "minus")):
        for c, ch in enumerate(("GW", "EM")[: rec.shape[1]]):
            vals[f"{side}_{ch}"] = energy(f(rec[o, c]), t) / e0
    out = {k: vals.get(k, 0.0) for k in names}
    out["EM_total"] = out["plus_EM"] + out["minus_EM"]
    out["residual_trapped_or_not_arrived"] = 1 - sum(out[k] for k in names)
    return out


def analyse(res):
    t = res["A"]["t"]
    dt = t[1] - t[0]
    win = windows(res["A"])
    uu = np.arange(0, 60 * TAU, dt)
    out = {
        "label": "derived; physical strain / flux of S4 readouts; one cell; NED + phantom, axial",
        "design": "docs/state/units/X2_s4_matter.md, section S4b (written before computing)",
        "a_minus_2M": res["A"]["a"] - 2,
        "dx": DX,
        "t_end_M": float(t[-1]),
        "readings": {
            "strain_source": "source strain = packet p; physical H2 = -d/dt(run output); strain = "
            "run output; physical flux ~ (d/dt run output)^2  [P2; used]",
            "moncrief_source": "H2 itself = packet p; its strain int p dt carries a DC offset "
            "(memory); physical flux ~ (run output)^2  [reported for comparison only]",
        },
    }

    # 1. strain check. Input g'(u), u = x0 - x - t, is -d/dt of input g, so int H2_C dt = -H2_B.
    hb = res["B"]["rec"][0, 0]
    hc = np.cumsum(res["C"]["rec"][0, 0]) * dt  # from t = 0; the output is zero before arrival
    chk = []
    for lo, hi in win:
        m = (t > lo) & (t < hi)
        rb, rc = math.sqrt(energy(hb[m], t[m])), math.sqrt(energy(hc[m], t[m]))
        chk.append(
            {
                "hrss_B": rb,
                "hrss_int_C": rc,
                "rel_diff": rc / rb - 1,
                "rel_L2_as_written": math.sqrt(energy(hc[m] - hb[m], t[m])) / rb,
                "rel_L2_sign_corrected": math.sqrt(energy(hc[m] + hb[m], t[m])) / rb,
            }
        )
    end_val = float(abs(hc[-1]) / np.max(np.abs(hc)))
    out["strain_check"] = {
        "pulses": chk,
        "end_value_over_peak": end_val,
        "pass_as_written": bool(
            all(abs(c["rel_diff"]) <= 0.01 and c["rel_L2_as_written"] <= 0.02 for c in chk)
            and end_val <= 0.01
        ),
        "pass_sign_corrected_post_hoc": bool(
            all(abs(c["rel_diff"]) <= 0.01 and c["rel_L2_sign_corrected"] <= 0.02 for c in chk)
            and end_val <= 0.01
        ),
    }

    # 2. per pulse at the + observer: strain h_rss ratio and receiver-side EM energy fractions
    ref_hrss = math.sqrt(energy(packet("p", uu), uu))

    def per_pulse(key):
        rec = res[key]["rec"]
        rows = []
        for lo, hi in win:
            m = (t > lo) & (t < hi)
            g = rec[0, 0, m]
            row = {"hrss_ratio_strain": math.sqrt(energy(g, t[m])) / ref_hrss}
            if rec.shape[1] > 1:
                e = rec[0, 1, m]
                gd, ed = np.gradient(g, dt), np.gradient(e, dt)
                row["em_fraction_receiver_strain_source"] = energy(ed, t[m]) / (
                    energy(gd, t[m]) + energy(ed, t[m])
                )
                row["em_fraction_receiver_moncrief_source"] = energy(e, t[m]) / (
                    energy(g, t[m]) + energy(e, t[m])
                )
            rows.append(row)
        return rows

    # 3. budgets relative to the initial packet energy, whole run, both observers
    def budgets(key, kind):
        pk = packet(kind, uu)
        rec = res[key]["rec"]
        e0s = energy(np.gradient(pk, dt), uu)
        e0m = energy(pk, uu)
        return {
            "strain_source": _budget(rec, t, e0s, "psidot2"),
            "moncrief_source": _budget(rec, t, e0m, "psi2"),
        }

    s4 = json.loads((OUT / "s4_axial_ned.json").read_text())
    ned = s4["comparison_gw150914_cell"]["ned"]["dx0.1"]
    for key, name in (("A", "A_ned"), ("D", "D_coupling_off")):
        rows = per_pulse(key)
        fake = {
            "pulses": [
                {
                    "hrss_ratio": rows[k]["hrss_ratio_strain"],
                    "arrival_after_free_flight_M": ned["pulses"][k]["arrival_after_free_flight_M"],
                }
                for k in (0, 1)
            ],
            "spacing_M": ned["spacing_M"],
        }
        out[name] = {"pulses": rows, "budget_vs_initial_energy": budgets(key, "p")}
        out[name]["pmw"] = an.pmw_eta_max(fake)
    out["B_smooth_packet_budget"] = budgets("B", "g")
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    scratch = Path(os.environ.get("X2_SCRATCH", "/tmp"))
    npz = scratch / "x2_s4b_waveforms.npz"
    a = cell_a()
    x_pk = an.gw_peak_x("af", a)
    if "--from-npz" in argv:
        z = np.load(npz)
        res = {
            k: {"t": z[f"{k}_t"], "rec": z[f"{k}_rec"], "a": a, "x_pk": x_pk, "x0": x_pk - 25.0}
            for k in "ABCD"
        }
    else:
        jobs = [("A", "ned", "p"), ("B", "ned", "g"), ("C", "ned", "gdot"), ("D", "af", "p")]
        with ProcessPoolExecutor(max_workers=min(len(jobs), os.cpu_count() or 1)) as ex:
            res = dict(ex.map(run, jobs))
        np.savez_compressed(
            npz,
            **{f"{k}_t": v["t"] for k, v in res.items()},
            **{f"{k}_rec": v["rec"] for k, v in res.items()},
        )
    out = analyse(res)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s4b_units.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
