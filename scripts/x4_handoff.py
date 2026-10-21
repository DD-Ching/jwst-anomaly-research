"""X4 handoff: a self-contained reproduction pack of the U fitting-factor problem (S1k C1, C2).

Unit file: docs/state/units/X4.md (handoff round). `simulated` / `derived`; no data. Needs lalsuite
(the X4 environment). Writes only to results/x4/handoff/; reads the stored S1k files and never runs
a search. Truths are used only to check the score and to compute the recovery error.

  python scripts/x4_handoff.py pack        # -> results/x4/handoff/case.json, targets.npz
  python scripts/x4_handoff.py pycbc       # -> results/x4/handoff/pycbc_check.json (PyCBC env)
"""

from __future__ import annotations

import importlib.util
import json
import platform
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


OUT_DIR = _DIR.parent / "results" / "x4" / "handoff"


def _modules():
    k = _load("x4_s1k")
    return k, k.e, k.e.b


def points():
    """Truth and stored final solution of C1 and C2 (box parameters) and their targets."""
    k, e, _ = _modules()
    truth = json.loads((e.OUT / "s1k_truth.json").read_text())
    rows = {}
    for r in json.loads((e.OUT / "s1k_validate.json").read_text())["rows"]:
        rows.setdefault(r["name"], r)  # the frozen first run of each case
    out = {}
    for n in ("C1", "C2"):
        spec, x_truth = k.c_spec(n)
        assert np.allclose(x_truth, truth[n]["x_box"])
        t, _ = k.target(n)
        out[n] = {
            "t": t,
            "spec": spec,
            "x_truth": np.asarray(x_truth),
            "x_found": np.asarray(rows[n]["x"]),
            "stored_eps": rows[n]["eps_final"],
        }
    return out


def pack():
    import lal
    import lalsimulation as ls
    import scipy

    k, e, b = _modules()
    pts = points()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    case = {
        "label": "derived; X4 U fitting-factor problem, reproduction pack (simulated targets)",
        "quantity": (
            "FF(x) = max over time shift (inverse FFT, pad 64 for final scores, 32 in searches), "
            "polarisation angle psi and amplitude > 0 of "
            "Re<t, a h+ + c hx> / (||t|| ||a h+ + c hx||), "
            "noise-weighted <a, b> = sum_f a conj(b) / S_n(f) over 20-1024 Hz; "
            "no free overall phase "
            "(the reference phase is a parameter). Analytic: max_{a,c} = sqrt(v^T G^-1 v) / ||t||, "
            "v = (Re<t,h+>, Re<t,hx>), G = Gram matrix of (h+, hx). Objective = 1 - FF."
        ),
        "code": "scripts/x4_s1d_de.py: stat(), objective(), hpc(), to_cart(); BOUNDS",
        "waveform": (
            "IMRPhenomXPHM via lalsimulation.SimInspiralChooseFDWaveform, "
            "distance 100 Mpc, f_ref = f_low = 20 Hz"
        ),
        "parameters": {
            "names": [
                "ln Mc",
                "q = m1/m2",
                "|chi1|",
                "cos tilt1",
                "azimuth1",
                "|chi2|",
                "cos tilt2",
                "azimuth2",
                "inclination",
                "phi_ref",
            ],
            "bounds": [list(map(float, bd)) for bd in e.BOUNDS],
            "periodic": [4, 7, 9],
            "spin_representation": (
                "Cartesian chi = |chi| (sin t cos az, sin t sin az, cos t) at f_ref, frame of L"
            ),
        },
        "noise": "lalsimulation.SimNoisePSDaLIGODesignSensitivityP1200087, single detector",
        "band_hz": [b.F_LOW, b.F_HIGH],
        "df_hz": b.DF,
        "accuracy_rule": (
            "injection recovery: eps = 1 - FF(found) <= 5e-4 and <= 0.1 (1 - FF_case), "
            "FF_case 0.915149"
        ),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "lal": lal.__version__,
            "lalsimulation": ls.__version__,
        },
        "reproduce": [
            "uv venv lalenv --python 3.12; uv pip install --python lalenv/bin/python "
            "lalsuite==7.26.16 numpy==2.5.4 scipy==1.18.1 mpmath dynesty",
            "PYTHONPATH=src lalenv/bin/python scripts/x4_handoff.py pack",
            "the failing search: set DEADLINE_UTC in scripts/x4_s1k.py to a future time, "
            "copy results/x4/s1k_validate.json, then "
            "PYTHONPATH=src lalenv/bin/python scripts/x4_s1k.py one_v C1 "
            "(about 35 min on 2 processes)",
        ],
        "cases": {},
    }
    npz = {"freqs_hz": b._FREQS, "weight_1_over_psd": b.weight()}
    for n, p in pts.items():
        t = p["t"]
        rec = {
            "truth_box": p["x_truth"].tolist(),
            "found_box": p["x_found"].tolist(),
            "stored_eps": p["stored_eps"],
            "injection": {
                "mc": p["spec"][0][0],
                "q": p["spec"][0][1],
                "chi1": list(p["spec"][0][2]),
                "chi2": list(p["spec"][0][3]),
                "inc": p["spec"][0][4],
                "phi_ref": p["spec"][0][5],
                "psi": p["spec"][0][6],
            },
        }
        for lab, x in (("truth", p["x_truth"]), ("found", p["x_found"])):
            hp, hc = e.hpc(*e.to_cart(x))
            rec[f"FF_{lab}"] = {f"pad{pad}": e.stat(t, hp, hc, pad=pad) for pad in (32, 64, 256)}
            mc, q, s1, s2, inc, phi = e.to_cart(x)
            rec[f"{lab}_physical"] = {
                "mc": mc,
                "q": q,
                "chi1": list(s1),
                "chi2": list(s2),
                "inc": inc,
                "phi": phi,
            }
        rec["eps_found_pad64"] = 1.0 - rec["FF_found"]["pad64"]
        case["cases"][n] = rec
        npz[f"{n}_target"] = t
    (OUT_DIR / "case.json").write_text(json.dumps(case, indent=1, default=float) + "\n")
    np.savez_compressed(OUT_DIR / "targets.npz", **npz)
    return case


def pycbc_check(pad=64):
    """Score layer re-computed with PyCBC only (no X4 code): its IMRPhenomXPHM generator, matched
    filter and the no-phase sky-location statistic (arXiv:1709.09181), against case.json."""
    import pycbc
    import scipy
    from pycbc.filter import matched_filter, overlap_cplx, sigma
    from pycbc.filter.matchedfilter import compute_max_snr_over_sky_loc_stat_no_phase
    from pycbc.types import FrequencySeries
    from pycbc.waveform import get_fd_waveform

    case = json.loads((OUT_DIR / "case.json").read_text())
    z = np.load(OUT_DIR / "targets.npz")
    df, (flo, fhi) = case["df_hz"], case["band_hz"]
    w = z["weight_1_over_psd"]
    n = len(w)
    m = n * pad  # zero-padded length: finer time lags, as the X4 score's pad
    psd_arr = np.ones(m)
    psd_arr[:n][w > 0] = 1.0 / w[w > 0]
    psd = FrequencySeries(psd_arr, delta_f=df)

    def fs(a):
        out = np.zeros(m, complex)
        out[: min(len(a), n)] = a[:n]
        return FrequencySeries(out, delta_f=df)

    rows = {}
    for name, c in case["cases"].items():
        d = fs(z[f"{name}_target"])
        sd = sigma(d, psd=psd, low_frequency_cutoff=flo, high_frequency_cutoff=fhi)
        for lab in ("truth", "found"):
            ph = c[f"{lab}_physical"]
            q, mc = ph["q"], ph["mc"]
            mt = mc * (1 + q) ** 1.2 / q**0.6
            hp, hc = get_fd_waveform(
                approximant="IMRPhenomXPHM",
                mass1=mt * q / (1 + q),
                mass2=mt / (1 + q),
                spin1x=ph["chi1"][0],
                spin1y=ph["chi1"][1],
                spin1z=ph["chi1"][2],
                spin2x=ph["chi2"][0],
                spin2y=ph["chi2"][1],
                spin2z=ph["chi2"][2],
                inclination=ph["inc"],
                coa_phase=ph["phi"],
                distance=100.0,
                delta_f=df,
                f_lower=flo,
                f_ref=flo,
                f_final=fhi,
            )
            hp, hc = fs(np.asarray(hp)), fs(np.asarray(hc))
            kw = {"psd": psd, "low_frequency_cutoff": flo, "high_frequency_cutoff": fhi}
            rp, rc = matched_filter(hp, d, **kw), matched_filter(hc, d, **kw)
            corr = overlap_cplx(hp, hc, normalized=True, **kw).real
            stat = compute_max_snr_over_sky_loc_stat_no_phase(
                rp.numpy(), rc.numpy(), corr, hpnorm=1.0, hcnorm=1.0
            )
            ff = float(np.max(np.abs(np.asarray(stat)))) / sd
            ours = c[f"FF_{lab}"][f"pad{pad}"]
            rows[f"{name}_{lab}"] = {"FF_pycbc": ff, "FF_ours_pad64": ours, "diff": ff - ours}
    out = {
        "label": "derived; PyCBC re-computation of the X4 score at fixed points (no search)",
        "versions": {
            "pycbc": pycbc.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "pad": pad,
        "rows": rows,
        "max_abs_diff": max(abs(r["diff"]) for r in rows.values()),
        "rule": "agree if max_abs_diff <= 1e-4",
    }
    out["agree"] = out["max_abs_diff"] <= 1e-4
    (OUT_DIR / "pycbc_check.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "pack":
        c = pack()
        print(
            json.dumps(
                {
                    n: {k: v for k, v in r.items() if k.startswith(("FF_", "eps"))}
                    for n, r in c["cases"].items()
                },
                default=float,
            )
        )
    elif argv and argv[0] == "pycbc":
        print(json.dumps(pycbc_check(), indent=1, default=float))
    else:
        print(__doc__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
