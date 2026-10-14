"""X3 G1: batch vetting of the held-out GAPS flags (the /vet-candidate order, Gaia light curves).

Unit file: docs/state/units/X3.md (G1 section 7). Ordinary explanations are tested in the
methodology's order; each flag gets the first rule it fails. A flag that fails none is
"unexplained", never a lens claim.

POST HOC: the rules were written after the held-out run and after seeing summary statistics of the
flags (fitted t_E and u0, G range, colour verdicts); D2 was added after viewing the contact sheets
(review of #228). Every rule is applied to the dev W3 injections as well, so its cost is in the
efficiency (``efficiency``): event rules on the injected light curve, host rules on the host.

Event rules (light curve):
   A1: the dip and the spikes are each within one Gaia visit (epochs within VISIT_DAYS): a single
   scan's photometry (scan angle, crowding, cosmic ray) rather than an event lasting days.
   D1: N1neg t_E outside 1-2,000 d or u0 > 1.9 (caustic grazing only).
   D2: the fitted dip is sampled in fewer than two visits (W3's defining feature is not observed).
Host rules (catalogue):
   A2: G < 13, where Gaia's per-transit photometry (gating, saturation) is outside the dev
   calibration.
   B1: a Gaia DR3 neighbour within NEIGHBOUR_ARCSEC at most 2 mag fainter.
   C1: a DR3 variability class (vari_classifier_result) or a SIMBAD variable otype (SIMBAD_VAR).
   A real event could itself earn a DR3 class; that path is not calibrated (stated condition).
Order: A1, A2, B1, C1, D1, D2.

Outputs: results/x3/g1_heldout_vetting.ecsv.gz (one row per flag), g1_vet_efficiency.json,
g1_heldout_limit.json, and contact sheets under the data root (not committed).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

from jwst_anomaly import gaia_mulens as gm

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


xb = _load("x3_bulk")
g1 = xb.g1
x3 = g1.x3
w3 = g1.w3

VISIT_DAYS = 1.0
NEIGHBOUR_ARCSEC = 2.0
BRIGHT_G = 13.0
TE_LIMITS = (1.0, 2000.0)
U0_MAX = 1.9
# SIMBAD otypes of variable stars and variable extragalactic sources (exact codes; "?" = candidate)
SIMBAD_VAR = frozenset(
    "V* Ir* Or* RI* Er* Fl* FU* RC* RC? Ro* a2* Psr Pu* RR* RR? Ce* Ce? cC* dS* RV* RV? WV* "
    "bC* sr* sr? SX* gD* LP* LP? Mi* Mi? EB* EB? El* El? BY* BY? RS* RS? Y*O Y*? TT* TT? Ae* "
    "Ae? Be* Be? No* No? CV* CV? XB* XB? AGN AGN? QSO QSO? Bla Bla? BLL BLL? Sy1 Sy2 SN*".split()
)
VAR_QUERY = (
    "SELECT g.source_id, g.phot_variable_flag, v.best_class_name, v.best_class_score "
    "FROM gaiadr3.gaia_source AS g LEFT OUTER JOIN gaiadr3.vari_classifier_result AS v "
    "ON g.source_id = v.source_id WHERE g.source_id IN ({ids})"
)
TE_RANGE = (30.0, 300.0)


def _key(ids) -> str:
    return hashlib.sha1(",".join(sorted(ids)).encode()).hexdigest()[:12]


def flags_table() -> Table:
    tabs = [Table.read(p) for p in sorted(x3.out_dir().glob("heldout_*.ecsv"))]
    for t in tabs:
        t.meta.clear()
    t = vstack(tabs)
    return t[np.asarray(t["fit_flag"], bool)]


def visits(t: np.ndarray) -> int:
    if t.size == 0:
        return 0
    return int(1 + np.sum(np.diff(np.sort(t)) > VISIT_DAYS))


def morph(lc, p: dict) -> dict:
    a = w3.magnification("N1neg", lc, p)
    X = np.vstack([a, np.ones_like(a)]).T
    fs, fb = np.linalg.lstsq(X * np.sqrt(lc.w)[:, None], lc.f * np.sqrt(lc.w), rcond=None)[0]
    model = fs * a + fb
    base = fs + fb
    dip = (model < base - 3 * lc.sf) & (lc.f < base - 3 * lc.sf)
    spike = (model > base + 3 * lc.sf) & (lc.f > base + 3 * lc.sf)
    return {
        "n_dip": int(dip.sum()),
        "dip_visits": visits(lc.t[dip]),
        "n_spike": int(spike.sum()),
        "spike_visits": visits(lc.t[spike]),
        "depth_frac": float((base - model.min()) / base) if base > 0 else np.nan,
        "fs_frac": float(fs / base) if base > 0 else np.nan,
        "model": model,
    }


def event_rule(m: dict, p: dict) -> str:
    """First failed event rule ("" if none), in the A1, D1, D2 sub-order."""
    if m["dip_visits"] <= 1 and m["spike_visits"] <= 1:
        return "A1 single-visit"
    if not (TE_LIMITS[0] <= p["tE"] <= TE_LIMITS[1]) or p["u0"] > U0_MAX:
        return "D1 degenerate fit"
    if m["dip_visits"] < 2:
        return "D2 no multi-visit dip"
    return ""


def host_rule(h: dict) -> str:
    if h["G"] < BRIGHT_G:
        return "A2 bright"
    if h["neighbour_dG"] <= 2.0:
        return "B1 neighbour"
    if h["vari_class"] or h["simbad_var"]:
        return "C1 known variable"
    return ""


def verdict(row: dict) -> str:
    p = {"tE": row["fit_tE"], "u0": row["fit_u0"]}
    ev = event_rule(row, p)
    if ev.startswith("A1"):
        return ev
    return host_rule(row) or ev or "unexplained"


def host_info(sids: list[str], ra, dec, gmag) -> dict[str, dict]:
    """Catalogue facts per host: DR3 class, SIMBAD otypes (2"), Gaia neighbours (CDS XMatch)."""
    import astropy.units as u
    from astroquery.xmatch import XMatch

    key = _key(sids)
    p_v = x3.out_dir() / f"vet_vari_{key}.csv"
    p_s = x3.out_dir() / f"vet_simbad_{key}.ecsv"
    p_g = x3.out_dir() / f"vet_gaia_nb_{key}.ecsv"
    if not p_v.exists():
        lines, header = [], None
        for k in range(0, len(sids), 500):
            rows = gm._tap_csv(VAR_QUERY.format(ids=",".join(sids[k : k + 500]))).strip()
            rows = rows.splitlines()
            header = rows[0]
            lines.extend(rows[1:])
        gm._write_atomic(p_v, ("\n".join([header, *lines]) + "\n").encode())
    up = Table({"sid": np.array(sids), "ra": np.asarray(ra, float), "dec": np.asarray(dec, float)})
    if not p_s.exists():
        XMatch.query(up, "simbad", 2.0 * u.arcsec, colRA1="ra", colDec1="dec").write(p_s)
    if not p_g.exists():
        res = XMatch.query(
            up, "vizier:I/355/gaiadr3", NEIGHBOUR_ARCSEC * u.arcsec, colRA1="ra", colDec1="dec"
        )
        res.write(p_g)
    cls = {}
    for r in Table.read(p_v, format="ascii.csv"):
        cls[str(r["source_id"])] = (
            "" if np.ma.is_masked(r["best_class_name"]) else str(r["best_class_name"])
        )
    missing = set(sids) - set(cls)
    if missing:
        raise RuntimeError(f"{len(missing)} ids missing from {p_v.name}")
    types: dict[str, list[str]] = {}
    for r in Table.read(p_s):
        types.setdefault(str(r["sid"]), []).append(str(r["otype"]))
    nb: dict[str, list[float]] = {}
    for r in Table.read(p_g):
        if str(r["Source"]) == str(r["sid"]):
            continue
        g = float(r["Gmag"]) if not np.ma.is_masked(r["Gmag"]) else 99.0
        nb.setdefault(str(r["sid"]), []).append(g)
    out = {}
    for s, g in zip(sids, gmag, strict=True):
        ty = types.get(s, [])
        out[s] = {
            "G": float(g),
            "vari_class": cls[s],
            "simbad_types": ";".join(ty),
            "simbad_var": any(t in SIMBAD_VAR for t in ty),
            "neighbour_dG": float(min(nb[s]) - g) if s in nb else 99.0,
        }
    return out


def contact_sheet(rows: list[dict], path: Path, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(rows)
    ncol = 6
    nrow = max(1, -(-n // ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.2 * nrow), squeeze=False)
    for ax in axes.flat:
        ax.set_visible(False)
    for ax, r in zip(axes.flat, rows, strict=False):
        ax.set_visible(True)
        lc, model = r["lc"], r["model"]
        o = np.argsort(lc.t)
        ax.errorbar(lc.t - 2457000, lc.f, lc.sf, fmt=".", ms=3, lw=0.5, color="k")
        ax.plot(lc.t[o] - 2457000, model[o], color="C3", lw=0.8)
        ax.set_title(
            f"{r['event_id'][-9:]} G{r['G']:.1f} {r['verdict'][:2]} dBIC{r['dbic']:.0f}", fontsize=7
        )
        ax.tick_params(labelsize=5)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def run(final_only: bool) -> Table:
    t = flags_table()
    ids = Table.read(x3.out_dir() / "heldout_ids_g19.csv", format="ascii.csv")
    info = {
        str(s): (float(r), float(d), float(m))
        for s, r, d, m in zip(
            ids["source_id"], ids["ra"], ids["dec"], ids["phot_g_mean_mag"], strict=True
        )
    }
    if final_only:
        t = t[np.asarray(t["final"], bool)]
    sids = [str(s) for s in t["event_id"]]
    hosts = host_info(
        sids, [info[s][0] for s in sids], [info[s][1] for s in sids], [info[s][2] for s in sids]
    )
    raw = x3.paths.data_root() / "raw" / "gaia_dr3_gaps_heldout_flags"
    rows = []
    for r, sid in zip(t, sids, strict=True):
        ra, dec, _g = info[sid]
        lc = x3.wg.gaia_light_curve(
            gm.parse_epoch_csv((raw / f"{sid}.csv").read_text()), sid, (ra, dec)
        )
        p = {k: float(r[f"fit_{k}"]) for k in ("t0", "tE", "u0", "rho")}
        m = morph(lc, p)
        row = {
            "event_id": sid,
            "ra": ra,
            "dec": dec,
            "dbic": float(r["dbic_N1neg"]),
            "fit_tE": p["tE"],
            "fit_u0": p["u0"],
            "achrom": str(r["achrom"]),
            **{k: v for k, v in m.items() if k != "model"},
            **hosts[sid],
        }
        row["verdict"] = verdict(row)
        rows.append({**row, "lc": lc, "model": m["model"]})
    keys = [k for k in rows[0] if k not in ("lc", "model")]
    out = Table(rows=[[r[k] for k in keys] for r in rows], names=keys)
    out.meta.update(
        provenance="observed light curves and catalogues; verdicts derived (ordinary tests only)",
        tests=__doc__.split("Outputs:")[0].strip(),
    )
    w3.write_ecsv_gz(out, x3.results_dir() / "g1_heldout_vetting.ecsv.gz")
    sheet = x3.out_dir() / "heldout_flags_sheet.png"
    order = sorted(rows, key=lambda r: (r["verdict"] != "unexplained", r["verdict"], r["dbic"]))
    contact_sheet(order, sheet, "X3 G1 held-out final flags (observed G flux, N1neg fit in red)")
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    print(json.dumps(counts, indent=1))
    print(f"contact sheet: {sheet}")
    return out


def efficiency(seed: int = 77) -> dict:
    """W3 efficiency of screen + vetting on dev injections into hosts that v1 does not veto.

    The exposure (held-out sources not v1-vetoed, >= MIN_POINTS) carries the v1 veto, so v1-vetoed
    hosts are left out of the denominator here. The injections are rebuilt deterministically with
    the seed of ``x3_g1.py run --set inject`` and matched to the stored fit parameters.
    """
    tab = Table.read(x3.out_dir() / "g1_inject.ecsv")
    kind = g1._strcol(tab, "kind")
    w = tab[kind == "W3_base"]
    stored = {str(r["event_id"]): r for r in w}
    dev = {s: (pos, g) for s, pos, g in g1.dev_sources()}
    sids = sorted(stored)
    hosts = host_info(
        sids, [dev[s][0][0] for s in sids], [dev[s][0][1] for s in sids], [dev[s][1] for s in sids]
    )
    by_te: dict[str, list[int]] = {}
    lost: dict[str, int] = {}
    n_all = n_found = 0
    for kind_, sid, lc, _b, _v, _sc, prm in g1.build_jobs("inject", seed):
        if kind_ != "W3_base":
            continue
        r = stored[sid]
        if abs(float(r["tE"]) - prm["tE"]) > 1e-9:
            raise RuntimeError(f"injection {sid} does not match the stored table")
        if bool(r["v1_vetoed"]):
            continue
        c = by_te.setdefault(f"{prm['tE']:g}", [0, 0])
        c[1] += 1
        why = "" if bool(r["final"]) else "screen"
        if not why:
            p = {k: float(r[f"fit_{k}"]) for k in ("t0", "tE", "u0", "rho")}
            ev = event_rule(morph(lc, p), p)
            why = ev if ev.startswith("A1") else (host_rule(hosts[sid]) or ev)
        if TE_RANGE[0] <= prm["tE"] <= TE_RANGE[1]:
            n_all += 1
            if why:
                lost[why] = lost.get(why, 0) + 1
            else:
                n_found += 1
        if not why:
            c[0] += 1
    rec = {
        "label": "derived: simulated W3 injected into observed dev light curves; screen + vetting",
        "denominator": "injections into hosts not vetoed by v1",
        "found_by_tE": by_te,
        "te_range": list(TE_RANGE),
        "n": n_all,
        "found": n_found,
        "efficiency": n_found / n_all,
        "efficiency_wilson95": g1.wilson(n_found, n_all),
        "lost_by_rule": lost,
        "host_rule_rate": float(np.mean([bool(host_rule(h)) for h in hosts.values()])),
    }
    (x3.results_dir() / "g1_vet_efficiency.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


POISSON_UP95 = {0: 3.0, 1: 4.74, 2: 6.3, 3: 7.75, 4: 9.15, 5: 10.51}


def rate_limit() -> dict:
    """95 % W3 rate limit per exposed held-out star and year after vetting (conditional)."""
    hs = json.loads((x3.results_dir() / "g1_heldout_summary.json").read_text())
    eff = json.loads((x3.results_dir() / "g1_vet_efficiency.json").read_text())["efficiency"]
    sel = json.loads((x3.results_dir() / "g1_selection_bound.json").read_text())
    vet = Table.read(x3.results_dir() / "g1_heldout_vetting.ecsv.gz", format="ascii.ecsv")
    n_ev = int(np.sum(np.array([str(v) for v in vet["verdict"]]) == "unexplained"))
    exp = hs["exposed_sources"] * xb.DR3_SPAN_YR
    keep = sel["retained_fraction_te30_300"]["1"]
    up = POISSON_UP95[n_ev]
    out = {
        "label": "derived 95 % upper limit; conditional (X3.md G1 sections 4, 7); not a detection",
        "te_range_d": list(TE_RANGE),
        "efficiency_after_vetting": eff,
        "exposed_star_years": exp,
        "n_unexplained": n_ev,
        "poisson_upper_95": up,
        "per_star_per_yr_photocentre_stable": up / (exp * eff),
        "retained_fraction_blended_prior": keep,
        "per_star_per_yr_blended_prior": up / (exp * eff * keep),
        "conditions": [
            "injection priors: u0 ~ U(0, 2), rho 0.01/0.1, f_s ~ U(0.5, 1), t0 uniform in span",
            "photocentre-stable: the event's own photocentre motion leaves RUWE <= 1.4",
            "a real event that earns a DR3 variability class is not modelled (C1)",
            "selection effects listed as not bounded in g1_selection_bound.json",
        ],
    }
    (x3.results_dir() / "g1_heldout_limit.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("cmd", choices=("flags", "efficiency", "limit"))
    ap.add_argument("--all-fit-flags", action="store_true", help="also vet flags removed earlier")
    a = ap.parse_args(argv)
    if a.cmd == "flags":
        run(final_only=not a.all_fit_flags)
    elif a.cmd == "efficiency":
        print(json.dumps(efficiency(), indent=1))
    else:
        print(json.dumps(rate_limit(), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
