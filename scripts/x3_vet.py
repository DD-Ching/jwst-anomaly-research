"""X3 G1: batch vetting of the held-out GAPS flags (the /vet-candidate order, Gaia light curves).

Unit file: docs/state/units/X3.md (G1 section 7). Ordinary explanations are tested in the
methodology's order; each flag gets the first test it fails. A flag that fails none is "unexplained
under A-D", never a lens claim.

POST HOC: the tests were written after the held-out run and after seeing summary statistics of the
flags (fitted t_E and u0, G range, colour verdicts), before viewing any flag's light curve. Their
cost on W3 injections is measured (g1_vet_rule_loss.json) and subtracted from the efficiency.

A. Processing / instrument.
   A1: the dip (or spike) is carried by one Gaia visit (epochs within VISIT_DAYS): a single scan's
   photometry (scan angle, crowding, cosmic ray) rather than an event lasting days.
   A2: G < 13, where Gaia's per-transit photometry (gating, saturation) is outside the dev
   calibration.
B. Catalogue. B1: a Gaia DR3 neighbour within NEIGHBOUR_ARCSEC at most 2 mag fainter.
C. Known populations. C1: a DR3 variability class (vari_classifier_result) or a SIMBAD variable
   type.
D. Fit degeneracy. D1: N1neg t_E outside 1-2,000 d or u0 > 1.9 (no dip, caustic grazing only).

Outputs: results/x3/g1_heldout_vetting.ecsv.gz (one row per flag), g1_heldout_limit.json, and
contact sheets under the data root (not committed).
"""

from __future__ import annotations

import argparse
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
SIMBAD_VAR = (
    "V*",
    "Ce",
    "cC",
    "dS",
    "RR",
    "EB",
    "LP",
    "Mi",
    "sr",
    "RG",
    "Ir",
    "Or",
    "WR",
    "Y*O",
    "BS",
    "Be",
    "RC",
    "No",
    "SN",
    "Pu",
    "s*",
    "AB",
    "C*",
    "S*",
    "WD",
    "El",
    "Er",
    "gD",
    "RV",
    "a2",
    "BY",
    "RS",
)
VAR_QUERY = (
    "SELECT g.source_id, g.phot_variable_flag, g.ruwe, g.phot_bp_rp_excess_factor, "
    "v.best_class_name, v.best_class_score FROM gaiadr3.gaia_source AS g "
    "LEFT OUTER JOIN gaiadr3.vari_classifier_result AS v "
    "ON g.source_id = v.source_id WHERE g.source_id IN ({ids})"
)


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


def morphology(sid: str, pos, r) -> dict:
    text = (x3.paths.data_root() / "raw" / "gaia_dr3_gaps_heldout_flags" / f"{sid}.csv").read_text()
    lc = x3.wg.gaia_light_curve(gm.parse_epoch_csv(text), sid, pos)
    p = {k: float(r[f"fit_{k}"]) for k in ("t0", "tE", "u0", "rho")}
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
        "lc": lc,
        "model": model,
    }


def catalogue(ids: list[str]) -> dict:
    path = x3.out_dir() / "heldout_flags_vari.csv"
    if not path.exists():
        text = gm._tap_csv(VAR_QUERY.format(ids=",".join(ids)))
        gm._write_atomic(path, text.encode())
    t = Table.read(path, format="ascii.csv")
    out = {}
    for r in t:
        cls = "" if np.ma.is_masked(r["best_class_name"]) else str(r["best_class_name"])
        out[str(r["source_id"])] = {
            "vari_class": cls,
            "phot_variable_flag": str(r["phot_variable_flag"]),
        }
    return out


def xmatch(tab: Table) -> tuple[dict, dict]:
    """SIMBAD main type within 2" and Gaia DR3 neighbours within NEIGHBOUR_ARCSEC (CDS XMatch)."""
    import astropy.units as u
    from astroquery.xmatch import XMatch

    path_s = x3.out_dir() / "heldout_flags_simbad.ecsv"
    path_g = x3.out_dir() / "heldout_flags_gaia_nb.ecsv"
    up = Table({"sid": tab["event_id"], "ra": tab["ra"], "dec": tab["dec"]})
    if not path_s.exists():
        XMatch.query(up, "simbad", 2.0 * u.arcsec, colRA1="ra", colDec1="dec").write(path_s)
    if not path_g.exists():
        res = XMatch.query(
            up, "vizier:I/355/gaiadr3", NEIGHBOUR_ARCSEC * u.arcsec, colRA1="ra", colDec1="dec"
        )
        res.write(path_g)
    s, g = Table.read(path_s), Table.read(path_g)
    simbad = {}
    for r in s:
        simbad.setdefault(str(r["sid"]), []).append(str(r["main_type"]))
    nb = {}
    for r in g:
        if str(r["Source"]) == str(r["sid"]):
            continue
        nb.setdefault(str(r["sid"]), []).append(
            float(r["Gmag"]) if not np.ma.is_masked(r["Gmag"]) else 99.0
        )
    return simbad, nb


def verdict(row: dict) -> str:
    if row["dip_visits"] <= 1 and row["spike_visits"] <= 1:
        return "A1 single-visit"
    if row["G"] < BRIGHT_G:
        return "A2 bright"
    if row["neighbour_dG"] <= 2.0:
        return "B1 neighbour"
    if row["vari_class"] or row["simbad_var"]:
        return "C1 known variable"
    if not (TE_LIMITS[0] <= row["fit_tE"] <= TE_LIMITS[1]) or row["fit_u0"] > U0_MAX:
        return "D1 degenerate fit"
    return "unexplained under A-D"


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
    t["ra"] = [info[str(s)][0] for s in t["event_id"]]
    t["dec"] = [info[str(s)][1] for s in t["event_id"]]
    sids = [str(s) for s in t["event_id"]]
    cat = catalogue(sids)
    simbad, nb = xmatch(t)
    rows = []
    for r in t:
        sid = str(r["event_id"])
        ra, dec, gmag = info[sid]
        m = morphology(sid, (ra, dec), r)
        types = simbad.get(sid, [])
        nbg = nb.get(sid, [])
        row = {
            "event_id": sid,
            "ra": ra,
            "dec": dec,
            "G": gmag,
            "dbic": float(r["dbic_N1neg"]),
            "fit_tE": float(r["fit_tE"]),
            "fit_u0": float(r["fit_u0"]),
            "achrom": str(r["achrom"]),
            **{k: v for k, v in m.items() if k not in ("lc", "model")},
            "vari_class": cat.get(sid, {}).get("vari_class", ""),
            "simbad_types": ";".join(types),
            "simbad_var": any(any(ty.startswith(v) for v in SIMBAD_VAR) for ty in types),
            "neighbour_dG": float(min(nbg) - gmag) if nbg else 99.0,
        }
        row["verdict"] = verdict(row)
        rows.append({**row, "lc": m["lc"], "model": m["model"]})
    keys = [k for k in rows[0] if k not in ("lc", "model")]
    out = Table(rows=[[r[k] for k in keys] for r in rows], names=keys)
    out.meta.update(
        provenance="observed light curves and catalogues; verdicts derived (ordinary tests only)",
        tests=__doc__.split("Outputs:")[0].strip(),
    )
    w3.write_ecsv_gz(out, x3.results_dir() / "g1_heldout_vetting.ecsv.gz")
    sheet = x3.out_dir() / "heldout_flags_sheet.png"
    order = sorted(
        rows, key=lambda r: (r["verdict"] != "unexplained under A-D", r["verdict"], r["dbic"])
    )
    contact_sheet(order, sheet, "X3 G1 held-out final flags (observed G flux, N1neg fit in red)")
    unexpl = [r for r in order if r["verdict"] == "unexplained under A-D"]
    if unexpl:
        contact_sheet(
            unexpl, x3.out_dir() / "heldout_unexplained_sheet.png", "unexplained under A-D"
        )
    counts = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    print(json.dumps(counts, indent=1))
    print(f"contact sheet: {sheet}")
    return out


def rule_loss(seed: int = 77) -> dict:
    """Cost of the post-hoc rules A1 and D1 on simulated W3 injections that pass the frozen screen.

    The injections are rebuilt with the seed of `x3_g1.py run --set inject` (deterministic) and
    matched to the stored fit parameters. Catalogue rules (A2, B1, C1) depend on the host, not
    the event.
    """
    tab = Table.read(x3.out_dir() / "g1_inject.ecsv")
    kind = g1._strcol(tab, "kind")
    fin = tab[(kind == "W3_base") & np.asarray(tab["final"], bool)]
    par = {str(r["event_id"]): r for r in fin}
    res = {"n": 0, "A1": 0, "D1": 0, "A1_by_tE": {}}
    for kind_, sid, lc, _bands, _vet, _sc, prm in g1.build_jobs("inject", seed):
        if kind_ != "W3_base" or sid not in par:
            continue
        r = par[sid]
        if abs(float(r["tE"]) - prm["tE"]) > 1e-9:
            continue
        p = {k: float(r[f"fit_{k}"]) for k in ("t0", "tE", "u0", "rho")}
        a = w3.magnification("N1neg", lc, p)
        X = np.vstack([a, np.ones_like(a)]).T
        fs, fb = np.linalg.lstsq(X * np.sqrt(lc.w)[:, None], lc.f * np.sqrt(lc.w), rcond=None)[0]
        model = fs * a + fb
        base = fs + fb
        dip = (model < base - 3 * lc.sf) & (lc.f < base - 3 * lc.sf)
        spike = (model > base + 3 * lc.sf) & (lc.f > base + 3 * lc.sf)
        a1 = visits(lc.t[dip]) <= 1 and visits(lc.t[spike]) <= 1
        d1 = not (TE_LIMITS[0] <= p["tE"] <= TE_LIMITS[1]) or p["u0"] > U0_MAX
        res["n"] += 1
        res["A1"] += int(a1)
        res["D1"] += int(d1 and not a1)
        c = res["A1_by_tE"].setdefault(f"{prm['tE']:g}", [0, 0])
        c[0] += int(a1)
        c[1] += 1
    res["label"] = (
        "derived: post-hoc vetting rules A1/D1 applied to simulated W3 injections "
        "passing the frozen screen"
    )
    (x3.results_dir() / "g1_vet_rule_loss.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def rate_limit() -> dict:
    """95 % W3 rate limit per held-out star and year after vetting (conditional; X3.md G1).

    Efficiency: frozen-screen W3 passes at t_E = 30-300 d minus those the post-hoc rules A1 and D1
    remove (g1_vet_rule_loss.json; D1 losses at any t_E subtracted: conservative). Two counts:
    0 events (the three "unexplained under A-D" flags judged not W3-like on inspection) and 3.
    """
    hs = json.loads((x3.results_dir() / "g1_heldout_summary.json").read_text())
    cal = json.loads((x3.results_dir() / "g1_summary.json").read_text())
    loss = json.loads((x3.results_dir() / "g1_vet_rule_loss.json").read_text())
    sel = json.loads((x3.results_dir() / "g1_selection_bound.json").read_text())
    k, n = cal["inject:W3_base"]["eff_te30_300"]
    a1 = sum(loss["A1_by_tE"][t][0] for t in ("30", "100", "300"))
    eff = (k - a1 - loss["D1"]) / n
    exp = hs["exposed_sources"] * xb.DR3_SPAN_YR
    keep = sel["retained_fraction_mean"]["1"]
    out = {
        "label": "derived 95 % upper limits; conditional (X3.md G1 sections 4, 7); not a detection",
        "efficiency_te30_300_after_vetting": eff,
        "exposed_star_years": exp,
        "retained_fraction_blended_prior": keep,
    }
    for n_ev, up in ((0, 3.0), (3, 7.75)):
        out[f"n_events_{n_ev}"] = {
            "poisson_upper_95": up,
            "per_star_per_yr_photocentre_stable": up / (exp * eff),
            "per_star_per_yr_blended_prior": up / (exp * eff * keep),
        }
    (x3.results_dir() / "g1_heldout_limit.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--all-fit-flags", action="store_true", help="also vet flags removed by later vetoes"
    )
    ap.add_argument("--limit", action="store_true", help="only recompute the rate limit")
    ap.add_argument("--rule-loss", action="store_true", help="recompute the rules' cost on W3")
    a = ap.parse_args(argv)
    if a.rule_loss:
        print(json.dumps(rule_loss(), indent=1))
    elif not a.limit:
        run(final_only=not a.all_fit_flags)
    print(json.dumps(rate_limit(), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
