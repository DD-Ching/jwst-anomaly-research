"""Weight-map coverage term: gates R0 / R1, corrected limits and flags (PLAN.md). Usage (repo root):
run.py GAIA_ROOT OUT FIELD...   (GAIA_ROOT/<field>/gaia.ecsv from persistence_scope/gaia_esa.py)

Writes OUT/<field>/{summary.json, limits_<variant>.ecsv, efficiency_<variant>.ecsv, wht_rel.npz}.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import dimming_screen as ds  # noqa: E402

COMMITTED = {
    "macs0416": "results/dimming/macs0416_regen",
    "abell2744": "results/dimming/abell2744_calib",
    "nexus": "results/dimming/nexus_calib",
}
WHATIF = "results/dimming/edge_coverage/runs/{}/limits_veto_on_masked.ecsv"

gaia_root, out_root, fields = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
p = ds.Params()


def no_se(fi, ei, det, p, inflation):
    r, _ = ds._classify_injected(fi, ei, det, p)
    rec = ds.flagged(r)
    if inflation is not None:
        r2, _ = ds._classify_injected(fi, ei * np.asarray(inflation)[None, None, :], det, p)
        rec &= ds.flagged(r2)
    return rec


def limits(field, lc, err, ordn, infl, fn):
    sel = ds.monitored(lc, err, ordn, p)
    static_veto = ds.catalogue_veto({**ordn, "single_epoch": np.zeros_like(ordn["single_epoch"])})
    orig = ds._recovered
    ds._recovered = fn or orig
    try:
        eff = ds.efficiency_table(
            lc,
            err,
            static_veto,
            sel,
            ds._times_yr(field),
            ds.T_E_GRID,
            ds.RHO_GRID,
            ds.DEPTHS,
            p,
            seed=0,
            forced_inflation=infl,
        )
    finally:
        ds._recovered = orig
    first = eff[
        (eff["model"] == "w3") & (eff["t_e_yr"] == eff["t_e_yr"][0]) & (eff["rho"] == eff["rho"][0])
    ]
    covd = np.isfinite(lc["flux"][:, :, 0]).sum(axis=1)
    area = ds.area_deg2(lc["ra"], lc["dec"], covd)
    lim = ds.rate_limits(eff, area, int(np.sum(first["n_sources"])))
    lim.meta.update(field=field["name"], area_deg2=area, calibrated=True)
    return eff, lim, sel


def same(a, b):
    return bool(
        np.allclose(
            np.asarray(a["exposure_source_yr"]),
            np.asarray(b["exposure_source_yr"]),
            rtol=1e-9,
            equal_nan=True,
        )
    )


for name in fields:
    t0 = time.time()
    out = out_root / name
    out.mkdir(parents=True, exist_ok=True)
    field = ds.load_field(ds.DEFAULT_CONFIG, name)
    field["name"] = name
    gaia = Table.read(gaia_root / name / "gaia.ecsv")
    infl = np.array(
        json.loads(Path(f"{COMMITTED[name]}/inject_summary.json").read_text())["forced_inflation"]
    )
    summ = {"field": name}
    flags = {}
    for tag, cache in (("off", None), ("on", out / "wht_rel.npz")):
        cats, lc, scale, err, res, ordn = ds._screen(field, p, gaia, cache)
        eff, lim, sel = limits(field, lc, err, ordn, infl, None)
        eff.write(out / f"efficiency_veto_on_{tag}.ecsv", overwrite=True)
        lim.write(out / f"limits_veto_on_{tag}.ecsv", overwrite=True)
        uid = np.array([f"{name}-d{k:05d}" for k in range(len(lc["ra"]))])
        fl = ds.flagged(res) & ~ds.catalogue_veto(ordn) & sel
        flags[tag] = set(uid[fl].tolist())
        summ[f"term_{tag}"] = {
            "n_monitored": int(sel.sum()),
            "n_injection_sources": int(np.sum(lim["n_monitored"][:1])),
            "weight_edge_pairs": int(lc["weight_edge"].sum()),
            "weight_edge_sources": int(lc["weight_edge"].any(axis=1).sum()),
            "n_flagged_compact_no_veto": int(fl.sum()),
            "area_deg2": float(lim.meta["area_deg2"]),
        }
        if tag == "off":
            summ["gate_R0"] = same(lim, Table.read(f"{COMMITTED[name]}/limits.ecsv"))
        else:
            summ["gate_R1"] = same(lim, Table.read(WHATIF.format(name)))
            eff2, lim2, _ = limits(field, lc, err, ordn, infl, no_se)
            eff2.write(out / "efficiency_veto_off_on.ecsv", overwrite=True)
            lim2.write(out / "limits_veto_off_on.ecsv", overwrite=True)
    summ["flags_removed_by_term"] = sorted(flags["off"] - flags["on"])
    summ["flags_added_by_term"] = sorted(flags["on"] - flags["off"])
    summ["wall_s"] = round(time.time() - t0)
    (out / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ), flush=True)
