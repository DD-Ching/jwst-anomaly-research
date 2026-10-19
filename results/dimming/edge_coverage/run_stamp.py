"""Post hoc sensitivity of run.py (gate V2 RESULT): injection-set pairs masked by stamp_edge
(gate_v2/pairs_<field>.ecsv) instead of the grid edge_pair; other pairs keep edge_pair.
Usage (repo root): run_stamp.py GAIA_ROOT OUT FIELD...

Writes OUT/<field>/{edge_summary.json, limits_<variant>.ecsv, efficiency_<variant>.ecsv}.
"""

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import dimming_screen as ds  # noqa: E402

from jwst_anomaly import pipeline  # noqa: E402

_spec = importlib.util.spec_from_file_location("edge", Path(__file__).with_name("edge.py"))
ed = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ed)

STEP = 32
COMMITTED = {
    "macs0416": "results/dimming/macs0416_regen",
    "abell2744": "results/dimming/abell2744_calib",
    "nexus": "results/dimming/nexus_calib",
}
NEXUS_EDGE_TARGETS = "results/dimming/nexus_o014/targets.ecsv"

gaia_root, out_root, fields = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
p = ds.Params()


def no_se(fi, ei, det, p, inflation):
    r, _ = ds._classify_injected(fi, ei, det, p)
    rec = ds.flagged(r)
    if inflation is not None:
        r2, _ = ds._classify_injected(fi, ei * np.asarray(inflation)[None, None, :], det, p)
        rec &= ds.flagged(r2)
    return rec


for name in fields:
    t0 = time.time()
    out = out_root / name
    out.mkdir(parents=True, exist_ok=True)
    field = ds.load_field(ds.DEFAULT_CONFIG, name)
    gaia = Table.read(gaia_root / name / "gaia.ecsv")
    cats, lc, scale, err, res, ordn = ds._screen(field, p, gaia)
    ne = len(field["epochs"])
    rel = np.full((len(lc["ra"]), ne), np.nan)
    cache = out / "wht_rel.npy"  # session cache (not tracked)
    for k, e in enumerate(field["epochs"]):
        if "F200W" not in e["bands"] or cache.exists():
            continue
        grid, wcs = ed.wht_grid(pipeline.l3_image_uri(field["cloud"], ds.obs_id(e, "F200W")), STEP)
        rel[:, k] = ed.wht_rel(grid, wcs, lc["ra"], lc["dec"], STEP)
        print(name, "epoch", k, "grid", grid.shape, f"{time.time() - t0:.0f}s", flush=True)
    if cache.exists():
        rel = np.load(cache)
    np.save(cache, rel)
    covered = np.isfinite(lc["flux"][:, :, 0])
    has_img = np.isfinite(rel)  # epochs without an F200W image are not tested
    edge = ed.edge_pair(covered & has_img, rel)
    pt = Table.read(f"results/dimming/edge_coverage/gate_v2/pairs_{name}.ecsv")
    pi = np.array([int(str(u).rsplit("-d", 1)[1]) for u in pt["uid"]], int)
    edge[pi, np.asarray(pt["epoch"], int)] = np.asarray(pt["stamp_edge"], bool)
    sel = ds.monitored(lc, err, ordn, p)
    inj = sel & ~ds.flagged(res)
    summ = {
        "field": name,
        "step_px": STEP,
        "n_master": len(lc["ra"]),
        "n_injection_sources": int(inj.sum()),
        "covered_pairs_injection": int((covered & has_img)[inj].sum()),
        "edge_pairs_injection": int(edge[inj].sum()),
        "f_edge_injection": float(edge[inj].sum() / max((covered & has_img)[inj].sum(), 1)),
        "injection_sources_with_edge_pair": int(edge[inj].any(axis=1).sum()),
        "covered_pairs_all": int((covered & has_img).sum()),
        "edge_pairs_all": int(edge.sum()),
        "edge_pairs_per_epoch_injection": [int(x) for x in edge[inj].sum(axis=0)],
    }
    if name == "nexus":  # gate V
        tg = Table.read(NEXUS_EDGE_TARGETS)
        tg = tg[np.asarray(tg["cut_edge"]) | np.asarray(tg["cut_low_weight"])]
        uid = np.array([f"{name}-d{i:05d}" for i in range(len(lc["ra"]))])
        called = []
        for u in np.asarray(tg["uid"], str):
            i = int(np.flatnonzero(uid == u)[0])
            called.append(bool(np.any(edge[i] & ~lc["det"][i, :, 0])))
        summ["gate_v_targets"] = [str(u) for u in tg["uid"]]
        summ["gate_v_called"] = called
        summ["gate_v_n_called"] = int(sum(called))
        summ["gate_v_passed"] = ed.gate_v(int(sum(called)))
    # what-if: edge pairs uncovered in both bands (ASSUMPTION: F444W footprint follows F200W)
    lcm = {**lc, "flux": lc["flux"].copy(), "det": lc["det"].copy()}
    errm = err.copy()
    lcm["flux"][edge] = np.nan
    lcm["det"][edge] = False
    errm[edge] = np.nan
    resm = ds.classify(ds.mask_shallow_nondetections(lcm["flux"], errm, lcm["det"], p), errm, p)
    selm = ds.monitored(lcm, errm, ordn, p)
    summ["n_injection_sources_masked"] = int((selm & ~ds.flagged(resm)).sum())
    inj_meta = json.loads(Path(f"{COMMITTED[name]}/inject_summary.json").read_text())
    infl = np.array(inj_meta["forced_inflation"], float)
    static_veto = ds.catalogue_veto({**ordn, "single_epoch": np.zeros_like(ordn["single_epoch"])})
    times = ds._times_yr(field)
    orig = ds._recovered
    for tag, (L, E, S), fn in (
        ("veto_on", (lc, err, sel), orig),
        ("veto_off", (lc, err, sel), no_se),
        ("veto_on_masked", (lcm, errm, selm), orig),
        ("veto_off_masked", (lcm, errm, selm), no_se),
    ):
        ds._recovered = fn
        eff = ds.efficiency_table(
            L,
            E,
            static_veto,
            S,
            times,
            ds.T_E_GRID,
            ds.RHO_GRID,
            ds.DEPTHS,
            p,
            seed=0,
            forced_inflation=infl,
        )
        first = eff[
            (eff["model"] == "w3")
            & (eff["t_e_yr"] == eff["t_e_yr"][0])
            & (eff["rho"] == eff["rho"][0])
        ]
        covd = np.isfinite(L["flux"][:, :, 0]).sum(axis=1)
        lim = ds.rate_limits(
            eff, ds.area_deg2(L["ra"], L["dec"], covd), int(np.sum(first["n_sources"]))
        )
        eff.write(out / f"efficiency_{tag}.ecsv", overwrite=True)
        lim.write(out / f"limits_{tag}.ecsv", overwrite=True)
    ds._recovered = orig
    summ["wall_s"] = round(time.time() - t0)
    (out / "edge_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summ.items() if not k.startswith("gate_v_t")}), flush=True)
