"""J2244-0550 HST test job lists (PLAN.md sections 3, 5 and 6; v2 / v3).

    python scripts/w2_j2244_hst_jobs.py build <scratch>          positives, development, validation
    python scripts/w2_j2244_hst_jobs.py threshold <scratch> <dev_results.json>
                                                                  T into validation + injection jobs
    python scripts/w2_j2244_hst_jobs.py target <scratch> <dev_results.json>
                                                                  target + target-injection jobs

Run in the isolated environment (it imports w2_j2244_hst for the target photometric zero-points).
"""

import json
import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w2_j2244_hst_rules as rules  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G_TARGET = 18.0  # mean Gaia G of the two J2244-0550 images (17.80, 18.21)
TARGET = {"name": "J2244-0550", "ra": 341.0473, "dec": -5.8471}
TARGET_OBS = {"F160W": ("if8a5k030", "17455"), "F814W": ("jf56c1010", "17308")}
POS = {  # PLAN section 3 pools (development: J0818, J0743)
    "F160W": [
        "SDSSJ0818+0601",
        "ULASJ0743+2457",
        "SDSSJ0806+2006",
        "SDSSJ1620+1203",
        "VDESJ2325-5229",
    ],
    "F814W": [
        "SDSSJ0818+0601",
        "ULASJ0743+2457",
        "SDSSJ0114+0722",
        "J0140-1152",
        "ULASJ0820+0812",
        "SDSSJ0946+1835",
        "SDSSJ1640+1932",
        "WGD2021-4115",
        "SDSSJ2146-0047",
        "144320-012537",
        "DESI-011.9124+03.3756",
        "022418+024703",
    ],
}
PROGS = {"F160W": ("17455", "17199"), "F814W": ("17308",)}
DEV = ("SDSSJ0818+0601", "ULASJ0743+2457")


def positives(S: str) -> dict:
    from astropy.table import Table

    c = Table.read(ROOT / "results/w2_hst_blend/classes.ecsv")
    cov = Table.read(ROOT / "results/w2_hst_blend/hst_coverage.ecsv")
    cov = cov[cov["contains"]]
    jobs = {"dev": [], "val": []}
    for band in ("F160W", "F814W"):
        order = []
        for n in POS[band]:
            r = c[c["name"] == n][0]
            rows = sorted(
                (
                    x
                    for x in cov
                    if x["system_id"] == r["system_id"]
                    and x["proposal_id"] in PROGS[band]
                    and band in x["filters"]
                    and not x["obs_id"].startswith("hst_")
                    and x["t_exptime"] > 100
                    and (
                        x["instrument"].startswith("WFC3/IR")
                        if band == "F160W"
                        else x["instrument"] == "ACS/WFC"
                    )
                ),
                key=lambda x: x["obs_id"],
            )
            x = rows[0]
            order.append(
                (
                    str(r["system_id"]),
                    n,
                    str(x["obs_id"]),
                    str(x["proposal_id"]),
                    float(r["ra"]),
                    float(r["dec"]),
                )
            )
        for sid, n, oid, prog, ra, dec in sorted(order):
            role = "dev" if n in DEV else "val"
            jobs[role].append(
                {
                    "kind": "positive",
                    "name": n,
                    "system_id": sid,
                    "band": band,
                    "obs_id": oid,
                    "programme": prog,
                    "ra": ra,
                    "dec": dec,
                    "library": f"{S}/data/lib_{prog}_{band}.pkl",
                    "seed": 1,
                    "png": f"png/{role}_{band}_{n}.png",
                }
            )
    return jobs


def negs(S: str, band: str, prog: str, ks: list[int], role: str) -> list[dict]:
    lib = pickle.load(open(f"{S}/data/lib_{prog}_{band}.pkl", "rb"))
    idx = sorted(
        [i for i, s in enumerate(lib) if abs(s["g"] - G_TARGET) <= 0.75],
        key=lambda i: (abs(lib[i]["g"] - G_TARGET), lib[i]["source_id"]),
    )
    out = []
    for k in ks:
        if 2 * k + 1 >= len(idx):
            break
        a, b = idx[2 * k], idx[2 * k + 1]
        out.append(
            {
                "kind": "negative",
                "name": f"pair{k}",
                "band": band,
                "programme": prog,
                "i": a,
                "j": b,
                "sep": 1.02,
                "pa": 37.0 * k,
                "g": [lib[a]["g"], lib[b]["g"]],
                "sat": [lib[a]["sat"], lib[b]["sat"]],
                "stars": [lib[a]["source_id"], lib[b]["source_id"]],
                "library": f"{S}/data/lib_{prog}_{band}.pkl",
                "seed": 1,
                "png": f"png/{role}_{band}_pair{k}.png",
            }
        )
    return out


def build(S: str) -> None:
    pos = positives(S)
    dev = (
        pos["dev"]
        + negs(S, "F160W", "17455", [0, 1, 2], "dev")
        + negs(S, "F814W", "17308", [0, 1, 2], "dev")
    )
    val = pos["val"] + negs(S, "F160W", "17455", [3, 4, 5, 6, 7, 8], "val")
    val += negs(S, "F814W", "17308", [3, 4, 5, 6, 7, 8], "val")
    json.dump(dev, open(f"{S}/data/dev_jobs.json", "w"), indent=1)
    json.dump(val, open(f"{S}/data/val_jobs_unthresholded.json", "w"), indent=1)
    print("dev", len(dev), "val", len(val))


def thresholds(dev_results: str) -> dict:
    d = json.load(open(dev_results))
    return {
        b: rules.threshold(
            [r["dchi2"] for r in d if r["band"] == b and r["kind"] == "negative" and r.get("ok")]
        )
        for b in ("F160W", "F814W")
    }


def zero_points() -> dict:
    import w2_j2244_hst as m
    from astropy.io import fits

    zp = {}
    for band, (oid, _) in TARGET_OBS.items():
        with fits.open(m.url_of(oid), use_fsspec=True, fsspec_kwargs=m.FSSPEC) as h:
            zp[band] = float(
                m.zp_ab(float(h["SCI"].header["PHOTFLAM"]), float(h["SCI"].header["PHOTPLAM"]))
            )
    return zp


def injection_jobs(base_jobs: list[dict], zp: dict, names=("pair3", "pair4")) -> list[dict]:
    """PLAN section 5: Sersic injections into base stamps (negative pairs or the target)."""
    out = []
    k = 0
    for band in ("F160W", "F814W"):
        base = [j for j in base_jobs if j["band"] == band and j["name"] in names]
        for b in base:
            for r_eff in rules.P.inj_reff:
                for mag in rules.P.inj_mags[band]:
                    j = dict(b)
                    j["kind"] = "injection" if b["kind"] == "negative" else "target_injection"
                    j["name"] = f"{b['name']}_r{r_eff}_m{mag}"
                    j["png"] = None
                    j["zp_ab"] = zp[band]
                    j["inject"] = {
                        "mag": mag,
                        "r_eff": r_eff,
                        "frac": rules.P.inj_frac,
                        "seed": 1000 + k,
                    }
                    out.append(j)
                    k += 1
    return out


def threshold_cmd(S: str, dev_results: str) -> None:
    T = thresholds(dev_results)
    val = json.load(open(f"{S}/data/val_jobs_unthresholded.json"))
    for j in val:
        j["T"] = T[j["band"]]
    json.dump(val, open(f"{S}/data/val_jobs.json", "w"), indent=1)
    inj = injection_jobs(val, zero_points())
    for j in inj:
        j["T"] = T[j["band"]]
    json.dump(inj, open(f"{S}/data/inj_jobs.json", "w"), indent=1)
    print("T", T, "val", len(val), "inj", len(inj))


def target_cmd(S: str, dev_results: str, bands: list[str]) -> None:
    T = thresholds(dev_results)
    zp = zero_points()
    tj = []
    for band in bands:
        oid, prog = TARGET_OBS[band]
        tj.append(
            {
                "kind": "target",
                "name": TARGET["name"],
                "band": band,
                "obs_id": oid,
                "programme": prog,
                "ra": TARGET["ra"],
                "dec": TARGET["dec"],
                "library": f"{S}/data/lib_{prog}_{band}.pkl",
                "seed": 1,
                "png": f"png/target_{band}.png",
                "T": T[band],
                "zp_ab": zp[band],
            }
        )
    inj = injection_jobs(tj, zp, names=(TARGET["name"],))
    for j in inj:
        j["T"] = T[j["band"]]
    json.dump(tj, open(f"{S}/data/target_jobs.json", "w"), indent=1)
    json.dump(inj, open(f"{S}/data/target_inj_jobs.json", "w"), indent=1)
    print("target", len(tj), "target injections", len(inj))


if __name__ == "__main__":
    cmd, S = sys.argv[1], sys.argv[2]
    if cmd == "build":
        build(S)
    elif cmd == "threshold":
        threshold_cmd(S, sys.argv[3])
    elif cmd == "target":
        target_cmd(S, sys.argv[3], sys.argv[4].split(","))
