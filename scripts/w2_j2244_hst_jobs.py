"""J2244-0550 HST test job lists: development and validation (PLAN.md section 3 and v2)."""

import json
import pickle
import sys

S = sys.argv[1]
G_TARGET = 18.0  # mean Gaia G of the two J2244-0550 images (17.80, 18.21)
pos = json.load(open(f"{S}/data/pos_jobs.json"))


def negs(band, prog, ks, role):
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
    print(band, prog, role, "window stars", len(idx), "pairs", len(out))
    return out


dev = (
    pos["dev"] + negs("F160W", "17455", [0, 1, 2], "dev") + negs("F814W", "17308", [0, 1, 2], "dev")
)
val = (
    pos["val"]
    + negs("F160W", "17455", [3, 4, 5, 6, 7, 8], "val")
    + negs("F814W", "17308", [3, 4, 5, 6, 7, 8], "val")
)
json.dump(dev, open(f"{S}/data/dev_jobs.json", "w"), indent=1)
json.dump(val, open(f"{S}/data/val_jobs_unthresholded.json", "w"), indent=1)
print(
    [
        (j["name"], j["band"], [round(g, 2) for g in j["g"]], j["sat"])
        for j in dev + val
        if j["kind"] == "negative"
    ]
)


def injection_jobs(val: list[dict], zp: dict) -> list[dict]:
    """PLAN section 5: Sersic injections into validation negative pairs 3 and 4 of each band."""
    import w2_j2244_hst_rules as rules

    out = []
    k = 0
    for band in ("F160W", "F814W"):
        base = [
            j
            for j in val
            if j["kind"] == "negative" and j["band"] == band and j["name"] in ("pair3", "pair4")
        ]
        for b in base:
            for r_eff in rules.P.inj_reff:
                for mag in rules.P.inj_mags[band]:
                    j = dict(b)
                    j["kind"] = "injection"
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
