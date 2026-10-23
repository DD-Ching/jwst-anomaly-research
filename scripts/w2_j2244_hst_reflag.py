"""v3 (review fix): re-derive the F814W library saturation flags with the per-exposure time
from the SCI-header NCOMBINE, from the stored star stamps (no image re-download). Stars that the
corrected rule would reject are dropped; stars the wrong rule rejected are missing."""

import json
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w2_j2244_hst as m  # noqa: E402


def main(lib_in: str, lib_out: str, summary_json: str) -> None:
    from astropy.io import fits

    lib = pickle.load(open(lib_in, "rb"))
    half = m.BAND["F814W"]["star_half"]
    sr = m.STAR_RULES["F814W"]
    per_exp = {}
    for oid in sorted({s["obs_id"] for s in lib}):
        with fits.open(m.url_of(oid), use_fsspec=True, fsspec_kwargs=m.FSSPEC) as h:
            ncomb = int(h["SCI"].header.get("NCOMBINE", h[0].header.get("NCOMBINE", 1)))
            per_exp[oid] = float(h[0].header["EXPTIME"]) / max(ncomb, 1)
    out, changes = [], {"sat_to_unsat": 0, "unsat_to_sat": 0, "dropped": 0, "kept": 0}
    for s in lib:
        st, wht = s["stamp"], s["wht"]
        medw = np.median(wht)
        core = wht[half - 3 : half + 4, half - 3 : half + 4]
        sat = bool(np.any(core <= 0.5 * medw))
        if not sat:
            pk = st[half, half]
            sat = bool(pk * per_exp[s["obs_id"]] > sr["sat"])
            if not sat and pk < st.max():
                changes["dropped"] += 1
                continue
        if sat and s["g"] > sr["sat_g_max"]:
            changes["dropped"] += 1
            continue
        if sat != s["sat"]:
            changes["sat_to_unsat" if s["sat"] else "unsat_to_sat"] += 1
        out.append({**s, "sat": sat, "sat_v2": s["sat"]})
        changes["kept"] += 1
    pickle.dump(out, open(lib_out, "wb"))
    summary = {
        "label": "derived (v3 re-flag of the stored GO 17308 F814W library; per-exposure "
        "times from SCI NCOMBINE)",
        "per_exposure_s": sorted(set(per_exp.values())),
        "changes": changes,
        "n_unsaturated_v2": sum(not s["sat"] for s in lib),
        "n_unsaturated_v3": sum(not s["sat"] for s in out),
        "g_unsaturated_v3": [
            round(min(s["g"] for s in out if not s["sat"]), 2),
            round(max(s["g"] for s in out if not s["sat"]), 2),
        ],
        "stars": [
            {"gaia_dr3": s["source_id"], "g": s["g"], "sat": s["sat"], "sat_v2": s["sat_v2"]}
            for s in out
        ],
    }
    Path(summary_json).write_text(json.dumps(summary, indent=1) + "\n")
    print(
        summary["per_exposure_s"],
        changes,
        summary["n_unsaturated_v2"],
        summary["n_unsaturated_v3"],
        summary["g_unsaturated_v3"],
    )


if __name__ == "__main__":
    main(*sys.argv[1:4])
