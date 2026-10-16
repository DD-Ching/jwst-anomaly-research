"""Post hoc (SMACS 0723 arc_diff RESULT.md): h061's F444W z vs controls by distance to the star.

One-off post hoc analysis for h061 (hard-coded positions, obs_ids and the sibling ``h061_spike``
directory); writes ``near_star.json`` and ``along_line.json`` into the given directory.

    python scripts/h061_near_star.py results/transients/smacs_arcdiff/h061_ctrl
"""

import json
import sys

import numpy as np
from astropy.table import Table

D = sys.argv[1]
ra0, dec0 = 110.82204720838848, -73.44918269902084
star = (110.82627304494622, -73.44867513014545)  # smoothed F444W peak, 2022 ERO mosaic (derived)


def sep(p):
    return float(np.hypot((p[0] - star[0]) * np.cos(np.radians(p[1])), p[1] - star[1]) * 3600)


t = Table.read(f"{D}/lightcurves.ecsv")
pos = {
    ln.split()[0]: (float(ln.split()[1]), float(ln.split()[2]))
    for ln in open(f"{D}/pick.txt").read().split("--pos ")[1:]
}
pos["h061"] = (ra0, dec0)
out = {
    "star_radec": star,
    "star_how": "maximum of the 3-px Gaussian-smoothed 2022 ERO F444W box (12 arcsec)",
    "pairs": {},
}
for _band, o1, o2 in [
    ("F444W", "jw02736-o001_t001_nircam_clear-f444w", "jw04043-o001_t001_nircam_clear-f444w"),
    ("F444W", "jw02736-o001_t001_nircam_clear-f444w", "jw06882-o057_t057_nircam_clear-f444w"),
]:
    z = {}
    for n in pos:
        a = t[(t["name"] == n) & (t["obs_id"] == o1)]
        b = t[(t["name"] == n) & (t["obs_id"] == o2)]
        if len(a) and len(b):
            z[n] = float(
                (b["flux_njy"][0] - a["flux_njy"][0]) / np.hypot(a["err_njy"][0], b["err_njy"][0])
            )
    zh = z["h061"]
    rows = {}
    for lab, lo, hi in [("near_le7", 0, 7.0), ("far_gt7", 7.0, 1e9)]:
        # c57 lies 0.03 arcsec from h061 (the same knot) and is excluded
        zs = [v for n, v in z.items() if n not in ("h061", "c57") and lo < sep(pos[n]) <= hi]
        rows[lab] = {
            "n": len(zs),
            "median_z": float(np.median(zs)),
            "robust_std_z": float(1.4826 * np.median(np.abs(np.array(zs) - np.median(zs)))),
            "frac_abs_z_ge_h061": float(np.mean(np.abs(zs) >= abs(zh))),
            "max_abs_z": float(np.max(np.abs(zs))),
        }
    out["pairs"][f"{o1} -> {o2}"] = {"h061_z": zh, "h061_star_sep_arcsec": sep((ra0, dec0)), **rows}
with open(f"{D}/near_star.json", "w") as fh:
    fh.write(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=1))


# along the 2024 stripe line (stripe_phot.py line_profile.ecsv): 2024 - 2022 F444W at 0.4-2.0 arcsec
lp = Table.read(f"{D}/../h061_spike/line_profile.ecsv")
x24 = lp[lp["program"] == "04043"]
x22 = lp[lp["program"] == "02736"]
a = np.asarray(x24["along_arcsec"])
dlt = np.asarray(x24["flux_njy"]) - np.asarray(x22["flux_njy"])
m = (np.abs(a) >= 0.4) & (np.abs(a) <= 2.0) & np.isfinite(dlt)
med = float(np.median(dlt[m]))
line = {
    "what": "2024 GO-4043 minus 2022 ERO F444W raw aperture flux along the stripe, 0.4-2.0 arcsec",
    "n": int(m.sum()),
    "median_njy": med,
    "robust_std_njy": float(1.4826 * np.median(np.abs(dlt[m] - med))),
    "max_abs_njy": float(np.max(np.abs(dlt[m]))),
    "frac_abs_ge_h061": float(np.mean(np.abs(dlt[m]) >= abs(float(dlt[a == 0][0])))),
    "h061_delta_njy": float(dlt[a == 0][0]),
}
with open(f"{D}/along_line.json", "w") as fh:
    fh.write(json.dumps(line, indent=2) + "\n")
print(json.dumps(line, indent=1))
