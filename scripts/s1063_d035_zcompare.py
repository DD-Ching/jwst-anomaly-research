"""Post hoc (after the frozen D-035 run): CATS v4.1 arcs.txt redshifts and positions against the
Chan+2020 list for images matched within 1 arcsec. Writes zcompare_posthoc.json.

    python scripts/s1063_d035_zcompare.py
"""

import json
import sys

sys.path.insert(0, "scripts")
import lens_consistency as lc
import numpy as np
import s1063_d035 as d
from astropy.coordinates import SkyCoord

from jwst_anomaly import lensmodel

f = lc.model_files(d.MODEL, ("arcs.dat",))
cats = lensmodel.load_lenstool_images(f["arcs.dat"])
chan = d.chan_images()
a = SkyCoord(chan["ra"], chan["dec"], unit="deg")
b = SkyCoord(cats["ra"], cats["dec"], unit="deg")
idx, sep, _ = a.match_to_catalog_sky(b)
rows = []
for i in range(len(chan)):
    if sep[i].arcsec <= 1.0:
        rows.append(
            {
                "chan_id": chan["image_id"][i],
                "cats_id": str(cats["image_id"][idx[i]]),
                "sep_arcsec": round(float(sep[i].arcsec), 3),
                "chan_z": float(chan["z"][i]),
                "chan_z_spec": bool(chan["z_spec"][i]),
                "cats_z": float(cats["z"][idx[i]]),
            }
        )
s = [r for r in rows if r["chan_z_spec"]]
out = {
    "provenance": "derived (post hoc)",
    "n_matched": len(rows),
    "n_matched_zspec": len(s),
    "n_cats_z_zero": sum(r["cats_z"] == 0 for r in rows),
    "n_cats_images": len(cats),
    "n_cats_z_positive": int((np.asarray(cats["z"]) > 0).sum()),
    "rows": rows,
}
open("results/lens_consistency/abells1063_d035/zcompare_posthoc.json", "w").write(
    json.dumps(out, indent=1) + "\n"
)
for r in s:
    print(r)
print({k: v for k, v in out.items() if k != "rows"})
