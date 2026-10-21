"""Post hoc (after the frozen D-035 run): CATS v4.1 ``arcs.txt`` redshifts and positions against
the Chan+2020 list, one-to-one matches within ``Params.match_arcsec``. Writes
``zcompare_posthoc.json`` next to the run's ``summary.json``.

    python scripts/s1063_d035_zcompare.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402
import s1063_d035 as d  # noqa: E402

from jwst_anomaly import lensmodel, schema  # noqa: E402


def compare(chan, cats, max_sep: float) -> dict:
    idx, sep, m = d.match_cats(chan, cats, max_sep)
    rows = [
        {
            "chan_id": str(chan["image_id"][i]),
            "cats_id": str(cats["image_id"][idx[i]]),
            "sep_arcsec": round(float(sep[i]), 3),
            "chan_z": float(chan["z"][i]),
            "chan_z_spec": bool(chan["z_spec"][i]),
            "cats_z": float(cats["z"][idx[i]]),
        }
        for i in np.where(m)[0]
    ]
    spec = [r for r in rows if r["chan_z_spec"]]
    dz = [abs(r["chan_z"] - r["cats_z"]) for r in spec]
    return {
        "provenance": schema.Provenance.DERIVED.value,
        "post_hoc": True,
        "match_arcsec": max_sep,
        "n_matched": len(rows),
        "n_matched_zspec": len(spec),
        "n_systems_matched_zspec": len({lensmodel.image_family(r["chan_id"]) for r in spec}),
        "median_sep_arcsec_all": float(np.median([r["sep_arcsec"] for r in rows])),
        "median_sep_arcsec_zspec": float(np.median([r["sep_arcsec"] for r in spec])),
        "max_abs_dz_zspec": float(max(dz)),
        "n_cats_images": len(cats),
        "n_cats_z_positive": int((np.asarray(cats["z"]) > 0).sum()),
        "rows": rows,
    }


def main() -> int:
    f = lc.model_files(d.MODEL, ("arcs.dat",))
    out = compare(
        d.chan_images(), lensmodel.load_lenstool_images(f["arcs.dat"]), d.Params().match_arcsec
    )
    with open(d.OUT / "zcompare_posthoc.json", "w", encoding="utf-8") as fh:
        fh.write(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
