"""Post hoc: c1d continuum level vs the 4043 F444W imaging flux, for every extracted source."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def ratios(c1d: Path, cat: Path, band=(4.2, 4.6)) -> dict:
    from astropy.io import fits
    from astropy.table import Table

    t = Table.read(cat)
    flux_by_label = dict(
        zip(np.asarray(t["label"]), np.asarray(t["aper_total_flux"], dtype=float), strict=True)
    )
    with fits.open(c1d) as h:
        d = h[1].data
        out = []
        for row in d:
            w, f = np.asarray(row["WAVELENGTH"], float), np.asarray(row["FLUX"], float)
            ok = (
                np.isfinite(w)
                & np.isfinite(f)
                & (np.asarray(row["DQ"]) == 0)
                & (np.asarray(row["N_INPUT"]) > 0)
                & (w > band[0])
                & (w < band[1])
            )
            img = flux_by_label.get(int(row["SOURCE_ID"]), np.nan)
            if ok.sum() > 10 and np.isfinite(img) and img > 0:
                out.append((int(row["SOURCE_ID"]), float(np.median(f[ok]) / img)))
    r = np.array([v for _, v in out])
    return {
        "file": c1d.name,
        "band_um": list(band),
        "percentiles": [16, 84],
        "n": len(r),
        "ratio_median": round(float(np.median(r)), 3),
        "ratio_p16_p84": [round(float(np.percentile(r, p)), 3) for p in (16, 84)],
        "per_source": {str(k): round(v, 3) for k, v in out},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--c1d", nargs="+", type=Path, required=True)
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--source-id", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    res = {}
    for p in args.c1d:
        r = ratios(p, args.cat)
        r["source"] = r["per_source"].get(str(args.source_id))
        vals = np.array(list(r["per_source"].values()))
        r["source_percentile"] = (
            None if r["source"] is None else round(float(np.mean(vals < r["source"]) * 100), 1)
        )
        del r["per_source"]
        res[p.stem] = r
    args.out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
