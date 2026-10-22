"""MACS1423 CANUCS reference run: post hoc checks (not frozen rules).

For the predicted extra images and the `absent` controls: the 0.3" aperture sum at the position on
the F200W BCG-subtracted mosaic, the source whose 3-sigma ellipse contains the position (rule 1),
its centroid separation, flux and ellipticity, and the separation to the nearest catalogued image.
Writes ``posthoc.json``; ``--sheet`` also writes a contact sheet PNG (not committed).

    python scripts/m1423_canucsref_posthoc.py --cat <photometry-cat> --image <f200w bgsub-sci>
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import m0417_arcref as ar
import m0417_canucsref as cr
import m0417_extra as mx
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from m0417_canucsref_posthoc import _clean, aperture_sum
from m1423_canucsref import AP_ARCSEC, FIELD, OUT, PIX

STAMP_ARCSEC = 3.0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--image", type=Path, required=True)
    ap.add_argument("--sheet", type=Path, help="optional contact sheet PNG path")
    args = ap.parse_args(argv)
    cat = cr.to_cat(Table.read(args.cat))
    with fits.open(args.image) as hl:
        w, d = WCS(hl[0].header), np.asarray(hl[0].data, float)
    pred = Table.read(OUT / f"predicted_{FIELD}.ecsv")
    ctrl = Table.read(OUT / f"control_{FIELD}.ecsv")
    imgs = mx.load_images(FIELD)
    ic = SkyCoord(np.asarray(imgs["ra"], float), np.asarray(imgs["dec"], float), unit="deg")
    cs = SkyCoord(cat["ra"], cat["dec"], unit="deg")
    items = [(f"{r['system']} {r['image_class']}", r) for r in pred]
    items += [(f"{r['system']} {r['images']}", r) for r in ctrl if r["class_local"] == "absent"]
    rows = []
    for name, r in items:
        pos = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
        x, y = (float(v) for v in w.world_to_pixel_values(float(r["ra"]), float(r["dec"])))
        k = ar.counterpart(float(r["ra"]), float(r["dec"]), cat)
        si = ic.separation(pos).arcsec
        j = int(np.argmin(si))
        rows.append(
            {
                "item": name,
                "ra": float(r["ra"]),
                "dec": float(r["dec"]),
                "f_pred_njy": float(r["f_pred"]) * 1e9,
                "sigma_loc_njy": float(r["sigma_loc_njy"]),
                "ap03_njy": aperture_sum(d, x, y, AP_ARCSEC / PIX),
                "containing_source": int(cat["label"][k]) if k >= 0 else None,
                "containing_sep_arcsec": float(cs[k].separation(pos).arcsec) if k >= 0 else None,
                "containing_flux_njy": float(cat["flux"][k]) * 1e9 if k >= 0 else None,
                "containing_flux_over_pred": (
                    float(cat["flux"][k] / r["f_pred"]) if k >= 0 else None
                ),
                "containing_ell": float(cat["ell"][k]) if k >= 0 else None,
                "search_arcsec": float(r["search_arcsec"]),
                "nearest_catalogued_image": str(imgs["image_id"][j]),
                "nearest_catalogued_sep_arcsec": float(si[j]),
            }
        )
    res = _clean(
        {
            "label": "post hoc (not frozen rules); observed pixels and catalogue, derived sums",
            "aperture_radius_arcsec": AP_ARCSEC,
            "items": rows,
        }
    )
    (OUT / "posthoc.json").write_text(json.dumps(res, indent=1, allow_nan=False) + "\n")
    print(json.dumps(res, indent=1))
    if args.sheet:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        h = int(STAMP_ARCSEC / PIX)
        fig, axs = plt.subplots(1, len(rows), figsize=(3.2 * len(rows), 3.4))
        for a, q in zip(np.atleast_1d(axs), rows, strict=True):
            x, y = (int(round(float(v))) for v in w.world_to_pixel_values(q["ra"], q["dec"]))
            s = d[y - h : y + h + 1, x - h : x + h + 1]
            v = np.nanpercentile(s, [5, 99.5])
            a.imshow(s, origin="lower", cmap="gray_r", vmin=v[0], vmax=v[1])
            a.add_patch(plt.Circle((h, h), AP_ARCSEC / PIX, fill=False, color="r"))
            a.add_patch(plt.Circle((h, h), q["search_arcsec"] / PIX, fill=False, color="c"))
            a.set_title(q["item"], fontsize=9)
            a.set_xticks([])
            a.set_yticks([])
        fig.tight_layout()
        fig.savefig(args.sheet, dpi=90)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
