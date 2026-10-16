"""Contact sheet for ``results/arcs/rxj2129-dja``: F200W (program 2767) stamps of the Lenstool run's
strong-shear background anti rows and four aligned rows, with the observed axis (red, catalogue PA)
and the predicted z = 2 stretch axis (cyan) drawn through the WCS. Stamps are S3 range reads,
written to ``--work`` and not stored.

    python scripts/arcs_rxj2129_dja_sheet.py --work <scratch dir>
"""

from __future__ import annotations

import argparse
from pathlib import Path

import astropy.units as u
import matplotlib
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from jwst_anomaly.cutouts import make_cutouts

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

I2D = "s3://stpubdata/jwst/public/jw02767/L3/t/o002/jw02767-o002_t001_nircam_clear-f200w_i2d.fits"
OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "rxj2129-dja"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)
    t = Table.read(OUT / "rxj2129-lenstool-maps" / "arcs.ecsv")
    b = t[t["strong_shear"] & (t["population"] == "background")]
    rows = list(b[b["orientation_class"] == "anti"]) + list(
        b[b["orientation_class"] == "aligned"][:4]
    )
    targets = Table(
        {
            "source_uid": [f"{r['orientation_class'][:4]}_{r['label']}" for r in rows],
            "ra": [r["ra"] for r in rows],
            "dec": [r["dec"] for r in rows],
        }
    )
    cuts = make_cutouts(I2D, targets, size_arcsec=3.0, out_dir=args.work)
    fig, axs = plt.subplots(2, 7, figsize=(21, 6.5))
    for ax, r, cr in zip(axs.flat[: len(rows)], rows, cuts, strict=True):
        ax.set_axis_off()
        if not cr["path"]:
            continue
        d = fits.getdata(cr["path"], "SCI")
        lo, hi = np.nanpercentile(d, [5, 99.5])
        ax.imshow(d, origin="lower", cmap="gray_r", vmin=lo, vmax=hi)
        w = WCS(fits.getheader(cr["path"], "SCI"))
        c0 = SkyCoord(r["ra"], r["dec"], unit="deg")
        x0, y0 = (float(v) for v in w.world_to_pixel(c0))
        for pa, col in ((r["pa_obs"], "r"), (r["pa_pred_z2"], "c")):
            x1, y1 = (
                float(v)
                for v in w.world_to_pixel(c0.directional_offset_by(pa * u.deg, 0.9 * u.arcsec))
            )
            ax.plot([2 * x0 - x1, x1], [2 * y0 - y1, y1], col, lw=1)
        ax.set_title(
            f"{r['orientation_class']} {r['label']} z{r['z_phot']:.1f} g{r['g_min']:.2f}",
            fontsize=8,
        )
    fig.savefig(OUT / "sheet_anti_aligned.png", dpi=70, bbox_inches="tight")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
