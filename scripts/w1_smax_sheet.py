"""F200W cutouts around W1 shear-screen maxima, for looking at them (D-053).

    python scripts/w1_smax_sheet.py --out sheet.png \
        FIELD=SUMMARY_JSON=I2D_URI[=GAIA_ECSV] [FIELD=...]

``SUMMARY_JSON`` is an ``inject_shear.py`` per-field summary (its ``screen.s_max_radec``);
``I2D_URI`` a level-3 ``_i2d.fits`` (S3 byte-range read). Rings mark 1″ and the 10″ aperture;
cyan crosses are Gaia stars of the spike veto. Display only (``observed`` pixels).
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.wcs import WCS  # noqa: E402
from astropy.wcs.utils import proj_plane_pixel_scales  # noqa: E402

from jwst_anomaly import cutouts  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("fields", nargs="+", help="FIELD=SUMMARY_JSON=I2D_URI[=GAIA_ECSV]")
    ap.add_argument("--size-arcsec", type=float, default=24.0)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    fig, axs = plt.subplots(1, len(args.fields), figsize=(5 * len(args.fields), 5.4), squeeze=False)
    for ax, spec in zip(axs[0], args.fields, strict=True):
        name, summary, uri, *gaia = spec.split("=")
        s = json.loads(Path(summary).read_text())["screen"]
        ra, dec = s["s_max_radec"]
        tg = Table({"source_uid": [f"{name}_smax"], "ra": [ra], "dec": [dec]})
        with tempfile.TemporaryDirectory() as tmp:
            cut = cutouts.make_cutouts(
                uri,
                tg,
                size_arcsec=args.size_arcsec,
                out_dir=Path(tmp),
                storage_options={"anon": True},
            )
            with fits.open(cut["path"][0]) as h:
                img = h["SCI"].data.astype(float)
                w = WCS(h["SCI"].header)
        lo, hi = np.nanpercentile(img, [5, 99.5])
        ax.imshow(np.arcsinh((img - lo) / (hi - lo) * 10), origin="lower", cmap="gray_r")
        px = 1 / (proj_plane_pixel_scales(w)[0] * 3600)
        cy, cx = (np.array(img.shape) - 1) / 2
        for r in (1, 10):
            ax.add_patch(plt.Circle((cx, cy), r * px, fill=False, color="tab:red", lw=1))
        if gaia:
            g = Table.read(gaia[0])
            x, y = w.world_to_pixel_values(g["ra"], g["dec"])
            ok = (x > 0) & (y > 0) & (x < img.shape[1]) & (y < img.shape[0])
            ax.plot(x[ok], y[ok], "c+", ms=12)
        ax.set_title(
            f"{name} S_max {s['s_max']:.2f} (p_rot {s['p_random_max']:.2f})\n"
            f"({ra:.5f}, {dec:.5f}); F200W {args.size_arcsec:g}″; rings 1″/10″",
            fontsize=9,
        )
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(args.out, dpi=90)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
