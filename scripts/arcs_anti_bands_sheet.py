"""Contact sheet of the 3 single-component anti rows in F150W, F200W and F444W / F410M (by eye).

Re-reads 8″ S3 cutouts of the same observation (temporary) and ``summary.json`` of
``scripts/arcs_anti_bands.py``; the red bar is the F200W catalogue axis (sky PA, north up):

    python scripts/arcs_anti_bands_sheet.py

Output: ``results/arcs/anti-bands/anti_bands.png`` (pixels observed; titles derived).
"""

from __future__ import annotations

import json
import sys
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
from scipy.ndimage import median_filter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arcs_anti_bands as aab  # noqa: E402

from jwst_anomaly import cutouts, paths  # noqa: E402

OUT = paths.repo_root() / "results" / "arcs" / "anti-bands"


def image_url(obs: str, band: str) -> str:
    """S3 URI of the i2d mosaic (same layout as ``arcs_bands.s3_url``)."""
    https = aab.ab.s3_url(f"{obs}_nircam_clear-{band}_i2d.fits")
    return https.replace(aab.ab.S3_HTTPS, "s3://stpubdata/")


def main() -> int:
    summ = json.loads((OUT / "summary.json").read_text())
    fig, ax = plt.subplots(len(aab.ROWS), 3, figsize=(9.6, 3.4 * len(aab.ROWS)), squeeze=False)
    with tempfile.TemporaryDirectory() as tmp:
        for i, (label, field, model, obs, sw, lw) in enumerate(aab.ROWS):
            stored = Table.read(
                paths.repo_root() / "results" / "arcs" / field / model / "arcs.ecsv"
            )
            r = stored[np.asarray(stored["label"]) == label][0]
            info = summ["rows"][str(label)]
            for j, (band, key) in enumerate(((sw, "SW"), ("f200w", None), (lw, "LW"))):
                axis = ax[i][j]
                tg = Table({"source_uid": [str(label)], "ra": [r["ra"]], "dec": [r["dec"]]})
                c = cutouts.make_cutouts(
                    image_url(obs, band), tg, size_arcsec=8.0, out_dir=Path(tmp) / f"{label}{band}"
                )
                with fits.open(c["path"][0]) as h:
                    ext = next(x for x in h if x.data is not None)
                    img, w = ext.data.astype(float), WCS(ext.header)
                x, y = (round(float(v)) for v in w.world_to_pixel_values(r["ra"], r["dec"]))
                k = round(2.0 / (proj_plane_pixel_scales(w.celestial)[0] * 3600))  # 2″ half-size
                s = img[max(y - k, 0) : y + k + 1, max(x - k, 0) : x + k + 1]
                s = s - median_filter(np.nan_to_num(s), 2 * k + 1)  # 4″ window
                lo, hi = np.nanpercentile(s, [5, 99.7])
                axis.imshow(s, origin="lower", cmap="gray", vmin=lo, vmax=hi)
                # F200W axis as a 1.2" bar, placed through the WCS (sky PA east of north)
                pa, half = np.deg2rad(float(r["pa_obs"])), 0.6 / 3600
                dra = half * np.sin(pa) / np.cos(np.deg2rad(float(r["dec"])))
                ends = [
                    (r["ra"] - dra, r["dec"] - half * np.cos(pa)),
                    (r["ra"] + dra, r["dec"] + half * np.cos(pa)),
                ]
                px = [w.world_to_pixel_values(a, d) for a, d in ends]
                ox, oy = max(x - k, 0), max(y - k, 0)
                axis.plot(
                    [q[0] - ox for q in px], [q[1] - oy for q in px], color="r", lw=1, alpha=0.6
                )
                if key is None:
                    t = f"{field} {label} F200W\npa {r['pa_obs']:.1f} e {r['ellipticity']:.2f}"
                else:
                    b = info["bands"][key]
                    t = f"{band.upper()} {b['verdict']}\ndpa {b['dpa']} e {b['e']} S/N {b['snr']}"
                axis.set_title(t, fontsize=8)
                axis.set_xticks([])
                axis.set_yticks([])
    fig.tight_layout()
    fig.savefig(OUT / "anti_bands.png", dpi=80)
    return 0


if __name__ == "__main__":
    sys.exit(main())
