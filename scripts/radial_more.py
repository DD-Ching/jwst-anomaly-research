"""Run the ``radial`` screen (D-031, W1 signature) on fields staged for the W1 shear screen.

Each field uses ``inject_shear.field_spec(<field>)`` (model, F200W catalogue, photo-z, Gaia veto
centre), loaded by ``inject_radial.load_field`` (CANUCS photo-z converted to an eazy-like zout,
Gaia DR3 spike stars queried at the centre), then ``exotic_screens.cmd_radial`` with the CLI
defaults and the spec's ``max_radius``. Results go to ``OUT/<field>/<model>/radial.json`` and the
compact table ``OUT/summary.json``. Plan: results/radial_more/PLAN.md.

    python scripts/radial_more.py --fields macs0417 macs1423 --out results/radial_more

``--sheet`` draws, per field, the F200W stamp (S3 byte range) around the strongest centre with
the arcs whose lines pass it (``OUT/sheet_radial_top.png``). ``--reproduce-s1063`` runs Abell
S1063 with ``inject_radial.FIELDS`` (no photo-z, no Gaia stars), the headline run of
docs/fields/abells1063.md.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import exotic_screens as es  # noqa: E402
import inject_shear as ish  # noqa: E402

ir = ish.ir
KEYS = (
    "n_elongated",
    "n_spike_segments_dropped",
    "n_not_background_dropped",
    "n_anti",
    "n_anti_not_model_radial",
    "n_peaks",
    "n_dark_centre_peaks",
    "random_peaks_mean",
    "random_peaks_p95",
    "max_lines",
    "p_random_max",
)
P_FLAG = 0.05  # a centre is flagged when p_random < 0.05 (the screen's rule, D-031; ASSUMPTION)


def flagged(summary: dict) -> list[dict]:
    """Centres the frozen rule flags: p_random < P_FLAG."""
    return [p for p in summary["peaks"] if p.get("p_random") is not None and p["p_random"] < P_FLAG]


def run(name: str, spec: dict, out: Path, tag: str) -> dict:
    fdir = out / tag
    fdir.mkdir(parents=True, exist_ok=True)
    _, _, extra, cat, photoz = ir.load_field(name, fdir, spec=spec)
    args = SimpleNamespace(**es.radial_defaults())
    args.max_radius = spec["max_radius"]
    args.model, args.catalog, args.photoz, args.out = spec["model"], cat, photoz, fdir
    args.spike_stars = fdir / f"{name}_gaia_stars.ecsv" if extra is not None else None
    s = es.cmd_radial(args)
    row = {"field": tag, "model": spec["model"], **{k: s[k] for k in KEYS}}
    row["n_spike_stars"] = None if s["spike_stars"] is None else s["spike_stars"]["n"]
    row["flags"] = flagged(s)
    return row


def contributors(cand, model, centre, args) -> np.ndarray:
    """Arcs whose major-axis line passes the centre under ``exotic_screens.line_counts``."""
    x, y = model.to_frame(cand["ra"], cand["dec"])
    gx, gy = np.array([[centre["x"]]]), np.array([[centre["y"]]])
    hit = []
    for xi, yi, pa in zip(np.asarray(x), np.asarray(y), np.asarray(cand["pa_obs"]), strict=True):
        hit.append(
            bool(
                es.line_counts([xi], [yi], [pa], gx, gy, args.line_tol_arcsec, args.max_len_arcsec)[
                    0, 0
                ]
            )
        )
    return np.array(hit, bool)


def sheet(fields: list[str], out: Path, size: float = 36.0) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.table import Table
    from astropy.wcs import WCS
    from w1_streak_mask import i2d_uri

    from jwst_anomaly.cutouts import make_cutouts

    lc = es.lc
    args = SimpleNamespace(**es.radial_defaults())
    fig, axes = plt.subplots(1, len(fields), figsize=(4.2 * len(fields), 4.6), squeeze=False)
    for ax, tag in zip(axes[0], fields, strict=True):
        name = tag.removesuffix("_reproduce")
        spec = ir.FIELDS[name] if tag.endswith("_reproduce") else ish.field_spec(name)
        d = out / tag / spec["model"]
        summ = json.loads((d / "radial.json").read_text())
        if not summ["peaks"]:  # no centre: an empty panel
            ax.set_title(f"{tag}: no centre", fontsize=9), ax.set_axis_off()
            continue
        top = max(summ["peaks"], key=lambda p: p["n_lines"])
        cand = Table.read(d / "radial_sources.ecsv")
        model, _, _ = lc.load_model(spec["model"])
        lc.apply_frame_offset(spec["model"], model)
        use = contributors(cand, model, top, args)
        tgt = Table({"source_uid": [f"{tag}_top"], "ra": [top["ra"]], "dec": [top["dec"]]})
        cdir = out / "stamps"
        row = make_cutouts(i2d_uri(spec["catalog"]), tgt, size_arcsec=size, out_dir=cdir)[0]
        with fits.open(row["path"]) as h:
            img, w = h["SCI"].data, WCS(h["SCI"].header)
        lo, hi = np.nanpercentile(img, [5, 99.5])
        ax.imshow(np.arcsinh((img - lo) / (hi - lo) * 10), origin="lower", cmap="gray_r")
        # each arc's major axis, 1.5" either side, drawn through the stamp WCS (mosaics are rotated)
        ra0, dec0 = np.asarray(cand["ra"][use], float), np.asarray(cand["dec"][use], float)
        t = np.deg2rad(np.asarray(cand["pa_obs"][use], float))  # E of N
        half = 1.5 / 3600.0
        dra = half * np.sin(t) / np.cos(np.deg2rad(dec0))
        x1, y1 = w.world_to_pixel_values(ra0 + dra, dec0 + half * np.cos(t))
        x2, y2 = w.world_to_pixel_values(ra0 - dra, dec0 - half * np.cos(t))
        for a, b, c, e in zip(x1, y1, x2, y2, strict=True):
            ax.plot([a, c], [b, e], color="tab:orange", lw=1.5)
        cx, cy = w.world_to_pixel_values(top["ra"], top["dec"])
        ax.plot(cx, cy, "+", color="tab:blue", ms=14, mew=1.5)
        ax.set_xlim(0, img.shape[1]), ax.set_ylim(0, img.shape[0])
        pr = top["p_random"]
        ptxt = "n/a" if pr is None else f"{pr:.3f}"
        ax.set_title(
            f"{tag}: {top['n_lines']} lines, p {ptxt}\n{sum(use)} arcs; dark {top['dark_centre']}",
            fontsize=9,
        )
        ax.set_xticks([]), ax.set_yticks([])
    fig.tight_layout()
    path = out / "sheet_radial_top.png"
    fig.savefig(path, dpi=80)
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fields", nargs="*", default=[], choices=ish.FIELD_NAMES)
    ap.add_argument("--reproduce-s1063", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--sheet", nargs="*", help="field tags to draw (after their runs)")
    a = ap.parse_args(argv)
    if a.sheet:
        print(sheet(a.sheet, a.out))
        return 0
    rows = []
    if a.reproduce_s1063:
        rows.append(run("abells1063", ir.FIELDS["abells1063"], a.out, "abells1063_reproduce"))
    rows += [run(f, ish.field_spec(f), a.out, f) for f in a.fields]
    path = a.out / "summary.json"
    old = json.loads(path.read_text()) if path.exists() else {}
    old.update({r["field"]: r for r in rows})
    path.write_text(json.dumps(old, indent=1) + "\n")
    for r in rows:
        print(
            json.dumps(
                {
                    k: r[k]
                    for k in (
                        "field",
                        "n_anti_not_model_radial",
                        "max_lines",
                        "p_random_max",
                        "n_peaks",
                    )
                }
            ),
            len(r["flags"]),
            "flags",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
