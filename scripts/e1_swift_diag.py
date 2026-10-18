"""Post hoc diagnostics of the E1 Swift BAT positive-control channel (BAT-GBM), not part of the
frozen plan (``results/e1_swift/PLAN.md``). Each writes one JSON under ``results/e1_swift/``.

  python scripts/e1_swift_diag.py peryear|sunframe|lagcurve|nocounterparts
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import json  # noqa: E402

import astropy.units as u  # noqa: E402
import e1_swift as sw  # noqa: E402
import numpy as np  # noqa: E402
from astropy.coordinates import GeocentricTrueEcliptic, SkyCoord, get_sun  # noqa: E402
from astropy.table import Table  # noqa: E402
from astropy.time import Time  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402


def peryear():

    ev = sw.load(sw.fetch(False))
    s = {k: en.Sample.from_table(v) for k, v in ev.items()}
    p = sw.P

    def per_year(b, g):
        i, j = en.pairs_within(b.mjd, g.mjd, 7 * en.DAY, False)
        lag = np.abs(g.mjd[j] - b.mjd[i]) * en.DAY
        sep = en.sep_deg(b.ra[i], b.dec[i], g.ra[j], g.dec[j])
        cls = en.sep_class(sep, b.sigma[i], g.sigma[j], p)
        m = (lag >= 3600) & (cls == 1)
        return np.bincount(s["BAT"].year[i[m]] - 2004, minlength=24)  # anchor year (original order)

    obs = per_year(s["BAT"], s["GBM"])
    rng = np.random.default_rng(5)
    null = np.stack(
        [
            per_year(en.scramble_jit(s["BAT"], rng, p), en.scramble_jit(s["GBM"], rng, p))
            for _ in range(300)
        ]
    )
    yrs = np.arange(2004, 2028)

    # downtime proxy: longest gap between consecutive events per year, per catalogue (days)
    def maxgap(t, y):
        tt = np.sort(t.mjd[t.year == y])
        return float(np.max(np.diff(tt))) if len(tt) > 2 else float("nan")

    rows = []
    for k, y in enumerate(yrs):
        if null[:, k].mean() < 1:
            continue
        rows.append(
            {
                "year": int(y),
                "obs": int(obs[k]),
                "jit_mean": round(float(null[:, k].mean()), 1),
                "jit_sd": round(float(null[:, k].std()), 1),
                "bat_n": int((s["BAT"].year == y).sum()),
                "gbm_n": int((s["GBM"].year == y).sum()),
                "bat_max_gap_d": round(maxgap(s["BAT"], y), 1),
                "gbm_max_gap_d": round(maxgap(s["GBM"], y), 1),
            }
        )
    out = {
        "provenance": "derived",
        "note": "post hoc diagnostic; BAT-GBM wide pairs at 1 h-7 d by BAT anchor year "
        "(scrambled anchor years are the original ones: the per-year permutation keeps the "
        "year); 300 jit scrambles",
        "total_obs": int(obs.sum()),
        "total_jit_mean": round(float(null.sum(1).mean()), 1),
        "rows": rows,
    }
    with open(sw.OUT / "control_deficit_by_year.json", "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    for r in rows:
        print(r)
    print(out["total_obs"], out["total_jit_mean"])


def sunframe():

    ev = sw.load(sw.fetch(False))
    s = {k: en.Sample.from_table(v) for k, v in ev.items()}
    p = sw.P
    # ecliptic coordinates and Sun longitude at the true times, once
    ecl = {}
    for k in ("BAT", "GBM"):
        c = SkyCoord(s[k].ra * u.deg, s[k].dec * u.deg).transform_to(GeocentricTrueEcliptic())
        ecl[k] = (c.lon.deg, c.lat.deg)
    grid = np.arange(53000, 61500, 1.0)
    sun_lon = get_sun(Time(grid, format="mjd")).transform_to(GeocentricTrueEcliptic()).lon.deg
    sun_lon = np.rad2deg(np.unwrap(np.deg2rad(sun_lon)))

    def sunl(m):
        return np.interp(m, grid, sun_lon)

    def sun_scramble(k, rng):
        j = en.scramble_jit(s[k], rng, p)  # times only (its RA is the hour-angle-kept one)
        lon = np.mod(ecl[k][0] + sunl(j.mjd) - sunl(s[k].mjd), 360.0)
        c = SkyCoord(lon * u.deg, ecl[k][1] * u.deg, frame=GeocentricTrueEcliptic()).icrs
        return en.Sample(k, j.mjd, c.ra.deg, c.dec.deg, s[k].sigma, s[k].year)

    obs = en.count_channel(s["BAT"], s["GBM"], p, same=False)[:5, :2]
    rng = np.random.default_rng(9)
    null = np.stack(
        [
            en.count_channel(sun_scramble("BAT", rng), sun_scramble("GBM", rng), p, same=False)[
                :5, :2
            ]
            for _ in range(300)
        ]
    )
    lab = en.lag_labels(p)
    rows = [
        {
            "lag": lab[k],
            "cls": c,
            "obs": int(obs[k, j]),
            "sun_null_mean": round(float(null[:, k, j].mean()), 2),
            "sun_null_sd": round(float(null[:, k, j].std()), 2),
            "z": round(
                float((obs[k, j] - null[:, k, j].mean()) / max(null[:, k, j].std(), 1e-9)), 2
            ),
        }
        for k in range(5)
        for j, c in enumerate(("same", "wide"))
    ]
    out = {
        "provenance": "derived",
        "note": "post hoc diagnostic, not part of the frozen plan; 300 scrambles; "
        "per-year permutation + jit shifts as the jit null, then each BAT and GBM event keeps "
        "its ecliptic latitude and its ecliptic longitude relative to the Sun",
        "rows": rows,
    }
    with open(str(sw.OUT / "control_sun_frame_null.json"), "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    for r in rows:
        print(r)


def lagcurve():

    ev = sw.load(sw.fetch(False))
    s = {k: en.Sample.from_table(v) for k, v in ev.items()}
    p = sw.P
    edges = np.array([1 / 24, 1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8, 16, 30]) * en.DAY

    def h(b, g):
        i, j = en.pairs_within(b.mjd, g.mjd, edges[-1], False)
        return np.histogram(np.abs(g.mjd[j] - b.mjd[i]) * en.DAY, edges)[0]

    # also the signed version: GBM after BAT (+) vs before (-)
    def hs(b, g):
        i, j = en.pairs_within(b.mjd, g.mjd, edges[-1], False)
        d = (g.mjd[j] - b.mjd[i]) * en.DAY
        return np.histogram(d[d > 0], edges)[0], np.histogram(-d[d < 0], edges)[0]

    obs = h(s["BAT"], s["GBM"])
    oa, ob = hs(s["BAT"], s["GBM"])
    rng = np.random.default_rng(13)
    null = np.stack(
        [
            h(en.scramble_jit(s["BAT"], rng, p), en.scramble_jit(s["GBM"], rng, p))
            for _ in range(300)
        ]
    )
    rows = [
        {
            "lag_lo_d": round(float(edges[k] / en.DAY), 4),
            "lag_hi_d": round(float(edges[k + 1] / en.DAY), 4),
            "obs": int(obs[k]),
            "gbm_after_bat": int(oa[k]),
            "gbm_before_bat": int(ob[k]),
            "jit_mean": round(float(null[:, k].mean()), 1),
            "jit_sd": round(float(null[:, k].std()), 1),
            "ratio": round(float(obs[k] / null[:, k].mean()), 3),
        }
        for k in range(len(edges) - 1)
    ]
    out = {
        "provenance": "derived",
        "note": "post hoc diagnostic; all-separation BAT-GBM pairs; 300 jit scrambles",
        "rows": rows,
    }
    with open(str(sw.OUT / "control_lag_curve.json"), "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    for r in rows:
        print(r)


def nocounterparts():

    ev = sw.load(sw.fetch(False))
    cp = Table.read(str(sw.OUT / "same_dir_pairs_BAT-GBM.ecsv"))
    b = ev["BAT"][~np.isin(ev["BAT"]["name"], cp["a"])]
    g = ev["GBM"][~np.isin(ev["GBM"]["name"], cp["b"])]
    sb, sg = en.Sample.from_table(b), en.Sample.from_table(g)
    p = sw.P
    edges = np.array([1 / 24, 1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8, 16]) * en.DAY

    def h(x, y):
        i, j = en.pairs_within(x.mjd, y.mjd, edges[-1], False)
        return np.histogram(np.abs(y.mjd[j] - x.mjd[i]) * en.DAY, edges)[0]

    obs = h(sb, sg)
    rng = np.random.default_rng(17)
    null = np.stack(
        [h(en.scramble_jit(sb, rng, p), en.scramble_jit(sg, rng, p)) for _ in range(300)]
    )
    rows = [
        {
            "lag_lo_d": round(float(edges[k] / en.DAY), 4),
            "lag_hi_d": round(float(edges[k + 1] / en.DAY), 4),
            "obs": int(obs[k]),
            "jit_mean": round(float(null[:, k].mean()), 1),
            "jit_sd": round(float(null[:, k].std()), 1),
            "ratio": round(float(obs[k] / null[:, k].mean()), 3),
        }
        for k in range(len(edges) - 1)
    ]
    tot = {
        "obs_1h_4d": int(obs[:6].sum()),
        "jit_mean_1h_4d": round(float(null[:, :6].sum(1).mean()), 1),
        "jit_sd_1h_4d": round(float(null[:, :6].sum(1).std()), 1),
    }
    out = {
        "provenance": "derived",
        "note": "post hoc diagnostic; counterparts removed: "
        f"{len(ev['BAT']) - len(b)} BAT and {len(ev['GBM']) - len(g)} GBM events; "
        "all-separation pairs; 300 jit scrambles",
        "rows": rows,
        "sum": tot,
    }
    with open(str(sw.OUT / "control_lag_curve_no_counterparts.json"), "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(out["note"])
    [print(r) for r in rows]
    print(tot)


if __name__ == "__main__":
    {
        "peryear": peryear,
        "sunframe": sunframe,
        "lagcurve": lagcurve,
        "nocounterparts": nocounterparts,
    }[sys.argv[1]]()
