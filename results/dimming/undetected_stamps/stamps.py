"""Reason for every NaN-flux exposure of the NEXUS monitored pairs in e4-e7 and of flag d30295
(results/dimming/undetected_stamps/PLAN.md, step 1). Range reads of a 31 x 31 px SCI + DQ stamp.

Usage (repo root): stamps.py OUT_DIR  -> OUT_DIR/reasons.ecsv, OUT_DIR/stamps.png
"""

import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import persistence_check as pc  # noqa: E402

RUN = Path("results/dimming/persistence_run/nexus")
EPOCHS = (4, 5, 6, 7)
DQ_BITS = {0: "DO_NOT_USE", 1: "SATURATED", 2: "JUMP_DET", 3: "DROPOUT", 9: "DEAD", 10: "HOT",
           13: "LOW_QE", 31: "REFERENCE_PIXEL"}  # fmt: skip


def reason(sci: np.ndarray, x: float, y: float, nx: int, ny: int, r_px: float, x0: int, y0: int):
    """PLAN.md step 1 classes for one exposure. ``sci`` is the stamp whose (0, 0) is detector
    pixel (x0, y0); (x, y) are detector pixels; (nx, ny) the array size."""
    if not (r_px <= x <= nx - 1 - r_px and r_px <= y <= ny - 1 - r_px):
        return "edge"
    yy, xx = np.mgrid[0 : sci.shape[0], 0 : sci.shape[1]]
    r = np.hypot(xx + x0 - x, yy + y0 - y)
    if not np.isfinite(sci[r <= r_px]).all():
        return "aperture_no_data"
    if ((r >= 2 * r_px + 1) & (r <= 4 * r_px) & np.isfinite(sci)).sum() < 10:
        return "annulus"
    return "not_reproduced"


def _uri(name: str) -> str:
    return f"s3://stpubdata/jwst/public/jw{name[2:7]}/{name[:13]}/{name}"


def stamp(job):
    uid, epoch, verdict, name, x, y = job
    long = "LONG" in name.upper()
    r_px = 2.0 if long else 3.0
    xi, yi = int(round(x)), int(round(y))
    with pc._open(_uri(name)) as hdul:
        ny, nx = hdul["SCI"].header["NAXIS2"], hdul["SCI"].header["NAXIS1"]
        y0, y1 = max(0, yi - pc.HALF), min(ny, yi + pc.HALF + 1)
        x0, x1 = max(0, xi - pc.HALF), min(nx, xi + pc.HALF + 1)
        sci = np.asarray(hdul["SCI"].section[y0:y1, x0:x1], dtype=float)
        dq = np.asarray(hdul["DQ"].section[y0:y1, x0:x1]).astype(np.int64)
    why = reason(sci, x, y, nx, ny, r_px, x0, y0)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    ap = np.hypot(xx - x, yy - y) <= r_px
    bad = ap & ~np.isfinite(sci)
    bits = int(np.bitwise_or.reduce(dq[bad])) if bad.any() else 0
    return {
        "uid": uid,
        "epoch": epoch,
        "verdict": verdict,
        "file": name,
        "x": x,
        "y": y,
        "reason": why,
        "n_ap": int(ap.sum()),
        "n_ap_nan": int(bad.sum()),
        "dq_bits_nan": ",".join(v for b, v in DQ_BITS.items() if bits >> b & 1),
        "stamp_nan_frac": float(np.mean(~np.isfinite(sci))),
    }, sci


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    v = Table.read(RUN / "verdicts.ecsv")
    jobs = []
    for k in (2, *EPOCHS):
        e = Table.read(RUN / f"exposures_e{k}.ecsv")
        vk = v[v["epoch"] == k]
        vk = vk[vk["uid"] == "nexus-d30295"] if k == 2 else vk[vk["kind"] == "monitored"]
        verd = dict(zip(vk["uid"], vk["verdict"], strict=True))
        for r in e[np.isin(e["uid"], list(verd)) & ~np.isfinite(np.asarray(e["flux"], float))]:
            jobs.append((r["uid"], k, verd[r["uid"]], r["file"], float(r["x"]), float(r["y"])))
    with ThreadPoolExecutor(12) as ex:
        res = list(ex.map(stamp, jobs))
    t = Table(rows=[r for r, _ in res])
    t.meta.update(provenance="derived", source="results/dimming/undetected_stamps/stamps.py")
    t.write(out / "reasons.ecsv", overwrite=True)
    sheet(t, [s for _, s in res], out / "stamps.png")
    print(dict(zip(*np.unique(t["reason"], return_counts=True), strict=True)))


def sheet(t: Table, stamps: list, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sel = [i for i, r in enumerate(t) if r["verdict"] == "undetected"]
    n = len(sel)
    cols = 6
    fig, axs = plt.subplots(-(-n // cols), cols, figsize=(2 * cols, 2.2 * -(-n // cols)))
    for ax in np.ravel(axs):
        ax.axis("off")
    for ax, i in zip(np.ravel(axs), sel, strict=False):
        s = stamps[i]
        fin = s[np.isfinite(s)]
        lo, hi = np.percentile(fin, [5, 99.5]) if fin.size else (0, 1)
        ax.imshow(np.ma.masked_invalid(s), origin="lower", vmin=lo, vmax=hi, cmap="gray")
        r = t[i]
        ax.set_title(f"{r['uid'][6:]} e{r['epoch']} {r['file'][26:32]}\n{r['reason']}", fontsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=80)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
