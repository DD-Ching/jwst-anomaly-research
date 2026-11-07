"""M2 image-embedding pilot: Zoobot / DINOv2 vs a pixel baseline on SMACS F200W arc injections.

Plan: ``results/m2_embeddings/PLAN.md`` (frozen before any embedding was computed). The F200W
``_i2d.fits`` SCI array is streamed into memory from ``s3://stpubdata`` (never written to disk);
hosts come from the matching pipeline catalog. ``decide`` encodes the plan's outcome rule.

    python scripts/m2_embeddings.py --out results/m2_embeddings
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
from astropy.stats import sigma_clipped_stats

S3 = "s3://stpubdata/jwst/public/jw02736/L3/t/o001/jw02736-o001_t001_nircam_clear-f200w_{kind}"
# Frozen plan settings (ASSUMPTIONs).
SIZE = 64
EDGE = 40
MAG_RANGE = (20.0, 27.0)  # amended before any embedding (PLAN.md deviation)
N_CLEAN = 1000
N_INJ = 50
HOST_SEED = 0
INJ_SEEDS = (1, 2, 3)
STRENGTHS = (0.3, 1.0)
R_RANGE_ARCSEC = (0.35, 0.70)
EXTENT_DEG = (90.0, 180.0)
SIGMA_PX = 1.5
HOST_R_ARCSEC = 0.3
PIX_ARCSEC = 0.0312  # F200W i2d pixel scale (header CDELT; checked at run time)
ASINH_SOFT = 0.1
NORM_PCTL = 99.5
N_PCA = 50
N_TREES = 300
TOP_K = 50
MARGIN = 0.10
ENCODERS = {
    "zoobot": "hf_hub:mwalmsley/zoobot-encoder-convnext_nano",
    "dinov2": "vit_small_patch14_dinov2.lvd142m",
}


def decide(enc: list[float], base: list[float], margin: float = MARGIN) -> str:
    """Frozen per-encoder, per-strength outcome: ``beats`` / ``worse`` / ``tie``."""
    if len(enc) != len(base) or not enc:
        raise ValueError("need equal, non-empty per-seed recall lists")
    me, mb = float(np.mean(enc)), float(np.mean(base))
    tol = 1e-9  # recalls are multiples of 1/50; compare margins without float round-off
    if me >= mb + margin - tol and all(e >= b for e, b in zip(enc, base, strict=True)):
        return "beats"
    if me <= mb - margin + tol:
        return "worse"
    return "tie"


def arc_model(
    rng: np.random.Generator, size: int, pix_arcsec: float, total: float
) -> tuple[np.ndarray, dict[str, float]]:
    """Noiseless arc (ring segment, Gaussian radial profile) of total flux ``total``."""
    r = rng.uniform(*R_RANGE_ARCSEC) / pix_arcsec
    ext = np.deg2rad(rng.uniform(*EXTENT_DEG))
    phi0 = rng.uniform(0, 2 * np.pi)
    c = (size - 1) / 2
    yy, xx = np.mgrid[:size, :size]
    rr = np.hypot(yy - c, xx - c)
    dphi = np.angle(np.exp(1j * (np.arctan2(yy - c, xx - c) - phi0)))
    img = np.exp(-0.5 * ((rr - r) / SIGMA_PX) ** 2) * (np.abs(dphi) <= ext / 2)
    s = img.sum()
    img = img * (total / s) if s > 0 else img
    return img, {"r_px": float(r), "extent_deg": float(np.rad2deg(ext)), "phi0": float(phi0)}


def background(cut: np.ndarray) -> float:
    return float(sigma_clipped_stats(cut)[1])


def host_flux(cut: np.ndarray, bkg: float, pix_arcsec: float) -> float:
    c = (cut.shape[0] - 1) / 2
    yy, xx = np.mgrid[: cut.shape[0], : cut.shape[1]]
    near = np.hypot(yy - c, xx - c) <= HOST_R_ARCSEC / pix_arcsec
    return float(np.clip(cut[near] - bkg, 0, None).sum())


def preprocess(cut: np.ndarray, bkg: float, scale: float) -> np.ndarray:
    return np.arcsinh((cut - bkg) / scale / ASINH_SOFT).astype(np.float32)


def scale_of(cut: np.ndarray, bkg: float) -> float:
    s = float(np.percentile(np.abs(cut - bkg), NORM_PCTL))
    return s if s > 0 else 1.0


def select_hosts(cat: Any, shape: tuple[int, int]) -> np.ndarray:
    ny, nx = shape
    x = np.asarray(cat["xcentroid"], float)
    y = np.asarray(cat["ycentroid"], float)
    mag = np.asarray(np.ma.filled(cat["isophotal_abmag"], np.nan), float)
    ok = (
        np.isfinite(mag)
        & (mag >= MAG_RANGE[0])
        & (mag <= MAG_RANGE[1])
        & (x >= EDGE)
        & (y >= EDGE)
        & (x <= nx - 1 - EDGE)
        & (y <= ny - 1 - EDGE)
    )
    return np.flatnonzero(ok)


def cut_at(sci: np.ndarray, x: float, y: float) -> np.ndarray:
    ix, iy = int(round(x)), int(round(y))
    h = SIZE // 2
    return np.asarray(sci[iy - h : iy + h, ix - h : ix + h], dtype=np.float64)


def embed(name: str, imgs: np.ndarray, batch: int = 32) -> np.ndarray:
    """Pooled pretrained features for preprocessed greyscale cutouts (N, 64, 64)."""
    import timm
    import torch
    import torch.nn.functional as F

    model = timm.create_model(ENCODERS[name], pretrained=True, num_classes=0).eval()
    out = []
    x01 = (imgs - imgs.min(axis=(1, 2), keepdims=True)) / (
        np.ptp(imgs, axis=(1, 2), keepdims=True) + 1e-12
    )
    mean = torch.tensor([0.485, 0.456, 0.406]).view(1, 3, 1, 1)
    std = torch.tensor([0.229, 0.224, 0.225]).view(1, 3, 1, 1)
    with torch.inference_mode():
        for i in range(0, len(x01), batch):
            t = torch.from_numpy(x01[i : i + batch]).float().unsqueeze(1)
            t = F.interpolate(t, size=(224, 224), mode="bilinear", align_corners=False)
            t = t.repeat(1, 3, 1, 1)
            if name == "dinov2":
                t = (t - mean) / std
            out.append(model(t).numpy())
    return np.concatenate(out)


def rank_recall(feats: np.ndarray, is_inj: np.ndarray, seed: int) -> float:
    from sklearn.decomposition import PCA
    from sklearn.ensemble import IsolationForest
    from sklearn.preprocessing import StandardScaler

    z = StandardScaler().fit_transform(feats)
    z = PCA(n_components=min(N_PCA, z.shape[1]), random_state=seed).fit_transform(z)
    iso = IsolationForest(n_estimators=N_TREES, random_state=seed).fit(z)
    score = -iso.score_samples(z)
    top = np.argsort(-score, kind="stable")[:TOP_K]
    return float(is_inj[top].sum() / is_inj.sum())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    import fsspec
    import torch
    from astropy.io import fits
    from astropy.table import Table

    torch.set_num_threads(4)
    t0 = time.time()
    with fsspec.open(S3.format(kind="cat.ecsv"), "rt", anon=True) as f:
        cat = Table.read(f.read(), format="ascii.ecsv")
    with fsspec.open(S3.format(kind="i2d.fits"), "rb", anon=True) as f:
        with fits.open(f, lazy_load_hdus=True) as hdul:
            hdr = hdul["SCI"].header
            sci = np.asarray(hdul["SCI"].data, dtype=np.float32)
    pix = float(abs(hdr.get("CDELT1", np.nan)) * 3600) if "CDELT1" in hdr else np.nan
    if not np.isfinite(pix) and "CD1_1" in hdr:
        pix = float(np.hypot(hdr["CD1_1"], hdr["CD2_1"]) * 3600)
    if not np.isfinite(pix) and "PIXAR_A2" in hdr:
        pix = float(np.sqrt(hdr["PIXAR_A2"]))
    t_io = time.time() - t0
    idx = select_hosts(cat, sci.shape)
    rng = np.random.default_rng(HOST_SEED)
    cand = rng.permutation(idx)
    hosts, cuts = [], []
    for i in cand:
        c = cut_at(sci, float(cat["xcentroid"][i]), float(cat["ycentroid"][i]))
        if c.shape == (SIZE, SIZE) and np.isfinite(c).all():
            hosts.append(int(i))
            cuts.append(c)
        if len(hosts) == N_CLEAN + N_INJ:
            break
    del sci
    cuts_a = np.stack(cuts)
    n_avail = len(hosts)
    clean, inj_hosts = cuts_a[:N_CLEAN], cuts_a[N_CLEAN:]
    bk = np.array([background(c) for c in cuts_a])
    sc = np.array([scale_of(c, b) for c, b in zip(cuts_a, bk, strict=True)])
    hf = np.array([host_flux(c, b, pix) for c, b in zip(cuts_a, bk, strict=True)])
    pre_clean = np.stack([preprocess(c, bk[i], sc[i]) for i, c in enumerate(clean)])
    is_inj = np.r_[np.zeros(N_CLEAN, bool), np.ones(len(inj_hosts), bool)]

    sets: dict[tuple[float, int], np.ndarray] = {}
    arcs: list[dict[str, Any]] = []
    for f_ in STRENGTHS:
        for seed in INJ_SEEDS:
            r = np.random.default_rng(seed)
            pre = []
            for j, c in enumerate(inj_hosts):
                k = N_CLEAN + j
                m, par = arc_model(r, SIZE, pix, f_ * hf[k])
                pre.append(preprocess(c + m, bk[k], sc[k]))
                arcs.append(
                    {"strength": f_, "seed": seed, "host_label": int(cat["label"][hosts[k]]), **par}
                )
            sets[(f_, seed)] = np.stack(pre)

    recall: dict[str, dict[str, list[float]]] = {}
    timing: dict[str, float] = {}
    # Pixel baseline
    recall["pixel"] = {}
    for f_ in STRENGTHS:
        recall["pixel"][str(f_)] = [
            rank_recall(np.r_[pre_clean, sets[(f_, s)]].reshape(-1, SIZE * SIZE), is_inj, s)
            for s in INJ_SEEDS
        ]
    for enc in ENCODERS:
        t1 = time.time()
        e_clean = embed(enc, pre_clean)
        recall[enc] = {}
        for f_ in STRENGTHS:
            recall[enc][str(f_)] = [
                rank_recall(np.r_[e_clean, embed(enc, sets[(f_, s)])], is_inj, s) for s in INJ_SEEDS
            ]
        timing[enc] = time.time() - t1
    outcome = {
        enc: {str(f_): decide(recall[enc][str(f_)], recall["pixel"][str(f_)]) for f_ in STRENGTHS}
        for enc in ENCODERS
    }
    h_emb = {enc: all(v == "beats" for v in outcome[enc].values()) for enc in ENCODERS}
    summary = {
        "provenance": "simulated (injections into observed cutouts); recall is derived",
        "plan": "results/m2_embeddings/PLAN.md",
        "image": S3.format(kind="i2d.fits"),
        "pix_arcsec": pix,
        "n_hosts_eligible": int(len(idx)),
        "n_clean": int(len(clean)),
        "n_injected": int(len(inj_hosts)),
        "n_hosts_used": n_avail,
        "random_recall_at_k": TOP_K / (len(clean) + len(inj_hosts)),
        "params": {
            "size_px": SIZE,
            "mag_range": MAG_RANGE,
            "strengths": STRENGTHS,
            "inj_seeds": INJ_SEEDS,
            "r_range_arcsec": R_RANGE_ARCSEC,
            "extent_deg": EXTENT_DEG,
            "sigma_px": SIGMA_PX,
            "n_pca": N_PCA,
            "n_trees": N_TREES,
            "top_k": TOP_K,
            "margin": MARGIN,
            "encoders": ENCODERS,
        },
        "recall_at_k": recall,
        "mean_recall_at_k": {
            m: {f_: float(np.mean(v)) for f_, v in d.items()} for m, d in recall.items()
        },
        "outcome": outcome,
        "h_emb_supported": h_emb,
        "timing_s": {"io": t_io, **timing},
        "versions": {
            "torch": torch.__version__,
            "timm": __import__("timm").__version__,
            "sklearn": __import__("sklearn").__version__,
        },
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    Table(rows=[[a[k] for k in arcs[0]] for a in arcs], names=list(arcs[0])).write(
        args.out / "injections.ecsv", overwrite=True
    )
    # Contact sheet: first 8 injected hosts at each strength (seed 1) plus 8 clean, preprocessed.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [pre_clean[:8], sets[(STRENGTHS[0], 1)][:8], sets[(STRENGTHS[1], 1)][:8]]
    fig, axs = plt.subplots(3, 8, figsize=(12, 4.8))
    for r_, row in enumerate(rows):
        for c_, im in enumerate(row):
            axs[r_, c_].imshow(im, origin="lower", cmap="gray")
            axs[r_, c_].axis("off")
    for r_, lab in enumerate(["clean", f"f={STRENGTHS[0]}", f"f={STRENGTHS[1]}"]):
        axs[r_, 0].set_title(lab, fontsize=8, loc="left")
    fig.tight_layout()
    fig.savefig(args.out / "contact_sheet.png", dpi=70)
    print(json.dumps({"mean": summary["mean_recall_at_k"], "outcome": outcome}, indent=1))


if __name__ == "__main__":
    main()
