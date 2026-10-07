# Roadmap

| Milestone | Goal | Exit criterion |
|---|---|---|
| **M0 Bootstrap** (done 2026-10-07, pending merge) | Durable repo structure, interfaces, CI, agent harness, reuse decisions | Skeleton + 9 bootstrap PRs merged |
| **M1 Catalog-level slice** (first run 2026-10-07) | Rank sources in program 2736 pipeline catalogs; cutouts and cross-check for the top k; candidate report; then quality gating and star/galaxy separation | One reproducible run on real data with honest results and limitations recorded; top-k artifact fraction tracked |
| **M2 Better representations** | Consistent matched-aperture photometry (DJA v7.4 catalogs, checked against CEERS DR1.0), then image embeddings (Zoobot via timm, DINOv2) | Documented gain over the M1 baseline on injection-recovery and top-k artifact fraction |
| **M3 Lensing-focused search** | Lens/arc-specific features (e.g. tangential alignment in cluster fields), AnomalyMatch-style lens finding, consistency checks against published SMACS 0723 models (Mahler+2022, RELICS, Caminha+2022) | Vetted lensing candidate list with tests recorded |
| **M4 Scale-out** | Many programs and fields via cloud access, batch orchestration, candidate DB at scale | Survey-scale run with provenance |

Direction changes when evidence says so; record why in DECISIONS.md. Tool choices for M2–M4 are in
DECISIONS.md D-010 and docs/landscape.md.
