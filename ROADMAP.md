# Roadmap

| Milestone | Goal | Exit criterion |
|---|---|---|
| **M0 Bootstrap** (done 2026-10-07, pending merge) | Durable repo structure, interfaces, CI, agent harness, reuse decisions | Skeleton + 9 bootstrap PRs merged |
| **M1 Catalog-level slice** (done 2026-10-08: clean SMACS galaxy top 20, D-011 to D-021) | Rank sources in program 2736 pipeline catalogs; cutouts and cross-check for the top k; candidate report; then quality gating and star/galaxy separation | One reproducible run on real data with honest results and limitations recorded; top-k artifact fraction tracked |
| **M2 Better representations** | Consistent matched-aperture photometry (DJA v7.4 catalogs, checked against CEERS DR1.0), then image embeddings (Zoobot via timm, DINOv2) | Documented gain over the M1 baseline on injection-recovery and top-k artifact fraction |
| **M3 Lensing-violation search** (current priority, D-023) | Test observations against published lens models (SMACS Mahler+2022 first): predicted vs observed image positions, parities and flux ratios; arcs whose orientation or curvature disagrees with the model shear; images where no multiple images are predicted. Exotic-lens signatures (demagnification, configurations impossible for positive mass) come from verified literature. Time-domain: two epochs (VENUS) for caustic-crossing transients | Lensing-consistency tests run on every cluster field; each inconsistency vetted against ordinary explanations (substructure, microlensing, model error, artifacts) |
| **M4 Scale-out** (started 2026-10-08: Abell 2744, El Gordo, Sunrise in parallel) | Many programs and fields via cloud access, batch orchestration, candidate DB at scale | Survey-scale run with provenance |

Direction changes when evidence says so; record why in DECISIONS.md. Tool choices for M2–M4 are in
DECISIONS.md D-010 and docs/landscape.md.
