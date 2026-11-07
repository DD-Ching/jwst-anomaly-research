# M2 image-embedding pilot on SMACS 0723 F200W (plan, frozen before any embedding is computed)

TASKS Next M2 "Image embeddings"; D-022 representations bullet ("any embedding must beat the sklearn baseline on
injection-recovery"). Script: `scripts/m2_embeddings.py`; outcome function `decide` unit-tested before the run.

## Question and hypothesis
- **Question:** do pretrained encoders (Zoobot `convnext_nano` greyscale-trained encoder, DINOv2 ViT-S/14), used
  through `timm` without fine-tuning, rank image-level lens-like injections higher than a simple pixel baseline?
- **Hypothesis (H-emb):** Zoobot recall@50 exceeds the pixel baseline's at both injection strengths.
- **Ordinary explanations for any gain or loss:** the normalisation and resize (natural-image encoders see a
  64 px cutout upsampled to 224 px), noiseless injections (optimistic for every method alike), host selection.
- **What changes the plan:** if no encoder beats the baseline, embeddings are not adopted for M2 ranking and the
  next M2 step is a fine-tuned or JWST-trained encoder only on a new input (D-022 Revisit); if one beats it, the
  next step is a held-out field (CEERS t021) and a real-source ranking with a contact sheet.

## Data (observed)
- `jw02736-o001_t001_nircam_clear-f200w_i2d.fits` SCI (10277 × 4762 float32, 196 MB) streamed into memory from
  `s3://stpubdata` (anonymous); not stored on disk (cloud disk rule). Hosts from the matching pipeline `_cat.ecsv`.
- **Hosts (ASSUMPTION):** `is_extended` sources with 21 ≤ `isophotal_abmag` ≤ 25, 64 × 64 px cutout fully finite,
  centre ≥ 40 px from the image edge; 1050 drawn with `numpy` seed 0: 1000 clean + 50 injection hosts (disjoint).

## Injection (simulated)
- An arc: circular-ring segment of radius r ~ U(0.35″, 0.70″) around the cutout centre, angular extent
  ~ U(90°, 180°), random orientation, Gaussian radial profile σ = 1.5 px, uniform along the arc; total flux
  = f × host flux (sum of background-subtracted pixels within 0.3″). Strengths f ∈ {0.3, 1.0}. Noiseless model added to
  the real cutout (ASSUMPTION). Three injection seeds (1, 2, 3); each draws new arc parameters on the same 50 hosts.

## Representation and ranking
- Cutout preprocessing (all methods): subtract the sigma-clipped median, divide by the 99.5th percentile of the
  absolute values of the clean cutout of that host, arcsinh(x / 0.1) (ASSUMPTION).
- **Pixel baseline:** flattened 64 × 64 preprocessed pixels.
- **Zoobot:** `hf_hub:mwalmsley/zoobot-encoder-convnext_nano`, pooled features; input resized to 224 px, greyscale
  repeated to 3 channels, scaled to [0, 1].
- **DINOv2:** `vit_small_patch14_dinov2.lvd142m` (timm), CLS-pooled features, 224 px, ImageNet normalisation.
- Each representation → standardise → PCA 50 → IsolationForest (300 trees, `random_state` = injection seed) on
  the 1000 clean + 50 injected rows; score = −`score_samples`.
- **Metric:** recall@50 = injected rows in the top 50 of 1050, per method, strength and seed (random 0.048).

## Frozen rule (`decide`, per encoder and per strength)
- `beats`: mean recall@50 over the 3 seeds ≥ baseline mean + 0.10 **and** recall ≥ baseline in each seed.
- `worse`: mean ≤ baseline mean − 0.10.  Otherwise `tie`.
- **H-emb supported** for an encoder only if it `beats` at both strengths.
