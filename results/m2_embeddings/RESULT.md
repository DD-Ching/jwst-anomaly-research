# M2 image-embedding pilot (SMACS 0723 F200W): no encoder meets the frozen rule at both strengths; H-emb not supported

Plan: `PLAN.md` (frozen before any embedding; one pre-run deviation below). Script: `scripts/m2_embeddings.py`.
Numbers: `summary.json`; arc parameters per host, strength and seed: `injections.ecsv`. Contact sheet:
`contact_sheet.png` (looked at; the same 8 injection hosts without and with the seed-1 arcs: the arcs are visible by
eye at both strengths on 7 / 8; on the bright point source in the last column they are not; cause not measured).

## Data (observed) and sample (derived)
- `jw02736-o001_t001_nircam_clear-f200w_i2d.fits` SCI streamed from `s3://stpubdata` (not stored); pixel scale
  0.0312″. Eligible hosts 2655 (20–27 mag, inside the edge margin); 1000 clean + 50 injection hosts drawn (seed 0).
- **Deviation (pre-run):** the plan's first host rule (`is_extended`, 21–25 mag) gave too few hosts (count in the
  PLAN deviation note, not stored in `summary.json`); amended to any source at 20–27 mag before any embedding.

## Result (simulated injections, derived recall@50 of 1050; random 0.048)
| Method | f = 0.3 (seeds 1, 2, 3) | mean | f = 1.0 (seeds 1, 2, 3) | mean | Outcome (0.3 / 1.0) |
|---|---|---|---|---|---|
| Pixel PCA + IF (baseline) | 0.02, 0.04, 0.04 | 0.033 | 0.08, 0.10, 0.06 | 0.08 | — |
| Zoobot `convnext_nano` | 0.18, 0.18, 0.16 | 0.173 | 0.14, 0.20, 0.18 | 0.173 | beats / tie |
| DINOv2 ViT-S/14 | 0.04, 0.04, 0.06 | 0.047 | 0.08, 0.06, 0.04 | 0.06 | tie / tie |

- **Frozen rule:** H-emb needs `beats` at both strengths; Zoobot misses at f = 1.0 (mean margin 0.093 < 0.10), so
  **H-emb is not supported**; DINOv2 ties the baseline at both strengths.
- Every method stays at recall ≤ 0.20 per seed although the arcs are visible by eye: an unsupervised Isolation
  Forest on whole-cutout features from this mixed population (stars, cluster members, arcs, faint galaxies) does not
  rank added arcs near the top. This describes this configuration only (global pooled features, PCA 50, IF); it is
  not a test of a supervised or fine-tuned lens finder.
- Post hoc (not part of the rule): Zoobot's recall does not rise from f = 0.3 to 1.0 (both means 0.173), so its
  gain over the baseline does not scale with arc brightness in this test; a cause is not isolated. The design is
  paired: each seed draws the same arc geometry at both strengths (only the flux differs), so the two strengths are
  not independent tests.
- Post hoc (`summary.json` `posthoc_host_scaled`): the frozen encoders' input min-max scales each injected cutout by
  its own range, so an arc also changes the host's contrast. Scaling injected cutouts by their host's clean range
  instead gives Zoobot 0.16 / 0.16 and DINOv2 0.047 / 0.08, so Zoobot's gain at f = 0.3 does not come from that
  rescaling.

## Limitations
- Noiseless arcs (optimistic for all methods alike); a single field and band; 50 injections per seed, so one
  injected source changes recall by 0.02; global pooled features; no catalog-feature comparison, because image-level
  injections would need re-extraction (the catalog baseline cannot see them without it).
- Encoder input is min-max scaled per image (frozen); the pixel baseline is not (see the post hoc check above).
- The full SCI array (196 MB) is read into memory once, not stored: 1050 byte-range cutouts of 64 full-width rows
  would move more bytes (estimate in the script docstring, not stored).
- Runtime on 4 CPU cores (stored; both variants): Zoobot 36 s, DINOv2 58 s, each for the 1000 clean cutouts plus 50 × 2 strengths × 3 seeds injected ones, twice (frozen and post hoc).

## Next
- Embeddings are not adopted for M2 ranking. A lens-specific route stays as D-010 lists it (AnomalyMatch, then a
  fine-tuned Zoobot), which needs labelled JWST lens examples; a JWST-trained encoder is the D-010 Revisit trigger.
