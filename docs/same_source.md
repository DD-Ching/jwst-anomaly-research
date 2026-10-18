# Recognising one source at several positions ("same-source" pairing)

Owner instruction 2026-10-18. Inventory, tool evaluation and the limits on what a later search may allow. Result of
the first feasibility test: `results/same_source/a2744_feasibility/RESULT.md`.

## What the repository can and cannot do (2026-10-18)
| Capability | Status | Where |
|---|---|---|
| Align two epochs of one sky region | working, calibrated per box | `arc_diff.py`, `batch_*.py` |
| Find brightening / fading point sources between epochs | working, injection-calibrated | `arc_diff.py`, `fade_calib.py` |
| Predict counter-images of a given position with a published lens model | working (dPIE port, map models) | `lensmodel.find_images`, `lens_consistency.py images`, `arcs_counter.py` |
| Test catalogued multiple images for flux / colour / position anomalies | working | `lens_consistency.py`, D-029–D-045 |
| **Task 1: propose which sources are the same source, without a model** | **not available**: SED pairing at 0.3–3″ (`orphan_pairs.py`, dark-lens screen) is null-dominated; at 2–60″ it fails the 2026-10-18 test (26,145 SED-matched pairs vs 11 true ones in Abell 2744) | `same_source_pairs.py` |
| **Task 2: test a given pair against ordinary lensing** | partly: positions and fluxes with published models; the morphological (flip / stretch) test is paused, failing its acceptance test (D-032) and not discriminating true pairs from decoys (AUC 0.53, 6 + 6 pairs) | `parity_test.py` |

Image registration and lens-model root finding are not same-source searches: they assume which pixels or which
system belong together.

## External tools evaluated (2026-10-18, isolated venv in the session scratch; not added to `pyproject.toml`)
| Tool | Licence, version | Installs / runs | What it solves |
|---|---|---|---|
| lenstronomy | BSD-3, 1.14.2 | yes (Python 3.12, ~35 s with PyAutoLens) | task 2 for galaxy-scale systems: forward-models given images with a lens and source; no catalogue pair proposal |
| PyAutoLens | MIT, 2026.10.17.1 | yes | task 2, galaxy (and group) scale; automated pipelines; no pair proposal |
| Lenstruction | MIT, last commit 2022-11, pinned to lenstronomy 1.10.3 | imports with 1.14.2 except `data_process` (needs photutils) | task 2 for cluster arcs: reconstructs one source from images already identified as one system, with local lens corrections; needs the image list and a cluster model as input |
- None proposes same-source pairs from a catalogue (task 1). They would refine task 2 for a pair already in hand; the
  repository's published-model checks cover the ordinary-lensing test at cluster scale, so no tool was adopted this
  round. Revisit with a specific pair whose ordinary explanation needs a source-plane reconstruction.

## Limits for any later search (hypothesis constraints)
- **Identity of data:** sky coordinates, observation date, band, programme and product; never pixel IDs.
- **Differences a model may allow, each declared before the run with its range:** rotation, mirror (parity), stretch
  (affine singular values), magnification, dust and band-dependent contamination, redshift (photo-z interval), time
  delay, source evolution, position shift, blending, non-detection, several foreground deflectors.
- **What is not allowed:** free per-pixel warps, per-line or per-epoch shifts until two data sets agree, or adding
  freedoms after seeing a score. Each extra freedom needs a decoy test showing it does not make unrelated sources
  match (the 2026-10-18 test: redshift freedom alone gives ~2,400 unrelated SED pairs per true pair).
- **Competing explanations kept for every pair:** two different sources each lensed; group members at one redshift;
  knots of one galaxy; instrument artefacts (persistence, spikes, edges).
- **Geometry:** no assumption that an exotic lens is a flat ring seen in projection; no exotic shape template is used.
- **Ordinary lensing:** a failed ordinary fit is not evidence of a wormhole, and an ordinary model is not required for
  keeping a candidate; a model-free test (e.g. time delays, parity, flux ratios) must be named instead.
- **Missing data:** without a spectrum or enough epochs the corresponding test is "not testable", never "passed".

## Next step if reopened
- A ranking-based proposal (relative SED similarity, which put the true counterpart first in 25 of 39 queries) needs a
  calibrated threshold and a new held-out field (e.g. MACS0416 CANUCS); tuning on Abell 2744's holdout is not allowed.
