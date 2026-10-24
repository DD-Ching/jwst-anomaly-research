# MARS MACS1149 map vs the Schuldt et al. 2024 new multiple images (held-out stage; plan, frozen)

Frozen 2026-10-24 (commit time of this file), before any residual is computed. Extends D-034 addendum ap
(`../mars_hff_heldout/`), which could not test MACS1149 because no post-MARS image list was stored. Script
`scripts/mars_m1149_heldout.py`; tests `tests/test_mars_m1149_heldout.py` (written before the run).

## Inputs
- model_prediction: MARS free-form deflection maps of MACS J1149.5+2223 (Cha & Jee 2023; Zenodo
  10.5281/zenodo.7804575, md5 checked; arcsec reading), loaded exactly as in `mars_hff_heldout.py`.
- observed: MARS catalogue rows of the field, used only to decide which images MARS saw.
- observed: `schuldt24_new_images.ecsv`, the 18 rows of Schuldt et al. 2024 (A&A 689, A42; arXiv:2405.12287v2)
  Table 2, "new multiple images, in addition to those of Grillo et al. 2016": SXa and families 201, 202, 203, 204,
  301, 401 (17 images, z 2.95–5.98, MUSE spectroscopic, QF 3 or 9). Published after the MARS catalogue (2023).
  ASSUMPTION: listed z are fixed; positions are in the HST (CLASH / HFF) frame MARS uses.

## Question, ordinary explanations, what changes the plan
- Question: does the MARS map, fitted before 2024, reproduce spectroscopic multiple images it never saw on a third
  field (MACS0416 passed, Abell 370 failed, addendum ap)?
- Pass: MARS is a validated second model on MACS1149, and the next stage is a two-model extra-image comparison
  against the CANUCS v1 maps there (gate K 0.990″, `../macs1149_canucsref/`).
- Fail or not testable: no MARS two-model stage on MACS1149. A large residual is **not** a flag here. Ordinary
  causes first: model incompleteness away from the HST-era constraints (these families are at 4–6 z and partly
  outside the core), coarse map pixels near critical curves, a frame offset, the listed z, MUSE position errors
  (MUSE-only images 203, 204, 401), the mean-source estimator.

## Frozen rules
Rules 1–5 of `../mars_hff_heldout/PLAN.md`, unchanged and with the same functions (`frame_offset`, `split_seen`,
`anchors`, `heldout_residuals`, `gate_h`, `outcome` in `scripts/mars_hff_heldout.py`) and thresholds (pair 1.0″,
≥ 5 pairs for an offset, seen 0.5″, gate H rms ≤ 1.0″ and ≥ 80 % within 3.0″, ≥ 5 tested images), with the
Schuldt list in place of the CANUCS list. Training gate K rms for comparison: 0.742″ (addendum ao). SXa is a single
image without a family partner in the list, so rule 3 gives it no anchor (expected untested).

## Stop rules
- No rule change after the first residual is computed; any departure is labelled deviation or post hoc.
- No flag, flux or extra-image reading in this stage. No image pixels read.
- Data: the 2.7 MB zip lives in the session scratchpad and is deleted after the run; only `summary.json` and the
  residual table are stored.
