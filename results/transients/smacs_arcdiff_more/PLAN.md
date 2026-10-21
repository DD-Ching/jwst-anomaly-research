# SMACS 0723 arc difference imaging, remaining arc boxes (ERO o001 → VENUS o057) — plan (written before the run)

**Why.** TASKS.md (SMACS 0723 arc `arc_diff`): "Next: the rest of the SMACS arcs (more boxes)". Pair 1 of
`results/transients/smacs_arcdiff/` (jw02736-o001 2022-06-07 → jw06882-o057 2026-06-05, F150W + F444W) passed its
known-case control (SN Atalanta A/B) and gave F150W recovery 0.82 / 0.60 at 27.5 / 28.0 AB on the system 5 box.
The other 45 Mahler+22 ICLv2 images lie outside that box.

**Boxes** (24″ × 24″; greedy cover of every ICLv2 image (JWST frame) not ≥ 1″ inside the system 5 box, each image
≥ 1″ inside its box; computed from `arcs.dat`, not from the data):

| box | RA | Dec | images |
|---|---|---|---|
| b1 | 110.812643 | −73.457179 | 4.1 4.2 5.3 8.1 9.1 9.2 10.3 11.1 11.2 12.3 15.2 15.3 17.2 18.3 19.3 21.2 21.3 |
| b2 | 110.830002 | −73.450431 | 1.1 2.1 3.1 6.1 6.2 12.1 13.1 15.1 17.3 (holds Atalanta A/B; overlaps the system 5 box) |
| b3 | 110.832266 | −73.457832 | 1.2 1.3 2.2 2.3 3.3 6.3 17.1 |
| b4 | 110.802563 | −73.449931 | 4.3 7.1 7.2 7.3 25.1 25.2 25.3 |
| b5 | 110.789807 | −73.458128 | 8.2 26.1 26.2 |
| b6 | 110.813362 | −73.464486 | 13.3 |
| b7 | 110.807082 | −73.451355 | 21.1 |

Command per box: as pair 1 (`--kernel-px 4 --kernel-stamps --inject 100 --min-bands 2 --procs 4`), `--known`
only on b2. Then `arc_diff_morph.py --all-hits` and `smacs_arcdiff_images.py`. Order b1, b2, b3, b4, then b5–b7
as time allows; boxes not run are reported as not run.

**Validation first (ASSUMPTIONs, fixed here).**
- b2: Atalanta A and B must be hits passing morphology (the known case, now off-centre in a new box). If not, b2
  is a tool failure, not a null.
- Each box: injection recovery ≥ 0.50 at 27.5 AB in both bands, else the box gives no depth statement (hits are
  still recorded).

**Hypothesis (to test).** No point-like flux change on the boxes' arcs beyond instrument, processing and known
populations. Ordinary explanations first: kernel residuals on bright cores, spikes, persistence, registration
dipoles, pipeline-version differences, ordinary caustic microlensing, supernovae.

**Frozen rule.** A non-control morph `pass` within **2.0″** of an ICLv2 image (pair 1's post hoc radius, adopted
here prospectively because 1.0″ failed the known case) goes to multi-epoch forced photometry
(`hit_lightcurve.py`); a change ≥ 3σ in either band there → `/vet-candidate`. Everything else is recorded only.

**What changes the plan.** A `/vet-candidate` survivor → needs-human PR. A clean set of validated boxes → depth
statement per box and the SMACS arc line closes for this pair.

## Addendum 1 (2026-10-21T16:21Z, after b1, before its control test): b1 h019
- b1 h019 (morph `pass`, 1.51″ from 17.1; arc_diff sign −) has forced-photometry changes 2022 → 2026-06 of
  +73.6 ± 16.1 nJy F150W (4.6σ) and +33.4 ± 5.1 nJy F444W (6.5σ) on a ≈ 1.7 / 1.0 μJy source, so the frozen rule
  sends it to `/vet-candidate`. A uniform +3–4 % with the opposite sign to the kernel-matched difference points first
  to an epoch-to-epoch flux-scale (zero-point) difference, which the random-aperture errors do not include.
- **Test (frozen here):** `hit_lightcurve_controls.py pick` on the 2022 F150W mosaic, flux range 500–5,000 nJy,
  n 60, half-box 16″; `hit_lightcurve.py` at h019 + controls; `compare` for the o001 → o057 pairs in F150W and
  F444W. z_cal = (z_h019 − median z_ctrl) / (1.4826 × MAD z_ctrl) per band. **explained** |z_cal| < 2 in both
  bands; **open** |z_cal| ≥ 3 in either (→ counter-image check of system 17); else **undetermined**. ASSUMPTIONs.

## Addendum 2 (2026-10-21T16:23Z, before b2 h032's photometry)
- b2 passes its known-case gate (Atalanta A = h021, B = h023, both `pass`). b2 h032 is a morph `pass` 0.09″ from
  image 1.1. It gets the forced photometry of the frozen rule and, if that shows ≥ 3σ, the Addendum 1 control test
  with the same thresholds (controls picked around h032). Every later rule hit gets the same two steps.

## Addendum 3 (2026-10-21T16:28Z, after the b2 h032 control test, before these steps)
- b2 h032 control test (Addendum 1 thresholds): F150W z_cal 0.44, **F444W z_cal 5.65** (ratio 1.153, control tail
  fraction 0.029) → **open**. On the stamps the aperture sits on the steep side of the elongated image 1.1, where a
  small WCS disagreement between the two mosaics changes the aperture flux. Other epochs: F444W 2024 ratio 1.045,
  2026-07 1.104; F300M changes 1,111 → 1,368 nJy in one month (2026-06 → 2026-07).
- **Step R (frozen):** `reg_offsets.py --pair 6882:<o057 cat> --ref 2736:<o001 cat> --bands f150w f444w
  --offsets-at <h032>` (its own frozen model choice), then `hit_lightcurve.py --offsets` at h032 + the same
  controls and `compare`. **explained by registration** if the corrected F444W |z_cal| < 2; **survives** if ≥ 3;
  else undetermined.
- **Step C:** forced photometry at the ICLv2 counter-images 1.2 and 1.3 (shifted catalogue positions), same pair and
  aperture, reported only (time delays not computed; a real change in 1.1 need not be simultaneous in them).
