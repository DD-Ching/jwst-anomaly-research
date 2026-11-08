# Match / search radii scaled to the image-plane rms (plan, frozen before the run)

TASKS Supporting 1, `forced_check` item: "match and search radii scaled to the image-plane rms when it exceeds 1″".

## Question
In the stored `validate` image-plane tables (`results/lens_consistency/*/imageplane*.ecsv`, column `dtheta_arcsec`,
offline), is an rms above 1″ a broad model scatter (then scaling the 1.5″ match radius by the rms is justified) or
driven by one image (then the rms is the wrong scale and a robust one is needed)?

## Ordinary explanations to keep apart
- One catalogued radial / central image far from its prediction (MACS0417 1b.4, MACS1423 system 1; D-034 readings).
- A model that fails gate K outright (rms of several arcsec): scaling cannot rescue it and must not try.

## Frozen rule (functions decide; thresholds are ASSUMPTIONs)
- `rms_scaled_radii(rms, match=1.5, search=1.0, threshold=1.0, cap=3.0)`: rms ≤ threshold → (match, search)
  unchanged; else both × rms / threshold, each capped at `cap` (= `MAX_SEARCH_ARCSEC`).
- `rms_outlier_driven(dtheta, threshold=1.0)`: True when rms > threshold and the rms with the single largest offset
  removed is ≤ threshold.
- Per table with finite offsets (n ≥ 3): rms, leave-max-out rms, robust scale 1.4826 × median |dθ|, unmatched counts
  at 1.5″ and at the scaled match radius.
- Reading, over tables with 1″ < rms ≤ 2″ (scaled match ≤ the 3″ cap): **"rms scale outlier-driven"** if
  `rms_outlier_driven` holds in at least half of them; else **"rms scale broad"**. Tables with rms > 2″ are
  "gate-K failure, not scalable" and are only counted.
- Changes the plan: "broad" → wire the scaling into `images` (default off) next; "outlier-driven" → do not wire it;
  the next step is a robust scale with its own gate.
