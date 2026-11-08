# `arcs`: does the RELICS GLAFIC v3 model treat RXC J0600 581 as a cluster member? (exploratory, 2026-11-08)

Plan: `PLAN.md` (frozen and committed before any κ value was read). Script `scripts/arcs_kappa_members_581.py`;
peak, gate, class and label functions unit-tested before the run (`tests/test_arcs_kappa_members_581.py`, committed
before the run). Every number below is in `summary.json`. κ maps (RELICS HLSP, DOI 10.17909/T9SP45; sha256 in
`summary.json`) and the NED `objdir` table were read in the scratchpad and deleted. Redshifts `observed` (NED);
κ `model_prediction`; peak flags and labels `derived`; thresholds ASSUMPTIONs (`summary.json` `params`).

## Frozen outcome: **member (model input)** for the GLAFIC v3 model
| Map | Pixel | Members (n) / peak recall | Background (n) / peak rate | Gate G | 581: κ max, distance, max / median | Label |
|---|---|---|---|---|---|---|
| GLAFIC v3 (deciding) | 0.1″ | 41 / 0.854 | 60 / 0.0 | passes | 5.24, 0.223″, 6.43 → `peak` | **member (model input)** |
| Lenstool v1 (secondary) | 0.2″ | 41 / 0.805 | 60 / 0.0 | passes | 2.31, 1.645″, 1.12 → `no_peak` | background (model input) |

116 NED spectroscopic galaxies (G / G_Lens, `zflag` SLS or S1L) fall inside both maps; 15 foreground ones are not used.

## Reading (narrowest statement)
- In both models a galaxy-scale κ peak marks a spectroscopic cluster member in ≥ 80 % of cases and never marks a
  spectroscopic background galaxy (0 / 60). The peak is a usable proxy for "the model treats this galaxy as a member".
- GLAFIC v3 has such a peak on 581: the GLAFIC modellers' input treats 581 as a cluster member. Its GLAFIC axis
  prediction at 581 is therefore not a test of a background arc, and the GLAFIC run's `anti` class for 581 is not a
  lensing-model residual.
- Lenstool v1 has no peak there: it treats 581 as background (or below its member cut), and under Lenstool 581 stays
  `anti` in 100 / 100 realisations (addendum 13). The two models disagree on the input, not only on the shear.
- This is about model inputs, not about 581's true redshift. Spectroscopy still decides membership (none public,
  addendum 15); 581's photo-z are 0.70 (DJA eazy) and 0.537 (RELICS BPZ, lower bound below the member cut; addendum 10, `../photoz-contam/summary.json`).
- Not a lensing anomaly claim and not a candidate.

## Looked at (post hoc figure, not a rule)
`kappa_581.png` (log κ, 8″ boxes centred on 581, made after the run): GLAFIC shows an elongated galaxy halo on 581
along its disc and a second halo about 3″ north; Lenstool shows only the northern halo. Consistent with the table.

## Deviations (labelled)
- PLAN.md says the pre-plan NED count was 78 spectroscopic G / G_Lens rows; the correct count of the same table is 116
  (SLS 69 + 39, S1L 8; the split by flag is not stored). A counting slip in the plan text; the rules did not depend on it.
- After the first run the label text was changed from "GLAFIC member / background" to "member / background (model
  input)", because the same function labels the Lenstool map; rules and numbers unchanged (re-run gave the same
  gates and peaks).
- `kappa_581.png` is post hoc.
