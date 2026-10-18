# W1 shear screen at lower masses: filter pilot (plan frozen 2026-10-18 11:20Z, before any injection or screen)

**Question.** The forty-field W1 shear limit (`results/w1_shear_primer/RESULT.md`) reaches 0.638× Takahashi &
Asada (2013) at 2 × 10¹³ M☉ but only 9.32× at 2 × 10¹² M☉, and the screen is blind at ≤ 2 × 10¹¹ M☉. The PRIMER
plan names lower masses as the next step. Before spending a forty-field grid at new masses, this pilot asks whether
the frozen filter (Schirmer, R = 10″, r_min = 1″; D-053) is the right one for small lenses.

**Reasoning (derived, order of magnitude; ASSUMPTION: weak-shear point mass).** A negative point mass leaves no
images inside θ_E, and its tangential shear falls as (θ_E/θ)². At 2 × 10¹² M☉ and Z = 0.3, θ_E(z_s = 2) = 3.63″
(`theta_e_zs2_arcsec` in the CEERS summaries), so most of the signal sits at θ_E–2θ_E ≈ 1.8–7″ for 5 × 10¹¹–2 × 10¹²
M☉. The Schirmer weight of a 10″ aperture peaks near 4″ and falls inward, so a filter that weights small radii more
(the point-mass filter, or a smaller aperture) may recover more small lenses.

**Hypothesis being limited.** As `results/w1_volume/PLAN.md` (W1 lenses uniform in comoving volume). This is a
method pilot: recoveries of simulated lenses only.

## Inputs (frozen)
- **Fields:** `ceers-o001` and `goodsn-o004`, the two blank pointings with the most 2 × 10¹² M☉ recoveries in the
  committed Z = 0.3 summaries (14 / 100 and 18 / 100). Chosen by baseline efficiency (stated selection).
- Inputs staged with `w1_stage_inputs.py` (pinned catalogues and zouts). **No streak masks** in the pilot (the
  SCI-mosaic streaming is skipped to save time; deviation from the headline runs, applied to every variant alike).
  V0 is therefore re-run here rather than read from the committed (masked) summaries.
- `inject_shear.py --lens-z 0.3 --mass 2e12 1e12 5e11 --n-inject 100 --seed 50`, other options at their defaults
  except the filter settings below.

## Variants (frozen)
| Variant | `--filter` | `--aperture-arcsec` | `--r-min-arcsec` |
|---|---|---|---|
| V0 (frozen screen) | schirmer | 10 | 1 |
| V1 | pointmass | 10 | 1 |
| V2 | schirmer | 5 | 1 |

## Decision rule (frozen; thresholds are ASSUMPTIONs)
- Per variant, N = recoveries at 1 × 10¹² plus 2 × 10¹² M☉, per field.
- **Adopt** the variant with the largest summed N for the lower-mass grid if it exceeds V0's N by ≥ 50 % in **both**
  fields. Then the next step is that variant's real-data screen on all forty fields, with its E/B flags vetted
  under the frozen rule (`results/w1_shear_more` addendum A, V1 ∪ V2) before any limit.
- **Otherwise** V0 stays, and the lower-mass step is a V0 grid at 1 × 10¹² M☉ over the forty fields (cost about
  one PRIMER run).
- 5 × 10¹¹ M☉ is reported but not used in the rule (expected near zero).
- **Real-data screens at V1 / V2 are a new look at the two fields.** Their S_max and E/B outcome are recorded
  (`screen` block of each summary) but not vetted in this pilot; they enter no limit.

## What would change the plan
- If no variant recovers ≥ 5 / 100 at 1 × 10¹² M☉ in either field, the shear screen cannot open lower masses on
  JWST depth, and W1 below 2 × 10¹² M☉ belongs to the radial screen or deeper imaging.
