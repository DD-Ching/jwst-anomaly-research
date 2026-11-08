# `arcs` MACS0308 2157 MUSE redshift, second plan (written before any data read)

Follows `../muse-2157/RESULT.md` (gate M failed as frozen). Same cube (ESO `ADP.2020-01-04T23:56:39.902`), streamed
once, nothing stored. The earlier readings are known (control peak contrast 1.79 off-peak, target 1.07), so the
thresholds below are not blind; they are frozen here before this run.

## Hypothesis and ordinary explanations
- H: 2157 (anti row, `robust` in 100 / 100 Zitrin realisations) is a cluster member at z ≈ 0.356 (then the anti
  orientation is no lensing residual) or a background source (then it stays an open anti row).
- Ordinary readings first: member galaxy; foreground galaxy; background source with no line in range.

## Changes from the first plan
1. Apertures (0.6″) recentred on the flux-weighted white-light centroid within 1.0″ of the catalogued position
   (ASSUMPTION: 1.0″ search radius), iterated 3 times.
2. Noise: the cube's STAT extension (variance) is streamed too; per pixel the noise is the larger of
   sqrt(sum STAT) and the random-aperture MAD (conservative; ASSUMPTION). Both are reported.
3. Gate M2 (validation, central galaxy): |z_best − 0.356| < 0.01 **and** peak contrast ≥ 1.5 (best score over the
   highest score more than 0.01 away; ASSUMPTION 1.5).
4. Target reading, only if M2 passes: contrast ≥ 1.5 → member (|z − 0.356| < 0.02), background (z > 0.40), else
   foreground; contrast < 1.5 → `undetermined`.
5. Line-flux limit: [O II] 3726.0 / 3728.8 doublet (1:1, Gaussian σ 1.1 Å observed, ASSUMPTION) injected into the
   target aperture spectrum at 40 random z in 0.40–1.45 per flux level. Recovered: |z_best − z_inj| < 0.002 and
   [O II] window S/N ≥ 5. F90 = lowest flux level with ≥ 90 % recovered (aperture flux, no aperture correction).

## What changes the plan
- M2 fails → the absorption route is closed for this cube; only the [O II] limit is reported.
- Target `member` → 2157 leaves the anti list; `background` → stays open, z recorded; `undetermined` → limit only.

## Amendment 1 (2026-11-08, before any data read)
The unit test showed that a lone injected [O II] doublet is z-degenerate under the scan (other lines at other z
score the same: z_best 0.2656, contrast 1.01 for a z 0.7 injection, not stored). Item 5's recovery is therefore
replaced by: **detected = [O II] window S/N ≥ 5 at z_inj**; F90 uses detection. The z-identified fraction
(|z_best − z_inj| < 0.002) is reported as a diagnostic only.
