# Plan: host redshift of AT 2024jfz (frozen 2026-10-19 before any catalogue row near the position was read)

**Why.** Handoff of #360: the counter-image line is closed except z_s 1.55–1.70 (undetermined) and z_s outside
1.35–2.5 (not tested). A host redshift decides which applies.

**Hypotheses.**
- H_bg (hypothesis): the host is a lensed background galaxy (z > 0.39).
- Ordinary: the host is a cluster member or foreground galaxy (z ≤ 0.40: one image, line closed); or the
  transient has no catalogued host (hostless; then the redshift stays unknown).

**Sources (pinned in SOURCES.md).** Caminha+22 MUSE (J/A+A/666/L9, `offplane_potentials.MUSE["smacs0723-iclv2"]`,
moved to the Mahler+22 frame); DJA v7.4 SMACS 0723 phot + eazy zout (sha256 pinned).

**Frozen rule (`host_z.classify`, unit-tested before the run).** Candidates: catalogue sources within
R_MAX = 1.0″ of the TNS position (ASSUMPTION). Host = the DJA source with the smallest normalised distance
d / r_half (DJA `flux_radius`, pixels × 0.04″, floor 0.1″); MUSE z_spec with a trusted flag within 0.5″ of the
host's DJA position overrides the photo-z. Outcome on the host's z interval (z_spec ± 0.01, else eazy z160–z840):
- `member_or_fg` if the upper edge ≤ 0.42;
- `tested_window` if the interval lies inside 1.35–1.50 or 1.75–2.50 (counter-images absent there);
- `undetermined_window` if it overlaps 1.55–1.70;
- `untested` otherwise (any other z, or no host within R_MAX).
**What changes the plan:** `undetermined_window` → a PSF-matched difference or same-instrument epoch is the next
step; `member_or_fg` or `tested_window` → close the line; `untested` → list what would decide it.
