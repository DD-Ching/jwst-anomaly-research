# Bullet Cluster: gate K with the WLS family source (deviation, frozen after gate K failed)

Frozen 2026-10-24 (commit time of this file), after `PLAN.md`'s gate K failed (`RESULT.md`: rms 1.196″; K1b and
K1e get one predicted image under the mean source). This is the Abell 2744 precedent (D-034 addendum k,
`../abell2744_uncoverref/PLAN_1b.md`) applied unchanged: **post-failure deviation**, labelled so wherever quoted.

## Rule (only change from `PLAN.md`)
- Each family's source is D-041's magnification-weighted least-squares source (`a2744_uncover.wls_source` via
  `a2744_uncover.imageplane`) instead of the mean back-traced source. Same 135 gold images, same z, same maps, same
  0.25″ grid, same gate K function and thresholds (rms ≤ 1.0″ **and** ≥ 80 % within 3″).
- Pass → the maps go on to gate C2 / extra-image stages **with the WLS source**, as on Abell 2744 (own-constraint
  validation only). Fail → the route stops on the Bullet Cluster.
- Reported, not gates: families with fewer or more predicted than catalogued images.
