# W1 faint-streak veto for the CANUCS flanking fields (plan, frozen 2026-10-18 before any mask, screen or injection)

Cloud routine, claim PR #324. TASKS "JWST focus" 1 (W1), handoff of `results/w1_shear_ncf/RESULT.md`: a row/diagonal
streak veto, frozen and run through the injections, then the two flagged NCF re-evaluated.

**Hypothesis (H).** As `results/w1_volume/PLAN.md`: W1 lenses give a radial aperture-mass peak.
**Ordinary explanation tested.** Faint linear artefacts the D-097 mask misses (findsat_mrt threshold 5): the
diagonal streak through the macs0416-ncf S_max aperture and the horizontal stripes through the macs0417-ncf
aperture (seen on `results/w1_shear_ncf/sheet_smax_flags.png`). D-097's "Revisit if" names exactly this case: a lower
threshold, recalibrated through injections.

**Veto v2 (frozen; every number an ASSUMPTION).** `scripts/w1_streak_mask.py` unchanged except the TrailFinder
threshold, which becomes a `--mrt-threshold` option (default 5, the D-097 value). Variants, tried in this order
on all five NCF F150W mosaics, the first that passes K1n adopted, no other variant tried:
- A: threshold 3;
- B: threshold 2.5.

**Known-case check K1n (computable).** With the 1″ dilation, the variant's mask covers the S_max position
(`screen.s_max_radec` of `results/w1_shear_ncf/summaries/z0.3_<field>.json`) of both macs0416-ncf and macs0417-ncf.
Neither variant passes → the faint-streak route fails, recorded as a failed approach; the headline stays at
fourteen fields.

**Cost check (computable).** A variant that vetoes > 25 % of used shape rows in any NCF is rejected as too costly
(go to the next variant). The cost inside the limit is measured by the injections, not by this check.

**Run (if a variant is adopted).** All five NCF re-run with the v2 masks (`inject_shear.py --streak-masks`, the
frozen `results/w1_shear_ncf` grid: Z ∈ {0.1, 0.3, 0.6, 1.0, 1.5}, masses 2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉, 100
lenses per cell, seed 50, R 0.425). An NCF whose v2 screen passes the E/B rule stays out of the limit and goes to
`w1_peak_vet.py` / `/vet-candidate`. **Combined limit:** `w1_volume.py` over the eleven stored cluster summaries
(D-097 masks, unchanged) plus every NCF whose v2 screen does not pass, reported beside the fourteen-field headline.
The cluster fields keep the D-097 masks: each field's screen and injections share one veto, which is what the
volume integral needs.

**Reported regardless.** Per NCF: status-2 trails, share of used rows vetoed, S_max, p_rot, max |S_×|, recoveries
per cell. Every mask sheet is looked at before any conclusion.

**What changes the plan.** Both flagged NCF admitted and the limit within 2× of Takahashi & Asada at 2 × 10¹³ M☉
→ blank deep fields next (`results/w1_shear_ncf/PLAN.md`). A v2 screen that passes the rule on a previously null
NCF → that field is flagged and vetted like any other.
