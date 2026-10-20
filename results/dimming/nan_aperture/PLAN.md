# W3 dimming: NaN-tolerant aperture in the D-039 persistence test (frozen 2026-10-20 06:15 UTC)

Handoff 1 of `results/dimming/undetected_stamps/RESULT.md`. Written before any re-run with the new aperture.

**Seen before this plan (observed / derived, committed):** 13 / 13 monitored NEXUS `undetected` pairs are visible
sources whose every covering aperture has NaN pixels; 66 / 71 `aperture_no_data` exposures in e4–e7 have ≤ 2 NaN
aperture pixels, all DQ DO_NOT_USE. Measurable held-out ≤ 4 pairs are `on_sky` 33 / 45 (0.59–0.84) against 33 / 89
over all pairs (0.28–0.48).

**Hypothesis.** (H_fill) Measuring through ≤ 2 isolated no-data aperture pixels recovers most of these exposures
without biasing flux, and raises the ≤ 4-exposure `on_sky` efficiency. Ordinary alternatives: the filled flux is
biased (the NaN pixel is often the PSF core, e.g. a saturated or jump-flagged peak), or the gained exposures are
first-of-visit (`n_prior` = 0) or low S/N, so efficiency barely moves.

**Method (ASSUMPTIONs).** `measure_pixel(..., max_nan=k)`: when the position is interior (the unchanged edge rule)
and 1 ≤ n_nan ≤ k aperture pixels are NaN in SCI, each is replaced, in SCI and ERR, by its point reflection through
the rounded centre pixel when that is finite; otherwise the measurement stays NaN. Then the unchanged
`aperture_flux`. Default `k = 0` keeps D-039 exactly. Runs here use `k = 2`. The same rule applies to the
earlier-exposure (illuminator) measurements. `n_ap_nan` is written per exposure.

**Measurement.**
1. G0 (code identity): for every exposure with finite flux under `k = 0` (committed
   `persistence_run/nexus/exposures_e*.ecsv`), the `k = 2` flux is identical (|Δ| = 0) and the verdicts of pairs
   with no NaN exposure are unchanged.
2. G1 (fill accuracy, known cases): `validate.py` range-reads ≤ 150 finite-aperture NEXUS exposures from the e5 run
   (frozen sample: first 150 rows by `file`, `uid` with S/N ≥ 5), masks 1 and 2 non-centre aperture pixels (one random
   pixel, then two random pixels that are not each other's reflection; seed 0) and measures with `k = 2`. Pass: median |f_fill / f_true − 1|
   ≤ 0.05 for both and the S/N ≥ 5 call unchanged in ≥ 95 %.
3. Re-run monitored + targets with `k = 2` for NEXUS, frozen order e5, e6, e7, e4, e3, e1, e0 (held out first), into
   `results/dimming/nan_aperture/nexus/`, until the run budget ends; a partial set is not chosen by outcome.
4. `score.py`: `on_sky` / n by class (≤ 4 / > 4) and split (train e0, e1, e3, e4 / held out e5–e7), monitored
   pairs with mosaic S/N ≥ 10, Wilson 95 %, side by side with the committed `k = 0` numbers.

**Decision rules (ASSUMPTION).**
- D1 (adopt): G0 and G1 pass, held-out ≤ 4 Wilson lower bound > 0.412 (the current swap value), and monitored
  `persistence` stays ≤ 2 % (Wilson lower bound). Then the next PR re-runs the 37 targets in all three fields with
  `k = 2` and sets the ≤ 4 swap efficiency from the new train Wilson lower bound, checked on held out as in
  `expcount_eff` rule 1.
- D2 (keep `k = 0`): any gate fails, or the held-out lower bound ≤ 0.412. The limits keep 0.412.
- Nothing here changes a limit: limits change only in the follow-up PR.

**What would change the plan:** G0 failing (the change touches more than NaN apertures), a filled exposure turning
a monitored pair into `persistence` (the fill invents an illuminator), or G1 bias > 5 % (try area correction instead,
under a new plan).

**Amendment 2026-10-20 06:15 UTC (before G1 or any score):** the rounded centre pixel is its own reflection, so a
NaN centre pixel stays unfilled (the rule above, unchanged). G1 therefore masks non-centre pixels (step 2 edited);
the count of exposures left NaN because the centre pixel is NaN is reported.
