# Changelog

Newest first. Results, failed approaches and the handoff state; not a diary.

## 2026-10-19: W3 ZTF v16: inventory of 25 partitions, frozen queue of 8 (3 pass); 5,995 groups searched, 2 survivors, both read as data problems
- **Inventory** (metadata only, 48.7 min): 54,080 catalogue-eligible unused groups in 25 partitions; 44,354 in the 17
  sampling-qualified ones. The frozen queue of 8 (by unused groups per GB) holds 11,909 searchable.
- **Queue:** 3 of 8 pass the frozen v13 validation (1627, 1624, 1625; inventory sampling worst cell 0.593–0.66). 5 fail on
  short-event cells (worst cell 0.508–0.578). Controls are within limits everywhere (fade ≤ 9.7 %).
- **Search:** 5,995 groups (11,919 series), 2 final survivors.
  - 1624: consistent with a low-response column (official LOW RESPONSIVITY mask, dark column at the source; visual
    reading over 4 automatically undetermined NDs).
  - 1625: consistent with a catalogue miss in poor-seeing images (source visible).
  - Neither is an image-confirmed dimming.
- **Cumulative searched with v13:** 17,773 groups, with no confirmed target event. Groups removed by vetting are not
  image-proven artefacts.
- **Speed:** queue of 8 areas in 14 min, 190 stage-seconds and 239 MB per 1,000 searched groups. Downloads
  3.56 GB (light curves from parquet metadata).
- **Still open:** 4 unqueued partitions with worst cell ≥ 0.59 (2900, 8124, 8126, 8127) hold 17,702 searchable groups
  in 5.76 GB, all sharing fields with earlier samples (`results/w3_ztf/v16/RESULT.md`).

## 2026-10-19: W3 dimming: D-039 persistence-test `on_sky` efficiency on the full control set; NEXUS fails rule 1
- **Why:** #384 handoff. Control-completion addendum frozen first (`results/dimming/persistence_run/PLAN.md`,
  23:12 UTC; NEXUS epochs in a fixed order until the time budget ends).
- **Result (`results/dimming/persistence_run/RESULT.md` "Control completion"):** 438 monitored pairs, `persistence`
  0. MACS0416 and Abell 2744 pass rule 1 (247 / 248 `on_sky`); NEXUS (e0, e1, e3, e4) fails it (129 / 174 at S/N
  20–50). Post hoc: NEXUS pairs with ≤ 4 exposures are `on_sky` 49 / 96, with > 4 94 / 94 (D-052 addendum i).
- **Failed approach:** `pkill -f drive.sh` killed the calling shell again; stop workers by PID (scripts/CLAUDE.md).
- **Handoff:** (1) freeze a per-(source, epoch) efficiency by exposure count and check it on NEXUS e5–e7;
  (2) `forced` + depth / footprint test on the 20 MACS0416 / Abell 2744 `on_sky` targets; (3) the swap in `inject`.

## 2026-10-19: W3 dimming: weight-map coverage term in `screen` / `inject`; corrected headline limits
- **Why:** #382 handoff (decision `over_counts`). Plan and gates R0 / R1 frozen first
  (`results/dimming/wht_coverage/PLAN.md`); `weight_edge` / `mask_weight_edges` / `wht_rel` unit-tested.
- **Result (`results/dimming/wht_coverage/RESULT.md`):** R0 (term off = committed limits) and R1 (term on = the
  validated gate V2 what-if) pass in all three fields. Injection sources 235 → 202; the term removes 10 compact
  unvetoed screen flags and adds none. Corrected headline (`limits_combined.ecsv`, D-052 addendum h): per-source
  limits ×1.07–1.33 weaker at t_E ≥ 0.03 yr, ×0.82–0.92 at 0.01 yr; D-039 swap gain ×1.80–3.52 with the term on.
- **Handoff:** the D-039 swap can now use the corrected chain: measure the persistence test's `on_sky` efficiency
  against mosaic S/N on monitored sources (D-052 addendum f), then fold the swap into `inject`.
## 2026-10-19: W3 ZTF v15: batch with the frozen v13 method; Npix 3837 searched (4,778 groups, 0 vetted survivors); Npix 183 fails validation
- **Choice:** by metadata and cost, Npix 3837 and 183. Npix 8170 has only about 738 unused eligible groups.
  - 3837 shares fields 826 and 1857 with v1–v5, so the sources are new but the exposures are partly shared.
- **Validation** (frozen rule, new seeds, fade × 1,000):
  - **3837 passes:** cells 0.62–0.83, fade 6.6 %, 0 survivors.
  - **183 fails:** 2 short-event cells at 0.48 and 0.455, mainly from sampling (post hoc trace); fade 15.4 %. Its
    1,803 search groups were not searched.
- **Search of 3837:** 4,778 groups (9,519 series), 45 flagged runs in 4 groups, 0 pass vetting, so no image check.
- **Cumulative searched new groups with the v13 rule:** 2,000 (2382) + 5,000 (8170) + 4,778 (3837) = 11,778.
  There is no confirmed dimming. No limit is derived.
- **Fixed:** `configure()` stacked paths on repeated calls (+ test).

## 2026-10-19: W3 dimming edge coverage: gate V2 passed; the coverage proxy over-counts edges (`over_counts`)
- **Why:** #381 handoff 1. Gate V re-frozen (PLAN.md addendum, 21:11 UTC, before any stamp read) on per-epoch
  `forced` stamp cuts of every proxy-covered injection pair (608 pairs, three fields); rules `edge.stamp_edge`,
  `edge.gate_v2` unit-tested first.
- **Result (`results/dimming/edge_coverage/RESULT.md` "Gate V2"):** pooled sensitivity 0.915, false-edge 0.0455 →
  passed; the frozen decision reads `over_counts` (masked veto-off exposure ×0.81–0.87). Post hoc, masking by the stamp
  flags instead gives ×0.84–0.93, still `over_counts`. Disagreement sheet looked at: the grid rule is the conservative one.
- **Deviation:** gate V replaced by V2 after the first gate set was mis-specified (#381); labelled in PLAN and RESULT.
- **Handoff:** add a weight-map coverage term to `build_light_curves` coverage, re-run `inject` for MACS0416,
  Abell 2744 and NEXUS, and quote corrected headline limits (veto-on may drop by up to about 25 % at t_E = 3 yr)
  before using the D-039 swap gain.
## 2026-10-19: W3 ZTF v14: #370 pending checks (official masks), Npix 8170 passes validation, 5,000 groups searched, 1 survivor not a real dimming
- **#370 pending checks:**
  - 26 new cutouts, 17 science and 9 official masks; 656 rebuilt NDs stay pending.
  - The 3 newly checked auto "dimming" NDs sit on mask value 516 (bits 2 LOW RESPONSIVITY + 9 DEAD; names from the
    mask header comments).
  - The 6th group's NDs are undetermined.
  - No item needed priority follow-up.
- **Npix 8170** (fields share nothing with the development data): 500 calibration and 5,000 search groups. The frozen
  v13 rule passes validation:
  - cells 0.57–0.81;
  - fade 50/1,000 = 5.0 % (Wilson 3.8–6.5 %), against 8.5 % (17/200) in v13;
  - 0 calibration flags.
- **Search** of 5,000 groups (9,995 series): 8 flagged runs, **1 vetted survivor** (635216100004733, g band).
  - The image check (7 cutouts) shows the star in all 3 NDs: 2 shallow images and 1 at 0.88 × reference flux.
  - That is a near-limit catalogue drop, **not a real dimming**.
- No limit is derived (`results/w3_ztf/v14/RESULT.md`). W3 / ZTF is paused.
## 2026-10-19: W3 dimming: injection edge-coverage check: gate V failed as written; no decision
- **Why:** #378 handoff 2 (does the injection chain count mosaic-edge positions as covered?). Plan, rule
  (`edge.edge_pair`, F200W `WHT` on a 32-px grid, `wht_rel < 0.5`), gate and decision frozen and unit-tested first
  (`results/dimming/edge_coverage/PLAN.md`).
- **Result (`results/dimming/edge_coverage/RESULT.md`):** gate V needed ≥ 9 of the 11 NEXUS mosaic-edge targets
  called; 5 were, so no decision is read. Post hoc: the gate set was mis-specified (`forced` cut flags are OR'ed
  over epochs the proxy already calls uncovered); all 4 "no other epoch measurable" targets are called. The
  `veto_on` exposures reproduce the committed limits exactly.
- **Post hoc, not a result:** 9–18 % of proxy-covered injection pairs are edge pairs; masking them gives veto-off
  exposure ×0.81–0.87 and veto-on ×0.75–1.22 of the unmasked, the swap gain ×1.80–3.37 (unmasked ×1.99–3.36).
- **Handoff:** re-freeze gate V on per-epoch `forced` stamp flags in proxy-covered epochs (targets and controls,
  both directions), then read the decision; if `over_counts`, add a weight-map coverage term before the swap.
## 2026-10-19: W3 dimming: D-039 persistence test on all 37 targets: 32 `on_sky`, 0 `persistence`
- **Why:** #378 handoff 1, under the frozen `results/dimming/persistence_run/PLAN.md`.
- **Result (`results/dimming/persistence_run/RESULT.md` "Continuation"):** MACS0416 complete and Abell 2744 epoch 1
  done: their 20 targets are all `on_sky`, each only in one epoch per field (MACS0416 `jw01208-o004` CANUCS,
  Abell 2744 `jw02561-o002`), the NEXUS `o014` pattern that `forced` explained as depth / footprint (#378).
  Controls 194 / 194 `on_sky`, 0 `persistence`, all at mosaic S/N 11–33 (fainter unmeasured).
- **Failed approach:** both fields in parallel shared the stream (about 150 MB/s; 1,920 `_cal` × 118 MB); Abell 2744
  epoch 2 (controls only) did not finish. Exposure tables over 1 MB are committed gzipped.
- **Handoff:** (1) `forced` + epoch depth / footprint test on the 20 targets; (2) reconsider the swap, since the
  persistence verdict removed none of the 37 targets.

## 2026-10-19: W3 dimming: NEXUS `o014` targets are mosaic-edge and catalogue-depth effects (0 survivors)
- **Why:** #377 handoff 1. Plan, rule (`summarise.explanation`, unit-tested) and target builder frozen first
  (`results/dimming/nexus_o014/PLAN.md`).
- **Result (`results/dimming/nexus_o014/RESULT.md`):** `forced` on the 12 NEXUS `on_sky` single-epoch vanish flags:
  0 survivors; 5 footprint / image quality, 7 catalogue depth or deblend (forced F200W S/N 4.1–28.2 in covering
  epochs). 11 / 12 lie at or beyond the other epochs' mosaic edges (sheet looked at); `o014` is about 1.5× deeper
  than six of the other seven epochs (control forced errors).
- **Failed requirement:** 8 / 300 controls flagged against the plan's ≤ 1 % (mis-set: the committed NEXUS run had
  10 / 300), so the run is not a calibrated null; the footprint reading uses coverage only.
- **Handoff:** finish MACS0416 e2–3 and Abell 2744 e1–2 persistence; check whether injections count mosaic-edge
  positions as covered before the D-039 swap; set the NEXUS control false-flag requirement.

## 2026-10-19: W3 dimming: D-039 persistence test on the single-epoch vanish flags (partial)
- **Why:** #376 handoff. Plan and rules frozen first (`results/dimming/persistence_run/PLAN.md`, rules as tested
  functions in `summarise.py`).
- **Result (`results/dimming/persistence_run/RESULT.md`):** NEXUS 17 targets (all in epoch `o014`): 12 `on_sky`,
  4 `inconclusive`, 1 `undetected`, 0 `persistence`. Controls (partial): 61 / 61 monitored pairs `on_sky`, 0
  `persistence` (S/N 10–50). No swap decision from the partial set; no flag is interpreted.
- **Deviations / failed approaches:** NEXUS ran without controls (I/O: each `_cal` open streams several MB);
  MACS0416 and Abell 2744 unfinished (their 20 targets untested); `pkill -f` killed the calling shell.
- **Handoff:** (1) `forced` on the 12 NEXUS `on_sky` targets, `o014` depth/footprint test first (its F200W
  catalogue is 65 MB against 12–27 MB); (2) finish MACS0416 e2–3, Abell 2744 e1–2; (3) then the swap.

## 2026-10-19: W3 dimming: the D-039 persistence test in place of `single_epoch` is worth building (scoping)
- **Why:** #375 handoff, D-052 "Revisit if". Plan and decision rule frozen first
  (`results/dimming/persistence_scope/PLAN.md`); screen code and thresholds unchanged.
- **Result (`results/dimming/persistence_scope/RESULT.md`):** with the veto on, the committed limits of all three
  fields are reproduced exactly. Workload: 37 compact single-epoch vanish flags (12 / 8 / 17 in MACS0416 /
  Abell 2744 / NEXUS). Value (simulated what-if, veto off, ρ = 0.1): joint exposure ×2.0–3.4 by t_E
  (Abell 2744 ×17–53); 100 % dimming recovered 240 → 680, 20 → 1,238, 1,140 → 2,779. Rule met: build the swap.
- **Limit of the reading:** an upper bound; it assumes the persistence test calls every real full vanish `on_sky`
  (two clean exposure detections at S/N ≥ 5), which the faint bins will not reach.
- **Failed route:** VizieR timed out; the Gaia mask came from the ESA Gaia TAP (same cone; reproduction exact).
- **Handoff:** (1) persistence `on_sky` efficiency vs mosaic S/N on monitored sources; (2) `persistence_check.py`
  on the 37 flags; (3) then swap the veto in `efficiency_table` and re-run `combine`.

## 2026-10-19: W3 dimming: chromatic veto folded into `forced`; weekly network test failure is upstream
- **Change:** `dimming_screen.py forced` applies the D-052 addendum c chromatic veto (simultaneous F444W below
  F200W / 20 at 3σ in some epoch) before `survives_artefact_tests` (D-052 addendum e, unit-tested). No field
  re-run: injections do not pass through `forced` and an achromatic W3 event cannot meet the veto, so the joint
  headline stands as written.
- **Issue #374 (weekly network tests):** 3 of 1,604 failed, all connection errors (MAST `RemoteDisconnected` twice, NED
  `ConnectTimeout`); the run took 24 min against 5 min on 2026-10-12. MAST API timeouts were also seen in the
  13:12Z run. Treated as upstream; no code change. Re-running the job from the cloud session is refused (403).
- **Handoff:** D-039 persistence test in place of the single-epoch veto for vanish flags; JADES / new NEXUS epochs.

## 2026-10-19: W3 dimming: MACS0416 regenerated (reproduces D-052); joint three-field `combine` gives the per-deg² headline
- **Why:** #371 / #369 handoff. Plan and a computable reproduction gate frozen first
  (`results/dimming/macs0416_regen/PLAN.md`, `scripts/dimming_repro_gate.py` + unit test).
- **Result (`results/dimming/macs0416_regen/RESULT.md`):** screen identical to D-052 (854 / 358 / 10 compact);
  forced calibrated, 1 confirmed (the saturated star d03349, cutout tests remove it), 0 survivors, 0 / 300 controls
  flagged; exposures −3.7 % to +0.7 % of D-052 (gate passed).
- **Joint headline** (`limits_combined.ecsv`, ρ = 0.1): < 0.44 / 0.11 / 0.043 / 0.029 / 0.021 / 0.016 per source per yr
  and < 3.0e3 / 776 / 295 / 201 / 157 / 118 per deg² per yr at t_E = 0.01 / 0.03 / 0.1 / 0.3 / 1 / 3 yr.
- **Failed approach / deviation:** `fetch` via the MAST API timed out (3 tries); catalogues read from S3 `cloud_uri`,
  all 9 sha256-identical to the committed manifest. Forced noise scale differs from D-052's offline re-calibration
  (cause not investigated).
- **Handoff:** fold `chromatic_veto` into `forced` and re-run the three fields' injections; D-039 persistence test
  for vanish flags.

## 2026-10-19: W3 ZTF v13: edge-duration rule passes validation; 2,000 groups searched, 0 vetted survivors, no image-confirmed dimming
- **New rule:** one W3 vetting rule replaced. The intermediate-point count becomes edge duration in days, allowed up to
  0.25 × the run span. The 0.25 comes from the W3 model only (99th percentile at ρ 0.1, blend 0.7). It was frozen
  before any comparison (`results/w3_ztf/v13/`).
- **Development set** (same events and noise; not a validation):
  - recovered 1,355 → 1,517 of 2,000;
  - the failing cell 84 → 148;
  - 166 of 202 sharp-edge losses kept, 4 earlier recoveries lost;
  - fade control 4 → 9 of 200.
- **Validation** on 500 never-read groups with new seeds: **passed**.
  - cells 0.64–0.87;
  - fade 17/200 (8.5 %, limit 10 %), other controls 0;
  - 0 real survivors.
  - Same area and exposures as the development set, so not an independent-area validation.
- **Search of 2,000 groups:** 96 flagged runs in 6 groups, **0 pass vetting**. Image check of 5 groups (50 cutouts):
  all 17 "credible dimming" non-detections have a dark detector column within 2″ of the source, against 0 of 20
  detections (post hoc check). That is consistent with a data problem. One group is pending. **No image-confirmed dimming.**
- **Handoff:** W3 / ZTF is paused. The fade-control margin (8.5 % vs 10 %) is the known weak point of the new rule.
## 2026-10-19: W3 dimming: NEXUS survivors vetted at exposure level; 0 W3 events, three-field headline stands
- **Why:** #369 handoff (its headline was conditional on 5 by-eye survivors). Plan frozen first
  (`results/dimming/nexus_vet/PLAN.md`).
- **New:** `scripts/dimming_exposure_check.py`: fixed-aperture flux, JUMP / DO_NOT_USE counts and stamps on every
  level-2 `_cal` exposure covering a flagged point; `chromatic_veto` (post hoc, unit-tested).
- **Result (`results/dimming/nexus_vet/RESULT.md`, `vetting.json`):**
  - d53823 is in one of two dithers (93.6 µJy, 60 JUMP pixels; 0.0003 in the other).
  - d60156 and d79928 are covered by one exposure each, are JUMP-flagged, and fail the chromatic veto
    (simultaneous F444W < F200W / 20 at 3σ).
  - d73390's F200W variation is not reproduced on `_cal` (min / max 0.865 against 2.03 in the i2d forced
    photometry).
  - d79956 is a brightening in both bands, a real sky source (SIMBAD none; TNS 403, NED timeout).
  - Controls show the same class: 3 covered single-epoch control points are single-exposure, JUMP-flagged events.
- **Consequence:** 0 W3 events; the three-field headline holds: 0.020 per source per yr at t_E = 1 yr and 0.029 at
  0.3 yr (docs/exotic_limits.md, D-052 addendum c).
- **Deviations:** the chromatic veto is post hoc. The nearby-galaxy comparison for d73390 was not run. 6 faint flagged
  controls were not examined at exposure level.
- **Handoff:** MACS0416 `forced --workers` regeneration, then a joint `combine` for the per-deg² headline. Fold
  `chromatic_veto` into `dimming_screen.py forced` and re-run injections before quoting a limit that uses it in
  the chain.

## 2026-10-19: W3 dimming: NEXUS forced stage calibrated; 5 artefact-test survivors, none a W3 dimming by eye; headline ×2.6–6.2 tighter if they hold
- **Why:** handoff of #367 (NEXUS was the one uncalibrated field that could move the headline). Plan frozen first
  (`results/dimming/nexus_calib/PLAN.md`).
- **New:** `dimming_screen.py forced --workers N` reads the (band, epoch) images concurrently. Validated before use:
  Abell 2744 summary identical, 666 → 236 s (`workers_validation.json`).
- **Result (`results/dimming/nexus_calib/RESULT.md`):** screen reproduces D-052 (3,793 / 1,051); `forced` calibrated
  (noise 1.91 / 2.57, inflation 1.00 / 2.06); 32 confirmed, 5 survive the artefact tests, 10 / 300 controls flagged.
  Survivors looked at (`survivors_sheet.png`): three single-band single-epoch F200W points, an extended-nucleus
  F200W-only variation (WISE source), an SN-like brightening in both bands. None is a W3 dimming (by eye, post hoc).
- **Limits:** W3 efficiency 0.026–0.218 (ρ = 0.1). Three-field headline, conditional on the survivors: 0.020 per
  source per yr at t_E = 1 yr (was 0.053), 0.029 at 0.3 yr (was 0.14) (`headline.json`).
- **Deviation:** `/vet-candidate` on the survivors did not run in this cycle (time); flagged controls not examined.
- **Handoff:** vet the 5 survivors (TNS for d79956) and look at the 10 flagged controls; then MACS0416 `forced
  --workers` and `combine`; update docs/exotic_limits.md once the survivors are vetted.

## 2026-10-19: W3 ZTF v12: new-field batch area Npix 2382 fails validation (ρ 0.1, blend 0.7, t_E 30 d at 0.42); search not executed
- **Setup:** plan, dedupe list (4,072 earlier groups) and candidates were fixed before any read
  (`results/w3_ztf/v12/`).
- **Area check (metadata):** of 5 candidates, 4 qualify. 2382 (0.6965) and 8170 (0.6615) share no field with
  earlier samples. Their partitions need 4.4 GB together, so the budget rule kept 2382 only.
- **Npix 2382:** 3,000 new groups, 500 calibration and 2,500 search. 9 of 10 cells reach ≥ 0.5, but
  ρ 0.1 / blend 0.7 / t_E 30 d recovers 84/200.
  - Of the 87 flagged injections lost in that cell, 80 fail first on the frozen sharp-edge test (post hoc trace).
    Sampling is not the limit: 171/200 are flagged.
  - Controls ≤ 2 %; 0 calibration survivors.
- **Not executed:** the 2,500 search groups. 8170 was not read. No threshold changed and no retry.
- **Handoff:** W3 / ZTF is paused again.
  - In the one validated area, the failing cell lost most flagged injections at the sharp-edge test (80 of 87, post
    hoc). A link to cadence is a hypothesis; it was not tested.
  - Another area under the same rules may pass or fail; it would not change the method. Changing the test is a method
    change that needs the owner.
## 2026-10-19: W3 dimming: Abell 2744 forced stage calibrated (300 controls); still null; headline 5–11 % tighter
- **Why:** D-052 "Revisit if" holds (S3 byte-range jobs have worked since 10-08); Abell 2744 and NEXUS were never
  calibrated. Plan frozen first (`results/dimming/abell2744_calib/PLAN.md`); code unchanged.
- **Result (`results/dimming/abell2744_calib/RESULT.md`):** screen reproduces D-052 (530 / 32 flags); `forced`
  calibrated on the o002–o006 pair (247 / 256 controls); 0 of 32 flags confirmed, 0 controls flagged. Noise scale
  4.05 / 7.21 and F444W inflation 3.60 (higher than MACS0416; unexplained).
- **Limits:** calibrated W3 efficiency 0.000–0.013 (the uncalibrated values were 4–7× too high); headline
  per-source rate 0.053 per source per yr at t_E = 1 yr (was 0.056), 5–11 % tighter over t_E.
- **Handoff:** NEXUS `forced` with ≥ 300 controls next (the field that could move the headline), then a joint
  `combine` for the per-deg² headline.

## 2026-10-19: W3 ZTF v11: area chosen by observing conditions (Npix 3839) passes validation; 488 search groups, 0 flags
- **Records fixed first** (`results/w3_ztf/v10/RESULT.md` addendum):
  - v10's 711 vs v6's 703 differs only in the noise realisation (per-expid seeds) on the same events; the re-read inputs have the same
    row counts.
  - +161 compares the same events, code and noise with and without the psfcat fill.
  - The direction of the 261 dropped points is not verified.
- **Area choice:**
  - Metadata only (`results/w3_ztf/v11/obs_conditions.json`): the worst-cell necessary-sampling fraction is 0.687,
    0.6855 and 0.7025 for Npix 3837, 3838 and 3839, with seasonal gaps kept.
  - 3839 was chosen. Its new 1,000-group sample excludes the v1–v4 samples.
  - The series hold a median of 98 % of the official exposures.
- **Validation passed** (frozen v6 rule): cells 0.615–0.855; controls ≤ 4 %; 0 calibration flags.
- **Search:** one search of the 488 groups found **0 flagged runs**, so no image check was needed. This is not an
  independent-area validation (same fields as v1–v5), and no limit is derived.
- **Handoff:** the W3 / ZTF DR24 unit is paused. No follow-on batch without the owner.

## 2026-10-19: TASKS sync: two stale "Next" items were already done
- The halo-phot elongation gate and its `fluxratios --halo-phot --mu-samples 100` run (D-032 addenda 2026-10-13) and
  the no-injection control of the off-arc vetoed positions (`results/transients/veto_control/RESULT.md`, 2026-10-16:
  mixed) still read as open in TASKS.md; now struck through with their results.
- **Handoff:** every open lensing line waits on an owner decision (off-arc efficiency, #199, #237, #279, X4 budget)
  or is optional low-value (VENUS/Bullet `arcs` without photo-z, critical-curve curvature). No new unit started.

## 2026-10-19: Review-list h039 (SMACS) and h005 (Sunrise o052): frozen rule inconclusive; post hoc within control scatter
- **Why:** the two open by-eye artifact calls of `results/transients/review_list/RESULT.md`. Plan and rule frozen
  first (`results/transients/review_list/vet_h039_h005/PLAN.md`).
- **New:** `scripts/review_pair_rule.py` (+ tests; written after the run, implementing the frozen text: deviation).
- **Result (`vet_h039_h005/RESULT.md`, `outcome.json`):** both **inconclusive** (one band ≥ 2σ each): h039 F150W
  +37.2 ± 14.7 nJy (2.53σ), F444W 1.94σ; h005 F444W −46.2 ± 6.6 nJy (−7.01σ), F150W −0.59σ.
- **Post hoc:** 40 matched compact controls per entry: |z_cal| ≤ 1.17, control tail fractions ≥ 0.375, so neither
  changes beyond control scatter. Stamps looked at: no visible change. Closed post hoc; not candidates; no lensing
  prediction tested (both off the arc mask and catalogued images).
- **Handoff:** the review list has no open entries left. The arc-box line still waits on the owner's off-arc
  efficiency choice (TASKS, Two-epoch search).

## 2026-10-19: AT 2024jfz z_s 1.55–1.70: PSF-matched difference untestable (star core residual); post hoc no counter-image beyond 0.35″
- **Why:** handoff of #361 (the z_s 1.55–1.70 rows of #360 were undetermined: a cross-instrument residual of a
  static source). Plan and rule frozen first (`results/transients/smacs_at2024jfz/psfmatch/PLAN.md`);
  `outcome_psfmatch` unit-tested before the run.
- **New:** `scripts/smacs_at2024jfz_psfmatch.py` (+ tests): T resampled to E's grid, `arc_diff.kernel_fit` kernel
  with the predictions and TNS masked, D = E − K * T, injection through the whole chain.
- **Result (`psfmatch/RESULT.md`):** frozen outcome **untestable** in all four rows: control passes (D at TNS
  228.6 nJy), injections pass, but the subtraction gate fails (D −131.9 nJy at the source core; thresholds
  96.9–99.7 nJy).
- **Post hoc:** the source is star-like (spikes in T, by eye); outside 0.35″ of it the disc maxima are 9.9–45.8 nJy,
  "absent" under the same rule in every row. A counter-image within ≈ 0.35″ of the star stays untested.
- **Handoff:** the counter-image line of AT 2024jfz is done to the depth these epochs allow (z_s 1.35–2.50 under the
  best ICLv2 model; z_s 1.55–1.70 only post hoc). Remaining deciders: a spectrum at the position or a new epoch.
  Not a candidate.

## 2026-10-19: AT 2024jfz host redshift: no catalogued host within 1″; frozen outcome untested
- **Why:** handoff of #360 (a host z decides which counter-image window applies). Plan and rule frozen first
  (`results/transients/smacs_at2024jfz/host_z/PLAN.md`); `classify` / `pick_host` unit-tested.
- **Result (`host_z/RESULT.md`):** no DJA v7.4 source within 1.0″ and no Caminha+22 MUSE source within 3″ of the TNS
  position: outcome **untested**, host z unknown.
- **Post hoc:** the nearest DJA source (1.632″, d / r_half 7.79) is image 6.1 of system 6, eazy z 1.597–1.634, inside
  the undetermined window 1.55–1.70. Faint template flux at the position (31 nJy, #360) sits on a bright vertical
  band, so a blended host is not excluded.
- **Handoff:** if pursued, a PSF-matched difference at the z_s 1.55–1.70 track or a spectrum at the position. Not a
  candidate.
## 2026-10-19: W3 ZTF v10: psfcat in the W3 chain raises recovery by 161 of 2,000 same events; 8 of 10 cells stay below 0.5; search not executed
- Plan and fill set fixed before any psfcat read (`results/w3_ztf/v10/PLAN.md`). The fill set has 724 exposures
  where DR24 has no detection of any control star; 32 of 33 field / quadrant / filter combinations pass the v9
  overlap test (672_31_2 has 3 pairs and was not filled). In total, 790 files (183.6 MB) were downloaded once.
- The merge adds 32,643 r-band and 13,103 g-band epochs to the calibration series. DR24 is never overwritten, and a
  missing psfcat row stays unknown.
- Same v6 events (exact replay): recovered 711 → 872 (v6 stored 703). The main gain is at the MIN_BRACKET stage
  (289 → 171). Only 2 of 10 cells reach 0.5 (124 and 102 of 200).
- Controls (fade 16/200) and burden (0 survivors) pass, but the frozen rule fails, so the **539-group search was not
  executed**.
- Of the 1,128 events still lost:
  - 602 have no dim run in the umbra, 527 of them with fewer than 3 official exposures there;
  - 86 lack official exposures for the bracket;
  - 128 fail vetting;
  - 85 have ≥ 20 official r-band exposures on both sides, but fewer than 20 non-dim series epochs on a side.
- **Handoff:** paused. No threshold change, no new sample, no follow-on batch without the owner
  (`results/w3_ztf/v10/RESULT.md`).

## 2026-10-19: AT 2024jfz: no counter-image in NIRISS 4368 F150W at z_s 1.35–1.40; the transient is still visible at +116 d
- **Why:** handoff of #358 (the only archival epoch matching the counter-image timing, if z_s ≈ 1.35–1.40). Plan
  frozen first (`results/transients/smacs_at2024jfz/niriss4368/PLAN.md`); rule `outcome` unit-tested before the run.
- **New:** `scripts/smacs_at2024jfz_niriss.py` (+ tests): NIRISS 4368 epoch vs the ERO NIRCam F150W template, S3
  byte-range boxes, injection control.
- **Result (derived, `RESULT.md`):** frozen outcome **absent**: max Δ 9.5 nJy (1.0σ) over 119 grid points within
  0.5″ of the track, against T_ap 271 nJy; injections 6 / 6 recovered. AT 2024jfz itself is still in E at 254.7 ± 5.6
  nJy (+116.1 d; template 31.3 ± 3.2).
- **Post hoc:** under a single-peaked light curve, nearby counter-images are not seen for z_s 1.35–1.50 and
  1.75–2.50 (21 rows, max Δ ≤ 40.9 vs thresholds ≥ 93.5 nJy; expected counter-images injected and recovered in all
  25 rows). z_s 1.55–1.70 is undetermined (tested rule `classify_flag`): the flags sit on the wing of a static compact
  source (centroids agree to 0.014″; flux change fraction 0.0266 across instruments).
- **Review fixes:** post hoc injection gate, `classify_flag` replaces a by-eye reclassification, off-footprint and
  NaN rows return "untestable", the epoch is a checked constant.
- **Handoff:** the AT 2024jfz counter-image line is closed except z_s 1.55–1.70 (needs a same-instrument epoch or
  a PSF-matched difference) and z_s outside 1.35–2.5. A host redshift would decide which applies. Not a candidate.

## 2026-10-19: AT 2024jfz (SMACS 0723): counter-images not testable in the JWST epochs for z_s ≤ 2.5 under Mahler+22 ICLv2
- **Why:** open TASKS item (SMACS arc `arc_diff`, h000 = AT 2024jfz, 1.81″ from image 6.1): does the cluster model
  predict counter-images that existing epochs could test? Plan frozen first (`results/transients/smacs_at2024jfz/PLAN.md`).
- **New:** `scripts/smacs_at2024jfz.py` (+ tests): images by `lensmodel.find_images`, Fermat delays from line integrals
  of the analytic deflection (SIS and path-independence checks pass; real-model path difference 1.5 × 10⁻¹³), frozen
  rule `testable`.
- **Result (model_prediction, `RESULT.md`):** frozen outcome "not testable" in the four NIRCam epochs for z_s 0.6–2.5;
  "testable" at z_s 4.0 only as the unresolved 0.087″ merging pair at the 2024 epoch itself (source on the critical
  curve; no resolved-pair measurement done). At system 6's model z (1.71) the counter-image 4.5″ away (1.88″ from 6.2)
  trails by +31.0 d (μ ratio −0.89), the far one led by 5,464 d. Best model only (no chain).
- **Review fixes:** next-window formula mirrored in PLAN.md, fixed (deviation); the path gate now withholds the
  outcome; delay factor tested against #308; epoch list and separations written by the script.
- **Post hoc:** NIRISS 4368 F150W (MJD 60558.0) covers the near counter-image's arrival only if z_s ≈ 1.35–1.40; HST
  WFPC2 11103 matches the far image's timing only if z_s ≈ 1.75–1.80 (depth not checked); footprints not checked.
- **Handoff:** optional next step: footprint check plus forced photometry of the NIRISS 4368 epoch at the
  z_s 1.35–1.40 positions (frozen plan first); a host redshift would decide which z_s applies. Not a candidate.

## 2026-10-19: E1 BAT–GBM control deficit: long lags behave as the shared-event displacement predicts; frozen rule undetermined
- **Why:** handoff of #356 (optional: the BAT–GBM deficit, D-074 addendum 6, with the shared-event test of
  addendum 9). Plan frozen first (`results/e1_swift/selfpair/PLAN.md`); rule `outcome` unit-tested before the run.
- **Result (derived, `results/e1_swift/selfpair/RESULT.md`):** frozen outcome "undetermined". Without the 476
  BAT bursts that have a GBM burst within 100 s: all-separation z −0.33 / −0.66 (1 h–1 d, 1–7 d), from
  −4.78 / −5.04; random drop keeps −5.83 / −5.31. The rule fails on 100 s–1 h (z +4.2): 25 same-direction pairs,
  23 within 0.1° (e.g. GRB 221009A at −3,197 s, GRB 130925A), read as joint GRBs with triggers 100 s–1 h apart.
  **Post hoc** 1 h variant reads "explained" but is circular in that cell (deviation, exploratory).
- **Review fixes:** GBM scrambled first so variants are paired; lag-bin windows only; the 1 h pair listing
  stored (`joint_1h_only_pairs`); wording narrowed. Numbers are from the re-run.
- **Handoff:** E1 waits on the CHIME uptime series (owner). Optional: a frozen joint-event rule (position match
  over > 100 s) for EP–GBM (addendum 7) and BAT–GBM.

## 2026-10-19: E1 DR2–ICECAT control deficit explained: the null displaces the 276 shared events
- **Why:** handoff of #354 (DR2–ICECAT wide deficit at 1 h–7 d, z −4.5 / −6.5, cause open). Plan frozen first
  (`results/e1_icetracks/selfpair/PLAN.md`); rule `outcome` unit-tested before the run.
- **Result (derived, `results/e1_icetracks/selfpair/RESULT.md`):** frozen outcome "explained". Without the
  276 DR2 tracks that are ICECAT alerts' own events: wide z −0.09 / −1.24 / +1.01 (100 s–1 h, 1 h–1 d, 1–7 d);
  without 276 random tracks: −1.27 / −4.63 / −6.98. The null mean falls by 122.38 / 450.22 pairs, the data by
  60 / 240; the null-only parts (62.38 / 210.22) are within 22 % / 11 % of the planned displacement estimate
  (≈ 80 / 190). D-074 addendum 9.
- **Consequence:** with shared events removed the jit null matches this same-detector channel (|z| ≤ 1.24);
  tested channels (no shared events) are unchanged; in long-lag wide cells shared events raise the null mean,
  lowering sensitivity to an excess there. The BAT–GBM deficit (addendum 6) is consistent with the same
  mechanism (not recomputed).
- **Review fixes:** z for zero-spread cells via `en.z_score`, scrambles chunked over all cores with an
  initializer; numbers above are from the re-run (first run: same outcome, not stored).
- **Handoff:** E1 waits on the CHIME uptime series (owner). Optional: recompute the BAT–GBM deficit with the
  expected displacement; the EP–GBM joint-event rule (addendum 7).

## 2026-10-18: E1 + IceTracks-DR2 (seventh node set, northern tracks): 35-cell family null; control passes
- **Why:** handoff of #353 (last open E1 node set). Plan frozen first (`results/e1_icetracks/PLAN.md`); data are
  the D-079 pinned DR2 files (fetched in 13 s, reduced cache only), other catalogues from `e1_swift.ecsv`.
- **Result (derived, `results/e1_icetracks/RESULT.md`):** 7,583 tracks (Dec > −5°, log10 E ≥ 4). Control:
  270 DR2–ICECAT same-direction pairs at 0–10 s, one-to-one, vs 0.0013. Family min DR2–CHIME 1–7 d same (35 vs
  26.51, p 0.065, × 35 = 1); pooled global p 0.79: **null**. DR2–GBM same at 0–100 s: 0 vs 0.47.
- **Post hoc:** 6 same-neutrino DR2–ICECAT pairs at 0–10 s fall in "wide" (reconstruction differences). Review
  found a DR2–ICECAT wide deficit at 1 h–7 d (z −4.5, −6.5); six distinct null variants leave it (not downtime, DR2 end,
  or one side's scramble); cause open. Review also fixed seed reproducibility (set order), so all numbers are
  from the re-run.
- **Handoff:** E1 node sets are complete; E1 waits on the CHIME uptime series (owner). Optional: the DR2–ICECAT deficit (e.g. lag
  curve around alert times, alerts split by realtime vs archival); f_astro per cut (D-079 revisit).

## 2026-10-18: E1 + Einstein Probe (sixth node set, from GCN Circulars): 18-cell family null; control passes
- **Why:** handoff of #351 (E1 open node sets: Einstein Probe). No public EP event table exists, so the table is
  parsed from the GCN Circulars JSON archive (`scripts/e1_ep.py build`, rules unit-tested; 5 / 5 random rows read
  against their circulars). Plan frozen first (`results/e1_ep/PLAN.md`). ICECAT-1, CHIME Cat 2 and the pinned GWTC
  end before EP starts, so only GBM and BAT channels run.
- **Result (derived, `results/e1_ep/RESULT.md`):** 231 EP transients; control 27 EP–GRB same-direction pairs at
  0–1 h vs 0.63. Family min EP–GBM 1 h–1 d same (27 vs 14.53, p 2.2 × 10⁻³, × 18 = 0.039), pooled global p 0.017:
  **null**. EP–EP: 95 % limits 1.3 × 10⁻² to 9.7 × 10⁻² dependent wide pairs per EP transient.
- **Post hoc / deviation:** 5 of the 27 min-cell pairs are one source twice (EP250702a = GRB 250702B–E; EP250919a,
  EP261007a). D-074 addendum 6 asks to drop joint detections first; the post hoc re-count does (26 vs 12.85,
  × 18 = 0.015), still above the rule. `jitday` dropped (equals `jit` here). Failed approaches: `control_sum`
  mixed axes on the null ensemble; the parser took observation starts as onsets; both fixed and re-run.
- **Handoff:** E1 open node set: IceTracks-DR2. Optional: a joint-event rule by name/position (not |lag| < 1 h)
  for EP–GBM, and a time-convention check of EP onsets vs GBM for the long-lag same cells.

## 2026-10-18: E1 + Swift BAT GRBs (fifth node set): 35-cell family null; BAT–GBM control passes
- **Why:** TASKS Now 0 (E1) lists Swift as an open node set; every other open line waits on the owner or is
  optional (2026-10-18 15:08Z queue check). Plan frozen first (`results/e1_swift/PLAN.md`); rules
  `control_passes` / `outcome` unit-tested before the run (`scripts/e1_swift.py`, `tests/test_e1_swift.py`).
- **Data:** HEASARC `swiftgrbba` (current to 2026-10-07; `swiftgrb` stops in 2012), 1,635 on-board BAT triggers;
  the four D-074 catalogues re-pinned in `data/manifests/e1_swift.ecsv` (the live `fermigbrst` moved since D-074).
- **Result (derived, `results/e1_swift/RESULT.md`):** control 466 BAT–GBM same-direction pairs at 0–100 s vs 0.10
  expected (passes). Family: min cell BAT–CHIME 10–100 s wide, analytic p 0.047 (× 35 = 1), pooled global p 0.54:
  **null**. 95 % limits per anchor, e.g. ≤ 1.8 × 10⁻³ dependent wide BAT–BAT pairs per burst at ≤ 100 s,
  ≤ 3.3 × 10⁻³ BAT–ICECAT at ≤ 100 s (`limits.ecsv`).
- **Post hoc (control channel only):** 8 same-GRB pairs land in "wide" through the 0.05° σ floor; a BAT–GBM deficit
  at 1 h–4 d (1 h–1 d z −4.9) is in all 19 years, survives a Sun-frame null and vanishes when the 466 joint bursts
  are removed (−0.5σ); its mechanism is not determined (`scripts/e1_swift_diag.py`, four JSONs).
- **Handoff:** E1 open node sets: IceTracks-DR2 (1.6 M events; would need a ghost-pair-style reduction, D-079),
  Einstein Probe (no public event table found yet; not searched this run).

## 2026-10-18: halo-subtracted flux ratios on Abell 370 (CANUCS, per-system rms gate): 2 images testable, no anomaly
- **Why:** fifth cluster for the D-032 flux-ratio line (mission: flux ratios against published models). The four
  usable Lenstool models (SMACS, El Gordo, Abell 2744, MACS0416) were done; `abell370-canucs`' list is gated off
  (D-044). Plan frozen first (`results/halo_phot/abell370/PLAN.md`, with the pre-run deviation).
- **New:** `lens_consistency.systems_passing_rms` / `gated_systems` and `--gate-imageplane` / `--gate-max-rms` on
  `halo_phot.py` and `fluxratios` (+ tests): a gated list opens only for systems with every image solved and
  image-plane rms ≤ 0.5″ (ASSUMPTION). Abell 370: 27 / 45 systems, 64 / 115 images.
- **Result (derived):** 2 images tested (system 9; F150W and F444W), residuals 0.09–0.23 mag, no `flux_outlier`,
  no colour outlier; `flux_tension` untested (too few images). Injections (simulated): a −1 mag change recovered
  in 1 / 2 targets, +1 mag in 0 / 2. Stamps looked at. **No anomaly; almost no power** (D-032 addendum).
- **Handoff:** the halo-flux line now covers five clusters; no further cluster has a usable Lenstool par + list.
  More power needs a gate that keeps resolved images (a surface-brightness test, D-032) or new image-plane models.

## 2026-10-18: `arcs` on MACS J0308.9+2645 (RELICS maps, RELICS BPZ): the S1063 anti excess does not replicate (1 vs 0.31, P 0.27)
- **Why:** next field of the S1063 replication tally (handoff of #346). Plan frozen first
  (`results/arcs/macs0308-relics/PLAN.md`); F200W catalogue from program 5293 o005 (1,578 s; chosen over 6882 o015 by
  coverage only, 705 vs 640 S/N ≥ 10 rows within 50″); photo-z RELICS BPZ through the SPT0615 conversion
  (95 % interval as z160–z840, ASSUMPTION; JWST − BPZ offset +0.34″/−0.10″ removed). z_lens 0.356 (Acebron et al.
  2018, Table 1). Exploratory.
- **New:** models `macs0308-relics-glafic` (v2) / `macs0308-relics-zitrin` (LTM-Gauss v1), both pass `validate`;
  `scripts/arcs_macs0308.py` (+ tests, reuses `arcs_spt0615`), `scripts/arcs_macs0308_sheet.py` (F410M as second
  band; the program has no F444W).
- **Result (derived; `RESULT.md`):** GLAFIC 17 / 21 aligned (control passes), anti 1 vs 0.309, P 0.266 → **frozen
  outcome "does not replicate"**; Zitrin 15 / 18, anti 1 vs 0.39, P 0.32 (same source). The anti (2157, m 20.6,
  z_phot 0.76) lies where the predicted axis flips between z = 1 and 2 (82° vs 8° offset); its class rests on BPZ
  z < 1. Sheets looked at. Not a candidate.
- **Handoff:** S1063 tally: MACS0417 no, MACS0647 no result, MACS1423 explained, RXJ2129 no result, WHL0137 no,
  RXC J0600 no, SPT0615 no result, MACS0308 no. Four "no", one explained, three "no result": nothing supports the
  S1063 excess beyond S1063 itself. Remaining routes are optional (VENUS clusters without photo-z; Bullet).

## 2026-10-18: `arcs` on SPT-CL J0615−5746 (RELICS maps, RELICS BPZ photo-z): no result, both controls fail
- **Why:** one more independent field for the S1063 anti-excess tally. Plan frozen first
  (`results/arcs/spt0615-relics/PLAN.md`, which also records the field survey below). F200W catalogue from program
  4212 o001; photo-z from the RELICS BPZ catalogue (no DJA table exists; 95 % interval used as z160–z840,
  ASSUMPTION; JWST − HST offset +0.31″/−0.27″ removed). Exploratory.
- **New:** models `spt0615-relics-glafic` (v2) / `spt0615-relics-lenstool` (v1), both pass `validate`;
  `scripts/arcs_spt0615.py` (+ tests), `scripts/arcs_spt0615_sheet.py`.
- **Result (derived; `RESULT.md`):** strong-shear background rows 7 (GLAFIC: 3 aligned, 2 anti) and 4 (Lenstool:
  1 aligned, 3 anti). Both controls fail (43 %, 25 % < 70 %) → **frozen outcome "no result"**; the P values (0.0065,
  0.00034) decide nothing. Sheets looked at: Lenstool anti 2681 sits at g_min 4.84 (critical-curve flip region; 0.67
  in GLAFIC); 2669 has F200W striping across it; 2773 is faint. No candidate.
- **Field survey (so it is not repeated):** DJA v7 photo-z exists for only two unrun lensing clusters
  (RCS2 J0327, MACS J0138), neither with public maps found. The LAM Strong Lensing Cluster Atlas
  (data.lam.fr/sl-cluster-atlas; ANIS API `server/search/sl-cluster-atlas/{model,map_model}` with cookie
  `Yogsototh_opens_the_door=1`) serves 17 Lenstool model tarballs; of those clusters only the Bullet has deep JWST
  F200W (program 4598, no public photo-z found); VENUS (6882) visits are 1,245 s. SGAS HLSP models give only κ and
  |γ| (no shear direction).
- **Handoff:** S1063 tally: MACS0417 no, MACS0647 no result, MACS1423 explained, RXJ2129 no result, WHL0137 no,
  RXC J0600 no, SPT0615 no result. Remaining routes: MACS0308 (RELICS GLAFIC + Zitrin maps, RELICS BPZ, JWST ≤
  1,578 s); Bullet (LAM Lenstool model) needs a photo-z source.

## 2026-10-18: `arcs` on RXC J0600−2007, RELICS maps: the S1063 anti excess does not replicate (1 vs 0.20, P 0.18)
- **Why:** the S1063 replication tally needs fields never used to tune `arcs`. Correction to the 2026-10-17 handoff:
  RXC J0600 does have public maps (RELICS `rxc0600m20`: GLAFIC v3, Lenstool v1). Plan frozen first
  (`results/arcs/rxc0600-relics/PLAN.md`); F200W catalogue from program 6796 o010, DJA v7.2 photo-z.
- **New:** models `rxc0600-relics-glafic` / `rxc0600-relics-lenstool` (sha256-pinned; both pass `validate`);
  `scripts/arcs_rxc0600.py`, `scripts/arcs_rxc0600_sheet.py`.
- **Result (derived; `RESULT.md`):** GLAFIC 11 / 12 aligned (control passes), anti 1 vs 0.203, P 0.184 → **frozen
  outcome "does not replicate"**. Lenstool fails its control (12 / 18, 67 %); its 3 vs 0.499 (P 0.014) decides
  nothing under the rule. The GLAFIC anti (581, m 21.2, z_phot 0.70, offset 63.6°) is a bright edge-on-looking disc
  with the same axis in F200W and F444W (sheets looked at); in Lenstool it sits at g_min 1.19. Not a candidate.
- **Handoff:** S1063 tally: MACS0417 no, MACS0647 no result, MACS1423 explained, RXJ2129 no result, WHL0137 no,
  RXC J0600 no. Further fields with public maps and JWST (map survey 2026-10-18): SPT0615, MACS0308 (RELICS), Abell
  383, MACS1206, MACS1931 (CLASH), many VENUS (6882) clusters; none has a DJA photo-z catalogue yet.

## 2026-10-18: S2 with LS DR10 z < 22 galaxies: lensing control not detected (1.1σ); S2 paused (D-070 addendum 2)
- Frozen plan `results/s2_flat_kernel/dr10/PLAN.md`; `scripts/s2_flat_kernel.py --survey dr10 --mag-z-max 22` (new
  `--fit-mag-z-max`, photo-z gate `photoz_gate`). Same DES-SN5YR SNe as the DR9 run, so exploratory.
- 116,190 DR10 galaxies (2.5× DR9 z < 21). Gate Q passed (21–22 mag: σ_NMAD 0.035, outliers 3.1 %, 2,021 z_spec).
- **Gate C failed:** γ_L = −0.0145 ± 0.0129 (1.13σ, scramble p 0.25) vs DR9 −0.0183 ± 0.0110 (1.67σ); post hoc
  DR10 z < 21 gives −0.0098 ± 0.0115. The deeper column raises the error; λ_L 0.60 → 0.65.
- Flat term α = 2: γ_F = −0.0011 ± 0.0031; computed bound 0.0076 not adopted (control failed); the quoted S2 limit
  stays γ < 0.019. α = 1 kernel-degenerate (corr 0.991).
- **Failed approach:** a deeper LS photo-z column as the S2 positive control.
- **Handoff:** S2 paused; reopen only with DES Y3 Gold or a spec-z / deeper column and the full covariance.

## 2026-10-18: Upkeep: `docs/state/pending.md` synced with the open needs-human PRs
- Added rows for #199 (Abell 370 1565, parked), #237 (o052 / SMACS-VENUS grids; SN Atalanta) and #279 (five-field
  clump-pair census); #273 moved to a "merged since" note. Other rows' heads checked: unchanged.
- **Queue check (handoff):** no agent PR awaits merge; 18 await the owner. The recent lines are paused on owner
  calls: W3 ZTF (v4/v5), W1 metacal R (PSF-limited; reopen with a validated drizzled PSF or a star-rich field),
  MOA gb13 (#145), S1063 spectroscopy (116.28VK.001). Remaining TASKS items are optional or low-value (W1 1 × 10¹²
  M☉ grid ≈ 30–45× T&A; G-dependent external spike cap). No new search was started.

## 2026-10-19: W3 ZTF v9: official psfcat catalogues fill the field 672 gap for 3 calibration stars (overlap test passed)
- psfcat (IRSA IBE, per exposure and quadrant) exists and is readable for all 72 files requested. Its `mag` is
  instrumental, so MAGZP is added. Labelled deviation: the first run omitted it.
- Overlap with DR24: 13 pairs, pooled median Δ −0.015 mag. Per-star offsets −0.015, −0.015 and +0.020 (range −0.016
  to +0.030). Passes the frozen test. This shows consistency only: DR24 is built from the same PSF photometry.
- Gap fill (60 of 63 candidate exposures): strict usable gap observations rise from 1 per star to 28, 60 and 31.
  All 55 missing targets lie outside the footprint (52) or inside its 10″ margin (3).
- ZFPS not needed for these stars; difference images exist for 5 of the 6 checked.
- One cell's #350 cases would need ≤ 533 files (≤ 114 MB, about 13 min).
- **Handoff:** paused until the owner authorizes a psfcat data path in the W3 chain and a re-run of the calibration
  (`results/w3_ztf/v9/RESULT.md`).

## 2026-10-18: W3 ZTF v8: field 672 gap exposures absent from DR24 in both HATS and API; image fill fails its control; paused
- Fixed sample of 13 r-band objectids (Npix 7134). HATS and the IRSA light-curve API agree exactly on detections in
  the gap. In field 672, IRSA lists 110 covering exposures in the gap (64 with infobits 0) against 10 DR24
  detections. The release notes say DR24 does not include every field and exposure. No newer release is on IRSA's S3.
  The forced-photometry service needs an account, so it was not used (`results/w3_ztf/v8/`).
- Implementation fix: the exposure mid-time now uses `exptime` (299–300 s private exposures missed the match by
  134.5 s), and single-exposure matching is fixed. Regression test added. No stored result changes: every non-30-s
  exposure has infobits ≠ 0.
- Image trial (3 stars, 36 cutouts, 2.51 MB): the star is present in all 24 gap exposures (0.62–1.16 × its controls).
  The method fails its control: 0.54–0.75 × the DR24 flux at controls, no field stars found.
- Filling the 49 bracket cases of one cell would need about 530 star × exposures (about 30 min, about 37 MB).
- **Handoff:** paused; reopen conditions in `results/w3_ztf/v8/RESULT.md`.

## 2026-10-18: W3 ZTF v7: the v4→v6 recovery drop is mostly missing catalogue epochs (2022-10 to 2023-10), not a bug
- Post hoc stage trace (`results/w3_ztf/v7/`). It reproduces every stored count in both areas.
- Gap 756 of 2,000 injections: about 654–704 from pre-injection sampling eligibility (depending on the rate used),
  about 52–102 from pre-screen and vetting given eligibility (20 common-mode failures in the new area, 0 in the old).
- The 72 unexplained ρ 0.01, t_E 100 d cases all stop at MIN_BRACKET. 49 of them lack a post-event baseline that
  IRSA lists as official good exposures. Most new-area exposures from 2022-10 to 2023-10 have no DR24 catalogue
  detection; an IRSA image is listed for 1 of 2 distinct exposures checked.
- No implementation error changes the outcome; 2,600 unmatched 2023 detections are left open.
- **Handoff:** route paused. Reopen only with a later release covering those epochs, or owner-authorized image
  forced photometry.

## 2026-10-18: W3 ZTF v6: cross-sky validation failed (10 of 10 cells below 0.5); search group not read
- New area Npix 7134 (RA 184, Dec +33), chosen by a metadata rule; no ZTF field shared with v1–v5. Sets were fixed
  before reading: 533 calibration groups, 539 search groups, a control group (5,212 with light curves) (`results/w3_ztf/v6/PLAN.md`).
- The 10 calibrated cells recover 48–90/200 (0.24–0.45; 0.60–0.87 in v4). Controls pass vetting in ≤ 8/200. Real
  calibration light curves give 0 survivors. **Failed**, so the round stopped. The search group was never loaded
  (`source_usage.jsonl`).
- Post hoc (only ρ 0.01, blend 1 diagnosed): at t_E 30 d, 82 of 112 unflagged injections have fewer than 3 umbra
  nights (median 54.8 nights per year here). At t_E 100 d, 72 unflagged injections with ≥ 3 umbra nights are
  unexplained. The v4 cadence was not measured on the same basis.
  Reads: 386.7 MB.
- **Handoff:** the v4 capability did not reproduce in Npix 7134; route paused (`results/w3_ztf/v6/RESULT.md`).

## 2026-10-18: W3 ZTF v5: exploratory holdout search within the calibrated cells, 0 flags in 486 groups
- #339's failed verdict stands. This is a narrower, exploratory search with the #339 code unchanged; the plan was
  frozen first (`results/w3_ztf/v5/PLAN.md`). Capability is stated for the 10 injection cells at ≥ 0.5 only.
- Searched 486 holdout groups: 971 series, 1,438,841 series epochs. Result: 0 flagged runs, so no image check (0 MB).
  The holdout was not strictly blind: its light curves were downloaded in v1–v3 and used in aggregate comparison
  counts (`independence.json`).
- Post hoc functional check: the development half also gives 0 flags (1,520 under the v3 mapping).
- Reads: 1,024.8 MB compressed (`cost.json`).
- **Handoff:** no verified dimming; the route stays paused. No next batch without the owner (`results/w3_ztf/v5/RESULT.md`).

## 2026-10-18: W1 metacal R v3: empirical-PSF gate failed (2 usable Gaia stars); post hoc ePSF lowers R to 0.71 ± 0.17
- Frozen plan `results/w1_r_realmetacal/PLAN_v3.md`: ePSF = median of Gaia DR3 star stamps (17 ≤ G ≤ 21, RUWE < 1.4,
  5/6-parameter) from the same GOODS-N F200W mosaic; gate E needs ≥ 5 stars and FWHM 0.060–0.090″.
- **Gate E failed** (`results/w1_r_realmetacal/v3/RESULT.md`): 6 Gaia sources in the footprint, 2 usable (FWHM
  0.0782″ passes); v3 gives no R. Stamps looked at: clean stars with six-armed spikes, which the analytic PSF lacks.
- Post hoc, two-star ePSF: the simulated gate fails V1 (0.969 vs 1.634, tolerance 0.647); real tiles R = 0.707 ±
  0.169 vs analytic 1.006 ± 0.134 in the same run (19 tiles). The borrowed R = 0.414 stays; limits unchanged.
- **Failed approach:** Gaia-star ePSF in a high-latitude deep field (too few stars). A unit test caught a flat-core
  bug (plateau corner); fixed, outcome unchanged.
- **Handoff:** the metacal R route is PSF-model-limited at about ±0.3 and does not reach the 0.8 bar; pause it.
  Reopen only with a validated PSF (STPSF per exposure, drizzled) or a star-rich field.
## 2026-10-18: W3 ZTF v4: new-source validation failed (one W3 cell); footprint fix removes all flags; 0 dimming in 38 checked non-detections
- `scripts/w3_ztf.py` v4 (`--pix 3838`): non-detections only ≥ 10″ inside each exposure's quadrant footprint
  (IRSA corners, gnomonic after review), ambiguous groups dropped, one entry per expid; plan frozen before reading (`results/w3_ztf/v4/PLAN.md`).
- 1,000 new groups (Npix 3838; same fields and windows as v1–v3): 0 flagged runs (2,565 under the v3 mapping);
  controls pass vetting ≤ 9/200; W3 recovery 0.85/0.865/0.71/0.405 in the pass cells → **failed** (ρ 0.1,
  t_E 100 d; sharp-edge test, post hoc diagnostic). t_E 3 d 0.40–0.45.
- Images of the 5 flags the fix removed (post hoc, 44 cutouts, 3.10 MB): of 38 non-detections, 21 lie in a
  quadrant without the objectid (11 photometrically present), 17 are not covered or at an edge, 0 show dimming.
  Archive reads ≤ 1.42 GB (file sizes, `cost.json`).
- **Handoff:** paused; no expansion without the owner (`results/w3_ztf/v4/RESULT.md`).
## 2026-10-18: W1 real-image metacal R v2: 1.056 ± 0.156 on 18 GOODS-N tiles, still "conservative" (D-105 addendum)
- Frozen plan `results/w1_r_realmetacal/PLAN_v2.md`: noise-filled mask with a 16-px source pad, PSF σ floor
  (0.97 px), unique matching, a tile gate (raw PSF σ below the floor or < ½ median detections), all tiles, tile
  bootstrap. Simulated gate re-run: passes (unchanged from v1).
- **Result (`results/w1_r_realmetacal/v2/RESULT.md`):** 22 tiles, 4 gated out; R = 1.056 ± 0.156 (804 sources),
  R − 2σ = 0.744 < 0.8 → "conservative"; borrowed R = 0.414 and the blank-field limits unchanged. No fixnoise 1.075 ±
  0.052 (not the estimator); PSF_d × 1.15 0.875 ± 0.162. R_11 < R_22 in every variant and in the simulation (a
  chain property, not investigated). Review: the frozen 0.97 px floor is PSF_d's σ, but σ is measured on PSF_r
  (1.37 px) images, so it acts only via the gate; post hoc with 1.37 px, R = 1.114 ± 0.168, R − 2σ 0.778 (same call).
- **Failed / gap:** the noise fill does not stop the ringing; the gated tiles all ring after deconvolution
  (`v2/sheet.png`, looked at), one without any masked pixel, so a sharp data feature drives it.
- **Handoff:** the PSF model is now the limiting term: an empirical drizzled PSF (stars of the mosaic or STPSF on
  the drizzle) as PSF_d, then a second field (CEERS F200W) if the bar is still missed.

## 2026-10-18: W1 R by metacalibration of real GOODS-N F200W images: 0.95 ± 0.24, borrowed R = 0.414 conservative (D-105)
- `scripts/w1_r_realmetacal.py` (plan frozen first, `results/w1_r_realmetacal/PLAN.md`): GalSim metacalibration
  (deconvolve an analytic PSF, shear ±0.01, reconvolve with a 0.10″ Gaussian, white-noise fixnoise) of six 1024²
  tiles of the JADES 1181 o001 F200W i2d streamed from MAST S3; shapes from the `w1_r_sim.measure` emulation.
- **Gate (simulated):** V1 metacal 1.442 ± 0.226 vs direct shearing 1.732 ± 0.064 (tolerance 0.642) and V2 sign
  both pass; weak gate. Post hoc four scenes: metacal reads lower than direct (1.519 vs 1.844), so it tends to
  understate R.
- **Result (`results/w1_r_realmetacal/RESULT.md`):** R = 0.953 ± 0.239 (323 sources), frozen rule "conservative"
  (R − 2σ 0.475; a post hoc bootstrap over tiles gives ±0.184); no-fixnoise 0.979 ± 0.067 (anisotropic: R_11 0.866, R_22 1.093); PSF × 1.15 0.886 ± 0.237. The
  GOODS-N / CEERS / PRIMER limits stand at R = 0.414; the error is too large to adopt R ≈ 1.
- **Failed / gap:** one tile rings after deconvolution (a sharp feature; by eye on `sheet.png`), 94 detections vs
  1,755–2,370, and its PSF σ estimate is 1/√12 px (a one-pixel artefact); post hoc without it R = 1.004 ± 0.270.
  The tool needs a bad-pixel mask, apodised masked pixels, a PSF σ floor and unique matching.
- **Handoff:** to tighten R: empirical drizzled PSF (stars of the mosaic), a mask before deconvolution, more tiles and
  a second field; then a frozen plan to re-run the blank-field limits with the measured R if its 2σ bound clears ~0.8.
## 2026-10-18: W3 ZTF v3: counting fixed; for three survivors the non-detections are coverage gaps and quadrant changes, not dimming
- `scripts/w3_ztf.py` v3: official exposure metadata (IRSA TAP), each star counted once per exposure, target
  excluded, in-play denominator; regression tests. Development half only (not blind); holdout unread.
- Survivors 9 → 8 (same two windows); W3 recovery at ρ = 0.1, t_E = 100 d 0.24 → 0.30; frozen conditions still fail.
- Images (IBE cutouts, 23 exposures of 3 of the 8 survivors, 1.34 MB): of 11 catalogue non-detections, 6 have no
  image, an invalid exposure or an edge, and 5 show the star present in a neighbouring quadrant that has no objectid
  for it; none is a real drop. The other 5 survivors were not imaged (`results/w3_ztf/v3/RESULT.md`).
- **Handoff:** paused. A new independent validation would add footprint containment from metadata corners,
  cross-quadrant merging and image checks of survivors (costs in RESULT); not started.

## 2026-10-18: W3 feasibility trial on ZTF DR24 light curves: failed its frozen validation (v1 and bug-fixed v2); no holdout search
- Plan, sample (1,000 isolated stars, metadata only) and rules committed before any light curve was read
  (`results/w3_ztf/PLAN.md`); data: IRSA ZTF DR24 HATS on S3 (one 953-MB light-curve partition read in 61 s); model:
  W3, negative-mass point lens (n = 1, ε < 0; finite-source spikes), a restricted hypothesis, not wormholes in general.
- **Result (`results/w3_ztf/RESULT.md`):** v1 failed (W3 recovery 0.215 at ρ = 0.1, t_E = 100 d; 6 real survivors).
  A review found code bugs (time bins splitting exposures, night boundary, control models, missing injection noise,
  seeds); v2 fixed them with unchanged thresholds and failed again on the same development half: W3 recovery
  0.56–0.71 for t_E 10–100 d at ρ = 0.01 and 0.565–0.62 for t_E 10–30 d at ρ = 0.1 (blend 1), but 0.24 at ρ = 0.1,
  t_E = 100 d; controls pass vetting at 0–1.5 %; 9 real survivors, all in two sample-wide non-detection windows.
  The holdout half was not searched.
- **Why:** detections-only data (non-detections only by inference), a common-mode test that misses sample-wide
  episodes (its dim fraction can go negative: a known, unfixed v2 bug), and slow large-ρ t_E edges that look like
  ordinary fades. **Handoff:** paused.
## 2026-10-18: W1 lower masses: filter pilot keeps the frozen Schirmer 10″ screen (D-103)
- Frozen plan `results/w1_low_mass/PLAN.md`: three filters (V0 Schirmer 10″, V1 point-mass 10″, V2 Schirmer 5″) ×
  masses 2 × 10¹², 1 × 10¹², 5 × 10¹¹ M☉ × 100 injected lenses at Z = 0.3 on ceers-o001 and goodsn-o004 (no streak masks).
- **Result (`results/w1_low_mass/RESULT.md`):** recoveries at 1e12 + 2e12: V0 15 / 21, V1 5 / 9, V2 14 / 12; no
  variant reaches the frozen ≥ 1.5× V0 in both fields, so V0 stays. V0 at 1 × 10¹² M☉: 4 and 3 / 100; 5 × 10¹¹: ≤ 1.
  V0's real-data S_max reproduces the committed Z = 0.3 summaries in both pointings. V1 / V2 real-data screens (new
  looks): no E/B pass.
- **Handoff:** the plan's next step is a V0 grid at 1 × 10¹² M☉ over the forty fields; expected ≈ 30–45× Takahashi &
  Asada (rough), so low value per cost (T4). Other W1 options: an R measurement on real blank-field images; a
  G-dependent external spike cap (removes false flags only).

## 2026-10-18: W1 PRIMER UDS x3y0 explained (post hoc) as a G 10.07 star's optical artefact in one visit
- `scripts/w1_cal_trail.py` (plan frozen first: `results/w1_shear_primer/vet/x3y0_cal_plan.md`) resamples every
  `_cal` exposure covering a flag onto the mosaic grid and tests the line through it. Four F200W exposures of
  observation 001 cover the x3y0 peak (visit 001 and visit 011, two dithers each).
- **Frozen outcome: undecided.** 2 of 4 exposures are above 5σ, from two different pointings, which matches no
  branch of the plan. The per-exposure calls also disagree with the sheets: the line is clear in both visit-011
  dithers (one scores 4.50 because stray light enters its control band) and absent in visit 001. The band-median contrast is a
  weak per-exposure statistic.
- **Post hoc:** the line is radial from Gaia DR3 2489648647202619264 (G 10.07, 57.6″ from the peak) and passes
  0.17–0.19″ from it. Both visit-011 dithers put it on the same sky row; a detector-fixed line would move about 17
  rows. Visit 001 shows nothing at the same positions (chunk peaks 0.009–0.028 vs 0.40–0.89 MJy/sr). x3y0 is
  explained as bright-star contamination, not a candidate. The forty-field limit is unchanged (flagged tiles were
  already excluded).
- **Method gap found:** the screen's Gaia spike veto caps external stars at 60″ (`EXTERNAL_SPIKE_CAP`). This star's
  uncapped law gives 291″, and its arm is traced to about 71″ (not stored). Bright (G ≲ 11) stars near tile edges can therefore
  still make E/B flags.
- **Handoff:** (1) W1 at lower masses (a new frozen plan); (2) an R measurement on real blank-field images;
  (3) optional: a G-dependent external spike cap tested by injection on the PRIMER tiles (it would only remove
  false flags, not change the limit's tiles).

## 2026-10-18: MACS1423 h000 HST archive light curve: no change before (−181 d) or after (+180 to +409 d); identity still open
- `scripts/hst_forced.py` (HST counterpart of `hit_lightcurve.py`: S3 byte-range drizzled products, 0.2″ aperture,
  local random-aperture errors, field peaks; `--offsets`) on 43 HST visits 2006–2024; visits registered to the
  JWST 1208 F150W frame (labelled deviation: the host gradient made per-visit WCS offsets of up to 0.27″ look like
  variability). Plan, result and light curve: `results/transients/macs1423_h000_hst/`.
- **Result:** four new epochs tested against same-instrument baselines (2022-08 before; 2023-08, 2023-12, 2024-04
  after). No change at the position: −181 d vs −380 d F435W/F606W 3σ < 19–22 nJy (and vs 2013 < 41–63); after, F435W
  < 28 nJy (+180 d), UVIS F606W < 41 (+290 d, vs 2014), F110W < 68 / 97 nJy (+290 / +409 d, vs 2013). Consistent with
  one episode of < ~470 observed days; a source that varied the year before or stayed bright after is weakened;
  supernova types, other one-off transients and magnification of a faint persistent source are not separated (no
  colours over time). No low-cost archive step remains for this question. Review fixes: field zero points matched
  in the registered frame and de-duplicated; single-visit IR errors scaled ×1.5; three baselines added (labelled).
- #332 wording narrowed: the achromaticity test is a band-consistency test under stated assumptions, not a general
  supernova-versus-lensing classifier.

## 2026-10-18: W1 shear in PRIMER: forty-field limit 0.638× Takahashi & Asada at 2 × 10¹³ M☉ (R = 0.414 ASSUMPTION); UDS x3y0 flag unresolved
- Staged: the four DJA PRIMER zouts (tarball names end in `-fix`; the planned names return 403), pinned and stacked
  per survey (`w1_stage_inputs.stack_zouts`, content sha256 pinned). Stacking every tile's survey pair instead of
  only tiles that straddle two DJA mosaics is an implementation choice (deviation from the plan's wording; the
  0.3″ nearest match differs only where DJA mosaics overlap).
- Streak masks: the SCI extensions are 1.9–3.3 GB, not the 31–69 GB file sizes, so `w1_streak_mask.py` streams each
  mosaic once (concurrent range reads, ≈ 40 s each), normalises it whole and cuts each tile (`tile_slices`).
  Contact sheets (`masks/*.png`): wide mask bands on several tiles; visible diagonal satellite trails in UDS x0y2 /
  x3y2 are partly unmasked (the frozen rule, not tuned).
- **Gate P (`gate_p.json`):** all 20 tiles pass (photo-z fraction 0.592–0.873).
- **Screens (`screen_flags.json`, Z = 0.3 summaries):** four tiles pass the E/B rule: COSMOS x0y0 (S 4.44), UDS x0y2
  (4.17), x2y0 (4.85), x3y0 (6.47); all four peaks lie within 20″ of a tile edge (7 of 20 maxima do). Frozen V1 ∪ V2
  vetting (`vet/`): x0y0, x0y2 and x2y0 explained (spike / streak). **UDS x3y0 is not explained by the frozen rule**;
  on `sheet_smax_flags.png` the peak sits on a bright horizontal line lined with catalogue detections (by eye a
  streak; not a frozen test). Removing the row band (V1r, not part of the frozen outcome) drops S to 2.99 (p 0.99);
  5 % of rotated control bands do as much (`vet/vet_primer-uds-x3y0.json`). Not a candidate; unvetted.
- **Run 3 (2026-10-18 08:10–08:52Z, takeover of the released claim):** inputs re-staged and streak masks regenerated;
  every regenerated mask JSON equals the committed one (wall time aside), and every new summary's `screen` block and
  mask `data_sha256` equal the committed Z = 0.3 ones (`assert` before copying), so the restaged chain reproduces run 2.
  Injections at Z ∈ {0.6, 1.0, 1.5} on the sixteen unflagged tiles (flagged tiles not injected, D-098, as CEERS):
  25 of 48 done (the `summaries/z{0.6,1.0,1.5}_*.json` files). ≈ 230 s per tile-Z job with four in parallel (not 6 min per Z); the streak-mask tool
  needs the `streaks` extra (`acstools`). UDS x3y0 (by eye on the sheet): one straight horizontal line crosses the whole
  40″ sheet through the peak, with fainter parallel streaks above it; the tile's mask lists near-horizontal status-2
  trails (θ 0.85°, 178.4°, 178.5°), so the trail family is found, but the line through the peak is not fully vetoed
  (vetoed fraction not measured). Still a flag, not a candidate.
- **Run 4 (2026-10-18 09:12–09:50Z, takeover):** the 23 remaining UDS tile-Z injection jobs, then R = 0.231 at Z = 0.3
  on the sixteen; regenerated masks and every new `screen` block equal the committed ones (checked in the session).
  **Limit (`results/w1_shear_primer/RESULT.md`, `limits_40fields_unflagged.json`):** forty fields, 446.2 arcmin²,
  n₉₅ < 6.38 × 10⁻⁵ h³ Mpc⁻³ at 2 × 10¹³ M☉ (0.638× Takahashi & Asada; 1.50× / 9.32× at 8 × 10¹² / 2 × 10¹²);
  V_eff +71 %. The plan's trigger holds: next is lower masses, not more area. Screened area is 0.937 of the
  covered tile area (`area_check.json`). R = 0.231 keeps 0.15–0.28 of the PRIMER recoveries (`r_sensitivity.json`,
  sixteen-tile totals), so the sub-T&A reading needs R ≳ 0.414 (D-101 supports it for blank fields).
- UDS x3y0 visit check (post hoc, `vet/x3y0_visit_check.json`): at the peak c1002 is effectively visit o001
  (max |diff| 0.053 MJy/sr), so the visit test cannot decide; by eye the line ends ≈ 15″ east of the peak.
- **Handoff:** (1) x3y0: the o001 L2 `_cal` exposures at the peak (a line in one exposure is a trail); (2) W1 at lower
  masses (a new frozen plan: smaller Schirmer aperture or stacking, injections at ≤ 2 × 10¹² M☉); (3) a PRIMER or
  real-image R measurement, since the result depends on R.
## 2026-10-18: Confirmation round 2: SMACS h005 within control scatter; an achromaticity test cannot decide G165 h002
- **SMACS 1 h005 (`results/transients/confirm_r2/RESULT.md`):** the control-photometry stall was one MAST portal
  query per position hanging on an HTTP response; `hit_lightcurve.py --group-arcsec` shares one query. With 22
  controls (entry and duplicates removed) the decisive F150W tail is 0.591: consistent with matched-control scatter.
- **Discriminating test (D-102):** band consistency (achromaticity, under stated assumptions; not a general classifier) of B1-PLCKG165 h002's fade (intrinsic transient vs magnification by a
  point mass, Ellis wormhole or negative-mass lens). Frozen rule: not decided (16″ pick: all bands < 20 usable
  controls; 30″ deviation: F090W 14 < 20). Simulated power against SN templates 0.15–0.62 (≤ 0.29 for type II), so
  the two epochs are weak for this question; paused until a calibrated bright-epoch repeat or a spectrum exists.
- **Prior work (`docs/value_comparison.md`):** JADES (2406.05060) and the Dragon-arc papers already search brightening
  and fading; Stein et al. find look-alikes and lens systems, not the same source (2110.00023 is the lens-mining
  paper, 2110.13151 the similarity search). Our measured added value is per-box calibration; coverage beyond other
  teams' searches is not shown.
- Wording correction to #327: the 26,145 SED-matched pairs minus the 10 published ones are unverified, not disproved.
- **Handoff:** no new search; nothing confirmed or exotic.

## 2026-10-18: W1 shear in PRIMER: plan frozen, mosaics tiled, catalogues pinned (no screen yet)
- `results/w1_shear_primer/PLAN.md` (frozen before any screen). MAST holds two UDS combined mosaics (c1001, c1002;
  129.2 / 127.9 arcmin²) and two COSMOS mosaics (o003, o004; 92.4 / 87.4 arcmin², overlap 38.7) in F200W
  (`footprints_*.json`); the single-visit products (837 s) lie inside them and are not used.
- The mosaics exceed the screen's 200″ radius, so they are cut into non-overlapping 280″ tiles
  (`w1_footprints.tile_grid`, new; rows kept by `inject_radial.in_box` through the spec key `box_half`): 12 UDS and
  8 COSMOS tiles, share-weighted 204.8 + 129.8 arcmin², about 1.7× the current headline area if all pass gate P.
  Tile sky covered only by the other mosaic is dropped (conservative; UDS keeps ≈ 80 % of its union, not stored).
- Four F200W catalogues pinned (`w1_stage_inputs.py`). The `_i2d.fits` are 31–69 GB, so the whole-SCI streak-mask
  stream is not used; a tile-range read is needed.
- **Handoff:** next run: stage the four DJA PRIMER zouts (streamed out of 345–536 MB tarballs; a tile may need two),
  register the 20 tiles in `inject_shear.py`, tile-range streak masks, then gate P, screens, injections and the
  combined limit as the plan says.
## 2026-10-18: W1 blank-field responsivity R from image simulations: borrowed R = 0.414 is conservative
- `scripts/w1_r_sim.py`: GalSim Sérsic galaxies sheared by ±0.02, jwst-like 2 px kernel, 3σ segmentation and
  photutils `SourceCatalog` moments on the smoothed data, then `intrinsic_ellipticity`; plan frozen first
  (`results/w1_r_metacal/PLAN.md`). Checks: noiseless R 0.998; sign and round-source checks pass.
- **Result (`results/w1_r_metacal/RESULT.md`, simulated):** weighted by 2,259 GOODS-N o001 rows, R_sim =
  1.017 ± 0.006 (n = 1) / 0.965 ± 0.007 (n = 4) with the field's PSF σ (post hoc), 1.34–1.37 with the plan's stamp
  PSF σ, 0.940–2.057 over every post hoc variant (resolved-cell weights, jwst-style threshold). All pass the frozen
  rule (≥ 0.414 − 2σ): the GOODS-N / CEERS limits are conservative with respect to R; no limit changed.
- The cluster-fit R (0.41–0.60) is about half the simulated measurement response; the cause (lens model, photo-z
  default, member dilution) is untested.
- D-101. **Handoff:** a real-image metacalibration on F200W i2d cutouts (S3 byte ranges, STPSF), and the simulation on a
  cluster field's rows to locate the factor-two gap. If the real-image R stays ≈ 1, the blank-field limits can be
  re-run with it (a new frozen plan).
## 2026-10-18: Same-source pairing feasibility fails its frozen rule; VENUS batch 4: 13 boxes calibrated, nothing confirmed
- **Inventory and tools (`docs/same_source.md`, D-100):** the repository aligns epochs, finds brightening/fading
  sources and predicts counter-images with published models; it has no model-free step that proposes same-source
  pairs. lenstronomy 1.14.2, PyAutoLens 2026.10.17.1 and Lenstruction (2022, pinned to lenstronomy 1.10.3) were tried
  in an isolated venv: all test or reconstruct a chosen pair/system; none proposes pairs. None adopted.
- **Feasibility test (`results/same_source/a2744_feasibility/`, plan and pass function committed first):** identity-blind
  SED ranking at 2–60″ on Abell 2744 (Bergamini+23 images for scoring only). Holdout recall@10 0.86, but true pairs
  pass the absolute SED rule only 17 % / 32 % (dev / holdout): **fails** (P2). Morphology (affine, either handedness)
  AUC 0.53 vs matched decoys. Post hoc: 26,145 SED-matched pairs field-wide, of which 10 match published same-system pairs (the rest unverified, not shown to be false).
- **Identity check (`results/transients/identity_r2/identity_check.json`):** PLCK G165 h002, MACS1423 h000, Bullet
  h019: no record found in the sources searched; TNS 403 (result unknown). Leads: a Legacy Survey galaxy 0.97″ from
  G165 h002 (possible host, unchecked); JWST GO 9811 (SN-ECHOES) lists a planned Cycle 5 G165 epoch (footprint not checked); HST 17504 (Kelly SNAP)
  imaged MACS1423 in 2023-12 – 2024-04 (not checked).
- **VENUS batch 4 (`results/transients/batch_venus4/RESULT.md`):** 20 boxes, 13 calibrated (2.08 arcmin²), 7
  undeterminable (2 timed out); 204 hits → 70 → 18 review entries, 10 reviewed, 8 pending (calibration-failed AS1063 1).
  Two rule-level real changes (A370 1 h000, SMACS 1 h005), checked against matched compact controls: A370 1 h000 consistent with matched-control scatter in the decisive F150W pair (tail 0.455, 33 controls; flux ratio 0.998); SMACS 1 h005 not decided (control photometry stalled twice); A370 2 h004 and A370 1 h006 undetermined; 8 AS1063 1 entries pending.
- **Handoff:** no next batch without the owner; held undeterminable boxes 33.
## 2026-10-18: W1 shear in CEERS (program 1345, no lens model): twenty-four-field volume limit 1.09× Takahashi & Asada
- Ten F200W pointings (o052 left out for overlap; `scripts/w1_footprints.py`), DJA v7.4 photo-z streamed out of the
  397 MB tarball, R fixed at 0.414 (D-099 addendum); plan frozen first (`results/w1_shear_ceers/PLAN.md`). All ten
  pass gate P (0.695–0.813).
- **Screens:** eight null; o004 (S 4.96) and o069 (S 4.07) pass the E/B rule. On the sheet o004 sits on a detector-gap
  edge lined with detections next to a bright star's spikes, o069 on a horizontal spike-like line; both explained by
  the frozen V1 ∪ V2 rule (S 1.54 / 1.85). Not candidates; out of the headline (D-098).
- **Result (`results/w1_shear_ceers/RESULT.md`):** twenty-four fields, 199.9 arcmin²: n₉₅ < 1.50 × 10⁻³ / 2.53 × 10⁻⁴ /
  1.09 × 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉, i.e. 15.0× / 2.53× / 1.09× weaker than T&A (from 24.6× /
  4.29× / 1.82×). V_eff at 2 × 10¹³ M☉ +68 %.
- **Caveat:** with R = 0.231 the CEERS recoveries fall to 0.19–0.28 of the adopted ones (eight-pointing totals, Z = 0.3).
- **Failed step (process):** deleting the mask FITS files mid-grid made the next runs stop with FileNotFoundError (no
  silent loss); the regenerated masks are byte-identical to the sha256 in the summaries.
- **Handoff:** the borrowed R is now the dominant systematic and the headline is within reach of T&A: calibrate R for
  blank fields (image simulations of the pipeline shapes) before adding PRIMER-UDS / COSMOS area.

## 2026-10-18: W1 shear in GOODS-N (JADES 1181, no lens model): sixteen-field volume limit 1.82× Takahashi & Asada
- Six F200W pointings with small mutual overlap, DJA v7.3 photo-z, R fixed at 0.414 (D-099); plan frozen first
  (`results/w1_shear_goodsn/PLAN.md`). Gate P (photo-z share ≥ 0.5) drops o009 / o010 (0.22 / 0.02).
- **Screens:** o001, o005 null; o003 (S 4.77) and o004 (S 4.35) pass the E/B rule. Both sit on mosaic edges on the
  sheet (a spike line of detections; an edge column of detections), and both are explained by the frozen V1 ∪ V2
  rule (S 3.80 / −0.18). Not candidates.
- **Result (`results/w1_shear_goodsn/RESULT.md`):** sixteen fields (D-098: flagged pointings out), 119.8 arcmin²:
  n₉₅ < 2.46 × 10⁻³ / 4.29 × 10⁻⁴ / 1.82 × 10⁻⁴ h³ Mpc⁻³ at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉. That is 24.6× / 4.29× /
  1.82× weaker than T&A, down from 30× / 5.4× / 2.35×. Eighteen fields (the plan's rule, not adopted) would give 1.54×.
- **Caveat:** with R = 0.231 the GOODS-N recoveries fall to 0.21–0.26 of the adopted ones (four-pointing totals;
  0.13–0.33 per pointing and mass), so the gain depends on the borrowed R.
- **Handoff:** more blank F200W area with photo-z (PRIMER-UDS / COSMOS, CEERS: MAST catalogues plus DJA zout), or
  calibrate R for blank fields first (the dominant systematic now). A photo-z catalogue covering o009 / o010 would add
  ≈ 20 arcmin².

## 2026-10-18: W1 row/column-segment stripe veto for the CANUCS flanking fields: failed on cost, headline unchanged
- `w1_streak_mask.py --seg-len` (`segment_stripes`: field-wide segment medians, 41-row detrend, 5σ; rows and
  columns) OR-ed onto the D-097 trail mask; `w1_row_veto_eval.py` decides K1r and the cost check. Plan frozen first
  (`results/w1_row_veto/PLAN.md`). `w1_stage_inputs.py` now stages the NCF photo-z and pins the F150W catalogues.
- **Result (`results/w1_row_veto/RESULT.md`):** K1r passes for both variants, but the screens keep only 0.68–0.75
  (L = 80) and 0.56–0.67 (L = 40) of the rows the D-097 screens used (rule ≥ 0.75 in all five): not adopted, no
  injections, the W1 headline stays at fourteen fields (2.35× T&A).
- Post hoc: with either union mask macs0417-ncf no longer passes (S 3.15 / 3.21); abell370-ncf passes under L = 80
  only (S 4.62, p 0.0), so S_max moves a lot when a third of the rows go. Sheets: the mask also hits bright-source
  halos, spikes and detector edges.
- **Failed approach:** an unmasked segment-median stripe detector. **Handoff:** mask source footprints before the
  segment medians and require multi-segment stripes (new frozen plan), or move to other blank fields (PRIMER, CEERS).

## 2026-10-18: W1 faint-streak veto (lower findsat_mrt threshold) for the CANUCS flanking fields: failed, headline unchanged
- `w1_streak_mask.py --mrt-threshold` (default 5 unchanged); plan frozen first (`results/w1_streak_rowdiag/PLAN.md`).
- **Result (`results/w1_streak_rowdiag/RESULT.md`):** K1n failed for both pre-declared variants (threshold 3, 2.5):
  more candidates, but the persistence check keeps the same two status-2 trails, masks identical pixel for pixel;
  neither flagged NCF peak is covered. No injection run; the W1 headline stays at fourteen fields (2.35× T&A).
- Sheets looked at: the macs0416-ncf "diagonal streak" lines up with a bright star's diffraction spikes (the V2 case);
  macs0417-ncf sits on detector-wide horizontal row striping. Neither is a single long trail.
- **Failed approach:** lowering the TrailFinder threshold (D-097 "Revisit if"). **Handoff:** a field-wide row-segment
  median veto (V1r made field-wide), frozen and run through the injections, and a check of why the Gaia spike veto
  misses the macs0416-ncf star; or skip to other blank fields (PRIMER, CEERS) under the D-097 mask.
## 2026-10-18: Confirmation round 1: A370 box-2 changes not source-specific / not supported; MACS0417 h002 is SN Athena (#323)
- **A370 h003:** consistent with 54 compact-source controls, which share its F444W drop between 1208 and VENUS
  (median ratio 0.852 vs 0.887): no evidence of a source-specific change. **A370 h010:** within 80 matched arc-pixel controls (tail 0.10): no support
  for a real change. Batch 3 has no change supported beyond matched-control scatter.
- **MACS0417 h002 = SN Athena**, found in October 2025 by another team (JWST DD 12781, GO 4568, HST 18470 follow-up).
- **HST 18470 visits 1–4 (public, to 2026-10-07):** no new point source at the 12781 or our predicted counter-image
  position (+1.4σ at both).
- **Decision:** not worth requesting new observations now; the predicted reappearance is covered by approved
  programmes (`results/transients/confirm_r1/DECISION.md`).

## 2026-10-18: W1 shear in the CANUCS flanking fields (no lens model): fourteen-field volume limit 2.35× from Takahashi & Asada
- `inject_radial.NullLensModel` (κ = g = 0) and five NCF entries in `inject_shear.py` (F150W: the NCF have no
  F200W; CANUCS NCF photo-z; fixed R). Plan frozen first (`results/w1_shear_ncf/PLAN.md`).
- **Gate R passed** (`gate_r.json`): F150W R 0.425 ± 0.035 / 0.618 ± 0.037 on the MACS0417 / MACS1423 cluster
  pointings vs F200W 0.488 / 0.604; adopted 0.425 (D-098; a first run without the Gaia spike veto gave 0.421, re-run).
- **Screens:** three NCF null; macs0416-ncf (S 4.64, p 0.005) and macs0417-ncf (S 6.02, p < 0.005) pass the E/B
  rule. Sheet first: a diagonal streak the trail mask missed and horizontal stripes cross the two apertures.
  macs0416-ncf explained under the frozen V1/V2 rule (S 3.53, p 0.575); macs0417-ncf only post hoc (V1 on rows:
  S 1.35; the same band rotated about the peak reaches that in 9 % of 200 angles: weak). Not candidates.
- **Result (`results/w1_shear_ncf/RESULT.md`):** fourteen fields (both flagged NCF out), 97.5 arcmin²: n₉₅ <
  3.03 × 10⁻³ / 5.35 × 10⁻⁴ / 2.35 × 10⁻⁴ h³ Mpc⁻³, 30× / 5.4× / 2.35× weaker than Takahashi & Asada (from 65× /
  8.5× / 3.6×). Sixteen fields would give 1.9× at 2 × 10¹³ M☉ (rests on the post hoc explanation).
- **Failed approach:** the frozen `w1_peak_vet.py` V1 looks for stripe columns only; it cannot see horizontal or
  diagonal streaks. **Handoff:** freeze a row/diagonal streak veto and run it through the injections, then more blank
  fields (TASKS).

## 2026-10-17: W1 shear null as a volume density over lens redshift: 3.6× weaker than Takahashi & Asada at 2 × 10¹³ M☉
- `inject_shear.py --lens-z` (W1 lens off the cluster redshift; default path unchanged) and `scripts/w1_volume.py`
  (ε integrated over comoving volume); plan frozen first (`results/w1_volume/PLAN.md`, addendum A before the grid).
- Gate Z (Abell 2744 at its cluster z reproduces 0 / 0 / 14 / 72 / 146) passed on the second attempt; inputs equal
  the stored streak-veto run on all eleven fields (`inputs_check.json`).
- **Result (`results/w1_volume/RESULT.md`):** eleven fields, Z 0.1–1.5: n₉₅ < 3.65 × 10⁻⁴ / 8.5 × 10⁻⁴ / 6.5 × 10⁻³
  h³ Mpc⁻³ at 2 × 10¹³ / 8 × 10¹² / 2 × 10¹² M☉, i.e. 3.6× / 8.5× / 65× weaker than Takahashi & Asada (2013).
- **Failed approach:** a three-mass gate run cannot reproduce a five-mass run (seeds spawned per `--mass` entry).
- **Deviations:** RXJ2129 out (its MAST F200W catalogue `jw02767-o002_t001` now 404 on S3 and the MAST API); 100
  lenses per cell. **Handoff:** the W1 shear line does not reach the literature bound; reaching it at 2 × 10¹³ M☉
  needs ~3.6× the area at similar efficiency. Single-plane assumption untested (a multi-plane check is the next
  step only if the line continues). Restore RXJ2129 when MAST serves its catalogue again.

## 2026-10-17: Radial screen on MACS0717 with ASTRODEEP-FF photo-z: no flag
- `scripts/radial_more.py --fields macs0717` (unchanged code), plan frozen first (`results/radial_macs0717/PLAN.md`).
- **Result (`results/radial_macs0717/RESULT.md`):** 411 elongated, 77 dropped as not behind the lens, 37 `anti` not
  model-radial, 4 centres, max 5 lines, p_random_max 0.69: no flag (2026-10-08 run without photo-z: 68 / max 5 /
  p >= 0.70). Top-centre sheet looked at: faint segments between two bright members, no dark centre.
- #319 follow-ups also landed there: zout pinned by a column-value hash and verified before replacing the file.
- **Handoff:** the TASKS item is closed. Photo-z cover ~17 % of MACS0717 rows, so this run is only partly cleaner;
  next `radial` work stays the collinearity streak rule (TASKS JWST focus).

## 2026-10-17: W1 shear on MACS0717 with ASTRODEEP-FF photo-z: D-053 rule not passed, nearly blind; headline unchanged
- **Why:** MACS0717 was the one screened HFF cluster without W1 (no DJA photo-z). HFF-DeepSpace downloads return
  HTTP 401; ASTRODEEP-FF (Di Criscienzo+17, VizieR J/A+A/607/A30) is public: `scripts/astrodeep_zout.py` → zout.
- Plan frozen first (`results/w1_shear_macs0717/PLAN.md`); same code, seeds and streak veto as D-097.
- **Result (`results/w1_shear_macs0717/RESULT.md`):** S_max 3.14, p_rot 0.28, below max |S_×| 3.28: D-053 rule not
  passed. R = 0.23; recoveries 0 / 0 / 0 / 7 / 9 of 200. 83 % of the shape rows lack photo-z (not stored), so MACS0717 stays out
  of the headline (twelve fields, D-097); counted anyway, limits move ≤ 1.5 %. S_max sheet looked at.
- **Failed approach:** HFF-DeepSpace v3.9 tarballs (cosmos.phy.tufts.edu) need credentials (401, 2026-10-17).
- **Handoff:** W1 field additions are at diminishing returns (T4); the low-R fields (MACS0717, MACS0647) would need
  a member-aware responsivity calibration (hypothesis: members without photo-z dilute R; untested). Open
  bookkeeping (none changes the null): `inject_radial.combine` admits any field with a photo-z file, so MACS0717's
  exclusion is prose only; a photo-z coverage threshold there is a frozen-rule change and needs its own plan. The
  62 spectroscopic rows with `relFlag` 0 are kept (possible blends; keep or drop is open). The zout pin is a
  column-value hash (`astrodeep_zout.content_sha256`), so an astropy upgrade does not reject identical content.

## 2026-10-17: W5 counts at the `radial` centres: no deficit; the deficit test has no power at JWST depth
- `scripts/w5_radial_counts.py`; plan frozen first (`results/w5_radial_counts/PLAN.md`, `RESULT.md`). 30 centres,
  four fields, radial arcs removed: Q = N_obs / control 1.91 / 1.76 / 1.70 at 3″ / 6″ / 9″, deficit p ≥ 0.9995 (null).
- **Failed approach:** at m_lim 27.3–28.2 the counts are shallow, so the negative-mass model predicts an *excess*
  (Q_H 1.6–2.9 at 3″ for 1 × 10¹² M☉), not a deficit; power 0–0.0035. The recipe's "~10 % missing" holds for DR10 depth only.
- Post hoc: without the rows within 2″ of the arcs, Q 1.34 / 1.23 / 1.26 (p_high 0.104 / 0.058 / 0.0075), nearly
  flat in R, unlike H's profile; read as a local overdensity at the centres. Not a candidate. **Handoff:** W5 at the
  `radial` centres closed; an excess test needs a density-matched control and a new plan (low priority).

## 2026-10-17: W1 stack injections, mass points 5 × 10¹¹ and 1 × 10¹² M☉: f = 1 limit 1 × 10¹² M☉
- Same code, sample, seed and noise-aware rule (`results/w1_stack_radial/inject/PLAN_mass.md`,
  `mass_points/summary.json`): 1 × 10¹² M☉ at all 19 centres excluded (share above observed 0.997; one MACS0417 centre undefined at that mass); 5 × 10¹¹ not
  excluded at any f (≤ 0.618). **Handoff:** the W1 radial-centre stack line is done.

## 2026-10-17: W1 lenses injected at the `radial` centres through the stack: limit ≥ 2 × 10¹² M☉ (post hoc rule)
- `scripts/w1_stack_inject.py` (`ShearInjector.lensed_sample` split out of `trial`); plan frozen first
  (`results/w1_stack_radial/inject/PLAN.md`). Gates Z (θ_E = 0 reproduces S) and R (stack responds) passed.
- **Failed approach:** the frozen rule adds lenses to the observed shapes only, so any positive shift "excludes"
  (2 × 10¹⁰ M☉ moves the stack 0.51 → 0.53 and is "excluded"); no limit quoted from it.
- **Post hoc noise-aware rule** (ΔS added to rotation-null draws; `RESULT.md`): excluded ≥ 2 × 10¹² M☉ at all 19
  centres, ≥ 8 × 10¹² at half; ≤ 2 × 10¹¹ untestable (share above observed ≤ 0.28). The p 0.025 excess stays an
  exploratory note. **Handoff:** W1 stack line done unless a mass point between 2 × 10¹¹ and 2 × 10¹² is wanted.

## 2026-10-17: W1 stack at `radial` centres, second model for MACS0417: S unchanged where both models cover
- `scripts/w1_stack_radial.py --shear-model` (centres and arc exclusion fixed, only the cluster-shear removal swaps).
- RELICS Lenstool vs CANUCS maps on the 7 MACS0417 centres valid under both: mean S −0.12 vs −0.18, largest change
  0.30 (`results/w1_stack_radial/model2/common_centres.json`). The MACS0417 excess (CANUCS, 14 centres, p 0.061) sits
  in the 7 centres outside the RELICS area. **Handoff:** injections at the centres through the stack for a limit;
  the cluster-model-error explanation stays open only outside the RELICS area.

## 2026-10-17: W1 shear stacked at the `radial` centres: excess mostly circular; residual p 0.027 (exploratory note)
- **Why:** the `radial` centres (arc-orientation signature) had never been checked against the W1 shear statistic of
  other sources (TASKS JWST focus: "stack S around `radial` ... centres").
- **New:** `scripts/w1_stack_radial.py` (S, S_× at the centres via `exotic_screens.ApertureMass`, rotation null,
  `--exclude-radial-arcs`, `--exclude-near-arcsec`); plan frozen first (`results/w1_stack_radial/PLAN.md`).
- **Result (`results/w1_stack_radial/RESULT.md`):** 19 valid centres (S1063 5, MACS0417 14); planned stack mean S 0.72,
  p 0.0035. Post hoc vet: without the 85 arcs that define the centres, 0.51, p 0.027 (exploratory note under the frozen
  rule); without their 2″ neighbours too, 0.48, p 0.071 (diagnostic only). Not a candidate.
- **Failed approach:** the planned sign check (cluster shear kept, off-centre apertures) barely responds (−0.007 to
  −0.15); use a synthetic or injected case. **Handoff:** next is the same stack with a second lens model per field
  (model-error explanation), then injections at the centres for a limit.

## 2026-10-17: Radial (W1) centre-level axial rule for streaks: gate K failed, not adopted
- **Why:** handoff of `results/radial_streak`: an arc-level mask veto thins the null and misses the faint MACS0647
  streak; a centre-level rule applied alike to data and null avoids the thinning.
- **Rule tested (frozen in `results/radial_axial/PLAN.md`):** axial concentration A = |mean exp(2iφ)| of the
  contributing arcs' azimuths about the centre; A ≥ 0.9 → the centre does not count.
- **Result (`results/radial_axial/RESULT.md`):** MACS0647's 6-line centre (reproduced, p 0.075) has A = 0.731, so
  gate K fails and the rule is not adopted (injections not run). Sheet looked at: five of six arcs sit on the faint
  streak, but the grid centre lies up to 1″ off it, so their azimuths spread over ≈ 20°, and one arc is off it.
- **Failed approach:** azimuth-about-the-grid-point as a streak test. Post hoc: 5 of 6 contributors lie within 0.5″
  of one line. **Handoff:** a collinearity rule needs a new frozen plan, a different known streak (or injected
  streaks) as its gate, and W1 injections; `radial` conclusions unchanged.

## 2026-10-17: VENUS batch 3: 20 new boxes, 10 calibrated; two rule-level changes in a failed A370 box, neither confirmed (#311)
- **Coverage:** 20 new 24″ boxes in 10 fields, ranked to prefer extra public epochs at similar depth (owner
  instruction); 10 passed calibration (1.6 arcmin²), 10 failed (2 false hits, 8 vanish recovery); 0 timed out.
- **Screen:** 136 hits → 22 unflagged → 9 review entries, all reviewed: 3 no significant change, 4 undetermined, 2 real
  change by rule, both in calibration-failed A370 box 2 (`docs/candidates/batch_venus3_A370_box2.md`):
  - h003 is a compact z 0.36 galaxy; the F444W drop is not confirmed by an epoch 18 d later (5324);
  - h010 sits on the Dragon arc next to 33 catalogued microlensing detections; no point source is seen.
- **Tools:** `batch_pairs.py --prefer-extra-epochs`; `batch_run.py` per-task timeout and stop time (portable without
  `os.wait4`). Results: `results/transients/batch_venus3/RESULT.md`.
## 2026-10-17: Radial (W1) screen with the W1 streak veto: gate failed, not adopted
- **Why:** handoff of the radial_more run: MACS0647's strongest `radial` centre is built from streak segments.
- **New:** `radial_more.py --streak-masks` (`exotic_screens.drop_streak_arcs`: candidate arcs within 1″ of a
  `findsat_mrt` trail dropped after selection; the null uses the same arcs). Masks regenerated; pixel counts match
  `results/w1_streak_veto` in 12 / 12 fields; the no-veto runs reproduce `results/radial_more`.
- **Result (`results/radial_streak/RESULT.md`):** gate K fails: the mask covers 2 of the 6 arcs; the centre stays
  with 5 lines, three on the faint streak the mask misses (sheet looked at), and the thinner null puts it at p 0.03,
  a veto-made flag in a no-result field (not a candidate). Twelve-field diagnostic: no flag changes in the ten
  control-passing fields; Abell 370 loses 22 of 100 arcs, SMACS max lines 4 → 3.
- **Failed approach:** `findsat_mrt` masks as a `radial` streak veto. **Handoff:** radial conclusions unchanged; a
  radial streak veto needs a faint-streak detector validated on the MACS0647 streak plus injections.

## 2026-10-17: Radial (W1) screen on MACS0417, MACS1423, WHL0137, RXJ2129, MACS0647: no flag
- **Why:** the five W1 shear-only fields had never run the `radial` screen (exotic W1 signature, D-031/D-049).
- **New:** `scripts/radial_more.py` (the shear fields' inputs through `exotic_screens.cmd_radial`, Gaia spike veto;
  `--sheet` draws the strongest centre with its arcs).
- **Known case:** Abell S1063 headline run reproduced exactly (46 not model-radial, 10 centres, max 4, p 0.435).
- **Result (`results/radial_more/RESULT.md`):** 0 centres with p_random < 0.05. MACS0417 max 5 lines (p 0.17), MACS1423
  3 (0.915), WHL0137 3 (0.705, repeats the D-036 run with the Gaia veto). RXJ2129 and MACS0647 give no result (models
  fail the `arcs` control); MACS0647's 6-line centre (p 0.075) is five segments of one linear streak (sheet, by eye).
- **Handoff:** `radial` now covers all twelve W1 fields; its D-049 limits are unchanged (no new injections; the shear
  limits are 9–14× stronger). Next radial item: the low-weight veto, and a streak veto (W1 masks) for MACS0647-like
  streak segments.

## 2026-10-17: MACS0417 h002 delays re-checked with a fixed source position (#309)
- #308 evaluated the Fermat potential with each pixel's own β; fixed β0 plus root-refined images gives 748 d (best fit)
  and 708–778 d (10-realisation sample range) for the 7.6″ counter-image (#308: 650–830 d, partly from a non-image grid
  candidate). Self-delay 0; order-1 vs order-3 interpolation within 3.6 d.
- Conclusion unchanged: not testable with current epochs; next useful epoch 2027-09-14 to 2027-11-22
  (`results/transients/batch_venus2/followup/macs0417_h002/recheck/RESULT.md`).

## 2026-10-17: W1 calibrated streak veto (STScI findsat_mrt masks, through injections): adopted; twelve-field limits 4.4e3 / 7.2e2 / 3.5e2 deg⁻²
- **Why:** handoff of `results/w1_shear_rxj_m0647/RESULT.md`: a streak veto run through injections (S1063, MACS1423,
  MACS0647 peaks on streaks); the post hoc V1 test could not be used for a limit.
- **New:** `scripts/w1_streak_mask.py` (streams the F200W SCI, `acstools.findsat_mrt` persistent trails),
  `inject_shear.py --streak-masks` (`shear_screen` gets a `veto` row flag), `scripts/w1_stage_inputs.py` (stages
  the twelve fields' inputs), extra `streaks` in pyproject.
- **Result (`results/w1_streak_veto/RESULT.md`):** K1 and K2 hold, so the veto is adopted. With it, MACS1423 gives S_max 3.41 and p_rot 0.22
  and no longer passes. No field passes. MACS0647 loses 25 % of its used rows (named). Twelve fields, 70.5 arcmin²:
  Σ(W1) < 4.4 × 10³ / 7.2 × 10² / 3.5 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉.
- **Failed approaches:** zeroing > 5σ pixels before the trail finder removes bright streaks (PLAN addendum A).
  Dropping vetoed rows before the screen changed the spike veto too, and S1063 then passed spuriously at p_rot
  0.045. Code review found it; the run was superseded.
- **Second review (final diff):** fixes that change no number (deterministic mask files, veto checks, a
  veto-path test); traced items in RESULT "Second code-review".
- **Handoff:** W1 next: faint parallel streaks (S1063 S_max aperture) are missed: a lower findsat threshold,
  calibrated through injections; fewer false trails on MACS0647; the open review items in RESULT.
## 2026-10-17: Batch-2 follow-up: AS1063 h002 edge, MS0451 h011 detector-fixed, MACS0417 h002 counter-images not testable yet (#308)
- **AS1063 1 h002:** artifact (edge, 12.54″ from the box centre); no significant change in any mosaic either.
- **MS0451 1 h011:** artifact. The signal sits on one NRCALONG pixel in all 18 exposures of visits 1–2 (S/N 9–28),
  absent at the sky position in other dithers; mechanism (persistence or hot pixel) not distinguished.
- **MACS0417 2 h002:** public RELICS glafic/lenstool models reproduce the catalogued images; counter-images lag by
  ≈ 650–830 d (glafic best fit and 10 realisations), so current epochs cannot test them. Next useful epoch ≈ 2027-07 to 2028-01.
- Results: `results/transients/batch_venus2/followup/RESULT.md`.

## 2026-10-17: W1 shear screen on RXJ2129 and MACS0647: null; eleven-field limits 4.5e3 / 7.6e2 / 3.8e2 deg⁻²
- **Why:** TASKS "JWST focus" 1 (W1): the next fields of the previous handoff (DJA v7.0 photo-z).
- **Result (`results/w1_shear_rxj_m0647/RESULT.md`):** neither passes D-053 (S_max 2.91 / 2.97, p_rot 0.375 / 0.325).
  R 0.551 (RXJ2129, in the earlier range) and 0.231 (MACS0647, lowest of all fields; cause not measured). Sheet looked
  at: both MACS0647 top centres sit on two crossing linear streaks (by eye), the third field with this pattern.
- **Injections:** RXJ2129 8 / 53 / 92 of 200 at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉; MACS0647 0 / 1 / 1 (close to blind).
- **Limits:** eleven fields, 64.0 arcmin²: Σ(W1) < 4.5 × 10³ / 7.6 × 10² / 3.8 × 10² deg⁻² (was 4.9e3 / 8.2e2 / 4.1e2);
  without MACS0647 < 0.2 % different. `--combine-only` reproduces the nine-field file exactly.
- **Handoff:** W1 next: the calibrated streak veto (now S1063, MACS1423, MACS0647), run through injections; then MACS1423
  back in; then more photo-z fields (MACS0717 once DJA has photo-z).

## 2026-10-17: W1 shear screen on MACS0417, MACS1423, WHL0137: one hit, carried by sources on an image stripe; nine-field limits 4.9e3 / 8.2e2 / 4.1e2 deg⁻²
- **Why:** TASKS "JWST focus" 1 (W1): more fields with photo-z (MACS0717 still has none in the DJA v7 index).
- **Result (`results/w1_shear_more/RESULT.md`):** R 0.488 / 0.604 / 0.533. MACS0417 and WHL0137 not passed (S_max
  3.56, p_rot 0.065; 2.75, 0.72; both maxima at mosaic edges, sheet looked at). MACS1423 **passed** the D-053 rule
  (S_max 3.88, p_rot 0.01, > B max 3.64); post hoc streak/star test frozen before computing (addendum A, `scripts/w1_peak_vet.py`): the
  15 sources on a vertical image stripe through the aperture carry it (S 3.88 → −0.07 without them; stripe origin
  not determined). Not a candidate.
- **Limits (deviation from the plan):** nine fields (MACS1423 left out: hit plus R 0.60 above the earlier range), 55.1 arcmin²: Σ(W1) < 4.9 × 10³ /
  8.2 × 10² / 4.1 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉; ten fields 4.3 × 10³ / 7.0 × 10² / 3.5 × 10² only
  once a streak veto is run through injections.
- **New:** shear-only field entries in `inject_shear.SHEAR_FIELD_OVERRIDES` (`FIELD_NAMES`); `scripts/w1_peak_vet.py`.
- **Handoff:** W1 next: a calibrated stripe/streak veto in the screen (now two fields, S1063 and MACS1423), re-run
  injections, then the ten-field limit; RXJ2129 and MACS0647 (DJA v7.0 photo-z) are the next fields.

## 2026-10-17: W1 shear screen on SMACS 0723 and El Gordo: null; seven-field limits 5.3e3 / 9.6e2 / 5.0e2 deg⁻²
- **Why:** D-053 "Revisit if": the DJA photo-z tarballs (404 on 2026-10-08) answer HTTP 200 again; sha256 match.
- **Result (`results/w1_shear_smacs_elgordo/RESULT.md`):** R 0.438 / 0.486 (cluster shear seen); S_max 2.91 / 2.95,
  p_rot 0.77 / 0.62, below the rotation-null median: null. Injections 11 / 62 / 106 (SMACS) and 0 / 14 / 48 (El
  Gordo) of 200 at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉. Seven fields, 43.7 arcmin²: Σ(W1) < 5.3 × 10³ / 9.6 × 10² /
  5.0 × 10² deg⁻² (was 5.8 × 10³ / 1.1 × 10³ / 5.7 × 10²). S_max cutouts looked at (`sheet_smax.png`).
- **New:** `SHEAR_FIELD_OVERRIDES` for both fields (Gaia veto); `DEFAULT_FIELDS` now seven; `--combine-only` reads screen
  options missing from older summaries as their legacy values (Abell 2744 re-run reproduces exactly);
  `scripts/w1_smax_sheet.py`.
- **Handoff:** W1 next: MACS0717 needs a photo-z; streak veto (S1063, and the SMACS S_max aperture has spike
  streaks); PSF-weighted moments for the floor (reuse-check GalSim HSM first).

## 2026-10-17: W1 STPSF + detection-kernel/jitter/drizzle smoothing vs the star χ: rule S3 fails; PSF-model route stops
- **Why:** `results/w1_stpsf` "Limitations": STPSF over-predicts the stars' area-matched |χ| by 0.05–0.07; the
  jwst catalogue measures morphology on kernel-smoothed data (FWHM 2 px default; ASSUMPTION for this catalogue).
- **New:** `scripts/w1_psf_smooth.py` (five pre-listed variants); STPSF data streamed (not stored).
- **Result (`results/w1_psf_smooth/RESULT.md`):** the unsmoothed variant reproduces #301. Smoothing cuts the real-
  component excess by 0.016–0.023 (stars 623, 7644) and 0.06–0.07 (star 8291), but S3 (all three stars within 0.05
  per component) fails in every variant: star 8291's imaginary component stays 0.10 off in all of them (orientation, which
  isotropic smoothing cannot change); star 7644 passes only with the added 0.010″ jitter (−0.048). W1 limits unchanged.
- **Failed approach:** a smoothed STPSF model as the per-source χ_PSF forward model.
- **Handoff:** next W1 floor step, if any, is PSF-weighted pixel moments (reuse-check GalSim HSM first); otherwise
  accept the 4.21 floor and spend W1 time on new fields.

## 2026-10-17: W1 model PSF (STPSF) vs the star χ_PSF and Abell 2744's B-mode floor: axis close, size level-dependent; not adopted
- **Why:** TASKS "JWST focus" 1 (W1), `results/w1_psf` "Next". Rules S1/S2 frozen in `results/w1_stpsf/PLAN.md`.
- **New:** `scripts/w1_stpsf_floor.py` (`psf`, `stars`, `screen`); STPSF data streamed (88 MB, not stored).
- **Result (`results/w1_stpsf/RESULT.md`):** model χ at sky PA 125.4–128.3° (stars 133°). Pixel-frame
  orientations: 2 of 3 stars within 1.2° of the model, one 16.7° off. The faint isophotes follow the horizontal
  spike pair (figure looked at). S1 (deviation: grid median) is consistent at 2 of 6 isophote levels (|χ| 0.18–0.21);
  model |χ| ranges 0.008–0.209 over the levels and is not monotonic. Area-matched stars (post hoc): model |χ|
  exceeds each star's by 0.05–0.07. S2 with a constant χ per level (deviation): b_max 3.57–3.97; the two passing
  levels cost 11–12 % of the sources, and at 3e-4 the E-mode maximum exceeds the null's 95th percentile. Not
  carried to injections: a constant correction would over-correct fainter sources. W1 limits unchanged.
- **Failed approach:** a constant PSF χ from the model (level-dependent, over-predicts the stars).
- **Handoff:** next W1 floor step is a per-source forward model (source ⊗ drizzled STPSF PSF at the catalogue
  isophote, validated on the three stars first) or PSF-weighted moments.

## 2026-10-17: VENUS batch 2: 20 new boxes, 12 calibrated; three real changes (ordinary-transient hypotheses), identities open (#304)
- **Coverage:** 20 new 24″ boxes (3.2 arcmin²) in 11 fields, metadata-only selection with the batch-1 rule; 12 passed the
  per-box calibration (1.92 arcmin²), 8 undeterminable (`results/transients/batch_venus2/RESULT.md`).
- **Screen:** 203 hits → 56 unflagged → 11 review entries; 10 reviewed, 1 pending. Final classes: 3 real changes
  (MACS1423 fading, MACS0417 brightening on lensed image 1.1, Bullet faint brightening), 3 undetermined, 3 artifacts,
  1 no significant change. TNS not queried (403). Vetting: `docs/candidates/batch_venus2_real_changes.md`.
- **Fix:** `batch_pairs.overlaps()` lets touching boxes through (old test rejected 13 of 40 batch-1 neighbours; `overlap_fix_check.json`).
- **Failed approach / cost:** PLCK G165 2 hit the per-stamp fallback (stamp pixels > `max_pix`): 38 min injection,
  screen stopped at 83 min (post hoc 12:45Z cutoff); undeterminable by calibration anyway.

## 2026-10-17: h002 identity follow-up and h003 limited re-check (#300)
- **h002:** real fade; an ordinary transient is reasonable; identity and type not confirmed. No reported transient
  matched in queried sources; TNS not queried (403). SN H0pe images 24.9–37.9″ away; NS_0481 (3.77″) is a separate
  galaxy. HST 2016 F160W 3σ 28.5 nJy (most likely appeared after 2016; indicative, different instrument and aperture). Only F090W declines significantly over 23.8 d
  (`results/transients/batch_venus/followup/h002_identity_followup.json`).
- **h003:** large differences are same band, position and aperture, also within programme 1176; 31 matched arc-pixel
  controls give tail fraction 0.39 on the decisive F356W pair: consistent with matched arc-pixel scatter (pre-fixed
  rule met); no mechanism isolated, so it stays undetermined; closing it is the owner's call (`followup/h003_recheck.json`).
- **Wording:** dated revisions in `RESULT.md` and the h002 record; review queue empty ≠ all objects explained.
## 2026-10-17: W1 PSF anisotropy from stars vs Abell 2744's B-mode floor: 4.21 → 3.85, misses the rule; not adopted
- **Why:** TASKS "JWST focus" 1 (W1), `results/w1_floor` "Next". Frozen rules R1/R2 in `results/w1_psf/PLAN.md`.
- **New:** `exotic_screens.psf_moments` and `psf_chi` in `intrinsic_ellipticity` / `shear_sources` /
  `shear_screen`, `inject_radial.lensed_cov` / `paint_lens`, `inject_shear.py --psf-chi` (default isotropic,
  unchanged); `scripts/w1_psf_floor.py`.
- **Result (derived; `results/w1_psf/RESULT.md`):** 3 PSF stars (Gaia-matched), median χ_PSF with |χ| 0.157 at
  PA ≈ 133°. R1 (star map at the 10″ scale) fails: median nearest star 92″. Constant correction: b_max 4.21 → 3.85
  at (88, −42) (rule ≤ 3.79: fails), used sources 4,951 → 4,496, field-mean ε₂ −0.026 → −0.006; E-mode still null.
  W1 limits unchanged.
- **Failed approaches:** a star-based PSF map (too few stars); a trace-preserving PSF model and round-PSF
  injection painting (code review; first-run numbers withdrawn).
- **Handoff:** next W1 floor step is a model PSF per detector (STPSF; reuse-check first).

## 2026-10-17: PLCK G165 h002 vetting record corrected: the z = 1.32 galaxy is 3.77″ away, not 0.06″
- **Error:** the #296 record read NED's separation (0.063, in arcmin) as arcsec and wrote that h002 "sits on" the
  z = 1.3236 galaxy `[FPP2024] NS_0481`. The separation is 3.77″ (`results/transients/batch_venus/h002_identity.json`).
- **Checked (literature subagent, 2026-10-17):** NS_0481's z is spectroscopic (Frye et al. 2024 Table 2, NIRSpec
  DDT 4446); Frye et al. report only SN H0pe in G165, whose three images are 25–38″ away (not stored); TNS returned
  HTTP 403. Host and redshift of h002: not determined.
- **Verdict unchanged:** known population, supernova-like transient (hypothesis; the light curve carries it).

## 2026-10-17: W1 B-mode floor diagnostic on Abell 2744: no streak; carried by many sources near bright galaxies
- **Why:** TASKS "JWST focus" 1 (W1): find what sets the highest floor (Abell 2744 max |S_×| 4.21 at (88, −42)).
  Claimed as a streak veto; the sheet showed no streak, so it became a diagnostic (deviation; nothing adopted).
- **Result (derived; `results/w1_floor/RESULT.md`):** no streak in the box (sheet looked at); top 1/3/5 of 44
  sources give 0.12/0.31/0.47 of S_×; aperture mean ε 0.044 (field mean |ε| 0.331); 1 of 200 rotation-null maxima
  reaches 4.21 (indicative). A bright-neighbour veto lowers b_max to 3.86 (k 3) or 3.89 / 3.47 (k 5, Δm 2 / 3)
  only with a 31–74 % source loss and flagged rows are not over-represented at the peak: not adopted.
- **Failed approaches:** collinear-chain streak proxy (21 % of sources flagged by chance at 10°/10°/20″);
  bright-neighbour rule as a floor-specific veto.
- **Handoff:** next W1 floor test is a spatially varying systematic (star-based PSF map or a bright-galaxy
  light-gradient term per source), then injections under a frozen adoption rule.

## 2026-10-17: W1 shear screen aperture-coverage cut (c = 0.5): not adopted; limits move ≤ 1.6 %
- **Why:** TASKS "JWST focus" 1 (W1): "drop edge apertures with < 50 % coverage" to lower the B-mode floor. Plan
  and computable rules first (`results/w1_edge_cut/PLAN.md`, `scripts/w1_edge_cut.py decide`).
- **New:** `exotic_screens.aperture_coverage`, `--min-coverage` on `shear` and `inject_shear.py` (default 0,
  unchanged behaviour); `ApertureMass(valid=)` blanks centres in the map, the rotation null and injection trials.
- **Result (derived; `RESULT.md`):** c = 0 reproduces the committed recoveries in all five fields. At c = 0.5 only
  MACS0416's b_max moves (4.04 → 4.01); Abell 2744's 4.21 floor sits at full coverage. Five-field limits
  5,765 → 5,672 / 1,082 → 1,088 / 567 → 570 deg⁻² (2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉): the frozen adoption rule
  (lower at each mass) fails. Real fields stay null. **Failed approach:** edge apertures do not set the floor.
- **Handoff:** the floor is elsewhere; next W1 gains are a PSF-anisotropy model from stars, blend rejection or a
  streak veto, tested against Abell 2744's 4.21 centre (88, −42) first.

## 2026-10-17: bounded VENUS batch, 20 new boxes with fixed `arc_diff`: one real fade, explained as an ordinary transient (#296)
- **Data:** 20 unsearched 24″ boxes (3.2 arcmin²) in 10 clusters; earlier NIRCam epoch → VENUS, F150W + F444W.
  - Chosen from metadata only: coverage, finite fraction, depth.
  - Depth ties broken alphabetically (post hoc label).
- **Calibration per box (fixed rule):** 12 passed (1.92 arcmin²); 8 undeterminable.
  - Pooled over the passing boxes, vanish recovery is 0.87 / 0.82 / 0.79 at 26.0 / 27.0 / 27.5 AB.
  - Partial (50 %) fades: 0.79 / 0.75 / 0.52.
  - Unchanged-control hits: 8/207.
- **Screen:** 113 hits, 13 unflagged, 5 review entries, all 5 reviewed (capacity 10, 0 pending).
  - Three show no significant change in forced photometry. h003 is undetermined after the full-rule post-review
    check: mixed-sign differences across programmes on a bright arc.
  - **PLCK G165 h002** is a real fade: a point source in all PEARLS bands in March to April 2023, gone in 2025 and
    2026. That light curve is an ordinary transient's (supernova-like).
  - Vetted (`docs/candidates/batch_venus_PLCKG165_h002.md`): known population (hypothesis). Its TNS identity is
    pending. (Corrected 2026-10-17: the z = 1.32 galaxy NS_0481 is 3.77″ away, not under it.)
  - Not an unexplained signal and not wormhole evidence.
- **Time:** 22 min of computation on 4 workers; about 60 min for the whole unit (not stored).
- **#293 report wording corrected:**
  - vanish 51/60, not 51/51;
  - 180 changes + 60 controls;
  - per-source vs totals;
  - throughput vs per-box time.

## 2026-10-17: SFFT vs `arc_diff` (bounded tool comparison): keep `arc_diff` (D-095, #293)
- **The JWST B-spline SFFT** (Hu & Wang 2024) is GPU-only in sfft 1.7.3, so it was not tested. SFFT's automatic
  masking needs SExtractor, so it was not tested either.
- **Ordinary SFFT** (CPU, unmasked) replaced only `kernel_match` on B1 + B2 with the same pre-matching injections.
  - After two rounds of post-review adapter corrections (union NaN mask over the kernel footprint, propagated variance,
    sign accounting), against `arc_diff`:
    - fades on the review list: 83/120 vs 92/120 (partial 32 vs 41; vanish 51/60 each, not the same sources);
    - brightenings: 42/60 vs 45/60;
    - unchanged false hits: 5/60 vs 1/60.
  - About ×1.4 faster (indicative; loads differ), with 1.49–2.51 × the peak memory.
  - **Verdict under the fixed rule:** keep `arc_diff`.
- **Failed first adapter (kept in `runs_v1/`):** its non-union NaN mask and unpropagated variance cost SFFT about a
  third of the fades.
- **Batch denominator:**
  - 56 VENUS targets, 18 with another public NIRCam F150W epoch, about 31 pairs (an upper bound: F444W and overlap
    not checked);
  - about 4 min per box per core; a throughput of about 1 box per minute with 4 boxes in parallel;
  - about 0.3–1 review entries per box.
## 2026-10-17: W1 shear screen on Abell S1063: null; five-field W1 limits 5.8 × 10³ / 1.1 × 10³ / 5.7 × 10² deg⁻²
- **Why:** TASKS "JWST focus" 1 (W1): S1063 has a pinned model and DJA photo-z but had only the radial screen. Plan
  first (`results/w1_shear_s1063/PLAN.md`).
- **Result (derived; `RESULT.md`):** R 0.472 ± 0.034; S_max 3.61 (p_rot 0.07) > max |S_×| 3.33, but fails p < 0.05:
  **null** under D-053. Both top-peak apertures are crossed by long straight streaks in F200W (sheet looked at).
- **Injections:** 0 / 0 / 17 / 69 / 91 of 200. The four earlier fields re-run here reproduce their counts and limits
  exactly (post hoc, `rerun_4fields/`). Five fields (35.8 arcmin²): Σ(W1) < 5.8 × 10³ / 1.1 × 10³ / 5.7 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² /
  2 × 10¹³ M☉ (was 7.7 × 10³ / 1.3 × 10³ / 6.5 × 10²).
- **Code:** `inject_shear.SHEAR_FIELD_OVERRIDES` (S1063 photo-z and Gaia veto for the shear run only);
  `inject_radial.load_field(spec=)` and `combine(specs=)` so the photo-z check uses the shear inputs. Tests added.
- **Handoff:** next W1 gains: a streak veto (long linear features beyond the Gaia spike veto) to lower the E/B floor;
  MACS0717 when photo-z exists; SMACS/El Gordo photo-z (404).

## 2026-10-17: `arc_diff` per-stamp kernel fits share one design (D-094): B2 run ×14 (screen), same results (#291)
- **Diagnosis:** each of B2's 1,076 F150W stamps rebuilt an ~86,800 × 82 design and ran 3 SVD least-squares on it
  (≈1.4 s per stamp).
- **Change:** `kernel_fit_stamps` builds the design and Gram once per band. Each fit subtracts its own stamp and
  clipped rows, or sums them directly when they dominate, and refines its solution until it converges.
  - Leave-one-out, per-fit clipping, iterations and precision are unchanged.
  - It is the default; `--kernel-stamp-fit refit` keeps the old path.
- **Failed first version (review):** plain Gram downdating lost digits for a bright-star stamp (synthetic: up to
  1.2e-4 relative; review measurement, not stored). Fixed with the guard and refinement, and checked on every
  stamp.
- **Validation (tolerances fixed first):**
  - all B1 and B2 stamp kernels agree to ≤ 1.2e-13, with the same clipping;
  - B2 S ≤ 3e-8;
  - hits, flags, morphology, review list and all 60 injection outcomes identical;
  - no threshold crossings or near-threshold hits.
- **Speed (B2, no cache; similar, not identical load, timed in separate batches):**
  - whole run 1,508 → 106 s (screen, ×14) and 1,341 → 157–176 s (injection run, ×8–9);
  - matching ×16–17; stamp kernels ×28–30;
  - peak memory unchanged.
- **Wording corrected in `fade_calib/RESULT.md` (owner):**
  - the direction of the injection-registration effect is not verified;
  - the 8 looked-at rejections do not show that no fade was missed.
## 2026-10-17: SMACS ERO → VENUS `arc_diff` with tile-median registration: shift 2.4 px off, outcome unchanged
- **Why:** pending item of #288: pair 1's F150W whole-box shift (2.11, −1.82) px was never compared with
  `--reg-tile-arcsec 6`. Plan first (`results/transients/smacs_arcdiff/o001_o057_tile/PLAN.md`).
- **Result:** tile medians F150W (0.20, −0.40), F444W (0.28, −0.09) px: the whole-box peak was pulled (2.4 / 1.0 px),
  as in El Gordo. Re-screen: Atalanta A/B recovered as `known`, 0 `vet`, no new pass within 2″ of an image; the two
  new `off_image` passes are core residuals on the sheet (looked at). Injection recovery at 27.5 AB 0.85 / 0.82
  (was 0.82 / 0.76).
- **Handoff:** item closed; `results/transients/smacs_arcdiff/RESULT.md` pair 1 stands, h061 stays open as recorded.

## 2026-10-17: MOA-II manifest lists every streamed field tar; retrieval dates no longer hard-coded
- **Why:** TASKS follow-up from the gb10 review: `data/manifests/moa_ii.ecsv` had only the metadata tar and gb22,
  and `write_manifest` stamped every row `retrieved_utc` 2026-10-08 (the day of the HTTP HEAD reads, not of the
  streaming).
- **New:** `w3_moa.py` pre-screen chunks record `retrieved_utc` (when their last byte range arrived);
  `field_retrieved_utc` turns them into a day or `first/last` range; the manifest gains `head_utc` (2026-10-08, the
  size and Last-Modified read). Test `test_field_retrieved_utc_from_chunk_stamps`.
- **Result:** manifest regenerated from the tracked chunks: 22 streamed tars with `range_sha256` (gb1–gb22; gb22
  also keeps its whole-file sha256 row). Their `retrieved_utc` is "not recorded": the older chunks carry no stamp,
  and git dates of the chunk files (all 2026-10-15) are not retrieval dates, so none is inferred.
- **Handoff:** no science change. A re-stream of any field fills its dates.

## 2026-10-17: `arc_diff` fading calibration (pre-matching injections), registration fix, El Gordo B1–B3 search (#288)
- **Scope (owner, 2026-10-16):**
  - calibrate fades, then a measured speed-up, then a small batch only if the calibration passes;
  - plan first: `results/transients/fade_calib/PLAN.md` with addenda 3–5.
- **Tools:**
  - `arc_diff.py --inject-pre`: sources go into both raw epochs before registration and matching. Epoch 2 follows
    the chain's injection-free registration.
  - `--inject-pre-psf empirical`: per-epoch PSFs.
  - **The one chain fix:** `--reg-tile-arcsec 6`, the global shift as a tile median. El Gordo's whole-box shift was
    pinned at ±3 px.
  - `--cache-dir`: raw-data products only.
  - `scripts/fade_calib.py`: `make`, `eval`, `boxes`.
- **Failed approach, kept for the record:**
  - The first injection placed epoch-2 sources at their WCS position, so after the chain's registration shift each
    unchanged source became a 0.9 px dipole.
  - Held-out seeds 20 + 21 "failed" (false hits 0.25–0.28) because of that artefact.
  - Found by `/code-review`; corrected in addendum 4. Fresh held-out seeds 30 + 31 were used after it.
- **Calibration (box 1, held-out seeds 30 + 31; n = 12–14 per cell):**
  - vanish 26.0 AB on the review list with the sign: 14 / 14 [0.93, 1.00];
  - unchanged false hits: 1 / 40 (empirical PSF) and 4 / 40 (Gaussian, exactly at the 0.10 limit); passes §4;
  - corrected dev seeds: 6 / 40.
  - Partial fades at 27.5 AB: detected 5 / 14 (empirical) and 8 / 14 (Gaussian); on the review list 4 / 14 and
    6 / 14.
  - Known limitation: injected sources keep a 0.02–0.15 px F150W misregistration
    (`inject_shift_residual.json`); not corrected, no rerun.
  - These are small-sample figures for one box, and an injected-source false-hit fraction is not a real-sky
    false-alarm rate.
- **Batch:**
  - B1 (#286) re-screened: review list 4 → 1 (the vetted h024 position); #286 stands.
  - B2 and B3 were chosen by coverage before any difference. Review lists 0 and 0, and all 8 rejected hits were
    looked at.
  - Local vanish 26.0 AB efficiency: 5 / 7 and 6 / 7.
  - Nothing pending; no exclusion.
- **Speed-up (`speedup.json`):** a warm cache gives ×1.9 (empirical) and ×1.3 (Gaussian) on box 1.
  - The outputs are byte-identical, and the first use is not slower.
  - There is no gain in B2 and B3, where the per-stamp kernel fits take 1,466–1,948 s (1,065–1,199 stamps against
    217); that step is left for a decision.
- **Pending (not done here):**
  - the SMACS ERO → VENUS whole-box shift (2.11, −1.82) px, unchecked;
  - the kernel-stamp cost.
## 2026-10-17: `arcs` on WHL0137 (Sunrise), RELICS maps: the S1063 anti excess does not replicate (1 vs 0.15, P 0.14)
- **Why:** the S1063 replication tally needs fields never used to tune `arcs`; WHL0137 has an F200W catalogue, DJA
  v7.5 photo-z and RELICS maps (the HEAD probe of DJA v7 names found only Sunrise, SMACS, MRG0138 and RXCJ0600 among
  untested clusters; the last two lack public maps). Plan frozen first (`results/arcs/whl0137-relics/PLAN.md`).
- **New:** `whl0137-relics-glafic` model (sha256-pinned RELICS v1 maps); `arcs_rxj2129.decide` factored out of its
  `main` and reused by `scripts/arcs_whl0137.py`; `scripts/arcs_whl0137_sheet.py`.
- **Result (derived; `RESULT.md`):** Lenstool 13 / 18 aligned (control passes, 72 %), anti 1 vs 0.150, P 0.139;
  GLAFIC 14 / 19, anti 1 vs 0.232 → **frozen outcome "does not replicate"**. The anti (2194, m 24.4, g 0.24) is an
  edge-on-looking disc with the same axis in F200W and F444W (sheet looked at): intrinsic shape; not a candidate.
- **Handoff:** S1063 replication tally: MACS0417 no, MACS0647 no result, MACS1423 explained, RXJ2129 no result,
  WHL0137 no. No untested public field with maps + F200W + DJA photo-z is left in the probe list (SOURCES.md);
  S1063 stays a single-field exploratory result.

## 2026-10-17: MACS1423 clump-pair census re-run with the fixed #273 code: unchanged (P 0.259); citable (#287)
- **Why:** the 2026-10-16 flag (below): #282 ran the census with the pre-fix `arcs_vet_deblend.fit_models`. Plan
  first: `results/arcs/macs1423-canucs/vetting/PLAN_census_rerun.md` (only the code changes).
- **Result (derived; `RESULT.md` "Step 2 re-run"):** every sample number unchanged (45 / 55 gate-passed, f̂ 0.289,
  expected anti 2.608, P(≥ 4) 0.259, P(≥ 3) 0.494); outcome `explained by clump orientation (exploratory)`. Only
  1818 changes, `single` → `undecided` (old fit pinned at S/N 3 × 10⁻¹¹); its gate still fails, so P does not move.
  No `gate_passed` flips. A second run gave identical output; sheet regenerated and looked at.
- **Handoff:** the MACS1423 census may be cited. #279 (five-field census) has its own re-run, awaiting the owner.

## 2026-10-16: El Gordo arc difference imaging, PEARLS 2022 → VENUS 2026 (one-box pilot): nothing open (#286)
- **Scope (owner):**
  - one box and one epoch pair;
  - the #277 recovery ≥ 0.50 rule is waived for this pilot only (the failure record stands; no wider scan);
  - exploratory, not blind.
- **Plan first** (`results/transients/elgordo_arcdiff/PLAN.md`):
  - 24″ box at (15.7300, −49.2530) with 11 Caminha+2023 images;
  - F150W + F444W, jwst 3.0.0 in both epochs;
  - the #275 / #277 settings, unchanged.
  - `arc_chain.py` gains an `elgordo` image list.
- **Sensitivity:**
  - full-chain brightenings 23 / 27 = 0.85 at 26.0 AB (control gate ≥ 0.80 passed), 0.78 at 27.5, 0.54 at 28.0;
  - fadings not calibrated;
  - bright-core hits (32 / 50) not reviewed.
- **Screen:** 50 hits, 5 morph `pass`, 4 review entries (all fadings, all in the arc mask).
- **Vetting:**
  - All four lie 2.6–3.9″ from one bright point source, and three are on its diffraction spikes; PA_V3 differs by
    0.84°, so the
    spikes move about 0.05″. Tool: `scripts/elgordo_spike_check.py`.
  - h024 and h031: no forced-photometry change ≥ 2σ.
  - h049: on-spike controls 6 / 8 fade alike → spike shift.
  - h043: no significant F150W change. Its F444W −4.6 % is within the program-to-program scatter of 20 bright
    compact controls (control median z −9.8; tail fraction 0.45).
- **Result:** no real change without an ordinary or instrumental explanation; nothing to `/vet-candidate`; no SN or
  variable found. No exclusion or abundance limit.
- **Deviations** (all labelled in PLAN / RESULT):
  - the spike source fix;
  - the control box widened from 16″ to 30″ (17 < 20 controls);
  - the addendum wording corrected before measuring.
- **Flagged, not re-run:** the MACS1423 census (#282) and #279 used the pre-fix `fit_models`. TASKS says not to
  cite them until they are re-run.
## 2026-10-16: `arcs` PA convention: the deblend knot axis `pa_sep` and the catalogue `pa_obs` agree (median 3.0°)
- **Why:** the MACS1423 and RXJ2129 handoffs (below) list the `pa_sep` vs `pa_obs` convention as unchecked. Plan
  frozen first (`results/arcs/pa-convention/PLAN.md`); `scripts/arcs_pa_convention.py` reads the stored S1063 and
  MACS1423 census outputs, no download.
- **Result (derived; `RESULT.md`, `summary.json`):** 18 gate-passed `two-knot` rows: median offset 3.0° (identity)
  vs 35.2° (mirrored) and 87.0° (rotated 90°) → **frozen outcome `convention agrees`**. Secondary set (66 rows,
  ellipticity ≥ 0.2): 2.5° / 46.1° / 87.5°. The S1063 1316 (78.5°) and 1557 (75.2°) knot-axis readings stand.
- **Post hoc:** MACS1423 2197 (screen anti) has its knot pair 6.6° from the prediction but its catalogue axis
  74.7° from it, so its anti class does not come from the fitted pair; #282's census null does not use `pa_sep`.
- **Handoff:** the convention item is closed for S1063 and MACS1423. Not checked on the five fields of #279 (that
  branch's census JSON). The z_phot > 15 grid-edge item (RXJ2129 entry) is still open and needs its own plan.

## 2026-10-16: `arcs` on Abell S1063: 4 background anti vs 0.52 expected; 2004 at κ ≈ 1, 2184 compatible with the tested noise model; 1316, 1557 pending vetting; public spectra cannot settle them (#273)
- **New field for the arc-orientation test** (D-024 defaults; plan committed first): CATS v4.1 maps, GLIMPSE F200W,
  DJA v7.5 photo-z. Strong-shear background: 20 / 24 aligned, **4 anti vs 0.524 (P 0.0021)**.
- **Ordinary tests (new scripts):**
  - `arcs_bands.py` gains S1063: F444W controls 18 / 18; the anti have no F444W match;
  - `arcs_faint.py` (post hoc): all anti at m ≥ 27.4, but other fields' faint ends do not show it;
  - `arcs_bgsub_shapes.py`: ICL-subtracted moments, controls 18 / 20, 3 of 4 stay anti;
  - `arcs_second_models.py`: GLAFIC / Sharon / Keeton v4 κ, γ by byte range, controls 20 / 20, all 4 stay anti.
- **Per item:**
  - 2004 explained (κ ≈ 1);
  - 1557 likely an intrinsic shape or blend in weak shear;
  - **1316 and 2184 open, unvetted, not candidates**.
- **Deviation:** no known-image check (map model; as MACS1149).
- **Data:** the 75 MB DJA photo-z tarball was streamed to the scratchpad; the F200W catalogue sits in the gitignored
  cache; the team maps were read by byte range.
- **Vetting, second run (`vetting/`, plan first):** `scripts/arcs_vet_noise.py` injects aligned Gaussians with each
  row's flux, size and ellipticity at 60 empty positions per row and re-measures them. The gate passed (min 0.9667).
  1316 and 2184 both have **P_anti 0 / 60** → not noise-explicable. Sample: 3 observed vs 0.017 expected.
  - Post hoc: at injected e = 0.2 the targets keep P_anti ≤ 2 / 60 (their own gate 60 / 60), but the all-row gate
    fails on large galaxies. At e = 0 the test cannot exclude a near-round source.
  - Independent HST axes (HFF-DeepSpace SExtractor) are also anti for both, at S/N ≈ 1.
  - Cutouts looked at: 1316 shows two knots along its axis (blend untested); 2184 is a blob 2″ from the giant arc.
  - All-VizieR query within 1.5″: no spectroscopic redshift (not stored).
- **Vetting, step 2 (`vetting/PLAN_deblend.md` + addendum first; `scripts/arcs_vet_deblend.py`):**
  - 1316 (frozen): two knots 0.208″ apart (ΔBIC 134, gate 20 / 20), the separation is the anti axis (78.5°), knot
    colours consistent (p 0.146) → clump pair, stays open. Post hoc: at its fitted flux (⅓ of the catalogue flux)
    the gate gives 0.75 < 0.8, so the two-knot call is not validated there.
  - 2184 (frozen): arc-fragment test invalid (control fails without an error floor). Post hoc 10 % floor: control
    passes, 2184's SED differs from the giant arc (p ≤ 9 × 10⁻⁶).
  - **Failed assumption found:** step 1 injected catalogue fluxes, but the 0.3″ cores are −0.47 to +3.09 mag
    fainter (1316 +1.24, 2184 +2.49). Re-run with each row's core flux (`--flux core`): gate passes; **1316 0 / 60,
    2184 7 / 60 (noise-explicable)**; sample 3 vs 0.333 expected, P 0.00039 (excess remains; 1316 and 1557 at 0).
- **Handoff:** 2184 explained by measurement noise at core flux (post hoc), not a candidate. 1316 open, unvetted,
  not a candidate. Next: 1557 through the same tests; a redshift for 1316 (MUSE or a second photo-z code); use
  `--flux core` on the other `arcs` fields (`vetting/RESULT.md` "Next").
- **Vetting, step 3 (`vetting/PLAN_1557.md` first; `arcs_vet_deblend.py --target 1557`):** 1557 is `two-knot`
  (ΔBIC 227, 0.147″), and this time the gate ran at the source's own flux and separation and passed (0.85 / 0.95);
  separation 75° from the prediction, knot colours consistent (p 0.019) → **clump pair, stays open**.
  - Post hoc (`scripts/arcs_vet_clumps.py`): a PSF-sized third point 0.29″ away on the same axis, inside the same
    catalogue segment, same colour (p ≥ 0.0069).
  - Post hoc model: a randomly oriented point pair is anti with probability 0.144 / 0.137 at g 0.43 / 0.45 and
    0.084–0.113 at g 1.6–1.9; the shape-noise null behind P 0.0021 / 0.00039 does not model clump pairs.
- **Handoff:** both remaining anti rows are clump pairs; the excess is measured against a null without clumps.
  Next: a clump-pair census of the 24 strong-shear rows with the same frozen rule (plan first; exploratory), then
  the expected anti count with clump rows. MUSE redshifts need a download decision (owner).
- **Vetting, step 4 (`vetting/PLAN_census.md` first; `scripts/arcs_vet_census.py`; exploratory):** the step-3
  deblend with an own-flux gate on all 24 strong-shear rows: gate passed 21, `two-knot` 5 (1061, 1316, 1452, 1557,
  2449; f̂ 0.238). Two-knot rows get the clump-orientation anti probability: expected anti 0.52 → 1.24, P(≥ 4)
  0.00082 → **0.026** → frozen outcome **excess remains under the clump-aware null**. P ≥ 0.05 only for k = 3 / 2 or
  f at its 95 % upper end (0.093).
- **Handoff:** 1316 and 1557 open, unvetted, not candidates. Next: z_spec from the ESO MUSE cube (needs the owner's
  download decision), else a second photo-z code or the same census on another field's anti rows.
- **Owner-corrected summary (2026-10-16; records kept):**
  - P 0.026 is exploratory and depends on a simplified ordinary model: the clump fraction is uncertain (0.093 at
    f = 0.472), and 2004's κ ≈ 1 degeneracy is not modelled.
  - Ordinary morphology, noise and model error are not excluded in general; the models' errors are not fully
    independent.
  - 2184 is compatible with the tested noise model.
- **Step 5, spectroscopy feasibility (`vetting/PLAN_spectra.md` first; `scripts/arcs_vet_spectra.py`):**
  - No catalogued redshift within 2.4″ (all-VizieR 6″ cone). Richard+2021 does not cover AS1063.
  - ESO MUSE cubes are 5.7–5.8 GB each, but SODA cutouts are 12–17 MB.
  - Lines, after the review fixes (below):
    - 1557: no usable position reaches the threshold.
    - 1316: one does, at 6000.3 Å, and it is classified as not significant (below).
    - Positions that are not usable: 37–42 % of the range and 59 % of 1557's [O II] window.
  - The median 5σ limits (1.34 × 10⁻¹⁸ for [O II] of 1557, 4.4 h; 3.37 × 10⁻¹⁸ for C III] of 1316, 2.2 h, 1.45″)
    are sensitivity estimates. They depend on the line width, channel noise, continuum and aperture assumptions,
    and they do not exclude a line anywhere in the window. They correspond to rest EW of order 240 Å / 250 Å.
  - At the present resolution (≥ 0.62″) and S/N, MUSE cannot reliably separate the 0.15–0.21″ clumps.
  - The result is undetermined. Implementation deviation: the sky flag was replaced before use.
- **Conclusion:** the spectra examined did not confirm a redshift; ordinary explanations are not excluded; the
  pause is a resource decision.
- **Owner decision 2026-10-16:**
  - Wait for 116.28VK.001.
  - 1316 and 1557 stay as pending-vetting flags with insufficient data.
  - No stacking, no whole-cube download, no new spectroscopy tool unit.
  - Reassess only for a reliable target redshift, a usable 1-D spectrum or materially better data; the arrival of
    data is not an automatic download or re-run.
  - Other authorized searches continue.
- **Final code review → fixes and a limited re-run (owner-approved; `vetting/PLAN_fixes.md` first):**
  - Fixes: usable-channel noise rescale; per-channel validity; deblend flux start and bound from the positive
    pixels; census null drops the rows the count drops. Each has a unit test (`tests/test_arcs_vet_fixes.py`).
  - Re-run used the same two MUSE cutouts and the same JWST stamps. Before-fix outputs are kept as
    `*_before_fixes.*`, with the diff in `fixes_diff.json`.
  - **Changed:**
    - 1316 rescale 2.44 → 1.73; C III] limit 4.75 → 3.37 × 10⁻¹⁸; C IV 5.51 → 3.89 × 10⁻¹⁸;
    - the 6000.3 Å peak S/N 4.77 → 6.72;
    - census P(≥ 2) without 2184 0.359 → 0.354;
    - injection-gate recalls (1316 post hoc 0.75 → 0.95 two-knot, still failing on the single side at 0.75; 1557
      0.90 / 0.95 → 1.00).
  - **Unchanged:** every class, the 21 / 24 census gates (the flux bound did not cause the 1404 / 2004 / 2184 gate
    failures), P(≥ 4) 0.026, P(≥ 3) 0.116, and `clumps_1557.json`.
  - **6000 Å peak (pre-registered checks):** valid and above threshold, but:
    - the narrow-band peak is 2.79″ away (sky PA 224°, at the edge of the valid map);
    - the 1.5″ W offset aperture is stronger (S/N 22.9);
    - the planned chance metric gives 1.83 per window search.
    - → **not significant under the calibration** (frozen class), i.e. not a line of 1316; no redshift.
    - Deviation, reported only: excluding the feature itself, no chance peak this strong elsewhere in the window,
      so the feature is real emission centred off 1316.
    - Decision to wait: unchanged.
## 2026-10-16: `arcs` on RXJ2129 with the DJA v7.0 catalogue: no result again (controls 64 % and 67 %)
- **Why:** the RXJ2129 handoff (below): 12 strong-shear rows were too few. Plan frozen first
  (`results/arcs/rxj2129-dja/PLAN.md`, `ea2feb6`); `scripts/arcs_rxj2129_dja.py` converts the DJA SEP catalogue
  (`ir` detection image) to the pipeline columns and reuses `arcs_rxj2129.main` and its outcome rule unchanged.
  Same field and models, so a re-test, not an independent field.
- **Tool check (derived; `convention.json`):** PA = θ_image + 90° matches the pipeline `sky_orientation` on 68
  matched rows (median |ΔPA| 2.8°, gate ≤ 10°). The contact-sheet axes are now drawn through the WCS.
- **Result (`RESULT.md`, `summary.json`):** 89 strong-shear background rows; aligned 57 (Lenstool) and 60 (GLAFIC),
  both under 70 %: **frozen outcome "no result"**. Post hoc (`posthoc.json`): 17 rows have z_phot > 15 (eazy grid
  edge) and hold 5 of the 10 anti; the F200W sheet (looked at, stored) shows most anti rows without a clear F200W
  counterpart or on linear structures. Deviation: `mag` is `mag_auto` (the planned `flux_iso` formula was in
  detection-image units); it enters no class.
- **Found for later:** `arcs` treats z_phot at the eazy grid edge (> 15) as background. No stored anti row on main
  is affected (all `arcs.ecsv` checked, not stored); a z_phot validity cut would need its own plan.
- **Handoff:** both public RXJ2129 catalogues used so far give "no result". The S1063 tally
  stays: MACS0417 no, MACS0647 no result, MACS1423 explained, RXJ2129 no result. Next field options need public
  maps and F200W (MRG0138, RXCJ0600 lack public maps; SOURCES.md).

## 2026-10-16: `arcs` on RXJ2129: no result (both JWST-updated models fail the control, 8 / 12 aligned)
- **Why:** the MACS1423 handoff (below): RXJ2129 as the next independent field for the S1063 anti-excess
  replication question. Plan frozen first (`results/arcs/rxj2129/PLAN.md`, `bbe0306`); script
  `scripts/arcs_rxj2129.py` reads the Lenstool (primary) and GLAFIC (second) maps out of the 4.30 GB Zenodo zip by
  HTTP range requests and reuses `arcs_macs0647.outcome` unchanged.
- **Result (derived; `results/arcs/rxj2129/RESULT.md`, `summary.json`):** both maps pass `validate` (median κ
  relative difference 2.9 × 10⁻⁵ and 1.2 × 10⁻⁵); 12 strong-shear background rows, 8 aligned (67 %) under both
  models, so both controls fail: **frozen outcome "no result"**. Post hoc: Lenstool anti 2 vs 0.22 shape-noise
  expectation, GLAFIC anti 1 vs 0.176; these carry no weight under the rule. Cutouts looked at (not stored).
- **Failed approach:** the MAST download API returned "internal error" for the program 2767 catalogue; the public
  S3 bucket served it. A range read of the Zenodo zip dropped once on a TLS EOF; the script now retries each block.
- **Handoff:** the S1063 replication tally is MACS0417 no, MACS0647 no result, MACS1423 explained, RXJ2129 no
  result (12 strong-shear rows on a 2,061 s catalogue). A deeper RXJ2129 catalogue (DJA `rxj2129-grizli-v7.0-ir.cat`,
  all NIRCam bands) would add rows but changes the shape source, so it needs its own plan. PA conventions of the
  deblend `pa_sep` vs catalogue `pa_obs` stay unchecked.

## 2026-10-16: `arcs` on MACS1423, vetting step 2: the clump-pair census explains the anti count (P 0.259, exploratory)
- **Why:** the step-1 handoff (below): the S1063 clump-pair census on the MACS1423 rows, anti classes at each
  row's photo-z. Plan frozen first: `results/arcs/macs1423-canucs/vetting/PLAN_census.md`; script
  `scripts/arcs_macs1423_census.py` reuses #273's census functions unchanged (`43fbbbd`).
- **Result (derived; `vetting/RESULT.md` step 2, `census.json`):** 45 / 55 rows pass the deblend gate; 13 two-knot
  (f̂ 0.289). 2197 and 2322 are two-knot, 4234 single, 1818 gate-failed. Expected anti 2.608 (clump-aware) vs 0.946
  (shape noise); **P(≥ 4) 0.259 → "anti excess explained by clump orientation (exploratory)"**; k = 3: P 0.494.
  Cutouts looked at (`census_anti.png`).
- **Reading:** two independent fields with JWST-era models (MACS0417, MACS1423) do not replicate the S1063 anti
  excess under their nulls; S1063 stays a single-field exploratory result (P 0.031 clump-aware, #273). Not a candidate.
- **Handoff:** the replication question needs a field where S1063-like anti rows survive the clump-aware null; the
  next independent field is RXJ2129 (six JWST-updated models, Dhanasingham et al. 2026; plan first). PA conventions
  of the deblend `pa_sep` vs catalogue `pa_obs` are unchecked.

## 2026-10-16: `arcs` on MACS1423 (CANUCS DR1 maps): the frozen rule says the S1063 anti excess replicates at screen level; after ICL subtraction it is marginal (post hoc P 0.071)
- **Why:** the MACS0647 handoff said MACS1423 lacked a public model; CANUCS DR1 has a JWST-era Lenstool model for it
  (Desprez, preliminary v1, 5 spectroscopic systems). New entry `macs1423-canucs-maps` (`validate`: κ reproduced to
  3.0 × 10⁻⁵). Plan and `outcome()` committed before the run: `results/arcs/macs1423-canucs/PLAN.md`,
  `scripts/arcs_macs1423.py`.
- **Screen (derived; `results/arcs/macs1423-canucs/RESULT.md`, `summary.json`):** control 44 / 55 aligned (80 %);
  **anti 4 vs 0.949, P 0.016 → frozen outcome "replicates"** (screen level). Cutouts looked at: 1818 on the BCG's
  ICL, 2197 one knot of a clumpy chain, 2322 a tiny faint blob, 4234 a bright galaxy at the g_min threshold.
- **Vetting step 1 (plan first; #273's scripts at `43fbbbd`, unchanged; `vetting/RESULT.md`):** ICL subtraction
  keeps 2197, 2322, 4234 anti and moves 1818 to `mixed` (controls 40 / 43); measurement-noise injections (gate
  100 %) give P_anti ≤ 0.017 per row → frozen outcome "excess remains" vs measurement noise. **Deviation, post hoc
  (`vetting/posthoc.json`):** the scripts classify at z_s = 2 only, so the frozen count includes 1758 (aligned at its
  z_phot 6.38); the 3 screen anti that stay anti vs the shape-noise expectation give P(≥ 3) 0.071.
- **Reading:** not a candidate. Two independent fields with JWST-era models now give one non-replication (MACS0417,
  P 0.195) and one screen-level replication that becomes marginal after ICL subtraction. The `arcs` shape-noise null
  covers intrinsic shapes only; clumpy fragments (2197 by eye) are the next ordinary test.
- **Found for later (`SOURCES.md`):** RXJ2129 has six JWST-updated models (Dhanasingham et al. 2026, Zenodo 4.3 GB
  zip, range-readable); MRG0138 and RXCJ0600 have JWST-era models without public maps.
- **Handoff:** done in the entry above (clump-pair census); then RXJ2129 as a third independent field (plan first; check the
  GLAFIC ZS = 3 header).

## 2026-10-16: `arcs` on MACS0647 as a second independent test of the S1063 anti excess: no result (both controls fail)
- **Why:** after MACS0417 (#280), the next field with JWST F200W, a DJA photo-z catalogue and public lens models.
  Plan and outcome rule frozen and committed before the run: `results/arcs/macs0647-clash/PLAN.md`,
  `scripts/arcs_macs0647.py` (`outcome()`, unit-tested). New model entries `macs0647-clash-nfw` (primary) and
  `macs0647-clash-ltm` (CLASH Zitrin v2 maps; `validate` reproduces κ and z = 2 μ to ≤ 3 × 10⁻⁵).
- **Result (derived; `results/arcs/macs0647-clash/RESULT.md`, `summary.json`):** JWST 1433 epoch 1 F200W: PIEMD + eNFW
  10 / 15 background rows aligned (67 %), LTM-Gauss 7 / 14 (50 %): **both controls fail → no result**. Post hoc:
  the primary failure is one row from the threshold and passes under one of four 0.3″ frame shifts (10 / 14; anti
  stays 2); epoch 2 (deviation, not deciding) gives 13 / 18 aligned, anti 2 vs 0.425 (P 0.068).
- **Reading:** all three HST-era map models tried in the two replication fields failed the alignment control (one
  marginally); the one JWST-era model passed. The S1063 excess stays a single-field exploratory result (one non-replication, one
  field that cannot test it).
- **Data:** CLASH maps (70.6 MB each) cached and deleted after the run; the DJA photo-z tarball was streamed.
- **Handoff:** a further replication needs a JWST-era model on a field with F200W and photo-z. The remaining DJA
  photo-z cluster fields (RXJ2129, MRG0138, RXCJ0600, rj0018p1626) have no JWST-era model pinned here (not searched
  further in this run); MACS1423 (CANUCS photo-z) still lacks a public model.

## 2026-10-16: `arcs` on MACS0417 as an independent test of the S1063 anti excess: does not replicate
- **Why:** the S1063 anti excess (#273, #279) was found post hoc in one field. Plan before any run:
  `results/arcs/macs0417-relics/PLAN.md` (D-024 rules unchanged, control ≥ 70 % aligned, replication if P < 0.05).
- **New model entries:** `macs0417-relics-lenstool` (RELICS v1 maps; reproduce the published κ and z = 2 μ to
  ~2 × 10⁻⁵) and `macs0417-canucs-maps` (CANUCS DR1 preliminary maps; κ from the deflections matches to 3 × 10⁻⁵).
- **Result (derived; `results/arcs/macs0417-relics/RESULT.md`):**
  - run 1, RELICS v1 (HST-era): control fails, 39 / 61 aligned (64 %); its 8 anti are not interpreted;
  - run 2, CANUCS v1 (chosen after run 1 failed, addendum A written before running it): 37 / 49 aligned,
    **anti 3 vs 1.52 expected, P 0.195 → does not replicate**. Cutouts looked at: a low-surface-brightness smudge,
    a compact knot, a source on a detector stripe.
- **Reading:** the S1063 excess stays a single-field, exploratory result. An HST-era model can fail the alignment
  control in a JWST-depth catalogue (most run 1 misalignment at g_min ≥ 1).
- **Data:** three CANUCS maps (85–89 MB each) were cached in the session and deleted after the run.
- **Handoff:** no further `arcs` field is queued for the S1063 question; MACS1423 has CANUCS photo-z but no public
  model was found in this run.

## 2026-10-16: No-injection control at vetoed SMACS injections: mixed; vetoed and passed positions both sit on pre-existing negative structure
- **Why:** handoff of the review-list unit: off the arc mask, morphology vetoes most recovered SMACS injections; was
  that pre-existing residuals (H1) or the injection itself (H2)? Plan before measuring:
  `results/transients/veto_control/PLAN.md`; code `scripts/arc_veto_control.py` (unit-tested). No rule changed.
- **Result (derived; `results/transients/veto_control/RESULT.md`, `compare.json`):** off the mask (where all
  vetoes are), on the re-run no-injection map, 19 / 19 rosette or dipole positions and 15 / 19 passed ones carry a
  negative region (|S| ≥ 3) touching 0.8″. Frozen rule: **mixed** (neither H1 nor H2 as defined). Post hoc: for 11
  of the 19 vetoed (8 of 19 passed) that region is under the position itself. By eye (sheet looked at): vetoed stars
  mostly sit on compact ± core residuals, spike stripes or a bright neighbour's residual field.
- **Deviations:** the whole-box row (H1) is diluted by in-mask passes and not used; the re-run `hits.ecsv` matches
  the stored one to 5 × 10⁻¹² in σ (not stored), not byte for byte.
- **Handoff:** next, if wanted: location-aware lobes (outside the region under the star) on injected and
  no-injection maps. Any morphology-rule change needs the owner's go.

## 2026-10-16: Arc difference imaging with position as a label: known cases kept, review burden small, injection efficiency fails (0.34 < 0.50)
- **Why:** owner decision of 2026-10-16: position becomes a label and review priority, not a veto, and every other
  criterion stays. Plan and acceptance rules were written before any new count:
  `results/transients/review_list/PLAN.md`. Code: `scripts/arc_chain.py review`, `inject-review`.
- **Result (derived; `results/transients/review_list/RESULT.md`):**
  - **R1 (development controls, not independent):** 6 / 6 are on the review list (Atalanta A/B; Yan+2023 D21-W1,
    W2, S1, S2).
  - **B1 (burden):** 4, 1, 0 and 1 added entries on SMACS ERO → VENUS, Sunrise o052, Sunrise o120 and the MACS0416
    Warhol box. Of the 6 added entries, 3 are known real transients, 1 was vetted with no ≥ 2σ change, and 2 are
    artifacts by eye (h039 SMACS, h005 Sunrise o052; not vetted by photometry).
  - **R2 fails:** SMACS whole-box efficiency (hit, unflagged, morph pass) is 0.34 at 27.5 AB (0.50 in the mask, 0.25
    outside it). Off the arcs, morphology vetoes most of the recovered stars (8 rosette, 5 dipole, 3 extended of 29).
- **Decision per plan:** no new box was searched. The review lists and the failing case are saved. The two by-eye
  entries remain open.
- **Handoff:** a cheap check would show whether the off-arc vetoes come from neighbouring sources (the MACS0416 pattern of
  2026-10-15): a no-injection control at the vetoed stars' positions. Any change to the morphology rules needs the
  owner's go (owner, 2026-10-16: no further veto tuning without a decision).

## 2026-10-16: Full-chain known-case check of arc difference imaging: MACS0416 4 / 4 pass, SN Atalanta fails at the position stage
- **Why:** the owner asked that known real signals pass the whole chain (detection, morphology, position, multi-band)
  before more arcs are searched. Plan before any run: `results/transients/arc_chain/PLAN.md`. Post-hoc arc mask
  (image-based, frozen parameters) in new `scripts/arc_chain.py`, unit-tested.
- **Result (derived; `results/transients/arc_chain/RESULT.md`):**
  - **MACS0416 Yan+2023.** The 4 transients detected in o211 → o212 pass every stage; 3 of them are non-trivial
    tests of the mask, since D21-W1 is 0.19″ from a catalogued image. D21-S2 sits 1.43″ from its
    catalogued image, so a 1″ rule would have dropped it.
  - **SMACS Atalanta A/B.** Detected and morph `pass`, but outside the mask (0.76–0.86″) as well as outside 1″. The
    host light at those positions is below the mask threshold (0.30–1.78 σ).
  - **SMACS full-chain injection efficiency inside the mask:** 0.50 at 27.5 AB and 0.31 at 28.0 AB (3 seeds). The
    pre-morphology per-band figures were 0.82 / 0.76 (SMACS) and 1.00 / 0.97 (Sunrise).
- **Failed approach / plan error:** Earendel was listed as an independent mask check, but it is catalogued image 1.1e
  and seeds its own component.
- **Limits:**
  - the mask misses transients on host light outside the dilated, seeded components, by construction;
  - SMACS in-mask samples are small (26 and 29 stars);
  - there is no masked full-chain efficiency for MACS0416 or Sunrise.
- **Handoff:** owner decision before the next new arc data. The frozen rule counts the chain as validated on its
  independent cases, but a known real transient (Atalanta) still fails. Not started; it needs the owner's go:
  - report every morph `pass` regardless of position, with position as a label rather than a veto;
  - a host-segmentation position rule (tool work, so the owner decides).

## 2026-10-16: SMACS 0723 arc difference imaging: Atalanta recovered, no candidate; h061 undetermined (near-star spikes); the 2024 → 2026-07 pair fails its control
- **Why:** the Sunrise handoff names SMACS next; a catalog grid misses transients blended into arcs (D-027 addendum
  2026-10-15j). Plan before the run: `results/transients/smacs_arcdiff/PLAN.md` (24″ box on the system 5 arc, SN
  Atalanta as the known-case control, frozen "on an image" rule = `classify` in new `scripts/smacs_arcdiff_images.py`).
- **ERO o001 → VENUS o057 (F150W + F444W, derived):** Atalanta A, B recovered (S +49.5 to +58.2, both morph passes);
  injections 0.82 / 0.76 at 27.5 AB per band. 73 hits, 0 `vet` under the frozen 1″ rule, **but that radius
  fails the known case** (Atalanta is 1.59–1.69″ from its catalogued images; found in review). Post hoc 2.0″:
  h063 (pass, 1.37″ from 6.1) shows no ≥ 2σ change in multi-epoch forced photometry (F150W −16.0 ± 10.4 nJy).
  **h061 open:** on image 12.1 (box-edge veto; S depends on the box), flat 2022 → 2026-06, but F444W +35.6 ± 3.4 nJy
  in the 2024 GO-4043 mosaic, where a stripe crosses it by eye (spike geometry not computed). No candidate so far.
- **GO-4043 → DD 12774 (F115W + F444W; added, not in the PLAN):** Atalanta not a hit (F115W S ≈ 2.4; forced F115W
  13.7 ± 16.7 nJy; band depth and SN colour not separated), injections 0.26 / 0.34 at 27.5 AB: no null from this pair. Its one pass, h000, is a point
  source only in the 2024 mosaics (F115W 2688 nJy in 2024 vs 49 ± 17 nJy in 2026-07; F090W–F444W 7–49 nJy in every other epoch; `hit_lightcurve.py`), 0.081″
  from TNS AT 2024jfz: known transient, not a candidate.
- **h061 follow-up (PLAN addendum 1, frozen before the run):** new `scripts/stripe_phot.py` subtracts the 2024
  stripe's across-profile (PA 65.5°) in every F444W mosaic; 2024 stays 2.39σ above the other epochs →
  `undetermined` under the frozen rule. Its injection check turned out tautological (the source never enters
  its own profile; stated in RESULT). Post hoc: h061 is 4.70″ from a bright star whose spikes rotate per visit;
  4 / 17 compact controls within 7″ of it reach |z| ≥ h061's 14.2 in the same pair (0 / 42 beyond, but they fade
  where h061 brightens), and 7 / 34 positions along the same 2024 spike reach |Δ| ≥ its 35.6 nJy. **Not
  distinguishable from bright-star spike systematics with these tests; not a candidate.** Reopen on a PSF model
  of the star (spikes) or a counter-image change (RESULT.md).
- **Run time:** each 24″ two-band `arc_diff` run took ≈ 1 min at `--procs 4` with `OMP_NUM_THREADS=1` (not stored).
- **Handoff:** more SMACS boxes (systems 1–4, 6–18); a counter-image check for AT 2024jfz (1.8″ from image 6.1);
  El Gordo arcs (Sunrise o010 → o120 done, #274).
## 2026-10-16: Sunrise arc difference imaging, o010 → o120 (F277W + F356W): 0 candidates on the arc
- **Rules:** fixed before the run in `results/transients/sunrise_arcdiff/PLAN.md` addenda 1–2:
  - settings as for o052;
  - a band counts if its recovery at 27.5 AB is ≥ 0.8;
  - "on the arc" means within 1.0″ of a Scofield+2025 image;
  - a candidate must pass morphology, be on the arc and appear in ≥ 2 bands.
- **Deviation:** the four-band run gave no output in 90 min, so long-wavelength only. Short-wavelength-only changes
  are not covered.
- **Result (derived; `o010_o120/`):**
  - injection recovery at 27.5 AB is 1.00 (F277W) and 0.97 (F356W);
  - 12 hits: 10 rosette, 1 dipole, 1 edge, 0 morphology passes; all ≥ 3.87″ from every catalogued arc image;
  - looked at on the sheet: point-spread-function rotation residuals, a disk and the map edge; the Earendel control
    is clean.
  - **0 candidates.**
- **Also:** the TASKS line for the 7 undecided Abell 2744 orphan pairs now states the owner's scope. A redshift
  difference argues against lensing, equal redshifts would not confirm it, and not pursuing them now is a resource
  decision.

## 2026-10-16: Sunrise arc difference imaging, o010 → VENUS o052: 0 candidates on the arc; o120 pair still to run
- **Why:** the Sunrise catalog grids are null, but a catalog grid misses transients blended into an arc (PLATE-23a,
  D-027 addendum 2026-10-15j); `arc_diff.py` had never run on the Sunrise arc. Plan before the run:
  `results/transients/sunrise_arcdiff/PLAN.md`.
- **Run (derived):** 24″ box on Earendel (holds the arc's catalogued images), `--kernel-px 4 --kernel-stamps --inject 100`, F150W + F444W
  (jwst 2.0.1 vs 3.0.0); injections on bright pixels recovered 0.86 / 0.83 at 27.5 AB, 0.14 / 0.28 at 28.0.
- **Result:** 9 hits; `arc_diff_morph.py` passes 4. On the sheet (`o010_o052/sheet.jpg`, new
  `scripts/arc_diff_sheet.py`) all 4 are off the arc: bright-core residuals (h001, h002), an edge-on disk (h007), a
  compact-core dipole at the box border (h005); all hits ≥ 3.98″ from every Scofield+2025 multiple image (post hoc rule). **0 candidates on the
  arc's images.** Earendel is not a hit (S −5.2 F150W, −1.5 F444W; stamp clean).
- **Failed approach (times not stored):** two `--kernel-stamps` runs without `OMP_NUM_THREADS=1` gave no output in
  ≈ 20 min (known, CHANGELOG 2026-10-15). The o120 run at `--procs 2` did not finish its first two SW bands in
  ≈ 45 min and was stopped.
- **Deviations:** `--inject 100` instead of the PLAN's `--inject-grid`; o120 postponed; "on the arc" defined post hoc
  (RESULT.md).
- **Handoff:** run o010 → o120 alone with `--procs 4` (the next highest-value step on this arc); then the SMACS 0723
  and El Gordo arcs.

## 2026-10-16: Abell 2744 orphans: per-pair vetting of the 11 escalated pairs: 4 ordinary, 7 insufficient data, 0 anomalies
- **Plan and timing:** `results/orphan_pairs/abell2744/vetting/PLAN.md` was written before measuring. Its addendum
  (E2 controls, `classify` / `e2_valid` as unit-tested functions) came after the first orphan measurements, so it is
  post hoc for the orphans.
- **Tool:** `scripts/orphan_vet.py` reads UNCOVER boxes by S3 byte range, and `--classify` runs and stores every
  VizieR query. Nothing large was downloaded.
- **Fixes after review:**
  - the "LW stack" had summed a misaligned band; stacking now requires member pixel positions to agree to ≤ 0.1 px,
    which leaves F277W only, plus F444W for #7 and #10;
  - E4 is now a stored flux + 3σ upper bound;
  - #3's E1 strip excess fell from 3.2σ to 1.2σ.
- **E2 validation** on 11 CATS pairs ≤ 4″ apart: in every control one member reaches 3σ in at most 1 band; 1 of 7 measured
  |sb_z| > 3 (rule ≤ 10 %). So E2 is not valid and supports or excludes nothing.
- **E5/E6:**
  - no public spectroscopic redshift lies within 1.5″ of any member (Price+25, Wang+24 `zsp`);
  - same-z neighbours number 0–3 per pair.
- **#7, member a:** UNCOVER DR3 detects it 0.006″ from our aperture. It is below 3σ only in our conservative
  photometry.
- **Classes** (`classification.json`):
  - ordinary: #1 (visual artifact flag), and #5, #6 and #14 (E1 strip excess 3.3–4.8σ);
  - insufficient: #2, #3, #7, #9, #10, #11 and #13;
  - anomaly: none.
- **Separate from the #269 population null and its ~0.3–1 % per-deflector efficiency:** nothing here is a candidate.
- **Handoff:** the deciding data are spectra of both members. The no-new-data route needs a validated PSF-matched
  morphology tool. Neither is in this unit.

## 2026-10-16: Process: "Before writing results" rules in scripts/CLAUDE.md
- **Why:** reviews of #262, #266 and #269 each found 8–12 problems, mostly of three kinds:
  - frozen rules applied loosely;
  - numbers not traceable to a committed file;
  - overclaiming wording.
- **Change:** `scripts/CLAUDE.md` gains a short section with five rules:
  - computable frozen rules, with functions and tests where possible;
  - every number traceable to a committed JSON/ECSV, or labelled "not stored", and ranges cover only the items
    classified;
  - narrowest wording for null / not affected / consistent with chance / equals / low completeness and for any causal
    reading;
  - deviations labelled;
  - ruff format and fix before ruff check.
- `.claude/skills/research-cycle/SKILL.md` step 8 points to the new section.
- No tool, result or threshold changed.
- **Open:** the Abell 2744 orphan vetting (#271) is paused unfinished; its state is in
  `results/orphan_pairs/abell2744/vetting/STATUS.md` on its branch.

## 2026-10-16: Orphan-pair dark-deflector screen on Abell 2744: null under the frozen rule; 11 top orphans unvetted
- **Plan before the run** (`results/orphan_pairs/abell2744/PLAN.md`, committed first): the D-048/D-055 rules
  unchanged. An excess needs P < 0.05 under both null (e) and null (f). It is the first run of this statistic on this
  field, which is exploratory, not blind.
- **Tool change:** `orphan_pairs.py` gains `FIELDS["abell2744"]` (DJA v7.2-fix phot and zout, already pinned; the zout
  is a plain file). `inject_pairs.py --fields abell2744` is allowed.
- **Checks:**
  - alignment 0.02″ against the published images;
  - **no positive control** (0 `published` pairs at 0.3–3″);
  - injection efficiencies match the deep fields within noise (point 0.30 % / 0.95 % per deflector at θ_E 0.3″ /
    0.7″; W1 15 % at 1.5″), so completeness is low here as everywhere; the cluster magnification of the injected
    images is not modelled, and the pooled `limits.json` is not valid for a cluster field;
  - artifacts in the top 15: a mosaic edge (#12), a star spike (#1), one off-footprint (#8).
- **Result (derived):** 9,385 sources, 174 SED-matched close pairs, 51 orphans against 43.8 / 42.4 expected; P 0.16 /
  0.11. **Null under the frozen rule.** The excess over null (a) persists under the redshift nulls (c) and (d) (P 0.016 /
  0.021) and disappears under the conditioned nulls (e) and (f).
- **Top orphans:** 11 of the top 15 meet the frozen escalation rule. They are unvetted and not candidates until
  vetted; the count test cannot classify individual pairs.
  - The per-pair `/vet-candidate` step is postponed (a deviation from the plan).
  - A first sheet pass with the D-051 categories finds 5 companions or group members, 2 possible knots, 1
    artefact-affected pair and 3 featureless faint pairs.
- **Data:** the 233.5 MB phot file was downloaded to the cloud scratchpad and deleted after use.
- **Code fixes after review:**
  - `inject_pairs.py` refuses to pool cluster fields into `limits.json`;
  - `orphan_pairs.lens_check` no longer falls back to the placeholder DJA `mu` outside the CATS maps.
- **Handoff:** the open ordinary check is `/vet-candidate` on the 3 featureless pairs (#9, #11, #13) and the 2 possible
  knots (#2, #3). MACS0417 and MACS1423 cluster pointings need published image lists.

## 2026-10-16: Abell 2744 `arc_diff` stripe, spike and offset-core vetoes on the annulus cut: they remove what they target, but other artifacts pass; not adopted
- **Hypothesis before measuring (#265 handoff):** the 31 extra real passes of the v2 annulus cut are detector
  stripes, spike fragments, compact-core residuals offset from their peak and faint compact fades; vetoes aimed at
  the first three keep the injected point brightenings on smooth light.
- **Vetoes (written before the run; ASSUMPTIONs), applied only to hits the v2 cut calls extended:**
  - *line coherence* (stripes and spikes are lines in S): for position angles 0–175° in 5° steps, the signed median
    of S along each half-line from 0.25″ to 1.0″ from the hit; the hit is vetoed if, at some angle, **both**
    half-line medians have the hit's sign and the smaller is ≥ 0.25 |S_hit| and ≥ 1;
  - *offset core*: an epoch-1 local maximum 0.15–0.45″ from the hit with contrast ≥ 20 annulus σ sends the hit back
    to the peer cut (as a compact hit).
- **Decision rule (written before the run):** adopt if the 28-mag chain efficiency is ≥ 0.30 with 27.5 mag ≥ 0.60,
  real peer-stage passes on the 4 tiles are ≤ 1.5 × the old 31 (≤ 46), **and** the sheet of newly passing real
  hits shows at most 2 recognisable stripes, spikes or core residuals.
- **Tool (derived):** `scripts/arc_diff_extended.py` (`line_coherence`, `ext_pass_vetoed`, columns `line_frac`,
  `line_theta`, `core_off_contrast`); `arc_diff_inject_chain.py --extended --vetoes` (stages `peers_veto` /
  `confirm_veto`); `arc_diff_extended_eval.py` adds a `veto` column. Unit-tested on a synthetic stripe, an isolated
  point and an offset core.
- **Run:** the 8 injection maps and the 4 real-tile maps re-made with the committed `arc_diff.py` commands (12 runs,
  4 at a time with `OMP_NUM_THREADS=1`, ≈ 6 min; S3 byte ranges, maps in the scratchpad only). Stage counts up to
  `peers` reproduce #264 exactly.
- **Result (simulated stars in real maps; real passes `derived`):**

  | version | 27.5 mag | 28.0 mag | real peer-stage passes (o001-confirmed) |
  |---|---|---|---|
  | peer cut (#261) | 0.637 | 0.188 | 31 (16) |
  | annulus v2 (#265) | 0.675 | 0.375 | 62 (45) |
  | **annulus v2 + vetoes** | **0.650** | **0.338** | **46 (31)** |

  The efficiency and count parts of the rule pass (46 is exactly the 1.5 × limit).
- **Looked at both sheets** (`extended/sheet_vetoes_removed.png`, `sheet_vetoes_newly_passing.png`):
  - the 12 hits the vetoes remove from v2 are all artifacts: 4 detector stripes (h270, h299, h331, h412), a spike
    (h321), 2 bright-star spike fans (h612, h638) and 5 compact-core residuals (h312, h344, h562, h607, h657);
  - the 19 that still pass and the peer cut vetoed include **≥ 4 recognisable artifacts**: residuals in a bright
    star's halo (h292, h317), mosaic-edge structure (h303, h564; h559 next to a vertical stripe) and h183, the edge-on
    galaxy dipole whose offset-core contrast is just under 20. The rest are faint fades (S −5 to −7) on galaxy
    outskirts; no positive point brightening on an arc.
  - **Sheet rule fails → not adopted.** The vetoes stay in the code (off by default) for a later version.
- **Consequence:** #261's "0 lensing candidates" still holds only at ≳ 36 nJy (#264); no new candidate.
- **Handoff:** what is left are bright-residual neighbourhoods (|S| ≫ 5 within ~1.5″) and coverage edges. A third
  version would need both vetoes and a new sheet check; each round has cost ~0.5 h for ×1.8 at 28 mag, so moving
  on from this tiling (the 1176 registration is done, #266) is the better use of time unless a fainter limit
  becomes decisive.
## 2026-10-16: PEARLS 1176 o221 vs UNCOVER registration on Abell 2744: a 0.027″ shift, no 0.3″ offset; h182 stays closed
- **Plan before measuring** (`results/transients/a2744_grid_alt/reg1176/PLAN.md`, committed first):
  - controls from the official level-3 catalogues (compact, isolated, AB 21–25.5), with a held-out half B by crc32 parity;
  - rule 1: M1 (affine) instead of M0 (shift) if it cuts B's RMS by ≥ 20 % with slopes ≥ 3σ;
  - rule 2: a band passes if |B median| ≤ 0.010″ and B RMS ≤ 1.5 × the ALT-vs-UNCOVER baseline.
- **Result (derived; `reg_offsets.py`, `reg1176/RESULT.md`).** In all 7 shared bands, 1176 sits +0.018…+0.027″ in RA
  cosδ and +0.006…−0.014″ in Dec from UNCOVER at h182.
  - Rule 1 picks M1 in F150W, F200W and F277W (slopes ≤ 0.007″/arcmin) and M0 elsewhere.
  - Every band passes rule 2. F115W is weak: 7 calibration sources and 10 held-out.
  - No control is off by > 0.11″. Near h182 the local shifts are within about 0.01″ of the global one.
  - The ALT baseline has its own +0.013″ RA offset; the plan had assumed none.
- **Image check** (`reg_imgcheck.py`, F356W, 80″ box): the FITS-header WCS used by `hit_lightcurve.py` agrees with
  the catalogue centroids to ≤ 0.008″ (median). **#263's by-eye "≈ 0.3″ offset" on c63 is not a 1176 registration
  offset.** The cause is unexplained; that control came from an unsaved first draw.
- **Re-measurement** (`hit_lightcurve.py --offsets`, `--only-program`):
  - **Reproduction:** the uncorrected rerun reproduces #263's 632 1176 fluxes exactly.
  - **Size:** the shift (≈ 0.4 px LW, ≈ 0.9 px SW) moves on-arc fluxes by a median of 1–3σ.
  - **h182:** z_cal in the 1176 pairs is now −0.5…+1.2, so it **stays closed**. This is exploratory: h182 and its
    controls had been seen before.
  - **Control spread:** it falls 27–36 % in F200W and F410M but rises 15 % in F356W, so registration does not explain
    the F356W on-arc spread.
- **Possibly affected at a smaller level:** the UNCOVER→ALT pairs of #261, #263–#265, through ALT's own uncorrected
  0.013″ shift. They were not re-measured.
- **Review fixes (post hoc, labelled in RESULT §5):**
  - rule 1 is now applied per axis as written;
  - a clip floor stops a fit from collapsing to one source (F115W had);
  - the pixel scale is corrected.
  No threshold changed; no new candidate.

## 2026-10-16: X2 S6a archival brief (Markdown + PDF) and README results entry; S6a stays paused
- **Added:** `docs/results/s6a_5d_rotating_wormhole_modes.md` and `.pdf` (3 A4 pages), a short brief of the S6a
  result built only from `X2_s6a_summary.md`, `X2_s6a_repro.md` and `X2_s6a_prior_work.md`; links pinned to a57b824
  (S6a files unchanged on main 3a80e18). No computation re-run, no new physics.
- **PDF route:** pandoc `--mathml --section-divs` to standalone HTML with embedded A4 CSS, printed with the Playwright
  Chromium (`--headless --print-to-pdf --no-pdf-header-footer`); pages checked as images, links checked (HTTP 200).
- **README:** new "Research results" section with the one S6a entry (scope, PDF, brief, reproduction package; no
  human expert review yet).
- **Handoff:** S6a remains paused; the next step is still a human expert review, the owner's call.

## 2026-10-16: Abell 2744 `arc_diff`: a local-annulus cut for hits on extended light doubles 28-mag efficiency but lets artifacts through; not adopted
- **Hypothesis before measuring (#264 handoff):** the peer cut vetoes faint point brightenings on smooth light because
  its brightness-bin spread is set by compact sources; a cut that keeps the peer test only at compact epoch-1 peaks
  and otherwise tests S against the local S scatter recovers them without new false passes. Decision rule
  (ASSUMPTION, written before the run): adopt if the 28-mag chain efficiency rises to ≥ 0.30 with 27.5 mag ≥ 0.60
  **and** real morphology passes surviving the peer stage on the same 4 tiles at most double.
- **Tool (new, derived):** `scripts/arc_diff_extended.py` (compact = epoch-1 local maximum within 0.15″ standing out
  from the 0.3–0.9″ annulus; otherwise z_ann = (S_hit − annulus median S) / annulus robust σ, other hits masked,
  cut |z_ann| ≥ 3); `arc_diff_inject_chain.py --extended` counts it as stages `peers_ext` / `confirm_ext`;
  `scripts/arc_diff_extended_eval.py` applies it to the real passes. Unit-tested on synthetic maps.
- **Run:** the 8 #264 injection runs and the 4 real tiles re-made with `arc_diff.py` (commands in the committed
  `summary.json` files; S3 byte ranges, maps in the scratchpad only). The injected stars and hits reproduce #264 exactly.
- **v1 (pre-written; compact = peak ≥ 2 × annulus median):** 28 mag 0.19 → **0.39**, 27.5 mag 0.64 → 0.68; real
  peer-stage passes 31 → **89** (×2.9) → rule fails. Looked at the 58 added passes: mostly fades on isolated faint
  compact sources. The ratio test fails because the mosaics are not background-subtracted (sky pedestal: median
  epoch-1 contrast of real passes 1.26), so faint compact sources counted as extended.
- **v2 (post hoc fix; compact = peak ≥ 5 annulus σ above the annulus median):** 28 mag **0.375**, 27.5 mag 0.675; real
  31 → **62** (×2.0, on the rule's edge; o001-confirmed 16 → 45). Looked at all 31 added passes
  (`arcdiff_tiles/extended/sheet_newly_passing.png`): detector stripes (h270, h299, h331, h412), spike fragments (h183,
  h321, h638), faint compact fades, knots on galaxies, and 7 positive hits that are all compact-core residuals of
  bright galaxies (the hit is > 0.15″ from the peak, so no compactness is measured). **None is a point brightening on
  an arc. Not adopted:** v2 is post hoc and the sheet shows that the peer cut was also doing artifact rejection.
- **Failed approach:** the first runs of this cycle used 4 `arc_diff.py` processes without `OMP_NUM_THREADS=1`; load
  15 on 4 cores, none finished in 16 CPU-min each; with one thread each a run takes ≈ 1.5 min (scripts/CLAUDE.md
  rule, now confirmed for `arc_diff.py`).
- **Limit of the comparison:** the real "old" verdict is rebuilt as |z_peers| ≥ 3 or no peers from `passes.ecsv`;
  over all 667 passes this keeps 185, not the 196 #261 reports (the injections use the stored `stored_pass`).
- **Consequence:** #261's "0 lensing candidates" still holds only at ≳ 36 nJy (#264). Outputs:
  `arcdiff_tiles/extended/` (`inject_summary.json`, `real_passes_4tiles.ecsv`, the sheet).
- **Handoff:** a replacement needs explicit stripe and spike vetoes (row/column-coherent S; position along a bright
  star's spike angles) and a compactness rule for hits offset from their source's peak before the annulus cut;
  re-calibrate on these injections and re-check the sheet. Alternatively accept 36 nJy as the tiling's limit and move
  on: the PEARLS 1176 registration offset (#263) is still unmeasured.

## 2026-10-16: Abell 2744 UNCOVER vs ALT `arc_diff` chain completeness: the peer cut halves 28-mag efficiency
- **Question:** what fraction of point brightenings on arcs and galaxy light survives the #261 tiling chain (hit,
  not `near_bright_core`, morphology, peer-population |z| ≥ 3, same-sign S ≥ 3 in UNCOVER→ALT o001)? No decision
  rule was written beforehand; this is a calibration run, and it changes the plan only if one stage dominates the losses.
- **Run (simulated stars in real maps, derived):** `arc_diff.py --inject-grid 40` (27.5 / 28.0 AB = 36 / 23 nJy, on
  pixels above the epoch-1 90th percentile, 1.5″ apart) on 4 tiles (t−1+0, t+0+0, t+0−1, t+1−1), UNCOVER→ALT o002
  and →o001 with the same stars (same epoch 1 and seed), tiling settings of #261; ≈ 2 min per run, 4 in parallel.
  New `scripts/arc_diff_inject_chain.py` runs morphology and peers and counts each stage
  (`results/transients/a2744_grid_alt/arcdiff_tiles/inject/summary.json`).
- **Result (80 stars per magnitude):**

  | mag | hit | unflagged | morph pass | peer pass | o001 confirm | chain efficiency |
  |---|---|---|---|---|---|---|
  | 27.5 | 77 | 65 | 54 | 51 | 51 | **0.64** |
  | 28.0 | 63 | 47 | 31 | 15 | 15 | **0.19** |

  The o001 confirmation costs nothing. Morphology removes 11 / 65 and 16 / 47 unflagged stars. The peer cut removes
  3 / 54 at 27.5 but **16 / 31 at 28.0**: it is the main avoidable loss for faint changes.
- **Looked at (contact sheet `inject/vetoed_28mag.png`, every vetoed unflagged 28-mag star):** the peer-vetoed stars
  are clean isolated positive points on smooth galaxy outskirts and ICL (z_stored 1.2–2.9). The peer MAD of their
  brightness bin (2–7, set by compact sources) is too wide for them. Morphology vetoes are a neighbour's rosette
  within 0.8″ (13), a mosaic edge (4) or an extended residual (3), as the morphology calibration on MACS0416 found.
- **Failed approach:** `scripts/arc_diff_peers.py` first rebuilt the uncommitted #261 peer step from its description
  (local maxima of epoch 1 within ±0.1 dex). It does not reproduce the stored `passes.ecsv` values (cut agreement 0.35–0.54,
  median |Δz| 1.2–2.6, 2–23× more peers), so the chain uses `z_stored`: S against the stored peer median and
  MAD of the 3 stored passes of the same tile nearest in `log_peak`.
- **Consequence for #261 (no new candidate claimed):** its "0 lensing candidates" holds at ≳ 36 nJy (0.64
  efficiency) on these tiles, not at 23 nJy (0.19). Of its 667 morphology passes, 471 fell to the peer cut; |z| ≥ 2
  would keep 307 (+111), and 24 / 31 of the 28-mag stars.
- **Limits:** Gaussian stars (no spikes); 4 of 33 tiles; one band (F356W) as in the tiling; `z_stored` borrows the
  peer statistics of real passes, not of the star's own position; the o001 check reads the tile map, not the
  10″ boxes of #261.
- **Handoff:** replace the peer cut for hits on smooth extended light (e.g. use peers only when the epoch-1 position is
  a compact peak, and otherwise the local S scatter in an annulus), re-calibrate on these injections, then re-screen
  the 471 peer-vetoed passes. Measure the PEARLS 1176 registration offset (from #263) before any on-arc photometry.

## 2026-10-16: Abell 2744 h182 against arc-pixel controls: within on-arc systematics; closed, not a candidate
- **Hypothesis before measuring:** h182's 2022→2023 red-band brightening (+45 % F356W; compact-control z_cal +4.9,
  #261) is a PSF/registration effect on steep extended light that compact-peak controls do not sample. Decision rule
  (ASSUMPTION, written before the run): h182 stays open only if |z_cal| ≥ 3 against arc-matched controls in ≥ 2
  independent program pairs; otherwise it closes as on-arc systematics.
- **Tool (new):** `hit_lightcurve_controls.py pick-arc`: grid positions (0.2″ step) whose 0.15″ aperture flux **and**
  annulus median both lie within ×1.5 of the hit's own (UNCOVER F356W), ≥ 0.6″ apart and from the hit. Unit-tested
  on a synthetic band with knots (only on-band knots are picked).
- **Run (derived):** 80 picks in a 32″ square box centred on h182 (corners ≤ 22.6″); 78 measured in the same 64
  mosaics as h182 (2 fell in a separate position group without a coverage query, as c08 in #261). Positions and map:
  `vet/ctrl182arc/positions.txt`, `controls_map.png` (looked at: knots on the 301 arcs, including two on h182's own
  arc, and many on the steep outskirts of compact cluster members). 5,056 measurements, 0 skipped. Command in the
  script docstring.
- **Result (derived; `vet/ctrl182arc/control_z.ecsv`):** control robust std 1.6–3.8 per pair (the compact controls'
  spread was smaller). h182 z_cal: F356W UNCOVER→ALT o002 **+2.3** (9 % of controls as extreme), →ALT o001 +1.4,
  ALT o002→o001 −1.4, UNCOVER→PEARLS +0.7; F444W UNCOVER→PEARLS +1.6, →SAPPHIRES +0.3; F410M +1.4, F277W +0.9,
  F200W +0.4, F150W −0.2. No pair reaches 3: **closed under the pre-written rule; not a candidate.** The rule's
  "independent" pairs share the UNCOVER epoch; with no pair at 3 the outcome does not depend on that.
- **Pick variance:** a first run, before a review fix aligned the annulus edge with `aperture_photometry`, drew a
  different random 80 (47 shared): F356W UNCOVER→ALT o002 z_cal +2.8 (6 % tail), every other pair ≤ 1.6. Same verdict.
- **Looked at:** tail controls in that first run changed as much as h182 (c63 +8.7σ raw, c18 +7.9σ, c12 −9.3σ). c63's
  stamps show a galaxy-outskirt position, and its PEARLS 1176 cutouts are offset by ≈ 0.3″ against the cross
  relative to the other programs: a mosaic registration difference, which moves flux on steep gradients.
- **Sensitivity:** with this control spread a real change on h182's knot needs raw z ≈ 11 (≈ +31 nJy, +56 %) in
  F356W UNCOVER→ALT o002 to reach z_cal 3.
- **Blocked:** TNS for h003 (www.wis-tns.org answers 403 from its own load balancer, not the proxy); h003 stays an
  ordinary SN-like transient with that check undone.
- **Handoff:** the `--inject-grid` completeness run on 2–3 Abell 2744 arc tiles is next; the PEARLS 1176
  registration offset should be measured (one shift per mosaic) before any further on-arc photometry with 1176.

## 2026-10-16: Abell 2744 UNCOVER vs ALT F356W `arc_diff` tiling and multi-epoch vetting of 6 priority hits: one SN-like transient, one variable, 0 lensing candidates
- **Hypothesis before measuring:** difference imaging of the whole UNCOVER o001 (2022-11) ∩ ALT 3516 o002 (2023-11/12)
  F356W overlap recovers PLATE-23a and the grid's real changes g1 and g2, and every other morphology pass is a
  star/spike residual, an ALT single-visit artifact, a mosaic-edge effect or a registration residual on a core.
  Decision rule (ASSUMPTION, written before the run): same sign at |S| ≥ 3 in UNCOVER→ALT o001 (9 days after o002)
  **and** |S| < 3 in ALT o002→o001.
- **Run (derived):** `arc_diff.py --band F356W <UNCOVER t003> <ALT o002> --half-arcsec 30 --min-bands 1 --kernel-px 4
  --kernel-stamps --maps` on 34 tiles of 60″ (`arcdiff_tiles/tiles.txt`; S3 byte ranges, nothing stored), 4 jobs,
  ≈ 4 CPU-min per tile. 33 tiles ran; t−2+3 (footprint corner) crashed on an empty array (no common data; not
  investigated). `arc_diff_morph.py --min-rosette-bands 1 --min-dipole-bands 1 --min-ring-bands 1` (single band):
  9,196 hits, 3,665 without the bright-core flag; extended 1,275, rosette 1,149, edge 443, dipole 116, ring 15,
  **pass 667** (`arcdiff_tiles/passes.ecsv`).
- **Population test (derived, post hoc):** S of each pass against S at all compact epoch-1 peaks of similar brightness
  in its tile. Compact sources show a systematic UNCOVER→ALT deficit (peer median S down to −6.7): the MACS0416
  compact-core class again (PSF difference between programs). **Failed cut:** |z| ≥ 5 drops g1 (z +4.7) and g2
  (z −3.1); the cut used is |z| ≥ 3 (or no peers), calibrated on those two and PLATE-23a (z +29): 196 survive.
- **Failed approach:** the pre-written rule rejects **PLATE-23a** (ALT o002→o001 S −7.2: the SN fades in 9 days) and
  g1 (+5.4). Revised rule: same sign, |S| ≥ 3 in UNCOVER→ALT o001 only (small 10″ `arc_diff` boxes per pass):
  **76 pass** (`arcdiff_tiles/survivors.ecsv`, sheets `sheet_survivors_{a..d}.png`: epoch-1 image and S).
- **Looked at, all 76 (by eye, not vetted):** 20 spikes/edges, 16 compact-core residuals, 7 UNCOVER pixel artifacts
  (bad-pixel blocks visible in the UNCOVER image), the 3 known changes, 24 faint compact fades or brightenings
  (mostly isolated faint sources, ±5–9σ; transient or variable population not separated), and **6 priority hits**:
  - t−1+0 h182 (3.552148, −30.376850): +8σ point on a bright arc, 4.6″ from Bergamini+23 image 301.1;
    UNCOVER→ALT o001 +7.3;
  - t−1+0 h123 (3.552406, −30.373976) and h287 (3.554220, −30.372327): −10σ / −6σ on compact knots of an arc,
    6.9″ / 5.7″ from 301.1; o001 −7.8 / −4.2 (knots fainter in both ALT visits);
  - t−1+0 h294 (3.550096, −30.378202): +6σ blob on blank sky 7.1″ from 301.1; o001 +3.2, o002→o001 −5.7;
  - t+0−1 h003 (3.577785, −30.394141): +48σ new compact source 19″ from image 64.1; o001 +41.5, o002→o001 −15.9
    (fading between the ALT visits; SN-like). Not in the grid's 45 groups (a host blend, like PLATE-23a);
  - t+1+1 h011 (3.592116, −30.366745): +14σ, growing (o001 +40, o002→o001 +30); variable source 51″ from images.
- **Vetting of the 6 priority hits (derived; `scripts/hit_lightcurve.py`, new):** forced photometry (r 0.15″, annulus
  0.35–0.6″, error = robust spread of 300 random apertures in the same box) at each fixed position in **every** public
  NIRCam level-3 mosaic covering it (MAST query; 46–64 mosaics, 2022-10 to 2024-10: DD 2756, UNCOVER, PEARLS 1176,
  MegaScience 4111, 3538, 2883, ALT, SAPPHIRES 6434; S3 byte ranges, nothing stored). Outputs and stamp sheets
  (looked at, all 6): `arcdiff_tiles/vet/`. Controls: 53 compact peaks (20–150 nJy) within 16″ of the 301.1 hits
  measured in the same mosaics (`vet/ctrl301/`, `scripts/hit_lightcurve_controls.py`); each hit's epoch-pair z is
  calibrated by the controls' median and robust std (`vet/ctrl301/control_z.ecsv`).
  - **h003: SN-like transient (ordinary).** Host flat at ≈ 25–35 nJy (2022-10 → 2023-11-06, 6 programs); a point source
    appears by 2023-11-20 (MegaScience: F090W 307 ± 3, F182M 187 ± 5, F360M 79 ± 3, F460M 61 ± 8 nJy; blue), rose
    ≥ 170 nJy in F182M within ≤ 10 days, fades through ALT (F090W 216 → 144 nJy in 9 days) and is gone by 2024-10
    (F090W 34 ± 4). Seen in ≥ 15 bands and 3 programs: not an artifact. 19″ from image 64.1; not a lensing anomaly.
    Known-transient lists (TNS, the PLATE catalog in prep) not checked.
  - **h011: variable compact source (ordinary; AGN-like).** ±10–20 % changes on 9–13-day scales in 2022 and 2023
    (F277W 171 → 136, F356W 147 → 124 → 171 → 185, F360M 198 → 166 nJy), no rise-and-fall. 51″ from any image.
  - **h294: single-visit artifact.** Only ALT o002 F356W (17 ± 3 nJy); nothing in the same visit's F070W/F090W, in ALT
    o001 nine days later, or in any of 44 other mosaics.
  - **h287: consistent with arc systematics.** Calibrated |z| ≤ 3.3 both ways; shows the size of on-arc excursions
    (|z_cal| up to 3.3) that compact controls do not capture.
  - **h123: marginal fade.** Red bands 25–40 % fainter after 2022-11 (F444W z_cal −3.0 in 2023-07, −3.3 / −2.5 in
    2024-10), blue bands unchanged; within the on-arc excursion size of h287. Not significant.
  - **h182: inconclusive brightening, ordinary explanations open.** Arc knot 4.5″ from 301.1b; red bands 20–45 %
    brighter in 2023 than in 2022-11 (z_cal +4.9 ALT o002, +3.7 ALT o001, +3.0 F444W / +2.8 F410M PEARLS 2023-07),
    −2.6 between the two ALT visits, back to the 2022 level in 2024-10 (SAPPHIRES F444W/F200W z_cal ≤ 1.3); no
    visible new point source in the stamps. Bergamini+23 gives μ 71 (z 3) to −359 (z 4) here (system 301 z 5.29;
    model_prediction): near a critical curve, so microlensing of a compact knot is the leading ordinary
    explanation, with the on-arc photometry systematic (h287) as the other. Microlensing changes one image only, so
    the lens model makes no counter-image prediction to test. Not a candidate.
- **0 lensing candidates.** Limitations: controls are compact peaks, not arc pixels (on-arc systematics are larger);
  the pair choice and calibration are post hoc; one control (c08, > 20″ from the hits) got no mosaics under `--query-only` and is absent; no injection through the tiling chain; the 24 faint open hits of the
  tiling are not classified; program-dependent PSFs not modelled.
- **Handoff:** h182 is the only open item (low priority): an arc-pixel control set (`arc_diff_vet.py` arc noise) on the
  same mosaics would decide whether +45 % exceeds on-arc systematics; a TNS check for h003. Then the `--inject-grid`
  completeness run on 2–3 tiles.
## 2026-10-16: S0 check: message-and-reply window in an unstable wormhole: literature already answers it; stop, no new unit
- Literature and code check only (≤ 2 h; no solver, no evolution). Record: `docs/state/units/S0-2026-10-16.md`.
  The S0-2026-10-15b radiation-fluid eigenvalue computation is cancelled (owner accepted #256).
- Published (Ellis + ghost scalar, two separate universes): self-gravitating ordinary pulses cross, then collapse
  the throat (Shinkai–Hayward 2002; Doroshkevich et al. 2009; Calhoun–Fay–Kain 2022). Test light rays return
  (CFK22 Fig. 3) or make a round trip (a possible path, González–Guzmán–Sarbach 2009 Fig. 4) before the black hole forms. The window is
  of order r_throat/c (e-folding 0.846 r_throat/c) and grows only logarithmically for weaker perturbations.
- Every published scheme for longer traversability uses fine-tuned **active phantom injection** (SH02 §III D; Xu–Chew–Yeom 2025), and is only
  temporary. No paper computes a finite-energy reply, receiver redshift, message capacity or a same-exterior geometry,
  so no delay advantage over a 20–100 ly path can be stated.
- Code: Shinkai's `wh_dual` is public (2002, page updated 2007) but has **no licence**; not downloaded. The other codes are not stated as public.
- Verdict: stop. The remaining computable gap (a two-leg self-gravitating exchange with redshift) would not change the
  feasibility judgement. Reopen conditions are in the record.

## 2026-10-16: X3 reopen check: DR4 epoch-astrometry fields confirmed (draft data model), precision not yet published; X3 stays paused
- Documents only (≤ 30 min; no Gaia data read; no X3 computation). Record: `docs/state/units/X3.md`, "DR4 reopen
  check".
- The official draft DR4 data model (2026-06-26) defines `centroid_pos_al`, `scan_pos_angle` and
  `centroid_pos_error_al` per CCD transit in `epoch_astrometry`, for all sources in the main astrometric
  processing.
- No official page or document read gives per-transit AL precision against G. The data model is still a draft.
- Verdict: partially confirmed. Waiting for the DR4 release documentation.

## 2026-10-15: Process: `merge=union` for CHANGELOG.md and trial.md; operations notes on background jobs
- `.gitattributes`: `CHANGELOG.md` and `docs/state/trial.md` use `merge=union`.
  - Checked locally with two test branches: both entries kept, no conflict markers, no duplicated lines.
  - Limits: a blank separator line can be lost; both versions of a line edited on both sides are kept silently.
    It applies to local `git merge` only; GitHub's web mergeability check may not use it.
- `docs/operations.md` §11: the session's `pkill -f` / wait-loop self-matches, per-row result writing, and
  background-wait timeouts, each with the fix used. No research content changed.

## 2026-10-15: Abell 2744 UNCOVER vs ALT F356W: grid misses PLATE-23a, `arc_diff` recovers it; 0 lensing candidates
- **Hypothesis before measuring:** the single-band grid UNCOVER 2561 o001 (2022-11) vs ALT 3516 o002 (2023-11/12)
  recovers PLATE-23a (arXiv:2608.20781; F356W < 29.7 → ≈ 26.1 mag); every other hit is a spike, star core, persistence
  or a single-visit artifact. Plan change: a missed known case means the grid has a blind spot that must be named.
- **PLATE-23a position checked:** the printed "03h33m09.92s" is degrees (3.552756°, −30.400340°); a point source
  appears there in ALT and is absent in UNCOVER.
- **Run (derived):** `transient_grid.py --min-bands 1 --inject 500 --inject-frac 1.0 0.5`; UNCOVER cut with
  `i2d_region.py` (207 MB, 16 s). 6,440 positions, 45 hit groups; injections −0.75 mag 95.6 %, −0.44 mag 50 %.
  Outputs `results/transients/a2744_grid_alt/` (`grid.ecsv` > 1 MB not committed). Control grids: UNCOVER vs ALT o001
  and ALT o002 vs o001 (9 days; subdirectories).
- **Failed approach / limitation:** the grid **misses PLATE-23a**: no catalog entry (blended into its host 2.2″ away).
  `arc_diff.py --min-bands 1` on a 24″ box finds it at +50σ, 0.06″ off (`arcdiff_plate23a/`).
- **Vetting (contact sheets looked at, all 45):** 16 spikes/wings, 9 star cores, 2 dipoles, 3 UNCOVER streaks, 1 edge;
  9 ALT o002-only pixel blocks, all absent in ALT o001 a week later; 3 UNCOVER F356W-only blocks, absent in the same
  visit's F277W/F410M/F444W. Two real changes, both ≥ 45″ from any Bergamini+23 multiple image (not lensing anomalies):
  - **g2** (3.613638, −30.416380): on sky in all 7 UNCOVER bands (12–43 nJy, r 0.25″), band centroids agree ≤ 0.03″
    across exposures hours to days apart (not a moving object), gone in both ALT visits (F356W −4 ± 3 nJy): a faded
    transient, SN-like, class not determined; not in the PLATE paper (its catalog is in prep). The GLASS 1324,
    DD 2756 and UNCOVER o002/o003/o006 F356W mosaics do not cover g1 or g2; SAPPHIRES 6434 F356W mosaics were not
    read (WCS conversion failed in the byte-range script) and MegaScience was not checked.
  - **g1** (3.562169, −30.400185): compact, F356W r 0.4″ 67 ± 2 → 111 ± 5 / 112 ± 6 nJy (+0.55 mag) in both ALT visits;
    variable AGN or a transient on a faint host; not vetted.
- **0 lensing candidates.** Data (≈ 1.8 GB of mosaics) streamed and deleted after the run.
- **Handoff:** next, a difference-image pass over the whole UNCOVER∩ALT footprint (`arc_diff.py` tiles,
  `--min-bands 1`) to close the host-blend blind spot; vet g1 with SAPPHIRES 2024–25 F444W/F200W if they cover it.

## 2026-10-15: Sunrise and SMACS two-epoch grids re-run under the MACS0416 defaults: 0 candidates
- **Hypothesis before measuring:** the first Sunrise o010/o120 run and the #237 runs (Sunrise o010/o052,
  SMACS o001/VENUS o057) used per-band recentring and the same-sign rule, which missed MACS0416 SN01 (#238). Under
  `--recentre-in origin` and mixed-sign hits, a colour-changing or neighbour-confused transient could appear.
  Known case: SN Atalanta (#237) must be recovered. Plan change: any new row that is not a star, spike or
  persistence goes to `/vet-candidate`.
- **Run (derived):** `transient_grid.py`, same bands, catalogs and injections as before; outputs
  `results/transients/{sunrise_grid_v2,sunrise_grid_o052_v2,smacs_grid_v2}/`.
  - SMACS/VENUS: 16 groups; all 15 earlier groups re-found, SN Atalanta A/B at −3.3/−4.7 mag, ≥ 14σ in both bands;
    2 new rows (g04188, g04192): bright-star proper-motion dipoles. Injection −0.75 mag 97 %, −0.44 mag 58 %.
  - Sunrise o010/o052: 10 groups; 9 of 10 earlier groups re-found (the lost one was a spike on a galaxy); new
    g03622 lies on a bright star's diffraction spike that moves between epochs. Injection 98 % / 68 %.
  - Sunrise o010/o120: 38 groups; all 36 earlier groups re-found; new g00650 (star spike) and g02509 (streak),
    both moving with the roll angle. Injection 99 % / 58 %.
- **Contact sheets looked at** for every new row. **0 candidates.** One mixed-sign row in total (SMACS g00612, a
  proper-motion star).
- **Failed approach:** my first contact sheets cut each epoch in its own pixel frame; Sunrise o010 and o120 differ
  in roll, so the difference stamps showed a false 2″ shift. Resampling epoch 2 onto epoch 1's WCS fixed it (the
  catalogs agree to ≤ 0.02″ in every field region).
- Data: 16 `_i2d` mosaics (≈ 9 GB) and their catalogs streamed from MAST, deleted after each run; `grid.ecsv`
  (> 1 MB) not committed.
- **Handoff:** the grid line is null under both rule sets in four clusters. Next per TASKS.md: Abell 2744 UNCOVER
  vs ALT F356W with PLATE-23a as the known case, or arc differencing on the Sunrise arc.

## 2026-10-15: MACS0416 arc tiling o211→o213, 28 / 28 tiles: 27 morphology passes, all ordinary, 0 candidates
- **Hypothesis before measuring:** o211→o213 shows the same artifact classes as o211→o212 (spikes, rosettes,
  persistence, the o211 SW compact-core deficit of #251), so its passes are either o211-epoch features already seen in
  o211→o212 or o213 detector effects. Plan change: a pass with no o212 counterpart that is on the sky, in ≥ 2 bands
  with LW support and not persistence goes to `/vet-candidate`.
- **Run (derived):** `arc_diff.py --kernel-px 4 --kernel-stamps --maps` on the 28 tiles of
  `tiles_o211_o212/tiles.txt`, F150W–F444W, `--procs 1` × 4 jobs with OMP=1 (2–16 min per tile; slower while the
  persistence check ran). Outputs: `results/transients/macs0416_arcdiff/tiles_o211_o213/<tile>/`, morph records
  `morph/tile_<tile>_o211_o213.json` (current `arc_diff_morph.py` defaults). 23 tiles: 1,433 hits, 330 without the
  bright-core flag; verdicts extended 164, edge 66, rosette 59, sign_flip 12, dipole 5, ring 1, **pass 23**.
- **The 23 passes, looked at (contact sheets `vet_o213/sheet_passes_{a,b,c}.png`: epoch-1 image, difference, S in
  F150W/F200W/F444W):**
  - 8 coincide (≤ 0.3″) with an o211→o212 hit (h007/h008 t−3−3, h023 t−3−1, h017 t+0−1, h020/h031 t−2+0,
    h021/h022 t+0+2, h029 t+0−3; t−1−2 h087 at 0.11″): o211-epoch features, the SW-positive core class of #251;
  - 5 in t−2−2 (h028/h029/h035/h039/h040): one negative diagonal streak on a diffraction spike (broken, so the
    extent test misses it);
  - h045 (t−2−1) and h087 (t−1−2): PSF residuals on bright stars (SW and LW of opposite sign); h037 (t−1−1): along
    the edge of a bright elongated galaxy, LW only;
  - h147 (t−1−4, SW +6.8/+8.4): **persistence** (`persistence_check.py` on o211/o212/o213 `_cal`: 5 detections,
    all suspect, earlier o213 NRCB4 exposures saturated that pixel; `vet_o213/persistence_*.ecsv`);
  - h014 (t−3−3, F150W/F200W/F277W +6.8/+5.1/+4.4): `on_sky` (47 clean detections), a faint compact source present
    in every epoch; its per-exposure `_cal` fluxes (not recentred) do not show an o213 brightening (F150W
    1.52 ± 0.09 o211, 1.28 ± 0.09 o212, 1.33 ± 0.04 o213; F200W 0.98 ± 0.05, 0.83 ± 0.04, 1.07 ± 0.07), so the
    difference-image signal is not confirmed; **open, low priority** (`cal_dither_phot.py` would settle it);
  - h053 (t+1−1), h035 (t+0+2), h021/h023 (t+0−3): SW-dominant positives on compact cores or knots without an
    o212 hit, the same class as the o211 SW deficit; not separately vetted.
- **Reading:** no candidate in 23 / 28 tiles. The o211 SW compact-core deficit dominates the remaining passes, as in
  o211→o212; detector-dependent PSF stays degenerate with epoch (o211 is all NRCB1).
- **Last 5 tiles (t+1+2, t+2−2, t+2+0, t+2+1, t+2+2; same command and defaults):** 160 hits, 21 without the
  bright-core flag; extended 0, edge 7, rosette 4, sign_flip 3, dipole 3, **pass 4** (sheet `vet_o213/sheet_passes_d.png`):
  t+1+2 h006 and t+2−2 h011 coincide (≤ 0.03″) with o211→o212 hits (h009, h014): o211-epoch SW-core features;
  t+2−2 h019: adjacent +/− pair on a compact core, SW only, LW null: same class; t+2+2 h002 (F150W −35σ, F200W −7σ,
  LW null; a point source in the o211 F150W image only): **persistence** (`persistence_check.py`, 64 `_cal`
  exposures: one detection, in o211 NRCB4 dither 4, and dither 3 had put 26 saturated pixels on that pixel; the other
  o211 F150W dither has the position on the array edge; o212/o213 flat; `vet_o213/persistence_*_h002.ecsv`).
- **h014 closed:** `cal_dither_phot.py` (16 knots within 14″, `knots_o211_o213.json`; `vet_o213/dither_h014.json`):
  o211→o213 excess in the r 0.08″ core aperture (3 + 3 exposures) F150W +3.0 ± 4.0 %, F200W −8.6 ± 8.4 %; no
  band-consistent brightening. The r 0.16/0.24″ values rest on 2 o213 exposures (masked pixels) with errors of
  ±11–80 % and opposite signs between epoch pairs, so their large `z_vs_knots` (knot scatter only) is not evidence.
  Not a candidate.
- **Totals, 28 tiles (derived):** 1,593 hits, 351 without the bright-core flag; extended 164, edge 73, rosette 63,
  sign_flip 15, dipole 8, ring 1, pass 27; every pass looked at; **0 candidates**. Limitation: the o211 SW
  compact-core class is not separately vetted (detector-dependent PSF stays degenerate with epoch).
- **Handoff:** o211→o213 tiling closed. Next per TASKS.md: re-run Sunrise/SMACS under the new `arc_diff` defaults.

## 2026-10-15: MACS0416 h044 per-dither `_cal` photometry: o211 deficit in every dither but not band-consistent; not a candidate (closed)
- **Hypothesis before measuring (from the previous handoff):** a deficit in some o211 exposures only is processing
  (resampling, outlier rejection); added point-source flux on the sky shows in every o211 exposure, in both SW bands
  at the i2d level (+4–7 % over the knots), and does not shrink with aperture radius. Plan change: if neither holds,
  h044 closes as not a candidate.
- **New `scripts/cal_dither_phot.py` (offline tests):** for every level-2 `_cal` exposure covering the target and the
  knots (`knots_o211_c004.json`, 150–1000 nJy, 17 knots): per-exposure knot offset (level-2 WCS lacks tweakreg; found
  ≤ 0.5 px here), recentred aperture fluxes at 0.08/0.16/0.24″, a per-exposure knot normalisation (leave-one-out for
  knots), epoch excess ranked against the knots, and sky-centroid shifts relative to the knots. A first version
  without the offset step and with a 1.5 px recentring limit gave a spurious +27 % o211 excess (h044's core sits
  ~1.2 px from the hit position, so recentring failed in some epochs): **failed approach, fixed**.
- **Data (derived):** o211 (4 dithers, all NRCB1), o212 (4, NRCB2), o213 (4–6, NRCB4/B2), CANUCS o004 (3–7), F150W
  and F200W; `results/transients/macs0416_arcdiff/vet_sw/dither_h044_h030.json` (per-exposure `t` included; epochs
  are labelled `<program>o<obs>`, e.g. `01176o211`, `01208o004`; pairs `a_b` are chronological, excess = b / a − 1;
  `n_measured` gives the exposures behind each number: the F200W r 0.16″ o211→o212 value rests on 3 + 3).
- **Result:** h044's o211 deficit is the same in every measured o211 dither (F200W r 0.16″: t = 0.93, 0.96, 0.93;
  the 4th dither's aperture has a masked (DO_NOT_USE or NaN) pixel; vs 1.00–1.05 in o212/o213/o004), so it is **not** a
  single-exposure or mosaic artifact. It is not consistent with
  added point flux: F200W +3.6 ± 2.8 % (r 0.08″, z 0.9), +9.0 ± 1.8 % (0.16″, z 3.1), +6.2 ± 2.2 % (0.24″, z 2.2) for
  o211→o212; F150W +3.7 % (0.08″, z 1.0) but −1.5 % and +1.7 % at 0.16/0.24″ (a hot source would be brighter in F150W;
  LW was already null; the F150W 0.16/0.24″ and F200W 0.24″ values rest on 3 + 2 exposures, so their errors are
  rough). o212→o213 control ≈ 0 in both bands. Centroid: h044 moves 9 mas west from o211 to later
  epochs in F200W (z −2.6), 4–5 mas in F150W (z −1.2), but also 12 mas south from o212 to o004 (z −2.5, no o211),
  so its centroid is less stable than the knots in general (compact core on a host gradient).
- **/code-review fixes (two passes):** sources whose recentring fails are dropped per exposure (no fallback centroid
  or flux; every h044/h030 number quoted above is unchanged); epoch labels carry the zero-padded program and pairs
  are ordered by time (the o004 pairs' keys and signs changed accordingly); per-target `n_measured`; NaN written as
  null (strict JSON); empty runs exit cleanly; target/knot id clash refused; one shared NaN-safe `robust_z`. `arc_diff_knots.py`: knots whose aperture touches a NaN pixel are excluded and targets get an
  `aperture_has_nan` flag, a NaN target fraction has no rank, and skipped targets are recorded. The stored
  `knots_*.json` predate these two rules (their maps were deleted); their h044/h030 rows were finite.
- **Reading:** h044 is **not a candidate (closed)**: a ≤ 3σ tail of the 17-knot distribution, picked as the most
  extreme of the o211→o212 hits over 28 tiles, with an aperture- and band-inconsistent excess. Untestable ordinary
  remainder: every o211 exposure is on NRCB1 and no later one is, so detector-dependent PSF is degenerate with epoch
  here. h030 stays ordinary (F150W excess follows o213, not o211). **All 9 open o211→o212 passes are now closed.**
- **Handoff:** next, o211→o213 arc tiling over the 28 tiles (`arc_diff.py --maps`, then `arc_diff_morph.py`), with
  `cal_dither_phot.py` as the per-dither check for SW-only survivors; HST pre-JWST epoch for h044 only if the
  owner wants it.

## 2026-10-15: MACS0416 h044/h030 knot test on o211→c004 and the o212→o213 control: h030 ordinary, h044 open
- **Hypothesis before measuring:** if the o211 SW core deficit explains h044/h030, both sit inside the
  compact-knot distribution in o211→c004 (an independent program), and neither stands out in o212→o213 (no o211).
  Plan change: a target that is the top knot in o211 pairs only stays open; one that also stands out without o211,
  or sits inside the knots, is ordinary.
- **New `scripts/arc_diff_knots.py` (tested offline):** the previous run's ad-hoc knot measurement as a script
  (epoch-1 peaks, footprint 0.25″, `arc_diff_vet.py` apertures, targets ranked against knots of 150–500 nJy,
  excluding knots within 0.2″ of the target; robust σ). **Validation:** the h044/h030 aperture fluxes of the
  rerun boxes equal the stored `vet_t32h044.json` / `vet_t33h030.json` values (e.g. h044 o211→c004 F150W 15.42 /
  281.7 nJy). Below 150 nJy, F150W picks up hundreds of noise peaks on extended light, so 150 nJy (the previous
  run's floor) is kept (ASSUMPTION).
- **o211→c004 (derived):** SW knot medians +1.7 % (F150W, n = 16) and +1.5 % (F200W, n = 11), LW +0.5 %: the
  o211 SW deficit is confirmed against CANUCS. h044 is the top knot in both SW bands (F150W +5.5 %, z 2.7;
  F200W +8.75 %, z 4.4); h030 is inside (z 1.6 / 1.3). LW ordinary for both.
- **o212→o213 control (derived):** F200W knot median +0.1 % (the F200W deficit is o211-specific), but F150W
  knots still give +1.3 % (+3.1 % for the two brightest), so part of the F150W offset is a general kernel-matching
  systematic. h044 is ordinary (F150W z −0.3, F200W z −1.7). h030 is the top F150W knot here too (+4.4 %, z 4.5).
- **o211→o213 (derived; finished at the end of the run):** SW knot medians +2.0 % / +1.8 %. h044 is the top knot
  in both SW bands again (F150W +8.0 %, z 4.0; F200W +6.4 %, z 2.2). h030 is the top F150W knot (+6.7 %, z 3.1),
  F200W z 1.8; its F150W excess tracks o213 (also in o212→o213), not o211.
- **Looked at the cutouts** (`sheet_h044_h030_c004.png`): both residuals sit on compact cores and are
  dipole-like (adjacent positive and negative lobes in F150W/F200W), the pattern of a core PSF or sub-pixel
  registration mismatch, not an isolated point source.
- **Reading:** h030 is **ordinary** (inside the knots against an independent program; its F150W excess appears
  without o211; dipole residual). h044 stays **open and unvetted, not a candidate**: an o211-only SW excess of
  about 5–7 % over the knot median (≈ 13–18 nJy F200W) on a compact core with a dipole-like residual. A hot point
  source would also fit SW-only (Rayleigh–Jeans: F356W ≈ F150W / 5.6, within the LW noise) — an ordinary
  lensed-star explanation is untested, not excluded. h044 is the top SW knot in all three o211 pairs and
  ordinary in the control. Limits: one box; aperture ratios, not PSF fits; a 4-band `--kernel-stamps` box takes
  7–22 min on 4 cores.
- **Handoff:** (1) per-dither o211 `_cal` photometry of h044's
  core against two nearby knots (a dither-dependent deficit is processing); (2) a sub-pixel-shift test: refit
  h044's core with a free 0.01–0.03″ offset and see if the net flux excess survives; (3) HST pre-JWST epoch.

## 2026-10-15: MACS0416 h044/h030: o211 compact knots are systematically fainter in SW; explains about half of h044
- **Hypothesis before measuring:** the o211-faint SW signature of h044/h030 is an o211-wide compact-core deficit
  (PSF or processing), so other compact knots in the same box show the same fractional residual. Plan change: if
  h044/h030 lie inside the knot distribution, they are ordinary; if they are clear outliers, go to per-dither tests.
- **Method (derived):** `arc_diff.py --maps` on the west box (64.0225 −24.0815, ±10″) for o211→o212 (F150W, F200W,
  F277W, F356W; `--kernel-px 4 --kernel-stamps`); difference / epoch-1 flux in the `arc_diff_vet.py` apertures at
  every epoch-1 local peak (0.25″ footprint) above 100 nJy. Table and plot:
  `results/transients/macs0416_arcdiff/vet_sw/knots_o211_o212.json`, `knots_frac_o211_o212.png`.
- **Result:** bright compact knots are fainter in o211 in SW only. Knots with an epoch-1 aperture flux ≥ 500 nJy
  give +1.3 to +3.6 % (F150W, up to 30.8σ) and +0.4 to +2.4 % (F200W); 150–500 nJy knots have a median of +1.8 %
  (F150W, n = 18) and +1.5 % (F200W, n = 13). LW knots scatter around 0. At the hit positions h044 is +6.4 % / +6.6 %
  (F150W / F200W) and +4.4 % / +5.4 % on its own epoch-1 peak 0.05″ away; h030 is +3.4 % / +4.8 %. The largest other
  knots are +3.6 % / +2.8 %. Only the F200W values stand out (h044 and h030 are the two largest, about 2.5 robust σ
  above the knot median).
- **Reading:** an ordinary o211 SW core deficit (cause not identified: PSF, focus or the kernel fit led by extended
  light) accounts for roughly half of h044's signal and most of h030's. The excess (≈ 2–3 %, 7–9 nJy) is **open and
  unvetted; not a candidate.** Limits: one epoch pair (the o211→o213 run was stopped for CPU time), one box, the
  ratio is of a 0.08″ aperture, not of a PSF-fitted flux.
- **Handoff:** (1) the same knot test on o211→c004 and o211→o213 (is the deficit o211-specific?) and on o212→o213 as a
  control; (2) per-dither o211 `_cal` differences registered to the i2d; (3) HST Frontier Fields pre-JWST epoch.
  CPU note: two `--kernel-stamps` runs at once on 4 cores took > 20 min each (load 15); run them one at a time.

## 2026-10-15: MACS0416 o211→o212 morph passes vetted: 7 ordinary; h044 and h030 are an open SW signature (unvetted)
- **Hypothesis before measuring:** the 6 SW-dominant passes (t−3−3 h028/h030, t−3−2 h044/h051, t−2−2 h035, t−1+0
  h104) are o212 afterimages or matching residuals, and t−1−4 h124/h159/h246 are o212 LW core-PSF mismatches. A real
  point transient would keep its sign and flux against independent epochs. Plan change: a hit that keeps a ≥ 4σ
  forced flux in pairs without o212 goes to `/vet-candidate`.
- **Method (derived):** `arc_diff.py --maps` in three boxes (west 64.0225 −24.0815, south 64.0326 −24.0889, north
  64.033681 −24.069529) for o211→o212, o212→c004, o211→c004, o212→o213, o211→o213 (CANUCS o004 does not cover h028);
  `arc_diff_vet.py` on each hit. JSONs, box summaries and the west-box sheet: `results/transients/macs0416_arcdiff/vet_sw/`.
- **Ordinary (not candidates):**
  - h124, h159, h246 (t−1−4): an **o212-only LW dip** (F277W/F356W −8 to −19 nJy in o211→o212, +7 to +14 in
    o212→c004 and o212→o213, ≤ 4 nJy in o211→o213) on three separate compact cores at once: an o212 LW matching
    artifact; h246 is arc-coherent in three pairs.
  - h051: a c004-only F200W deficit (−11.5 ± 3.8 in both pairs with c004; o211→o212 +2.6 ± 2.6); its o211→o212 S sits
    on a SW streak (sheet).
  - h035: forced fluxes ≤ 1.7σ against arc noise in every pair; c004/o213 pairs show rings around the core.
  - h104: ≤ 2.2σ in o211→o212 in the re-centred box (S +3.4 / +4.1); o212→o213 arc-coherent.
  - h028: o211→o212 F150W +3.6 ± 1.5 only; o212→o213 rings around the compact source.
- **Open, unvetted (not candidates):** **h044** (64.022499, −24.081711) and **h030** (64.023039, −24.081869), 1.9″ apart:
  o211 is fainter than **all three** later epochs in F150W and F200W with no LW change and no arc coherence.
  h044: o211→o212 +18.0 ± 3.9 / +17.2 ± 3.2 nJy, o211→c004 +15.4 ± 6.1 / +22.6 ± 3.4, o211→o213 +22.5 ± 5.1 /
  +16.6 ± 3.7; o212→c004 and o212→o213 within 1.7σ of 0. h030: o211→o213 +14.7 ± 3.2 / +8.8 ± 2.1, same pattern.
  The sheet shows sharp positive spots in every pair starting at o211, and sharp pixel-scale dipoles across the
  o211 SW stamps. Leading ordinary explanation (hypothesis): an o211-side SW deficit (detector or processing),
  shared by both. Alternative: a rise between 2022-10-07 and 2022-12-29 that then stays flat for 43 d.
- **Failed approach:** comparing per-exposure `_cal` aperture fluxes (`persistence_check.py` output) between o211
  and o212 on galaxies gives the opposite sign to the mosaic for h044/h035/h051; archive `_cal` WCS lack the level-3
  alignment, so the aperture lands on different parts of the galaxy. `persistence_check.py` returns `on_sky` for all
  nine because the host galaxy is detected in every exposure; it cannot test persistence on a galaxy.
- **Handoff:** vet h044/h030 first: per-dither o211 `_cal` differences in a box registered to the i2d (is the deficit
  in every o211 dither, one detector, or one exposure?), HST Frontier Fields F160W/F814W as a pre-JWST epoch, and
  whether other compact knots within 10″ show the same o211 offset. Then `/vet-candidate` if it survives.


## 2026-10-15: MACS0416 `arc_diff_morph.py` ring and sign_flip verdicts; all 28 o211→o212 tiles have morph records
- **Hypothesis before measuring:** t−1−4's passes are PSF-mismatch rings around galaxies; a ring test (and a
  SW/LW sign test) removes most of them while losing < 0.1 of injected point sources and none of Yan+2023. Plan
  change: if either rule vetoes a Yan+2023 transient or > 0.1 of injections, drop or loosen it.
- **t−1−4 re-run with maps (derived):** 46 unflagged passes without the new rules. Ring (fill ≥ 0.5 in ≥ 2 detected
  bands) and sign_flip (±5) leave 20; injections in the same tile (`--inject-grid 40` at 27.5/28.0 and at 29.0/29.5,
  simulated): of 56 unflagged recovered stars, ring removes 1 and sign_flip 2. Chain pass fraction in this
  artifact-heavy tile: 10 / 20, 10 / 20, 9 / 20, 5 / 20 at 27.5, 28.0, 29.0, 29.5 mag (rosettes take 8 of the 25
  faint unflagged). `results/transients/macs0416_arcdiff/inject_morph/summary_t-1-4_ring.json`.
- **Validation:** Yan+2023 11 / 11 pass on fresh o211/o212 Warhol maps; o211/o213 unchanged (D21-S2 rosette). A sign
  level of 3 was rejected: it vetoes D21-S1 (F200W −3.2, LW +7). Rejudge of the stored tiles: 10 former passes,
  already judged ordinary by eye (#245), become sign_flip (h081 and h074 among them); none else changes.
- **Looked at t−1−4's 20 remaining passes (contact sheet):** rings centred on a neighbour or cut by the window,
  streaks (h143, h253, h310) and three compact single-sign LW residuals on epoch-1 compact cores (h124, h159,
  h246); the last three are unvetted (likely core PSF mismatch; tile-wide SW rings around every compact source).
- **The 6 tiles without morph records** (t−3−3, t−3−2, t−3−1, t−2−2, t−2−1, t−1+0) re-run with maps: 6 passes, all
  SW-dominant with weak or no LW: t−3−3 h028/h030, t−3−2 h044 (SW +11.7/+11.3, LW < 2.2) and h051 (on a SW
  streak), t−2−2 h035 (F200W −3.0 against F150W +8.2), t−1+0 h104 (S/N 5.9). **Unvetted, not candidates**; sheet
  `results/transients/macs0416_arcdiff/morph/unrecorded_tiles_passes.png`.
- **Limits:** Gaussian injection PSFs; the ring centre is the brightest epoch-1 pixel within 1.2″, so rings around
  a fainter neighbour escape the test; stored JSONs of the 21 earlier tiles have no `ring_fill`.
- **Handoff:** next, vet the 6 SW-dominant passes (persistence on `_cal` frames first, then o213 / CANUCS o004 with
  `arc_diff_vet.py`) and h124/h159/h246; then o211→o213 over the 28 tiles.

## 2026-10-15: MACS0416 morphology veto calibrated by whole-chain injections; rosette rule now needs strong lobes
- **Hypothesis before measuring:** point transients on arcs and galaxies pass `arc_diff_morph.py` (≥ 0.9 of the
  recovered, unflagged ones) unless they sit in an artifact-dominated area. Plan change: if a veto removes more than
  0.1 of them, find which rule and fix it before quoting a limit.
- **New `arc_diff.py --inject-grid N`:** N Gaussian stars (27.5 and 28.0 AB in turn) at common sky positions on
  pixels above the epoch-1 90th percentile, 1.5″ apart, added to epoch 2 of every band before detection; the hits,
  maps and `injected.csv` carry them, so `arc_diff_morph.py --known injected.csv` tests the whole chain.
- **o211→o212, 11 tiles incl. the Warhol box (403 stars; simulated in real maps, derived):** 402 recovered as hits;
  326 without the bright-core flag. With the original rosette rule, 0.82 of those pass morphology (maps kept: 207 /
  251); the losses are rosettes 40, extended 10, dipole 3, edge 1. A no-injection control of t−1−3 shows the
  vetoing lobes are there without the star (lobe count = control + 1): neighbouring galaxy and BCG residuals within
  0.8″ in LW, not the star's own kernel ringing.
- **Fix:** a rosette counts only lobes whose peak is ≥ 0.3 |S_hit| (`--rosette-lobe-ratio`, default 0.3; the per-band
  `lobe_peaks` are now stored). Pass fraction 0.82 → 0.92 (231 / 251); it changes **none of 244 real-hit verdicts**
  on the same maps (0.5 would change 4). Yan+2023: 11 / 11 still pass (plain o211/o212 Warhol maps). Chain
  efficiency (hit, unflagged, morph pass) at 27.5 / 28.0 mag: 120 / 152 = 0.79 and 111 / 150 = 0.74 (maps kept).
  Summary: `results/transients/macs0416_arcdiff/inject_morph/summary_o211_o212.json`; ring sheet
  `t-1-4_real_passes_rings.png` there.
- **Looked at the remaining vetoed injections (contact sheet):** inside bright-star rosettes, on t−1−4's horizontal
  streaks and its ring residuals, on a compact core with its own −14σ SW residual. A veto there is the stated cost.
- **Gap found in the #245 record:** 7 of 28 o211→o212 tiles have no morph JSON (t−1−4, t−3−3, t−3−2, t−3−1, t−2−2,
  t−2−1, t−1+0). On the injection runs' maps, t−1−4 has 40 real morph passes: all sit on ~1″ ring residuals
  (opposite signs in SW and LW, or SW-only fading), an artifact class the rules do not cover; **not candidates**.
  t−3−3 has 2 passes, not yet looked at. So "28 / 28 tiles, 0 candidates" holds only for the 21 tiles with morph
  records until these are looked at.
- **Limits:** Gaussian PSFs (no spikes, so no self-rosette); only 11 tiles; two magnitudes; maps for 3 early tiles
  were deleted, so their ratio-0.3 pass rate is not measured.
- **Handoff:** next, a ring-residual rule (or a SW/LW sign-consistency rule) validated on t−1−4 and Yan+2023; then
  look at t−3−3 and the other 5 unrecorded tiles; then o211→o213 over the 28 tiles with `--maps` and the new default.

## 2026-10-15: `arc_diff_morph.py` review fixes; re-judged o211→o212 tiles: 4 new passes, all ordinary, still 0 candidates (#245)
- **Hypothesis before measuring:** the open review findings (top-band-only extent, lobes and edge counted in bands
  without the hit) can hide or invent vetoes; fixing them should change few verdicts. Plan change: any new pass not
  explained by an ordinary cause goes to `/vet-candidate`.
- **Fixes:** `classify` uses only detected bands (|S_hit| ≥ 3; the top band if none) for the edge, extent and rosette
  tests, and takes the extent over all of them. The dipole ratio uses the opposite region's local peak (within
  0.8″), not a far streak peak. Border targets become edges rather than being dropped, and NaN pixels are now
  `not_covered`. Per-band products are cached and computed band by band (memory), and the thresholds now live in
  `Params`. `--rejudge` re-classifies stored JSONs and refuses changed measurement thresholds. 7 new tests, each
  failing on the old code.
- **Re-judged 25 morph JSONs (stored per-band values; maps are gone):** Yan+2023 still 11 / 11 pass (o211/o212). Of
  the 13 changed verdicts, 9 became extended (t+0−3 h064/h067 on the o212 trail, t+0+2 h050 and six t−2−3 rosettes). The
  **4 new passes** are all ordinary (`results/transients/macs0416_arcdiff/vet_rejudge/`):
  - t+2+2 h007 (SW +6.2/+9.3): one o212 SW dither, near an NRCB3 corner. The same detector pixel (1977, 1996) shows S/N
    2.6, 3.8, 3.2 and 16.3 in the four o212 F150W dithers, which see four different sky positions. In F200W it then
    decays 8.8 → 4.9 → 2.8 → 0.6. That is detector-fixed persistence; a control pixel stays below 1.7σ.
  - t−2+0 h036 (SW −6.2/−5.1): the o211 F150W flux is only in the two dithers within 30 px of the NRCB3 right edge
    (S/N 4.8 and 2.4); the interior dithers give 0.6 and 0.8. That is a detector-edge artifact.
  - t+1+0 h081 and t−2−3 h074 change sign between detected bands (PSF and registration residuals).
- **Limits of the re-judge:** stored `dipole_ratio` values keep the old definition (in-disk maximum), and 8 early files
  have no dipole measurement; a fresh `--maps` run could add dipoles (fewer passes), never passes. t−2−3 keeps its
  pre-fix extents.
- **Failed approach:** `persistence_check.py` called with two directories after one `--cal-dir` (argparse error; repeat
  the flag). The script also only tests exposures where S/N ≥ 5, so a fixed-pixel scan of every exposure decided
  h007.
- **Handoff:** #245 should now be mergeable. Next: injection-calibrate the morphology veto through the whole chain,
  then o211→o213 over the 28 tiles (`--maps`, OMP=1).

## 2026-10-15: MACS0416 arc tiling o211→o212 (28 / 28 tiles): 0 candidates; residual-morphology triage `arc_diff_morph.py`
- **Hypothesis before measuring:** outside the Warhol box, the adopted `arc_diff.py` settings (#242) leave few
  unflagged hits per 20″ tile, and those are point-like. Plan change: if outer tiles are dominated by a new artifact
  class, build and validate a veto for it before screening further.
- **Tiles:** 28 boxes of 20″ (18″ step) holding the 106 CATS v4.1 multiple images outside the #242 box (shifted
  +0.208″, −0.025″ to the JWST frame); list in `results/transients/macs0416_arcdiff/tiles_o211_o212/tiles.txt`.
  o211→o212 done on 16 / 28 (first 10 with 30 injections per magnitude: 28.0 mag recovered 0.43–1.0 per band; the
  rest without injections). 1,353 hits, **390 without the bright-core flag** (t-1-4 alone 199, t-2-3 60).
- **Looked at t-2-3:** the unflagged hits are a diffraction spike of a bright star rotated between the epochs (a
  negative diagonal streak in every band), PSF rosettes around compact sources below the bright-core level, and
  residuals at the mosaic edge.
- **New `scripts/arc_diff_morph.py`** (unit tests on synthetic maps): per hit, the length of the S region containing
  it (S smoothed by σ = 2 px, ≥ 3σ), the number and signs of |S| ≥ 3 lobes within 0.8″, and the distance to missing
  data. Verdicts edge (< 0.5″), extended (> 3.0″ in the top band), rosette (≥ 4 lobes, both signs, ≥ 2 bands)
  (ASSUMPTIONs).
- **Validation on known cases (derived):** at 1.0″ the extent test flagged two published o211/o212 transients
  (D21-S1, D21-S2; 1.1–1.6″ on their arc), so the threshold was raised to 3.0″ (streaks measure 14–23″). At 3.0″:
  all 11 Yan+2023 positions pass in o211/o212; in o211/o213 one (D21-S2, outside its discovery pair) is a rosette.
  Not yet calibrated by injections through the chain (follow-up).
- **Morphology on 9 tiles with maps:** t-2-3 60 → 3 pass; over t+0±, t-1±, t-2-3, 15 passes remain, unvetted and
  not yet vetted. Looked at t-1-3's 7: low-S spots (mostly ≲ 6σ, two bands) on compact sources, in the BCG halo noise
  or beside the mosaic edge; none is point-like in every band. JSON: `results/transients/macs0416_arcdiff/morph/`.
- **Failed approach:** a 1.0″ extent cut (kills arc transients); rosette in the top band only (missed h074).
- **Handoff:** finish o211→o212 on the 12 remaining tiles with `--maps`, run `arc_diff_morph.py` on every tile, look
  at the passes (contact sheets), vet survivors with `arc_diff_vet.py` (CANUCS o004 as fourth epoch); calibrate the
  morphology veto with injections; then o211→o213. t-1-4 (199 unflagged) needs a look first.
- **Failed approach (run 016oUw9Y, 08:08–08:50Z):** two `arc_diff.py --procs 2` tile jobs in parallel without
  `OMP_NUM_THREADS=1` pushed the load to about 15 on 4 cores, and no tile finished in 30 min. Always export
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1` (scripts/CLAUDE.md, parallel topology). The `--maps` files from the first
  16 tiles were in an ephemeral scratchpad and are gone, so the morphology passes of t+0-1, t-1-1, t-1-2, t-1-3
  and t-2-3 (15) need those tiles re-run with `--maps` before their contact sheets can be made.
- **t+0+2 done (17 / 28; 456 s per tile at 2 jobs × `--procs 2`, OMP=1):** 51 hits, 10 without the bright-core flag; morphology
  edge 5, rosette 1, pass 4. Contact sheet looked at: h015 and h033 are +/− dipoles on one compact elongated galaxy
  core (registration/PSF residual), h028 a dipole on a compact core in SW only, h035 a weak SW-only spot on a knot of
  a clumpy galaxy (nothing in LW). None is single-signed and point-like across bands; not candidates (not formally
  vetted). Remaining o211→o212 tiles: t-2+0, t+1±, t+2± (11) plus `--maps` re-runs of the 7 map-less earlier
  tiles and the 5 with passes (≈ 70 min at this throughput: one run's worth). Suggest adding a dipole test
  (opposite-sign lobe within 2 × FWHM, about equal |S|) to `arc_diff_morph.py`; validate it on Yan+2023 first.
- **Run 01QcfDQb (10:09–10:50Z): 28 / 28 tiles.** t−2+0, t+1±0/±1/+2/−2 and t+2−1/−2/+0/+1/+2 done (`--procs 1`, 4 jobs, OMP=1;
  2–15 CPU-min per tile, so the 456 s figure was not a constant). 598 hits, 50 without the bright-core flag.
- **`arc_diff_morph.py` dipole verdict (new):** a hit is a dipole when an opposite-sign region (|S| ≥ 3, ≥ 4 px)
  within 0.3″ reaches ≥ 0.5 × |S_hit| in ≥ 2 bands (ASSUMPTIONs; unit tests). **Validation on Yan+2023 (Warhol box re-run with `--maps`):** the first version
  flagged 3 / 11 (Dc2-W4, D31-W7, D21-S2) from bands where |S_hit| < 3, so a nearby noise lobe gave a large
  ratio (failed approach); a band now counts only when |S_hit| ≥ 3, and 11 / 11 pass. Tile verdicts unchanged by
  the fix. On the eleven tiles: dipole 11, rosette 8, edge 13, extended 3, **pass 15**.
- **Contact sheets of the 15 passes looked at:** t−2+0 h029 is a positive core inside a negative ring on a bright compact source (PSF residual). **t+2+2 h002 (64.04865, −24.05870 tile; F150W +20.6, F200W +5.6)
  and h003 (F150W +15.0, F200W +7.3) are compact positive spots on blank sky in o211 with nothing in LW (|S| ≤ 3.5):
  unvetted, the open items of this unit.** SW-only epoch-2 spots without a source match the persistence class seen
  in Abell 2744 (#246), but that is a hypothesis until `persistence_check.py` runs on them. t+2+2 h009 is +5.2 in
  F277W and +5.8 in F444W but −1.3 in F356W (inconsistent). t+2+1 h063 t+2+1 h063 fades in SW (−6.0, −6.5; F200W dipole) on a knot beside a galaxy but rises in LW (+3.3, +2.6); t+2+0 h035 is an F150W dipole with nothing in F200W (+1.8) and F277W −1.9; h087, h098 (t+1±), h025 (t+1+2) and t+2−1 h009 change sign between
  bands (e.g. F150W −7.1, F277W +6.6); h092 is a compact source with an F150W dipole; h100 is a weak SW-only
  spot (+5.3) on a galaxy knot. t+1+2 h009, h014, h019 are positive SW residuals centred on compact galaxy cores:
  **in t+1+2 the S at the 10 brightest compact peaks is > +3 for 9 / 10 in F150W and 10 / 10 in F200W (median +12,
  +9; t+2+0 F150W 7 / 7)**, so these three are a tile-wide SW core-matching residual, not single brightenings (derived). Not
  candidates (not formally vetted).
- **`--maps` re-runs of t+0−1, t−1−1/−2/−3, t−2−3 (hits identical to the stored tables):** morphology with the dipole
  verdict leaves 10 passes there (was 15; 5 are now dipoles); their contact sheets are not yet looked at.
- **Failed approach:** none new; the tile cost varies a lot, so queue dense tiles first.
- **Handoff:** **first `persistence_check.py` and `arc_diff_vet.py` (CANUCS o004 as third epoch) on t+2+2 h002 and h003**, then
  h063 (t+2+1); then o211→o213 over the tiles, then the
  `--maps` re-runs of t+0−1, t−1−1/−2/−3, t−2−3 for their contact sheets, and the Warhol box with `--maps` to
  validate the dipole verdict on Yan+2023. Add a **core-population test** to `arc_diff_morph.py` (a hit's S
  against the S distribution at the tile's compact-source peaks in that band) and validate it on Yan+2023 first.
- **Run 01UyLcH7 (12:09–12:45Z): the open items vetted** (`results/transients/macs0416_arcdiff/vet_t22/`, derived):
  - Hypothesis before measuring: t+2+2 h002/h003 are o212 afterimages. Ordinary alternative: a real short
    transient, which should show in more than one dither. Plan change: an on-sky verdict sends it to `/vet-candidate`.
  - `persistence_check.py` on 396 `_cal` exposures of o211 and o212: h002 and h003 → **persistence**. There is one
    detection each in F150W and F200W, in dither 1 of 06101/08101 on NRCB3, after F115W 04101 exposures that put 28
    saturated pixels there. h003 is the grid's g40 pixel (2033, 1926). The `_cal` contact sheet shows the saturated
    core and the fading blob.
  - With CANUCS o004 (`arc_diff_vet.py`): h002 and h003 are present only in o212. o212→c004 F150W gives S −20.1 / −14.7,
    and o211→c004 gives |S| ≤ 1.9.
  - h063, re-differenced in a box centred on it: F200W −6.7 ± 1.6 nJy (7 % of the knot), F150W −2.5 ± 2.0. The LW
    change is arc-coherent. The sheet shows a −/+ pair 0.2″ apart in F200W only, which is a registration dipole.
    Its S depends on the box (F150W −6.0 in tile t+2+1, −3.1 centred). Not a candidate.
  - **Code review found a bug in the extent test (fixed):** `arc_diff_morph.py` rescaled the smoothed S map as if S
    were white noise. S is a matched-filter output, so the rescaled noise was 2.1–3.7σ, and 18–42 % of pixels
    passed the "3σ" streak level. The map is now normalised by its own robust spread, and a correlated-noise test
    was added. Only "extended" verdicts can change. The four tiles with short extents or a uniform 20″ extent
    (t+0−1, t+0−3, t+1−2, t−1−3) were re-run with `--maps`, and their morph JSONs were replaced:
    - t+0−1 h031/h041: now rosettes;
    - t+0−3 h064/h067 and t+1−2 h041: now pass, but the contact sheets put them on one straight positive trail that
      crosses all five bands in o212 (a satellite or asteroid trail). They pass only because the top band, F200W,
      shows the trail broken up;
    - t−1−3 h016: now passes, SW-only (+5.6 / +6.0). `persistence_check.py` → **persistence**. It appears only in
      dither 4 of each o212 filter, after dither 3 put a saturated star on those NRCB2 pixels, with weak LW too, as
      in Abell 2744 g4;
    - t+0−3's other 19 hits remain "extended" (a 19.9″ F150W trail);
    - t−2−3's 28 "extended" verdicts come from the old code. That tile's diffraction spike was already confirmed by
      eye.
  - **o211→o212 over the MACS0416 arcs: 0 candidates.** This is not a limit: the morphology veto is not
    injection-calibrated (D-027 addendum 2026-10-15f).
  - **Failed approach:** the extent test uses only the top band, so a trail that is clean in F150W but broken up in
    F200W passes. Next, take the extent over every band with |S_hit| ≥ 3, re-validated on Yan+2023.
- **Handoff:** first fix the open code-review findings listed in #245's body (multi-band extent, the |S_hit| guard for
  lobes and edge, the border case, per-band caching, `Params`), re-run t−2−3 morphology and re-validate on Yan+2023; then
  o211→o213 over the 28 tiles (`--maps`, OMP=1, 4 jobs `--procs 1`), with `persistence_check.py` on
  every single-epoch compact pass. Calibrate `arc_diff_morph.py` with injections through the whole chain. Add the core-population test.

## 2026-10-15: Abell 2744 two-epoch grid, UNCOVER o001 vs DD 2756 o002: 46 hit groups, all ordinary, 0 candidates
- **Hypothesis before computing:** a catalogued position in the Abell 2744 core changed by ≥ 0.3 mag at ≥ 5σ in two
  bands between UNCOVER (2022-11-02/04) and DD 2756 o002 (2022-12-06). This pair's first epoch is independent of the
  o003/o002 run. Ordinary explanations: spikes (PA 40.8° vs 59.5°), persistence, registration, edges. Plan change: any
  survivor goes to `/vet-candidate`.
- `i2d_region.py` cut UNCOVER to the o002 footprint (SW 373 MB and LW 90 MB per band written; the SW sources are 5.2–7.7 GB each).
  `transient_grid.py` was unchanged: 4,182 positions, noise scales 1.22–1.44. Injections: −0.75 mag 94.6 %,
  −0.44 mag 82.6 %.
- **46 hit groups, by contact sheets (6 bands × e1 / e2 resampled / difference):**
  - 24 star spikes; 21 of them come from one bright star at (3.5996, −30.4291);
  - 9 registration dipoles: the 0.5″ aperture changes by ≤ 0.06 mag in every band;
  - 4 edges and 3 background;
  - 3 persistence (`persistence_check.py`: g4, g7 on UNCOVER, g39 on o002; every detection suspect, none clean);
  - 2 detector-fixed artifacts: the o003/o002 run's g16/g17, recovered;
  - 1 processing artifact: g11, a star at F115W −0.82 mag only, from a suppressed core in the UNCOVER F115W mosaic.
  - **0 candidates.** `results/transients/a2744_grid_uncover/`.
- **Failed approach / lesson:** g4 looked like the best fading source: compact in all six UNCOVER bands and gone in
  o002. Every one of its detections is in dither 4, on pixels where the bright star's saturated core sat in dither 3.
  This is SW + LW persistence, so a multi-band detection does not rule out persistence. D-027 addendum 2026-10-15e.
- Data: 6 o002 `_i2d` (4.6 GB) plus the UNCOVER cuts, streamed and deleted.
- **Handoff:** the next Abell 2744 pair is UNCOVER F356W vs ALT 3516 F356W, single band (`--min-bands 1` plus
  `persistence_check.py`), validated on PLATE-23a once its position is checked. Or move the grid to another cluster
  with two NIRCam epochs.

## 2026-10-15: S0-2026-10-15b re-check: the proposed static radiation-fluid background does not exist; two over-strong claims withdrawn (derived)
- **Question (owner, ≤ 1 h):** is the proposed next background self-consistent, and is a solver worth building?
- **Background obstacle** (sympy-verified):
  - A static perfect fluid in −N²dt² + dl² + r²dΩ² obeys p′ = −(ρ + p)N′/N. With p = ρ/3 that gives ρ = C/N⁴.
  - With N finite and non-zero at both asymptotic ends (no horizon), ρ → 0 forces ρ ≡ 0.
  - It cannot be truncated either: without a surface layer, p must be continuous.
  - So no Ellis + static isotropic radiation background exists, and no solver will be built. This holds only for
    this closure (not for directed, injected or anisotropic radiation).
- **Withdrawn:**
  - The general spherical no-go: Birkhoff needs a vacuum gap, and Ellis's phantom scalar has monopole waves.
  - The claim that a first-order drive must be non-perturbative. h/Ω² is the mode's relative jitter, not
    necessarily the background's. Whether the background stays perturbative is undetermined until h is derived
    from a source, so a direct drive is not excluded.
- **Restated:**
  - The thin-shell toy is the static stability of already-trapped radiation, comparing different exotic-matter
    theories: the EoS intercept moves with f.
  - The 35–130 % and 10⁵² W figures are model estimates only.
  - "Not found" describes four queries only.
- **Paused:** a stability calculation is proposed only with a specified radiation distribution, stress tensor,
  boundary conditions, confinement or supply mechanism, and a fixed matter theory.

## 2026-10-15: S0-2026-10-15b Kapitza screen (partly superseded by the re-check above): an ordinary high-frequency drive can act on a wormhole throat only through its averaged stress, at about 10⁵² W or under an unproven confinement (derived, toy)
- **Question (owner):** can Kapitza-style dynamic stabilisation, driven by an explicit source, suppress known
  wormhole instabilities? Record: `docs/state/units/S0-2026-10-15b.md`. ≤ 2 h; X4/R1 paused.
- **Literature** (checked in the bodies):
  - Kapitza criterion h > √2 γΩ (1407.4803 §5; 2507.02736).
  - Floquet SYK (2404.08394) is a boundary-coupled AdS₂ model whose drive heats and closes the wormhole; it does
    not transfer to 4D.
  - Ellis radial mode β ≈ 1.18 c/r (0806.0608 Table I).
  - No prior paper stabilises a wormhole by driving or applies Isaacson averaging to one.
- **Derived:**
  - Spherical external drives cannot reach the l = 0 throat mode (shell theorem; 4D only).
  - Non-spherical fast drives act at second order, through the Isaacson / Brill–Hartle stress, which is positive
    and satisfies the NEC.
  - A first-order drive would need a stiffness modulation of about 20/a², 6–7 times the Ellis potential depth.
- **Toy (exact thin-shell junction conditions):**
  - Confined 2D radiation can stabilise the throat when the exotic EoS slope is < ½. For throats a₀ ≥ 4M it needs
    radiation energy 0.35–1.3 × |σ₀| and raises the exotic negative energy by 35–130 %; a₀ = 2.5M needs more.
  - It destabilises for slopes > ½.
  - Unconfined, it needs about 2–8 × 10⁵² W whatever the throat size.
- **Not done:** the field-supported (Ellis + radiation fluid) eigenvalue, about 3–4 h, which decides the sign for
  an explicit Lagrangian; S6a's 5D pair.
- **Verdict:** pause. A concrete shaker exists, but it needs more exotic matter plus 10⁵² W or a hypothetical
  confinement, and it does not make a communication channel.

## 2026-10-15: X4 S1f: the validated unlensed search finds no ordinary solution that makes R1 fail the screen; R1 stays undetermined (simulated)
- **Gate revision (owner, recorded before running):**
  - An exactly verified ordinary solution with ρ_req > 30 would settle R1 without proving the whole tool
    reliable.
  - Not finding one supports nothing.
- **Check:** the S1e unlensed procedure (validated on V6 and V7), seeds 21–23, on the S1c family-D R1 target, with
  the saved solutions re-scored.
  - Best FF 0.91351, at q 2.65 with a nearly anti-aligned |χ1| = 0.99 primary, viewed edge-on.
  - An independent explicit-ψ rebuild at pad 64 and pad 256 agrees to 1.2e-7.
  - ρ_req 7.4 against the 30 needed (the screen needs FF > 0.994987); intrinsic SNR needed 38 (≤ 100 holds).
  - The best ordinary fit of any kind is still 0.945 (lensed, S1c).
- **Result:** no counter-example, so R1 stays undetermined. No pilot proposed; S2 not started.
  - Missing: a validated positive-lens search, and a consistent unlensed optimum (seeds give 0.892 and 0.914).
  - 0.6 h used; X4 has about 7.2 h left.

## 2026-10-15: X4 S1e: joint lens + intrinsic search makes the unlensed fits reliable; one lensed check still fails, so R1 stays undetermined (simulated)
- **Question (owner):** with a joint search, does the tool pass its injection checks, and what is R1's verdict?
  Diagnostics and plan frozen first (`docs/state/units/X4.md`, S1e).
- **Record corrections:**
  - Bilby can run XPHM with phase marginalisation turned off. It was not used for task fit and cost; the
    10⁵–10⁶-call figure is an estimate.
  - "No existing tool" is narrowed to "no ready-to-use complete pipeline found in this round's check".
- **Diagnostics** (truth used only here):
  - The formal parameter path scores every truth at ≤ 3e-16.
  - The cached lens grid alone leaves ε 9e-4 – 1.3e-3 even with the true intrinsics, so it can only locate the
    region.
  - Exact continuous refinement from there reaches 6e-6.
- **Validation** (exact lens, pad 64):
  - U: ε ≤ 2.1e-6 on both seeds of V6 and V7.
  - P: V5 2.4e-5 / 3.7e-5; V9 1.5e-4 and **1.5e-3**. The failing seed's DE was still improving at its
    400-generation cap, and its refinement was still gaining at its 10-round cap.
  - The tool fails the frozen rule, so **R1 was not fitted** and stays undetermined. No physics conclusion.
- **Handoff:**
  - Minimal failing case: `scripts/x4_s1e_joint.py one V9 11`.
  - Continuing likely needs about 1,000–1,500 P generations (an estimate: about 4–5 h on 4 cores, at least about
    2.5–4 h on any number of cores without a code change) or advice on a precession parametrisation.
  - Owner's option: an R1 U fit testing outcome 1 only (about 0.5 h).
  - No pilot proposed; S2 not started. 1.8 h used; X4 has about 7.8 h left.

## 2026-10-15: X4 S1d: the reference global search (SciPy DE) fails its injection checks; R1 stays undetermined (simulated)
- **Question (owner):** does R1's negative-lens XPHM waveform keep a reliably identifiable difference, enough for a
  real-data pilot? Reuse check and plan frozen first (`docs/state/units/X4.md`, S1d).
- **Reuse:**
  - SciPy differential evolution is the reference search.
  - Polarisation and amplitude are maximised analytically, as in precessing searches; checked against a ψ scan.
  - Bilby/dynesty was not used: it does inference, not max-match, at an estimated (not measured) 10⁵–10⁶
    calls. It can run XPHM with phase marginalisation turned off (wording corrected in S1e).
  - PyCBC has no max-match optimiser for this objective (and is GPLv3). GLoW computes lens factors, which we already
    have exactly.
- **Result:** the tool fails the frozen rule (ε ≤ 5e-4 on both seeds).
  - Unlensed: ε 5e-6 – 1.5e-3, much better than S1c's 0.026, but every run hit its generation cap.
  - Lensed: ε 0.045–0.091 in 3 of 4 runs. The lens-grid scan from unlensed intrinsics picks a weak lens, and the
    node-repeat stop ends the search there.
  - So the R1 case was not fitted. R1 stays **undetermined**; this is not evidence for or against the
    negative-mass model.
- **Handoff:** no pilot proposed; S2 not started. Remedy (owner's call, about 4–4.5 h): a joint lens + intrinsic DE
  on the cached lens grid, and longer unlensed runs. About 1.3 h used; X4 has about 9.6 h left.

## 2026-10-15: Abell 2744 two-epoch grid (DD 2756 o003/o002): 22 hit groups, all ordinary, 0 candidates
- **Hypothesis before computing:** a catalogued position in the Abell 2744 core changed by ≥ 0.3 mag at ≥ 5σ in two
  bands between 2022-10-20 and 2022-12-06 (47 d; e.g. the DD 2756 lensed SN or a caustic transient). Ordinary:
  spikes, persistence, proper motion, edges, background. Plan change: any survivor goes to `/vet-candidate`.
- **Pair choice (failed approach first):** UNCOVER o001 (2022-11) vs o006 (2024-07) looked ideal by MAST centre
  coordinates, but o006's mosaic covers Dec −30.43 to −30.55 and misses the core. The ALT and MegaScience epochs share
  only F356W or F090W with earlier data. The DD 2756 visits have six shared bands over the core.
- `transient_grid.py`, unchanged (6 bands, 3,431 positions, noise scales 1.18–1.34). Injections: −0.75 mag 99 %,
  −0.44 mag 89 %. **22 hit groups:** 12 star spikes/wings (the PA changes 30°→60°), 1 star wing (g13, 20.5 mag
  point source, no proper motion beyond the 13/22 mas frame offset), 5 persistence, 2 detector-fixed artifact, 1 galaxy-edge background, 1
  mosaic edge. **0 candidates.**
  `results/transients/a2744_grid_2756/` (`hit_verdicts.ecsv`).
- **Persistence:** g0/g3/g4/g6/g7 are diffuse SW-only blobs in epoch 1. They are absent in the simultaneous LW
  exposures, in epoch 2 and in UNCOVER 13 days later. `persistence_check.py` gives "persistence" for all five.
  g16/g17 (epoch 2, SW only): the tool gives on_sky / inconclusive, but every detection of both sits on one NRCB3 pixel
  (~1855, 1976), in dithers 2 and 1, with S/N decaying (F115W 64 → 34, F150W 11, F200W 5). That makes them a
  detector-fixed feature at two sky positions (by hand). Gap: the tool checks each position alone.
- **Contact sheets looked at:** six bands × (e1, e2 resampled onto e1, difference), 2″ stamps; UNCOVER F150W, F200W
  and F277W stamps at the eight artifact positions. First sheets without resampling were misleading: the 30° roll
  made stars look displaced.
- **Validation limit:** there is no published Abell 2744 transient with coordinates in this window (the DD 2756 SN
  has no public position), so recovery rests on injections and on the grid's MACS0416 validation (D-027 addendum).
- New: `scripts/i2d_region.py` (byte-range SCI/ERR cut of a huge mosaic to another mosaic's footprint; unit tests),
  for an UNCOVER pair. Data: 12 `_i2d` (9.6 GB) streamed and deleted.
- **Handoff:** UNCOVER o001 (2022-11-02) vs DD 2756 o002 (2022-12-06), with `i2d_region.py` cutting UNCOVER SW to
  the 2756 footprint. UNCOVER F356W vs ALT 3516 F356W (2023-11/12) is single-band, so it needs `--min-bands 1` plus
  the persistence test; PLATE-23a (Pang+2026) would validate it once its position is checked.

## 2026-10-15: X4 S1c: R1 keeps its residual against validated aligned-spin fits; undetermined with IMRPhenomXPHM (simulated)
- **Question (owner):** does R1's large difference (30 M☉ behind a 300 M☉ negative lens, y 0.3) survive a validated
  ordinary fit and a more complete waveform? Plan frozen first (`docs/state/units/X4.md`, S1c). 2.0 h, the cap.
- **Tool:**
  - One fitter serves validation and the case fit: Nelder–Mead restarted to convergence, plus a frozen 1,271-node
    lens-grid scan.
  - The FFT length is the next fast length (FF change ≤ 1.6e-7, 5.6× faster). A batched-FFT grid scan was measured
    slower and not used.
  - S1b's failed injection is now recovered (ε 3.5e-6); every aligned-spin check passes.
  - The XPHM unlensed check fails (ε 0.026 > 0.0125). The XPHM lensed check passes from one start only.
- **Result:**
  - Aligned-spin, dominant-mode templates: the residual is kept, ρ_req 5.8 (observed SNR), intrinsic SNR 30. It is
    bound-conditional (q 15, |χ| 0.99) and provisional: the frozen saved-solution starts were cut for budget. The best mimic is a strong positive lens (2,400 M☉, y 0.1).
  - IMRPhenomXPHM: **undetermined**. The best ordinary fit found is FF 0.945 (ρ_req 9.2), from precession with
    near-maximal spins plus a positive lens (1,400 M☉, y 0.22). The starts disagree and the fits did not converge.
- **Disclosed:**
  - family-A saved starts were cut for budget;
  - family D used only 4 of the 25 saved starts, and two of those fits were stopped at the cap;
  - the results file was post-processed after the stop (label, grid maps);
  - the convergence flag excludes the intrinsic restart status.
- **Wording tightened:** "absorbed" means only ρ_req > 30; the difference is smaller, not gone. It says nothing
  about the 16 untested S1 cases.
- **Handoff:** no real-data pilot proposed; S2 not started. Next (owner's call): a validated precessing fitter for
  R1.

## 2026-10-15: MACS0416 arc-difference hit h108 is an arc-segment residual; arc-coherence vetting test
- **Hypothesis before measuring:** h108 (o211→o213, #242) is a point source on a faint thin arc that brightened. It
  would be one PSF-sized spot, with the same sign in every band, and additive over the epochs. A matching residual
  would also appear on the neighbouring knots of the arc.
- There is no jw01176 o214; the #242 handoff meant o213. **CANUCS jw01208 o004 (2023-01-11)** is an independent fourth
  epoch between o212 (2022-12-29) and o213 (2023-02-10), covering h108 with 8 exposures. Five pairs differenced.
- `scripts/arc_diff_vet.py` (new; unit tests on synthetic arcs): forced aperture fluxes in nJy with blank-sky and
  **arc noise** (apertures of similar epoch-1 arc flux), S at the hit, and **arc coherence**: the fraction of
  neighbouring knots (0.25–1.0″) with S of the same sign and ≥ 0.5 × |S_hit| (ASSUMPTIONs).
- **Result (derived):** the fluxes are additive (F356W o211→o212 + o212→o213 = +5.6 vs +5.0 nJy direct); the largest
  change is o211→c004, about 13 nJy (28.6 AB) in F277W–F444W. But the arc's knots follow h108 in 4 of 5 pairs, including both
  pairs without o211: c004→o213 F356W knots −5.5, −5.3, −4.2 (h108 −3.7). In o212→o213, F277W rises (+4.2 to +6.6)
  while F356W falls (−3.8 to −5.2). **Verdict: artifact (arc-segment PSF/background matching residual); not a
  candidate.** Record: `docs/candidates/macs0416-arcdiff-o211-o213-h108.md`; JSON in
  `results/transients/macs0416_arcdiff/vet_h108/`.
- **Validation:** the Yan+2023 arc transients are never flagged arc-coherent in their discovery pair (0 / 9
  testable); flags appear only in the other pair, where the transient is weak and the Warhol arc carries other
  transients. Caveat: coherence can veto a weak transient on an arc with several real ones.
- **Contact sheets looked at:** 1.5″ and 0.6″ stamps, five pairs × five bands.
- **Failed approach:** the first per-band "max |S| within 0.2″" numbers (mixed signs) were misleading; S at the exact
  position, and aperture fluxes against arc noise, are the measurements to use.
- **Handoff:** tile the rest of MACS0416's arcs with `arc_diff.py` (adopted settings), add CANUCS o004 as a fourth
  epoch, and run `arc_diff_vet.py` on every unflagged hit. A per-segment local kernel may remove the arc residual.

## 2026-10-15: MACS0416 arc differencing: fitted pixel kernels remove the arc ridges and the compact-core residual
- **Hypothesis before computing:** the arc ridges left by #241 are a PSF-matching failure that a least-squares pixel
  kernel (Bramich 2008 style) removes without losing the published arc transients. Plan change: if recovery of the
  Yan+2023 transients drops, the kernel is not adopted.
- `scripts/arc_diff.py`: `--kernel-px 4` (9 × 9 kernel + constant, fitted on the whole box, 4σ clipping),
  `--kernel-stamps` (second kernel for compact-source stamps, 97th-percentile peaks, leave one stamp out),
  `near_bright_core` flag (0.8″, 99.9th percentile; ASSUMPTIONs). Unit tests on synthetic scenes.
- **Results (20″ box, five bands; derived):**

  | PSF match | Ep1→Ep2 published hits | Ep1→Ep3 (pair 31) | compact-peak residual, SW | unflagged hits o212 / o213 |
  |---|---|---|---|---|
  | Gaussian + scale (#241) | 4 / 5 | — | +20–43 % | — |
  | whole-box kernel | 4 / 5 | 3 / 4 | +11–28 % | 30 / 28 |
  | stamp-only kernel | 1 / 5 | 2 / 4 | ≈ 1 % | 2 / 9 |
  | stamps ∪ arc pixels, one kernel | 2 / 5 | 2 / 4 | +2–10 % | 2 / 9 |
  | **whole-box + stamp kernels (adopted, after review fixes)** | **4 / 5** | **2 / 4** | **≤ 3 %** | **4 / 8** |

  Adopted run: injections on bright pixels 28.0 mag 0.79–1.0, 28.5 mag 0.38–0.96 (F444W lowest); the bright-core
  flag costs about 30 % of them (recovered and unflagged in that band: 0.64–0.73 at 27.5–28.0 mag). Missed: D21-W3,
  D31-W5 (≤ 4.8σ in every band) and, after the review fixes, D31-S4 (it was a hit before them). D23-S3 is a hit in
  o211/o213 (negative). D21-W1 is negative in o211/o213 (brighter in epoch 1 than in epoch 3; observed).
  `results/transients/macs0416_arcdiff/summary_o211_o21{2,3}_kernel.json`, `hits_*_kernel.ecsv`.
- **Contact sheets looked at:** the whole-box kernel removes the Warhol and Spock arc ridges in every band; the
  residual it leaves is positive cores of compact sources in F150W/F200W (stars, knots, galaxy cores), which a
  stamp kernel fixes. All 4 unflagged o211/o212 hits are published transients. Of the 3 unexplained o211/o213
  hits, 2 are a rosette wing (h080) and the box edge (h109); **h108** (64.03682, −24.06626; F200W/F277W, 5.7σ) sits
  on the arc noted in #241 with no visible point source. **Unvetted; not a candidate.**
- **Failed approaches:** a stamp-only kernel (it fits cores but not the arc light, and loses the arc transients);
  one kernel on stamps ∪ arc pixels (a compromise worse on both); leave-one-stamp-out alone changed nothing for
  the published transients (they are not in stamps), so self-absorption was not the cause.
- Review fixes (code-review): leave-one-out fits with too few rows fall back to the whole-box kernel; stamp
  pixels are never differenced with a kernel trained on themselves; the final clip is refitted; stamp-region
  variances use the stamp kernel.
- **Not covered:** PSF rotation rosettes around bright stars (flagged, not removed); only one 20″ box.
- **Handoff:** vet o211/o213 h108 (third epoch o214, cutouts per band); then tile the rest of MACS0416's arcs with
  the adopted settings.

## 2026-10-15: X4 S1b: against precessing, higher-mode binaries the representative passing cases fall below the screening threshold (simulated)
- **Question (owner):** after the ordinary explanations get a fuller fitting chance, is there still a representative
  case worth taking to real data? Plan and criteria frozen first (`docs/state/units/X4.md`, S1b).
- **S1 status corrected:** "18 pass" is a preliminary result of the searches run so far, not a proven optimum.
- **Numerical stage** (4 exploratory cases chosen from S1; multi-start Nelder–Mead with restart convergence;
  fit-error injections):
  - R2 and R3 keep a residual against aligned-spin dominant-mode models. Strictly, they are undetermined: the
    best P block descent was not converged (review finding). Step D was run on them anyway, as a disclosed
    deviation.
  - R1 (the large-difference case) is undetermined: the positive-lens fit fails its own injection test
    (ε 0.13).
  - The R4 control is absorbed.
- **IMRPhenomXPHM** (target and templates; no free overall phase): R2 (ρ_req 30.8) and R3 (35.8) are absorbed, which here means only ρ_req > 30: the difference is smaller, not gone.
  - R2 through an edge-on inclination (higher modes); R3 through in-plane spins at inclination 49–84°. q stays
    1.2–1.7.
  - R3 was already absorbed in the frozen run (31.4).
  - R2 passed the frozen run (28.9) and is absorbed only after post-hoc extra starts. More search can only raise
    ρ_req.
- **Post-hoc, labelled:** convergence continuations and 8 extra R2 starts, after the frozen runs hit their
  evaluation caps.
- **Handoff:** no real-data pilot proposed; S2 not started. Open: R1-type cases need a P optimiser that passes its
  injection test, then XPHM (owner's call). Not covered: one target orientation, real noise, priors, other lenses.
## 2026-10-15: MACS0416 arc difference imaging (o211/o212): tool validated on compact arc transients; arc ridges remain
- **Hypothesis before computing:** two-epoch difference images reach the caustic-crossing stars on the Warhol and
  Spock arcs that the catalog grid cannot (D-027 addendum). Ordinary: misregistration, PSF and zero-point differences,
  star spikes, persistence; known transients (Yan et al. 2023). Plan change: if the published arc transients are not
  recovered, fix the tool before any search.
- New `scripts/arc_diff.py`: S3 byte-range box (both epochs, SCI + ERR), resample, global + tiled registration
  (3″), Gaussian blur of the sharper epoch and a flux scale fitted on the brightest 10 % of pixels, 0.6″ median
  background, matched-filter S ≥ 5 in ≥ 2 bands (ASSUMPTIONs). Unit tests on synthetic images.
- **Validation (20″ box, F150W/F200W/F277W/F356W/F444W):** 4 of 5 published Ep1→Ep2 arc transients are hits
  (D21-W1 13σ, W2 6σ, S1 9σ, S2 7σ); W3 peaks at 4.9σ. Injected point sources on bright pixels: 28.0 mag 0.89–1.0,
  28.5 mag 0.58–0.95 (LW lowest). Flux scale epoch 2/epoch 1 = 0.977–0.986 in every band (derived).
  `results/transients/macs0416_arcdiff/summary_o211_o212.json`.
- **Failed approaches (contact sheets looked at, three iterations):** global shift + global-rms blur left dipole and
  positive ridges along every thin arc (blurring an epoch also lowers its noise, so the global rms over-blurs);
  tiled registration alone did not remove them; scoring the blur on bright pixels cleaned the Warhol region in LW,
  but the Spock arc and the arc near (64.0369, −24.0664) still show a positive ridge in F200W and F356W.
  D31-S4 at 6.8σ in this Ep1→Ep2 difference sits on that ridge, so it is not a detection.
- **No search result:** 149 hits, unvetted, mostly arc ridges and bright-star residuals. No candidate claimed.
- Data: ten 20″ boxes streamed by byte range (no mosaic stored).
- **Handoff:** fix the arc ridge (fitted PSF-matching kernel, e.g. a pixel kernel by least squares on stars, or a
  local flux scale along arcs), re-validate on the five Ep1→Ep2 transients, then o211/o213 (W5–W7, S4).

## 2026-10-15: MACS0416 grid, epoch 1 against epoch 3 (o211/o213): no new transient
- **Hypothesis before computing:** a two-band change between 2022-10-07 and 2023-02-10 that the published search
  (Yan et al. 2023, F356W-led) did not report. Ordinary: as in the o211/o212 entry below. Plan change: any survivor
  goes to `/vet-candidate`.
- Same command and defaults as #238: 8,345 positions (197 merged after refinement), 3,470 measured in all four
  bands, 104 hits in 50 groups: 47 bright-star spikes, wings or footprint edges, 2 SW-only persistence
  afterimages (`persistence_check.py`: every detection suspect), SN01 (`on_sky`, 26 clean exposures; mixed sign).
  **0 new candidates.** Injection −0.75 mag 95 %, −0.44 mag 40 %. Contact sheets looked at.
  `results/transients/macs0416_grid_o213/`.
- **Blind spot, larger here:** neither the o211 nor the o213 catalogs have an entry within 3.2″ of the Warhol arc,
  so its five published transients (D21-W1/W2, D31-W5/W6/W7) have no grid position; in o211/o212 only the o212
  catalogs supplied them. Arc transients need image differencing, not the catalog grid.
- Data: o211 and o213 mosaics (6.4 GB) streamed and deleted.
- **Handoff:** o212/o213 is the last PEARLS pair (expect the same limits); the higher-value step is PSF-matched
  difference imaging of the Warhol and Spock arcs, validated on the 12 published arc transients.

## 2026-10-15: MACS0416 two-epoch grid, validated on published transients: two grid fixes, no new transient
- **Hypothesis before computing:** the grid recovers the published PEARLS epoch-1→2 transients (Yan et al. 2023:
  6 caustic stars, SN01, SN02); any further two-band change in o211/o212 is a candidate. Ordinary: spikes,
  persistence, SNe and caustic stars (the expected populations). Plan change if the known cases are missed.
- **Validation failed first: 0 of 8 recovered.** Causes (D-027 addendum 2026-10-15): a 0.20″ frame offset of the
  published positions; per-band recentring pulled SN01's LW apertures onto a neighbour; SN01 changed colour
  (fades in F115W, brightens in F356W), which the same-sign rule rejects.
- New: `transient_grid.py --recentre-in origin` (default; centroid once where the catalog saw the source) and
  mixed-sign hits (column `mixed`), `--catalog-bands`, re-merge after refinement, with unit tests. SN01 recovered; the new rules add only that row.
- **Result:** 174 hits in 59 groups: 50 star spikes and wings, 6 edges, 2 persistence (`persistence_check.py`),
  SN01. **0 new candidates.** Injection −0.75 mag 95 %, −0.44 mag 46 %. Contact sheets looked at.
  `results/transients/macs0416_grid/`.
- **Failed approach / limit:** the caustic stars (F356W 27.1–28.6 on arcs) are below the grid's depth (|σ| ≤ 3.3),
  and SN02 has no catalog position of its own. The grid cannot test transients on arcs; image differencing can.
- Data: 8 `_i2d` mosaics (6.0 GB) streamed and deleted.
- **Handoff:** Ep1/Ep3 (o213) and Ep2/Ep3 with the same command; arc-region image differencing (PSF-matched) is
  the step that would reach the caustic stars; re-run Sunrise/SMACS under the new defaults (cheap, ~10 s each,
  needs the mosaics again).

## 2026-10-15: X4 S1: free aligned-spin binaries mostly do not absorb a negative-mass lens's GW signature (simulated)
- **Question:** with masses, spins, amplitude, time and phase free (IMRPhenomXAS, aLIGO design curve LIGO-P1200087),
  can an unlensed binary or a positive point lens fit a binary lensed by a negative point mass? Configuration
  and screening threshold were frozen first (`docs/state/units/X4.md`, S1).
- **Tool:** validated (FF 0.9998–0.99995 for matching models; positive-lens optimiser floor 0.9963). Exact lens
  factors with a tested geometric-optics fallback at very large w, cached on disk.
- **Result (36 cases, all converged):**
  - Main run: 19 pass the screen. 17 fail: 13 are too demagnified to be seen (intrinsic SNR ≥ 124), and 4 are
    absorbed by a positive lens (200 M☉ binary behind 100 M☉, or behind 300 M☉ at y = 1.5).
  - The frozen 6-parameter positive-lens refinement keeps all 19.
  - With mass ratio up to 15, 18 of 19 pass, for 30–60 M☉ binaries behind 100–1,000 M☉ lenses and 200 M☉ behind
    300–1,000 M☉, using the best fit over all stages. Best ordinary FF 0.63–0.98 (unlensed) and 0.75–0.99
    (positive lens); observed SNR needed 4.5–20, intrinsic 20–95. The result is bound-conditional: fits run to
    q = 15.
  - Extreme mass ratios absorb part of the signature (fits run to q = 15), so higher modes and precession are
    needed next.
  - IMRPhenomD fits agree within 0.0065.
- **Review fixes:**
  - The first IMRPhenomD cross-check had silently used IMRPhenomXAS; it was rerun.
  - The frozen refinement was run.
  - Bounds and time-resolution effects are disclosed.
- **Not covered:** precession, higher modes, eccentricity, other lens models (SIS, embedded, binary), a detector
  network, real noise, priors and selection. Screening threshold only, not a detection threshold.
- **Handoff:** S1 gate met, conditionally (aligned spins, dominant mode). The S2a pilot on GW231123 (≤ 3 h, after a
  higher-mode and precessing check) is proposed; it waits for the owner.

## 2026-10-15: Sunrise two-epoch grid on every catalogued position: no transient (D-027 addendum)
- **Hypothesis before computing:** a source in the o010/o120 overlap changed by ≥ 0.3 mag at ≥ 5σ in two bands, but
  the catalog stage and its two-band rule missed it (e.g. a blue transient seen only in F090W and F115W).
  Ordinary: spikes, persistence, moving objects, background structure. Plan change: any survivor goes to
  `/vet-candidate`.
- New: `scripts/transient_grid.py` (whole-mosaic forced photometry on the union of both epochs' catalogs, noise
  scale from all positions, injections with the real zero point) and unit tests; `results/transients/sunrise_grid/`.
- **Result:** 6,401 positions (2,331 measured in all four bands), 66 hits in 36 groups: 34 star spikes, streaks or
  wings, and the known persistence afterimages `n0022` and `n0150` (recovered with no catalog-stage input).
  **0 candidates.** Injection at nominal amplitude: −0.75 mag 100 % (198/198), −0.44 mag 53 % (before by-eye
  vetting). Contact sheets looked at.
- **Failed approaches:** a difference of two negative fluxes (annulus brighter than the aperture in an extended
  galaxy) passed as appeared/disappeared; now one epoch needs positive flux. A zero point recomputed on the injected positions absorbed the injected brightening (5 %
  recovery); fixed by using the real grid's zero point and scale.
- Data: eight `_i2d` mosaics (6.4 GB) streamed into the gitignored data root and deleted after the run.
- **Handoff:** the same grid for Sunrise o010/o052 and SMACS/VENUS (F150W, F444W), with a spike veto first.

## 2026-10-15: S0 screen after X3: GW wave-optics lensing by negative-mass objects qualifies (X4)
- Three directions (`docs/state/units/S0-2026-10-15.md`):
  - A, GW lensing by repulsive lenses: **qualifies**.
  - B, VASCO vanishing sources: not qualified; no mechanism links it to the hypothesis.
  - C, the quantum-processor wormhole: not qualified; it is a simulation of a toy model.
- **Minimal step (A; derived and simulated, no strain):**
  - Closed-form point-lens amplification factor for either sign, checked against the diffraction integral
    (≤ 3e-7 at the tested points) and against the W3 magnification (< 3 %).
  - A negative-mass lens is a low-pass filter on GWs.
  - With intrinsic parameters fixed, it separates from the best positive point mass at observed SNR
    ρ ≈ 3.1–16.5 for M_Lz ≥ 300 M☉ and y ≤ 1 (4.8–46.5 at y = 1.5), but not for 30 M☉ lenses.
  - The separable cases are demagnified, so their detectable volume falls by ≳ 100×.
  - Review fixes: exact residual ρ²(1 − O²); the Ellis wormhole is a separate model; Kobrin, Schuster & Yao
    cited for C; X4 gates sharpened.
  - Real-data hook: GW231123, the GWTC-4.0 lensing outlier, has been fitted only with positive lenses.
- **Handoff:** X4 S1 (intrinsic-parameter degeneracy with an IMR waveform, ≤ 4 h) is proposed, not started.

## 2026-10-14: MACS1149 system 16: z 2.5–2.9 from the pair; a possible far image (CANUCS 5106117)
- **Hypothesis before computing:** system 16's far image is detectable along its redshift track and fixes z, or
  no z keeps every predicted image consistent in both models. Ordinary: a faint far image, chance field sources.
  Plan change if no z in 2–6 is allowed in both models.
- New: `scripts/m1149_sys16.py` (CATS v4.1 and Sharon v4cor z scan, pair rms, forced photometry of every extra
  image in three bands, colour matches at one offset) and unit tests; output
  `results/lens_consistency/macs1149-cats/sys16_track.json`.
- **Result** (D-037 addendum 2026-10-14, docs/fields/macs1149.md): the pair is reproduced at z 2.10–2.90 (CATS) and
  2.50–3.30 (Sharon); joint 2.5–2.9, overlapping the CANUCS photo-z of 16.2. The far image is predicted ≤ 2.5σ
  there (point-source μ). One colour-matched source near both models' far image at z ≈ 2.5–2.65, CANUCS 5106117
  (z_phot 2.94, 2.49–3.16), is a possible third image; the three other matches have photo-z 1.23, 1.68 and 4.59.
  Stamps looked at. **Null:** no system 16 image is unmet.
- **Failed approach:** the far-image track as a redshift ruler; the point-source |μ| of the fold pair (≥ 39) makes
  the predicted far image ≤ 2.5σ, so photometry did not constrain z. First-pass colours measured each band at its
  own best offset (caught in review; now one F277W offset for all bands, which also found S1).
- **Handoff:** MACS1149 system 16 is settled up to spectroscopy (16.1/16.2 or 5106117). Next field items in TASKS:
  the Sunrise rel_weight plateau, duplicate merged sources (Abell 2744 `3999`/`7179`).

## 2026-10-14: MACS1149's 11 untested counter-images: none unmet in both models
- **Hypothesis before computing:** one of the 11 CATS forced rows without a reference flux (systems 13, 14, nine SN
  Refsdal host knots) is a predicted image both models need and the data do not show. Ordinary: model copies of
  catalogued images, resolved host light, images below the noise. Plan change if a gated row is `absent` and Sharon
  v4cor predicts it too.
- New: `scripts/m1149_refs.py` (sibling references without the |μ| cap in F150W/F277W/F444W; Sharon v4cor solved per
  system, D-041 cross-model rule) and unit tests; tracked input and output in `results/lens_consistency/macs1149-cats/`.
- **Result** (D-037 addendum 2026-10-14, docs/fields/macs1149.md): 0 `absent` in any band, gated or ungated. 13 is
  CATS's misplaced 13.2; 1002 and 1222 are CATS-only extras; 1007 is both-model with μ 18.3 vs 0.1 (untestable); 14
  is undetectable (2.4σ); six both-model rows lie in the resolved host images 1.2/1.3 (ungated `recovered` ratio
  1.4–1.9, or untestable). Stamps looked at.
  **Null:** no MACS1149 counter-image is unmet in both models.
- **Failed approach (caught in review):** reading "Sharon's nearest image lies on a catalogued image" as "CATS
  misplaced it" without checking that CATS reproduces that image; fixed by splitting Sharon's images first.
- **Handoff:** no MACS1149 counter-image is unmet in both models (six rows untestable, no injection limit). Next field items in TASKS: system 16's far-image
  track (fixes its z), the Sunrise rel_weight plateau, duplicate merged sources.

## 2026-10-14: X3 G2 prerequisite: astrometry mostly repeats the blend effect; narrow lens-shift gain
- **Question before computing:** can Gaia epoch astrometry plus photometry tell a W3 negative-mass lens from an
  ordinary blended variable, with the same blend freedom in both? The criterion was committed first:
  Δχ² ≥ 25 measurable, 9–25 marginal, < 9 not measurable.
- **Analytic result (`docs/state/units/X3.md`, G2 prerequisite):** the blend-ratio photocentre shift (∝ 1/F along
  a fixed line) is reproduced exactly by any one-variable blend (residual ~1e-24). Only the lens's own
  image-centroid shift is new.
- **Cases (`scripts/x3_g2_ident.py`, `results/x3/g2_ident.json`; simulated):**
  - The image shift is measurable against that ordinary model for G ≲ 17, θ_E ≳ 3 mas and t_E ≈ 100 d, under the
    assumed per-transit σ (doubling σ divides Δχ² by 4). Base
    case Δχ² 70; G 15: 486; G 19: 10; θ_E 1 mas: 8; t_E 30 or 300 d: 6–7.
  - It adds sign-of-mass information at most at G 15 or θ_E 10 mas (Δχ² ≤ 42, ≤ 67; upper bounds, because the
    positive-mass comparison keeps W3's t₀, t_E and u₀).
  - Against two variable blended stars, astrometry alone adds nothing.
- **DR4 (ESA DR4 content page):** epoch astrometry for all sources in the main astrometric processing. Whether
  it carries AL coordinates, scan angles and uncertainties, and its faint-source coverage, are unknown.
- **G1 corrections (owner):**
  - the 196 flags were excluded by rules, not individually proven ordinary;
  - the 92–99.6 % RUWE figure and the blended-prior limit are model estimates;
  - the claim that DR4 astrometry would make the photocentre motion a signature is withdrawn.
- **Handoff:** no DR4 development until ESA documents the epoch-astrometry content. If it fits (AL, scan angle,
  σ for G ≲ 17), build a joint photometry + image-shift likelihood on dev injections against the ordinary model.

## 2026-10-14: X3 G1, W3 screen of held-out Gaia GAPS: 196 flags, no candidate, conditional limit
- **Hypothesis before computing:** negative-mass lenses would show the W3 dip in unselected Gaia light curves.
  Ordinary: eclipsing binaries, dust dimmers, LPVs, single-scan photometry, crowding. Plan change if any flag
  survives ordinary vetting with W3 morphology.
- **Counts:** the dev split `MOD(source_id, 10) = 0` is 20 %, not 10 %, because GAPS ids are even. At G ≤ 19 the
  held-out set is 423,598 sources. The split is kept as frozen.
- **New scripts:**
  - `x3_g1.py`: full screen with a repeat-dip veto (no period search) and a per-band blend colour test; ECL, YSO
    and LPV references; selection bound.
  - `x3_bulk.py`: streams the bulk epoch-photometry files; verified identical to DataLink on 309 dev sources.
  - `x3_vet.py`: batch vetting.
- **Continue conditions** were frozen before the injection and reference results were read; the dev real count and
  the bulk check were already known. All were met. Results on dev: W3 efficiency 0.32 (t_E 30–300 d), PSPL
  0 / 1,995, dev real flags 0 / 1,995.
- **Limits of the screen (measured):**
  - The colour test removes 39 % of W3 and 48 % of dust dips.
  - Blended W3 events inflate their own RUWE, so an estimated 92–99.6 % would be vetoed by v1. This is a model
    estimate under assumed blend and error priors, not a measured Gaia rate or a strict bound (owner correction).
- **Held-out run (once):**
  - 395,764 exposed sources, 196 final flags.
  - The flags were excluded by rules; no ordinary cause was proven for each flag individually (owner correction).
    The rules were post hoc, and every rule was costed on the dev injections. First failed rule: 124
    single-visit, 8 bright, 29 neighbour, 27 known variables, 5 degenerate fits, 3 with no multi-visit dip.
  - **No candidate.**
  - Vetted efficiency 0.276 [0.253, 0.300].
  - Γ₉₅ < 9.7 × 10⁻⁶ per star per year for photocentre-stable W3 (t_E 30–300 d); 4.6 × 10⁻⁴ under the blended
    prior (a model estimate with the same label as the RUWE estimate).
- **Review fixes before merge:**
  - D2 coded and costed;
  - host rules costed;
  - v1 no longer counted twice;
  - t_E-matched retained fraction;
  - SIMBAD otype match;
  - id-keyed caches.
- **Failed approach:** the global `vari_classifier_result` reference query timed out (HTTP 408); per-HEALPix
  source_id ranges worked. The M31 field has too few dev YSOs, so Orion was used for them (different cadence).
- **Handoff:** X3 G2 waits for DR4, and its design on dev (epoch astrometry, joint colour fit) comes first.

## 2026-10-14: the El Gordo and Sunrise follow-up sources are ordinary
- **Hypothesis before computing:** one of El Gordo `2828`/`1601` or Sunrise `1756`/`1907`/`1578` could be a lensed
  image with an unmet model prediction, or a source no ordinary class explains. Ordinary: spikes, persistence,
  catalogue effects, dusty disks, stars, emission-line galaxies. Plan change if a counter-image is `absent` with
  injection `recovered`, or a source survives every ordinary test.
- New: `scripts/field_followups.py` (DJA matched photometry, spike geometry on the mosaic axes, five-band stamp
  photometry, Caminha+23 counter-images with D-029 forced photometry, RELICS magnifications and counter-images,
  two-epoch proper motion) and unit tests; `results/field_followups/summary.json`.
- **Result** (details: docs/fields/elgordo.md, sunrise.md "Follow-up sources"):
  - El Gordo 2828 is the horizontal diffraction spike of Gaia DR3 4931899405421624704 (G 15.9), on its mosaic
    row 6.0″ away. An artefact.
  - El Gordo 1601 is a real faint source whose SW bands are ~2 mag brighter than its LW bands (an emission-line
    galaxy is a hypothesis). Its eight Caminha+23 counter-images at z 2–6 have no `absent` verdict.
  - Sunrise 1578 is an edge-on dusty disk; 1907 is a dusty clumpy disk (μ ≤ 6, singly imaged).
  - Sunrise 1756 is a point-like source on the sky in 46 / 48 exposures, steady over 164 d, with no detected
    proper motion (≤ 48 mas/yr, 2σ) and a stellar-like SED. As a z ≥ 2 source, 2–3 of the 4 RELICS models at
    each z predict counter-images that forced photometry with injection finds `absent`. So it is read as
    foreground (likely a star). DJA misses it (undetermined; a stripe 0.6″ away is not star 1817's spike).
  - **Null:** no unmet lens-model prediction, no candidate.
- **Failed approach:** the DJA v7.0/v7.5 `eazypy.zout.fits` URLs return 403 (they ship only in tarballs), so the
  lens checks scan z_s ∈ {1, 2, 3, 6} instead of using photo-z.
- **Handoff:** the field follow-up lists of Abell 2744, El Gordo and Sunrise are closed. Next from TASKS
  "Field follow-ups": Sunrise `1869` (spike), the rel_weight plateau, duplicate merged sources.

## 2026-10-14: the five Abell 2744 follow-up sources are ordinary
- **Hypothesis before computing:** one of the five unvetted M1 follow-ups (docs/fields/abell2744.md) could be a
  lensed image no published model explains, or a configuration around an invisible deflector. Ordinary: known
  compact populations (little red dots, brown dwarfs), clumpy dusty galaxies, rejected candidate images, noise.
  Plan change if a source is a lensed image with an unmet model prediction, or `f277w_264` shows images of one SED.
- New: `scripts/a2744_followups.py` (published-image match, Bergamini+23 magnification and counter-images per
  photo-z percentile, `--knots` F200W-peak photometry in seven bands) and unit tests.
- **Result** (details and references in docs/fields/abell2744.md "Follow-up sources"):
  - 5904 and 7987 are catalogued little red dots (NED: z_spec 4.46; Kocevski et al. 2025 sample).
  - 7298 is Furtak+23's candidate image c65.13 of system 65 (z 3.519), but its NIRSpec z_spec is 3.987: excluded.
  - 264 (dusty star-forming galaxy, McKay et al. 2025): of the three F200W peaks, one is the red core
    (F150W/F200W 0.05 ± 0.03 vs 0.36–0.45), and the two blue clumps sit on one side with no counter-image.
  - 4731 is a faint compact galaxy at the noise level.
  - Bergamini+23 (`model_prediction`): four are singly imaged (μ 1.6–6.2). 4731 has μ ≈ 22 near system 700.1 in
    the model-limited infall region; its brightest predicted counter-image (F200W 27.8) has a 27.0 catalogue
    source 0.37″ away and the other three are predicted at 28.6–29.4. **Null:** no unmet prediction, no candidate.
- **Failed approach:** a 0.1″ solver grid over ±244″ did not finish in 15 min; 0.25″ over ±214″ took ~3 min.
- **Handoff:** the M1 Abell 2744 follow-up list is closed. Next field follow-ups in TASKS: El Gordo `2828`/`1601`,
  Sunrise `1756`/`1907`/`1578`. Optional: 4731's counter-images in UNCOVER v2.0.

## 2026-10-14: off-plane potentials in the Abell 370 and SMACS 0723 models: null
- **Hypothesis before computing:** a CANUCS (Abell 370) or Mahler+22 (SMACS 0723) potential with a MUSE z_spec off
  the cluster plane might, on its own plane, break the image fit or predict a robust extra image. Ordinary: real
  foreground/background structures absorbed by a single-plane fit, z_spec of a lensed knot on a member, unrefitted
  parameters. Plan change if Δχ² ≫ 10 % at every plausible σ, or a robust new extra image.
- New: `offplane_potentials.MUSE` entries for Richard+21 `a370.dat` and Caminha+22 `muse.dat`, with an optional
  per-catalogue `frame_offset_arcsec` (SMACS: Mahler+22's frame is 0.34″ from Caminha+22's; measured on 24 members,
  MAD 0.009″); unit test.
- **Result** (`model_prediction`, unrefitted; D-046 addendum d, docs/fields/abell370.md):
  - Abell 370: 36 / 284 off-plane. 23 foreground potentials at z 0.30–0.35 lower χ² 1003.0 → 985.3; 11 background
    at z 0.41–0.47 +2; 273 (a z 2.80 lensed source on a member) −3; 180 (z disputed: MUSE 0.465 vs GLASS 0.38) +12
    to +46. `images` is gated off for this source-plane-fitted model.
  - SMACS 0723: 2 / 149 off-plane: 610 is the BCG potential carrying a lensed source's z 0.94 (misassociated); 489
    is a real z 0.746 background galaxy. χ² 31.18 → 31.33; `images --plane` unchanged (75 images, no class change).
  - **Null.** The D-042 off-plane listing is now done for all five Lenstool models in the repo.
- **Failed approach:** matching SMACS without the frame offset left only 0.31–0.48″ matches near the 0.5″ radius.
- **Handoff:** next lens-model item from TASKS "Supporting"; a multi-plane refit only if a flag comes to depend on
  an off-plane potential.

## 2026-10-14: off-plane potentials in the Abell 2744 and El Gordo models: El Gordo's z 0.63 group, no tension
- **Hypothesis before computing:** a Bergamini+23 or Caminha+23 potential whose galaxy has a MUSE z_spec off the
  cluster plane might, on its own plane, predict a robust extra image or break the image fit. Ordinary: a known
  foreground structure absorbed by the single-plane fit, misassociated redshifts, unrefitted parameters. Plan
  change if Δχ² ≫ 10 at every plausible σ or a robust new extra image appears.
- New: `offplane_potentials.py` reads the authors' MUSE catalogues on CDS (Richard+21 for Abell 2744, Caminha+23
  for El Gordo; quality-flag trust, model frame) when no `--field` is given; unit test.
- **Result** (`model_prediction`, unrefitted; details docs/fields/elgordo.md, abell2744.md; D-046 addendum c):
  - Abell 2744 Bergamini+23: 72 / 180 with z_spec, **0 off-plane**.
  - El Gordo Caminha+23: 188 / 265 with z_spec, 20 off-plane, all at z 0.6305–0.6355: the paper's own
    foreground group of 20 ellipticals (positive control). Moved to z 0.6325: χ² 82.5 → 98.2 at best σ (× 1.10;
    100.3 fitted, 110.5 lensing-equivalent, 127.8 luminosity-rescaled; 312.9 removed). The paper's refitted
    model without the group has χ² 103, so the unrefitted multi-plane cost is smaller than dropping the group.
  - `images --plane`: the existing extras persist; one new system-17 extra at F200W 29.2 (undetectable). **Null.**
- **Handoff:** the D-042 off-plane listing is now done for MACS0416, Abell 2744 and El Gordo. Remaining: Abell 370
  CANUCS (CANUCS DR1 z_spec, `--field abell370`), SMACS 0723 (Caminha+22 MUSE, J/A+A/666/L9). A multi-plane refit
  stays out of scope until a flag depends on an off-plane potential.

## 2026-10-14: MACS0416 CANUCS potentials with off-plane spectroscopic redshifts: no tension
- **Hypothesis before computing:** a CANUCS potential whose galaxy has a z_spec off the cluster plane (D-042 rule)
  might carry an extra image, or the catalogued images might reject its true redshift. Ordinary: compiled z_spec
  belonging to a blended lensed image or neighbour, single-plane fits absorbing a foreground galaxy, unrefitted
  parameters. Plan change if a potential at its z_spec costs Δχ² ≫ 10 at every plausible σ, or predicts a robust
  new extra image.
- New: `scripts/offplane_potentials.py` (potential × CANUCS DR1 `Z_SPEC` within 0.5″; off-plane at |Δz|/(1+z_cl) >
  0.02; photo-z agreement |Z_ML − z_spec|/(1+z_spec) ≤ 0.15; all ASSUMPTIONs) and a unit test.
- **Result** (`validate`/`images --plane`, other parameters **not refitted**, so Δχ² is an upper bound; 303 images,
  single plane χ² 330.81; details docs/fields/macs0416.md):
  - 137 / 222 potentials have a z_spec; 12 are off-plane. 3 of them (4677 z 0.94, 8785 z 2.09, 3719 z 2.54 = MUSE
    system 45 image) disagree with their photo-z (0.31–0.39): misassociated, kept at the cluster.
  - 7 near-cluster interlopers (Δz/(1+z) 0.04–0.09, σ 40–68 km/s): χ² 330.98 (+0.17).
  - **O5** (σ 136 km/s, z_spec 0.114 Shipley+18, Z_ML 0.071; 66″ south-west, lenses D16 and D17 directly): χ² 404
    (fitted σ), 370 (150 km/s, best of 67–200), 443 (123, lensing-equivalent), 1255 (67, luminosity-rescaled).
    The Δχ² ≈ +39 at 150 km/s is spread over D16/D17 (beside O5) and D14/D20/D62 (multi-plane deflection of
    the rays). Only a refit can say how much survives; a ~12 % rise unrefitted is no tension.
  - `images --plane O5=0.1137[:150]` (F277W catalogue, no forced photometry): predicted extras change only at
    folds (|μ| 3–49) in D20, D25, D77, D115 within ≈1″ of catalogue sources; no robust new extra image. No
    existing MACS0416 flag (systems 27, 45, 51, 55) depends on O5.
  - **Null.** The 8757 explanation of system 51 (D-046 addendum) stays the only off-plane potential that matters.
- **Handoff:** next lens-model item: the same listing for Abell 2744 (UNCOVER v2.0 members are photometric; use
  UNCOVER/MUSE z_spec) and El Gordo (Caminha+23); a local refit (O5 σ, r_cut and the nearest halo) if a field's
  off-plane potential carries a flagged image.

## 2026-10-14: `--plane NAME=Z[:SIGMA]` for `lens_consistency.py validate/images`; MACS0416 system 51 check
- **Hypothesis before computing:** the D-046 explanation of MACS0416 system 51's missing fourth image (CANUCS
  potential 8757 is a z_spec 0.268 foreground galaxy; at its own redshift with σ ≤ 60 km/s the fourth image goes)
  would fail if moving 8757 broke the fit to the CANUCS catalogued images. Plan change if Δχ² ≫ 10 or a nearby
  system's images move by ≫ σ_pos.
- New: `--plane` (D-046 addendum 2026-10-14): `split_planes` on the CLI, provenance under `multi_plane`, separate
  output directory, no single-plane map comparisons for a multi-plane model. Unit tests for parsing, the move and
  the output directory.
- **Result (`model_prediction`, other parameters not refitted):** single plane χ² 330.81 over 303 images (Lenstool
  344.30). 8757 at z 0.268: χ² 330.55 (fitted σ), 332.71 (70 km/s), 336.44 (60 km/s). CANUCS D52a–c, 1.3″ from
  8757, absorb most of it (first images +0.28–0.45″ each at 60 km/s). **No tension**: the foreground explanation of the
  system-51 fourth image stands.
- **Handoff:** `images --plane` is wired but only exercised through unit tests (no field run yet). Next lens-model
  item: other model potentials with spectroscopic redshifts away from the cluster (D-042 rule) across the CANUCS
  and Abell 2744 models, then `images --plane` for any that predicts an extra image; or the next field's second
  model (Sharon v4cor, MACS1149).

## 2026-10-14: Abell 2744 both-model extra images get sibling reference fluxes: none unmet
- **Hypothesis before computing:** an extra image that Bergamini+23 and UNCOVER v2.0 both predict, left
  `no_reference` by the D-034 |μ| ≤ 50 cap and gates, would become a both-model unmet prediction if its siblings
  predict it at ≥ 10σ and nothing is there. Ordinary: fold-pair μ near a critical curve (prediction ∝ 1/μ_sib),
  resolved or contaminated references, a target on another knot's image.
- New: `scripts/a2744_xmodel_refs.py` (smacs_ref_fluxes approach on the tracked cross-model tables; S3 byte-range
  F150W/F277W/F444W stamps with the 1″ halo filter; no lens maps). Targets are the `both_loose` `no_reference`
  rows. Each sibling's peak flux / |μ_sib| × |μ_target|, and the fold pair's summed flux / summed |μ| (`pred_snr_pair`).
  `forced_class` only when every pair sibling passes the D-034 gates; otherwise `no_reference`, with the
  all-sibling verdict kept as `ungated_class` (a diagnostic). Most siblings fail the gates (resolved knots), so the
  verdicts below are ungated unless stated.
- **Result [derived; positions and μ model_prediction]:**
  - 3.2: on catalogued 3.1c (0.55″ / 0.75″); 5–8× brighter than predicted (`confused`): the galaxy's third image, seen.
  - 26.2: `recovered` in both models (9–19σ, flux 1.6–2.5× the pair prediction), 0.37″ from 26.1c.
  - 33.1 third, 62.1: predicted at only 0.1–0.7σ in every band and model; the pair's |μ| (45–216) would have to be
    overestimated 15–84× to reach 10σ (33.1 F277W/F444W: `undetectable` also with gated references). Not testable.
  - 4.2c: `absent` in both only ungated (pair 80–206σ, ratio 0.06–0.27); every sibling fails the gates (resolved knots);
    the injection test and the vetted verdict (docs/candidates/abell2744-family4-c.md) stand.
  - Cutouts viewed (all targets, three bands). **Null: no new flag.**
- `/code-review` fixes: the system's widened search radius (the stored one is NaN on `no_reference` rows), gated
  classes, `off_image`, the 1.5″ knot-match threshold applied, siblings a model leaves unmatched recorded (UNCOVER
  3.2a, 26.2a: its pair sum is then one image), sibling photometry once per system and band.
- **Handoff:** Abell 2744 counter-images are exhausted with two models. Next lens-model item: `--plane` for
  `images/validate`, or the next field's second model.

## 2026-10-14: D-041 cross-model rule on every Abell 2744 forced-photometry row: no image unmet in both models
- **Hypothesis before computing:** an extra image that both Bergamini+23 and UNCOVER v2.0 predict at ≥ 10σ, with
  nothing there, would be a lens-model residual that no single model's systematics explain. Ordinary: one-model
  predictions (D-041), fold pairs straddling a misplaced critical curve, member-buried images, fitted redshifts.
- New: `a2744_uncover.py --crossmodel --catalog --forced-image`. UNCOVER's extra images are force-photometered under
  the same D-034 rules and 1″ halo filter as Bergamini's (`predict_counter_images` from the WLS source). Each
  model's extras are matched to the other's extras within 1.5″ / 3.0″ (ASSUMPTIONs), and to its images of catalogued
  ones (an extra within 3″ of those is that image, misplaced).
- **Result (derived; positions and μ model_prediction):**
  - Extras: Bergamini 33 (12 / 15 also predicted by UNCOVER at 1.5″ / 3″); UNCOVER 45 (14 / 17).
  - Bergamini's forced classes reproduce the recorded ones (recovered 6, confused 2, absent 2).
  - **Unmet:** 4 `absent`, 0 `ambiguous`.
    - Bergamini 700.1 ×2: model-dependent (UNCOVER's nearest images are 3.6″ and 12.7″ away).
    - UNCOVER 3.1 (μ −59.5) and 34.1 (μ 28): fold pairs UNCOVER merges into one image. Bergamini predicts the
      partner 2.0″ / 2.7″ away, on the observed arc and on the recovered 34.1 image.
  - Cutouts viewed. Nothing is unmet in both models. Null; details in `docs/fields/abell2744.md`.
- `/code-review` fixes: verdicts now use what the other model predicts there (the first version labelled UNCOVER
  3.1 and 34.1 "both models (loose tier)", contradicting this entry); matching only against the other model's
  extras; NaN-free JSON; the Bergamini input must carry the same halo filter.
- **Failed approach:** the first UNCOVER map download was truncated (sha256 mismatch on alpha_y); a retry passed.
  The sha pin caught it.
- **Handoff:** five extra images are predicted by both models but untestable (`no_reference`): 3.2, 4.2c (vetted),
  26.2, 33.1 third and 62.1. Next: sibling reference fluxes for them without the |μ| cap
  (`smacs_ref_fluxes.py` approach). Other lens-model item: `--plane` for `images/validate`.

## 2026-10-14: UNCOVER v2.0 passes as the second Abell 2744 model; it also predicts 4.2c, and the vetted verdict stands
- **Route:** the UNCOVER DR4 lens maps (v2.0, Furtak & Zitrin, 2024-08-07) from Google Drive
  (`drive.usercontent.google.com`, sha256-pinned in `scripts/a2744_uncover.py`; deflection x/y 75 + 71 MB gz, 4561² at
  0.1″/pix, deleted after use) and UNCOVER's own image table (186 rows). The claim started as `bayes.dat` search
  radii and was redirected: D-045 already measured the chain spread of the flagged images (0.2–0.4″).
- **Validation gate (D-035 rule, ≤ 1.5 × the quoted 0.60″; ASSUMPTION): passes** with the magnification-weighted
  least-squares source per family (`wls_source`): own list 181 / 181 traced images re-imaged (63 families), rms **0.84″**,
  median 0.60″, p90 1.21″, max 2.31″, no family short of images; 28 images share their nearest predicted image with a
  sibling (knot or fold pairs; their offsets are optimistic). Bergamini+23's list through the same maps: rms 1.00″.
- **Failed approaches** (numbers from a first parse that split the 11 one-decimal-id systems, e.g. 83.1–83.3,
  into single images; fixed after `/code-review`): the mean back-traced source (the `imageplane_residuals` default) gives rms 4.8″: it falls
  outside the caustic for fold pairs, so 34.1 and 80.1 come back as single images 22–32″ away. The best single
  image's source gives 1.01″. Solver grid 0.1″ instead of 0.25″ changes nothing (≤ 1e-6″). The source definition
  was chosen after the mean failed (post hoc); WLS is the standard first-order image-plane χ² source.
- **4.2c, D-041 second-model test [model_prediction; observed photometry]:** UNCOVER also predicts 4.2c, μ 8.83,
  0.93″ from Bergamini's position. High-pass (1″) peak photometry within 1″ of it: F150W/F277W/F444W flux
  0.33/0.20/0.24 × the μ-scaled 4.2b expectation, at the west end of an arc-like feature 0.43–0.52″ E, and 0.32 mag
  bluer in F150W/F444W than 4.2b (cutouts viewed). Injections of 4.2b's knot at the prediction: recovered at 1× in
  all bands and at 0.5× in F277W/F444W. So both models predict an image not seen at ≥ 0.5× its flux; it is not
  model-dependent in the D-041 sense. The vetted verdict (docs/candidates/abell2744-family4-c.md: μ and position
  systematics next to member 34423, shared by any scaling-relation model) stands: UNCOVER, like Bergamini, puts the
  knots' c images 2.3″ apart against 0.45–0.49″ in a and b. Not a candidate.
- **Handoff:** UNCOVER v2.0 is usable as Abell 2744's second model (`a2744_uncover.py`). Next: the full D-041
  cross-model rule on every Bergamini forced-photometry row (37 / 63 UNCOVER families and 32 / 50 Bergamini families get extra predicted images, mostly
  third images of catalogued doubles or central demagnified images). Other lens-model item: `--plane` for
  `images/validate`.

## 2026-10-14: S0 screen picks X3 (W3 in unselected Gaia epoch photometry); X3 G0 passes on dev; Gaia DataLink fix
- **S0-2026-10-14** (`docs/state/units/S0-2026-10-14.md`; owner request to return to observational searches).
  - **A, pulsar-timing Shapiro sign test: not qualified.** A negative-mass lens would advance pulses
    (Δt = −(ε/c) ln b), so the sign of a timing bump is a real discriminator. But a NANOGrav-like array expects
    ≤ 0.08 events even for an all-dark-matter population (`scripts/s0_pta_shapiro.py`).
    - NANOGrav 15-yr searched positive amplitudes only, and not the dynamic Shapiro channel.
    - The densities a PTA could reach are already excluded by microlensing.
  - **B, FRB lensing by Ellis or negative-mass lenses: not qualified.** A single echo (Δt, R) cannot separate
    exotic from ordinary lenses in geometric optics (1302.7170).
  - **C qualifies as X3.** Gaia's Andromeda survey (GAPS) publishes G/BP/RP epoch photometry for all 1,257,319
    sources around M31 with no variability selection [observed, TAP]. DR4 is scheduled for 2 Dec 2026.
- **X3 G0** (`docs/state/units/X3.md`, `scripts/x3_gaps.py`, `results/x3/`). The dev split and the gates were
  frozen before any data.
  - **Screen v0 failed the false-alarm gate.** On 998 real dev light curves, 2.7 % flagged. 24 of 27 were
    blend-type astrometry (no 5-parameter solution, RUWE > 1.4, multi-peak) and 3 were known ECL/RS variables.
  - **Post-hoc v1 adds an astrometric-quality veto**, not the variability flag. It passed on a fresh dev sample
    of 1,000:
    - W3 efficiency **0.70** at t_E = 30–300 d (0.42 at 10 d);
    - real flags **0.1 %**;
    - PSPL→W3: 0;
    - vetoed: 7.5 %.
  - Simulated injections and dev data only. **No candidate, no held-out data touched.**
- **Fix.** The Gaia DataLink service now answers `DATA_STRUCTURE=INDIVIDUAL` with HTTP 500, which broke D-061's
  fetching. `gaia_mulens` now uses `DATAMODEL_STANDARD` (same CSV columns).
  - The manifest URIs in `w3_gaia.py` are updated.
  - The stored manifest updates on its next re-run.
- **Failed approaches:**
  - A single sync TAP dev list silently truncated at 50,000 rows in storage order.
  - The `random_index` thinning and the async jobs time out: anonymous jobs are capped at 60 s.
  - Fix: query per level-5 HEALPix `source_id` range.
- **Handoff.**
  - X3 next is G1: periodicity and achromaticity vetoes frozen on dev, then the single held-out GAPS screen. That
    screen is about 0.96 M sources, so it needs the bulk epoch-photometry files, not DataLink.
  - X3 used about 1.5 / 16 h.
## 2026-10-14: Sunrise counter-images in the four RELICS models: null (D-033 addendum 4)
- **Tool:** `scripts/sunrise_counter.py` (byte-range deflection blocks → `MapLensModel` → `find_images`, D-029 forced
  photometry on F200W with an injection check, D-041 cross-model rule); `results/sunrise/counter_images.json`.
- **Known cases:** system 2 re-predicted by all four models (rms 0.16–0.69″); system 1 only partly (1–5 of each knot's
  images). Leave-one-out control through the same chain: 16/20 `recovered`, 3 `confused`, 1 `undetectable`, 0 `absent`.
- **Result:** 0 flags. Glafic's `absent` arc-end images of 1.1 and 1.2 (|μ| 39, 79; predicted 16σ, 25σ; injections
  recovered) are predicted by 1–2 of 4 models and are its displaced 1.1c / 1.2c (5.1″, 3.3″; cutouts empty, arc ~2″
  away). The arc's counter-image across the cluster is
  predicted at S/N ≤ 1.7 in every model, at positions up to ~19″ apart. System 2 counter-images: `undetectable` or
  `confused`. Details: docs/fields/sunrise.md.
- **Handoff:** the Sunrise counter-image line is closed until a JWST-era model fits every knot. Remaining lens-model
  items: `--plane` for `images/validate`; `bayes.dat` search radii; Abell 2744 UNCOVER cross-check.

## 2026-10-14: X2 S6a writeup corrections (documents only; no computation re-run)
- **`X2_s6a_summary.md`, "Who did what".** Updated from the owner's cross-conversation record. It replaces the
  earlier "no recorded GPT output":
  - GPT read the paper, the code and the results. It helped with scope and verification design, and flagged three
    issues: min vs last step, the merger-grid differences and NaN.
  - The owner relayed these to Claude.
  - GPT did not re-derive the equations or re-run the numerics.
  - Still no domain-expert review.
- **`X2_s6a_repro.md` §7.**
  - The full rotation range, other sectors and time-domain evolution are now listed as possible strengthening,
    not as presumed submission requirements.
  - The claim stays limited to the two verified points. Sufficiency is for experts and editors to judge.
  - The Comment assessment is marked provisional, not an exclusion.
- **Literature-search scope** is now stated wherever it is summarised: 47 citers screened by title, relevant
  abstracts and some full texts checked. The search itself was not re-done.
## 2026-10-14: PSF-commuting parity fit makes the acceptance test worse (76 %); the PSF does not explain the flags; tool paused
- **Hypothesis before computing:** the 3 persistent Abell 2744 parity flags (1.1b, 7b, 600.3b; previous entry) come
  from the unlensed PSF, which the free affine map distorts and mirrors in A but not in B. If a statistic exact for
  an affine map under a PSF (P ⊗ T(A) = K ⊗ B, K(d) ∝ P(M d)) clears them and the acceptance reaches ≥ 95 %
  (ASSUMPTION), the tool can run on other fields; if the flags persist, the PSF route is falsified.
- New in `scripts/parity_test.py` (default; `--no-psf` keeps the plain statistic): `empirical_psf` (median stack of
  Gaia DR3 stars on the UNCOVER F200W mosaic, recentred, local pedestal subtracted, peak-rule rejection; cached in
  `results/parity/psf-jw02561-o001_t003_nircam_clear-f200w_i2d.json`; mirror asymmetry 8–15 % of the flux),
  `PsfPair` (the PSF-commuting statistic), `refine_psf` / `compare_psf` (refined from the plain optima, the
  predicted map and reflection starts of both), PSF-consistent synthetic pairs (A_syn = P ⊗ S, B_syn = P ⊗ T(S),
  S = smoothed A), and the stability rule (`unstable` unless 16- and 24-angle start densities agree). Test
  `test_psf_commuting_statistic`: on a compact source through a toy asymmetric PSF, χ² at the true map is 1147 for
  1185 pixels (plain statistic 34 483); handedness recovered both ways; a wrong prediction is flagged.
- **Result (derived; D-032 addendum; `results/parity/abell2744-bergamini23-psf.json`, ~20 min on 4 cores): worse.**
  First run 51 / 68 = 75 %. `/code-review` then found three defects (B's mask not applied to T(A) before its
  convolution, no PSF pedestal subtraction, a cache keyed on the image only); fixed and rerun: 27 `consistent`,
  13 `untestable`, 10 `parity_residual`, 6 `unstable`; acceptance (consistent + untestable + unresolved)
  52 / 68 = 76 % (plain: 87 %). 1.1b and 7b stay flagged and 600.3b becomes `unstable`, so the unlensed PSF does
  not explain them; 8 new flags (2.1a, 2.1b, 2.4d, 2.5c, 2.6d, 2.7c, 24b, 600.5b). Not candidates; no parity
  conclusion. Known remaining limit: K = P(M d) is truncated at 0.6″, exact only when M's smaller singular value
  ≥ 0.83.
- **Looked at the contact sheet** (1.1b, 2.1a, 2.4d, 7b, 2.3b): 2.4d is a long knotty arc with striped residuals in
  both families (the map is not affine over the 0.6″ disc); 7b keeps a core residual in the predicted family; 2.1a
  is a faint diffuse image beside an arc. |Δχ²| grew ×4–8 on every pair: convolving both sides correlates the noise,
  which the per-pixel χ² ignores, so the fixed DCHI2_MIN = 100 is effectively lower (ASSUMPTION: the main cause of
  the new flags). A goodness-of-fit gate (min χ²/n ≤ 2, applied afterwards to the stored tables) removes only 2 of
  first run's 12 flags and none of the plain run's 3. VizieR returned empty results for the Gaia cone during the
  rerun, so the cached Gaia rows of the first run were reused (`empirical_psf(stars=...)`).
- **Failed approach:** PSF commutation as the fix for the parity flags. **Next:** the tool is paused (three acceptance
  failures). Reopen only with a non-affine forward model (source-plane reconstruction through the published
  deflection field) or a Δχ² threshold calibrated for correlated noise on PSF-consistent synthetics.

## 2026-10-14: halo-subtracted flux ratios on MACS0416 (CANUCS model): five systems testable, no anomaly
- **Hypothesis before running:** a fourth cluster extends the D-032 flux-ratio line; no compact image deviates by
  > 0.75 mag (both routes) from the `macs0416-canucs` μ ratio, and `fluxratios --halo-phot` finds no 3σ tension.
  A flag would go to `/vet-candidate`.
- Unchanged code: `halo_phot.py` on all 303 images (111 systems; F150W + F444W, ±1 mag injections; 18 min on 4 cores,
  dominated by serial S3 stamp reads), then `fluxratios --halo-phot` F444W, both routes
  (`results/halo_phot/macs0416/`, `docs/fields/macs0416.md`).
- **Result (derived): null.** 5 systems testable (D4, D13, D21, D98, C3; F444W 11 images, F150W 7), max |r| 0.35 mag;
  `fluxratios` 14 / 11 images, σ_sys 0.25 / 0.26 mag, 0 `flux_tension`, 0 colour outliers. Injections (simulated):
  1 / 2 mag deficits flagged in 50 / 29 % (medfilt) and 36 / 36 % (annulus). Cutouts looked at: compact sources.
- **Limitation:** no bayes.dat for this model (no μ error per image); the q ≥ 0.80 gate uses the three-field star
  calibration. ~96 % of images fail the compactness gates, as in the other fields.
- **Next:** the halo-flux line is now null in four clusters; MACS1149 needs a Lenstool par (CATS is map-only for
  `fluxratios`). For arcs, a surface-brightness-conserving test (D-032) is the remaining flux tool.
## 2026-10-14: parity acceptance test fails (87 % vs 95 %); fits not converged, PSF suspected
- **Hypothesis before computing:** the 6 Abell 2744 parity flags (previous entry) come from neighbour light and
  point-like sources; with a segmentation mask and a resolution gate, ≥ 95 % of catalogued non-near-critical pairs
  (ASSUMPTION) are `consistent` or `untestable`. If not, the test is not used on any field.
- New in `scripts/parity_test.py`: `segment` (own vs neighbour light per stamp), `size_ratio` and the `unresolved`
  class, `mirror_starts` (72 reflection axes on the other family's best map), denser default starts, a process pool;
  `--no-mask` keeps the first version. Test: `test_segment_separates_neighbour_and_size_ratio`.
- **Result (derived; D-032 addendum): acceptance failed, 59 / 68 = 87 %.** 56 fitted pairs: 27 `consistent`,
  20 `untestable`, 6 `undetermined`, 3 `parity_residual` (1.1b, 7b, 600.3b); QSO1 and 42.1 now `unresolved`
  (`results/parity/abell2744-bergamini23.json`). Not candidates; no parity conclusion.
- **Failed approaches:** (1) the original cross-seed scanned the mirror axis in 15° steps; on sharp cores this left
  the predicted family Δχ² ~ 10³ worse (dipole residual at the core on the contact sheet). Reflection-axis starts
  rescued 2.2c, but (2) classes still depend on the start set: 600.4a went from Δχ² +32 030 to −52 264 to +3 329,
  2.5d from −111 to +183. (3) A toy asymmetric PSF (lobe + spike) reduced synthetic mirror-pair power but flipped no
  pair, so it does not yet explain the 3 persistent flags. All 3 are predicted-mirror pairs with bright compact cores.
- **Next:** a PSF-commuting fit, P ⊗ T(A) vs T(P) ⊗ B with an empirical F200W PSF stacked from mosaic stars (exact
  for affine T), plus a stability rule (class only when Δχ² keeps its sign across two start densities); then rerun
  this acceptance test. No field run before it passes.
## 2026-10-14: observational parity test built; validated on synthetics, fails on known Abell 2744 systems
- **Hypothesis before computing:** the relative parity of two images of one source, sign(μ_A μ_B), is visible as
  mirror vs non-mirror morphology; a resolved pair preferring the non-predicted handedness, far from critical curves,
  would be a lens-model residual. Ordinary: symmetric or unresolved sources, ICL/member light, PSF, non-affine
  structure, near-critical images (model-dependent parity). Plan: `docs/plans/parity-test.md`.
- New: `scripts/parity_test.py` (free affine fit of B by A in both handedness families; cross-seeded refits;
  per-pair power from synthetics drawn from a smoothed A; 8 threaded S3 byte-range readers) and
  `tests/test_parity_test.py` (handedness recovered both ways and flagged when the model is mirrored; symmetric
  source `untestable`).
- **Result (derived; D-032 addendum): the tool is not validated on real data.** Abell 2744 Bergamini+23, UNCOVER
  F200W, 1″ halo filter: 116 pairs, 31 `consistent`, 6 `parity_residual`, 4 `undetermined`, 27 `untestable`,
  37 `near_critical`, 6 `off_image`, 5 `faint` (`results/parity/abell2744-bergamini23.json`). The 6 flags include the point-source
  quasar A2744-QSO1 (53a/b) and system 1: tool failures, not residuals; not candidates; no parity conclusion.
- **Failed approaches:** (1) synthetics drawn from A itself share A's noise pattern and gave symmetric sources Δχ²
  ≈ ±100–160 of fake power; fixed by drawing from A smoothed by 2 px. (2) a single optimisation per family stalled
  in local minima (a symmetric source gave −30 in the wrong direction); fixed by grid starts plus cross-seeding.
  (3) a template chosen by peak S/N without halo subtraction picks images on the BCG halo (1.1b, 2.7d: looked at the
  cutouts); 6 of 15 decided pairs flagged. (4) 2.5″ stamps read row by row from S3 took minutes per image; 1.6″
  stamps and 8 readers take ~1 s.
- **Next:** before any field run counts, an acceptance test on catalogued non-near-critical pairs (≥ 95 %
  `consistent` or `untestable`), with member + ICL models subtracted (`member_sub.py`-type) and a fit mask from the
  source's own segmentation; then El Gordo, MACS0416/MACS1149, SMACS.
## 2026-10-14: member-galaxy subtraction for member-buried counter-image rows: one testable row, already explained
- **Hypothesis before computing:** the `confused` rows the median filter could not clear (D-029 addendum, the halo-filter
  entry) become decidable once the bright neighbours are modelled; an `absent` among them without a recorded
  explanation would be a lens-model residual. Ordinary: model position/topology errors near critical curves,
  galaxy-scale caustics, member-model residuals, uncatalogued light.
- New: `scripts/member_sub.py` (Sérsic per catalogued neighbour brighter than the prediction, search disc masked
  except the cores, empirical ring noise, injections, control discs rotated about the neighbour) and
  `tests/test_member_sub.py` (synthetic member + image: recovered at ratio 1.03, empty position `absent`, blend,
  faint neighbour). Runs: `images --forced-image <F277W i2d, S3>` (no halo filter) then
  `member_sub.py --table <images_predicted.ecsv> --catalog <F200W cat> --forced-image <F277W i2d> --classes
  confused,recovered,absent` for `macs1149-cats`, `macs0416-cats`, `abell2744-bergamini23` (catalogue sha256 `9953c869…`;
  no photo-z: the forced classes do not use it). Outputs: `results/member_sub/*.json`.
- **Result (derived; D-029 addendum): no new unexplained `absent`.** Of 31 rows, 1 testable (MACS1149 system 2,
  `absent`, 0.30 × prediction vs 0.32 × in the controls; explained by the CATS/Sharon topology disagreement; looked
  at the residual with the prediction injected), the rest `untestable` (control residuals 0.6–660 × the prediction),
  `blend` or `no_member`. MACS1149 21's forced `recovered` sits on a member: doubtful.
- **Failed approach:** initialising Sérsic amplitude with the peak pixel (Sersic2D's amplitude is I(r_eff), ~100×
  lower for n = 2.5) left 56× residuals on a synthetic member; fixed. First version masked only up to 0.3″ short of
  the core (0.95″ of a 2.25″ search disc), and selected neighbours within 2.5″ only; both widened to the search disc.
- **Next:** none on these rows unless a PSF-convolved member model is built for another reason (D-029 addendum);
  other Supporting 1 items (`--plane` option, orphan-pair injections).
## 2026-10-14: X2 S6a writeup: research summary, reproduction package and prior-work check (no new physics)
- **Summary** (`docs/state/units/X2_s6a_summary.md`). It is limited to the two fully verified points, J/A = 0.121
  and 0.142. Published results, published conjectures and what this work adds are kept apart. Contributors are
  listed: the owner, Claude, no recorded GPT output, and no human expert yet.
- **Reproduction package** (`X2_s6a_repro.md`):
  - equations (V₁₀ corrected), parameters, methods, criteria, a results table and run commands;
  - venue assessment: a Comment on 1309.2448 is poorly suited; a short article is possible only after a wider
    J/A range, odd parity and an expert check.
- **Prior-work check** (`X2_s6a_prior_work.md`).
  - Scope: all 47 INSPIRE citers of 1309.2448, author queries, arXiv and Crossref.
  - No prior computation of the post-merger pair was found. 2510.11406 conjectures it.
  - No erratum or typo note for V₁₀ was found.
  - "First" is not used.
- **Audits:**
  - `Shooter.root` returned the minimum step, not the last, while the docs said "last". The full step record was
    added and the shooting re-run (the static gate was recomputed whole). Every point reached 1e-9 (last = minimum), the roots are bit-identical,
    and no status changed.
  - At the merger, the eigenvalues at the bracket ends differ between grids (pair ±0.014i against ±0.079i), as
    expected near an exceptional point. J_EP moves by 4e-6 in J/A. The representative points are unaffected.
  - NaN is replaced by null, and the results are now strict JSON.
- **Handoff.** X2 is at about 13.9 / 16 h. The next step, a human expert review, is the owner's decision. Agents
  will not submit anything or contact authors.

## 2026-10-14: X2 S6a: in 5D equal-J rotating Ellis wormholes the two unstable modes continue past their merger as an oscillating growing pair
- **Equations re-derived** (sympy, 5D Ricci of 1309.2448 Eq. (pertlineel)).
  - The background equations and R_ll are reproduced exactly.
  - The gauge H₁ = H₃ − H₂ with H_φ = 0 is legal for complex Ω, via the δΦ = 0 gauge and the scalar equation.
  - The G system matches the printed one except V₁₀. The re-derived equations reproduce the published static
    eigenvalue (−21.80); the printed V₁₀ gives −146. This supports a typesetting error, but does not establish
    what the authors' code did.
  - The printed H_ω′ (Eq. hop) is 2× the derived one; it does not enter the G system.
- **Gates passed:**
  - background: constraint ≤ 3.4e-7, Smarr ≈ 2e-7, digitised Fig. 1 within its read-off precision (≤ 0.3 % relative);
  - units: 4πG = 1, from the extremal Myers–Perry point;
  - static mode −21.801758 (dense and shooting);
  - merger at J/A ∈ [0.100622, 0.100633], against 0.10063.
- **Result (even parity, paper's ansatz).** Past the merger the modes form a complex-conjugate pair, i.e. an
  oscillating growing mode.
  - **Verified** at J/A = 0.121 and 0.142, with all six frozen checks passed:
    - grid, box and outer radius;
    - complex eigenpair residual;
    - two un-reduced constraint equations converging at second order;
    - independent Evans shooting within ≤ 1e-7.
  - Growth is about 2.3 A^{−1/3}, half the static rate.
  - J/A = 0.214 is a **candidate**: the shooting-radius check gives 4e-6 against a limit of 1e-6.
  - Exploratory tracking reached 0.40.
  - Odd parity was screened only (no recurring unstable eigenvalue in 0.3 < |z| < 60); the near-extremal limit is
    unresolved.
- **Post-hoc method changes** (spectrum filter, shooting radius and stop rule, thresholds) are recorded with their
  before/after values. They were frozen before the verification run, not set in advance.
- **Failed approaches:**
  - a two-box filter at the same grid step lets grid modes through;
  - shooting from l ≈ 200 loses the determinant to cancellation;
  - `pkill -f` with the script name in the same shell killed the shell (again).
- **Handoff:** X2 at about 13.1 / 16 h, paused for the owner. No 4D, NED + SV or communication conclusion.
## 2026-10-14: halo-filtered counter-images in Abell 2744, MACS1149 and MACS0416: no new absent image
- **Hypothesis before computing:** as for SMACS system 17, rows that earlier `images` runs left `confused` or
  marginal on BCG/ICL light may turn `absent` (a ≥ 10σ prediction unmet) once the halo is median-filtered.
  Ordinary: member light the filter cannot remove, model copies and galaxy-scale caustics, fitted redshifts.
- Runs (`lens_consistency.py images --forced-image <F277W i2d, S3> [--halo-filter-arcsec 1.0 | 0.6]`; Abell 2744 with
  the D-029 catalogue and DJA photo-z, sha256 as `docs/fields/abell2744.md`; MACS1149/0416 with the F200W pipeline
  catalogues, sha256 `7f964bf8…` / `00cb6e18…`, no photo-z, which the forced classes do not use). Unfiltered runs
  reproduce the recorded classes (MACS1149 2/2/1, MACS0416 3/13/2/3; Abell 2744 3/3/2 with undetectable 4 vs 3, as the
  solver now predicts 178 images, not 176).
- **Result (derived): no new `absent` in any field.** Forced classes, unfiltered → 1.0″: Abell 2744 recovered 3 → 6,
  confused 3 → 2, absent 2 → 2; MACS1149 recovered 2 → 2, confused 2 → 4, absent 1 → 1; MACS0416 recovered 3 → 3,
  confused 13 → 11, absent 2 → 2. The `absent` rows are the explained ones (700.1 ×2, MACS1149 16, MACS0416 27 ×2).
- Looked at raw and 1.0″-filtered 5″ cutouts of the 16 `confused` rows with predicted S/N ≥ 10: member cores
  survive a 0.6–1.0″ filter (neighbour 4–900× the prediction). At the exact positions the filtered S/N is below the
  prediction in 14 of 16, but all have recorded ordinary explanations, and next to bright cores the annulus is biased
  (down to −8σ). **Not candidates.** Median filtering fixes smooth halos (SMACS 17), not member-buried images.
- Fix: the filter floor `3 × FORCED_R_AP` was 0.6000000000000001 and refused `--halo-filter-arcsec 0.6` (test added).
- **Next:** member subtraction (isophote/Sérsic models) in `forced_check` for member-buried rows; otherwise the
  other Supporting 1 items (`--plane` option, orphan-pair injections).

## 2026-10-14: SMACS systems 11, 16, 26: no reference flux makes their counter-images testable
- **Hypothesis before computing:** D-029 left systems 11, 16 and 26 `no_reference` (catalogued images at |μ| > 50
  or S/N < 5). Dropping the |μ| cap, a catalogued sibling may still predict an uncatalogued image at ≥ 10σ; if
  nothing is there in F200W and F444W under more than one model, that is a lens-model residual. Ordinary: μ of
  merging pairs near a critical curve is ill-determined, fitted redshifts, resolved or neighbour-contaminated
  references, model dependence.
- New: `scripts/smacs_ref_fluxes.py` (ICLv2 plus the RELICS Glafic v2 and Zitrin LTM v2 maps of
  `smacs_s17_vet.py`, each system at its ICLv2-fitted z: 1.65, 1.10, 4.58), F200W and F444W by S3 byte range; per
  predicted image the S/N each sibling predicts (all, and those passing the D-034 gates except the |μ| cap), the
  μ-error factor needed for 10σ, and the best forced S/N within 1″.
- **Result (`results/smacs_ref/ref_fluxes.json`; model_prediction vs derived): no testable image, so no `absent`.**
  - System 11: both catalogued images are undetected (S/N ≤ 3.7 at the peak, ≤ 2.7 at the position, both bands); the
    counter-image (μ 2.3–4.2 in the three models) is predicted at ≤ 1σ. Untestable at this depth.
  - System 16: the catalogued pair is one merging knot (μ 110–340 in ICLv2/Glafic); its counter-image is predicted at
    ≤ 0.6σ (a 17–28× μ error would be needed for 10σ). 16.1 and 16.2 share one 0.4″ peak, so they are one reference.
    Zitrin (pair μ 10–13) predicts a third image 2.5″ from the pair at 18σ from a gated reference: the 1″ search
    finds 41σ (×2.4, `recovered`) on the edge of a bright arc (looked at: the arc runs ~2″ N of the prediction, so
    the identification is uncertain); its μ 1.9 image is predicted at 3.9σ.
  - System 26: only ICLv2 predicts counter-images; Glafic and Zitrin give one image at the pair (μ 35 / 5.3; a
    merging pair the 0.25″ solver grid does not split) and no distant counter-image. ICLv2's
    μ 26.6 image is predicted at 10.3σ (F200W) / 10.7σ (F444W) only from 26.2, whose reference fails the gates: its
    0.4″ peak is 2.3–3.2× the flux at the catalogued spot (looked at: the wing of a bright compact neighbour
    1–2″ away; concentration 0.29–0.33). From 26.1 (diffuse, concentration 0.48–0.53) it is 3.3σ / 2.2σ, and the
    sky there shows 4.5σ / 2.9σ within 1″ (as D-029's 4.6σ scratch run). Its μ ≤ 5.8 images: ≤ 2.3σ.
- **Ordinary: the counter-images of these merging, high-μ pairs are 6–60× fainter than the pairs, below the
  forced-photometry depth in every model, except Zitrin's system 16 third image, which is recovered. The one other
  ≥ 10σ prediction (26, ICLv2) rests on a neighbour-contaminated reference and a single model. No `absent`; not
  candidates.** This closes the D-029 `no_reference` systems at the current depth.
- Review fixes (`/code-review`): wording for the recovered Zitrin image and the one-image RELICS pairs; NaN/zero μ
  and non-finite photometry guards; 4-decimal fluxes in the JSON. Not changed (recorded): the D-034 gates are
  restated, not shared with `forced_check`, and the RELICS loader duplicates `smacs_s17_vet.relics` (refactor
  when either changes); no halo filter (the systems are 20–36″ from the BCG).
- **Next:** `--halo-filter-arcsec` on the other fields' `images` runs near BCGs (Supporting 1).

## 2026-10-14: SMACS system 17 near the BCG: halo-subtracted counter-image absent, explained by its free redshift
- **Hypothesis before computing:** D-029's marginal recovery of system 17's μ −10.2 image (7.8σ, ratio 0.36) sits on
  the BCG's light gradient; with the halo removed, either a compact image remains (recovered) or the 7.8σ was halo
  light and a ≥ 10σ prediction is unmet. Ordinary: halo residuals, reference choice, model position error, the
  fitted (not spectroscopic) redshift, model dependence of central images.
- New: `lens_consistency.py images --halo-filter-arcsec` (median-filtered stamps, `halo_recovery` from ring
  injections). Reproduced the D-029 run first (24.4σ / 7.8σ / 0.36; the μ 5.4 image `confused` by a cluster galaxy).
  System 9's third image is now `no_reference` without the filter (the D-034 reference gates).
- **Result (derived, `results/smacs_s17/`):** filter 1.0″: μ −10.2 image `absent`, predicted 35.7σ, best 3.6σ, flux
  ratio 0.09, recovery 0.98; filter 0.6″: 25.5σ / 2.7σ. μ 5.4 image still `confused` (×81). Looked at the raw and
  filtered F200W cutouts: nothing within 1″; within 2″ only an unrelated compact source (7.9σ, ratio 0.22).
- **Vetted (`scripts/smacs_s17_vet.py`, model_prediction):** z is fitted (Mahler 2.12 +0.11/−0.09; Caminha
  1.84–1.92; Pascale 1.72; photo-z 2.34; literature check by a subagent, arXiv:2207.07101 Table 2, 2207.07567,
  2207.07102). In ICLv2 the central pair merges at z ≈ 1.92 (μ −36 / +31, best 3.3σ) and is gone at z ≤ 1.84. RELICS
  Zitrin LTM v2 (best z ≈ 1.5, spread 0.37″) predicts no central pair; RELICS Glafic v2 (z ≈ 1.5, 0.33″) predicts one
  0.9–1.5″ from ICLv2's, with nothing of the predicted flux within 1″ (3.0σ, 0.73 against 2.1–4.7 from the siblings;
  RELICS maps shifted by the ICLv2 frame offset, ASSUMPTION). **Ordinary: a near-caustic prediction that depends on the
  free redshift and the inner (BCG) mass profile. Not a candidate;** the absence is a constraint on those models.
- **Failed approach:** the annulus-only background next to the BCG (D-029): the gradient made the 7.8σ "recovery".
- Review fixes (`/code-review`): RELICS maps in the JWST frame, forced photometry at their images, maps deleted after
  use, a filter floor (≥ 0.6″), one shared median-filter helper. Not changed: `halo_recovery` does not rescale the
  prediction, because the reference passes the same filter (losses cancel for compact images; recovery 0.94–0.98 here).
- **Next:** reference fluxes for systems 11, 16 and 26 (Supporting 1); `--halo-filter-arcsec` on the other fields'
  `images` runs near BCGs.

## 2026-10-14: counter-images of uncatalogued bright arcs: five fields null (Abell 370 3881 is known system 7)
- **Hypothesis before computing:** every strongly sheared background arc that is not a catalogued image has the
  counter-images its model predicts at its photo-z; a robust predicted image (z16 and z84 agree within 3″) of ≥ 10σ
  with nothing detected would be a lens-model residual. Ordinary: photo-z error, knots of catalogued systems,
  cluster galaxies/ICL, model extrapolation, extended-source flux scaling.
- `scripts/arcs_counter.py` (new; `lensmodel.find_images`, the D-029 forced photometry, F200W i2d by S3 byte range)
  on the stored `results/arcs/*/arcs.ecsv`. Positive control (`--control`, catalogued arcs, Abell 2744): 18 of 20
  known counter-images not absent; of the two `absent`, 7247 ( a catalogued image detected at 86σ, flux ratio 0.27) shows
  the isophotal-to-aperture prediction over-predicts bright extended arcs by up to ×4, so a flag needs best S/N < 5.
- **Result (model_prediction vs derived, `results/arcs_counter/`), after the exclusion fix below:** Abell 2744 2 arcs /
  5 counter-images, Abell 370 12 / 42, El Gordo 4 / 9, MACS0416 4 / 7, MACS1149 3 / 7: **0 flags in five fields**
  (3 more images untested: their search circle would reach the arc).
- **Abell 370 3881 vetted (was the one robust flag: predicted 237σ, best 4.5σ, 64″ out):** it is image 7.5 of system 7
  (z_spec 2.751) in the HFF CATS v4 `arcs.txt` (0.0″); CANUCS lists no image within 7″ of it, and the A370 image list
  is gated off (D-035, D-044), so the 2″ exclusion never ran. Catalogued 7.4 lies 2.5″ from the flagged position and is
  detected (F200W source 2578, 24.05 mag against 23.6 predicted; looked at: compact and elongated like 3881). At
  z_spec, CANUCS puts the image 1.6″ from 7.4; the CATS model predicts it 3.3–6.4″ from the flag at z 2.6–3.4. The
  z_phot 2.99 (z16–z84 spread only 0.47″) shifted the prediction 2.6″, beyond the 1″ search. **Ordinary: photo-z
  error on a catalogued system.** Not a candidate.
- **Fix:** `known_images` excludes arcs near the model's own image list even where it is gated off for position tests,
  plus the HFF CATS list of the same cluster. 16 Abell 370 arcs, and one more in Abell 2744, are now excluded.
- **Controls (`--control`, catalogued arcs, all fields):** MACS0416 22 counter-images, 0 absent; MACS1149 8, one
  `absent` that is not a flag; El Gordo none testable. Abell 2744 20, flags 2528 (the system 26 knot with z_phot
  5.16 against z_spec 3.054, now in the control because it is within 2″ of 26.1a). Abell 370 8 arcs / 26, flags only 3881 (system 7 at z_phot 2.99,
  as above).
- **Review fixes (`/code-review`):** the arc no longer counts as its own counter-image's detection (catalogue match
  excludes the arc's label; images whose search circle comes within 2″ of the arc are untested); NaN magnitudes never
  win the catalogue match; `robust` needs both z16 and z84; masked redshifts are NaN; stale `counter.ecsv` is removed;
  the chain check runs after the outputs are written. Not fixed (recorded): a fixed 2″ self radius for long arcs,
  and `nearest_image` accepting any image within 4.2″ (only `--chain-samples`, not used here).
- **Lesson:** "robust" (z16 and z84 agree) does not guard against a wrong photo-z: both false flags (2528, 3881) were
  photo-z errors on known systems. A future flag needs a spectroscopic redshift or an independent-model image search
  wider than the 1″ floor before vetting.
- **Failed approach:** a 1″ exclusion radius. Abell 2744 arc 2528 lies 1.1″ from 26.1a; its photo-z 5.16 (z_spec
  3.054) moved the predictions up to 8″ and made a 23σ `absent` (looked at: F115W–F444W cutouts). At z = 3.054 its
  images land on 26.2b (0.6″) and 26.1c (1.5″). Exclusion is now 2″.
- `recovered` among unknowns is weak: a random 5σ source within 1–3″ at 1/3–3× the flux is common in deep fields
  (2528's z = 5.16 "counter-image" was recovered 4.7″ from 26.1c); only `absent` with no detection is informative.
- **Next:** the counter-image test of uncatalogued arcs is closed at this depth (null). Remaining Supporting 1 items:
  the BCG/ICL-subtracted test for system 17, reference fluxes for systems 11, 16 and 26.

## 2026-10-14: X2 S6 screen: rotation is not an established wormhole stabiliser; one checkable gap (5D complex modes)
- **Sources read.** 2403.08387 (4D, second-order slow rotation) and 1309.2448 (5D, full rotation); 2510.11406 v1
  and v2 compared. Follow-ups via INSPIRE; no errata and no public rotating-wormhole code.
- **Shown:**
  - In 4D (perturbative), the main radial mode weakens (C = 0: ω_I = −1.182 + 16.358 J²). That static value
    matches our S5 tool's 1.18199.
  - A second mode grows from the zero mode, and the two branches cross near J/A ≈ 0.013. That they merge and
    vanish is a conjecture.
  - In 5D, the purely imaginary modes merge and vanish at (|J|/A)_cr = 0.10063, but only real Ω² was searched.
- **Unsettled.** In the charged analogue the merged modes persist as a complex pair, with a growth rate that
  goes to zero near extremality (2510.11406 v2). For rotation the authors draw only an analogy ("long, close to
  extremal Kerr") and defer to a non-perturbative study.
- "Only real Ω² was searched" in 1309.2448 is our inference from its method.
- **Open:**
  - 4D full-rotation radial M_z = 0 (no calculation exists);
  - complex modes after the 5D merger;
  - non-radial sectors (charged: polar ℓ = 2 instability, 2607.18399).
- **Proposed S6a** (≤ 6 h, awaiting the owner): the complex spectrum of 5D equal-J rotating EB wormholes beyond
  J_cr, as ODEs, with a reproduction gate at J_cr.
- **S5 scope fix.** The wording is now "this background is radially unstable". ε ≲ 10⁻⁵⁸ is a linear-extrapolation
  requirement, not a computed merger initial condition or collapse time.

## 2026-10-13: X2 S5: the S4 echo background is radially unstable (e-folding 9.9 M, 139 e-folds per echo spacing)
- **Equations.** CPS framework (arXiv:1911.13103) in δΦ = 0 gauge: δr free, so the throat may move. Magnetic NED
  enters as a fixed function of the perturbed areal radius, since F is fixed by the Bianchi identity.
  - Background: Φ′² = 2a²/r⁴ (phantom), U = 6Ma²/(5r⁵), V = 4Ma²/(5r⁵).
  - The Ellis–Bronnikov master equation is reproduced symbolically.
- **Result.** At a − 2M = 3.575e-4: λ = 0.1011244 M⁻¹, i.e. 3.35 ms for GW150914.
  - Converged in grid, domain and initial data; the independent C-equation residual is 8.5e-6 and falling.
  - Ellis–Bronnikov reference 1.18199 (CPS 1.18207).
  - The mode sits at R ≈ 2.05M on each side of the throat.
  - The first echo is 133 e-folds, the spacing 139. Growth time is not collapse time.
- **Decision.** Echo searches on this background stop (owner rule).
- **Scope fix to S4b.** Both EM figures rest on the heuristic relative energy weight. Withdrawn: "instability
  makes echo bounds moot"; the growth time must be compared with the window instead.
- **Failed approaches.**
  - Evolving C on its own admits a constraint-violating pointwise mode (0.168 M⁻¹).
  - Single-sided shooting misplaced the mode at the throat.

## 2026-10-13: forward-modelled arc curvature: all 5 testable arcs bend as their models predict; El Gordo 2527 explained
- **Hypothesis before computing:** El Gordo 2527's bend away from −α (D-024 addendum 5) is a saddle effect that
  the Caminha+23 model reproduces; if the model bends it the other way, it is a lens-model residual.
- `scripts/arcs_curvature_model.py` (new): a round 0.1″ source at the arc's source position, ray-traced through the
  model, fitted with the same parabola; validated on an SIS in the tests.
- **Result (model_prediction vs derived, `results/arcs_curvature/*/forward.*`):** 5 of 18 curved arcs have a
  predicted bend ≥ 3σ; all 5 agree with their model, including El Gordo 2527 (model 12.8σ away from −α, as
  observed). The other 13 are not testable (predicted bend < 3σ). Null (D-024 addendum 6).
- **Next:** curvature is closed at this depth; critical curves at z_s = 1, 2, 4, or blend rejection in F200W
  segmentation for the `arcs` test (TASKS Supporting 1).

## 2026-10-13: arc curvature against the models' mass direction: not a lens-model test at these arc lengths
- **Hypothesis before computing:** a tangential arc bends around the mass (centre of curvature toward −α); a long,
  single arc that bends away would be a lens-model residual. Ordinary: saddle points between clumps, radial arcs,
  blends, intrinsically curved galaxies.
- `scripts/arcs_curvature.py` (new; F200W i2d S3 cutouts via `cutouts.make_cutouts`, halo median filter, parabola in
  the major-axis frame) on the stored `results/arcs/*/arcs.ecsv`; 43 long background arcs in strong shear.
- **Result (derived, `results/arcs_curvature/`):** 11 agree, 7 reverse, 24 flat; the aligned control is 9 / 7, so
  the sign carries no lensing signal. The circular-lens curvature c_pred = 1/(2r) is 0.011–0.024 arcsec⁻¹ against a
  measured median |c| of 0.056: intrinsic bends, knots and blends dominate. 5 testable arcs (c_pred ≥ 3σ_c): 3
  agree, 1 window-clipped (Abell 370 4613), El Gordo 2527 reverse in the bimodal merger, where c_pred does not apply.
  Looked at all five contact sheets. Not a flag (D-024 addendum 5).
- **Failed approach:** the circular-centre sign test; the formal σ_c from pixel residuals also understates the
  scatter (spurious 3–7σ on faint 1″ blobs: Abell 370 4401, 4429).
- **Next:** forward-model each arc (ray-trace an elliptical source at its source-plane position through the model,
  fit the same parabola) for the arcs ≳ 3″, starting with El Gordo 2527; widen the window for Abell 370 4613.

## 2026-10-13: arc orientations in F444W vs F200W: the same axes, no new anti arc (five fields)
- **Hypothesis before computing:** lensing is achromatic, so a lens-set axis agrees between F200W and F444W; a blend
  of differently coloured objects need not. An F200W anti source that is aligned in F444W would be a blend; one that
  stays anti is a single source whose shape lensing does not explain.
- `scripts/arcs_bands.py` (new; reuses `lens_consistency` `load_shapes`, `orientation_table`, `anti_expected`) on the
  stored `results/arcs/*/arcs.ecsv` with the F444W pipeline catalogues of the same observations (manifest sha256s).
- **Result (derived, `results/arcs_bands/`):** background, strong shear, 90 of 154 usable in F444W (60 unmatched, 4 below e 0.3). Control: 81 of 83
  F200W-aligned stay aligned, 0 turn anti; median axis difference 1.7–3.9°. F444W anti 3 vs 1.30 expected
  (P = 0.14; null at the F200W e ≥ 0.5 cut): Abell 370 340, 4540 and El Gordo 742, each one source with the same axis in both bands (looked at:
  6″ F200W/F444W cutouts). 1565, 1788, 4236 and 1703 have no F444W source within 1″. Null (D-024 addendum 4).
- **Limitation:** an intrinsic shape is achromatic too, so the band test is a blend veto, not a lensing test.
- **Next:** critical-curve curvature against arcs; F444W shapes for the 60 unmatched sources need forced moments on
  the i2d (S3 byte range), only if a new anti flag appears.

## 2026-10-13: X2 S4b: physical strain and energy of the S4 readouts; η_max 1.52 stands; no polar/radial equations exist
- **Conversion** (Yang & Lan definitions; Nagar–Rezzolla §6 as the Schwarzschild limit):
  - H₂ is the Moncrief / Regge–Wheeler function, so the strain is ∫H₂ dt. The NED correction falls as Σ^{-3/2}.
  - GW energy flux ∝ H₂². The physical EM flux has the same weight in the symmetrised form.
- **Strain check on the solver.** The h_rss of ∫H₂ dt agrees to 0.5–0.7 % with no drift. The L2 criterion
  failed as written because of the d/du = −d/dt sign; with the sign fixed (post hoc) and a trapezoid integral it is
  0.75–1.1 %. This tests linearity and the integration only; the physical relation rests on the derivation.
- **What stands.** For P2's strain-shaped ringdown source, the run output is the strain and the physical energy is
  its ψ̇² measure. So S4's η_max = 1.52 and the receiver-side EM share of 67–70 % stand. The pre-written note had
  said ψ²; it was reversed **post hoc**, on the basis of P2 (budgets cannot decide between the readings).
- **New.** By t_end, 31 % of the initial energy has gone into EM (both universes); the budget closes to 96 %.
- **Failed approach.** The ψ² reading (53–56 %) assumes an H₂-shaped source with a memory offset. Kept for
  comparison only. Arguing from budget closure was also withdrawn, because both measures are conserved.
- **Coupling-off diagnostic** (V₂₂ unchanged): η_max = 1.24. Barrier shape adds about 10 % over "rw" (1.13), and
  coupling adds about 23 % more.
- **S4c.** No polar (ℓ ≥ 2) or radial (ℓ = 0) equations exist for SV with phantom + magnetic NED. 2607.10718's
  radial system uses a gauge that fails at throats. Nearby phantom and charged wormholes are unstable (ℓ = 0;
  polar ℓ = 2).
- **Next if resumed.** The radial ℓ = 0 stability of this background, in a throat-compatible gauge.

## 2026-10-13: arc orientations with a shape S/N cut: the five-field anti excess is not significant at S/N ≥ 50
- **Hypothesis before computing:** if the anti excess (7 vs 2.6, five fields) is a lens-model residual, it survives a
  shape-quality cut; if it comes from noisy shapes near critical curves, it disappears at high S/N.
- **Result (derived, stored `results/arcs/*/arcs.ecsv`):** pipeline isophotal S/N ≥ 50 leaves 3 anti vs 1.33 expected
  (P = 0.15) of 91 background sources in strong shear, aligned fraction 0.92; all three are intrinsic shapes looked at
  earlier. Abell 370 1565 and 1788 drop out at S/N 26.6 and 21.5. Null. Details: `docs/fields/abell370.md`.
- **Next:** blend rejection and F150W vs F444W shapes need the catalogues again; critical-curve curvature.

## 2026-10-13: arc orientations in Abell 370: a 4-vs-0.94 anti excess, all at critical curves or intrinsic shapes
- `lens_consistency.py arcs` (unchanged) with the CANUCS DR1 model and photo-z. Known images: 53/59 within 30°; the
  two beyond 60° are knots of the resolved spiral system 6.
- **Result (derived):** background 47/56 aligned, 4 anti vs 0.94 (P = 0.016). Looked at: 1565 and 1788 on critical
  curves (g_min 1.2, 4.2); 4540 is catalogued image 6.1 (resolved spiral, e 0.56, g 0.23); 340 compact, g 0.26.
  Five fields: 4 anti vs 2.13 (P = 0.17) away from critical curves. Null; D-024 addendum 3 restricts the shape-noise
  null to g_min < 1.
- **Next:** for g_min ≥ 1 anti sources, test the flip against the curve position (z scan across the photo-z range,
  chain spread where `bayes.dat` exists); MACS0717 (no photo-z) or blend rejection + S/N ≥ 50 shapes.

## 2026-10-13: arc orientations vs predicted shear in MACS1149 and MACS0416: null
- **Hypothesis before running:** background arcs align with the published models' stretch (as in SMACS, Abell 2744,
  El Gordo); a real, single background galaxy ≥ 60° off at every z would be a lens-model residual. Ordinary
  explanations: low-S/N or blended shapes, halo segments, members/foreground, near-critical direction flips.
- `lens_consistency.py arcs`, unchanged, with CANUCS DR1 photo-z; F200W pipeline catalogues (manifest sha256s match).
- **Known cases (MACS0416, CANUCS Lenstool):** 63/85 catalogued images within 30°, every e ≥ 0.6 image within 38°.
  MACS1149 (CATS maps) has no image list in `arcs`; its background alignment is the control.
- **Result (derived):** background in strong shear 20/24 (MACS1149) and 26/27 (MACS0416) aligned; anti 2 (0.61
  expected) and 0 (0.43). Both MACS1149 anti looked at (10″ F200W cutouts, background-subtracted): 4236 is diffuse
  light beside a bright galaxy (3.9σ centre peak despite catalogue S/N 112), 1703 sits on a critical curve. Null.
- **Failed approach:** the first contact sheet drew PAs as if north-up; i2d frames are rotated, so lines must go
  through the cutout WCS.
- **Next:** MACS0717 and Abell 370 `arcs` (no DJA; CANUCS covers Abell 370), or blend rejection + S/N ≥ 50 shapes;
  critical-curve curvature test.

## 2026-10-13: arc orientations vs predicted shear in Abell 2744 and El Gordo: no anti-tangential arc unexplained
- **Hypothesis before running:** as in SMACS (D-024), background arcs align with the published models' stretch; a
  real, single background source ≥ 60° off at every z would be a lens-model residual. Ordinary explanations: blends,
  low-S/N shapes, unlensed members or foreground, near-critical direction flips, intrinsic shape.
- `lens_consistency.py arcs` gains a shape-noise expectation for background `anti` sources (σ_ε = 0.26, ASSUMPTION)
  and `--chain-samples` (anti fraction over `bayes.dat` rows) + tests.
- **Known cases (validation):** catalogued images 35/44 (Abell 2744) and 14/17 (El Gordo) within 30°; every
  e ≥ 0.6 image within 31°. The six beyond 60° are compact, low-e or near-critical (looked at).
- **Result (derived):** background in strong shear 33/35 and 10/12 aligned; anti 0 (0.45 expected) and 1 (0.17).
  El Gordo 742: band-consistent thin shape, 60.8° off, anti in 55 % of chain rows: marginal, ordinary. The three
  Abell 2744 anti without photo-z are a z 0.54 spiral-arm knot and two low-S/N segments (cutouts looked at). Null.
- **Next:** critical curves at z_s = 1, 2, 4 and arc curvature against them; or blend rejection + S/N ≥ 50 shapes
  (D-024 "Revisit if") before running `arcs` on the HFF CATS map models.

## 2026-10-13: X2 S4: matter-consistent axial GW equations (Yang & Lan) replace the test-field proxy at one cell; η_max 1.52
- **Input.** arXiv:2603.20594 (axial perturbations of SV under fluid and NED + scalar sources) and arXiv:2112.13198
  (phantom scalar + magnetic NED source). Neither was assessed before.
  - Checked: background Einstein tensor, closed-form potentials, the Regge–Wheeler limit, G⁰₀ − G²₂ = 3Ma²/Σ⁵ > 0
    (no ghost EM channel), a regular throat and radiative ends.
  - Doubts recorded: the "non-Hermitian" asymmetry is a normalisation artifact; the transition-region tables jump
    between modes; the fluid closure is an assumption.
- **Tool.** `scripts/x2_axial_ned.py`, a two-channel leapfrog. It reproduces the published a = 2.20M QNMs (NED
  GR-led 0.246142 − 0.000116i against 0.246140 − 0.000111i; AF control 0.308249 − 0.000263i), converged.
- **Comparison** at GW150914's cell a − 2M = 3.575e-4 (design frozen first, `docs/state/units/X2_s4_matter.md`):
  - About two-thirds of the echo energy leaves in the EM channel.
  - η_max,NED = 1.52, against 0.79 (s = 1) and 1.13 ("rw") from the same code: outside the 0.70–1.24 proxy
    bracket.
  - Timing agrees across the three models tested (one a): spacing 4L + 3.1 to 3.6 %, and a barrier delay of +27
    to +32 M on the first echo. Pulses fall at 0.445 s and 0.910 s, inside the PMW.
- **Failed approach.** A 60 M peak-search grid put every peak at the grid edge for a − 2M ≈ 1e-4 to 1e-3. It was
  fixed with a guard and a test before the comparison.
- **Missing.** The polar sector (where the phantom scalar enters), η, spin and stability. X2 is paused for owner
  review. Next if resumed: AF at the same cell (≤ 0.3 h); then a literature check for the polar equations.

## 2026-10-13: Sunrise system 2 in the four RELICS HST-era models: one source at z ≈ 3; no anomaly
- **Hypothesis before running:** system 2's misfit in M5-G1 (previous entry) is that model's; the RELICS v1 models
  should map 2a–c to one source at some z. If none did, the system 2 identification would be the open question.
- New `scripts/sunrise_relics.py` (byte-range patches of all four RELICS models, ~15 s) + offline SIS tests.
- **Validation:** |μ| from finite differences matches the published z = 6.2 maps (median 0.998–1.000 at 7–12
  images per model); system 1 (known z = 6.2) returns z = 5.8–6.3 in every model.
- **System 2 (derived):** best z = 2.6–3.1 in all four codes, source-plane rms 0.11–0.23″ (M5-G1 ≥ 0.62″), image-plane
  offsets 0.05–0.88″. The published tension is M5-G1's; system 2 closed as ordinary. Details: docs/fields/sunrise.md.
- Not checked: whether RELICS used system 2 as a constraint (if so, partly in-sample).
- **Next:** Sunrise lens-model work is done to the published products' precision; return to TASKS.md "Now".

## 2026-10-13: Sunrise deflection maps read as 0.02″ pixels at z_ref ≈ 3; system 2 fits at no redshift (published tension)
- **Hypothesis before running:** the empirical 0.0223″ per unit is 0.02″ × r(6.2)/r(z_ref) for a round z_ref; if
  so, the z fit recovers system 1's z = 6.2 and gives system 2's redshift. A system 2 that fits at no z would be a
  lens-model residual, to be checked against the published modelling first.
- New `scripts/sunrise_zfit.py` (offline, from `results/sunrise/scofield25_positions.json`) + SIS tests.
- **Unit (derived):** z_ref = 3.02 (16–84 % 2.86–3.16) from the κ match of all 15 images; the paper says the maps are
  scaled to z = 6.2. Reading used as an ASSUMPTION (a global factor cannot be proven a convention).
- **Known case:** system 1's source-plane rms is minimal at z = 5.1; with a level of 1.5× that minimum (ASSUMPTION)
  z = 4.2–6.5 is accepted, so z = 6.2 is recovered (best z biased low by ~1).
- **System 2:** source-plane rms ≥ 0.62″ for every z ≤ 20 (3.2× the level); no redshift reproduces it. Cutouts
  (looked at): one thin arc between two bright cluster galaxies with three similar knots. Scofield+2025 already
  report that models including system 2 are unphysical with ≈ 1″ scatter: a known tension, plausibly missing
  galaxy-scale mass at members Q and R, not a new anomaly. Details: docs/fields/sunrise.md.
- **Failed approach:** the linear image-plane rms |A⁻¹Δβ| diverges for system 1's near-critical images (3.3″ rms at
  z = 6.2), so it cannot drive the z fit; it is reported only.
- **Next:** system 2 in the RELICS HST-era models (D-033), which used it as a constraint.

## 2026-10-13: WHL0137 (Sunrise) image positions in the Scofield+2025 maps: unit check failed; provisional null
- **Hypothesis before running:** the JWST images map to common sources in the published M5-G1 deflection maps to the
  published 0.06″ source-plane rms; a held-out image (1.1a, b, e; 1.7b; system 2) off by ≫ that would be a lens-model
  residual to vet. Checks first: κ and μ from the deflection against the published κ and μ maps.
- New `scripts/sunrise_scofield.py` (byte-range patches; maps are 1.1 GB each) + offline SIS test.
- **Failed validation:** the deflection unit is undocumented and is neither arcsec nor 0.02″ pixels; κ needs
  0.0223″ per unit (ratio 44.3–45.2 across images), the source plane prefers 0.0218″. Not understood (D-033 addendum).
- **Conditional result (derived):** at that scale, constraint rms 0.12–0.17″; held-out 1.1a/b/e and 1.7b within
  0.07–0.24″ of their sources, like the constraints. System 2 offsets 0.9–1.0″ at an assumed z = 6.2 (its z is
  unknown). No position anomaly, provisionally. Details: docs/fields/sunrise.md.
- **Failed approach:** 8 concurrent range requests to Zenodo dropped some (no retry); fixed with back-off and 6 workers.
  Zenodo serves ~180 KB/s and ~6 s per range request from the cloud proxy.
- **Next:** settle the unit (Zenodo record or the authors' MrMARTIAN convention), then predict counter-images and
  fit system 2's redshift; or test the RELICS HST-era models (pinned, D-033) on the JWST image list (held out).

## 2026-10-13: X2 S1–S3: SV-throat echo scattering tool and a GW150914 demonstration (phenomenological; η bounded to order unity)
- Owner review fixed the comparison gate. The basis is Miani et al. Table II (95 % upper limit on the in-window
  h_rss, and the merger-ringdown h_rss), not Table I's 50 % detection-probability sensitivity. The model-to-strain
  link is stated as assumptions P1–P4: a test-field barrier as a GW proxy, a free source amplitude η, spectral
  weights and timing. Design frozen before computing (`docs/state/units/X2_s2_design.md`).
- S1 (`scripts/x2_scatter.py`):
  - peaks found from the actual potentials (s = 1 exactly at R = 3M; s = 0 at R ≈ 3.016M);
  - flux to 4.3e-11; two routes agree to 1.4e-6; Schwarzschild 16πM² limit reproduced;
  - resonance spacing vs 4L: passes 1 % at a − 2M = 1e-4, **fails** at 1e-3 (+1.7 / +3.1 %). Post hoc reading:
    a constant 14–25 M reflection delay, the 4L approximation's bias (owner: not a code failure).
- S2, GW150914 (PMW 0.2–1.2 s), applicable cells only (Δt 0.32–0.58 s, a − 2M ≈ 2.2e-4 to 7.1e-4 M; flat weight
  only 0.45–0.58 s):
  - η_max = 0.70–0.73 (scalar), 0.77–0.80 (EM), 1.18–1.24 (RW-form) with the ringdown weight; 0.99–1.01 flat;
  - between the frozen result classes: η ~ 1 is disfavoured under the scalar and EM proxies (EM with the
    ringdown weight only), not under RW-form, so there is no proxy-independent statement.
- S3 (time domain, double barrier): echo amplitudes match the frequency domain to 2.8 % and 4.5 % (pass); the
  pulse spacing is 4L + 4.1 % (**fails** 1 %; measured post hoc). Shifting S2's pulses by that delay (≈ 16 ms per
  round trip) removes only the Δt = 0.58 s cell; the η_max ranges are unchanged.
- Missing irreplaceable input: SV spin-2 perturbation equations with their matter source, and η. X2 paused for
  owner review; next if resumed: the multi-event η_max map (≤ 2 h).

## 2026-10-13: halo fluxes in the fluxratios tension test: Abell 2744's one flag is A2744-QSO1, ordinary
- **Hypothesis before running:** with halo-subtracted fluxes the floor drops from 0.36 mag to ≈ 0.15 mag, so the
  3σ tension test reaches ≈ 0.5 mag; a flag in both routes goes to vetting, and a null tightens the limit.
- New `fluxratios --halo-phot` (D-032 addendum) → `results/halo_phot/abell2744_all/fluxratios_{band}_{route}.json`;
  `--mu-samples 100`, `--inject-mag 1 1.5 2`.
- **Result (derived):** F444W 15 / 11 images (medfilt / annulus), σ_sys 0.17 / 0.13 mag; F150W 5 / 2, no flag. System 53 is
  the only `flux_tension` in either route: 53a, |z| 4.1 (medfilt), and the 53a/53c pair, |z| 4.9 (annulus; 53b fails
  its S/N there). Injections (simulated) flagged at 1 / 1.5 / 2 mag fainter:
  71 / 57 / 43 % and 89 / 89 / 33 % (losses to S/N grow with the deficit).
- **Vetting (docs/fields/abell2744.md):** system 53 is A2744-QSO1, the triply imaged broad-line "little red dot"
  (Furtak+2025). The offset of 0.55–0.68 mag is achromatic (F150W agrees) and two models agree on the prediction. It is
  1–2σ under the published 20–40 % μ systematics, and the images sample the AGN 855 d and 7043 d apart. Ordinary;
  not a candidate. Contact sheet looked at: 53a and 53c are isolated point-like sources, and 53b is 1.3″ from a bright member.
- **Limitation:** a floor fitted on four other systems understates per-image cluster-model μ systematics, so
  |z| > 3 is a screen. The halo test reaches only five Abell 2744 systems; El Gordo and SMACS have no star-like
  image at S/N ≥ 10.
- **Next:** the lensing-consistency flux-ratio line is exhausted on the current catalogues (no anomaly in three
  fields). Move to the next ordinary lens check with reach: WHL0137 (Sunrise) map model or a fourth cluster with a
  Lenstool chain, or surface-brightness conservation on arcs (D-032).

## 2026-10-13: S0 addendum: X2 (SV-throat echoes vs published post-merger limits) authorized; B conditions stated
- Owner follow-up review; `docs/state/units/S0-2026-10-13.md` addendum.
- 2606.08351 is horizon (black-bounce) greybody, not two-sided throat transmission; that citation is withdrawn.
- GrayHawk-v2: Mathematica only and no licence; its method and benchmark targets are reusable, its code is not.
- No public wormhole echo code exists.
- Minimal check (`scripts/x2_echo_delay.py`, tests; derived): the SV double-peak echo spacing grows as a power
  law in 1 − 2M/a. Miani et al.'s 50–350 ms window corresponds to a − 2M ≈ 4.5e-4 to 2.2e-2 M (60 M☉), so their
  published per-event h_rss thresholds can be mapped without new injections. Direction A therefore qualifies as
  X2 (16 h; premises: non-rotating SV, test-field proxy).
- B (Dai–Stojkovic vs Krasnikov): the ℓ = 0 part is exactly degenerate by Birkhoff in our vacuum region. A time
  variation needs ℓ ≥ 1 transfer through the throat, which neither paper computes. B stays not qualified.

## 2026-10-13: halo_phot elongation gate: El Gordo 1c's flag was an arc; few images stay testable; no anomaly
- **Hypothesis before running:** the El Gordo 1c flag (previous entry) is an arc knot that passed the concentration
  gate; an axis-ratio gate removes it and leaves star-like images. A flag surviving the gate would go to vetting.
- New `axis_ratio` in `scripts/halo_phot.py`: Gaussian-windowed (σ 0.15″) second moments of the median-filtered
  light about an iterated centroid; gate q ≥ 0.80 in F150W and F444W (ASSUMPTION), other bands refused until
  star-calibrated. `scripts/halo_phot_stars.py` → `results/halo_phot/stars/{field}.json` (Gaia DR3, pm or
  parallax > 3σ, unsaturated core).
- **Validation:** Gaia stars (observed positions, derived shapes; 55 F150W / 58 F444W counted, S/N ≥ 10):
  q 5th / 50th percentile 0.91 / 0.95 and 0.94 / 0.97. Known thin arcs El Gordo 1a–c: q 0.57–0.67 in both bands.
  602.1a/b/c (compact, Abell 2744): F150W 0.77 / 0.89 / 0.83, F444W 0.85 / 0.94 / 0.92. Synthetic: point sources
  off a pixel centre q ≥ 0.9; white-noise point sources at S/N 10 q ≥ 0.90 (5th pct; noise reads round, so the gate
  is lenient for faint images); a 5″ thin arc with a knot at f(0.2″)/f(0.4″) = 0.606 gives q 0.22.
- **Re-run (derived; ±1 mag end-to-end injections simulated), tested F150W / F444W:** Abell 2744 2 / 11 (was 18 / 23),
  El Gordo 0 / 0 (was 4 / 2), SMACS 0 / 2 (was 3 / 6). Of the images passing concentration and S/N, most are
  elongated (Abell 2744 28 / 43 F150W). **No flag** in any field; El Gordo 1c is now `untested`. Abell 2744 F444W
  (30, 7, 53, 602.1): MAD σ 0.17 mag, max |r| 0.66; injections recovered −1 mag 12 / 15, +1 mag 7 / 15.
- **Failed approach (caught in /code-review):** moments about the peak pixel made a point source half a pixel off
  centre read q 0.68 on a 0.04″ grid; the first star calibration (q 5th pct 0.80 / 0.85) carried that artefact.
- **Consequence:** halo-subtracted flux ratios are testable only for ~4 Abell 2744 systems in F444W; El Gordo and
  SMACS have no star-like multiple image with S/N ≥ 10. The earlier three-field "scatter 0.10–0.27 mag" was on a
  sample that included arcs.
- **Next:** feed the F444W halo fluxes of the tested Abell 2744 systems into `fluxratios` (`--mu-samples 100`) for a
  tension test with the lower floor; for arcs, a surface-brightness-conserving test (D-032) is the right tool, not
  aperture ratios.

## 2026-10-13: halo_phot over every catalogued image of Abell 2744, El Gordo and SMACS: no anomaly
- Ran `scripts/halo_phot.py` unchanged on all systems (149 / 56 / 60 images), F150W + F444W `_i2d` stamps by S3
  byte ranges (UNCOVER 2561 t003; PEARLS 1176 o241_t012; ERO 2736 o001_t001), end-to-end injections at ±1 mag →
  `results/halo_phot/{abell2744_all,elgordo,smacs0723}/halo_phot.json`. Wall time 10.0 / 4.3 / 3.8 min on 4 cores
  (CPU ≈ 45–50 %: the stamps are read in sequence).
- **Hypothesis before running:** halo subtraction lowers the flux-ratio scatter below `mag_auto`'s σ_sys
  0.32–0.36 mag; a two-route flag on a system's worst image goes to vetting.
- **Tested (both routes, ≥ 2 gated images in the system), F150W / F444W:** Abell 2744 18 / 23, El Gordo 4 / 2,
  SMACS 3 / 6 images. Most images fail the concentration gate (resolved: 67 / 63 of 149 in Abell 2744) or S/N.
- **Scatter (derived, MAD σ of medfilt LOO residuals):** Abell 2744 0.10 / 0.15 mag (max |r| 0.22 / 0.66, 53a
  F444W), SMACS 0.13 / 0.27. On the compact subset this is 2–3× below the `mag_auto` floor (selection caveat: only
  compact images enter).
- **Injections (simulated, testable targets recovered as `flux_outlier`):** Abell 2744 −1 mag 17/18, 19/23;
  +1 mag 10/18, 10/23. El Gordo and SMACS: −1 mag 2/4, 2/2 and 3/6, 4/6; +1 mag 0/6 and 1/12. A dimmed image mostly drops out
  through the S/N or concentration gate (`untested`: 6/6, 9/12, 14/41 in Abell 2744), so the fainter-image limit is
  weak outside Abell 2744.
- **One flag, explained:** El Gordo 1a/1c F150W ±0.97 mag (both routes). Relative to 1a, 1c is 0.85–1.2 mag fainter
  than the `bayes.dat` chain predicts (m1c − m1a predicted −1.19 to −1.54 at 95 %, μ(1c) 13–20; observed −0.34);
  1b (S/N 5.9, gated out) agrees with 1a, and F444W shows the same ratio (achromatic). **Contact sheet:** 1b and
  1c are ~5″ thin tangential arcs and the 0.2″ aperture sits on a knot, so their aperture flux follows surface
  brightness, not μ. Ordinary (resolved arc). **No anomaly** in any of the three fields.
- **Failed approach:** the concentration gate f(0.2″)/f(0.4″) ≥ 0.6 passes thin arcs with a bright knot (1c: 0.607),
  because a line's light grows ∝ r. It needs an elongation gate (e.g. second moments, or f(0.2″)/f(1″) of the
  filtered light, calibrated on PSF stars per band) before the fluxes feed `fluxratios`.
- **Next:** add and validate that elongation gate (PSF stars, 602.1a/b, Abell 2744 7/30/6.1, El Gordo 1c as the
  known arc), re-run the three fields, read stamps in parallel, then feed halo fluxes into `fluxratios`
  (`--mu-samples 100`) for the tension test with the lower floor.

## 2026-10-13: Halo-subtracted photometry explains Abell 2744 602.1c's flux flag; 600.1 is untestable
- New `scripts/halo_phot.py` → `results/halo_phot/abell2744/halo_phot.json`; `tests/test_halo_phot.py`. Compact
  images on UNCOVER `jw02561-o001_t003` F150W/F444W `_i2d` stamps (S3 byte ranges, 6″, D-034 frame offset), two
  routes: stamp minus a 1″ median filter, and a local-plane annulus (0.5–0.7″); r = 0.2″ aperture centred on the
  filtered image's peak. Gates (ASSUMPTIONs, as D-032 / `exotic_screens`): concentration f(0.2″)/f(0.4″) ≥ 0.6 on
  the unfiltered light, |μ| ≤ 20, S/N ≥ 10 with 1.5 × ERR ⊕ blank-aperture scatter. `flux_outlier` needs both routes
  beyond 0.75 mag with one sign, on the system's worst image (as D-032); `one_route` and `untested` are never
  flags. The blank apertures' median (leftover halo) is subtracted from each flux.
- **Hypothesis before running:** the 600.1b and 602.1c flags (`mag_auto` −1.17 / −1.52 mag, D-032 addendum) are
  BCG/ICL light. Halo-subtracted residuals below 0.75 mag confirm it; a residual > 0.75 mag in both routes would go
  to `/vet-candidate`.
- **Validation (simulated):** a Gaussian on a synthetic exponential halo is recovered to +3 % (median filter) / +10 %
  (plane) when the halo's aperture light equals the source's; at 10× it leaves +0.3 mag, always brightening.
- **Ring checks on the real stamps (8 apertures at 1″, same peak search):** blank apertures keep 0.8–11 % of
  602.1c's flux (≤ 0.11 mag of leftover halo; ≤ 0.8 % off-halo on 602.1a/b); injected PSF-width Gaussians are
  recovered at 95–100 % for the compact images. Off-halo control: 602.1a/b residuals −0.04 to +0.18 mag.
- **Result (derived; LOO residual against 602.1a/b):** 602.1c F444W −0.27 (filter) / −0.32 (plane), statistical
  ± 0.01 / 0.02 (ERR, blank scatter and the whole halo correction; D-032's σ_sys 0.32–0.36 mag is *not* included and
  the class uses the fixed 0.75 mag threshold); was −1.52 in `mag_auto`. `consistent`. In F150W 602.1c fails the
  concentration gate on the halo (0.57). The flag is halo light. **600.1 is untestable:** a resolved spiral on the
  contact sheet, all four images fail the concentration gate except 600.1b, which then has no sibling; its
  `mag_auto` flag stays explained only by the contact sheet (D-032 addendum). 602.1d: |μ| = 25. **No anomaly.**
- **End-to-end injection (`--inject-mag`, real stamps; each image's own filtered light within 0.45″ scaled, so the
  real profile and halo stay; whole chain re-run):** 602.1a/b are flagged at ±1 mag in both bands, and at +0.75 mag
  in F444W (F150W stays `consistent`, at the threshold). 602.1c is flagged at −1 mag in both bands (brightening
  also lifts it over the F150W concentration gate) but not at +1 mag: its −0.27 mag offset leaves +0.73. So the null
  for 602.1 is: no image deviates by more than 0.75 mag from its siblings beyond the measured offsets, in F444W.
  600.1 and 602.1d stay untested at every amplitude. A Gaussian subtracted at the catalogue position failed as a
  dimming injection (every target left `untested`: it does not match the real profile).
- **Failed approach:** concentration on the median-filtered image let the spiral's nuclei pass (0.63–0.81); the
  filter strips extended light. It is now measured on the unfiltered light.
- **Final review (takeover session 0117Bja4):** each stamp is filtered once for both routes and the blank
  apertures, each injected ring stamp once; ring checks skip images the gates drop; all measurements, injections
  included, run in one process pool (one BLAS thread per worker; base + three amplitudes 47 s on 4 cores, was
  6.5 min); fluxes are divided by the recovery, and a route keeping < 50 % of an injected image is dropped; the
  error adds the whole halo correction, since the ring 1″ out can misjudge the leftover halo at the cusp; injection
  records keep the class before and after; `classify()` is tested; `--systems` is stripped and unknown ids fail.
  Injection caveat: scaling the filtered light also scales the pixel noise and any neighbour within 0.45″. Kept on purpose: the fixed 0.75 mag threshold (as D-032) and
  a NaN policy stricter than `lc.aperture_snr` (a masked pixel cannot be halo-filtered).
- **Limits:** one system tested; one UNCOVER epoch and reduction. Stamps are read in sequence, then measured in
  parallel. NaN pixels are filled with the stamp median before filtering (none here).
- **Next:** run `halo_phot.py` over every compact image of Abell 2744, El Gordo and SMACS (SMACS 1.1, 6.x), then
  fold the halo-subtracted fluxes into `fluxratios` as a photometry route so `--mu-samples` refits σ_sys.

## 2026-10-13: X1 closure scope corrected; S0 screen finds no qualifying unit; routine prompt 2026-10-13.1
- Owner review. X1's closure separates three things:
  - the cross-checked tool, within its shared premises (point spot, primary + secondary images, bolometric g⁴,
    common emission assumption);
  - the exploratory, assumption-dependent position comparison;
  - the S2 gate observable (flux contrast), which no data tested.
  No candidate, no new-physics evidence. Gates and numbers are unchanged.
- S0 screen (`docs/state/units/S0-2026-10-13.md`, about 1.2 h):
  - throat transmission is already published, and mapping LVK echo results needs injections beyond the 20 h cap;
  - a far-side mass is invisible to S2-like probes: monopole degenerate with M, quadrupole ≤ 1.1e-4, suppressed
    ~1e-10 in the orbit;
  - the lensing arrival order is settled by theorem and mimicked by substructure.
  None qualifies; the high-risk share idles.
- Rules trial-2026-10-13.2:
  - T1 outranks a stale routine-search claim when an authorized unit has an open stage;
  - T7 sets the S0 screen format;
  - "ordinary physics cannot mimic" now means comparison with explicit ordinary models over plausible ranges.
- Routine prompt 2026-10-13.1: one saveable unit step per run (tool, theory or data); the lens-only restriction
  removed. **Repo updated; the routine's own prompt awaits the owner's paste.**

## 2026-10-13: `fluxratios` with `bayes.dat` μ errors: model statistics do not explain the flux scatter
- `lens_consistency.py fluxratios --mu-samples N`: μ at N random chain rows (`chain_magnifications`, family
  redshifts included); each image's model error is half the 16–84 % range of its leave-one-out residual. `resid_z`
  adds a scatter floor (median |z| = 0.6745, a pair counted once) fitted on the *other* systems, so an anomaly
  cannot raise its own floor; `flux_tension` = a system's highest-|z| image with |z| > 3 (ASSUMPTION);
  `mu_unstable` (whole system) where > 16 % of rows exceed |μ| = 20 (none here). `--inject-mag` measures recovery
  by making each consistent image fainter in turn (flux and S/N scaled together). Tests: synthetic chain, μ spread, critical-curve rows, floor.
- **Hypothesis before running:** chain μ errors are ≪ the 0.5–0.6 mag rms, so the scatter is photometric
  (halo light, faint counterparts) and the flags do not change. If instead the model error were ≳ 0.3 mag, the
  fixed 0.75 mag threshold would be too strict near critical curves and the next step would be model systematics.
- **Result (derived; 100 rows each, seed 1; DJA files sha256-matched, downloaded to the session disk and deleted):**
  model error median / max 0.04 / 0.14 mag (Abell 2744, 34 images) and 0.07 / 0.21 mag (El Gordo, 23); scatter
  floor σ_sys 0.36 and 0.32 mag. The `flux_tension`s are exactly the old `flux_outlier`s: Abell 2744 42.1a/c,
  600.1b, 602.1c (|z| 6.2–8.5; halo light, Lyα counterpart) and El Gordo 18b/c (|z| 4.6; faint counterpart), all
  explained on the contact sheets (D-032). **No anomaly; the hypothesis holds.**
- **Injection-recovery (simulated, 30 / 21 target images; flux and S/N scaled together, so a faint image can drop
  below S/N 10 and count as missed):** an image made Δ mag fainter is flagged in 17 / 27 / 60 / 93 % (Abell 2744)
  and 24 / 24 / 52 / 62 % (El Gordo) of cases at Δ = 0.5 / 1 / 1.5 / 2 mag. A demagnification by ≳ 4 (Δ ≥ 1.5)
  is found in about half to most tested images; Δ ≲ 1 mag hides in the scatter. S/N losses at Δ = 2: 2 / 30
  (Abell 2744) and 8 / 21 (El Gordo). 14 / 30 and 18 / 21 targets are pair members, where a hit detects the system
  but cannot say which image is faint.
- SMACS 6.3 cannot be rechecked this way: the Mahler+22 release pins no chain.
- **Limits:** chain errors are statistical within one model; model-to-model differences are larger (D-044).
- **Next:** photometry, not the model, limits this test: ICL-aware totals (DJA `_phot_apcorr`, local-background
  Kron) for images on halos, then re-run with `--mu-samples`.

## 2026-10-13: Abell 2744 flux-ratio and colour test on DJA v7.2: no anomaly
- `lens_consistency.py --model abell2744-bergamini23 fluxratios` on DJA v7.2 `fix_phot` + `eazypy.zout` (233.5 +
  63.6 MB, sha256 matched; downloaded to the cloud session disk and deleted after use; reason in D-032 addendum).
- **Hypothesis before running:** no image deviates beyond the D-032 thresholds (0.75 mag flux, 0.3 mag colour;
  ASSUMPTIONs) once instrument, blending and counterpart explanations are tested; a flag that survives the contact
  sheet and a second photometry route would go to `/vet-candidate`.
- **Result (derived):** frame offset negligible (−0.010″, −0.010″). 65 / 149 images matched; 34 flux-tested (median
  |r| 0.17 mag, rms 0.60), 35 colour-tested (median 0.06, rms 0.28). 7 flags, all ordinary on the contact sheet:
  BCG/ICL halo light (600.1b, 602.1c), a faint Lyα counterpart whose sign flips between DJA and forced i2d
  photometry (42.1a/c), a resolved giant-arc aperture (3.1b/c colour) and a near-blend (34.1c). **No anomaly.**
  Details: docs/fields/abell2744.md.
- **Failed approach:** forced photometry on one UNCOVER tile (2561 t003) cannot settle continuum-faint Lyα images.
- **Next:** an ICL-aware total (local-background Kron or the DJA `_phot_apcorr` columns) for images on cluster-galaxy
  halos in all three fields; `bayes.dat` μ errors in `fluxratios` (then a χ² in place of fixed thresholds).

## 2026-10-13: X1 S4b: the tracer's image weights pass an independent check; X1 closed
- New `scripts/x1_s4b.py` → `results/x1/s4b_image_weights.json`; `tests/test_x1_s4b.py`; `x1_tracer.image_mag`
  extracted (no numeric change). Gate written in docs/state/units/X1.md before the run (commit 75bd7b7). Hours count
  as verification under T1.
- **Hypothesis before running:** the tracer's weight b / (|dΦ/db| cos α sin γ) is the point-source magnification.
  If not, S2/S3 weights are wrong and S3 is not reported as a finding.
- **Method:** the magnification by another route: 1/|det J| of the screen → perpendicular-offset map, built in 3-D
  from Hamilton-integrated swept angles and the photon direction at the emitter (cos α from p_r).
- **Result (derived):** far side a = 2.1 and 5 M (the gate): ≤ 5.0e-7 relative (tol 5e-3); a = 3, 7 M ≤ 4.9e-6; ordinary
  a = 0 near side (primary and secondary) ≤ 3.4e-5, the largest at the fold. Weak-field control: magnification 1 to
  1.5e-7. The tracer's weight is per unit source solid angle, i.e. R² × per unit cross-section. This constant cancels
  at fixed R, so no S1–S4 number changes.
- **X1 closed** (3.7 / 20 h). Conditional results and the reopen inputs: docs/state/units/X1.md "Closure". Both
  trial slots are now empty; a new unit needs an S0-passed direction (T3).
- **Next:** routine searches are under their weekly share (T1), but MOA-II has no field left except gb13 (waits on
  the owner, #145). Pick the next routine unit from TASKS.md "Continuing" item 1 and the `w3-survey` skill (which
  notes the KMTNet/OGLE EWS limits; `scripts/w3_gaia.py` covers only `vari_microlensing`).

## 2026-10-13: X1 S4: the model-free disk bound is consistent with S3, but S3's Δχ² rests on the traced centroid
- New `scripts/x1_s4.py` → `results/x1/s4_disk_bound.json`; `tests/test_x1_s4.py`. Gate written in
  docs/state/units/X1.md before the run. The S4 hours count as verification under T1.
- **Hypothesis before running:** if S3's far-side disfavour follows from b < b_c alone, a lower bound with each
  epoch's centroid free in the b_c disk keeps Δχ² ≥ 9 (ASSUMPTION). If it does not, S3 depends on the tracer and
  orbit model.
- **Result (derived; same flare as S3 and #172, exploratory):** the convex SOCP bound is ≤ S3's far χ² for every
  a and offset cap (consistent). The model-free Δχ² against S3's best ordinary fit (18.6) is +6.8 (a = 2.1, 3 M),
  0.0 (5 M), −7.5 (7 M) and −12.9 (9 M). No a reaches 9. S3's +23 to +37 therefore needs the traced far-side
  centroid (max 1.6–2.7 M, far inside b_c ≈ 5.2 M).
- A second route for the far-side swept angle (Hamilton's equations, DOP853) matches the tracer to ≤ 3.3e-11 on 12
  rays.
- **Failed approach:** einsteinpy 0.4.0 FANTASY works in Schwarzschild (≤ 1.2e-6) but its fixed step is unstable
  at the SV throat (H → 10⁴; Δφ not converged at step 0.0025).
- **Next:** the far-side image-weight check (magnification by finite differences vs the tracer), then close X1.

## 2026-10-13: W3 MOA-II gb14 (twenty-first field): 0 survivors; combined limit over twenty-one fields
- **Result:** 103 flags, 0 survive; combined Γ₉₅ ≈ 1.2 × 10⁻⁸–7.4 × 10⁻⁸ W3 events per star per year (was
  1.4–7.8 × 10⁻⁸; gb14 tightens every cell by 5–10 %). Details: docs/exotic_limits.md "gb14" and the combined table.
- Hypothesis before running: 0 survivors under CHAIN_VERSION 2026-10-10.1; a survivor with caustic spikes in the data
  goes to `/vet-candidate` before injections. Plan: stream the tar in tracked pre-screen chunks across routine slots,
  then fit, vet, inject (6 chunks), limit, combine.
- Run 1 (cloud, 22:29–22:47Z, `prescreen --procs 4`, `timeout 1080`): chunks 1–16 of 50 tracked (0 errors, CPU 99 %;
  streamed, nothing stored; 16 chunks in 18 min). Per-chunk `wall_time_s` rose 70 → 159 s within the run, the
  overlapping-window artifact noted in TASKS.md.
- Run 2 (cloud, 2026-10-13 00:09–00:46Z, two `prescreen --procs 4` calls, `timeout 1560` then `600`): chunks 17–45
  tracked (0 errors, CPU 99 % until the last minutes). The archive throughput fell from ~65 MB/s to ~0–4 MB/s with
  4 retries near 00:44Z, so chunk 46 did not finish before the timeout. **Next:** `prescreen --procs 4 --chunks 46-50`,
  then merge-prescreen, fit (4 chunks; `.[dev,cloud,mulens]`), merge-chunks, vet, inject (6 chunks), limit, combine.
- Run 3 (cloud, 01:09–01:13Z, `prescreen --procs 4 --chunks 46-50`, `timeout 1500`): throughput back at 68–89 MB/s;
  chunks 46–50 tracked, 0 errors. **Pre-screen complete: 182,302 light curves, 2,066 shape passes, 109 not at a
  shared epoch** (merge ran automatically). Fits (`fit --procs 4 --chunk K/4`, sequential, detached, 01:14–01:42Z):
  **109 fits, 0 errors, 103 flags** (28/24/25/26 per chunk; 317–552 s each). `vet --procs 4` (01:42–01:45Z):
  **0 of 103 survive** (`results/w3_moa/vetting_gb14.json`). Quiet χ²/dof 95th percentile 3.37; CDS XMatch worked
  (Gaia DR3 RP for 55 / 109). Leading failures: exotic_in_domain 97, feature_bracketed 88, exotic_chi2_cap 78,
  repeated_deficit 47. Contact sheet inspected: box and eclipse dips, slow variables, noisy flat curves; caustic
  spikes only in the models. **Next:** `inject --procs 4 --chunk K/6` (6 chunks, ~18 min each), `merge-inject --n 6`,
  `limit`, `combine`, then docs/exotic_limits.md and the PR.
- Run 4 (cloud, 02:08Z–, session 01BiV92m): fresh clone, so `merge-prescreen` and `merge-chunks --n 4` first
  (the first `inject` calls failed in 2 s on the missing `prescreen_gb14.ecsv`; no output written). Injection chunks
  (`inject --procs 4 --chunk K/6`, sequential, detached): 1/6 367 in 641 s, 0 errors, W3 333 → 48 flagged → 22
  recovered, PSPL 0 / 34; 2/6 367 in 479 s, 0 errors, W3 333 → 28 → 17, PSPL 0 / 34.
  3/6 367 in 497 s, 0 errors, W3 334 → 28 → 10, PSPL 0 / 33; 4/6 367 in 617 s, 0 errors, W3 334 → 41 → 23,
  PSPL 0 / 33. Chunks 1–4 together: W3 1,334 → 145 flagged → 72 recovered; PSPL 0 / 134 flagged. Chunk 5 was
  stopped at the 40-min mark (nothing written). **Next:** `merge-prescreen`, `merge-chunks --n 4` (fresh clone),
  `inject --procs 4 --chunk 5/6` and `6/6`, `merge-inject --n 6`, `limit`, `combine`, docs/exotic_limits.md, PR ready.
- Run 5 (cloud, 04:10–04:45Z, session 01B15J5c): `merge-prescreen`, `merge-chunks --n 4`, injection chunks 5/6 (366 in
  496 s, 0 errors, W3 333 → 31 flagged → 12 recovered, PSPL 0 / 33) and 6/6 (366 in 449 s, 0 errors, W3 333 → 30 → 13,
  PSPL 0 / 33). All six: W3 2,000 → 206 flagged → 97 recovered (47 %); PSPL 0 / 200 flagged. `limit` first failed:
  the tracked `results/w3_moa/vetting_gb14.json` is the summary, without the `chain` stamp `limit` checks, and the
  full record lived only in run 3's `data/derived/`. Re-ran `vet --procs 4` (~3 min): 0 survivors, funnel identical
  to run 3. Contact sheet re-inspected: same picture. gb14 limits 1.3 × 10⁻⁷–1.4 × 10⁻⁶ per cell; `combine` over 21
  fields. MOA-II is now complete except gb13 (#145, owner). **Next (W3):** gb13 after the owner's decision; another
  MOA-II-like field moves the combined limit by ≲ 10 %, so new W3 coverage should come from a new survey.

## 2026-10-13: Trial tally covers every run (trial-2026-10-13.1)
- Acceptance check of D-083:
  - The cloud routine run of 23:10Z claimed X1 with `rules: trial-2026-10-12.1 prompt: none`. The routine's own
    prompt text is not updated yet; the owner pastes it.
  - That run passed S1 (#176) and opened S2 (#177).
  - The 00:09Z run worked W3 gb14 and did not log hours, so T1 shares could not be computed.
- Fix: every run appends a tally line, routine searches included; gb14's run back-filled.

## 2026-10-13: X1 S3: the 2018-07-22 GRAVITY flare disfavours the far-side spot under the premises (exploratory)
- New `scripts/x1_s3.py` → `results/x1/s3_avg_fit.json`; `tests/test_x1_s3.py` (exact recovery of injected
  tracks). Input: CDS J/A+A/618/L10 `avg.dat`.
- Astrometric χ² on 20 values (grid then bounded Powell; far side also given η 0.6–1): ordinary 18.6–28.4 (best
  p ≈ 0.14); far a = 2.1–7 M 45–57; a = 9 M 33–37. Δχ² far − ordinary: +23 to +37 (a ≤ 7), +14 (a = 9).
- Static offsets beyond the pre-written 20 µas lower the far χ² slowly (not bracketed); Δχ² ≥ 26 up to 400 µas.
- Input: `data/manifests/x1_gravity_flare.ecsv` (sha256 pinned). Downloaded to `$JWST_ANOMALY_DATA/raw/x1/`
  and deleted after use in cloud sessions.
- Same flare as #172, so this is not independent of the S0 bound. Flux contrast was not testable (one flare, half
  an orbit). Details: `docs/state/units/X1.md`.
- **Next:** X1 S4, a model-free disk-bound fit as an independent check.

## 2026-10-13: X1 S2: no astrometric far-side discriminant; a fragile flux-contrast one at a = 2.1 M (gate passed formally)
- `scripts/x1_tracer.py`: `spot_track` takes sub-Keplerian η and light-travel delays per image. Validated: primary
  delay = Roemer 2R sin i to < 1 % at R = 400 M.
- New `scripts/x1_s2.py` → `results/x1/s2_discriminant.json` (derived, conditional geometry). i = 155°, P = 60 ± 3 min.
  Ordinary: a = 0, η 0.6–1. Far: a ∈ {2.1, 3, 5, 7, 9} M. Both get static light 0–0.5 at 0–20 µas.
- Loop radius: far 2–27 µas vs ordinary 14.5–48 µas. The gaps are ≤ 0.92 × the 9 µas precision (ASSUMPTION), so
  **astrometry does not discriminate.**
- Flux contrast ln(F_max/F_min): far ≤ 0.33 at a = 2.1 vs ordinary ≥ 0.43, i.e. 1.98 × the assumed 0.05 precision;
  0.43 × at a = 3; none above. Peak count: 1 everywhere.
- The pass is fragile. It closes at static fraction 0.56, needs orbital modulation separated from intrinsic flare
  evolution, and holds for a point spot only. Details and the frozen S3 prediction: `docs/state/units/X1.md`.
- Design choice: L is measured about the track's own mean (free centre), because GRAVITY's centre is degenerate
  with the orbit radius.
- **Next:** X1 S3 on `avg.dat` against the frozen prediction, then S4 (independent check of the a = 2.1 contrast).

## 2026-10-12: X1 S1: Simpson–Visser ray tracer validated (gate passed)
- New `scripts/x1_tracer.py` (quadrature in u = 1/R with regularised turning-point and throat endpoints, plus an
  independent ODE route) and `tests/test_x1_tracer.py`; `scripts/wh_farside.py` and its tests copied from #172.
- Checks (`results/x1/s1_checks.json`, derived): Schwarzschild deflection vs Darwin ≤ 4e-15; b_c(a) by ODE rays 4.8e-7;
  swept angles (far side, and near side above b_c) quadrature vs ODE ≤ 8e-11. Numbers and limits: `docs/state/units/X1.md`.
- Finding for S2: GRAVITY 2023 App. A's r_o = R + 1 M is the primary image (−0.05 M face-on in our tracer); the
  all-light centroid is ~0.6 M smaller, so the ordinary model must be traced, not (R + 1) θ_g. Far-side spots at
  R = 9 M give centroid radii ≈ 4–14 µas and 2.6–17 % of the near-side flux (a = 2.1–7 M).
- Failed approaches: weighting images by b|db/dΦ|/sin γ alone (it omits 1/|cos α_e| and zeroes the image at the
  fold where the ray is tangent at the emitter); a uniform b grid misses roots next to the fold (now refined).
- **Next:** X1 S2 (state file "Next action").

## 2026-10-12: Exploration trial adopted (D-083): split WIP counters, pending index, unit state files; X1 F0 done
- The owner approved #174. New files: `docs/state/trial.md` (rules T1–T6, version trial-2026-10-12.1, weekly
  tally, ledger), `docs/state/pending.md` (provisional conclusions of open needs-human PRs) and
  `docs/state/units/X1.md`.
- Changed: research-cycle step 1 (trial dispatch) and step 3 (two WIP counters), the CLAUDE.md WIP paragraph,
  TASKS.md (trial section; S1 "done; do not redo") and `docs/cloud-routine-prompt.md` (prompt version 2026-10-12.1).
- **The routine's own prompt text is not updated from a session.** The repository files above take effect on the
  routine's next run, because it reads them on every run.
- X1 F0 feasibility: see the unit file.

## 2026-10-12: W3 MOA-II gb10 (twentieth field): 0 survivors; combined limit over twenty fields
- Hypothesis before running: 0 survivors under CHAIN_VERSION 2026-10-10.1; a survivor with caustic spikes in the data
  goes to `/vet-candidate` before injections. gb10 is a 220.9 GB tar, streamed (nothing stored) in 52 tracked
  pre-screen chunks over three cloud runs (16:10–18:13Z), then fits, vetting and injections over four more.
- Pre-screen: **181,820 light curves, 1,469 shape passes, 108 not at a shared epoch**; 0 errors. Run 2 did 20 chunks
  in ~30 min, but CPU fell from 99 % to 51–88 % for chunks 42–49 (below the 70 % rule: an I/O stall, recovered at the
  next run). The tracked per-chunk `wall_time_s` / `mb_per_s` overstate the stall: each chunk's clock starts at its
  first finished segment, so the windows overlap (follow-up in TASKS.md).
- Fits: 108, 0 errors, **101 flags**. `vet --procs 4`: **0 of 101 survive** (`results/w3_moa/vetting_gb10.json`;
  reproduced exactly by run 7). Quiet χ²/dof 95th percentile 3.73; CDS XMatch worked (Gaia DR3 RP for 63 / 108).
  Contact sheet inspected: box and eclipse dips, variables, noisy flat curves; caustic spikes only in the models.
- Injections (six tracked chunks, 0 errors): **W3 2,000 → 209 flagged → 99 recovered** (47 % of flagged kept by
  vetting); PSPL controls 0 / 200 flagged. gb10 limits (`results/w3_moa/limits_gb10.ecsv`; N_s = 2.15 × 10⁷,
  Nunota et al. 2024 Table 1, 70/80 subfields): Γ₉₅ ≈ 1.7 × 10⁻⁷–3.3 × 10⁻⁶ per star per year; both t_E = 300 d
  cells rest on 1 recovery, t_E = 3 d, ρ = 0.01 on 2.
- **Combined over twenty fields:** Γ₉₅ ≈ 1.4 × 10⁻⁸–7.8 × 10⁻⁸ per star per year (nineteen fields: 1.5 × 10⁻⁸–
  8.0 × 10⁻⁸); gb10 holds 6.4 % of the stars and tightens the cells by 2–9 %; counting uncertainty 6–15 %; max field
  share 0.204 (gb9, t_E = 300 d, ρ = 0.01). Details: docs/exotic_limits.md "W3 in MOA-II".
- Failed approaches: a `while pgrep -f "chunk 3/4"` wait loop matches its own command line and never exits (wait on a
  PID or an output file); a fresh venv needs the `mulens` extra before `fit`; run 7 ran `vet --procs 4` beside
  `inject --procs 4` on 4 cores (oversubscribed, against scripts/CLAUDE.md; both finished, avoid it).
- **Next:** gb13 waits on the owner (#145); then gb14, the last large Nunota et al. 2024 field (TASKS.md).

## 2026-10-12: W3 MOA-II gb9 (nineteenth field): 0 survivors; combined limit over nineteen fields
- Hypothesis before running: 0 survivors under CHAIN_VERSION 2026-10-10.1; a survivor with caustic spikes in the data
  goes to `/vet-candidate` before injections. Run 1 (cloud, 07:35Z): tracked pre-screen chunks only.
- Run 1 (cloud, 07:35–07:47Z): pre-screen chunks 1–13 of 101 tracked (~65–120 s per chunk on 4 procs, CPU 99 %,
  35–64 MB/s; streamed, nothing stored). The run hit its time limit during chunk 14, which reruns.
- Run 2 (cloud, 08:12–08:44Z, `prescreen --procs 4 --chunks 14-101`): chunks 14–50 tracked (37 chunks, 0 errors,
  0 retries, mean 117 s per chunk wall, 30–82 MB/s); stopped by PID at the slot's end (chunk 51 reruns).
- Run 3 (cloud, 09:10–09:40Z, `prescreen --procs 4 --chunks 51-101`): chunks 51–70 tracked (20 chunks, 0 errors,
  at most 1 retry, mean 141 s per chunk wall, ~40–50 MB/s, CPU 99 %); stopped by `timeout` at the slot's end
  (chunk 71 reruns).
- Run 4 (cloud, 10:10–10:39Z, `prescreen --procs 4 --chunks 71-101`): chunks 71–89 tracked (19 chunks, 0 errors,
  at most 1 retry, mean 131 s per chunk wall, ~45–50 MB/s, CPU 99 %); stopped by `timeout` at the slot's end
  (chunk 90 reruns).
- Run 5 (cloud, 11:13–11:24Z, `prescreen --procs 4 --chunks 90-101`): chunks 90–101 tracked (12 chunks, 0 errors,
  0 retries, mean 60 s per chunk wall, ~65–78 MB/s, CPU 99 %); the pre-screen is complete: **237,123 light curves,
  2,363 shape passes, 111 not at a shared epoch** (to fit; gb5 had 126). Fits (`fit --procs 4 --chunk K/4`):
  chunk 1/4 28 fits, 0 errors, 28 flags, 466 s; chunk 2/4 28 / 0 / 28, 547 s; chunk 3/4 28 / 0 / 28, 585 s (run 5 ended 11:53Z).
- Run 6 (cloud, 12:20–13:25Z): `merge-prescreen` (a fresh clone needs it before `fit`), fit chunk 4/4: 27 fits,
  0 errors, 25 flags, 1,074 s; `merge-chunks --n 4`: 111 fits, **109 flags**. `vet --procs 4`: **0 of 109 flags
  survive** (`results/w3_moa/vetting_gb9.json`); 100 / 109 took 21 s, the last 9 held one worker for ~20 min.
  Quiet χ²/dof 95th percentile 3.42 (1,927 quiet curves); CDS XMatch worked (Gaia DR3 RP for 61 / 111). Contact sheet
  inspected: box and eclipse dips, variables, noisy flat curves; caustic spikes only in the models.
  Failed approach: `photometry.fetch_catalog` (urllib) got the 97 MB MOA `metadata.ipac.tar.gz` truncated twice
  (different sha256 each time); `curl -C - --retry 3` into `$JWST_ANOMALY_DATA/cache/external/b339ec176933_metadata.ipac.tar.gz`
  fetched it intact on the first try.
- Run 7 (cloud, 13:24–14:10Z; started before run 6's late push): `vet --procs 4` reproduced the tracked record
  exactly (same funnel, 0 / 109). The slow flag is gb9-R-10-3-499239 (ΔBIC −413, PSPL t_E ≈ 888 d, 25k epochs):
  ~16 min in the binary-lens fit (py-spy: VBBL over every epoch, since the ±3 t_E window spans the whole light
  curve); it passes binary source and binary lens (ΔBIC −378) and fails `exotic_feature_sampled` (11 epochs on
  2 nights inside the feature). `vet` now logs each flag's id, so a stall names its flag.
- Injection chunk 1/6 (`--procs 3` beside the slow vet): 367 in 1,234 s, 0 errors; W3 333 → 35 flagged →
  **18 recovered**; PSPL 0 / 34 flagged. Chunk 2/6 (`--procs 4`): 367 in 502 s, 0 errors; W3 333 → 23 → **8**;
  PSPL 0 / 34.
- Run 8 (cloud, 14:19Z–): `merge-prescreen`, `merge-chunks --n 4` (111 fits). Injection chunk 3/6 (`--procs 4`):
  367 in 1,097 s, 0 errors; W3 334 → 33 → **16**; PSPL 0 / 33. Chunk 4/6: 367 in 681 s, 0 errors; W3 334 → 29 →
  **12**; PSPL 0 / 33. Chunk 5/6 was started at 14:51Z and stopped at the slot's end (reruns). Chunks 1–4 so far:
  W3 1,334 → 120 flagged → 54 recovered; PSPL 0 / 134.
- Run 9 (cloud, 15:10–15:50Z): `merge-prescreen`, `merge-chunks --n 4`. Injection chunk 5/6: 366 in 768 s, 0 errors;
  W3 333 → 29 → **12**; PSPL 0 / 33. Chunk 6/6: 366 in 556 s, 0 errors. All six: **W3 2,000 → 178 flagged → 81
  recovered** (46 % of flagged kept by vetting); PSPL controls 0 / 200 flagged. `limit` needs the full derived
  vetting record, which a fresh clone lacks, so `vet --procs 4` ran a third time (770 s, the slow flag last):
  identical record, 0 / 109.
- gb9 limits (`results/w3_moa/limits_gb9.ecsv`; N_s = 3.33 × 10⁷, Nunota et al. 2024 Table 1): Γ₉₅ ≈ 1.4 × 10⁻⁷–
  1.05 × 10⁻⁶ per star per year; the t_E = 300 d, ρ = 0.1 cell rests on 2 recoveries.
- **Combined over nineteen fields:** Γ₉₅ ≈ 1.5 × 10⁻⁸–8.0 × 10⁻⁸ per star per year (eighteen fields:
  1.5 × 10⁻⁸–8.6 × 10⁻⁸); gb9 holds 11 % of the stars and tightens the cells by 6–21 %; max field share 0.21 (gb9,
  t_E = 300 d, ρ = 0.01, 8 recoveries) (`results/w3_moa/limits_combined.ecsv`, docs/exotic_limits.md).
- Failed approach: `limit` on a fresh clone stops on the missing derived `vetting_gb9.json`; the tracked record
  lacks the `chain` and `fit_errors` keys `limit` checks, so a field finished across runs re-runs `vet` before `limit`
  (budget ~15 min here for gb9's slow binary-lens flag).
- **Next:** gb13 after the owner decides #145; then gb10, gb14.
## 2026-10-12: W3 MOA-II gb6 (eighteenth field): 0 survivors; combined limit over eighteen fields
- Hypothesis before running: 0 survivors under CHAIN_VERSION 2026-10-10.1; a survivor with caustic spikes in the data
  goes to `/vet-candidate` before injections. gb6 is small (7.7 GB tar). Stacked on `claude/w3-moa-gb5` (#168).
- Run 1 (cloud, 05:35–05:42Z): pre-screen 2 tracked chunks in 111 s (0 errors): **35,948 light curves, 800 shape
  passes, 32 not at a shared epoch**. Fit chunk 1/1: 32 fits, 0 errors, 32 flags, 118 s. `vet`: **0 of 32 flags
  survive** (`results/w3_moa/vetting_gb6.json`); CDS XMatch worked (Gaia DR3 RP for 29 / 32). The quiet χ²/dof 95th
  percentile is 5.5, against 2.1–2.6 in gb4 and gb5 (a noisier small field; the χ² cap scales with it). Contact sheet
  inspected: box and eclipse dips, single low points, variables; caustic spikes only in the models.
- Run 2 (cloud, 06:13–06:40Z): `merge-prescreen`, `merge-chunks --n 1`, `vet` reproduced the tracked record (0 / 32,
  identical funnel). Injections in four tracked chunks of 550 (0 errors, ~270–355 s each, ~2 per s: gb6's small
  carrier set is fast): W3 2,000 → 135 flagged → **80 recovered** (59 %); PSPL controls 0 / 200 flagged.
- gb6 limits (`results/w3_moa/limits_gb6.ecsv`): Γ₉₅ ≈ 8.0 × 10⁻⁷–1.0 × 10⁻⁵ per star per year. gb6 is not in
  Nunota et al. 2024 Table 1, so N_s = 6.76 × 10⁶ comes from the Cut-0 object model (ASSUMPTION, as for gb22). The
  t_E = 3 d, ρ = 0.01 cell rests on 1 recovery.
- **Combined over eighteen fields:** Γ₉₅ ≈ 1.5 × 10⁻⁸–8.6 × 10⁻⁸ per star per year (seventeen fields:
  1.6 × 10⁻⁸–9.0 × 10⁻⁸); gb6 holds 2.4 % of the stars and moves the cells by 0–5 %; max field share 0.17
  (`results/w3_moa/limits_combined.ecsv`, docs/exotic_limits.md).
- Failed approach: a `setsid nohup ... &` queue started from a cloud Bash call died when the call returned (after
  `merge-chunks`, during `vet`). Use the harness's background task (`run_in_background`) for long queues there.
- **Next:** gb13 after the owner decides #145; then gb9, gb10, gb14. This PR is stacked on #168, #165 and #160:
  merge those first.

## 2026-10-11: W3 MOA-II gb5 (seventeenth field): 0 survivors; combined limit over seventeen fields
- Hypothesis before running: 0 survivors under CHAIN_VERSION 2026-10-10.1; a survivor with caustic spikes in the data
  goes to `/vet-candidate` before injections.
- gb5 is the largest field so far: **119 pre-screen chunks** of 4 GiB (~485 GB streamed, never stored; gb4 had 66).
  Run 1 (cloud, 22:11–22:41Z, `prescreen --procs 4`): chunks 1–21 tracked in `results/w3_moa/prescreen/`
  (46,883 tar members, 0 errors), ~80–90 s per chunk CPU-bound at 99 % (chunk 21: 216 s, slower archive reads).
  At this rate the pre-screen needs ~3 more routine runs; then fit (chunked by hand if > ~60 passes), vet, inject
  in chunks (`inject --chunk K/N`, from #165), limit, combine.
- Run 2 (cloud, 23:13–23:46Z, same command, `--chunks 22-119`): chunks 22–61 tracked (40 chunks, 0 errors,
  mean 120 s per chunk wall, 78 s when archive reads are fast); stopped by PID at the slot's end.
- Run 3 (cloud, 00:10–00:44Z, `--chunks 62-119`): chunks 62–93 tracked (32 chunks, 0 errors, mean 101 s per
  chunk wall); stopped by PID at the slot's end.
- Run 4 (cloud, 01:09–01:34Z, `--chunks 94-119`): chunks 94–119 tracked (26 chunks, 0 errors, ~60–90 s per chunk);
  the pre-screen is complete: **259,556 light curves, 1,802 shape passes, 126 not at a shared epoch** (to fit; gb4
  had 118). A fresh clone needs `pip install -e '.[mulens]'` for `fit`.
- Fit chunk 1/4 (cloud, 01:35–01:47Z): 32 fits, 0 errors, 31 flags, 699 s (`results/w3_moa/fits_gb5_chunk1of4.ecsv.gz`).
- Fit chunk 2/4 (cloud, 02:11–02:30Z): 32 fits, 0 errors, 27 flags, 1,131 s (`results/w3_moa/fits_gb5_chunk2of4.ecsv.gz`).
- Fit chunk 3/4 (cloud, 02:30–02:48Z): 31 fits, 0 errors, 27 flags, 1,060 s (`results/w3_moa/fits_gb5_chunk3of4.ecsv.gz`);
  chunk 4/4 was started and stopped unfinished at the slot's end.
- Fit chunk 4/4 (cloud, 03:10–03:27Z): 31 fits, 0 errors, 27 flags, 1,014 s. `merge-chunks --n 4` (126 fits),
  `vet`: **0 of 112 flags survive** (`results/w3_moa/vetting_gb5.json`); CDS XMatch worked (Gaia DR3 RP for 67 / 126).
  Contact sheet inspected: variables, box and eclipse dips, slow trends; caustic spikes only in the models.
- Merged `claude/w3-moa-gb4` (#165, itself stacked on #160) into this branch for `inject --chunk K/N`, so `combine`
  will cover seventeen fields. gb5 injections use 6 chunks (~370 each, fits one routine slot).
- Injection chunk 1/6 (cloud, 03:32–03:43Z): 367 in 612 s (0.60 per s), 0 errors; W3 333 → 25 flagged → 8 recovered;
  PSPL controls 0 / 34 flagged.
- Injection chunk 2/6 (cloud, 03:44–04:00Z): 367 in 985 s, 0 errors; W3 333 → 21 flagged → 11 recovered; PSPL 0 / 34.
- Injection chunk 3/6 (cloud, 04:16–04:34Z): 367 in 955 s, 0 errors; W3 334 → 26 flagged → 13 recovered; PSPL 0 / 33.
- Injection chunk 4/6 (cloud, 04:34–04:49Z): 367 in 881 s, 0 errors; W3 334 → 23 flagged → 12 recovered; PSPL 0 / 33.
- Injection chunks 5/6 and 6/6 (cloud, 05:10–05:31Z): 366 each in 485 s and 659 s, 0 errors; W3 333 → 26 → 9 and
  333 → 32 → 14 recovered; PSPL 0 / 33 each. Total: W3 2,000 → 153 flagged → 67 recovered (44 %); PSPL 0 / 200.
- `vet` re-run on the fresh clone reproduced the tracked record (0 / 112 survive, identical funnel).
- gb5 limits (`results/w3_moa/limits_gb5.ecsv`): Γ₉₅ ≈ 1.8 × 10⁻⁷–1.2 × 10⁻⁶ per star per year
  (N_s = 2.91 × 10⁷, 65 / 80 subfields); the weakest cells (t_E = 3 d, ρ = 0.01 and t_E = 300 d, ρ = 0.1) rest on
  2 recoveries each.
- **Combined over seventeen fields:** Γ₉₅ ≈ 1.6 × 10⁻⁸–9.0 × 10⁻⁸ per star per year (sixteen fields:
  1.7 × 10⁻⁸–9.8 × 10⁻⁸); max field share 0.17 (`results/w3_moa/limits_combined.ecsv`, docs/exotic_limits.md).
- **Next:** gb13 after the owner decides #145; then gb6, gb9, gb10, gb14 (chunked fits and injections). This PR is
  stacked on #165 and #160: merge those first.

## 2026-10-11: W3 MOA-II gb4 (sixteenth field): 0 survivors; combined limit over sixteen fields
- Hypothesis before running: 0 survivors under CHAIN_VERSION 2026-10-10.1; a survivor with caustic spikes in the data
  goes to `/vet-candidate` before injections.
- gb4 is the largest field so far: 66 pre-screen chunks of 4 GiB (~265 GB streamed, never stored). Run 1 (cloud,
  15:08–15:43Z, `prescreen --procs 4`): chunks 1–35 tracked in `results/w3_moa/prescreen/` (~112k light curves,
  ~57 per s, ~1 chunk per min, CPU-bound at 99 %).
- Run 2 (cloud, 16:10–16:50Z, `prescreen --procs 4 --chunks 36-66`): chunks 36–66 tracked (94,639 light curves in
  2,252 s, ~42 per s). Pre-screen complete: **205,886 light curves, 2,350 shape passes, 118 not at a shared epoch**
  (gb3: 75 of 167k), so the fit stage has ~1.6× gb3's work (~55 min at 4 procs).
- Run 3 (cloud, 17:12–17:45Z): `fit --procs 4 --chunk 1/2` tracked (`results/w3_moa/fits_gb4_chunk1of2.ecsv.gz`):
  59 fits in 1,727 s, 0 errors, 57 flags. `run-field` would fit all 118 passes as one chunk (`FITS_PER_CHUNK` = 120,
  ~55 min, past a routine slot), so the fit is split in two by hand.
- Run 4 (cloud, 18:08–18:31Z): `fit --procs 4 --chunk 2/2` tracked (`results/w3_moa/fits_gb4_chunk2of2.ecsv.gz`):
  59 fits in 932 s, 0 errors, 57 flags; `merge-chunks --n 2` (118 fits), `vet`: **0 of 114 flags survive**
  (`results/w3_moa/vetting_gb4.json`). CDS XMatch worked (Gaia DR3 RP for 67 / 118). Contact sheet inspected: box and
  eclipse dips, variables and slow trends; caustic spikes only in the models. The first fit attempt stopped on a
  truncated `metadata.ipac.tar.gz` download (sha256 mismatch, same content-length upstream); the retry was clean.
- Run 5 (cloud, 19:07–19:45Z): `inject` is now chunked and resumable (`inject --chunk K/N`, `merge-inject --n N`;
  tracked `results/w3_moa/injections_<F>_chunkKofN.ecsv.gz`, each carrying a SHA-256 of the seeded job list, so a
  re-merged pre-screen cannot mix runs; `run-field` resumes the chunks; a chunk reads only its own carriers, 240 of
  287). The seeded job list and the outcome per injection are unchanged; only the row order (now by job index) and the
  new `inj` column differ. Chunk 1/4 tracked: 550 injections in 1,352 s (0.41 per s), 0 errors; W3 500 → 47 flagged
  → 28 recovered; PSPL controls 0 / 50 flagged.
- Failed approach: `pkill -f` on the injection also killed the harness shell that ran it (scripts/CLAUDE.md warns);
  the first chunk-1 start was lost, ~3 min.
- Run 6 (cloud, 20:08–20:58Z): chunk 2/4 tracked: 550 injections on 246 carriers in 878 s (0.63 per s), 0 errors;
  W3 500 → 51 flagged → 23 recovered; PSPL controls 0 / 50 flagged. A fresh clone needs `merge-prescreen` before
  `inject` (it reads `data/derived/w3_moa/prescreen_gb4.ecsv`), and the `metadata.ipac.tar.gz` download was truncated
  once more (sha256 mismatch; the retry was clean).
- Failed approach: a wait loop `while pgrep -f "<inject command>"` never ends, because the pattern also matches the
  waiting shell's own command line; it held the run ~25 min past the chunk's end. Wait on the PID (`wait` or
  `kill -0 <pid>`) instead.
- Run 7 (cloud, 21:08–22:25Z): chunks 3/4 (1,746 s) and 4/4 (1,828 s), 0 errors; W3 500 → 52 / 57 flagged → 28 / 28
  recovered; PSPL controls 0 / 50 flagged each. `merge-inject --n 4` (2,200 injections on 287 carriers). `limit` needs
  the full `vetting_gb4.json`, which a fresh clone lacks (only the compact record is tracked), so `merge-chunks --n 2`
  and `vet` were re-run: **0 of 114 survive** again, and the regenerated compact record is byte-identical.
- gb4 limits (`results/w3_moa/limits_gb4.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 1.7 × 10⁻⁷–9.0 × 10⁻⁷ per star per year
  (N_s = 2.60 × 10⁷); vetting keeps 107 / 207 flagged W3 injections; the t_E = 300 d cells rest on 3 / 4 recoveries.
  **Combined, sixteen fields: Γ₉₅ ≈ 1.7 × 10⁻⁸–9.8 × 10⁻⁸** (was 1.9 × 10⁻⁸–1.1 × 10⁻⁷ over fifteen; max field share
  0.19). This branch merges `claude/w3-moa-gb3` (#160) so the combine includes gb3: merge #160 first.
- Failed approach: the venv lacked the `mulens` extra, so the first `inject` exited at once; install
  `.[dev,cloud,mulens]` on a fresh clone.
- `/code-review` (2026-10-12) fixes in the chunk resume: a tracked chunk is reused or merged only if it has no
  failed injection and matches the current job list (SHA-256) and baseline χ² calibration; `merge-inject` now checks
  the current pre-screen too, and `run-field` resumes with the chunks' own seed, sampling and carrier count instead
  of the defaults. The tracked gb4 chunks pass the new checks (`merge-inject --n 4`: 2,200 injections).
- **Next:** gb13 after the owner decides #145; then gb5, gb6, gb9, gb10, gb14 with chunked fits and injections
  (a fresh clone needs `merge-prescreen` and, before `limit`, `merge-chunks` + `vet`).

## 2026-10-11: W3 MOA-II gb3 (fifteenth field): 0 survivors; combined limit over fifteen fields
- Hypothesis before running: 0 survivors under CHAIN_VERSION 2026-10-10.1. A survivor with caustic spikes in the data
  would have gone to `/vet-candidate` before injections.
- Four cloud runs. Pre-screen (`prescreen --procs 4`, chunks 1–30, then 31–50): 167,431 light curves, 1,797 shape
  passes, 75 not at a shared epoch. `run-field --procs 4` (takeover 12:11Z) reused the pre-screen chunks, refit the 75
  passes in 2,070 s (0 errors, 72 flags), vetted (**0 survive**, `results/w3_moa/vetting_gb3.json`, same as the
  earlier vet) and ran 2,200 injections in 6,430 s (~0.34 per s; chain 8,885 s). CDS XMatch worked (Gaia DR3 RP for
  52 / 72). Contact sheet inspected: box and eclipse dips, variables, slow trends; caustic spikes only in the models.
- Limits (`results/w3_moa/limits_gb3.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 1.5 × 10⁻⁷–1.5 × 10⁻⁶ per star per year
  (N_s = 2.27 × 10⁷); the t_E = 300 d cells rest on 2 recoveries each. **Combined, fifteen fields: Γ₉₅ ≈ 1.9 × 10⁻⁸–
  1.1 × 10⁻⁷** (was 2.1 × 10⁻⁸–1.2 × 10⁻⁷; max field share 0.19).
- Failed approach: the earlier run's tracked fit chunk (98a6aac) carried the default `w3_microlensing` bounds
  (t_E ≤ 5,000 d, log ρ ≤ 0) in its header, not `MOA_FIT_BOUNDS`, so `run-field` did not reuse it and refit (1a4689c).
  Vetting is unchanged. Run fits through the `w3_moa.py` CLI, which sets the MOA bounds in `main`.
- Timing: injections ran at ~0.34 per s on 4 procs (gb8: 0.77 per s), so a 167k-curve field needs ~110 min for
  injections alone. `inject` still writes only at the end (gb8 entry).
- **Next:** gb13 after the owner decides #145; then gb4–gb6, gb9, gb10, gb14.

## 2026-10-11: W3 MOA-II gb8 (fourteenth field): 0 survivors; combined limit over fourteen fields
- Hypothesis before running: 0 survivors under CHAIN_VERSION 2026-10-10.1. A survivor with caustic spikes in the data
  would have gone to `/vet-candidate` before injections.
- `run-field --procs 4` (cloud), split across three runs: 120,726 light curves streamed in 1,068 s (19 tracked
  pre-screen chunks; 2,256 shape passes, 77 off shared epochs); 77 fits in 675 s, 0 errors, 73 flags; **0 survive**
  (`results/w3_moa/vetting_gb8.json`, reproduced by the third run's re-vet). CDS XMatch worked (Gaia DR3 RP for
  43 / 73). Contact sheet inspected: variables, box and eclipse dips, slow trends; caustic spikes only in the models.
- Injections (2,200 on 267 carriers) in 2,876 s (third run; the second run's were lost when it ended mid-run). Limits
  (`results/w3_moa/limits_gb8.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 1.7 × 10⁻⁷–1.0 × 10⁻⁶ per star per year
  (N_s = 2.23 × 10⁷); the t_E = 300 d cells rest on 3 recoveries each. **Combined, fourteen fields: Γ₉₅ ≈ 2.1 × 10⁻⁸–
  1.2 × 10⁻⁷** (was 2.4 × 10⁻⁸–1.4 × 10⁻⁷; max field share 0.19).
- Failed approach: a takeover run that starts injections near the end of its slot loses them, because `inject`
  writes only at the end. Start injections at the top of a run, or make `inject` chunked and resumable like `fit`.
- **Next:** gb13 after the owner decides #145; then gb3–gb6, gb9, gb10, gb14.

## 2026-10-11: W3 MOA-II gb1 (thirteenth field): 0 survivors; combined limit over thirteen fields
- Hypothesis before running: 0 survivors under CHAIN_VERSION 2026-10-10.1 and a finite limit in every cell. A
  survivor with caustic spikes in the data would have stopped injections and gone to `/vet-candidate`.
- `run-field --procs 4` (cloud): 70.1 GB streamed in 1,249 s (109,715 light curves; 1,619 shape passes, 51 off shared
  epochs); 51 fits in 590 s, 48 flags; **0 survive**; 2,200 injections in 2,982 s; chain 4,932 s. Contact sheet
  inspected: variables, box dips, one long trend; caustic spikes only in the models. CDS XMatch worked (RP for 30 / 48).
- Limits (`results/w3_moa/limits_gb1.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 1.7 × 10⁻⁷–1.7 × 10⁻⁶ per star per year
  (N_s = 2.10 × 10⁷); the t_E = 300 d cells rest on 2 recoveries each. **Combined, thirteen fields: Γ₉₅ ≈ 2.4 × 10⁻⁸–
  1.4 × 10⁻⁷** (was 2.8 × 10⁻⁸–1.5 × 10⁻⁷; max field share 0.21, gb1 at t_E = 3 d, ρ = 0.01).
- Timing: injections ran at ~0.5–1 per s (gb2: ~1.3 per s), so a 70-GB field needs ~85 min, past one routine slot.
- **Next:** gb8 (78 GB; start it at the top of a run or split stages across runs: pre-screen chunks and fits are
  tracked and reused), `combine`; gb13 waits on the owner (#145).

## 2026-10-11: W3 short events (t_E = 3 d): `eclipse_dip` losses are unsampled spikes; the `jackknife_nights` share is undecided (D-082)
- Hypothesis: sparse sampling leaves the caustic spikes on few nights, so genuine peak-dip-peak signals fail as
  ordinary dips. If failed injections had well-sampled spikes, the rule would be the loss.
- Frozen sample: 2,000 t_E = 3 d W3 injections, production chain, gb12 carriers, seed 3003 (846 s on 4 cores).
  dev / validation are assigned by a SHA-256 of each injection's inputs. Flagged → recovered: dev 57 → 21
  (0.37, 95 % 0.24–0.51); validation 61 → 27 (0.44), read for the baseline only.
- Dev: 11 / 12 `eclipse_dip` losses have ≤ 1 spike night, so they are indistinguishable from a flat dip in the
  data. The `jackknife_nights` losses (14) and the recoveries (21) overlap in spike nights, so the proxy cannot
  attribute that loss. Stage 2 must record the jackknife's own feature nights.
- **Failed approach (rule):** a split keyed on the row index is not frozen. `run_inject` writes rows in pool
  completion order, so slow (flagged) injections collect late (`/code-review`; flagged 72/46, p ≈ 0.02). Key a
  split on each injection's inputs. The first draft's "≤ 8 % rule-attributable" bound is withdrawn.
- Injection profile (`docs/audits/2026-10-10-acceleration-v2.md`, V2.1 priority 2): 97 % of injection CPU is
  fitting and vetting of the pre-screen passes (median 18.6 s each); the finite-source integrand is 46 % of that.
  Data reuse, vectorized injection and caching cannot gain more than ~3 %.
- **Next:** stage 2 of the short-event study (re-vet the dev flags with the jackknife details recorded; then a
  night-robust likelihood term; FP check on real flags and PSPL controls; one validation run); a compiled W3
  finite-source integrand with an identical-outcome test on fixed injections.

## 2026-10-10: W3 MOA-II gb2 (twelfth field): 0 survivors; combined limit over twelve fields
- Hypothesis before running: 0 survivors under CHAIN_VERSION 2026-10-10.1 and a finite limit in every cell. A
  survivor that held up on the contact sheet would have stopped injections and gone to `/vet-candidate`.
- `run-field --procs 4` (cloud): 60.5 GB streamed in 774 s (95,624 light curves; 2,189 shape passes, 61 off shared
  epochs); 61 fits in 291 s, 60 flags; **0 survive**; 2,200 injections; chain 1,754 s. Contact sheet inspected: box
  dips, variables, sparse few-night dips, caustic spikes only in the models. CDS XMatch worked (RP for 43 / 60).
- Limits (`results/w3_moa/limits_gb2.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.1 × 10⁻⁷–9.9 × 10⁻⁷ per star per
  year (N_s = 1.76 × 10⁷); every cell has ≥ 4 recoveries. **Combined, twelve fields: Γ₉₅ ≈ 2.8 × 10⁻⁸–1.5 × 10⁻⁷**
  (was 3.2 × 10⁻⁸–1.8 × 10⁻⁷; max field share ≤ 0.20).
- Failed approach (cloud I/O): `metadata.ipac.tar.gz` (97 MB) came back truncated twice through `urllib` (two
  different wrong sha256 values); one `curl -C - --retry 3` into `data/raw/moa/b339ec176933_metadata.ipac.tar.gz`
  gave the pinned file. Pre-fetch it that way before `run-field` in a fresh cloud session.
- **Next:** gb1 (70 GB) and gb8 (78 GB), `combine` after each; gb13 waits on the owner (#145).

## 2026-10-10: Neutrino Frontier benchmark 1 reproduced: SkyLLH TXS 0506+056 2014–15 box on IceTracks-DR2
- Hypothesis before running: SkyLLH 26.1.0 on IC86_IV with the Table 6 box gives n̂s ≈ 12.7, γ̂ ≈ 2.3
  (arXiv:2605.19040 section 5). A miss by > 1 event or > 0.1 in γ (ASSUMPTION) would mean our DR2 handling is wrong.
- `scripts/nf_txs_benchmark.py` (D-080): the printed window (MJD 56927.5–57112.5) gives n̂s = 11.1, γ̂ = 2.22, TS = 20.3
  under both IRF versions, 1.6 events short. Table 6 rounds T0 and ΔT to whole days, and an on-source event
  (MJD 57112.653, 0.38° from TXS, σ = 0.20°) lies 0.15 d past the printed stop. With both edges snapped to the
  on-source events (56927.860–57112.653: T0 = 57020.26, ΔT = 184.8 d) the fit gives **n̂s = 12.72, γ̂ = 2.26, TS = 23.0**
  (Dataverse 1.0 IRFs, as in the paper), or 12.75 / 2.26 with the ≥ 2.0 IRFs. Paper: 12.7 / 2.3 (SkyLLH), 12.56 / 2.26
  (internal tools). **Reproduced**, with a post hoc window rule (D-080): it shows the paper's box ends on these
  events; it is not a blind test.
- Fixed-window background (v2 IRFs, SkyLLH scrambles, 50,000 trials per window): 0 trials reach the observed TS, p < 2 × 10⁻⁵
  (half-χ²₂ Wilks: 2 × 10⁻⁵ / 5 × 10⁻⁶). This is a p-value for a fixed window, not the paper's free-window pre-trial
  p (4.18 × 10⁻³). It is not a new result.
- Failed approaches: multi-season time-dependent fits (SkyLLH supports one dataset only; the box is inside IC86_IV
  anyway). `pkill -f <script>` killed the waiting shell itself, as scripts/CLAUDE.md warns.
- Disk: about 1.4 GB of IRFs streamed, sha256-checked (`data/manifests/nf_skyllh_dr2.ecsv`) and deleted (`--cleanup`).
- **Next:** time-integrated benchmark (Table 8: TXS 8.8 / 2.0, NGC 1068 80.1 / 3.2; all 14 seasons, 3.4 GB of IRFs,
  cloud only). It is the first SkyLLH check without a free window choice.

## 2026-10-10: W3 MOA-II cloud re-vet of the lenient-reference fields: 0 survivors
- Hypothesis before running: a working Gaia DR3 RP reference only tightens `exotic_in_domain`, so 0 survivors stays.
  A new survivor would have reopened the combined limit.
- `merge-prescreen`, `merge-chunks --n 1`, `vet --procs 4` and `summary` for gb7, gb11, gb16 and gb19–gb22 (41–127 s
  per field; no injections needed: they never used the reference). CDS XMatch worked for every field (RP for
  192 / 276 flags). **0 survivors.** 274 flags fail the same first test as before. gb19-R-1-0-51315 and
  gb20-R-5-5-38429 fail a different ordinary test (another in-domain exotic solution); contact sheets inspected,
  no caustic spikes in the data. The combined limit (below) stands.
- **Next:** gb13 (stopped at 25 / 78 fits), then more Nunota et al. 2024 Table 1 fields, `combine` after each.

## 2026-10-10: W3 MOA-II combined limit over the eleven re-run fields
- Hypothesis before running: `combine` gives a finite limit in every cell, and no single field's few recoveries
  carry a cell (the gb18 t_E = 300 d, ρ = 0.1 cell has 0 / 200). A cell with share near 1 would have meant more
  injections before quoting it.
- `w3_moa.py combine` (`results/w3_moa/limits_combined.ecsv`, docs/exotic_limits.md "Combined limit"): 11 fields,
  0 survivors, Σ N_s ≈ 1.37 × 10⁸ stars, 1.18 × 10⁹ star-years. **Γ₉₅ ≈ 3.2 × 10⁻⁸–1.8 × 10⁻⁷ W3 events per star
  per year** (t_E 3–300 d, ρ 0.01/0.1); best 3.2 × 10⁻⁸ at t_E = 30 d, ρ = 0.01. Injection-counting uncertainty
  8–18 %. New column `max_field_share` (≤ 0.19 everywhere); a zero-efficiency field cell adds no exposure (test).
- Failed approach (coordination): a 60-min-stale claim (#142, no heartbeat during a ~32-min run) was taken over
  while its session was still alive; it merged 5 min later and the duplicate run was stopped. Heartbeat during
  long runs, and re-check the PR before starting a takeover run.
- **Next:** the cloud re-vet of the lenient-reference fields (gb7, gb11, gb16, gb19–gb22); then gb13 and more
  Nunota et al. 2024 fields (re-run `combine` after each).


## 2026-10-10: W3 MOA-II gb18 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors; gb18-R-9-4-24509 (the D-068 chain survivor, a red giant's slow dimming)
  should now fail inside the chain. A survivor that held up on the contact sheet would have changed the plan.
- `run-field --procs 4` in 1,890 s: all 10 pre-screen chunks reused (83,855 light curves, 1,497 shape passes, 70 off
  shared epochs); 70 fits 400 s, 67 flags; 2,200 injections 1,126 s. **0 survive.** gb18-R-9-4-24509 fails
  `slow_dip_seasons` (ΔBIC 8.6). Contact sheet inspected: box dips, slow dimmings, variables; caustic spikes only
  in the models. CDS XMatch worked from the cloud (Gaia DR3 RP for 39 / 67 flags).
- Limits (`results/w3_moa/limits_gb18.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.8 × 10⁻⁷–2.3 × 10⁻⁶ per star per year
  (N_s = 1.55 × 10⁷, 78 / 80 subfields); the t_E = 300 d, ρ = 0.1 cell recovered 0 / 200 (no limit there).
  Vetting keeps 83 / 162 flagged W3 injections (51 %); PSPL controls: 0 false W3 calls.
- All eleven D-068 fields are now re-run under the new chain: 0 survivors in each.
- Failed approach (tooling): a `pgrep -f "<pattern>"` wait loop matches its own shell command line and never exits;
  wait on the PID instead.
- **Next:** `combine` (check it handles an infinite-Γ cell), the cloud re-vet of the lenient-reference fields
  (gb7, gb11, gb16, gb19–gb22), then gb13.

## 2026-10-10: W3 MOA-II gb17 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors; gb17-R-6-1-3829 (the D-068 chain survivor, a ~30-d periodic variable)
  should now fail inside the chain. A survivor that held up on the contact sheet would have changed the plan.
- `run-field --procs 4` in 1,797 s: all 13 pre-screen chunks reused (D-078 addendum 3, 0 s instead of ~8 min;
  100,448 light curves, 1,437 shape passes, 85 off shared epochs); 85 fits 479 s (fit chunks still key on the whole
  `Params`), 75 flags; 2,200 injections 1,175 s. **0 survive.** gb17-R-6-1-3829 fails `exotic_chi2_cap`
  (6.67 > 5.46), as the D-068 addendum intended. Contact sheet inspected: box dips, slow dimmings, variables;
  caustic spikes only in the models.
- CDS XMatch worked from the cloud again (Gaia DR3 RP for 53 / 75 flags).
- Limits (`results/w3_moa/limits_gb17.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.3 × 10⁻⁷–1.9 × 10⁻⁶ per star per year
  (N_s = 1.80 × 10⁷, Nunota et al. 2024, 79 / 80 subfields; the 300-d cells rest on 2 recoveries). Vetting keeps
  86 / 156 flagged W3 injections (55 %); PSPL controls: 0 false W3 calls.
- **Next:** gb18 (the last re-run field), then `combine`; the fit-chunk reuse fix (D-078 audit) would save ~8 min here.

## 2026-10-10: W3 MOA-II gb15 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068; a survivor that held up on the contact sheet would have
  changed the plan (vetting record, owner notified).
- `run-field --procs 4` in 2,686 s (pre-screen re-streamed ~810 s, 8 chunks, 56–200 light curves/s as archive
  throughput fell from 78 to 22 MB/s; 83,145 light curves, 4,936 shape passes, 122 off shared epochs; fits 2 chunks,
  118 flags; 2,200 injections 1,182 s). **0 survive.** First failing test: repeated deficit 30, residual deficit 24,
  eclipse dip 18, bracketing 14, fit domain 12, χ² cap 10, other 11. Contact sheet inspected: smooth U/V dimmings,
  eclipse-like boxes, variables; caustic spikes only in the models.
- CDS XMatch worked from the cloud again (Gaia DR3 RP reference for 88 / 118 flags).
- Limits (`results/w3_moa/limits_gb15.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 3.2 × 10⁻⁷–3.4 × 10⁻⁶ per star per year
  (N_s = 1.04 × 10⁷, Nunota et al. 2024, only 62 / 80 subfields: the low end of 1.04–1.34 × 10⁷; both t_E = 300 d
  cells rest on 2 recoveries). Vetting keeps 114 / 171 flagged W3 injections (67 %); PSPL controls: 0 false W3 calls.
- The re-streamed gb15 pre-screen rows equal the tracked ones row for row (only header metadata differs), a third
  check of D-078 addendum 3; the tracked chunks are kept unchanged.
- **Next:** gb17, gb18 (pre-screen now reused, D-078 addendum 3), then `combine`.

## 2026-10-10: Acceleration V2 audit, first pass: MOA pre-screen chunks reused across vetting changes; injection-loss map
- Owner brief "Autonomous research acceleration V2": measured audit in `docs/audits/2026-10-10-acceleration-v2.md`
  (what the repository already covers, the gb12 injection-loss funnel, ranked next units).
- **Redundancy removed:** a vetting-only `Params` change invalidated every streamed MOA pre-screen chunk; gb12
  (342 s) and gb7 (487 s) re-streamed rows byte-identical to the tracked ones. Chunk validity now depends only on
  the pre-screen parameters, the tracked-row constants and `PRESCREEN_CODE` (D-078 addendum 3). gb15, gb17, gb18:
  all 31 tracked chunks are reused (0 s instead of ~6–8 min each).
- **Where W3 injections are lost (gb12):** Cut-0 45 % kept, shape cut 23 %, shared-epoch test 75 %, vetting 63 %.
  Faint sources (I_s > 19.5) carry most of the loss (near the noise). Testable design losses: the shared-epoch veto
  (25 % of shape passes) and short / long t_E vetting (`jackknife_nights` at 3 d, `feature_bracketed` at 300 d).
  Sparse sampling does make real peak-dip-peak signals fail as ordinary dips at t_E = 3 d.
- **Next:** gb15, gb17, gb18 with the reuse fix, then `combine`; a shared-epoch masking study on a dev/validation
  injection split with the real flags as the false-positive check.

## 2026-10-10: W3 MOA-II gb12 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068; a survivor that held up on the contact sheet would have
  changed the plan (vetting record, owner notified).
- `run-field --procs 4` in 1,528 s (pre-screen re-streamed 342 s at ~215 light curves/s, CPU 99 %; 73,386 light
  curves, 5,532 shape passes, 51 left after the shared-epoch test; 51 fits 203 s, 46 flags; 2,200 injections). **0 survive.**
  First failing test: repeated deficit 17, residual deficit 9, eclipse dip 9, bracketing 4, χ² cap 3, other 4.
  Contact sheet inspected: box dips, quasi-periodic variables, caustic spikes only in the models.
- The batched CDS XMatch (cat2 `vizier:I/355/gaiadr3`, 1″) worked from the cloud this time: Gaia DR3 RP reference
  for 32 / 46 flags, so the source-flux bound used real references, unlike the earlier cloud runs.
- Limits (`results/w3_moa/limits_gb12.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.3 × 10⁻⁷–1.5 × 10⁻⁶ per star per year
  (N_s = 1.61 × 10⁷, Nunota et al. 2024, 79 / 80 subfields). Vetting keeps 99 / 156 flagged W3 injections (63 %);
  PSPL controls: 0 false W3 calls.
- **Next:** gb15, gb17, gb18, then `combine`; audit which vetting tests remove injected W3 signals
  (37 % loss per field).

## 2026-10-10: W3 MOA-II gb7 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068 (gb7-R-8-6-94052 was the known jackknife case); a survivor
  that held up on the contact sheet would have changed the plan (vetting record, CDS-dependent tests where CDS works).
- `run-field --procs 4` in 2,670 s (pre-screen re-streamed 487 s at ~150 light curves/s, CPU 99 %; 75,328 light
  curves, 1,450 shape passes, 44 off shared epochs; 44 fits 289 s, 42 flags; 2,200 injections 1,297 s).
  **0 survive**; gb7-R-8-6-94052 still fails `jackknife_nights`. Contact sheet inspected.
- Limits (`results/w3_moa/limits_gb7.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.7 × 10⁻⁷–1.4 × 10⁻⁶ per star per
  year (N_s = 1.63 × 10⁷, Nunota et al. 2024, 78 / 80 subfields; the t_E = 300 d, ρ = 0.1 cell rests on 3
  recoveries). Vetting keeps 100 / 169 flagged W3 injections (59 %); pre-screen pass 6–12.5 %. PSPL controls: 0
  false W3 calls. Quiet χ²/dof 95th percentile 6.1.
- CDS XMatch failed again from the cloud ("Too many jobs", recorded in `vetting_gb7.json`): lenient default
  source-flux bound; 0 survivors stands.
- **Next:** gb12, gb15, gb17, gb18, then `combine`.

## 2026-10-10: W3 MOA-II gb11 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068; a survivor that held up on the contact sheet would have
  changed the plan (vetting record, CDS-dependent tests where CDS works).
- `run-field --procs 4` in 1,981 s (pre-screen re-streamed 334 s at ~190 light curves/s, CPU 99 %; 64,090 light
  curves, 1,579 shape passes, 44 off shared epochs; 44 fits 238 s, 41 flags; 2,200 injections 973 s). **0 survive**.
  First failing test: eclipse dip 13, repeated deficit 12, residual deficit 6, bracketing 4, robust errors 2,
  χ² cap 1, neighbour 1, smooth dip 1, slow dip seasons 1. Contact sheet inspected (eclipse-like few-day boxes, slow
  dimmings, a deep V-shaped ~15-d dip in gb11-R-2-7-83307 that the box model misfits; caustic spikes only in models).
- Limits (`results/w3_moa/limits_gb11.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 4.0 × 10⁻⁷–1.1 × 10⁻⁶ per star per
  year (N_s = 1.09 × 10⁷ from Nunota et al. 2024, 76 / 80 subfields; every cell ≥ 6 recoveries). Vetting keeps
  103 / 163 flagged W3 injections (63 %); pre-screen pass 5.5–11 %. PSPL controls: 0 false W3 calls. Quiet χ²/dof
  95th percentile 6.1.
- CDS XMatch failed again from the cloud ("Too many jobs"; recorded in `vetting_gb11.json`): lenient default
  source-flux bound; 0 survivors stands.
- **Next:** `run-field` for gb7, gb12, gb15, gb17, gb18, then `combine`.

## 2026-10-10: W3 MOA-II gb16 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068; a survivor that held up on the contact sheet would have
  changed the plan (vetting record, CDS-dependent tests where CDS works).
- `run-field --procs 4` in 2,148 s (pre-screen re-streamed 381 s at ~170 light curves/s, CPU 99 %; 64,353 light
  curves, 1,427 shape passes, 32 off shared epochs; 32 fits, 31 flags; 2,200 injections 929 s). **0 survive**.
  First failing test: repeated deficit 10, neighbour 6, residual deficit 4, exotic domain 3, eclipse dip 3,
  bracketing 2, feature sampled 1, χ² cap 1, jackknife 1. Contact sheet inspected (slow dimmings, short few-point
  dips, a deep ~20-d dip in gb16-R-9-6-97939 that the box model misfits; caustic spikes only in the models).
- Limits (`results/w3_moa/limits_gb16.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 2.9 × 10⁻⁷–2.3 × 10⁻⁶ per star per
  year (N_s = 1.50 × 10⁷ from Nunota et al. 2024, 79 / 80 subfields; both t_E = 300 d cells rest on 2 recoveries).
  Vetting keeps 71 / 122 flagged W3 injections (58 %); pre-screen pass 3–11 %. PSPL controls: 0 false W3 calls.
  Quiet χ²/dof 95th percentile 6.2.
- CDS XMatch failed again from the cloud (non-VOTable reply): lenient default source-flux bound; 0 survivors stands.
  `summary` now keeps this `gaia_rp_xmatch` status in the tracked `vetting_gbN.json` (TASKS follow-up).
- **Next:** `run-field` for gb11, gb7, gb12, gb15, gb17, gb18, then `combine`.

## 2026-10-10: W3 MOA-II gb19 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- Hypothesis before running: 0 survivors as in D-068; a survivor that held up on the contact sheet would have
  changed the plan (vetting record, CDS-dependent tests where CDS works).
- `run-field --procs 4` in 1,596 s (pre-screen re-streamed, rows identical to the tracked tables; 51 fits; 2,200
  injections 686 s). 48 flags, **0 survive**; gb19-R-4-4-31159 now fails `feature_bracketed`. Contact sheet
  inspected (variables, slow dimmings, season steps, eclipse-like boxes; no caustic spikes in the data).
- Limits (`results/w3_moa/limits_gb19.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 4.1 × 10⁻⁷–2.9 × 10⁻⁶ per star per
  year (t_E 3–300 d, N_s = 1.21 × 10⁷ from Nunota et al. 2024; the 2.9 × 10⁻⁶ at t_E = 300 d, ρ = 0.01 rests on 2
  recoveries). Vetting keeps 78 / 112 flagged W3 injections (70 %); pre-screen pass 2–8.5 % stays the limiting
  factor. PSPL controls: 0 false W3 calls. Quiet χ²/dof 95th percentile 6.4.
- CDS XMatch failed again from the cloud ("Too many jobs"): lenient default source-flux bound; 0 survivors stands.
- **Next:** `run-field` for gb16, gb11, gb7, gb12, gb15, gb17, gb18, then `combine`.

## 2026-10-10: E1 re-pinned to one GBM snapshot (8cc2…, 4,392 bursts): every E1 conclusion unchanged (D-078 addendum 2)
- Hypothesis before running: only bursts added or revised after the old pins differ, so every E1 channel stays
  null; a family p below 0.01 would have changed the plan.
- Fresh fetch: only `fermigbrst.vot` changed (sha256 8cc2823f…, 424,672 B, 4,392 rows; D-078 lists 4,391 rows for
  this digest). ICECAT-1, GWTC and CHIME Cat 2 match their pins. `data/manifests/e1_gw_events.ecsv` is removed:
  every E1 script (GW ones included) now reads the single `e1_events.ecsv`, so the pins cannot diverge again.
- Like-for-like check (observed counts only): the new snapshot minus bn261008763 and bn261007236 (triggered
  2026-10-07/08, last modified 2026-10-09) reproduces all 75 tracked D-074 observed counts. The other revision
  after 2026-10-09 (bn260930833) changes no cell. No trigger is later than 2026-10-08, so a trigger-date cut is a
  no-op. With the new bursts, 5 GBM–GBM cells gain 1–10 pairs (1d–7d all 17,523 → 17,533).
- Re-runs (same seeds and scramble counts as tracked; 4 cores): `e1_events.py` 456 s, antipodal (`--n 2000`),
  signed lag, CHIME flag, GW directional, GW signed (15–51 s each).

  | Result | old pin | new pin |
  |---|---|---|
  | D-074 five requested channels, family p (jit) | 0.278 | 0.273 |
  | D-074 wide + GW cells, pooled p (jit / jitday) | 0.0058 / 0.0033 | 0.0053 / 0.0029 |
  | D-074 all 75 cells, global p (perm / jit / jitday) | 0.0021 / 0.0036 / 0.0047 | 0.0020 / 0.0043 / 0.0041 |
  | Antipodal global p (addendum 2) | 0.557 | 0.531 |
  | Signed-lag global p (addendum 3) | 0.676 | 0.568 |
  | GW sky-map channels global p (addendum 4) | 0.829 | 0.833 |
  | Signed-lag GW cells global p (addendum 5) | 0.555 | 0.551 |
  | CHIME flag tests (CHIME only) | — | identical |

  The wide + GW driver is still CHIME–CHIME 1 h–1 d wide, which the calibrated rate-modulated null explains
  (D-074 addendum; CHIME-only, not re-run). Changes outside GBM cells are Monte Carlo noise from shifted
  random streams: one RNG per scramble serves all catalogues in order, so two extra GBM bursts reshuffle the
  CHIME and ICECAT draws too (|Δz| ≤ 0.36, |Δp| ≤ 0.01 in counts.ecsv). Injection outputs in limits.ecsv move more
  (40 trials per cell): `ul95_rate_per_anchor` ≤ 24 %, `eff` up to 19 % (CHIME–CHIME 1d–7d wide 0.999 → 0.806),
  `n50_detect` GBM–ICECAT 0s–10s wide 3 → 5.
- Fixed after `/code-review`: the null-ensemble cache name now carries the input digests (an old-snapshot cache
  was silently reused without `--refresh`); `e1_gw_events.ecsv` and its override removed (see above).
- Not done (follow-up): per-catalogue random streams in the nulls (`SeedSequence` keyed by catalogue) would make
  re-pins like-for-like outside the changed catalogue, but they change every tracked number, so they need their
  own PR.
- **Open risk / next (local session):** `main` pins a snapshot that only this ephemeral session held (cloud
  sessions cannot create releases, D-078). HEASARC serves 8cc2… until its next GBM update, so a local session
  should fetch it right away, check the digest, publish it to the derived-data store and wire `e1_events.fetch`
  to it. If the digest has moved on, run `e1_events.py --refresh` there and publish that snapshot instead.

## 2026-10-10: W3 MOA-II gb20 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, limits tracked
- `run-field --procs 4` in 2,144 s (pre-screen re-streamed 299 s, rows identical to the tracked tables; 51 fits;
  2,200 injections 1,075 s). 48 flags (as in D-068), **0 survive**; contact sheet inspected (eclipse-like boxes,
  slow dimmings, season steps, variables; no caustic spikes in the data).
- Limits (`results/w3_moa/limits_gb20.ecsv`, docs/exotic_limits.md): Γ₉₅ ≈ 5.0 × 10⁻⁷–6.5 × 10⁻⁶ per star per
  year (t_E 3–300 d, published N_s; the 6.5 × 10⁻⁶ at t_E = 300 d, ρ = 0.1 rests on 1 recovered injection).
  Quiet χ²/dof 95th percentile 6.7 (gb22 5.3, gb21 6.3; D-068 revisits above ~10). Vetting keeps 73 / 124 flagged W3 injections (59 %); the four new tests cost
  8. Pre-screen pass 2.5–11 % stays the limiting factor. PSPL controls: 0 false W3 calls.
- Failed approach (env): CDS XMatch from the cloud failed (gb22/gb21 too), so the source-flux bound used the lenient
  default; 0 survivors stands (docs/exotic_limits.md, note under the gb20 table).
- gb19 was planned in the same run but not started (40-min run budget).
- **Next:** `run-field` for gb19, gb16, gb11, gb7, gb12, gb15, gb17, gb18, then `combine`.

## 2026-10-10: First derived-data release `derived-data-20261010` (GW sky maps); pinned GBM snapshot not recoverable (D-078 addendum)
- Local session on the owner's machine (cloud sessions cannot create releases); details in the D-078 addendum.
  Asset `gw_skymaps_nside32.npz`, 4,337,802 B, sha256 ab0330588b9a380f2f49b0d4dae809ea6cd1ed6ac1d99718e23f66d00c8e1ad5
  (pre-release, not latest), pinned in `data/manifests/derived_data.ecsv`.
- 200 MB rule: the 704 MB of tarballs were streamed into memory and never saved; only reduced products stay in the
  data root (the 4.3 MB npz and four per-tarball caches, 28 MB). The re-reduction reproduces every pin.
- Fresh session (empty data root): 2.7 s wall instead of streaming (row added to the Part 0 speed table below).
- GBM: neither pinned `fermigbrst.vot` is on this machine; nothing downloaded or invented.
- Fixed after `/code-review` (8 findings, 8 addressed): `derived_publish.py` creates the tag at the recorded code commit
  (`--target`; this release's tag sits on c92249c, its code commit is 67fffba); the TASKS re-pin item now names
  the procedure (one snapshot, a copy for the GW manifest, maps first, a like-for-like window).
- **Next (cloud routine):** E1 GBM re-pin (TASKS "Part 0"); TXS SkyLLH benchmark; round 2. Next local session:
  publish the new GBM snapshot and wire `e1_events.fetch` to the durable copy.

## 2026-10-10: W3 MOA-II gb22 and gb21 re-run under CHAIN_VERSION 2026-10-10.1: 0 survivors, first new-chain limits
- `w3_moa.py --field gbN run-field` (~20 min per field on 4 cores). gb22: 30 flags, gb21: 36 flags, **0 survive**;
  contact sheets inspected (dips, season steps, scatter; no caustic spikes in the data). Pre-screen rows identical to
  the tracked tables (only the `Params` header changed).
- Limits (`results/w3_moa/limits_gb22.ecsv`, `limits_gb21.ecsv`; docs/exotic_limits.md "Re-run under
  CHAIN_VERSION 2026-10-10.1"): gb22 Γ₉₅ ≈ 1.4–10 × 10⁻⁶, gb21 5.3 × 10⁻⁷–2.8 × 10⁻⁶ per star per year
  (t_E 3–300 d; gb21 has a published N_s). The efficiency is
  pre-screen-limited (3–12 % pass); vetting keeps 62 % (gb22) / 56 % (gb21) of flagged injections, and the four
  new tests cost ~5 %. PSPL controls: 0 false W3 calls.
- Failed approach (ops): `pkill -f q.sh` in a watcher killed the watcher itself (its own command line matched); kill
  queue processes by PID.
- **Next:** `run-field` for gb20, gb19, gb16, gb11, gb7, gb12, gb15, gb17, gb18 (gb17/gb18 have the most light
  curves, ~35 min each), then `combine`; the pre-screen dominates the efficiency loss, so a better pre-screen is the
  lever for deeper limits.

## 2026-10-10: E-NF1b IceTracks-DR2 short-lag ghost pairs: null; DR2 cannot reach R_g ~ 10⁻³ (D-079)
- Streamed IceTracks-DR2 v3.1 (28 files, ~210 MB, 3 s, sha256 pinned; raw not kept): 1,643,355 events, 4,963.4 d
  good-run union, no duplicates. Known case: TXS 0506+056 2014–15 box 6 on vs 1.11 expected (p = 1.0 × 10⁻³).
- Forecast first (most optimistic: every event astrophysical): reachable R ≥ 0.055 (1 h–1 d), 0.14 (1–7 d),
  0.27 (7–30 d), 0.64 (30–180 d) for every northern cut, so DR2 cannot improve ICECAT-1 beyond 1 h and R_g ~ 10⁻³ is
  out of reach at every lag (TASKS target withdrawn).
- Pre-registered short-lag test (6 cells, 5,000 uptime-aware scrambles; 33 s wall on 4 cores: null 12.4 s,
  injections 3.8 s): **null**, min Bonferroni
  p = 0.48. 95 % CLs limits per northern track: R ≤ 0.002 (log10 E ≥ 4, ≤ 100 s), ≤ 0.005 (≥ 4.5, ≤ 100 s),
  ≤ 0.01–0.02 (100 s–1 h); R_g per astrophysical neutrino is R / f_astro (not estimated).
- Found: the 9 pairs at 0–10 s (log10 E ≥ 4) are all from IC40/IC59, where the northern rate above the cut is
  ×7–18 the IC79/IC86 rate (the energy proxy differs by configuration); the jitter null absorbs it (p = 0.08).
  Post hoc IC86-only: 0 pairs at ≤ 100 s. Any future DR2 energy cut must be per season.
- `/code-review` (10 findings, 9 fixed): good-run union instead of a sum; cache tied to the manifest digest;
  `--mjd-min` also restricts the null's uptime and labels full-release fields; (run, event) group by unique
  index; retry with back-off; injections on the process pool (×4); shared `empirical_p`/`DETECT_P`/floor;
  `wide_pair_counts` reuses `wide_pair_stats`.
- Failed approach: urllib gets 403 from Dataverse (use requests/curl); Dataverse md5 is of the original CSV, so the
  manifest pins the sha256 of the served `.tab`.
- **Next:** NF-H04 needs a quantitative prediction before more compute (TASKS); TXS SkyLLH benchmark stays
  open for other DR2 uses.

## 2026-10-10: Owner brief: Neutrino Frontier round 1 (E-NF1 ghost pairs: null) and Part 0 speed-ups (D-077, D-078)
- **Neutrino Frontier** (new program, `docs/neutrino_frontier/`). N-A (literature-blind) wrote four frameworks;
  N-B (adversarial, checked references) judged three KNOWN-REDUNDANT: stochastic distance (NF-H01), pseudo-Dirac
  partners (NF-H02; N-A's δm² window is already excluded, arXiv:2406.06476), and the ν–γ Shapiro-delay test
  (NF-H03). Only NF-H04, delayed "ghost" copies at unrelated positions, is CONDITIONALLY VIABLE, as phenomenology
  (ladder B). The data audit confirms IceTracks-DR2 (doi:10.7910/DVN/MMIIZA, 1.64 M tracks, CC0) and that every
  audited host is reachable; track catalogs carry no flavor.
- **E-NF1 (NF-H04), ICECAT-1 v4, 340 alerts:** wide pairs in 7 lag bins out to 180 d, counts and signalness
  products (14 cells), 20,000 scrambles. It reproduces D-074 exactly (36 and 138). **Null:** pooled global
  p = 0.058; the largest cell is 1 h–1 d weighted, z = 2.7, Bonferroni p = 0.094. The 36 pairs were inspected:
  they span all 13 years, 6 are same-run and two are alert triplets (ordinary day-scale clustering). 95 % CLs
  limits on R_g per astrophysical alert: ≤ 0.02 (Δt ≤ 10 s), ≤ 0.05 (10 s–1 h), ≤ 0.2 (1 h–1 d), ≤ 0.1
  (1–7 d), ≤ 0.5 (7–30 d); at 30–180 d only R_g = 1 is excluded. Label: statistical null; a parameter limit on an ad-hoc model.
- `/code-review` (10 findings, 9 fixed): manifest written only after a successful release; deterministic pin
  choice; refusal of an underpowered `--n`; CLs limit needs every larger R_g to pass; durable copy accepted by
  content digest; NaN distance maps for unmapped GW events; unused masks removed; ghost Dec and error from one
  event. Not changed: `wide_pair_stats` duplicates `count_channel`'s binning (it adds weights; D-074 counts
  reproduce exactly).
- Failed approaches (E-NF1): classical limits excluded R_g = 0 on a low fluctuation, so CLs replaced them.
  Injections tested against the observation plus ghosts gave limits that were too tight; they now use a null draw
  plus the ghost excess. Without subtracting the random pairs that extra ghost events add, sensitivity at large R_g
  was overstated.
- **Speed-ups (Part 0)**, measured on this 4-core cloud machine:

  | Item | Before | After | Gain |
  |---|---|---|---|
  | S-2 `classify_gw_gw` per scramble (≈ 644 real GW–GW pairs, 30 seeds, identical classes) | 808 ms | 1.79 ms (numba); 28 ms (numpy fallback) | ×452; ×29 |
  | S-2 GW–GW pairs/s | 797 | 360,419 | ×452 |
  | E1 GW `count_all` per scramble (4 channels) | ≈ 0.8 s | 8.1 ms | ≈ ×100 |
  | `GWMaps` build (skip the 109 unmapped events' distance maps) | 150 s (36 maps + 355 uniform) | 27–29 s (282 maps) | ×5 |
  | E1 GW directional, full run (`--n 1000` null + injections + map build, same inputs, 60 cells identical, global p 0.833) | 333 s | 52 s | ×6.4 |
  | S-4 GW sky-map tarballs (704 MB, parallel ranges into memory, per-tarball cache) | > 40 min, unfinished (2026-10-09) | 159 s + 122 s for GWTC-4.1/5.0 (1.7–2.2 MB/s); cached re-run 1 s | — |
  | E-NF1 null (ICECAT, 14 cells; jit pass counted only to 7 d after review) | 1,375 scrambles/s | 1,982 scrambles/s (4 cores) | ×1.4 |
  | S-1 GW sky maps in a fresh session (owner's machine, empty data root, `derived-data-20261010`) | ~5 min to > 40 min of streaming (675 s on the owner's machine) | 2.7 s wall | ≈ ×250 |
- S-1: the store (`derived_store.py`, `derived_publish.py`, durable-first fetch in `e1_gw_skymaps.py`) is built and
  tested. **Release creation from cloud sessions is refused** (HTTP 403, "not permitted for this session type"),
  so the first asset waits for the owner or a local session. Today's re-reduction reproduces the pinned maps
  digest exactly (c5c2e4ff…).
- Found: the live HEASARC GBM table (8cc2…, 4,391 rows) matches neither pinned snapshot (D-074 5b42…, #122
  0cb3…), so both GBM-dependent E1 scripts stop in a fresh session until the snapshot is in the durable store or
  re-pinned. `nf_ghost.py` now fetches only ICECAT-1.
- S-5: the claim rule is now a 10-minute heartbeat and a takeover only after 60 minutes with no commit and no
  heartbeat (docs/cloud-routine-prompt.md, research-cycle skill, docs/operations.md).
- **Next:** publish `derived-data-20261010` (GW maps plus the GBM snapshot) from a session that may create
  releases; E-NF1 on IceTracks-DR2 (stream per-season files, good-run null); TXS 2014–15 SkyLLH benchmark first.

## 2026-10-10: W3 MOA-II: the four D-068 chain gaps closed as one post-hoc change (D-068 addendum, CHAIN_VERSION 2026-10-10.1)
- Hypothesis tested: the chain's survivors gb17-R-6-1-3829 (periodic variable) and gb18-R-9-4-24509 (slow dip with
  season levels) got through four gaps; new tests `exotic_chi2_cap`, `residual_deficit`, `periodic_variable`,
  `slow_dip_seasons` (thresholds in `Params`, ASSUMPTIONs).
- Validation: 6 synthetic W3 events in white noise pass every new test; unit tests for each model.
- Real-data controls (`vet` re-run, local only; tracked summaries wait for `run-field`): gb17 75 flags, gb18 67
  flags, **0 survivors**. gb17-R-6-1-3829 now fails `exotic_chi2_cap` (6.67 > 5.46); gb18-R-9-4-24509 fails
  `slow_dip_seasons` (ΔBIC +8.6, σ 415 d, p 5.2: the vetting note's model, found independently). First failures
  at the new tests: `residual_deficit` 13 + 11, `exotic_chi2_cap` 6 + 3, `slow_dip_seasons` 0 + 2.
- Failed approach: residual tests over every epoch failed a synthetic W3 event fitted at another u0 optimum
  (χ²/dof 1.69 vs 1.06, all on the spike nights): caustic-spike epochs are now left out. `/code-review` fixes:
  bright-source wings not spikes, whole-curve residual scan (window-only scans lost the red-noise normalisation),
  F-test error scale, P ≥ 2.5 d (nightly Nyquist), untestable cases named, no duplicate season refit.
  Not changed: `fit_slow_dip_seasons` duplicates `fit_smooth_dip`'s optimiser loop.
- **Next:** `run-field` for every field under the new chain (injections + limits; ~20 min per field), then `combine`.
  Every 2026-10-09.1 vetting record is stale; no limit may be quoted until then.

## 2026-10-10: D1 kinematic D_s/D_ds on 76 lenses: null, max 2.49σ vs 5.78σ (D-073 addendum 4)
- Hypothesis: one lens's kinematics demand a D_s/D_ds no FLRW model gives its redshifts. Data: the hierArc
  kinematic likelihood pickles in TDCOSMO 2025 (TDCOSMO 8, SLACS KCWI 13, SLACS SDSS 41, SL2S 14), read numpy-only.
- `d1_distance.py kinematic` (~1 min): LOO pull vs the others' offset and intrinsic scatter (τ = 0.15 in ln). Max
  SDSSJ2302−0840 −2.49σ; shuffled-z null median 2.71σ (p = 0.79); upper prior bound ×10 changes nothing.
  Detectable ×0.21 / ×5.6 (median lens; KCWI ×0.30–0.45 / ×1.6–2.9). Plot `results/d1_distance/kinematic_lenses.png` checked: SL2S sits low as a sample
  (δ −0.36 ± 0.15), KCWI vs SDSS differ by 0.2 on the same lenses; sample systematics, no single-lens outlier.
- Failed approach: mean/sd summaries and a grid-edge refusal; the hierArc term has a power-law upper tail
  (error ∝ prediction), so quantiles and a stated grid bound (D_s/D_ds ≤ 40) replace them. `/code-review`: γ_pl grid
  0.1 → 0.025 steps (moved SDSSJ1538 2.56σ → 2.48σ), axes by name, ASSUMPTIONs into `Params.kin_*`, τ-bound flag.
  `_ZGRID` now reaches z = 5 (SL2S sources to z 3.35).
- **Next:** W3 MOA-II chain gaps and re-injection (D-068); D1 is exhausted on public data until per-lens λ_int-free
  kinematics (JWST/KCWI IFU) or new time-delay lenses appear.

## 2026-10-10: D1 handoff: hierArc kinematic pickles are readable with the existing safe unpickler
- RXJ1131's `*_const_processed.pkl` is a `DdtHistKin` dict with numpy-only globals; details and next step in TASKS
  "Now" 0 (D1). Not yet used for any result.

## 2026-10-10: D1 composite model on SDSS1206: null; model choice moves its pull 0.61σ → 0.54σ (D-073 addendum 3)
- Hypothesis: the D-073 null depends on the power-law mass model; a composite (stars + NFW) model could shift one lens
  out of line. Only SDSS1206 has a public composite D_dt chain in TDCOSMO 2025.
- `d1_distance.py tdcosmo --j1206 final_composite` (κ_ext included in the chain): J1206 0.54σ, max |pull| 1.14σ vs
  5.39σ, null. Validation first: TDCOSMO's final power-law chain gives 0.606σ vs this pipeline's 0.608σ.
- Failed approach: the pre-LOS composite pickle cannot be paired per sample (κ_pert 20,000 vs 400,000 rows).
- **Next:** new per-lens D_d via a numpy-only reader for the hierArc `*_const_processed.pkl`; W3 MOA-II chain gaps
  and re-injection (D-068).

## 2026-10-10: W3 MOA-II eleven fields: 582 flags, 0 candidates; limits withdrawn pending the corrected chain (D-068, #95)
- Lands the streamed MOA-II work: gb7, gb11, gb12, gb15–gb22 (712,780 Cut-0 light curves, 582 flags). Nine fields
  have 0 chain survivors; gb17-R-6-1-3829 (periodic variable) and gb18-R-9-4-24509 (red-giant slow dimming) passed
  the chain and are explained only by vetting notes (docs/candidates/). `derived`. No W3 candidate.
- Merged `claude/w3-moa-fields-alt2` (injections vetted against the injected source flux, `CHAIN_VERSION`, refusal
  without MulensModel, truncated-tar and size-pin fixes). Every old-chain `limits_*` and `efficiency_*` table is
  deleted: none may be quoted.
- `/code-review` fixes: the vetting record carries `CHAIN_VERSION` and `limit` refuses a mismatch (the eleven
  tracked records predate the stamp, so every field is re-vetted before a limit); a non-finite night-jackknife refit
  now fails the test (it returned −inf and passed; this could only have kept flags, so 0 survivors stands);
  `fit_step_ramp` without an admissible step time; ragged rows go to the line parser; chunk ranges are checked.
  Not changed (follow-ups): duplicate check over all streamed ids, `prefetch` floor above `--conns`, repeated chunk
  reads, per-field metadata pass, unnamed quiet-carrier thresholds in `is_quiet`.
- **Handoff / next:** add the four chain gaps in D-068 as one post-hoc change, re-inject every field, `combine`;
  gb13 stopped at 25 / 78 fits (run died 2026-10-09 13:02Z); then more Nunota Table 1 fields.

## 2026-10-10: E1 signed-lag GW cells with sky-map classes: null, global p = 0.56 (D-074 addendum 5)
- Hypothesis: GW events lead (or follow) GBM, ICECAT-1 or CHIME events in a given sky-map class (same, wide,
  antipodal) at lags no ordinary path explains: D = N_after − N_before.
- `scripts/e1_gw_signed.py --n 1000` (71 s): 45 cells (3 channels × 5 lags × 3 classes; GW–GW left out because D
  is antisymmetric there). Pooled global p = 0.56, min analytic p × 45 = 0.78. Largest: GW–CHIME 1 h–1 d antipodal,
  D = 10 vs −0.3 ± 4.3 (z = 2.4; it was z = 3.1 with 40 scrambles: the smoke-run null sd was too small).
- Control GW170817 → GRB 170817A gives D = +1 at 0–10 s same; it is left out of the family.
- Sensitivity: 3 one-sided pairs at ≤ 10 s, 3–10 at 10–100 s, 10–30 at 100 s–1 h, mostly > 30 at 1 h–7 d. Per-cell 95 %
  upper limits on extra after-pairs are in `results/e1_events/gw_signed.json`. Cells are tested against the null mean
  (catalogue edges make it non-zero, e.g. GW–CHIME 1–7 d wide −66.5 ± 52.7).
- **Next:** D1 composite-model D_dt chains; E1 IceTracks-DR2 / Swift / Einstein Probe when they can be fetched.

## 2026-10-10: E1 GW channels with sky maps: null, global p = 0.83 (D-074 addendum 4)
- Hypothesis: GW events have same-direction, wide or antipodal partners (GBM, ICECAT-1, CHIME, GW) at lags no
  ordinary path explains. Until now the GW channels were lag-only.
- `scripts/e1_gw_skymaps.py` streamed the four GWTC PE sky-map tarballs (704 MB, about 3.5 min, one process per `--tar` then a merge run; the
  2026-10-09 run got 0.15–1 MB/s and did not finish) into 282 nside-32 maps. The 109 unmapped GWTC-4.1 / 5.0 entries
  have no PE parameters in the GWOSC CSV either.
- `scripts/e1_gw_directional.py --n 1000` (266 s, 4 cores): 60 cells, pooled global p = 0.83, min analytic p × 60 = 1.
  Largest z: GW–ICECAT 1–7 d antipodal (14 vs 8.3 ± 2.9). Control GW170817 × GRB 170817A is `same` at 0–10 s; it
  is left out of the family and reported apart.
- Sensitivity: 3 injected pairs at ≤ 10 s, 3–10 at 10–100 s, 10–30 at 100 s–1 h, mostly > 30 at 1 h–7 d.
- Catalogue snapshot pinned separately (`data/manifests/e1_gw_events.ecsv`): the GBM TAP table grew since the
  D-074 pin, so `e1_events.fetch` refuses the old one. Rerunning older E1 scripts needs `--refresh` or that pin.
- **Next:** signed-lag GW cells with same/wide classes; then D1 composite-model chains.

## 2026-10-09: E1 signed-lag ("which event comes first") channels: null, global p = 0.68 (D-074 addendum 3)
- Hypothesis: one catalogue leads another at lags no ordinary path explains (sign asymmetry D = N_after − N_before).
- `scripts/e1_signed_lag.py`: 45 cells (6 cross channels × 5 lags × class), `jit` null (2,000 scrambles, 52 s).
  Pooled global p = 0.68. The largest |z| is GW170817 → GRB 170817A (D = −1, Skellam p = 0.05), the known ordinary
  ordering.
- Failed approach: a Gaussian tail on one-pair cells (z = −4.4 for GW170817 alone, Bonferroni 6 × 10⁻⁴). Replaced by a
  Skellam tail where the null sd < 1. The empirical p could not be used instead: its floor (1/2001 × 45) makes
  injections undetectable.
- Sensitivity: 3–10 one-sided injected pairs at ≤ 100 s.
- **Next:** E1 GW sky maps (GW cells get directions; redo antipodal and signed lags), then D1 composite-model chains.

## 2026-10-09: E1 antipodal lag channels: null, global p = 0.56 (D-074 addendum 2)
- Hypothesis: a dependent event appears near the antipode of a first event, at 0 s–7 d lags. Ordinary: the D-074
  observer-side common causes.
- `scripts/e1_antipodal.py`: 180° − sep ≤ max(10°, 3σ_comb); 30 cells (6 localized channels × 5 lag bins); `jit` null
  (2,000 scrambles, 19 s on 4 cores). Pooled global p = 0.56. Max z = 1.7 (CHIME–GBM 1 h–1 d, 258 vs 230 ± 16).
- Injection-calibrated: 3 pairs at ≤ 100 s detected at family-wise 3σ in most channels, 3–30 at 100 s–1 h, 10–100 at 1 h–7 d.
- Fixed before the run: a same-catalogue channel must count the injected sample with itself (`pairs_within`
  assumes ta is tb); counting the original against the injected copy undercounted.
- **Next:** E1 GW sky maps (GW–X same/wide/antipodal channels), signed-lag ("which comes first") channels.

## 2026-10-09: E1 CHIME 100 s–1 h wide cell under the calibrated rate-modulated null: z 3.7 → 2.5, family-wise null (D-074 addendum)
- Hypothesis: the 100 s–1 h CHIME–CHIME wide excess (396 vs 330, z = 3.7 `jit`) is lag dependence beyond a smooth
  detection-rate modulation. Ordinary: week-scale rate modulation, daily duty cycle, sub-day outages.
- `scripts/e1_chime_flag.py --cell 100s-1h` (2,000 scrambles, 30 s): 7-day null z = 2.52 (analytic p × 45 = 0.26);
  time of day kept z = 2.42; 3-day z = 0.48 (not primary, D-074 (d)). The 7-day null keeps ≥ 97 % of 300 injected pairs.
- Result: null at the family level; a local 2.5σ residual remains that only hour-scale uptime could decide.
- Failed approach (not run, reasoned): an hour-resolved running mean of the catalogue's own counts would absorb the
  tested pairs, as the in-day null did.
- **Next:** E1 GW sky maps (GW channels get directions), then "which event comes first" / antipodal channels.

## 2026-10-09: Hypothesis round 3: derived small effects in large-N data; no testable survivor; new rounds paused (D-076)
- Brief: a derived small amplitude in a large-N public quantity never analysed for it. **Results:**
  - R3-A (same-exposure pair correlation, about 2.4 × 10⁻⁵): fails as written (axiom A5 gives zero). The Holometer
    disfavours it at about 67 % (reviewer estimate).
  - R3-B (condensate scalar charge in neutron stars): upper half (s ≥ 8 × 10⁻³) excluded at 6.7σ by J1738+0333; self-consistent
    minimum disfavoured at about 2σ; only the no-feedback corner survives. α_NS < 2.2 × 10⁻³ (95 %, reviewer's combination).
  - R3-C (SN Ia rate vs "clock depth"): not excluded but untestable. About 5 × 10³ public SNe Ia give 0.5–2σ, and
    environment systematics are 10–40× the signal (reviewer estimates).
- **Rule:** three rounds (9 worlds) produced survivors only in round 1, and those were null in data. New rounds are
  paused; compute goes to data-limited open searches (D-076).

## 2026-10-09: D1 on TDCOSMO 2025 power-law chains: 8 lenses null with κ_ext; without it RXJ1131 is 5.7σ off (D-073 addendum 2)
- Hypothesis: one lens's D_dt disagrees with the H0 + Ωm of the others (a sightline-specific distance anomaly).
  Ordinary explanations: κ_ext, internal mass sheet / model choice, time-delay systematics. A flag with κ_ext applied
  would have gone to composite models and κ_ext vetting.
- `scripts/d1_distance.py tdcosmo`: power-law D_dt^model chains (6 H0LiCOW lenses + DES0408 + WGD2038), each lens's
  own TDCOSMO κ_ext PDF applied identically, H0 + Ωm LOO; 100 shuffled-z nulls; injections.
- **With κ_ext: max |pull| 1.12σ (RXJ1131), threshold 5.39σ; null.** Reach ×0.41–0.78 / ×1.36–1.85 (WGD2038 ×0.26 /
  ×3.2). LOO H0 72.2–74.4. κ_ext files checked uncorrelated with their chains (|r| < 0.003).
- **Without κ_ext: RXJ1131 −5.67σ, DES0408 +3.9σ**, both removed by the measured κ_ext (positive control: the test
  sees a ~7 % line-of-sight convergence).
- Safety: SDSS1206's pre-LOS file is a pickle; read with a numpy-only unpickler (arbitrary classes refused, tested).
- Data: whole TDCOSMO2025_public checkout (661 MB, cloud) streamed by git, used, deleted; 16 files pinned in
  `data/manifests/d1_distance.ecsv`. Failed approach: `git clone --no-checkout` then `git checkout <sha>` fetches the
  whole tree; check out paths only.
- **Next:** composite-model chains (TDCOSMO SDSS1206 `final_composite_*`, others where released) for a model-choice
  check; per-lens D_d needs a data-only reader for `*_const_processed.pkl`.

## 2026-10-09: Hypothesis round 2: three worlds with derived Lorentz invariance; nothing testable here survives (D-075)
- System A (no literature) was bound by the round-1 failure modes:
  - derive Lorentz invariance;
  - no long-range links;
  - no achromatic dimming;
  - only derived amplitudes.
- **Results:**
  - R2-A (every two records share a future, so horizons re-emit everything): **excluded**. Ergoregion instability;
    GWTC-4.0 null; supporting: GW250114 reflectivity < 0.35 % (preprint).
  - R2-B (gravity as phase inference, so saddle images mirror polarisation): internally inconsistent; B0218+357
    polarisation variations correlate with the wrong sign.
  - R2-C (four axioms; massless lightest neutrino, so Σm_ν = 58.8 meV): not novel (minimal seesaw). 2.2–3.4σ
    tension with DESI DR2 + CMB (reviewer estimates; DESI quotes 3.0σ); tracked externally.
- **Lesson (rule):** deriving Lorentz invariance from a symmetric primitive is easy. The theory's risk then sits in
  one precision-measured sector, where O(1) effects are already excluded. Round 3 should look for derived small
  effects in large-N quantities that nobody has analysed for them.

## 2026-10-09: E1 causal event network: no wide-separation dependence between GRBs, neutrinos, GW events and FRBs (D-074)
- Hypothesis (owner idea 4): events in different directions that depend on each other at lags no ordinary path
  explains. Ordinary explanations: observer-side common causes (uptime, exposure, Sun, follow-up chains, shared
  triggers, duplicates).
- Data (4.7 MB, pinned): Fermi GBM 4,390; ICECAT-1 340; GWTC 391 (no sky positions, so lag only); CHIME/FRB Cat 2
  3,641 sources. Lag × separation pair counts in 75 cells. Nulls keep declination and hour angle; GBM is shifted in
  whole orbits.
- **Five requested channels: family p = 0.28.** All wide/GW cells: pooled family p = 0.0058 (jit) / 0.0033 (jitday).
  The driver, CHIME–CHIME 1 h–1 d wide (z ≈ 3.9; 5.8 per cell in the post-hoc nominal-flag subset), disappears
  under a calibrated rate-modulated null (7-day running mean of daily counts, keeps 87 % of an injected signal):
  z = −0.86. Same-direction CHIME excess: one unflagged six-burst episode (FRB20230825D–I).
- Positive controls: GBM re-triggers recovered. GW170817 × GRB 170817A is the only GBM–GW pair within 10 s, but
  blind lag counting gives only p = 0.058 without GW sky maps.
- Limits (95 %, family-wise 3σ detection, anchors inside the partner's live time) in `results/e1_events/limits.ecsv`.
- **Failed approaches (rules):**
  - Observed and scrambled p computed with different rank formulas lost the trials factor beyond the ensemble.
    Pool the observation with the scrambles.
  - A "keep each event's day" null removes ~85 % of an injected dependent signal, so it is not a vetting test.
  - Resampling times with replacement self-pairs events.
  - The CHIME Cat 2 exposure file (216 MB, downloaded and deleted) has no time axis, so it cannot model uptime.

## 2026-10-09: D1 handoff: TDCOSMO 2025 per-lens files listed
- File names in TDCOSMO2025_public d7f38db recorded in SOURCES.md: D_dt chains for the 6 D-073 lenses plus DES0408 and
  WGD2038; new per-lens D_d only inside hierArc likelihoods. Next: TASKS "Now" 0.
- An equal-weight NE2001/YMW16 mixture cannot flag more than the per-model runs (its tail is their average); a wider
  ISM error only lowers the pulls. Neither can turn the FRB null into a flag, so both wait.

## 2026-10-09: D1 FRBs with YMW16 beside NE2001: still null (D-073 addendum)
- Question: is the D-073 FRB null (and its 3.26σ extreme) an artefact of NE2001? Ordinary explanation for any
  low-side outlier: Galactic-model error at low latitude. A flag under YMW16 alone would have needed vetting.
- `scripts/d1_distance.py frb --ism ymw16`: YMW16 to 30 kpc per RA/Dec, same chain and 5.81σ threshold. Known case
  FRB 20121102A: 287 pc cm⁻³ (published YMW16 value). NE2001 path re-run: identical pulls and flags.
- **Result: no flag under YMW16 or under either model.** Max low pull FRB 20220319D 4.09σ (NE2001: 3.26σ;
  DM_obs 111 below both models, 127 and 211); max high 20190520B 2.37σ; injections 94/94 per side.
  YMW16/NE2001 median 0.84, range 0.36–1.91; > 20 % apart for 54/94 bursts. Outputs `results/d1_distance/*_ymw16*`,
  `frb_ism_compare.ecsv`.
- **Failed route:** `pip install pygedm` fails in the cloud image (its NE2001 extension needs libf2c) and the
  package import fails on current SciPy (`integrate.simps`); built without `ne21c`, `ymw16` imported directly.
- **Next:** TASKS "Now" 0 (per-lens TDCOSMO D_d; a model-marginalised ISM term).

## 2026-10-09: S2 on DES-SN5YR: no flat-kernel term; γ < 0.019 mag per unit T/⟨T⟩ (D-070 addendum)
- Same chain as D-070 (`scripts/s2_flat_kernel.py --sample des`, `results/s2_flat_kernel/fit_des_*.json`) on the
  DES-SN5YR Dovekie Hubble diagram: DES-discovered SNe only, MU/MUERR (BEAMS-renormalised), P(Ia) ≥ 0.5 (ASSUMPTION),
  0.1 < z < 1.3: 1,518 SNe, 1,414 fitted after coverage and 5σ cuts; LS DR9 galaxies (47.2 k) as before.
- **α = 2: γ_F = +0.0009 ± 0.0025 (scramble p = 0.68); trimmed −0.0012 ± 0.0044; colour term +0.0008 ± 0.0018.**
  λ_F = 0.86 (0.81 trimmed); whole-chain recovery 0.57 (limit divided by it). **One-sided 95 %: γ < 0.019 mag per
  unit T/⟨T⟩** (Pantheon+: 0.025). Still above A3's kill threshold 0.005; excludes the upper half of its 0.01–0.04
  fiducial only if the galaxy-to-matter bias is ~1 (ASSUMPTION).
- α = 1: kernel-degenerate again (corr(X_L, X_F) = 0.99), no limit, as in D-070.
- **Positive control:** γ_L = −0.018 ± 0.011 (lensing sign, scramble p = 0.095; λ_L = 0.60). Pantheon+ gave
  −0.017 ± 0.013; naively combined −0.018 ± 0.008 (~2σ), but the samples overlap (DES 3YR spectroscopic SNe are in
  Pantheon+), so this is indicative only. The control is still not a detection with z < 21 LS DR9 galaxies.
- Not done (failed route): the DES Y3 Gold / deeper-galaxy column the D-070 revisit asked for; LS DR9 at z < 21 was
  kept so that the two samples share one calibrated chain.
- **Next:** a deeper galaxy column (DES Y3 Gold or LS DR10 z < 22 with photo-z quality cuts) so the lensing control
  detects; the full DES-SN5YR STAT+SYS covariance (inverse matrices in the repo) instead of diagonal errors.

## 2026-10-09: D1 distance self-consistency: no inconsistent sightline in strong lenses or localized FRBs (D-073)
- Hypothesis (owner idea 1): one sightline whose independent distance measures no ordinary model reconciles.
  Ordinary explanations: mass-sheet degeneracy, line-of-sight convergence, kinematic anisotropy, substructure, wrong
  redshifts; for FRBs, Galactic ISM model error and host DM.
- Lenses: H0LiCOW public posteriors. The H0-free R = D_dt/((1+z_d)D_d) = D_s/D_ds for B1608, RXJ1131, PG1115, J1206
  (prior-predictive ΛCDM / wCDM and leave-one-out), plus leave-one-out D_dt for those and HE0435, WFI2033. **Max pull
  1.29σ** against 5.53σ (18 trials). Shuffled redshifts give a median max pull of 3.8σ (C) / 15σ (D).
- Injection-calibrated reach (|z| crossing the threshold, baseline pull included): R tests ×2–8 up / ×0.2–0.4 down;
  D_dt ×1.3–1.85 / ×0.54–0.78. A ≲ 30 % mismatch on one sightline is not excluded.
- FRBs: 94 localized FRBs (FRBs/FRB repo) against the Macquart predictive distribution (grid convolution, tails to
  ≪ 1e-10). Nothing flagged at the one-sided 5.81σ; the lowest is FRB 20220319D at 3.26σ (DM below its NE2001 DM_ISM,
  a known low-latitude case). Injections flag 94/94 low and 94/94 high.
- **Failed approaches (rules):**
  - A Monte Carlo predictive tail floored at 1/N capped z at 4.26σ below a 5.93σ threshold, so the screen could
    never flag (found by review). Tails must reach past the threshold: test it.
  - Sensitivity factors from a Gaussian extrapolation ignored baseline pulls. Read them from injections.
  - TDCOSMO 2025's SDSS1206 D_d file duplicates H0LiCOW's samples and adds no lens.
- **Next:** see TASKS "Now" 0.

## 2026-10-09: S3 hybrid images: no copies in 6 COSMOGRAIL doubles; sensitivity only r ≈ 0.4–2.4 (D-072)
- Hypothesis (A1 P2): image j carries a faint copy at lag s_i − s_j; the copy into the trailing image can precede the
  leading image. Ordinary mimics: quasar red noise (correlates at all lags), delay errors, microlensing.
- Data: COSMOGRAIL XIX (CDS J/A+A/640/A105), 7 doubles with secure delays. SIS lag window for f ∈ [0, 1]
  (D_ls/D_s = 0.46–0.68). `scripts/s3_hybrid.py`, `results/s3_hybrid/`.
- Known-case gate: the main-term fit recovers the published delay (sign, ≠ 0, within 3σ + 4 d) in 6 of 7. J1620 fails
  with a wrong-sign minimum at τ_BA = −206 d (published +171.5 d) and is not screened.
- **Result:** no copy. Window maxima against same-width off-model windows (|lag| < 10 d excluded from both): lowest
  p = 0.02 (J1226 into-trail, a window on the visible main-image wing) of 12 windows; the windows overlap, so p is
  approximate and the run-wide chance of p ≤ 0.02 is ~0.2. Sensitivity r95 ≈ 0.38 (J1455 into-trail) to 2.4: only
  copies comparable to the main image are excluded, far from B-on-A1's ~10⁻². Not quoted as limits: the
  window-maximum injection efficiency is 0.67–4.7 (> 1 is the copy–main degeneracy; capped at 1 in r95), the
  null windows are not matched in |lag| (windows next to the wing get low p), and for HE0047 the delay gate
  (±14.5 d around 10.4 d) is weak.
- **Failed approaches (rules):** microlensing splines with 120–730-d knots absorb the quasar variability and fit any
  delay (6 of 7 published delays missed); use the delay recovery as the gate for every setting. Quasar variability
  is too slow for copies at lags ≲ 100 d: the copy is collinear with the main image (a synthetic 20-d DRW recovers
  r = 0.1 ± 0.03, so the loss is the data). /code-review: inject through the window maximum, not the fit at the
  injected lag; hold the observed window to the null's wing rule; inject a smoothed copy (the template's noise
  otherwise rides along); fit the main term at the validated delay.
- **Next:** S3 needs a sharp template: SN Refsdal / SN H0pe pre-explosion and post-peak HST/JWST imaging at the
  model lags (B-on-A1: a few × 10⁻²). Low priority; S1 and S2 (DES-SN5YR) first.

## 2026-10-09: S1 burst twins in Fermi GBM: null; < 7.2 × 10⁻³ twin pairs per eligible burst (D-071)
- Hypothesis (D-069 S1, from A2 P1): pairs of GRBs far apart on the sky whose light curves match at any delay.
  Ordinary mimics: single-envelope look-alikes, re-triggers, ordinary lensing (excluded by the position cut).
- Data: HEASARC `fermigbrst` (4,390 bursts) and bcat HDU 2 light curves (4,389; streamed, sha256 recorded).
  `scripts/s1_ingest.py`, `scripts/s1_twins.py`, `results/s1_twins/`.
- 1,414 multi-pulse bursts; 510,294 pairs with inconsistent positions (sep > 3·√(σ₁² + σ₂² + 2·3.7²) deg).
  Pulse-shuffle null ρ* = 0.964 (a pre-screen: the surrogates under-predict the real high-ρ tail ~12×). 6 pairs
  flagged, all bright single-envelope look-alikes failing the χ² twin test (p < 10⁻¹⁰⁰). **0 survivors.**
- Injections with a 15 % per-band mismatch measured from the data: ratio-1 efficiency 29.5 %. **95 %: < 7.2 × 10⁻³
  twin pairs per eligible burst (1.4 % of bursts have a twin) at flux ratio 1, 1.6 × 10⁻² averaged over ratios 1–0.1.**
  Sensitivity: 5.1–9.4 × 10⁻³ for 5–20 % mismatch.
- **Failed approaches (rules):**
  - bcat HDU 1 `PHTCNTS` fill values made 84 % of bursts look multi-pulse.
  - Self-copy injections shared the target's noise (recovery 83 % vs 37 %).
  - Injections whose noise matches the χ² errors exactly pass χ² by construction: inject a measured band mismatch.
  - A ρ threshold chosen so that k = 0 (ρ > 0.90) is post hoc and illustrative, not a limit.
- **Next:** CHIME/FRB Catalog 2 (the proxy resets CANFAR downloads); TTE light curves for short and faint bursts; an
  s ≠ 1 chain with energy-channel mapping; a generative pulse-model null.

## 2026-10-09: S2 flat-kernel SN residuals: no flat-kernel term; γ < 0.025 mag per unit T/⟨T⟩ at α = 2 (D-070)
- Hypothesis (A3 P2b): SN Ia residuals track a flat-kernel foreground column ∫(1+δ)^α dχ, with sign +γ (fainter),
  after the lensing-kernel column is fitted. Ordinary mimics: lensing magnification, grey dust (colour), host-group
  overlap, leverage of a few sightlines, disc coverage, shot noise in the column (regression dilution).
- Data: Pantheon+ (944 unique Hubble-flow SNe, z 0.1–1.3) and LS DR9 photo-z galaxies (countmap selection, dereddened
  z < 21 mag, 2′ discs, 32.3 k galaxies). `scripts/s2_flat_kernel.py`, `results/s2_flat_kernel/`.
- **α = 2 (A3's fiducial), 505 SNe** (shells expecting ≥ 1 galaxy, so z_s ≳ 0.2; coverage and 5σ cuts remove 15):
  γ_F = −0.0003 ± 0.0018 per unit X (scramble p = 0.83); 1 % trimmed +0.0024 ± 0.0043; no colour term
  (−0.0017 ± 0.0012). Dilution λ = 0.92 (0.75 trimmed). Whole-chain injection (γ = 0.01 on full counts, measured on
  half-thinned counts, 20 trials) recovers 1.04 (full) / 0.78 (trimmed) of the input; the limit is divided by 0.78.
  **One-sided 95 %: γ < 0.025 mag per unit galaxy-traced T/⟨T⟩.** A3's fiducial is 0.01–0.04 (B-on-A3 expected a
  ~0.015 reach), and its kill threshold is 0.005. Converting to a matter column (bias b ~ 1–2, ASSUMPTION) weakens the
  limit by up to ~b².
- Lensing column: γ_L = −0.017 ± 0.013 (the lensing sign, p = 0.21). The positive control is not detected, as
  expected at this N (Smith+2014: 1.4σ with 608 SNe).
- **Failed approaches (rules):**
  - The linear (α = 1) flat column is 0.996-correlated with the lensing column in 2′ disc counts. It cannot separate
    the kernels: the chain injection is not recovered, so it sets no limit.
  - N(N−1)/E in shells with E ≪ 1 (z_g ≈ 0.045 in a 2′ disc) gave X ≈ 20–54 from 3–5 galaxies. These sightlines
    drove a spurious "γ < 0.0093" in the first pass. Require E ≥ 1 per shell.
  - Scrambled-column injections recover γ exactly but cannot see regression dilution. Use the Poisson-only λ and the
    thinned-count chain injection. λ on a trimmed sample must trim the simulations too (found by /code-review).
  - Raw `mag_z` with no maskbits/brick_primary/DUP cuts was the first query; use the repository's countmap selection.
- **Next:** DES-SN5YR + DES Y3 Gold (Shah+2024's setup: lensing detected at 6σ, the positive control; deeper
  galaxies give more shells with E ≥ 1); the full Pantheon+ covariance in place of the diagonal errors.

## 2026-10-09: Owner direction: System A/B hypothesis rounds; round 1 (3 worlds, 3 reviews, 3 conditional survivors) (D-069)
- Owner ideas (distance self-consistency, missing light with paired excess, transient connectivity, causal event
  network, Ontological Reset). Distance self-consistency was never tested here; the deficit half of missing light is
  W3/W5.
- System A (no literature): A1 states + transitions; A2 no identity, one quantity "carry"; A3 bridged Markov
  substrate with a fixed pattern of distance disagreement. System B (verified references) failed 14 of 17
  predictions (plus A2's foundations). All three link mechanisms are disordered locality (arXiv:0903.5303) renamed, and achromatic dimming
  is the known cosmic-opacity test.
- **Survivors (conditional):** S1 burst twins (GRB/FRB, wide separation, any delay: untested by 2006.07095 and
  2204.06014); S2 flat-kernel SN residuals (new ~0.015 mag limit reachable); S3 hybrid images (limit only).
- B corrected A1's "blindness theorem": a timing dependence between events ≳ 0.1° apart is provably spacelike at Gpc distances with spectroscopic redshifts (B-est); the
  obstacle is chance coincidences and common causes.
- **Next:** S1 on the Fermi GBM catalogue and CHIME/FRB Catalog 2, with time-scrambled nulls and lensing/duplicates as
  positive controls (docs/hypotheses/round-1/summary.md).
## 2026-10-09: W3 MOA-II six fields null; limits withdrawn after review (D-068, #95)
- Cloud runs streamed gb11, gb16, gb19, gb20, gb21 and gb22 (296,618 Cut-0 light curves; no tar on disk): 234 flags,
  **0 survive** (vetting records for gb20-R-4-0-49379 and gb19-R-4-4-31159; gb11 contact sheet inspected: eclipse
  dips, repeated deficits, unsampled features). New tests: `exotic_in_domain`, `smooth_dip`, `feature_bracketed`,
  `step_ramp`. `derived`.
- **Limits withdrawn:** `/code-review` found that injections were vetted against the lenient 14.2 mag default
  reference flux, while real flags used their own magnitude, so efficiencies were biased high. Injections now use
  I_s; `CHAIN_VERSION` makes `limit`/`combine` refuse old tables. Other fixes: a truncated tar looped forever in
  `read_segment`; a 200 reply skipped the size pin; the Gaia XMatch is retried, and a failure marks flags incomplete.
- **Failed approaches (rules):** a fresh cloud venv without the `mulens` extra skipped the parallax and binary-lens
  tests silently (the stages now refuse to run); `pkill -f` killed the calling shell again (kill by PID).
- **Handoff / next:** a run that died at 03:16 left gb11's injections unfinished; it was taken over at 05:12. Re-run
  `w3_moa.py --field gbN run-field` for the six fields (injections + limit, ~20 min each), then `combine`, in a new
  `[field: W3 MOA-II limits]` PR; then more Nunota Table 1 fields.

## 2026-10-09: W5/W1 Euclid Q1 radial-shear survey: all Deep Fields (60 deg²) null; R calibrated (D-067)
- Hypothesis as D-066 (negative mass shears background galaxies radially). Ordinary radial patterns: PSF-anisotropy
  gradients, blends, tile edges; known clusters must give the opposite (tangential) sign.
- **R = 0.56 ± 0.16** (z_s = 1, scaled by √(χ²/dof); 0.50–0.67 for z_s = 0.8–1.2) from the 4 SZ clusters against
  NFW haloes of their M500: the assumed R = 0.5 is consistent.
- **All 344 Q1 tiles of EDF-F/S/N** fetched by `tileid` (4.78 M rows, 4.60 M of them galaxies; 25 min) and screened at θ_E = 1′, 2′, 4′:
  **0 flags** (field maxima S = 3.3–4.4 against thresholds 4.3–5.1); known clusters 8/9 negative (EDF-S minimum is
  ACT-CL J0405.9-4915). Injection efficiency 0.43–0.58 at 1′ (no limit: below the 0.5 gate), 1.0 at 2′–4′.
  **n₉₅ ≈ 0.049 / 0.051 deg⁻²** at θ_E = 2′ / 4′. With D-063, negative point masses are limited from 2′ to 32′.
- **Failed approaches (rules):** injecting into the observed moments before the PSF deconvolution boosted the
  injected shear ~1.6× over the calibrated R (found by /code-review before merge; the first run's 1′ limit 0.062 deg⁻²
  is withdrawn and the D-066 pilot efficiencies are optimistic): inject where R is calibrated. 16 IRSA threads → 504s,
  no gain; killing a fetch truncated a cache file (writes
  are atomic now); `pgrep -f`/`pkill -f` waiters matched their own shell again (wait on the log file).
- **Next:** per-tile star-ellipticity gradient test (data in `survey.json`); θ_E ≤ 1′ needs fainter or better shapes
  (Euclid DR1); DR1 scale-up when public; θ_E ≈ 1° needs a larger contiguous area than Q1.

## 2026-10-09: J0728+2607 decided with an empirical PSF: no light between the images to F814W ≈ 23 (D-064 addendum)
- The PSF is the median of 6 unsaturated Gaia stars in the same HST F814W cutout. The control lens J2308+3201 is
  detected (S/N 49). J0728's residual is centred on its images, and after the halo correction there is no light
  between them (S/N 1.6). Whole-chain injections are recovered to F814W = 23, ≥ 1.8 mag below any ordinary lens.
  J0130+0725 has 1 usable star, so it has no empirical-PSF result.
- Both D-064 HST pairs therefore have no visible deflector; binary quasars remain the untested ordinary explanation
  (spectra of both images, or time delays, are needed).
- **Failed approach:** the lens-vs-mirror aperture test. The control gives only 3.1σ, because real lens positions
  depart from the SIS flux-ratio rule.
## 2026-10-09: W5/W1 radial-shear screen on Euclid Q1 MER shapes: pilot null, first limits at θ_E = 1–2′ (D-066)
- Hypothesis: a negative-mass lens shears background galaxies radially; ordinary mass gives tangential shear.
- **Convention found wrong in the archive docs:** MER `position_angle` is PA east of north (25 VIS cutouts, median
  |Δ| 0.78°), not THETA_IMAGE "CCW/x". Four SZ clusters in EDF-S then show tangential shear (S = −1.5 to −5.0).
- Pilot: 0.3° discs in EDF-F and EDF-S, θ_E = 30″/1′/2′: **0 flags** (field-max p_random 0.07–0.99); off-grid,
  whole-chain injections (R = 0.5 ASSUMPTION) 0.88–1.0 at 1′, 1.0 at 2′, ≤ 0.18 at 30″. n₉₅ ≈ 8.1 deg⁻² (1′),
  12 deg⁻² (2′).
  The DR10 count floor was 6′.
- **Failed approaches (rules):** unquoted `position_angle` breaks IRSA's ADQL parser (quote it); IRSA returns query
  errors as a VOTable with HTTP 200 (check the body); the EDF-N 0.3° row query hung > 15 min; IBE cutouts come
  gzip-compressed; `pkill -f` on a pattern in your own command kills your shell.
- **Next:** calibrate R on the SZ clusters (shear vs. their M500 NFW prediction); scale to all Q1 (63 deg²; fetch
  in ≤ 0.2° discs in parallel, EDF-N too); star-ellipticity maps for PSF gradients; flag vetting via `/vet-candidate`.

## 2026-10-09: Efficiency rules 8–13 for routine cycles (owner text)
- docs/cloud-routine-prompt.md "EFFICIENCY RULES" gets the owner's rules 8–13: merge main before the final review;
  state-file conflict handling (own CHANGELOG block first, both entries kept); one CI wait per head SHA; verify
  outcomes on GitHub; agents and the owner share the DD-Ching account; stop when the next unit needs a human.
- **Owner action:** paste the updated prompt into the routine (`trig_01PNAmgcfqef8CvhPAY8ggbP`).

## 2026-10-09: D-064 addendum: archival HST of the untestable rejected pairs; J0130+0725 has no lens light to F814W ≈ 23
- Coverage: 2 of the 11 untestable pairs have HST F814W imaging (program 17308), and so do 2 controls.
- Method: a two-PSF fit, halo correction, and the residual flux between the images, with empirical sky-aperture
  errors. Validated on both controls (S/N 42, 16) and on whole-chain injections (recovered to F814W = 23).
- **J0130+0725:** no residual (S/N 0.45). The limit F814W ≈ 23.0 is ≈ 3.5 mag below a typical ordinary lens and
  ≥ 0.7 mag below a 2σ under-luminous one. A binary quasar remains the untested ordinary explanation; it needs
  spectra of both images. Not a candidate.
- **J0728+2607:** inconclusive (residuals at the image cores); it needs an empirical PSF.
- **Failed approaches:**
  - the HAP "combined_skycells" product has a WCS that does not describe its pixels;
  - without the halo correction, PSF mismatch alone gives S/N 18–24;
  - white-noise errors on drizzled pixels overstate S/N by ~3.5;
  - a 60-px cutout cannot hold the pair, because the catalogue position is one image.
- **Next:** spectra or two-epoch flux ratios for J0130+0725 (binary vs lens); an empirical PSF for J0728+2607;
  Euclid DR1 when public for the 9 pairs without HST imaging.
## 2026-10-09: W5 in Euclid Q1: deeper counts cannot open θ_E < 6′; shapes could (D-065)
- Hypothesis: Euclid Q1's deeper counts lower the D-063 floor (blind below θ_E ≈ 6′). Measured (IRSA TAP counts,
  three deep fields): VIS < 24.5 extended galaxies 6.8 × 10⁴ deg⁻², 1.81× DR10, same count slope. With the
  clustering-inflated scatter measured in DR10, the count S/N gains only ×1.0–1.4 (Z scatter at 4′ is 1.6× Poisson
  and clustering does not shrink with depth), on 63 instead of 340 deg². Count screen not built.
- Forecast (model_prediction): a radial-shear test on the same galaxies reaches θ_E ≈ 26″ (1.5–3 θ_E annulus, S/N ≥ 6,
  σ_γ = 0.3, ASSUMPTIONs); the shear sign separates a negative-mass lens from ordinary foreground mass.
- **Failed approach (rule):** IRSA TAP does not index plain RA/Dec ranges (0.25 deg² box > 5 min); use
  `CONTAINS(POINT, CIRCLE)` (40 s for 12 k rows).
- **Next:** Euclid Q1 radial-shear screen: PSF-anisotropy check on stars, synthetic shear injections, trial centres
  on a grid, cross-field null (EDF-N/F/S).

## 2026-10-08: D-056 amendment — the LS pair must be the catalogued pair; one entry per lens (W1/W2)
- `lenscats.pair_match` (moved from D-064's `w12_niq`) now gates the D-056 quasar pair test: a status from an LS image
  pair whose separation differs from the catalogued 2θ_E by > 0.5″ is undecided. Only 115252+004733 changed (LS pair
  4.18″ vs 3.34″ expected; the cutout shows an 18.1 mag lens galaxy with an unrelated faint pair). None of the 15
  pair-decided systems has a catalogued θ_E, so the pipeline checks 0 of 15; a one-off manual check of the 3 also in
  SQLS matches (1.88/2.99/2.03″ vs LS 2.00/3.01/2.01″). For quads (≥ 3 LS images) only a pair wider than 2θ_E + 0.5″
  is a mismatch (a fold/cusp pair is closer than 2θ_E); the check applies only to "deflector"/"none" statuses.
- Same lens listed twice beyond the 3″ merge (MG0414+0534, B2114+022, B2319+052; ~11″ apart): decided systems are
  grouped by designation within 30″ (ASSUMPTION; decimal-degree names give no key) over all covered sensitive
  systems, and a deflector at any copy (decided or not) explains the lens. Each had counted once as
  "deflector" and once as "none".
- New limits (typical, `derived`): quasar 3/15 < 0.52, radio 0/10 < 0.30, all **3/25 < 0.31**; conservative 0/5 < 0.60
  (were 0.48 / 0.23 / 0.27 / 0.50). No new unexplained system; the 3 CHITAH pairs stay open as before.
- **Rules:** a pair test must check that its pair is the catalogued one; name-match decided systems before counting N.

## 2026-10-08: Efficiency rules for routine cycles (owner text)
- docs/cloud-routine-prompt.md gets the owner's "EFFICIENCY RULES" after "MOVE FAST, SAFELY": result first, a review
  stopping rule (fix only result/provenance/reproducibility/guarded-file findings, list the rest, ≤ 3 rounds), verify
  the branch after forked skills (they can leave HEAD detached, which stranded two commits on #98), fail fast on
  data access, calibrate before flagging, a 35-minute time box, one unit per cycle.
- **Owner action:** paste the updated prompt into the routine (`trig_01PNAmgcfqef8CvhPAY8ggbP`).

## 2026-10-08: W5 count deficits in Legacy Surveys DR10: 340.5 deg² null, first W5 limit (D-063, #94)
- Two sessions (a cloud run's 10° pilot, then a worktree worker on the full 20° × 10° regions; see the coordination
  entry below). DR10 Tractor galaxies (r < 23.5) counted per `nest4096` HEALPix pixel on the Data Lab server
  (~170 MB, not a catalogue download); predicted deficit profile from `exotic_sim.count_ratio` with the measured
  number-count slope (`model_prediction`): ~10 % of galaxies inside θ_E missing, ratio 0.02 at θ/θ_E = 0.1.
- Screen: matched filter at θ_E = 2–32′ over desA (RA 20–40°) and desB (RA 50–70°), Dec −30 to −20; cross-region
  null calibrated on vetted peaks (galaxy clustering makes Poisson errors 1.1–5.8× too small). 247,361 peaks →
  **40 flags → 0 survivors** (NGC 1398 sky over-subtraction ×3, bright stars, depth/tile edges; the rest consistent
  with the null). Contact sheet inspected by the worker and the coordinator. `derived`.
- 5,232 injections through screen + vetting: efficiency 0.50–0.73 for θ_E = 8–32′, 0 below 4′. **95 % limit on the
  sky density of W5 lenses n₉₅ ≈ 0.012–0.018 deg⁻² at θ_E = 8–32′** (≈ 3 × 10¹⁰–5 × 10¹¹ M☉ at 1 kpc, 3 × 10¹³–5 ×
  10¹⁴ M☉ at 1 Mpc; geometry ASSUMPTION). Untestable below θ_E ≈ 6′ with DR10 counts.
- **Failed approaches (rules):** Data Lab ADQL rejects sub-selects, CASE, SIGN and GROUP BY on expressions (group by
  `nest4096` only); RA/Dec chunks split edge pixels — sum them, or every chunk border looks like a deficit; an
  unbounded bright-end count slope makes the profile blow up at x → 0 (fix 0.6); a null from the other region's raw
  peaks inherits its artefacts (vet them first); veto radii growing with θ_E removed 90 % of random positions at 32′
  (veto only mimics that can empty ≥ 10 % of the core); HyperLEDA returns sexagesimal unless `_RAJ2000`/`_DEJ2000`
  are requested; a large-galaxy veto is needed (sky over-subtraction around NGC 1398); `pkill -f` killed the calling
  shell; JSON writers must end the file with a newline (pre-commit end-of-file-fixer failed CI twice).
- **Next:** deeper counts (HSC, Euclid) for θ_E < 6′; larger contiguous area for θ_E ≈ 1°; review leftovers: number
  counts divided by the full box area (small bias), per-job rebuild of region state in `inject`.

## 2026-10-08: W1/W2 in rejected lensed-quasar pairs: LS DR10 cannot decide them (D-064)
- Hypothesis: a dark deflector hides among pairs rejected as lenses for lack of a lens galaxy (Lemon 2023 UQP/QSO
  pair; SQLS "no lensing object", "QSO pair", "binary"). 123 rejected and 106 control lenses (≤ 3″) went through the
  D-056 chain.
- **Control efficiency 0 / 5** (1.9–2.6″; lens light blended into the images, fitted as PSFs): an LS "none" at these
  separations carries no information. Rejected: 24 "none". Of these, 10 are colour-mismatched, 2 are catalogued
  binaries and 1 has two redshifts; 11 are untestable. No limit, no candidate.
- **Failed approaches:** several input-format traps (sexagesimal Hennawi coordinates, SQLS two-row pairs, Lemon
  `z2` semantics, time-stamped VizieR headers) and sample-definition bugs were caught in review and fixed; the list is
  in D-064 "Evidence".
- **Next:** HST/Euclid/HSC image models or spectra for the 11 colour-matched pairs; the LaTeX-only NIQ tables.

## 2026-10-08: W3 MOA-II: calibrated baseline test, LF-drawn injections, streaming pipeline; gb22 re-run null (D-068)
- Cloud runs (session that opened #95, taken over at 23:12 UTC after 35 min idle to land it). The variable-baseline
  threshold is now the 95th percentile of the field's quiet χ²/dof (5.31 in gb22; ASSUMPTION), injections are drawn
  from the luminosity function, 200 per cell.
- gb22: same 30 flags as D-062 (event IDs identical; contact sheet inspected in D-062), **0 survivors**.
  **Γ₉₅ ≈ 1.1–3.3 × 10⁻⁶ per star per year for t_E = 10–300 d, 4–7 × 10⁻⁶ at 3 d** (103 / 2,000 recovered,
  0 / 400 PSPL controls; `derived`, N_s ASSUMPTION). Supersedes D-062's numbers.
- Streaming per-field pipeline (`moa_stream`, `w3_moa.py --field`): no tar on disk; gb22 by HTTP in 81 s;
  258 light curves/s on gb21 (test only, no gb21 result yet).
- **Failed approach:** a parent process reading ranges for its workers was OOM-killed on gb21.
- **Handoff / next:** `w3_moa.py --field gbN run-field` for the Nunota et al. 2024 Table 1 fields, smallest first
  (gb21 12.3 GB, gb20 15.2 GB, gb19 16.3 GB …; sizes in docs/exotic_limits.md), then `combine`.

## 2026-10-08: Coordination and dispatch rules for concurrent sessions (owner text)
- Two sessions worked PR #94 (W5) at once: a cloud run started from its 20:39 skeleton, another session pushed the full
  screen at 21:14, and the run's push was rejected (its pilot is parked on `claude/w5-clustered-null`, findings on #94).
  Commit-age heuristics cannot see a session that is coding but has not pushed.
- docs/cloud-routine-prompt.md gets the owner's "COORDINATION AND DISPATCH" section (dispatch first; `claimed` label +
  claim comment with a 10-minute heartbeat, stale after 20; re-check before every push, `-alt` branch on collision;
  "D-069" until merge). It replaces the 15-minute commit-age rule; docs/operations.md and its label table follow.
- **Owner action:** paste the new prompt into the routine (`trig_01PNAmgcfqef8CvhPAY8ggbP`); this file is only the copy.

## 2026-10-08: Governance — owner decisions on scope, parallelism, cloud disk and layout
- Owner direction (2026-10-08, given by the owner in the session that opened this PR) with a one-time
  authorisation to edit CLAUDE.md, including "Owner decisions", in this single PR and to merge it once CI is green
  and `/code-review` findings are fixed: four owner-decision bullets added verbatim (scope: any public dataset; smart
  parallelism and `[field: <unit>]` claim PRs; cloud disk streaming; layout). "Prefer single-thread work" removed
  from Budget. Later edits to "Owner decisions" again need the owner.
- Layout: root CLAUDE.md keeps owner decisions, merge policy, budget and pointers; module rules moved to
  `src/jwst_anomaly/CLAUDE.md` (stage contracts, signature layer), `scripts/CLAUDE.md` (screen / vet / inject /
  limit conventions, parallel topology) and `data/manifests/CLAUDE.md` (data, manifests, cloud disk). New
  directory-scoped skills `scripts/.claude/skills/w3-survey` and `w12-lenscats` carry the recipes and every
  failed-approach rule recorded for those searches; a `w5-counts` skill follows when W5 (#94) lands.
- Inventory before the change: one CLAUDE.md (root, 97 lines) and three root skills (`research-cycle`,
  `reuse-check`, `vet-candidate`), kept at the root; `research-cycle` now points to the new files and the claim rule.

## 2026-10-08: W3 MOA-II gb22 limit corrected: injections re-run after the `s_min` review fix (D-062)
- #92's limit (Γ₉₅ ≈ 0.5–5 × 10⁻⁶, 171 / 600 recovered) came from injections run before the review fix that made
  the pre-screen `s_min` the weaker of the two significances; the merged code already has the fix. Re-run of the whole
  chain on that code: **101 / 600** W3 injections recovered, 0 / 100 PSPL controls; **Γ₉₅ ≈ 1.2–4.4 × 10⁻⁶ per star per
  year for t_E = 10–300 d, 0.6–1.1 × 10⁻⁵ at 3 d** (`derived`; N_s ASSUMPTION unchanged). The real-data result is
  unchanged and reproduced with main's code: 30 flags, 0 survivors. Tracked fit table regenerated.
- Largest efficiency loss: the variable-baseline test (73 of 188 vetted injections; 35 % of quiet carriers have
  χ²/dof > 2 alone). **Rule:** re-run injections after any pre-screen change before quoting a limit.

## 2026-10-08: W3 in MOA-II Cut-0 light curves (gb22): 30 flags, 0 survivors, first W3 rate limit (D-062)
- Cloud runs (two sessions merged on `claude/w3-moa`). Hypothesis: the MOA-II 9-year release publishes every Cut-0
  object (positive *or negative* difference-image detections) before any bump cut, so unlike OGLE (0 / 600) and
  Gaia (0 / 240) its selection should keep W3 dimming events. `gb22.tar` (3.5 GB, smallest field, sha256 pinned).
- Pre-screen on all 18,599 light curves: 1,058 shape passes, **30** after the shared-epoch tests (1,022 removed at
  field-wide shared epochs, 6 at chip-wide ones); 0 / 100 PSPL controls pass. 30 flags (ΔBIC < −10), **0 survive
  vetting** (funnel and the three eclipse-like dips in docs/exotic_limits.md "W3 in MOA-II"). `derived`.
- Injection-recovery through Cut-0 emulation, pre-screen, fit and vetting (600 W3, 100 PSPL): recovered 8–53 % of
  injections per cell. **Γ₉₅ ≈ 0.5–5 × 10⁻⁶ per monitored star per year** for t_E = 3–300 d, strongest at ~30 d
  (N_s = 3.5 × 10⁶ is an ASSUMPTION from Nunota et al.'s N_s-per-object ratio; 1.43× weaker at its low end).
- **Failed approaches:** fitting all 1,058 shape passes (chunk 1/8 flagged 131 / 133: dip-shaped variables at
  shared epochs prefer negative-flux models; superseded by the shared-epoch cut, table kept as
  `fits_gb22_chunk1of8.ecsv.gz`). The first vetting chain left 9 survivors; the contact sheet showed one- or
  two-night drops and flat dips without caustic spikes, so the eclipse model, chip-level shared-epoch test and
  ≥ 3-night rule were added (their cost to W3 is in the injection numbers). Waiting with `pgrep -f`/`pkill -f`
  matches the waiting shell itself; wait on a log line or output file.
- Reproduced in a later cloud run from the tracked fit table (`prescreen`, `merge-chunks --n 1`, `vet`, `sheet`).
- Final `/code-review` fixes: Pool workers get the MOA fit bounds and shared-epoch populations through an
  initializer (spawned workers on Windows inherited neither, so the shared-epoch cut was silently skipped there);
  `fit --chunk` writes its own file instead of replacing the table `vet` reads; `merge-chunks` also checks the
  pre-screen `Params`; the sheet labels a PAR winner's PSPL curve as such.
- **Next:** a field with a published N_s (Nunota et al. 2024 Table 1) to make the limit model-light; the Cut-0
  spurious-detection filters are not emulated, so ε (and the limit) are optimistic by an unknown factor.

## 2026-10-08: W3 in the Gaia DR3 microlensing candidates: the published selection rejects W3 too (D-061)
- Cloud run. Hypothesis: the Gaia DR3 candidates (363; Wyrzykowski et al. 2023) are selected less PSPL-shaped than
  Mróz et al., so they could limit W3. Tested the selection before fitting, as D-057 requires.
- New `gaia_mulens.GaiaDR3Microlensing` adapter (TAP + DataLink epoch photometry, 8 parallel batches of 12 ids,
  ~3 s per source on the server; Table D.1 sample labels from the pinned arXiv source) and `scripts/w3_gaia.py`.
- Emulated Sample A cuts pass 143 / 163 real Sample A events. W3 injections: **2 / 240 selected, 0 / 240 selected
  and flagged**; PSPL controls 17 / 120. W3 dimming fails skewness < 0 and the skew–Abbe cut, as designed for
  brightenings. The fits of all 363 give one flag (4053892503992268288, ΔBIC −40.3). It is an event truncated at the
  window end on a variable baseline (light curve inspected), not a candidate.
- **Failed approaches (rules):** the Extractor cuts (n points, duration > 135 d, max σ > 50) cannot be emulated from
  the paper. The guessed definitions fail 126 / 163 real events, so they are left out. A single-id DataLink request
  answers with bare CSV, not a zip. Truncated chunked replies happen, so retry them.
- **Next:** W3 needs light curves taken before any microlensing selection. Gaia DR3 has epoch photometry only for
  its variable sources (vari_summary; ~11.7 M per the DR3 release, recheck), selected by variability, not shape. Next, check whether W3 survives that
  classifier (`vari_classifier_result`) by injection, and whether a sky-region subset is small enough to screen.
  Other options: KMTNet public seasons; OGLE EWS (owner decision, terms).

## 2026-10-08: W2 deflector test at HST resolution from the Hubble Source Catalog: not decisive (D-060)
- Cloud run. Hypothesis: HST resolution decides the lensed quasars that Legacy Surveys left blended or too close
  (D-056). `scripts/w12_hsc_probe.py`: HSC v3 summary sources within 4″ of each of the 444 galaxy-scale quasar/radio
  systems (MAST catalogs API, 58 s); ≥ 2 point sources (CI < 1.3) are the images, an extended source (CI ≥ 1.5)
  near their centroid and > 0.2″ from an image is the deflector; sources in < 2 HSC images are dropped as likely
  artifacts (MAST's recommendation; #90 review) (ASSUMPTIONs in `Params`).
- 71 of 444 systems have HSC sources in ≥ 2 HSC images (91 with any). Validation on systems with a published lens redshift (a lens galaxy is known
  to exist): 13 deflector, 15 none, 38 undecided → **efficiency 13/28 = 0.46** (`derived`; 15/35 = 0.43 without the
  artifact cut). Misses include quads and doubles (13 of 16 "none" have 2 point images). In four inspected misses (H1413+117, HE1104−1805, SBS0909+532,
  HE2149−2745) the HSC rows within 4″ are only the quasar images: the lens galaxy is absent from the catalogue, not
  mis-typed (likely lost in the quasar PSF; hypothesis, no cutouts inspected). Without a lens redshift: 1 none
  (HS0810+2554), 4 undecided. The 3 open SuGOHI IX pairs (D-056) have no HSC sources.
- **Failed approach (rule):** HST *catalogue* photometry cannot decide a dark deflector in lensed quasars — a
  "none" is more likely a missed lens (0.54) than a dark one. No limit, no candidate; HS0810+2554 is not flagged.
- **Next:** PSF-subtracted HST image modelling (e.g. drizzled frames from MAST, quasar PSF + Sérsic fit) validated
  on the same known-lens set, or HSC PDR3 photometry; until then W2 in wide imaging stays at D-056.

## 2026-10-08: Warp literature check: still nothing testable
- Subagent search (arXiv API 2023–2026, INSPIRE citations of Clough et al. 2024): no imaging or lensing prediction
  for a distant observer; Lentz & Felton 2024 give order-of-magnitude EM fluxes for a bubble 100 lyr away but no
  template that separates it from ordinary transients (found by review); the collapse-burst waveform is not public; an O3 search for superluminal-source GW
  bursts (Kuwahara & Cannon 2023) is already null. D-047 stands; recorded in SOURCES.md and docs/exotic_lensing.md.

## 2026-10-08: W3 OGLE bulge (all 5,790 events): no candidate; the published selection rejects every W3 event (D-057)
- Worktree worker, unbounded π_E (predates D-058; bounding can only add flags). Bulge: 0 fit failures; best ordinary
  PSPL 5,377 / PAR 401 / FSPL 12; ΔBIC(best exotic) 5/50/95 % = −3.6 / 5.7 / 12.0; 127 flags < −10 (`derived`).
- Vetting, cumulative: 127 → 113 refit all ordinary → 80 robust errors → 73 variable baseline → 14 season offsets →
  9 season drifts → 7 binary source/lens, VSX/Gaia, arXiv → 1 `feature_coverage` → 1 jackknife → **0** two unrelated
  events (new `revet` tests: ≥ 3 epochs where the models differ by > 3σ; drop up to 3 influential epochs keeping ≥ 3 in
  the feature; two independent PSPL bumps). BLG519.21.110304's exotic spikes sat on the 2011 event and a 1-day bump in
  2015 (two PSPL bumps better by ΔBIC 24.4). Disk: 6 → 0.
  Contact sheet of the 7 late survivors inspected. The D-059 chunk survivors fail `feature_coverage` here.
- Injections: 600 W3 events (n = 1, ε < 0; t_E 3–300 d, ρ 0.01/0.1) on real bulge cadences + 300 PSPL controls. The
  fitter flags 42–97 %, but **0 / 600 pass the emulated Mróz selection** (controls 15–43 %; cuts failed most: one
  bump, PSPL χ², χ₃₊). The emulation passes 63.9 % of the real selected events (somewhat stricter). No rate limit is
  derivable; ε_W3/ε_PSPL < 0.12–0.33 (95 %). Wall time ≈ 5 CPU h (bulge fit 10.4 ks, vet 6.6 ks, inject 2.0 ks).
- **Failed approaches (rules):** exotic fits started only from the PSPL solution miss the W3 geometry — start with
  caustic spikes on pairs of maxima and on absolute t_E; ΔBIC alone is not a candidate test (season blends and
  feature coverage remove 120 of 127); a jackknife must keep ≥ 3 epochs in the exotic feature or it kills real short
  events; an exotic fit whose spikes sit on two bumps years apart needs the two-unrelated-events test; pin BLAS threads (`OMP_NUM_THREADS=1`) with process pools; a published
  PSPL-selected sample cannot constrain a non-PSPL signal — inject through its selection before fitting it.

## 2026-10-08: W3 OGLE bulge, chunks 2–3/12: 20 flags, no candidate; chunk tables tracked (D-059)
- Cloud run. Chunk 1/12's fit table lived only in the ephemeral session and is lost. `fit --chunk K/N` now also
  writes a deterministic gzipped copy to `results/w3_ogle/` (~0.22 MB per chunk); `merge-chunks --n 12` joins
  chunks 1..12 into the table `vet` reads and marks it the whole sample only when every chunk is present, fitted
  with the current `Params` and holds exactly its own events (none skipped) (D-059). A chunk takes ~9 min on 4 cores (not ~17).
- Chunk 2/12 (483 events, 0 errors): best ordinary PSPL 454 / PAR 28 / FSPL 1; ΔBIC(min exotic) 5/25/50/75/95 % =
  −4.7 / 3.8 / 5.8 / 6.6 / 11.2; 49 below 0, 9 flags below −10 (`derived`). `vet`: 8 fail; **BLG624.18.69573**
  (no EWS name; t_E ≈ 180–240 d, best ordinary PAR) passes every automated test (N1neg ΔBIC −19.0 vs PAR; season
  offsets −15.0, drifts −15.5; binary source / lens −10.1; 0 VSX / Gaia matches). Contact sheet and residuals
  inspected: the N1neg model puts its first caustic spike inside a season gap (no data at t − t0 ∈ [−50, 0] d; the
  models differ by > 5 % over [−89, −7] d); its Δχ² comes from 2 peak points (−9.4) and 16 post-peak points (−6.7).
  An unsampled caustic plus a sparsely sampled peak is not evidence: not a candidate (ASSUMPTION-level judgement).
- Chunk 3/12 (483 events, 0 errors): PSPL 447 / PAR 36 / FSPL 0; ΔBIC(min) 5/50/95 % = −2.6 / 5.7 / 11.1; 39 below
  0, 11 flags; **0 survive** `vet` (season offsets/drifts remove 10, the refit of all ordinary models the 11th; one
  VSX match). Contact sheet inspected; in several flags the exotic and ordinary curves also differ mainly in gaps
  (e.g. BLG597.28.9837 has no peak data).
- **Failed approach:** chaining chunk runs with `while pgrep -f 'chunk 2/12'` — the waiting shell matches its own
  pattern and never starts the next chunk.
- **Next:** a `gap_coverage` vetting test (require data where the exotic and best ordinary models differ, else the
  flag fails); refit chunk 1 and fit chunks 4–12 (two or three per run); `merge-chunks`, `vet`, `sheet`; the
  empirical ΔBIC null and xallarap fit for BLG667.04.62161 (scratch null-simulation design: PSPL best fit plus white
  noise, and plus season-wise circularly shifted residuals); then `inject` / `limit`.

## 2026-10-08: W3 OGLE bulge, chunk 1/12 (483 events): one marginal flag survives automated vetting
- Cloud run. `w3_microlensing.py fit --chunk K/N` fits events K−1, K−1+N, … so sessions fit disjoint, field-balanced
  parts of the bulge sample (measured ~10 s CPU per bulge event, ~4 h for all 5,790 on 4 cores; one chunk of 12 is
  ~17 min wall). Chunk 1/12: 483 events, 0 errors; best ordinary PSPL 449 / PAR 32 / FSPL 2; ΔBIC(min exotic)
  5/25/50/75/95 % = −3.5 / 3.7 / 5.6 / 6.5 / 12.4; 41 below 0, 12 flags below −10 (`derived`).
- `vet`: 11 of 12 fail; **BLG667.04.62161 (OGLE-2015-BLG-1250) passes every automated test**: E2neg ΔBIC −12.1 vs
  PSPL; season offsets −12.7, season drifts −12.4, binary source / binary lens no better than PSPL, 0 outliers,
  baseline χ²/dof 0.87, 0 arXiv records, no VSX / Gaia variable within 1″. Contact sheet inspected: faint source
  (peak flux ~2× baseline), large scatter, a flattened peak and a few low points ~+20…+50 d; the E2neg plateau fits
  those. ASSUMPTION-level judgement: marginal, not a candidate — with 483 trials and 12 flags below −10 from a
  heavy-tailed ΔBIC distribution, one −12 survivor is expected without any exotic lens; untested ordinary
  explanations: xallarap, per-season error underestimation near the peak, blending/difference-imaging systematics
  of a faint source, and a calibrated null (the injection/limit stage).
- **Next:** chunks 2–12 (one per run: `fit --chunk K/12`, then `vet` / `sheet`); for BLG667.04.62161, an
  empirical ΔBIC null from the same chunk's PSPL-simulated light curves (does −12 occur at rate ≥ 1/483?) and an
  xallarap fit before any further attention; then `inject` / `limit`. `limit` refuses vetting of a single chunk
  (#86 Codex), so a `merge-chunks` step (concatenate `fits_bulge2019_chunk*of12.ecsv`, then `vet`) comes first.

## 2026-10-08: W3 OGLE disk re-fit with bounded parallax: still null (D-058)
- Cloud run. ASSUMPTION `Params.pie_max = 5`: PAR fits with |π_E| > 5 are rejected in every fit and vetting refit
  (the unbounded fits reached π_E ~ 30–1,400). Disk (460 events, 376 s): best ordinary PSPL 408 / PAR 52 / FSPL 0;
  36 of 368 PAR fits sit on the bound; ΔBIC(min exotic) 5/50/95 % = −4.0 / 3.6 / 6.9.
- 7 flags (was 6; new GD1217.10.8703 at ΔBIC −10.0, `E2pos`), **0 survive** `vet` (all tests complete). The new flag
  is four post-peak points 0.15–0.35 below baseline that neither model fits; per-season offsets remove it (ΔBIC 0.9).
  Contact sheet inspected.
- Absorption check (simulated: 40 E2pos/N1neg injections, 700 uniform epochs, white noise σ = 0.05–0.3): flags 25
  bounded vs 25 unbounded; ΔBIC shifts ≤ 0.46. On this cadence unbounded parallax does not absorb exotic signals,
  so the bound is a physical prior, not a sensitivity gain. Not tested: seasonal gaps and correlated systematics.
- Environment: install the `mulens` extra (`-e ".[dev,cloud,mulens]"`); without MulensModel `fit` silently skips PAR.
- Follow-up (#84 Codex): checkpoint rows carry a `Params` hash; `fit` refits rows from other Params, so a chunked
  bulge run never mixes bounded and unbounded parallax fits.
- **Next:** bulge `fit` in chunks or locally, then `vet` / `sheet`; `inject` / `limit` with season-drift vetting.

## 2026-10-08: W1/W2 in published lens catalogues: no dark deflector; weak limits (D-056)
- Worktree worker: lenscat (32,838), Euclid Q1 Discovery Engine (2,584) and SuGOHI (3,961) merged into 35,862 systems
  (`src/jwst_anomaly/lenscats.py`, `scripts/w12_lenscats.py`; pinned by sha256). Deflector test in Legacy Surveys
  DR10 Tractor; footprint and depth from the DR10 brick summary, not from detected sources.
- 17,555 galaxy-scale systems in the footprint. Galaxy-finder (17,102) and sub-mm (110) systems cannot show a dark lens
  and give no limit; only lensed-quasar (325) and radio-interferometric (18) systems are tested (pair test: two point
  images, nothing bright enough near the expected deflector). Blended, too-close and faint-galaxy cases are undecided.
- 29 decided (13 with a deflector, 16 without); 13 of the 16 have a literature lens galaxy, and 3 SuGOHI IX CHITAH
  pairs (090434−005328, 091517+040747, 104122−005618) are open only in the typical variant — a lens below the local
  LS depth explains them conservatively. Not candidates. Cutout sheets inspected.
- 95 % limits on the dark-deflector fraction, test completeness assumed (not measured): typical f_dark < 0.48
  (quasar, k = 3, N = 16), < 0.23 (radio, k = 0, N = 13), < 0.27 (all, k = 3, N = 29); conservative < 0.50 (quasar,
  N = 6). Earlier 1.5 × 10⁻⁴ / 8.2 × 10⁻³ / 0.13 are withdrawn (PR #81 reviews). No W1 geometry (no image positions).
- **Failed approaches (rules):** a limit is valid only over systems where the test could have found the signal —
  prove it by injection; a recovery factor from deleting deflectors and re-running the same code is 1 by
  construction; coverage must come from footprint/depth products, never from "a source nearby"; evaluate exclusion
  flags on every system; Data Lab TAP takes no table uploads or q3c (batch box ORs); lenscat types cluster-survey
  entries as "galaxy", has AGEL declination and SPT position errors and rounded positions, and keeps rejected
  candidates; "no lens redshift" ≠ "no lens".

## 2026-10-08: W3 OGLE-IV disk sample: null (D-057)
- Cloud run. `jwst_anomaly.ogle` (Mróz et al. 2019/2020 adapter) and `scripts/w3_microlensing.py`
  (`fit` / `vet` / `sheet` / `inject` / `audit` / `limit` / `manifest` / `summary`). Manifest
  `data/manifests/ogle_mroz.ecsv`. `www.astrouw.edu.pl` reachable from the cloud (2026-10-08).
- **Disk (Mróz 2020, all 460 Table B1 events; 525 s on 4 cores):** best ordinary PSPL 406 / PAR 54 / FSPL 0.
  ΔBIC(min exotic) quantiles 5/50/95 % = −3.9 / 3.6 / 6.9; 6 flags below −10, all `E2pos`. `vet`: 0 survivors;
  every flag loses the exotic preference under per-season baseline offsets and/or drifts (ΔBIC −4.6 … +8.0), and
  BLG568.12.9169 is already −5.1 on refit. arXiv 0 mentions, no VSX / Gaia DR3 variable within 1″. Contact sheet
  inspected: two flags have sparse peak coverage; two have post-peak points below baseline (the Ellis
  demagnification shape) that one season's drift absorbs.
- **Limits:** none yet (injections use bulge light curves); the disk null is a flag count, not a rate limit.
- Review follow-up (#82 Codex): binary-lens α starts were passed in radians to MulensModel (degrees), so only
  0.5–5.8° was searched; fixed and disk re-vetted: still 0 survivors (BL BICs move by ≤ 18). Injection vetting now
  refits PAR with season trends as `vet` does; failed XMatch queries leave a flag unvetted; `limit` refuses a
  zero-event limit unless the bulge vetting is a complete null.
- **Caveats:** unbounded parallax fits reach π_E ~ 30–1,400 (unphysical); they can absorb an exotic signal and cut
  sensitivity. The season-drift test may also absorb real W3 dips: calibrate both with injections.
- Timing: the bulge fit is ~1.1 s/event on 4 cores (~1.8 h for 5,790); `fit` checkpoints to
  `fits_<key>.partial.jsonl`, but cloud disks are ephemeral: fit in chunks per run or locally.
- **Next:** bound π_E (ASSUMPTION, e.g. |π_E| < 5) and re-fit disk; bulge `fit` (local or chunked), `vet`, `sheet`;
  `inject` / `limit` with season-drift vetting inside the injection loop.

## 2026-10-08: W1 shear (aperture-mass) screen: four clusters null; limits 5–8× stronger than radial (D-053)
- Cloud run. `exotic_screens.py shear` builds D-050: PSF-deconvolved catalogue ε, cluster shear removed, spike
  segments vetoed, Schirmer 10″ aperture-mass S/N map, rotation null, B-mode check. `scripts/inject_shear.py` reuses
  D-049's W1 painting.
- **Validation on known signals.** Measured ε along the cluster model's g rises with |g| (MACS0416, Abell 2744), with
  responsivity R = 0.41–0.48, not 1. The screen now removes R g; injected images keep R of their lens shear.
- **Real fields** (Abell 2744, MACS0416, MACS1149, Abell 370; 30.2 arcmin²): **null**. The first run (no spike veto)
  had Abell 370 E 4.51 / B 4.68 and Abell 2744 E 3.92, all p ≤ 0.01 against the rotation null. Diffraction spikes
  were the cause: the veto (77–178 segments per field) brought them to 3.48 / 3.46 and 3.39. B-mode extremes still
  beat the rotation null in three fields, so an E peak must also beat the field's max |S_×| (post hoc,
  conservative). None passes both.
- Injections (200 per field and mass; spike rows never painted): 33 / 197 / 405 of 800 at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (radial 8 / 84 /
  156 of 1,600); 1 at 2 × 10¹¹, 0 at 2 × 10¹⁰. 95 % limits 7.7 × 10³ / 1.3 × 10³ / 6.5 × 10² deg⁻² (radial 6.1 × 10⁴ /
  7.0 × 10³ / 4.0 × 10³). docs/exotic_limits.md.
- **Failed approaches (rules):** subtracting the full model g from catalogue moments (leaves a W1-signed radial
  residual); any radial-alignment statistic without the spike veto; trusting the rotation null without a B-mode
  check; point-mass 1/x² and top-hat filters (lost to Schirmer 10″ in injections).
- Data: the DJA SMACS v7.4 and El Gordo v7.0 photo-z tarballs return 404 (2026-10-08), so these fields are out. MAST
  catalogues and CANUCS DR1 catalogues (~30 MB each) downloaded fine.
- Wall time: 58–168 s per field for 1,000 injections; the base screen takes 6–13 s. Review rounds changed the injection model (spike and
  near-core and unresolved rows not painted, lens change from raw moments); every re-run kept the four real fields
  null. The counts moved from 59/264/459 (first model) to 33/197/405 (final); not painting unresolved rows was the
  largest step (from 50/238/449; docs/exotic_limits.md caveats).
- Final review: painted images now face the spike veto too (409 → 405 at 2 × 10¹³ M☉; limits unchanged at two
  digits). `psf_sigma_px` (now shared with `inject_radial.py`) raises when no S/N > 50 row has a size, where it
  returned NaN. Left as is (maintainability only): the E/B summary is computed in both `cmd_shear` and
  `inject_shear.run_field`; the null recomputes the |e|² noise term per draw.
- **Handoff:** TASKS "JWST focus" 1 W1: lower the B-mode floor (PSF anisotropy, blends), more fields, stacking.
## 2026-10-08: Orphan-pair null (e) fixed, companion-aware null (f): no deep-field excess (D-055)
- Cloud run. Hypothesis: the D-051 flanking-field orphan excess (246 vs 211, P = 0.010) comes from null (e)'s
  colour cell (member i only, non-finite colours in the 0–0.3 bin) or from physical companions the 10–30″
  reference misses.
- `pair_cells` uses both members' colour bins (unordered), a non-finite bin and a ≥ 5-pair cell floor (S/N × size
  fallback); null (f) conditions on 3–6″ pairs. Nine fields re-run (≈ 1 min each, in parallel), orphans unchanged.
  Deep fields: (e) P = 0.13, (f) P = 0.35; clusters P ≥ 0.13. **The excess was a null-model artefact.** Table:
  docs/orphan_pairs.md.
- Injections re-run (6 fields in parallel, ~4 min): no-candidate limits reproduce D-051; background-aware limits
  tighten by ×0.76 (docs/exotic_limits.md). Orphan set unchanged, so D-051's inspected contact sheets stand.
- Not separated: the shares of the cell changes.
- **Handoff:** W2 next is the segmentation-map same_galaxy rule (TASKS "Now" 1); PR #78 (W1 shear) in flight.
## 2026-10-08: Survey-agnostic signature layer; any public dataset in scope (D-054)
- Owner direction: find evidence of traversable wormholes / negative-mass objects or warp-drive spacetimes in any
  public dataset. `jwst_anomaly.signatures` registers W1/W2/W3/W5 (prediction, injection, screens, ordinary
  mimics, limits) and defines the `LightCurveSurvey` / `CatalogueSurvey` adapters.
- Reuse-check: MulensModel (extra `mulens`) for ordinary microlensing fits; OGLE-IV Mróz et al. 2019/2020 samples
  first for W3 (with published efficiencies), then Gaia DR3 `vari_microlensing`, then KMTNet; lens catalogues
  (lenscat, Euclid Q1, SuGOHI) for W1/W2. No published survey light-curve search for negative-mass or Ellis
  events exists, so a W3 limit would be new (owner first). OGLE EWS use waits for the owner (its terms).
- **Handoff:** W3 re-analysis of the Mróz OGLE-IV samples (TASKS "Now" 1).

## 2026-10-08: W3 multi-epoch dimming / inverted-microlensing screen: null, limits (D-052)
- Worktree worker: `scripts/dimming_screen.py` cross-matches per-epoch level-3 catalogues (F200W + F444W), flags
  vanishing, achromatic-dimming and rise-dip-rise sources, vetoes catalogue effects (incl. a `bright_neighbour`
  veto: a source ≥ 100× brighter within 1.5″), and confirms with S3 byte-range forced photometry. Fields:
  NEXUS-Center (8 epochs), MACS0416 (5; F200W baseline 0.35 yr, F444W 3.26 yr), Abell 2744 (3). El Gordo rejected
  (one shared band); JADES deferred.
- 5,177 catalogue flags → 1,441 after catalogue tests → 6 confirmed by forced photometry (MACS0416) → 0 after the
  bright-neighbour veto (all six sit on or beside saturated stars; cutout sheet inspected: the star's wings and
  spikes rotate through the aperture with each epoch's mosaic orientation). No candidate.
- Injection (`exotic_sim.inject_light_curve`, n = 1, ε < 0) with per-copy vetoes: efficiency 0.04–0.24. Headline
  95 % limits (MACS0416 only, the one field with calibrated forced errors, 34 sources): rate < 0.056 per source
  per yr at t_E = 1 yr (umbra fraction τ < 0.18); all fields indicative < 0.015 (τ < 0.049). The per-copy
  single-epoch veto makes full vanishes of two-epoch sources undetectable (main efficiency loss).
- **Failed approaches (rules):** a catalogue non-detection is zero flux only if the source would have been ≥ 10σ;
  vanish only when every testable band vanishes; `is_extended == False` is not a point-source cut (use CI_70_30);
  inject multiplicatively with matched noise; Gaia masks must drop the source's own match before looking for a
  neighbour; forced fluxes need each epoch's pixel solid angle; never credit exposure to sources flagged before
  injection; don't run parallel or multi-target S3 cutout jobs through the proxy (s3fs "bucket does not exist").
- Caveats: forced counts are from the 06:39–07:03 UTC run (later S3 failures), re-calibrated offline; NEXUS and
  Abell 2744 forced errors uncalibrated; PEARLS (MAST lists 2.0.1, headers 3.0.0) + CANUCS processing differ.
- Wall time: screen 4–76 s, injections 20–136 s, forced photometry 3–32 min per field.
- Final /code-review fixes (cloud run): the forced-stage SIMBAD/NED label read non-existent `*_otype` columns
  (always empty; now `best_match_*` of `crossmatch.XMATCH_COLUMNS`); `calibrated` is False when an epoch image was
  unreadable; the dead saturated-star branch of `bright_neighbour` removed (`near_star` covers it). Counts and
  limits unchanged (no survivor reached the cross-match). A failed SIMBAD/NED service is now named in
  `forced.ecsv` meta (`n_<service>` == -1), `inject` also treats `unread_images` as uncalibrated, and
  `check_params` guards the saturated-star / near_star coupling (#77).

## 2026-10-08: Orphan-pair cutout footprints: frame-token polygons, deterministic visit order
- Cloud run. Review follow-ups to #73 that its final squash did not carry: `_in_region` parses `POLYGON ICRS …`
  (a frame token used to match nothing, so every cutout read "outside"); MAST observations are sorted by obs_id
  before "first covering visit" is chosen (archive row order no longer picks the cutout). Counts and limits do
  not depend on either (cutouts only). TASKS: null (e) colour cell asymmetry.
- **Handoff:** TASKS "Now" 1, the D-050 W1 shear screen, is next; it needs the SMACS `_cat.ecsv` and DJA zout.

## 2026-10-08: Dark-deflector (orphan-pair) screen in deep fields, with injection-recovery: null (D-051)
- Worktree worker: `orphan_pairs.py` now reads CANUCS and DJA catalogues through one column layout (D-048 counts
  reproduced bit for bit) and runs on six deep fields: the five CANUCS NIRCam flanking fields and DJA GOODS-N
  (108 arcmin²). `scripts/inject_pairs.py` paints `exotic_sim` point-mass, W2 Ellis and W1 pairs into the real
  catalogues and runs the unchanged screen.
- 355 orphans against 315.4 expected (strictest null, P = 0.015); the excess sits in the flanking fields (246 vs
  211, P = 0.010), GOODS-N matches every null. Read as unmodelled physical companions (lensing would need ~140 dark
  galaxy-mass deflectors per arcmin²); open, not a candidate. 90 top orphans inspected: 50 ordinary (knots,
  satellites, groups, artefacts), 40 faint chance-like pairs. Nothing for `/vet-candidate`.
- Per-deflector efficiency ≤ 0.8 % (point mass, W2) and 0.1–18.6 % (W1, θ_E 0.3–1.5″). Background-aware 95 % limits:
  W1 < 1.3 × 10⁴ deg⁻² at θ_E 1.5″ (4.5 × 10¹¹ M☉); W2 < 2.9 × 10⁵ deg⁻² at 0.7″ (throat ≈ 10 pc). This covers
  θ_E ≤ 1.5″, where `radial` (D-049) is blind. docs/exotic_limits.md, docs/orphan_pairs.md.
- **Failed approaches (rules):** DJA Kron apertures for the same-galaxy rule (2.6× CANUCS; use 3.3 × flux_radius);
  the standalone DJA GOODS-N zout (older catalogue; stream the tarball member); a global null (c) in deep fields
  (under-predicts orphans; use the conditioned null (e)); unretried S3 reads (spurious NoSuchBucket via the proxy).
- Downloads > 200 MB (GOODS-N catalogue 224 MB, photo-z tarball 371 MB streamed) stated in D-051.
- Wall time: 3–20 s per field search; 190–266 s per field of injections.

## 2026-10-08: W1 limits re-run: mass-parametrised, blend-aware, independent nulls (D-049 update)
- #70 merged an intermediate version. This re-run fixes the review findings: one lens mass with θ_E per source
  redshift (photo-z, else z_s = 2), overlapping image pairs painted as one blend, an independent 200-draw null per
  batch of 10 trials, screen grid/null/defaults shared with `cmd_radial` (byte-identical SMACS output), a stated
  S/N ≥ 5 detection floor, measured footprint-border excess (3–16 %).
- Recovered of 1,600 per mass (2 × 10¹⁰ / 2 × 10¹¹ / 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉): 0 / 0 / 8 / 84 / 156.
  Headline 95 % limits (six photo-z fields, 38.0 arcmin²): none below 10¹² M☉; < 6.1 × 10⁴, 7.0 × 10³, 4.0 × 10³
  deg⁻² at 2 × 10¹², 8 × 10¹², 2 × 10¹³ M☉ (all eight fields, optimistic: 3.7 × 10⁴, 3.8 × 10³, 2.1 × 10³).
- **Failed approaches (rules):** one θ_E for every source; two lines for overlapping images; one fixed null for all
  trials; a `nanmin(S/N)` detection floor; recovery at the 3-line peak instead of p_random; catalogue shapes
  without PSF deconvolution; headline limits including fields without photo-z.
- Wall time: 371–993 s per field (1,000 injections, 100 independent nulls; 4 parallel).
- **Handoff:** a W1-specific screen benchmarked with this harness (TASKS).

## 2026-10-08: W1 shear screen reuse-check (D-050)
- `/reuse-check` for the W1-specific screen: build a catalogue aperture-mass map (−M_ap: a negative-mass lens gives
  negative tangential shear) with scipy cKDTree. TreeCorr NG only gives a stacked ⟨N M_ap⟩(R) and has no Windows
  wheels; lenspack works on pixelised maps. Rough S/N 1–2 at θ_E = 1″ in a 5″ aperture, so gains are expected
  mainly at θ_E ≥ 2–3″.
- **Handoff:** implement D-050 in `exotic_screens.py` and benchmark it with `scripts/inject_radial.py` on SMACS
  first (smallest field; needs its MAST `_cat.ecsv` and DJA zout in the cache).

## 2026-10-08: PR #70 merged after review fixes; Abell 2744 radial doc brought to the post-D-034 result
- Cloud run. PR #70 (D-049) was conflicted with main: merged `origin/main`, then `/code-review` on the final diff.
  Main finding: MACS0717 and Abell S1063 have no photo-z, so the injector painted their members and foreground
  galaxies as W1 images (Abell S1063 had the highest efficiency of all fields). The headline limits now use the six
  photo-z fields (38.1 arcmin²): < 7.0 × 10⁴ / 6.8 × 10³ / 3.2 × 10³ deg⁻² at θ_E = 3″ / 6″ / 10″, about 27× weaker
  than Takahashi & Asada; the all-field values stay in `limits.json` as optimistic. `inject_radial` now reads the
  screen defaults from `exotic_screens.radial_defaults()`.
- docs/fields/abell2744.md: radial result with the D-034 spike veto (134 `anti` arcs, max 5 lines, p 0.505; null).
- PR #61 (n0153, `needs-human`) brought up to date with main; it still waits for the owner.
- **Handoff:** the W1-specific screen (TASKS "Now" 1) is next; it fits one full cycle.

## 2026-10-08: W1 injection-recovery through `radial`: the screen is blind to negative-mass lenses (D-049)
- Worktree worker: `scripts/inject_radial.py` paints `exotic_sim` W1 lenses (n = 1, ε < 0) into the real catalogues
  of all eight null `radial` fields and runs the unchanged screen (`exotic_screens.radial_candidates`, a pure
  refactor of `cmd_radial`'s selection). 200 lenses per field and θ_E; 51.2 arcmin² screened.
- Recovered: 0 / 1,600 at θ_E = 0.3″ and 1″; 10 at 3″; 80 at 6″; 181 at 10″. 95 % limits on W1 lens surface
  density (six photo-z fields): none below 3″; < 7.0 × 10⁴ deg⁻² at 3″ (|M| ≈ 1.4–4.3 × 10¹² M☉), < 3.2 × 10³ deg⁻²
  at 10″. MACS0717 and Abell S1063 (no photo-z, members get lensed) are excluded as optimistic. About 27×
  weaker than Takahashi & Asada spread over 0 < z < 1. docs/exotic_limits.md.
- Why: an image reaches e ≥ 0.5 only for β ≲ 2.3 θ_E, and `anti` against the cluster keeps a third, so a lens puts
  1–3 arcs into a screen whose null needs 5–8 lines.
- **Failed approaches (rules):** recovery is the screen's p_random, not the 3-line peak; cache the null draws once
  per field (0.2–1 s per lens instead of 15–150 s); deconvolve the PSF before applying the lens Jacobian.
- Wall time: 322–877 s per field for 1,000 lenses (5 θ_E).
- **Handoff:** a W1-specific screen (collinear radial pairs flanking an empty centre, orientation against the
  candidate centre, local null), benchmarked with this harness.

## 2026-10-08: El Gordo lens model reproduces the published magnification maps (issue #68)
- `validate` now compares Lenstool-par models with published signed μ maps (`mag_map_files`, fetched only by
  `validate`); `map_check` reports `parity_agree` for signed maps. El Gordo (CDS J/A+A/678/A3, z=2 and z=4, |μ|<10,
  37-px sub-grid): median |μ| ratio 1.00001 / 1.00001, p95 relative difference 0.04 % / 0.07 %, parity 100 %;
  χ²pos unchanged (82.53). At full resolution (250,000 points), >20 % mismatches sit only on critical curves or
  within ~1″ of member cores (0.4″ map pixels).
- El Gordo `bayes.dat` `Chi2` explained (D-045 amendment): the sampling run uses `forme -10`, an image-plane χ² with
  σ² = a·b from the image list (19 of 56 images at 1.24″). It reproduces `Chi2` row by row to 0.1 % and the
  ln(Lhood) offset (75.904) exactly. The chain is validated. **Failed approach:** a source-plane χ² with free
  per-family weights (held-out ρ 0.81) fit only partly; the definition came from Lenstool's source.
- **Handoff:** El Gordo μ, parities and the MCMC chain are validated. Next: posterior μ spreads in `fluxratios`.

## 2026-10-08: Exotic-lens predictions first: wormhole / negative-mass searchable, warp not (D-047)
- Owner focus (2026-10-08): search only for signatures of traversable wormholes / negative-mass lenses and warp-drive
  spacetimes, predictions first. Research worker; every citation fetched from arXiv / Crossref.
- Searchable: W1 negative-mass dark lens (radial pair beside an empty centre), W2 Ellis pair without deflector, W3
  inverted microlensing (umbra between caustic spikes), W5 count deficit. Not searchable: the 4 % Ellis gutter, µas
  shifts, and every warp signature (Alcubierre exterior is flat; no published imaging/photometric prediction for a
  distant observer). The warp branch is stopped.
- `jwst_anomaly.exotic_sim`: Kitamura+2013 power-law lens family (either sign of ε), finite-source light curves,
  `inject_images` / `inject_light_curve` (`simulated`) for injection-recovery.
- Review fixes: an exact finite-source integral, checked against inverse ray shooting, puts the negative-mass spike
  peaks at ×7.0 / 2.35 / 1.53 (ρ = 0.01 / 0.1 / 0.3; first version 9.2 / 2.6 / 1.7). Exact demagnification onset
  added (2/(n+1) is KNA13's large-n estimate); KNA13's n = 3 "~10 %" is rounding of their Fig. 2c (13–14 %). NaN
  epochs no longer read as an umbra.
- **Handoff:** injection-recovery for `radial` (W1) and the dark-lens search (W2) to turn nulls into limits; a
  dimming class for the transient screen (W3); counts around `radial` centres (W5).

## 2026-10-08: Multi-plane lens models (D-046)
- `LensModel.split_planes` moves named potentials (e.g. a foreground galaxy fitted as a member) to their own
  redshift; `MultiPlaneLensModel` solves the multi-plane lens equation, and `find_images` / `backtrace_images` /
  `imageplane_residuals` accept it.
- Reproduces the D-042 system-51 result with library code: σ ≤ 70 km/s at z 0.268 leaves exactly 3 images.
- **Handoff:** use `split_planes` in `/vet-candidate` for any extra or missing image near a non-member galaxy.

## 2026-10-08: Orphan image pairs, a blind dark-deflector screen: null (D-048)
- Worktree worker: `scripts/orphan_pairs.py` looks for SED-matched close pairs with no published system and no visible
  galaxy between them in the CANUCS DR1 catalogues of MACS0416, MACS1149 and Abell 370.
- Pair excesses come from same-redshift groups (null (c)); 38 orphans against 34.5 expected. The 35 top orphans are
  knots, group members or chance matches, all at |μ| ≈ 1–2.5. No candidate.
- **Failed approaches:** a midpoint-only lens rule (missed galaxies between members); no per-band S/N cut (25 % false
  SED matches).
- Wall time: 5–10 s per field per search, 30–60 s with cutouts.
- **Handoff:** injection-recovery so the null becomes a limit; segmentation-map same-galaxy test; more clusters.

## 2026-10-08: Lenstool MCMC posteriors; Abell 2744 multiplicity residual explained (D-045)
- `lensmodel.read_lenstool_bayes` / `posterior_par` and `lens_consistency.py posterior` rebuild published models at
  MCMC samples (potfile rescaling, sampled family redshifts). Validated on Abell 2744: best.par is a chain row, and
  two random rows give χ²pos 173.5 / 173.0 against the chain's 178.1 / 174.1.
- Abell 2744 (12 samples plus best.par): 3.2a/b, 34.1a/b and 700.1a/b stay merged in every model. CATS v4.1 splits
  34.1a/b, and 3.2a/b sit on the caustic in both models. That is model resolution at folds; 0 surviving. 4.2c moves
  only 0.2–0.4″ (μ 8.7–10.0), well below its galaxy-scale systematics (D-036).
- **Failed approach:** the El Gordo CDS chain's `Chi2` column (54–77) does not track our χ²pos (93–106), not even in
  rank, so that chain is not validated.
- Wall time: about 70 s per Abell 2744 sample on a 0.25″ grid; 15 min for 13 models.
- **Handoff:** feed the posterior μ spread into `fluxratios` (then recheck SMACS 6.3 once a SMACS chain is pinned) and
  the position spread into `forced_check` search radii.

## 2026-10-08: CANUCS Lenstool models pinned (D-044)
- `macs0416-canucs` reproduces Lenstool's image-plane χ²pos (330.8 against 344.30; rms 0.51″ over 303 images).
  `abell370-canucs` is a source-plane fit (image-plane rms 2.3″).
- They are the independent second model for vetting (the D-042 system-27 and system-51 checks used scratch code).
- `fluxratios` now refuses gated image lists and map models, as `images` does; its `--offset-arcsec` defaults to
  the model's pinned frame offset.
- **Handoff:** a two-plane option in `LensModel`; use `macs0416-canucs` in `images` runs as a direct cross-check.

## 2026-10-08: Abell 370 and Abell S1063 radial screens: null (D-043)
- Worktree worker: CANUCS 1208 for Abell 370 (CANUCS DR1 photo-z; frame offset −0.121″, −0.015″ pinned) and GLIMPSE
  3293 for S1063 (DJA v7.5 photo-z). Only `radial` applies (image lists gated, D-035).
- Abell 370 raised two flags (15 and 8 lines, p ≤ 0.01). Both are diffraction-spike chains from Gaia stars that are
  off the mosaic or saturated and absent from the catalogue. The Gaia-seeded veto is now in the code
  (`radial --spike-stars`, `scripts/gaia_stars.py`); with it, p 0.945. S1063: p 0.435, unchanged by the veto.
- **Failed approach:** spike vetoes seeded from the pipeline catalogue miss saturated and off-mosaic stars, whose
  spikes reach 37″.
- Wall time: under 2 min of pipeline per field, plus cutout vetting.
- **Handoff:** a low-weight veto or an aper50 S/N floor in `radial`; use `--spike-stars` on every field.

## 2026-10-08: MACS0416 counter-images and flux ratios: null; CANUCS model as second model; two-plane check (D-042)
- With the D-040 solver, the MACS0416 image list is open. `images` and `fluxratio` raised 21 flags; 0 survive. The most
  persistent was system 51's fourth image: two independent models predict it, and the 51.1–51.3 photometry says it
  should appear at 9–25σ, but it is not seen. It comes from a z 0.268 foreground galaxy modelled as a member. A
  scratch two-plane model with that galaxy at σ ≤ 60 km/s gives exactly the 3 observed images.
- System 27's two bright `absent` predictions (S/N 247–341) are CATS-only galaxy-scale caustics. The JWST-era CANUCS
  Lenstool model (222 potentials, 111 spectroscopic systems) reproduces system 27 with exactly its 3 images, and the
  cutouts show empty sky there.
- Flux-ratio flags: 45.2 is blended with a bright galaxy 0.5″ away; 38.1 is marginal.
- Wall time: about 3 min of pipeline plus about 5 min of vetting.
- **Handoff:** pin the CANUCS models (MACS0416, MACS1149, Abell 370) as `MODELS` entries for second-model vetting.

## 2026-10-08: Sunrise transient candidates `n0022` and `n0150` are detector persistence (D-039)
- New `scripts/persistence_check.py`: per-exposure photometry on level-2 `_cal` files (S3 byte ranges), plus the
  same detector pixel in earlier exposures on that detector. Validated on a synthetic afterimage and on two real
  sources (`n0153`, Earendel: 0 suspect detections).
- `n0022`: afterimages of a bright galaxy (2.0 % and 0.8 % of its flux) in o010 dithers 3 and 4; the dither geometry
  stacks both on one sky position. `n0150`: in all three epochs an afterimage (0.04–0.07 %) of a saturated star, one
  exposure later; its "motion" is the dither vector. Cutouts inspected (dithers 1–4 at fixed detector pixels).
- `n0153` is on sky in all 8 dithers and F150W is ×1.9 brighter in 2025 at the `_cal` level: still a variable
  candidate.
- **Lesson:** a single-epoch source in a mosaic can be an afterimage that two dithers place on one sky position.
  The archived `_cal` files carry no DQ flag there.
- **Handoff:** `/vet-candidate n0153`; persistence-check new single-epoch candidates before vetting.

## 2026-10-08: MACS0717 screens: null (D-041)
- Worktree worker on VENUS 6882 o029 (10 bands, the only public NIRCam association) with `macs0717-cats` (rms 3.21″).
- 51 flags, 0 surviving. 29 are the model's own copies of catalogued images it does not reproduce (1.6–5.6″ off), 6 are
  CATS-only extra images (Sharon v4cor and Keeton v4 predict none), and 5 have μ more than 2× model-dependent.
  System 65's flux ratio is a 0.6″ catalogue offset. `radial` p ≥ 0.70.
- **Failed approach:** a fixed 1.5″ match radius for a 3.2″-rms model makes most catalogued images "unpredicted",
  and their model copies then flag as absent or confused.
- Wall time: about 7.5 min of pipeline plus about 20 min of vetting.
- **Handoff:** in `forced_check`, rms-scaled match and search radii and automatic copy classification; Sharon v4cor and
  Keeton v4 as pinned `MapLensModel` entries.

## 2026-10-08: `find_images` fold refinement; frame offsets for map models (D-040)
- Cells on a critical curve near the source are subdivided into ≤ 0.02″ sub-cells. MACS0416 system 26 is now solved
  (rms 1.57″ → 0.76″), and its image list is open with offset (0.208, −0.025). No other model changes beyond 0.03″
  rms; Abell 2744 `validate` takes 41 s instead of 29 s.
- **Bug fixed:** `apply_frame_offset` was a no-op for map models (the maps are looked up by sky position), and
  `radial` / `arcs` skipped it for them. No earlier result used a map model with an offset. MACS0416's radial screen
  was re-run in the JWST frame: still null.
- **Handoff:** `images` / `fluxratio` on MACS0416 (CANUCS photo-z); radial on Abell 370 and Abell S1063.

## 2026-10-08: MACS0416: system 26 is a solver-grid artefact; radial null (D-038)
- Worktree worker on CANUCS 1208 (`jw01208-o004_t002`, 8 bands). The CATS system-26 residual (11″) is the 0.25″
  `find_images` grid missing a merging pair near the critical curve; the source lies 0.001–0.005″ from the caustic. On a 0.1″ grid the rms is 0.811″ (1.13×
  quoted). The model stays map-only until the solver refines its grid near high |μ|.
- `radial`: 98 centres against 99.6 random; the 7-line centre (p 0.29) is low-S/N noise segments. With CANUCS DR1
  photo-z (background cut): 77 against 85.0, p 0.225. 0 flags.
- **Failed approach:** jwst 3.0.0 isophotal S/N admits noise segments (66 of 137 anti arcs have aper50 S/N < 3). The
  radial screen needs an aper50 S/N floor.
- Wall time: about 4.5 min (lens and exotic scripts plus vetting stamps; no pipeline `run`).
- **Handoff:** adaptive grid refinement in `find_images`; an aper50 S/N floor in `radial`.

## 2026-10-08: CANUCS DR1 photo-z for MACS1149 (D-037 addendum)
- **Failed approach:** the worker checked only DJA, but docs/landscape.md already listed CANUCS DR1 (PSF-matched
  EAzY photo-z for A370, MACS0416, MACS0417, MACS1149 and MACS1423). Check docs/landscape.md before reporting that a
  field has no photo-z.
- MACS1149 system 16.2: z_phot 2.25 (95 % 0.23–2.33). This disfavours the CATS-fitted 4.419.
- Radial re-run with the background cut: 12 peaks against a null mean of 12.8; null. 6 of 68 matched lensed images
  get a blended z < 0.6.

## 2026-10-08: MACS1149 screens: null (D-037)
- Worktree worker on CANUCS 1208 (8 bands) with `macs1149-cats`. The image-plane rms is 0.673″ (gate passes) and
  the frame offset is under 0.02″.
- 3 flags, 0 surviving: system 16 (fitted z 4.419; at z 2.5 in CATS or 3.0 in Sharon v4cor the image lands on 16.2), system 2 (CATS topology error,
  checked against Sharon v4cor), Refsdal-host knot 1192 (next to the BCG; μ differs >2× between models).
  `fluxratio` 0 flags; `radial` p = 0.91.
- **Failed approach:** forced photometry with an annulus on a BCG core gives negative fluxes; use high-pass.
- No DJA mosaic for MACS1149 (v7).
- Wall time: about 5 min of pipeline plus about 10 min of vetting.
- **Handoff:** `macs1149-sharon` (or the CANUCS models) in `MODELS`; the D-037 rules in `forced_check`; system 16's
  redshift (spectroscopic z, or forced photometry on its far-image track); `images` with CANUCS photo-z.

## 2026-10-08: Abell 2744 screens: null (D-036)
- A worktree worker ran `validate`, `images --forced-image`, `fluxratio` and `radial` on Bergamini+2023b, with
  cutouts of every flag. 15 of 16 flags were ordinary. The survivor, 4.2c, was re-run under the D-034 rules: it
  is `no_reference` (resolved-knot references, BCG halo), not a candidate.
- Radial: 35 peaks against a random mean of 50.4, p ≥ 0.575.
- **Lesson:** a "252σ absent" image can come from a resolved-knot reference. D-034's compact-reference rule now
  catches this.
- Family 4's c images were vetted (4.1c 4–8× underluminous after BCG subtraction; 4.2c undetected). They are
  explained by μ systematics next to member 34423 (3.9–28.7 under ±30 %; CATS 7.3). **Rule:** a μ that moves by
  more than 2× under member perturbation is untestable.
- **Handoff:** HFF field runs (D-035 models), `bayes.dat` position spreads, UNCOVER v2.0 cross-check.

## 2026-10-08: Sunrise third epoch with calibrated significances; two open transient candidates (D-027)
- `transient_forced.py --controls`: noise scale from ordinary sources, applied before thresholding.
- Sunrise o010 against VENUS o052 (2.9 yr; F150W, F444W): scales 1.30/1.18; **0 of 57 catalog candidates pass**;
  Earendel steady (Δm ≤ 0.12 mag, < 1σ).
- Among the 200 controls, `n0022` (gone after 2022-07 in four SW bands), `n0150` (a different position in each
  epoch) and `n0153` (+0.5–0.6 mag in both bands) change. Not vetted; ordinary explanations (supernova, moving
  object, AGN, edge artefact) are untested. docs/fields/sunrise.md has the numbers.
- **Failed approach:** the catalog stage plus the two-band rule misses single-pair, blue transients.
- Follow-up #53: an empty control selection is flagged uncalibrated instead of being skipped silently.
- **Handoff:** `/vet-candidate` for `n0022` and `n0150` (level-2 `_cal` exposures per filter: is the source in
  every dither? Does `n0150` move within one visit?); then a grid of forced photometry (all sources, not only
  catalog candidates) for every epoch pair.

## 2026-10-08: Six HFF clusters as CATS map models (D-035)
- MACS0416, MACS1149, Abell 370, MACS0717, Abell S1063 and Abell 2744 (CATS v4/v4.1) are pinned map models. All six
  reproduce their published z = 2 magnification maps (median 9e-5 to 6e-3).
- Image lists pass the rms gate for MACS1149 and MACS0717. MACS0416 has a `params.txt` but fails the gate (system 26).
  Abell 370, Abell S1063 and Abell 2744 have no `params.txt`, and their placeholder redshifts give 9–12″
  residuals. All four stay map-only.
- **Failed approach:** using CATS `arcs.txt` redshifts as published.
- **Handoff:**
  - per-cluster JWST field runs (configs, catalogs, `radial` on all six, `images` / `fluxratio` on MACS1149 and
    MACS0717);
  - `/vet-candidate` for Abell 2744 4.2c (the field worker's survivor).

## 2026-10-08: El Gordo counter-images and radial screen: no candidate (D-034)
- Forced photometry: 3 `absent` images on the first pass, all ordinary on vetting (cutouts plus numbers):
  - 23: model position error;
  - 6: reference on a galaxy wing; photo-z-consistent counterparts at 2.6–3.1″;
  - 7: HST→JWST frame offset of 0.22″.
- Radial screen: max 4 lines, p = 0.64. A 6-line "centre" was a star's diffraction spikes.
- Four new rules (D-034): frame offset, compact and consistent references, a residual-scaled search radius, and a
  spike-segment veto. SMACS re-run: unchanged, 0 absent.
- Tally: El Gordo screened 17 uncatalogued predicted images and 37 anti arcs; flags 3 + 1; **surviving 0**.
  Wall time: about 50 s for `images`, 21 s for `radial`.
- **Handoff:** Abell 2744 (worker running), then the HFF/RELICS map fields via `MapLensModel` (#47).

## 2026-10-08: Map-based lens models; Sunrise radial screen null (D-033)
- `lensmodel.MapLensModel` evaluates published deflection maps. The first one is `whl0137-relics-lenstool`, which
  reproduces RELICS κ to 3.5e-5 and μ(z = 6.2) to 2.0e-5 (medians).
- Sunrise `exotic_screens radial`: 29 usable anti arcs, 3 centres against a null mean of 2.2, max 3 lines, p 0.885.
  **No candidate.** Wall time 17 s, with cached maps and catalogs.
- Cycle tally:
  - SMACS + Sunrise screened: 60 + 0 images, 34 + 29 arcs;
  - flags: SMACS fluxratio 2 (removed by the two-band compactness rule), radial 0 significant;
  - surviving vetting: 0.
- **Handoff:** more map fields (HFF: Abell 2744, MACS0416, MACS1149, Abell 370; RELICS clusters), El Gordo and
  Abell 2744 `images --forced-image` and screens.

## 2026-10-08: Flux-ratio and colour test of catalogued images: no anomaly in SMACS or El Gordo (D-032)
- New: `lens_consistency.py fluxratios`.
  - It compares each image's DJA `mag_auto` + 2.5 log|μ| and its F150W−F444W colour with the other images of its
    system (leave-one-out).
  - It drops blends, segments larger than 20,000 px and counterparts whose photo-z excludes the system redshift.
- Results (`derived`):
  - SMACS: 8 images in 4 systems flux-tested, rms 0.60 mag.
  - El Gordo (0.3″ after a +0.221″ frame shift): 23 images in 11 systems, rms 0.49 mag.
  - No colour outlier (rms 0.11 and 0.07 mag).
  - The flagged pairs SMACS 6.1/6.3 and El Gordo 18b/18c (and 7b/7c, at 0.5″ without the shift) fail vetting.
    Forced photometry brings SMACS 6 down to 0.61–0.78 mag. The El Gordo pairs are low-S/N or chromatic, so they
    are measurement failures. **No anomaly.**
- **Failed approaches**, each checked on cutouts or SEDs:
  - DJA aperture × `tot_corr` is not a total flux for arcs.
  - SMACS 1.1 looked 1.4 mag too bright, from host-halo light in a 32,864 px segment.
  - El Gordo 9a/9c differed by 1.2 mag in colour, because 9a's DJA counterpart is a z_phot 0.89 interloper.
  - Bare eazy 95 % intervals exclude good images.
  - A median that includes the image itself halves pair differences.
- Limit: DJA misses most arcs inside cluster light (SMACS: 38 of 60 images unmatched within 0.3″).
- **Handoff:**
  - El Gordo `images --forced-image` (TASKS "Now" 2);
  - Abell 2744 `fluxratios` with the 233 MB DJA catalogue (needs a DECISIONS entry);
  - `bayes.dat` μ uncertainties, then recheck SMACS 6;
  - BCG/ICL-subtracted totals, to test the core images.
## 2026-10-08: Exotic-lens screens; SMACS null (D-031)
- New `scripts/exotic_screens.py`:
  - `fluxratio`: two-band forced photometry, luminosity ratio against sibling images, compactness and chromatic
    gates;
  - `radial`: anti-tangential arcs whose axes converge on a dark centre, with a false-alarm rate from randomised
    position angles.
- SMACS:
  - fluxratio: 6 compact images consistent, 0 flags. A first, single-band compactness gate had flagged system 7,
    an ordinary knot-vs-whole-arc mismatch;
  - radial (background sources only): 12 centres against a null mean of 9.0 (p95 15); max 4 lines, p 0.945.
  - **No exotic candidate.**
- **Failed approaches (now rules):**
  - fixed-aperture flux ratios on resolved arcs: surface brightness is conserved, so the ratios scale with 1/|μ|
    (systems 5 and 10);
  - a single-band compactness gate (knots of a clumpy arc pass in F150W);
  - a uniform-angle null for the radial screen (the selected arcs point at the mass centre).
- **Handoff:** fan out per cluster (El Gordo, Abell 2744, Sunrise, then HFF/RELICS). Each runs `validate`,
  `images --forced-image`, `exotic_screens fluxratio` and `radial`, with contact sheets of all flags.

## 2026-10-08: `find_images` 3–5× faster with identical images
- Seeds are pre-filtered with boolean sign tests on the mapped grid corners, and Newton steps run for every seed in
  one `fields_xy` call.
- Benchmark: all catalogued systems at a 0.25″ grid, images identical (max |Δ| 0 arcsec, same counts):
  - SMACS: 23.2 → 6.7 s;
  - El Gordo: 44.3 → 13.0 s;
  - Abell 2744: 118 → 23.5 s (at ±190″).
- **Handoff:** the deflection grid itself (one-time and cached) is now the main cost. Published deflection maps
  (UNCOVER, RELICS, HFF) could replace it for fields without a Lenstool model (TASKS).

## 2026-10-08: Image-plane χ² reproduces Lenstool for SMACS, El Gordo and Abell 2744 (D-030)
- Merged #40 (counter-images, D-029) after its last commit, which GitHub had not attached to the PR, was picked up
  by a follow-up commit.
- `validate` now computes the exact image-plane χ² for every model; El Gordo and Abell 2744 are `MODELS` entries.
  Issue #41's numbers are now reproducible in the repository: χ² 30.87/30.91 (SMACS), 82.53/80.22 (El Gordo),
  146.64/146.60 (Abell 2744). No catalogued image is off by more than 3σ in any field.
- Parser fixes: letter-suffixed image ids, several ids per `z_m_limit`, 6-decimal `_kpc` rounding.
- **Failed approach:** the per-image error column for El Gordo (χ² 52.0); a uniform 0.621″ matches Lenstool.
- CDS reset connections through the proxy in this run; the cache was seeded from a `curl` copy with the same sha256.
- **Handoff:** run `lens_consistency.py images` with forced photometry on El Gordo and Abell 2744 (TASKS "Now" 2);
  then the exotic screens (docs/exotic_lensing.md).

## 2026-10-08: SMACS counter-images: no predicted image is absent (D-029)
- New: `lensmodel.find_images` (an image-plane solver on a cached deflection grid), `lens_consistency.py images`, and
  `--forced-image`, which runs forced photometry on S3 byte-range stamps.
- Results:
  - ICLv2 reproduces all 60 catalogued images within 0.04–0.91″, with 4 demagnified central images.
  - Of the 11 testable uncatalogued images: 3 are recovered (systems 9 and 8, and system 17 marginally: flux ratio
    0.36 on the BCG gradient), 1 is
    confused, 1 is undetectable, 6 have no reference flux, and **0 are absent**.
  - System 8's model z = 11.76 is contradicted by its F090W/F150W detections.
- **Failed approach:** pipeline-catalog flux references near cluster galaxies. They falsely made system 9 `missing`.
- The Mahler κ tarball now comes from raw.githubusercontent.com, because github.com/raw returns 403 in cloud runs.
- **Handoff:** TASKS "Now" 1 is done for counter-images. Next:
  - parser fixes (issue #41), then the same test on El Gordo and Abell 2744;
  - the exotic screens of docs/exotic_lensing.md (demagnified images, radial images around a dark centre).

## 2026-10-08: El Gordo and Abell 2744 lens models, preliminary: no image-position anomaly (issue #41)
- Preliminary validation on published products. The numbers come from scratch code and are not yet reproducible
  in the repository; they need a re-run with #40's solver.
  - `lensmodel.py` reproduces the Caminha+2023 El Gordo magnification maps (median |Δμ|/μ 4e-5).
  - Exact image-plane solves match Lenstool's χ²: El Gordo 82.5 vs 80.22 (rms 0.754″, paper 0.75″); Abell 2744
    (Bergamini+2023b) 146.64 vs 146.60 (rms 0.427″).
- Results:
  - No image-position anomaly in either field.
  - Predicted uncatalogued counter-images: El Gordo 16, inconclusive, because the MUSE Lyα images are too faint in
    continuum. Abell 2744 30, not yet testable, because the pipeline catalogs miss most arcs in the core.
  - docs/fields/elgordo.md, docs/fields/abell2744.md.
- **Failed approach:** source-plane back-trace χ² for image-plane-optimised models (El Gordo 121.6 against 80.22).
- **Bugs found on `main`:** the parser fails on letter-suffixed image IDs and on 6-decimal `_kpc` rounding (#41).
- Also merged #38 after bringing it up to date with `main`.
- **Handoff:** #40 (image solver, local session); then reproduce these numbers in the repository, the parser fixes,
  and the counter-image flux test on DJA photometry (TASKS "Now" 2).

## 2026-10-08: Sunrise two-epoch search, a null result; Earendel steady (D-027 amendment)
- WHL0137 2282 o010 against o120 (same pipeline): 372 catalog candidates → 58 in ≥ 2 bands → 52 after the bright-Gaia
  mask → 0 by recentred forced photometry in all four bands. 2 pass in one or two bands, both ordinary on visual check
  (an epoch-2 streak; galaxy outskirts at the footprint edge). docs/fields/sunrise.md.
- Earendel: |Δm| ≤ 0.18 mag (≤ 1.5σ) in F090W/F115W/F277W/F356W over 164 days.
- Noise calibration on 150 ordinary sources: ERR-based significances are 1.2–1.5× too large.
- New: `scripts/transient_combine.py` (band coincidence and bright-Gaia mask) and
  `transient_forced.py --recentre-arcsec`.
- **Failed approach:** a fixed aperture at a literature position. Earendel seemed to brighten by 0.76 mag at 5.7σ,
  because the aperture sat 0.15″ off the source and the PSF rotates by about 180° between the epochs.
- **Handoff:** lens-model checks per field (#35 merged; El Gordo and Abell 2744 model validation in issue #41).
  Time domain: rescale the thresholds by the control std, and add VENUS 6882 o052 (F150W, F444W) as a third epoch.

## 2026-10-08: Cloud runs use the GitHub MCP tools; conflicting PRs get no CI (D-028)
- The first routine run opened #38 with the session's GitHub MCP tools. The prompt, the research-cycle skill and
  docs/operations.md §3 now prefer them, merge with `mcp__github__merge_pull_request` after checking the policy,
  and fall back to REST plus the `merge-ready` label.
- **Failed approach:** #37, a REST squash-merge allow rule, was closed. The glob `pulls/*/merge` also matches other
  `pulls/...` writes, so the Bash merge API stays denied.
- **Lesson:** a PR that conflicts with `main` gets no `pull_request` CI at all (#38). Merge `origin/main` in first.
- The cloud merge gate checks author, head repository (no forks), branch, labels, every page of files and the CI jobs
  (not the skipped `claude` runs), and pins `expectedHeadSha`. D-028 records the decision.
- **Handoff:** TASKS "Now" 1 (counter-images, then parity and flux ratios; see the lens-model entry below). #38
  needs `origin/main` merged in before its CI and merge.

## 2026-10-08: Lens-model stage validated; SMACS arc orientations agree with ICLv2 (D-024)
- `lensmodel.py` is a Lenstool dPIE port (from PyAutoGalaxy, MIT) that evaluates a published `best.par`. It
  reproduces Mahler+2022 ICLv2:
  - κ map: median |Δκ| 1.6e-5;
  - back-trace χ² of the 60 catalogued images: 31.18, against Lenstool's 30.91.
- `scripts/lens_consistency.py` provides `validate` and `arcs`.
- SMACS arcs: 21 of 25 elongated strong-shear background sources, each tested at its own photo-z range, are
  aligned with the predicted stretch (p = 2.6e-7). The 6 anti candidates are all ordinary:
  - 2 segmentation blends;
  - 1 galaxy at z ≈ 0.77 with an intrinsic shape;
  - 3 noisy low-surface-brightness shapes.
  This is a null result.
- **Lesson:** moment orientations from pipeline segments pick up blends and low-S/N shapes, so filter them before
  calling a source anti-tangential.
- **Handoff:** TASKS "Now" 1, which leaves counter-images (predicted but missing, or observed but not predicted)
  and parity/flux ratios. Then El Gordo and Abell 2744; check their model profiles first.


## 2026-10-08: Hourly cloud routine
- Routine `trig_01PNAmgcfqef8CvhPAY8ggbP` runs `/research-cycle` loops hourly at :07 UTC (Opus 5.5, no connectors)
  with the prompt in docs/cloud-routine-prompt.md.
- First test run: clone, venv, skills and MAST downloads work.
  - GitHub GraphQL is blocked in cloud sessions, so `gh pr` / `gh issue` fail. The REST equivalents and a single
    background CI-wait loop are in docs/operations.md §3, and the research-cycle skill points to them.
  - The merge API is denied by `.claude/settings.json`, so cloud runs stop at merge-ready with the label
    `merge-ready`. The owner or a local session merges.
  - Coordination: runs last about 40 minutes and skip `local-wip` PRs and branches committed to in the last 15
    minutes. An unmerged PR's handoff is read from its branch.
- **Handoff:** unchanged (TASKS "Now (M3)"); the cloud routine continues it.

## 2026-10-08: Two-epoch transient search, a null result; exotic-lens signatures (D-027)
- `scripts/transient_search.py` produces catalog-level candidates; `scripts/transient_forced.py` checks them with
  forced aperture photometry on byte-range cutouts of both epochs.
- SMACS against VENUS: catalog 441 (F444W) and 628 (F150W) → 144 in both bands → 4 by forced photometry (6 with
  WCS-centred apertures) → 0 after a visual check. All are bright stars or a source next to a bright neighbour, where
  fixed-aperture systematics dominate. No credible transient at |Δm| ≥ 0.3 mag, ≥ 5σ.
- **Failed approach:** catalog-only comparison across pipeline versions, which picks up deblending differences.
- docs/exotic_lensing.md records the verified exotic-lens signatures (wormhole demagnification, radial images for
  negative mass), the degeneracies with ordinary lensing, and existing limits. Warp drives have no imaging
  prediction to test.
- **Handoff:** the lens-model stage (TASKS "Now (M3)" 1), then Sunrise's same-pipeline epoch pair.
## 2026-10-08: Three cluster fields; robustness rules from them (D-025, D-026, cycle 17)
- Parallel worktree workers added Abell 2744 (#32), El Gordo (#31) and Sunrise/WHL0137 (#30). PR #29 fixed downloads
  for MAST products whose listed size is stale after the 2026-10-01 reprocessing. Proposals are folded into
  D-026 and SOURCES.
- **D-025** (found on Abell 2744):
  - catalogue "stars" with DJA r50 too large for a point source (bright cluster-galaxy cores in Gaia) go back to the
    galaxy ranking, before the locus calibrates;
  - very red spiky sources are exempt from the D-020 host test.
  - Run `20261007T211229Z-883a5066`: 8 vetoed, the locus recalibrated, the star top 10 is now mostly point sources with spikes, and
    `5904` is kept.
- **Handoff:** TASKS "Now (M3)". The lens-model stage (worker, `claude/lens-model`) is in flight; then per-field
  model checks and the transient forced photometry.


## 2026-10-08: Owner priorities, lensing-violation search first (D-023)
- The owner merged PR #20 (`2804`, inconclusive) and reprioritised: M3 lens-model consistency first; more clusters in
  parallel (Abell 2744, El Gordo, Sunrise; MACS J0416 has no DJA v7 catalog); a two-epoch transient search; no
  pauses; one review per PR.
- **Handoff:** TASKS "Now (M3)" items 1–3. Worker PRs with label `batch-clusters` are in flight.


## 2026-10-08: CEERS control on DJA matched photometry (D-022, cycle 16)
- `ceers_t021` now uses the DJA v7.4 `ceers-full` catalog, a 250.5 MB download with the reason stated in D-022.
  The control field is built the same way as SMACS for colours, confirmation, screening and the top-k metric.
- Run `20261007T202419Z-65728e80`:
  - 3,069 of 7,582 sources match DJA, and 3,156 pass the gate (2,757 before);
  - the top 20 is 85% known objects with 0% cutout-flagged, and 4 sources were screened;
  - the stellar locus does not apply (2 catalogued stars).
- **Found:** `f200w_2842` (#4) is a single-band, ellipticity-0.95 streak, probably a spike from the bright source
  4.6″ away. A coincident DJA object confirmed it (D-014 "Revisit if").
- **Handoff:** TASKS "Spike streaks confirmed by a coincident DJA object", then the star-stratum NED matches.


## 2026-10-08: Low-weight screening of the top-k pool (D-021, cycle 15)
- `stages.cutouts.screen_low_weight` re-applies D-011's 0.5 weight threshold with the cutout's own core weight. It
  runs on every non-star stratum, and `screened.ecsv` now records each source's `reason`.
- Run `20261007T194733Z-ac5f1b59`: the stripe artifact `3034` (#3) is removed: short-wavelength only, not in DJA, core weight 0.468.
  The galaxy top 20 is now 40% known objects, with 1 flagged source (`2915`, kept on purpose) and 4 screened
  (one by D-021, two by D-019, one by D-020).
- **Handoff:** TASKS "Star-stratum sources with NED galaxy matches", then MIRI's `G_Lens` vetting.


## 2026-10-08: Host test for spiky sources (D-020, cycle 14)
- `cutouts.host_ratio` measures annulus light against the peak around the spike peak. D-019 screening now also
  removes spiky sources with no host light (< 0.004), so the two-star segment `940` (0.0017) leaves the top 20.
  `2915` (0.0076, a galaxy nucleus) stays.
- Run `20261007T190843Z-db3c83fe`: the galaxy top 20 is 40% known objects, with 1 spiky source (`2915`, legitimate) and 3 screened.
- **Next contamination:** `3034` (#3) is a diagonal stripe artifact in a low-weight region (cutout flag `low_weight`).
  It passed the D-011 gate.
- **Handoff:** TASKS "Stripe artifact at galaxy rank 3".


## 2026-10-08: Top-k composition metric; `940` re-diagnosed (cycle 13)
- Every stratum now reports its top-k composition: known, lens-related, catalogued star, cutout-flagged,
  spikes and screened. It is stored in `run_record.json` as `samples[].topk` and defined in
  docs/methodology.md. The crossmatch rows now keep `is_lens_related`.
- Baseline, run `20261007T183007Z-4342c5c2`: the galaxy top 20 is 35% known objects, 10% lens-related, 0% stars and 15% flagged.
- **`940` is not saturated.** Its F200W cutout has no no-data pixels. The pipeline merged two stars into one
  5,788 px segment whose centroid lies 0.20″ from the brighter star's peak, so aper50 at the centroid is faint
  (25.7 against isophotal 20.1). It is a single-band detection confirmed by the DJA match (0.18″) to the bright
  star.
- **Failed approaches** (saturation tests):
  - aper50 − isophotal magnitude: cluster galaxies reach 10–21 mag, so `940` (5.6) is not an outlier;
  - a no-data core within 0.3″: `940` has none, and only `1345` does.
- **Handoff:** TASKS "Star-pair segments in the galaxy ranking" (a host-vs-PSF profile test or a centroid–peak
  offset rule) and "Star-stratum sources with NED galaxy matches".


## 2026-10-08: Vetted `f200w_2915`, a catalogued quiescent galaxy with a compact core (cycle 12)
- Run `20261007T174424Z-540be92d` (clean `a27d628`), galaxy rank 21, kept ranked by D-019 as a spiky source with non-stellar colours.
- It is SIMBAD `[VBG2023] SMACS 1060` (Valentino+2023 atlas of colour-selected quiescent galaxies at z > 3). DJA
  z_phot is 2.84 (95% 2.68–3.10), mildly below that selection. Its F090W−F150W break is 2.5 mag. The spikes come from its compact core. D-019 behaved as
  intended.
- **Two epochs:** the position residuals disagree between bands (7 mas F444W, 34 mas F150W), as expected for centroid
  shifts of an extended source, and Δm is F444W only (+0.11; F150W +0.01). So it shows neither motion nor
  variability. `epoch_compare.py` now documents that its significance assumes a point-like target.
- **New resource:** program 4043 has NIRCam F444W grism spectra over part of SMACS 0723 (SOURCES).
- **Handoff:** TASKS "Saturated stars in the galaxy ranking", then the star-stratum NED matches.


## 2026-10-08: Spike screening of the galaxy top k (D-019, cycle 11)
- The galaxy strata get cutouts for 2k. Spike-flagged sources with stellar colours (`1571`, `2242`) are removed
  and backfilled. Spiky sources with non-stellar colours stay ranked and are noted: `940` (a saturated star with
  corrupted colours) and `2915` (a bright nucleus, possibly an AGN).
- **Failed approach:** screening every spiky source removed `2915`. Spikes also mark bright galactic nuclei.
- Run `20261007T170539Z-b369bfa2`: the top 20 now runs to rank 22, with no star+galaxy blends left.
- **Handoff:** TASKS "Vet `2915`" (a bright red nucleus), "Saturated stars in the galaxy ranking" and "Star-stratum
  sources with NED galaxy matches".


## 2026-10-08: Diffraction-spike flag on cutouts (D-018, cycle 10)
- `cutouts.spike_statistic`: hexagonal-harmonic power around the brightest peak near the target. NIRCam cutouts with
  `spike_s6` ≥ 3 get the `spikes` flag, and the report lists them.
- Run `20261007T162133Z-2511977c`: it flags exactly the three PSF-like blends in the galaxy top 20 (`940`, `2242`, `1571`) and all of
  the star-stratum top 10. `940` is really a saturated star: its pipeline aper50 is 5 mag fainter than DJA's
  measurement, and the stellar locus missed it because its colours are corrupted.
- **Failed approach:** centring on the catalog centroid, which misses blends whose star is offset.
- PR #20 (`2804`, needs-human) was revised after review: F410M−F444W is now an observed colour, and the line
  is a hypothesis.
- **Handoff:** TASKS "Spike-flagged sources in the galaxy ranking" (do flagged sources leave the top k?) and
  "Star-stratum sources with NED galaxy matches".

## 2026-10-08: Vetted `f200w_2804`, inconclusive, a little-red-dot-like source (cycle 9b)
- Run `20261007T151316Z-5706ae30` (clean `6fade35`), galaxy rank 6. Uncatalogued, unresolved, F444W 22.9.
- **Second epoch:** no proper motion over 4.0 years (VENUS 6882, 2026-06-05; 1.9 mas = 0.95σ, so
  |μ| ≲ 1.5 mas/yr) and no F444W variability (0.004 mag).
- F444W is 0.725 mag brighter than F410M in the same epoch. That is an observed colour; a line at 4.25–4.98 µm
  and F410M absorption are both hypotheses.
- It meets Kokorev+2024 little-red-dot colour cuts. A brown dwarf is disfavoured on two counts: the 1.5–2.8 µm
  colours are red rather than blue, and there is no motion.
- Verdict: `inconclusive: needs NIRSpec spectroscopy`. Record: `docs/candidates/jw02736-o001_t001_nircam_f200w_2804.md`.
- **Handoff:** the owner decides whether to pursue it (`needs-human`).


## 2026-10-08: Two-epoch vetting tool (D-017, cycle 9a)
- `scripts/epoch_compare.py`: proper motion and variability of one target from two public level-3 catalogs of one
  filter, with a local median frame tie. Program 6882 (VENUS, 2026-06-05) gives SMACS 0723 a 4-year second epoch in
  F150W, F182M, F210M, F300M, F410M and F444W (SOURCES "Second epochs").
- **Handoff:** the vetting record for `f200w_2804` follows in its own PR (inconclusive, `needs-human`).


## 2026-10-08: One-sided stellar-locus size test (D-016, cycle 8)
- Bright calibration stars have inflated r50 (5–8 px under 20.5 mag), so D-015's two-sided band around 2.86 px
  missed the unsaturated stellar sequence at r50 ≈ 2.0 px. The size test is now 0.5–1.2 × r50_psf.
- Run `20261007T144011Z-11342d69`: the locus adds 90 stars, not 16. The star stratum has 109; its top 10 is all PSF-like.
  `f200w_1874` is now a star.
- Point-like sources left in the galaxy top 20 are faint (mag > 24), blends, or `2804`, whose red
  F356W−F444W makes it a brown-dwarf-like candidate to vet. Details are in D-016.
- **Handoff:** TASKS Now 1 (vet `2804`), Now 2 (blended stars).


## 2026-10-08: Vetted four high-ranked SMACS galaxies; DJA photo-z for vetting (cycle 7)
- PR #16 (D-015) review fixes merged. Run `20261007T132538Z-96e911bb` (clean `725eac3`) reproduces the D-015 numbers.
- `scripts/vet_evidence.py --photoz`: the nearest DJA eazy entry, labelled model_prediction (SOURCES "Vetting
  (cycle 7)"). `/vet-candidate` now names the script and the real `candidates` subcommands.
- All four are ordinary; none is an artifact. Records are in `docs/candidates/`.

  | Source | Rank | Verdict |
  |---|---|---|
  | `829` | 3 | red, dusty galaxy at spectroscopic z = 2.74 (NED), a published F150W-dropout "z ≈ 11–20" candidate |
  | `1032` | 8 | lensed arc `[MJR2023] 028.2` (Mahler et al. 2023), a tangential pair with 028.1; blue colours likely from emission lines (hypothesis) |
  | `438` | 11 | uncatalogued; Balmer-break galaxy at z_phot 4.8 (hypothesis, 6-band photo-z) |
  | `1243` | 12 | catalogued red point-like high-z candidate, photometric z ≈ 5.6–5.75 |

- **Lesson:** after D-011/D-014/D-015 the top of the galaxy ranking is dominated by real but rare populations
  (lensed arcs, dusty and high-z galaxies), most already catalogued. The pipeline recovers known lensing
  features. Uncatalogued sources such as `438` are the discovery set.
- **Side finding:** the D-015 stellar locus missed the PSF-like star `f200w_1874` (F200W 21.4 mag, spikes)
  because its r50 (2.10 px) is 27% below r50_psf (TASKS Now 1).
- **Handoff:** TASKS Now 1 (stellar-locus lower bound) and Now 2 (blended stars).


## 2026-10-08: Stellar locus for stars missing from Gaia (D-015, cycle 6)
- `classify.stellar_locus`: a source counts as a star when its DJA detection-image r50 is within 20% of the
  catalogued stars' median and its colours are stellar. The DJA join now carries `r50_pix` and `mag_auto`.
- Run `20261007T130207Z-c8ffc470`: the locus added 16 stars to SMACS, and 3 of the 6 PSF-like sources left the galaxy top 20.
  `829`, `1032` and `438` are ranks 3, 8 and 11. Details are in D-015.
- **Failed approaches:**
  - pipeline CI and `semimajor_sigma` do not separate stars from galaxies, because spikes and
    saturation inflate them;
  - size alone admits compact galaxies, because their colours are too broad.
- **Handoff:** TASKS Now 1 (vet `829`, `1032`, `438`, `1243`). Three blended stars remain in the galaxy top
  20 (TASKS Now 2).


## 2026-10-07: Rank only confirmed detections (D-014, cycle 5)
- Quality gate: best-band S/N floor (`low_snr`) and single-band confirmation (`single_band`, confirmed by a
  DJA match). Too few survivors now fall back to the D-011 tests instead of ranking ungated.
- S/N is now inverted exactly from the pipeline's `abmag_err` everywhere: the features S/N ≥ 3 gate,
  the quality gate and the scripts. The linear approximation overstated S/N at low S/N.
- Run `20261007T124944Z-bbad4ab3`: spike and stripe detections in the SMACS galaxy top 20 fell from 5–6 to 1, and near-noise
  sources from about 6 to 0. `[YML2023] F150DB-C-4` and `[MJR2023] 028.2` are ranks 4 and 9. Evidence and
  counts are in D-014.
- **Failed approach:** an S/N floor in the reference band only removed red dropouts, including `829`
  (D-014).
- **Handoff:** TASKS Now 1 (stellar locus: about 5–6 faint PSF-like stars remain), Now 2 (vet `829`,
  `1032`, `438`, `1243`).

## 2026-10-07: Size-robust colours from DJA matched apertures (D-013, cycle 4)
- New `photometry.py`. The DJA v7.4 catalog is verified against its sha256 before it is cached, and is
  joined one-to-one within 0.2″. Colours now come from 0.5″ matched apertures.
- `features`: `ref_mag` and the morphology features keep the pipeline aperture's S/N gate, so sources
  without a DJA match keep their morphology. With matched photometry, a band counts as detected when DJA
  measured it at S/N ≥ 3. DAOFind sharpness/roundness apply only when 0 < CI_50_30 ≤ 1.8.
- **The size bias is gone, measured on the same 1,423 rows** (run `20261007T122054Z-c52935ec`, F200W−F277W, F200W aper50
  S/N > 10). Spearman between log area and colour:
  - aper50: +0.181 (p = 6.7e-12);
  - isophotal: −0.312 (p = 1.5e-33);
  - DJA 0.5″: −0.012 (p = 0.66).

  Re-measure with `python scripts/feature_size_bias.py --run-dir <outputs>/runs/20261007T122054Z-c52935ec --sample
  smacs0723_nircam --compare dja05`.
- Join coverage: 2,729 of 5,254 sources matched one-to-one. Nearest-neighbour matching gave 2,920, of
  which 191 were fragments sharing one DJA object.
- **Ranking:** the sources vetted in #13 dropped (`2925` to rank 1297, `2559` to 482, `1096` to 296). Two
  catalogued lens-related objects are in the galaxy top 20 (observed cross-matches):
  - rank 7 `f277w_829`, SIMBAD `[YML2023] F150DB-C-4` (also NED G_Lens);
  - rank 13 `f200w_1032`, NED `SMACS J0723-73:[MJR2023] 028.2` (G_Lens).
- **Top-20 purity is still mixed** (visual, unvetted):
  - 4 faint PSF-like stars missing from Gaia;
  - 5–6 detections on diffraction spikes or parallel stripes (a new artifact class the D-011 gate misses);
  - about 6 faint, compact or near-noise sources, one of them hot-pixel-like;
  - 3 arc-like sources;
  - 1 catalogued elongated galaxy.
- Failed approaches, both fixed after `/code-review` (15 findings, all addressed):
  - nearest-neighbour joining (fragments shared photometry);
  - gating morphology on DJA S/N, which erased the morphology of 44% of sources.
  The first run, `20261007T120528Z-378af5a1`, used both.
- **Handoff:** TASKS Now 1 (mask detections on diffraction spikes and stripes), Now 2 (S/N floor for
  ranking), Now 3 (stellar locus), Now 4 (vet the catalogued lensed images and arcs).

## 2026-10-07: First vetting: three arc-like galaxy candidates (cycle 3)
- New `scripts/vet_evidence.py` gathers vetting evidence:
  - per-band catalog rows;
  - radius and tangential alignment from both the catalog orientation and image moments (they agree
    to within 1–7°);
  - nearest multiple image in Mahler+2022 `arcs.dat`;
  - nearest star;
  - SIMBAD/NED/Gaia cross-match at 1″ and 3″;
  - six-band S3 cutouts.
- Vetted SMACS NIRCam `f200w_2925`, `f200w_2559` and `f200w_1096` (records in `docs/candidates/`; store status
  `explained`):
  - none is among the 62 Mahler+2022 constraint images (nearest 24–45″); this does not show they are
    singly imaged, which needs the model's critical curves (M3);
  - 2925 and 2559 are not tangentially aligned (45–60°); 1096 is at 19°, consistent with weak shear
    (`hypothesis`);
  - NED knows 2925 (z 1.98) and 1096 (z 1.36, with catalogued clumps) as background galaxies.
- **Systematic finding (derived):** the baseline's colour and shape features are size-biased.
  - aper50 F200W−F277W reddens with isophotal area (median +0.77 → +0.99, Spearman 0.21,
    p ≈ 1e-16), while isophotal colour does the opposite (+0.71 → +0.02). EE apertures differ in angular
    size between SW and LW.
  - Reproduce with `python scripts/feature_size_bias.py --run-dir <outputs>/runs/20261007T040938Z-01527ace
    --sample smacs0723_nircam`: sources detected in F200W and F277W with F200W aper50 S/N > 10, F200W
    isophotal-area bins [0, 50, 200, 1000, ∞) px.
  - DAOFind sharpness reaches +145σ on clumpy extended galaxies.

  This explains these three candidates (verdict `catalog effect`) and sets the next priority (TASKS Now 1).
- **Handoff:** TASKS Now 1, size-robust features (DJA matched-aperture colours; sharpness only for compact
  sources).

## 2026-10-07: Stars ranked as their own stratum (D-012)
- New `classify.classify_sources`. One bulk Gaia DR3 + SIMBAD cross-match over the whole catalog
  (about 8 s per sample) labels stars. The runner ranks them as `<sample>-stars` with their own top k;
  galaxies are ranked without them.
- Failed approach: the pipeline's `is_extended` flag labels every bright, saturated star in the top 20
  as extended.
- Real data, run `20261007T040938Z-01527ace`:
  - Stars found: SMACS 52, MIRI 12, CEERS 7. Of these, 40, 8 and 7 passed the gate and were ranked as
    `-stars`.
  - Galaxy-stratum top 20s: 0 cross-matched stars (SMACS had 7), 1/60 image-quality flags.
- Visual check (unvetted): the SMACS galaxy top 20 is now mostly interacting, clumpy or elongated
  galaxies, including three arc-like sources (`f200w_2925`, `f200w_2559`, `f200w_1096`). About 4
  faint PSF-like sources without Gaia counterparts remain.
- `/code-review` found 15 issues, all fixed. The main ones:
  - a classification leaving fewer than 2 galaxies could abort the run;
  - stars were dropped entirely when there were fewer than `min_stars`;
  - the offline suite made live CDS calls;
  - the rule used any match instead of the nearest one.
- **Handoff:** TASKS Now 1 (vet the arc-like sources against published lens models), then Now 2
  (stellar locus for faint stars).

## 2026-10-07: Quality gate before ranking (D-011)
- New `quality.assess_sources` with `cutouts.sample_weight_map`. A coarse WHT map (1″ cells,
  concurrent S3 byte ranges, about 10 s per 1.8 GB mosaic) gives each source's relative weight and
  edge distance. A CI_50_30 test
  catches sources sharper than the PSF. Ranking now sees only sources that pass.
- Real data, run `20261007T033832Z-cfcf6032`: top-20 sources with a cutout image-quality flag fell from
  20/60 to 2/60.
  - SMACS NIRCam: 6→1. MIRI: 10→0. CEERS: 4→1.
  - Gated out: 28% of SMACS NIRCam, 55% of MIRI, 28% of CEERS.
- `/code-review` found 15 issues in the first version, all fixed:
  - edge distance was snapped to cells, so `edge` never fired for LW or MIRI;
  - a gate passing too few sources could abort the run;
  - the gate ran even when unconfigured;
  - noise was labelled as artifacts.
- Visual check (unvetted): the SMACS NIRCam top 20 is now about 12 bright stars plus 8 galaxies, which
  include interacting pairs and the two arc-like sources (`f200w_2925`, `f200w_2559`). Low-weight
  noise and streaks are gone. Stars dominate, so star/galaxy separation is next (TASKS Now 1).
- The agent self-merge policy was adopted (owner decision; CLAUDE.md "Merge policy and version
  control"), and PRs #1–#10 were merged.
- **Handoff:** TASKS Now 1 (star/galaxy separation).

## 2026-10-07: M0 bootstrap batch landed; first real-data run (M1 slice)
- Nine parallel units became PRs #1–#9: archive, catalog, rank, cutouts, crossmatch, runner/CLI, OSS,
  agent harness, landscape. The integration branch merges all of them and adds a contact sheet to the
  runner. 297 offline tests pass.
- First end-to-end run on real public data: run `20261007T020124Z-4bfabaaa`, config `reference_sample_v1`,
  about 3 minutes on the owner's laptop.
  - Merged sources: SMACS 0723 NIRCam (6 bands) 5,254; MIRI (4 bands) 530; CEERS t021 (7 bands) 7,582.
  - The top 20 per sample get cutouts through S3 byte-range reads (no full i2d download) and a
    SIMBAD/NED/Gaia cross-check.
- Results. These come from the integrating agent's visual classification of the contact sheets
  plus the cross-match; nothing is vetted.
  - SMACS NIRCam top 20: mostly bright stars with diffraction spikes (7 confirmed by Gaia/SIMBAD,
    about 11 PSF-like), 2 image artifacts (a streak and low-weight noise), and galaxies or blends.
    Ranks 6 and 19 are elongated and arc-like; they are hypotheses to vet against published lens models.
  - MIRI top 20: 10/20 carry image-quality flags (edge or low weight). One known lens-related object
    (NED `G_Lens`) was recovered; it is a known object, not a discovery.
  - CEERS control: with F200W only, about 16/20 were artifacts (edges, hot-pixel-like compact sources,
    linear streaks). With all 7 bands, about 5/20 are artifacts and the rest are structured galaxies, so
    cross-band features matter.
  - Interpretation: the baseline mostly surfaces stars and instrument or processing artifacts, as
    expected. Quality gating and star/galaxy separation come next (TASKS Now 1–2).
- Injection-recovery (unit 3, NIRCam): only LOF clearly beats random on every injection type. The
  combined score reaches precision@50 of 0.06–0.20. Nothing is tuned toward one injection type until
  vetted labels exist (D-004).
- Data caveats: program 2736 catalogs were made with jwst 2.0.1 and CEERS with jwst 3.0.0, and deblending
  changed in between (D-010). The CEERS `s_ra/s_dec` is the MIRI target; use `s_region` instead.
- Process lessons:
  - `setup-uv` has no floating major tag, so pin the full version.
  - Workers collided in the shared scratchpad; each now gets its own subdirectory (CLAUDE.md).
  - Unit 5's `gh pr create` was blocked, so the coordinator opened #9 with an inline body;
    `--body-file` is denied by project settings.
- Guardrails:
  - `main` is protected: a PR and 4 CI checks are required, admins included.
  - A ruleset blocks force-push on every branch.
  - Private vulnerability reporting and Dependabot security updates are on.
  - Direct commits to `main` happened only in Step 0, before protection existed.
- **Handoff:** the owner merges #1–#9 in any order, then this integration PR (updated by merging `main`
  into it). The next cycle starts at TASKS "Now" item 1.

## 2026-10-07: M0 bootstrap started
- Created repository skeleton: interfaces (`schema.py`, stage stubs with fixed signatures), state files,
  CI, BSD-3-Clause license, verbatim agent charter (`docs/agent-charter.md`) imported by `CLAUDE.md`.
- Verified on MAST: program 2736 level-3 NIRCam/MIRI imaging is public; catalogs ~3 MB vs NIRCam i2d
  ~1.8 GB, hence the catalog-first slice (D-001). The S3 mirror key pattern is confirmed.
- Owner decisions: public repo; every change via PR and the owner merges; chat in Chinese.
