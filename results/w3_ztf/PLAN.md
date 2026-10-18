# W3 (negative-mass point lens) feasibility trial on ZTF DR24 light curves — plan

Owner instruction 2026-10-18 (after #333). Budget 11:05Z → **15:05Z** (session clock; commit order is the record).
Written before any light curve of the sample was read (only the objects table and partition metadata were read).

## 1. Not a repeat
- W3 was searched in OGLE-IV (D-057–D-059), Gaia DR3 (D-061) and MOA-II (D-062, D-068): bulge fields, difference
  fluxes or PSPL-selected lists. D-054 set ZTF aside for the **bulge** (ineffective there). Not tried: ZTF light
  curves at all, a high-latitude northern field, and a detections-only catalogue where a vanished star shows up as
  missing epochs. New here: the ZTF DR24 data (IRSA HATS, version date 2026-05-13, s3://ipac-irsa-ztf/ztf/enhanced/
  dr24/) and the non-detection handling. Reused: `exotic_sim.inject_light_curve` (D-047), the W3 vetting ideas
  (bracketed feature, ≥ 3 nights, shared-epoch test, smooth/gradual dips, repeats).

## 2. Model (one)
- **W3: negative-mass point lens, Kitamura–Nakajima–Asada n = 1, ε < 0** (`exotic_sim`, sign = −1; Safonova et al.
  umbra; D-047). The lensed flux of a point source is zero for impact u < 2 (umbra), with caustic spikes at u = 2 and
  A = (u² − 2)/(u√(u² − 4)) beyond; a uniform-disk source of radius ρ caps the spikes (× 7.0 at ρ = 0.01, × 2.35 at
  ρ = 0.1; finite, never the infinite point-source spike). Umbra duration 2 t_E √(4 − u0²).
- Why this one: the largest photometric signature of the D-047 family (the source disappears), which ZTF's ~1–3 %
  photometry and multi-day cadence can see; Ellis (n = 2) gives ≤ 4 % dips or PSPL-like bumps (W4, not searchable).
- Scope: a **restricted hypothesis about a negative-mass lens**, not "wormholes" in general.

## 3. Data and what it can and cannot show
- ZTF DR light curves contain **detections only** (PSF-catalogue matches); non-detections are absent. A missing epoch
  can be: the source fainter than the limit, a masked/flagged/edge position, or no exposure. So:
  - **coverage**: for each (field, rcid, filter), an exposure epoch is *good* if ≥ 50 % of comparison stars (same
    quadrant, reference mag 14.0–18.0, all objects of the partition, not only the sample) have a detection with
    catflags = 0 within 0.0005 d;
  - **non-detection (ND)**: a good epoch where the target has no detection (it is then fainter than that epoch's
    limit; no flux is known); the per-epoch limit is estimated as the 95th percentile of comparison-star magnitudes;
  - detections with catflags ≠ 0 are dropped (cloud/moon, bit 15, and others).
- So a full vanish is visible only as ND on good epochs; a **partial dimming** (blended source) only while above the
  limit. Efficiencies below are for this stage on (light curves + this coverage rule), **not** for the upstream
  catalogue matching.

## 4. Sample (metadata only; `sample.ecsv`, seed 7)
- Light-curve HATS partition Norder 5 / Npix 3839 (centre RA 272.65°, Dec +64.95°, b ≈ +29°; low extinction,
  uncrowded) and the objects partition Norder 4 / Npix 959 restricted to that area.
- Stars: r-band (fid 2) objects with 14.5 ≤ refmag ≤ 17.5, |refsharp| ≤ 0.3, transient = 0, ngoodobs ≥ 300; isolated:
  no other object (any filter) within 5″ brighter than refmag + 2.5. **Groups**: every objectid (any field, filter)
  within 1.5″ of the chosen r-band object belongs to its group; groups are disjoint (greedy by ngoodobs).
- At most 1,000 groups drawn at random (seed 7) from those passing; no variability statistic is used.
- Split by group: holdout if crc32(str(group key oid)) % 2 == 0, else development. All filters/fields of a star stay
  in one half.

## 5. Pre-screen (the search statistic)
- Per group and filter, a light curve of detections (mag, err) and NDs over good epochs, all fields merged by time.
- **Dim state D** at an epoch: ND, or a detection fainter than the median by ≥ 1.0 mag.
- **Flag**: a run of ≥ 3 consecutive good epochs in state D spanning ≥ 3 distinct nights, bracketed by ≥ 20 baseline
  (non-D) epochs on each side in the same filter.

## 6. Vetting (in order; each a function, tested)
1. Common mode: at the run's epochs, < 20 % of the quadrant's comparison stars are also in state D.
2. Sharp edges: ≤ 1 detection in the intermediate range (median + 0.3 to + 1.0 mag) within one run-length before and
   after the run (dust/R CrB fades and long eclipses recover gradually; W3 ingress/egress last ~ρ t_E).
3. One-off: no second D-run anywhere else in the same light curve (eclipses, dippers repeat; a single lens does not).
4. Achromatic coverage: if the other filter has good epochs inside the run, they are in state D too.
- Survivors are looked at (light curves, neighbours within 2′, flags); at most 5 detailed checks.

## 7. Injection tests (development half only; seeds 11–14)
- Carriers: development-half groups, r band. W3 injected with `exotic_sim.inject_light_curve` (n = 1, sign = −1,
  blend = 1 and 0.7), t_E ∈ {3, 10, 30, 100} d, u0 ~ U[0, 2), ρ ∈ {0.01, 0.1}, t0 uniform inside the light curve;
  a model epoch fainter than its epoch limit becomes ND; otherwise mag = −2.5 log10(flux/baseline) + carrier mag +
  noise from the carrier's magerr–mag relation.
- Ordinary controls injected the same way: (a) periodic eclipses (depth 1.5 mag, duration 0.1–2 d, period 1–50 d);
  (b) dust/R CrB-like fades (decline 3 mag in 5–20 d, exponential recovery 30–200 d); (c) ordinary positive-mass
  PSPL (u0 ~ U[0,1], t_E 3–100 d); (d) common-mode: the same dim applied to the target and 50 % of its quadrant's
  comparison stars at one random epoch run of 5 epochs.
- 200 injections per W3 cell, 200 per control class.
- **Pass (all on development):** W3 flag+vet recovery ≥ 0.5 in at least the t_E = 30 and 100 d cells for both ρ
  (blend 1); each control class passes vetting in ≤ 10 % of injections; real development light curves leave ≤ 5
  survivors after vetting. If any fails: no holdout search; the failure is reported.
- If passed: the pipeline is frozen and run once on the holdout half; ≤ 5 survivors checked in detail.

## 8. Limits
- A pass is feasibility for this field and stage only; no W3 abundance or exclusion is computed (sample too small).
- Cap: data reading ≤ 30 min; stop computation at 14:30Z.

## Deviation (before any search or injection result; after the data-condition check)
- Per-epoch limit: the 95th percentile of comparison-star magnitudes would be ≤ 18 by construction (they are chosen at
  14–18), so the limit is instead the magnitude where the comparison stars' magerr–mag relation at that epoch reaches
  0.2171 (S/N 5), clipped to 18–22.5 (20.5 if fewer than 10 stars); it is used only by the injections.
- Data-condition check (comparison stars only, 3,000 light curves): no NaN magnitudes (detections only); catflags 0
  in 94 %, 32768 (cloud/moon) in 6 %; median spacing between epochs 1.0–4.0 d (10th–90th percentile); fields 826,
  848, 1856, 1857 (not stored beyond this note).
