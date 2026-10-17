# W5 count deficit at the `radial` centres: plan (frozen before any run, 2026-10-17)

**Why.** D-047 pairs W5 (fewer background sources inside about θ_E of a negative-mass lens) with the `radial`
centres; TASKS JWST focus lists "W5: counts N(>S) around `radial` centres" as the open item. The W1 shear stack at
the same centres excludes ≥ 1 × 10¹² M☉ at every centre (`results/w1_stack_radial/inject/RESULT.md`); this test uses
a different observable (counts, not shapes) on the same centres.

**Hypothesis (H, hypothesis only).** Each `radial` centre holds an n = 1 negative-mass lens of mass M; inside
x = r/θ_E ≲ 1 images are demagnified, so background counts above the flux limit drop by
`exotic_sim.count_ratio` (`model_prediction`).

**Ordinary explanations of a deficit.** Masking by cluster members and ICL near the centres (the centres sit in
the cluster core), catalogue depth or edge changes, cluster magnification bias (depletion for steep counts), and
circularity (the arcs that define the centres are themselves counted). Of an excess: lensed arcs fragmented by the
deblender near the centres, galaxy clustering.

## Data and sample (as `scripts/w1_stack_radial.py`)
- Fields with `radial.json`: Abell S1063, MACS0417, MACS1423, WHL0137 (`w1_stack_radial.RADIAL`), loaded with
  `inject_radial.load_field` on `inject_shear.field_spec` (photo-z attached as for the shear screen).
- Counted rows: background (z_phot > z_l + 0.1 where a photo-z exists; rows without one are kept, ASSUMPTION), flux
  above the field's 20th-percentile background flux (completeness proxy, ASSUMPTION), and **not** a `radial` arc row
  (`w1_stack_radial.radial_arc_rows`; circularity removed up front).
- Count radii R = 3″, 6″, 9″ (≈ θ_E at 1 × 10¹², 4 × 10¹², 1 × 10¹³ M☉ for z_l 0.4–0.5, z_s 2).

## Statistic and null
- N_obs(R) = Σ over centres of counted rows within R. Controls: per centre, positions drawn uniformly in the annulus
  of cluster-centric radius r_c ± 10″ (model frame), kept when ≥ 1 counted row lies within 20″ (footprint proxy,
  ASSUMPTION) and the position is ≥ 2R from every centre. One control per centre per draw, 2000 draws, seed 2026.
- Ratio Q(R) = N_obs / mean control sum; one-sided p_low = share of control sums ≤ N_obs (deficit),
  p_high likewise for an excess.

## Forecast (power, before the measurement is read)
- Predicted Q_H(R, M) = area mean of `count_ratio` over r ≤ R for θ_E(M) at each field's z_l, z_s 2, with N(>S) the
  field's own empirical background counts (ASSUMPTION: field counts stand in for unlensed counts).
- Power(M, R) = share of control draws whose sum, scaled by Q_H, lies at or below the control 5th percentile.

## Decision rule (thresholds are ASSUMPTIONs)
- p_low < 0.01 at any R → vet the centres (ordinary explanations above), no interpretation.
- 0.01 ≤ p_low < 0.05 → exploratory note.
- Otherwise null; quote Q(R) with the control 5–95 % range, and the power table. If power < 0.5 at every
  M ≤ 1 × 10¹³ M☉ and R, record that counts at these centres cannot beat the shear stack and close the W5 item.
