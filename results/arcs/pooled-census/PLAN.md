# `arcs`: clump-pair census of the pooled anti rows (plan, 2026-11-08, before the run)

**Question.** The pooled anti count over the S1063 replication fields (`../pooled-tally/`, variant B: 6 anti vs
2.27 expected, P 0.029) is not a deciding-model error (`../anti-model-swap/`: `mixed`, 4 / 6 persist). Is it a
clump-pair effect, as in MACS1423 (`../macs1423-canucs/vetting/`: explained)? A two-knot source's long axis follows
the knot separation, which the shear only partly aligns, so the shape-noise null under-predicts anti for such rows.

**Hypothesis (ordinary).** The 6 variant-B anti rows are mostly two-knot sources; under the clump-aware null of #273
the pooled tail P rises to ≥ 0.05. **Alternative.** They are single-component sources, P stays < 0.05, and the excess
needs another ordinary explanation (photo-z contamination next, then σ_ε of these catalogues).

**Code.** `scripts/arcs_pooled_census.py`. It reuses main's `arcs_vet_census.py` unchanged (`deblend_row`,
`row_p_anti`, `clopper_pearson_upper`, `outcome`, `vc.clump_anti_fraction`, `vn.excess_p`) and
`arcs_pooled.background` for the row selection. No new rule.

**Rows.** The strong-shear background rows of each variant-B field's deciding run (the `FIELDS` of
`scripts/arcs_pooled.py` without MACS1423): MACS0417 49, WHL0137 18, RXC J0600 12, MACS0308 21, RXC J0018 4 (104 rows;
anti 3 + 1 + 1 + 1 + 0 = 6).

**Images (observed).** Each field's F200W level-3 mosaic that its screen catalogue came from, read as 8″ S3
byte-range cutouts (SCI, ERR, WHT); nothing is stored:
`jw01208-o006_t003` (MACS0417), `jw02282-o010_t001` (WHL0137), `jw06796-o010_t001` (RXC J0600),
`jw05293-o005_t001` (MACS0308), `jw06882-o001_t001` (RXC J0018), each `_nircam_clear-f200w_i2d.fits`.

**Rules (frozen, #273's; all ASSUMPTIONs).** Deblend `two-knot` if ΔBIC ≥ 10, both knot S/N ≥ 3, separation ≥ 0.10″;
own-flux gate with 20 + 20 injections, ≥ 80 % correct each; fit radius 0.4″, 1″ median filter; seed 20261108 + global
row index. Per-row P_anti: P_cp(g_min) for gate-passed `two-knot`, P_sn(g_min; σ_ε 0.26, e_min 0.5, 60°) for gate-passed
others, the f̂ mixture for gate-failed rows. g_min rounded to 2 decimals as `arcs_pooled.field_stats` (so the
shape-noise sum reproduces the pooled tally's 2.273; checked, the script stops otherwise).
- **f̂ (choice made here):** pooled over all gate-passed rows of the 5 fields (one null for one pooled test).
  Per-field f̂ is reported as a sensitivity, not used.
- **Primary:** Poisson-binomial P(≥ 6) over the 104 rows. **Secondary (reported, not used):** P(≥ 4) for the 4 rows
  that persist under the second models; variant A adds MACS1423's stored `census.ecsv` `p_anti_row` (55 rows, k 4) for
  P(≥ 10); the f̂ and 95 % upper mixtures.

**Outcome (frozen).** Fewer than 12 gate-passed rows, or a gate pass fraction below 0.5 → `undetermined`; otherwise
#273's `outcome(P_primary)`: P ≥ 0.05 → `anti excess explained by clump orientation (exploratory)`; P < 0.05 →
`excess remains under the clump-aware null`.

**What changes the plan.** Explained → the pooled `arcs` excess (variant B) is closed as clump orientation; S1063's
1316 / 2184 stay as they are; next on the line is nothing in the replication fields. Excess remains → photo-z
contamination of the anti rows (DJA / CANUCS z_phot posteriors vs the g_min interval) runs next.

**Look before concluding.** A contact sheet of the 6 anti rows (cutout with the fitted knots) is made and looked at.

**Labels.** Pixels observed; fits and classes derived; injections simulated; P values model_prediction. Exploratory:
the census idea is post hoc on S1063, and the pooled counts were seen before this plan.
