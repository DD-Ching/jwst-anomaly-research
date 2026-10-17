# W1 shear stacked at `radial` centres: plan (2026-10-17, frozen before any run)

**Question.** Do the `radial` screen's line-convergence centres (exotic W1 signature, D-031/D-049) carry the
radial (negative tangential) shear a W1 lens would also produce? The two screens use different observables
(arc orientations vs weak-lensing ellipticities of other sources), so agreement would be a cross-check; a null is
expected if the centres are chance line alignments (ordinary explanation).

**Inputs (as in `results/w1_shear_more` and `results/radial_more`; no streak veto, ASSUMPTION).** Fields with a
`radial` result whose models pass the `arcs` control: Abell S1063 (`results/radial_more/abells1063_reproduce`),
MACS0417, MACS1423, WHL0137. Centres: every `radial` peak (`peaks`, n_lines ≥ 3), pooled.

**Statistic.** S_c (and S_×,c as a control) from `exotic_screens.ApertureMass` with the shear-screen defaults
(10″ Schirmer, r_min 1″, min 5 sources), evaluated at the centres' model-frame positions, with cluster shear removed
as in `shear_screen`. Stack: mean S over all valid centres, pooled over fields.

**Null.** 2000 random rotations of every used source ellipticity (positions and |e| kept), the same rotation
for the centres of one field; p = fraction of draws with stacked mean ≥ observed (one-sided, W1 sign).

**Sign check (known case).** The same estimator with the model's cluster shear **kept** (not removed) at the
model-frame origin of each field must give S strongly negative (tangential shear of the cluster).

**Decision (ASSUMPTION thresholds).** p < 0.01 → the centres get a vetting step (sheet, sources, streaks)
before anything else; 0.01 ≤ p < 0.05 → exploratory note only; p ≥ 0.05 → null: the two exotic screens do not
reinforce each other at the current depth; record the stacked mean and its null 95th percentile as the reach.
