# Clump-pair census on five other `arcs` fields (plan, written before any run)

**Question.** The Abell S1063 census (#273, `results/arcs/abells1063-cats/vetting/RESULT.md` step 4) estimated the
fraction of compact two-knot sources among strong-shear background rows from that field alone (f̂ = 5 / 21 = 0.238,
95 % upper 0.472). Its clump-aware null left P(≥ 4 anti) = 0.026, but P ≥ 0.05 at the upper f. Does an f measured on
independent fields keep or remove the S1063 excess?

**Hypotheses.** H0: two-knot sources are as common elsewhere as in S1063 (pooled f ≈ 0.24): the S1063 excess stays
at P ≈ 0.02–0.03 (exploratory). H1: they are more common (pooled f ≳ 0.45): the clump-aware null explains it.
Ordinary explanations this does not test: blends of unrelated galaxies, κ ≈ 1 rows (2004), field-to-field depth.

**Data.** The five stored `results/arcs/<model>/arcs.ecsv` tables (rows with `strong_shear` and
`population == background`), and 8″ F200W S3 byte-range cutouts of the i2d mosaic each table's catalogue came from
(`arcs.json` "catalog"). No download of whole mosaics.

**Method (frozen; unchanged S1063 code).** For each row: `arcs_vet_deblend.fit_models` + `deblend_class` at fit
radius 0.4″, own-flux gate with 20 injections per kind (`arcs_vet_census.deblend_row`, seed 20261016 + row index).
Per field and pooled over the five fields: n gate-passed, n `two-knot`, f̂ and its Clopper–Pearson 95 % upper.

**Use in S1063 (frozen).** Re-evaluate the S1063 mixture null P(≥ 4 anti) with the pooled f̂ and its 95 % upper
(the S1063 per-row p_cp and p_sn from `census.json`). Outcome rule:
- pooled-f̂ mixture P ≥ 0.05 → "S1063 excess explained by clump pairs at the five-field rate (exploratory)";
- otherwise → "S1063 excess remains at the five-field clump rate (exploratory)".

**Validity checks.** Fields with < 12 gate-passed rows are reported but not interpreted alone. Rows whose cutout
lies outside the mosaic are dropped and counted. Contact sheets of every `two-knot` row are looked at before
writing a conclusion.

**Labels.** Pixels observed; classes derived; injections simulations; f̂, P model_prediction (ASSUMPTIONs as the
S1063 census: point clumps, isotropic separations).
