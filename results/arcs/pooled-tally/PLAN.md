# `arcs`: pooled anti count over the independent replication fields (exploratory)

**Status of this plan.** Written 2026-10-20 before running the pooled computation, but **after** the per-field
anti counts and expectations were read in the CHANGELOG handoffs, and a rough sum was formed while choosing the
action. The result is therefore labelled **exploratory (post hoc)**; it cannot count as a replication.

**Question.** Each replication field of the Abell S1063 anti excess "does not replicate" on its own, mostly
because single fields have little power. Do the fields together show more `anti` background sources than the
Poisson shape-noise expectation (`anti_expected_shape_noise_background`, σ_ε 0.26)? If yes, the likelier reading is
an ordinary one: the shape-noise null under-predicts `anti` in every field (clump pairs, deblending, photo-z blends,
model error), which would also weaken S1063 as a candidate. If no, the S1063 excess stays field-specific.

**Fields (frozen):** every field never used to tune `arcs` whose deciding run passed the background-alignment
control, with the deciding model recorded in its own RESULT.md: MACS0417 (CANUCS v1), MACS1423 (CANUCS DR1),
WHL0137 (RELICS Lenstool), RXC J0600 (RELICS GLAFIC v3), MACS0308 (RELICS GLAFIC v2), RXC J0018 (RELICS Lenstool).
"No result" fields (MACS0647, RXJ2129, SPT0615) are excluded because their control failed. The tuning fields and
S1063 itself are excluded.

**Rule (ASSUMPTION thresholds).** N = Σ n_anti, λ = Σ expected; P = Poisson P(≥ N | λ).
- Variant A: all six fields at screen level. Variant B: without MACS1423, whose excess the clump-pair census
  already explains.
- **Pooled excess** if P < 0.05 in **both** variants; **no pooled excess** if P ≥ 0.05 in both; otherwise
  **mixed**.
- Secondary (does not change the outcome): the σ_ε at which variant B reaches P = 0.05, and the split of pooled
  anti rows by g_min < 1 / ≥ 1 (critical-curve flips).
