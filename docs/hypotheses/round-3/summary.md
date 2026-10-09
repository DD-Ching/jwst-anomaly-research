# Round 3 (2026-10-09): derived small effects in large-N data

Everything here is a **hypothesis** or a reviewer estimate ([B-est], [calc], [my derivation]). System A wrote
without literature or repository access. The System B reviews cite references their authors opened; the lists are
at the end of each B file.

## Brief

Round 3 kept the round-1 rules (D-069) and added two lessons from round 2 (D-075):
- Find a small effect whose size is **derived** from the axioms.
- It must be visible only statistically, in a **large-N public quantity** that has never been analysed for it.

| Group | Seed | World | File | Review |
|---|---|---|---|---|
| R3-A | Correlation between records is the primitive | Universal random pair correlation of innovations in the same exposure, about 2.4 × 10⁻⁵ | [R3A-correlation.md](R3A-correlation.md) | [B-on-R3A.md](B-on-R3A.md) |
| R3-B | The owner's "spacetime response material" | "Phase-Locked Dilaton": a scalar charge carried only by condensates (superfluid neutron-star interiors) | [R3B-response.md](R3B-response.md) | [B-on-R3B.md](B-on-R3B.md) |
| R3-C | The owner's "transient connectivity" | Runaways need a partner with the same "clock depth" λ; the SN Ia rate falls in deeper wells | [R3C-transient.md](R3C-transient.md) | [B-on-R3C.md](B-on-R3C.md) |

## Verdicts

| World | Main prediction | B verdict | Main reason |
|---|---|---|---|
| R3-A | ZTF split-sample pair statistic Ŝ = 1/D ≈ 6 × 10⁻¹⁰ | **FAIL as written**. Conditional only for a repaired axiom; the cheapest test is a lab test outside scope | Axiom A5 gives zero. The Lorentzian signature is not derived. D is chosen (it spans 10⁸–10¹⁵). The Fermilab Holometer cross-correlation already disfavours the nominal amplitude at about 67 % (reviewer's own normalisation). ZTF would need about 3 TB of downloads and reaches about 1.8σ at nominal D, and every systematic biases Ŝ upward |
| R3-B | NS–WD dipole radiation with α_NS ≈ 2–9 × 10⁻³ | **FAIL (excluded)** | J1738+0333 alone excludes s ≥ 8 × 10⁻³ at 6.7σ (7.1σ in preprint 2605.01436). The combined bound is α_NS < 2.2 × 10⁻³ (95 %, reviewer's combination); the self-consistent minimum is disfavoured at about 2σ and the no-feedback corner survives. The model is Brans–Dicke with ω = 0 plus scalarisation |
| R3-C | SN Ia rate δ ln Γ = K Δλ: −6 % (10¹⁴ M☉ groups) to −20 % (cluster cores) | Not excluded; conditional PASS for P1/P3, but **untestable as posed** | λ has no Φ term (my derivation: λ = 1 − 1/γ). K ≈ −2.4 to −4 × 10³ falls at the weak edge, where its kill cannot fire. The published cluster/field SN Ia ratios (×0.6 to ×3.5; ZTF excesses ×3–8) scatter 10–40× the signal. The ZTF BTS Ia sample at z < 0.1 is about 5 × 10³, which gives 0.5–2σ, and the systematics floor exceeds the signal at any N |

## What this project can test

Nothing in round 3 earns search compute:
- **R3-A:** needs a lab test of independent noise sources, about 400 s per device pair.
- **R3-B:** already excluded. Its bound is recorded: α_NS < 2.2 × 10⁻³ (95 %).
- **R3-C:** the public sample is about 16× too small, and the environment systematics exceed the signal.

## Lessons (rules)

- **Derived amplitudes collide with existing precision.** Once forced to derive an amplitude, each world either
  predicted something already excluded (R3-B), something that needs a lab test (R3-A), or something buried under
  astrophysical systematics 10–40× larger (R3-C).
- **Round 1 kept bound-setting tests; rounds 2–3 kept nothing.** Three rounds of 9 worlds produced 3 conditional
  survivors, all in round 1, and all were null in data (S1–S3, D-070–D-072). Under the D-069 rules the
  invention/destruction loop is now returning nulls before any data analysis.
- **Recommendation:** pause new System A rounds. Put compute into the searches that are still open:
  - finish MOA-II (W3, 22 fields);
  - extend the round-1 tests that are limited by data, not ideas: CHIME/FRB baseband, Euclid DR1 for W5, Gaia DR4
    when released.

  Restart System A only with a new brief, for example from the owner.
