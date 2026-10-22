# Abell 370 gate C2 failure: post hoc diagnosis (plan)

Written before any computation in this directory. Everything here is **post hoc**: the frozen outcome of
`../RESULT.md` (C2 fails, arc-referenced 19 / 25, no flag read) stands whatever this finds.

- **Question:** why are 6 of the 25 testable arc-referenced controls `absent`? Is the failure one cause (one physical
  source counted several times) or several independent ones?
- **Observation that prompted it:** controls 2 (2.4+2.5), 102 (102.3+102.4) and 202 (202.4+202.5) all take their
  reference flux from the same CANUCS source 2101809 and land on the same counterpart 2101806.
- **Ordinary explanations to separate:** (a) knots of one source referenced to a whole-arc segment (reference flux
  too large for a knot); (b) reference segment much larger than the control's own segment (15.3: a3 1.86″ vs 0.17″);
  (c) real faintness / local noise; (d) a wrong model magnification ratio.
- **Method:** recompute C2 with the repository's `gate_c2` on the stored `control_abell370.ecsv` (must reproduce
  19 / 25), then the same with controls that share one (reference, counterpart) source pair across systems counted
  once (absent if any is absent). Same count on the four other CANUCS-reference fields.
- **What would change the plan:** if the shared-pair counting changes no class on the four other fields and lifts
  Abell 370 to ≥ 0.80, propose it as a rule to **freeze before** any new field (it cannot be validated on Abell 370
  itself, and no flag is read here). If it changes other fields' outcomes, it is rejected.
