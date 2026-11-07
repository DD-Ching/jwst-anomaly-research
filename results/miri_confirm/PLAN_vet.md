# Plan: vet the five new top-20 F770W singles (2026-11-07)

Sources: `f770w_221`, `_239`, `_176`, `_118`, `_240` of run `20261107T022631Z-036f249e`
(`configs/smacs_jwst3_miri_confirm.yaml`). Procedure: `/vet-candidate`, methodology order A → B → C.

**Hypothesis (ordinary, to test first).** Each entry is ranked high because of MIRI processing or catalog effects
(local-background annulus inside a large isophotal footprint, deblending, bright-neighbour light) or because it is a
known population (dusty/PAH-bright galaxy, AGN, lensed arc in SMACS 0723).

**What would change the plan.** Any source that survives A–C ("unexplained") goes to a needs-human candidate PR.
If ≥ 3 / 5 are explained by the same processing effect (ASSUMPTION threshold), that effect becomes a quality-gate
proposal for the next cycle instead of further per-source vetting.
