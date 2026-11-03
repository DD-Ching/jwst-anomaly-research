# fgm_oneway_bit — one-way bit through the FGM window: necessary conditions

The canonical record (model, Helstrom benchmark, FGM differences, two path accounts, the two bounds, critical
inequality, plain answer) is [`docs/results/fgm_oneway_bit.md`](../../docs/results/fgm_oneway_bit.md) (D-113
addendum, 2026-11-02). Everything here is analytic arithmetic and small integrals on #619's stored FGM point
(`results/fgm_window/window.json`) [der]; the source S★ is a design choice and the halo-to-B transfer τ is an
estimate (an upper bound τ ≤ 1 and an isotropic value). Nothing is a simulation or an observation.

| File | Content | Produced by |
|---|---|---|
| `bound.json` | inputs; source S★; Helstrom and pilot thresholds; thermal and mixing factors; path timing, emerging frequencies, phase sensitivity; halo fraction, spectral and reflection leaks; FGM time and angle weights; the two intervals, critical τ and leak share; crossover carrier | `python scripts/fgm_oneway_bit.py` |

Verdict for S★: the intervals do not overlap (B needs N_A ≥ 0.2554; the closure limit allows N_A ≤ 10⁻⁵⁴·⁴ with the
isotropic τ). Missing coefficient: τ; a lower bound τ ≥ 10⁻⁶¹·⁴ settles it. A second missing input is an in-situ
phase reference (δ ln|ΔV| ≤ 6.9 × 10⁻¹⁶).
