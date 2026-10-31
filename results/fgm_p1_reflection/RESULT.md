# fgm_p1_reflection — lower bound on P1's own reflected energy into mouth B

The canonical record (stage table, closure functional, setting S, vortex completions, early tail, mean field
versus single run, verdict) is
[`docs/results/fgm_p1_reflection_bound.md`](../../docs/results/fgm_p1_reflection_bound.md) (D-113 addendum,
2026-10-31). Everything here is analytic arithmetic and small integrals on the stored #619 numbers
(`results/fgm_window/window.json`) [der]; the setting S and the vortex transfer are labelled assumptions and
estimates in the record. Nothing is a simulation or an observation.

| File | Content | Produced by |
|---|---|---|
| `bound.json` | inputs read from `window.json`; reproduction of the critical fraction (1.90 × 10⁻⁷⁹); FGM angular weights (uniform, equatorial and minimal deposits); setting-S lower bound and break-even mean delay; single-run probabilities and the relative fluctuation; vortex transfer estimates per \|m\|; early-tail coefficient, crossing leads and lifetime extrapolation; narrow-band erratum value | `python scripts/fgm_p1_reflection.py` |

Verdict (owner's rule): mean field, P1 as frozen: **C** (the two-mouth transfer lacks P1's angular content about
the A–B axis, its spin and its preparation time); mean field in the added setting S: lower bound ≥ 7.2 × 10⁴⁶ times
the threshold; single run: P1 passes this necessary condition. P2's labels (frozen-rule A, post-hoc B) are unchanged.
