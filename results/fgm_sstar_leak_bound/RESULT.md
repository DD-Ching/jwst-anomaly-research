# fgm_sstar_leak_bound — direct exclusion test for S★

The canonical record (thresholds, method, path table, bound, verdict, calibration versus fluctuations, physical
reading) is [`docs/results/fgm_sstar_leak_bound.md`](../../docs/results/fgm_sstar_leak_bound.md) (D-113 addendum,
2026-11-02). Everything here is derivation, RN null-geodesic quadrature and small integrals on the stored #619 point
and the Round 3 source S★ [der/computed]; allowances are labelled in the record. Nothing is a simulation of the full
wave problem or an observation.

| File | Content | Produced by |
|---|---|---|
| `bound.json` | reproduced thresholds (0.2554, 0.3125, K_critical 1.2235); RN photon sphere and paths from the aperture edge to B (impact parameters, travel times, orbit gate); reciprocal field at the edge (phase mismatch, ray-tube amplitude, edge-saddle factor and its envelope-peak checks on a simplified two-edge integral, family-1 ratio); the lower bound on K under the stated approximations, with allowances and break-even offset; aperture scan; isotropic comparison; verdict | `python scripts/fgm_sstar_leak_bound.py` |

Conclusion: Under the specified geometric-optics, aperture and mean-field approximations, the calculation supports
excluding the fixed source S★. A physical lower bound for the full two-mouth propagation remains to be verified.
Within these approximations K ≥ 6.5 × 10³⁶ against 1.224 (ideal phase reference, η ≤ 1). The isotropic estimate
(10⁵³·⁹) is kept as a comparison only.
