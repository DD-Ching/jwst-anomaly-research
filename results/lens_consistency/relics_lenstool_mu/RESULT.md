# RELICS Lenstool v1 amplitude over the high-magnification region (H_K re-test; 2026-10-25)

Frozen plan `PLAN.md` (09620d66; gate T added in c3d60577, both before any map was read). Numbers come from
`summary.json`, written by `scripts/relics_lenstool_mu.py`. Maps are observed (published deflection and κ fields,
sha256-checked against `jwst_anomaly.lensmodel`); μ, R, classes and gates are derived. No imaging was read.

## Frozen result
**Gate T passes** for all eight references (median |κ_derived − κ_pub| 2e-6 to 1.1e-5 ≤ 0.05; the Jacobian and axis
convention reproduce each published κ). **Gate V passes.** **H_K is undetermined** (plan: an intermediate class).

| Cluster | Role | R (median, \|μ_ref\| ≥ 10) | 16–84 % | \|S\| arcsec² (fraction) | R at μ 5 / 20 | Class | f_obs | 1/R |
|---|---|---|---|---|---|---|---|---|
| MACS J0553 | validation | 1.635 | 1.23–2.18 | 2529 (0.046) | 1.651 / 1.631 | offset | 0.665 | 0.612 |
| Abell 3192 | validation | 1.277 | 1.06–1.42 | 3140 (0.083) | 1.276 / 1.296 | offset | 0.692 | 0.783 |
| RXC J2211 | validation | 1.070 | 0.95–1.15 | 3557 (0.095) | 1.077 / 1.069 | near 1 | 1.017 | 0.934 |
| RX J0152 | validation | 1.059 | 0.82–1.32 | 319 (0.008) | 1.079 / 1.064 | near 1 | 1.045 | 0.944 |
| RXC J0600 | test, batch | 1.384 | 0.98–2.32 | 7561 (0.201) | 1.312 / 1.366 | **offset** | — | — |
| MACS J0417 | test, batch (JWST-era CANUCS ref) | 1.375 | 0.85–1.64 | 2941 (0.078) | 1.336 / 1.397 | **offset** | — | — |
| WHL0137 | test, other maker | 1.038 | 1.00–1.08 | 1533 (0.118) | 1.038 / 1.038 | **near 1** | — | — |
| SPT0615 | test, other maker | 0.772 | 0.68–1.00 | 2407 (0.064) | 0.744 / 0.781 | **intermediate** | — | — |

Gate V: every validation cluster in its class, |1/R − f| 0.053 / 0.091 / 0.082 / 0.101 (≤ 0.12), S ≤ 0.10 of each
overlap (≤ 0.5). Reading: both batch maps are offset and WHL0137 is near 1, as H_K predicts; SPT0615's Lenstool map is
**smaller** than GLAFIC (R 0.77), neither near 1 nor offset, so the frozen outcome is undetermined, not supported.
Not rejected: no batch map is near 1 and no other-maker map is offset.

## Post hoc (not in the plan)
- **Bias.** 1/R is below f_obs in all four validation clusters (by 0.05–0.10): the median |α| ratio over S is not the
  image-fitted amplitude factor; the tolerance absorbed it.
- **Contact sheet** (looked at: ratio map with the |μ_ref| ≥ 10 contour, 8 panels; not committed, regenerable from
  the script's functions). S follows the critical curves. Inside S the ratio is not uniform for RXC J0600 (one side
  ≈ 2, the other ≈ 0.8) and MACS J0417 (a gradient along the major axis): there the two models differ in shape or
  centre, not by one factor, so "offset" for these two is a profile difference as much as a normalisation (16–84 %
  0.98–2.32 and 0.85–1.64). MACS J0553 is the cleanest offset.
- **SPT0615** (z_lens 0.972, D_LS/D_S(z 2) 0.375): an amplitude difference in the other direction. Not tested against
  images; no flag.

## Reading
- The batch association stays an association: three of the four test predictions hold, one is intermediate. Per the
  plan, the batch line is paused (no third region variant). The addendum-bd per-cluster amplitude check stays the only
  rule before a RELICS Lenstool v1 map is used; RXC J0600 and MACS J0417 Lenstool v1 are now known to differ from a
  second model by ≈ 1.38 in |α| over S and should not be used without that check.
- Hosts: archive.stsci.edu with the proxy + certifi bundle in `SSL_CERT_FILE`; ≈ 1 GB streamed to memory twice (run
  and contact sheet), nothing stored.
