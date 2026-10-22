# Draft inquiry (not posted; for the owner's review)

Status: draft only, 2026-10-22. Nothing has been posted, sent or filed. The repository link works only once
PR #452 is merged into main.

## Where to ask
- **Primary: the cogwheel issue tracker** (https://github.com/jroulet/cogwheel/issues).
  - The package metadata names this repository as its source.
  - The README names no other support channel (checked 2026-10-22).
  - The question is how to use cogwheel's waveform generator and refinement consistently with lal's XPHM.
- **Backup: the LALSuite service desk.**
  - This is the helpdesk address listed in LALSuite's `CONTRIBUTING.md`, under "More Information".
  - The address was read from that file on 2026-10-22; that e-mailing it opens an issue was not verified.
  - Use it only if the cogwheel maintainers point to lal's conventions.

## Draft (English)

**Title:** Setting up a consistent IMRPhenomXPHM comparison with cogwheel at aligned-spin, edge-on points

Hello, and thank you for cogwheel.

We use cogwheel 1.6.0 (lalsuite 7.26.16) to maximise a single-detector, zero-noise fitting factor for
IMRPhenomXPHM over a box with full precessing spins. The box allows q ≤ 15, |χ| ≤ 0.99, any tilt (cos tilt = ±1
included), any inclination ι and any φ_ref. Our own score calls `SimInspiralChooseFDWaveform` with multibanding
switched off.

**What agrees.**
- We compared h+ and h× at identical parameters, mapped with `inplane_spins_xy_to_xy_n`, at 37 parameter sets.
  The noise-weighted 1 − overlap is ≤ 7.3e-13, including a nearly aligned, edge-on set (cos tilts ±0.999, ι = π/2).
- At those sets, `lnlike_fft` with the ASD drift set to 1 reproduces our normalised overlap to ≤ 1.1e-12.

**Where our comparison breaks down.** The set is m1 36.343, m2 4.543 M☉; χ1 = (0, 0, +0.99), χ2 = (0, 0, −0.99);
ι = π/2; φ_ref = π; f_ref = f_min = 20 Hz.
- Our check of `WaveformGenerator.get_hplus_hcross` against lal's h+/h× gives 1 − overlap = 2.0 for at least one
  polarisation. Which polarisation it was is not recorded.
- At this ι, h× is expected to be small, so that number may be ill-defined.
- Our own mapping from the fitted (a h+ + c h×) to distance and ψ is also ill-conditioned here: it returns
  d ≈ 0.01 Mpc.
- We therefore do not know whether the waveforms differ at all, and we assume no cause.

**Questions.**
1. At such a point (exactly aligned spins, ι = π/2), what is the recommended way to compare cogwheel's h+/h× with
   `SimInspiralChooseFDWaveform`, and the supported parameter conventions for doing so (φ_ref, in-plane spin frame,
   ι)? Is any difference expected there by design?
2. For a full precessing search whose box includes this boundary, how would you set up
   `Posterior.refine_reference_waveform`, or another cogwheel route? In our tests the original point was not always
   among its unfolded images. The gap is large near ι = π/2 and below 1e-4 in lnL elsewhere. Is that expected?

**Reproduction.**
- Package: https://github.com/DD-Ching/jwst-anomaly-research/tree/main/results/x4/handoff, file
  `cogwheel_aligned_case.json`. It holds the parameters, both call sequences, the stored outputs, the versions and
  the command.
- Script: `scripts/x4_cog_exact.py recheck`.
- Outputs: `results/x4/cog_exact/review_recheck.json`, rows `tilts_pm1_chi099` and `tilts_pm0999`.

Thank you for any pointers.
