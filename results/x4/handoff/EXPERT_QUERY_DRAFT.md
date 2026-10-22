# Draft inquiry (not posted; for the owner's review)

Status: draft only, 2026-10-22. Nothing has been posted, sent or filed.

## Where to ask
- **Primary: the cogwheel issue tracker** (https://github.com/jroulet/cogwheel/issues).
  - The package metadata names this repository as its source.
  - The README names no other support channel (checked 2026-10-22).
  - The question concerns using cogwheel's waveform generator and refinement with lal's XPHM conventions.
- **Backup: the LALSuite service desk**, the helpdesk address in LALSuite's `CONTRIBUTING.md`
  (https://git.ligo.org/lscsoft/lalsuite/-/raw/master/CONTRIBUTING.md, "More Information"; e-mail creates an issue).
  Use it only if the cogwheel maintainers point to lal's XPHM conventions.

## Draft (English)

**Title:** IMRPhenomXPHM h+/h× from WaveformGenerator differ from SimInspiralChooseFDWaveform at one aligned-spin, edge-on point: supported way to get matching waveforms?

Hello, and thank you for cogwheel.

We are using cogwheel 1.6.0 (lalsuite 7.26.16) to maximise a single-detector, zero-noise fitting factor for IMRPhenomXPHM over a box that includes the full precessing parameters:
- q ≤ 15 and |χ| ≤ 0.99;
- any tilt, including cos tilt = ±1;
- any inclination and φ_ref.

Our score calls `SimInspiralChooseFDWaveform` directly, with multibanding switched off.

**What agrees.** We compared h+ and h× at 38 points with identical parameters, mapped with `inplane_spins_xy_to_xy_n`:
- the noise-weighted 1 − overlap is ≤ 7.3e-13;
- `lnlike_fft` (ASD drift set to 1) reproduces our normalised overlap to ≤ 1.1e-12.

**What does not.** At one legal point the two disagree. The point is:
- m1 36.343, m2 4.543 M☉;
- χ1 = (0, 0, +0.99), χ2 = (0, 0, −0.99);
- inclination π/2, φ_ref π, f_ref = f_min = 20 Hz.

There, `WaveformGenerator.get_hplus_hcross` and lal's h+/h× give 1 − overlap = 2.0 for at least one polarisation; which one is not recorded. The fitting factors at our optimum differ by up to 0.75.

We have not separated the roles of the exactly aligned spins, the edge-on inclination and φ_ref, and we do not assume a cause.

**Questions.**
1. At such a point, what is the supported parameter convention for obtaining h+/h× from cogwheel that match `SimInspiralChooseFDWaveform` (φ_ref, in-plane spin frame, inclination)? Or is a difference expected here by design?
2. For a full precessing search whose box includes this boundary (cos tilt = ±1, inclination π/2), how would you set up `Posterior.refine_reference_waveform`, or an alternative cogwheel route?
   - Should such points be handled differently?
   - Is the folding at θ_JN = π/2 a known limitation? We saw the original point not always among its unfolds there.

**Reproduction package** (parameters, both call sequences, stored outputs, versions and the command):
https://github.com/DD-Ching/jwst-anomaly-research/tree/main/results/x4/handoff, file `cogwheel_aligned_case.json`.
- Script: `scripts/x4_cog_exact.py recheck`.
- Stored outputs: `results/x4/cog_exact/review_recheck.json`, row `tilts_pm1_chi099`.

Thank you for any pointers.
