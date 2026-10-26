# g85 host redshift from DJA catalogues: plan (frozen 2026-10-26 before any catalogue row is read)

**Question.** What is the redshift of the galaxy under the F200W transient g85 (RA 101.99665, Dec +70.250721;
`docs/candidates/macs0647-1433-g85.md`)? A cluster member (z 0.591) and a background galaxy give different
luminosities and, for a background host, possible lensed repeats.

**Ordinary explanations held in view.** Supernova in a cluster member, a foreground or a background galaxy; AGN
variability in the host nucleus. Exotic readings are not tested here.

**Inputs (pinned in SOURCES.md, streamed, not stored).**
1. DJA v7.0 `macs0647-grizli-v7.0-fix.photoz.tar.gz` (41,494,640 B, sha256 `d56246b5…cdf0c`), member
   `macs0647-grizli-v7.0-fix.eazypy.zout.fits` (sha256 `ceb199c6…2c04e`). Columns used: `id`, `ra`, `dec`, `z_phot`,
   `z160`, `z500`, `z840`, `z_spec` if present.
2. DJA NIRSpec `dja_msaexp_emission_lines_v4.6.csv.gz` (sha256 `ee64bb46…c6ddf3`), grade 3 rows only.

**Rules (ASSUMPTIONs; implemented as tested functions in `scripts/macs0647_g85_host.py`).**
- H1 host: the nearest zout row within 0.5″ of g85 (the DJA-to-pipeline frame offset is not measured here; the
  Meena/Nishida lists agree to ~0.1″). Also report every zout row within 1.5″ (neighbours), not used for the verdict.
- H2 class from the host's z160 / z500 / z840 (z_cl = 0.591):
  - `unconstrained` if z840 − z160 > 1.5;
  - else `member` if |z500 − z_cl| ≤ 0.1 (1 + z_cl);
  - else `foreground` if z500 < z_cl, `background` if z500 > z_cl.
- H3 spec-z: any grade-3 NIRSpec row within 1.0″ of g85; it overrides H2 when present (class by the same
  member window).
- No host within 0.5″ → outcome `no_host`.

**What would change the plan.** `background` (or a spec-z behind the cluster): the next step is the CLASH
LTM-Gauss magnification and whether the host is multiply imaged (a lensed host predicts repeats). `member` /
`foreground`: g85 is a field / cluster SN candidate with no lens-model leverage; the line ends at TNS (owner).
`unconstrained` / `no_host`: the verdict stays as it is.

**Derived, reported but not a gate.** Absolute magnitude at the host z500, M ≈ m − DM(z) + 2.5 log10(1 + z)
(no K-correction, μ = 1, flat ΛCDM H0 70, Ωm 0.3; ASSUMPTIONs), with m = 26.28 AB (the vetting's 0.15″ value,
not aperture-corrected, so M is an upper bound on the brightness).
