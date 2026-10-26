# g85 host redshift from DJA catalogues: result (2026-10-26)

Plan `PLAN.md` (pushed in b1706213 before any row was read). Tool `scripts/macs0647_g85_host.py`,
rules tested in `tests/test_macs0647_g85_host.py`, both pushed in 9caaa1b7 before the run. Numbers: `summary.json`. Catalogue
values are `observed` (DJA v7.0 EAZY photo-z is itself a model fit, `derived` by DJA); separations `derived`.

## Result (frozen rules)
- **Outcome `no_host`.** No DJA v7.0 zout row (7,956 rows, member sha256 matched the pin) lies within 0.5″ of g85
  (H1), so H2 gives no class. DJA NIRSpec v4.6 has **no row at all** within 1.0″ (H3; gz sha256 matched the pin).
- The only zout row within 1.5″ is id 5954 at 1.079″: z_phot 0.636, z500 0.616 (16–84 % 0.568–0.656), no z_spec.
  Under H2 it would read `member` (|0.616 − 0.591| = 0.025 ≤ 0.159, not stored), but by the plan it is a neighbour, not the
  host, so the verdict does not use it.

## Reading (post hoc, not a gate)
- If id 5954 is the host, g85 sits about 1.1″ from its catalogue centre, which is about 7.3 kpc at z 0.62 (not
  stored; flat ΛCDM H0 70, Ωm 0.3). A supernova in the outskirts of a cluster member is then the simplest
  ordinary reading: at z 0.616 the vetting's 26.28 AB (no aperture correction, μ = 1, no K-correction) is
  M ≈ −16.0 (not stored), faint for a supernova near peak, so a late or a low-luminosity phase is required. A
  cluster member at z ≈ 0.6 is not magnified by the cluster, so this reading carries no lens-model leverage.
- The 0.5″ host radius may be too tight for a supernova offset from its host; that was the frozen ASSUMPTION, and
  the plan named no wider rule. The NED association in the vetting record ([CDK2017] 07064 at 0.02′ = 1.2″) is
  consistent with the same neighbour (not cross-checked by id).
- The g85 contact sheet (`../g85_sheet.png`, 1.6″ stamps) does not reach 1.1″ from the source in every panel, so
  the neighbour was not inspected by eye in this run.

## Limits
- The DJA-to-pipeline frame offset was not measured (the plan assumes ≲ 0.1″, as the Meena / Nishida comparison).
- One photo-z catalogue (DJA v7.0, 2023-07); no newer DJA MACS0647 photo-z was checked.
- Data: zout member (7.7 MB) and the v4.6 csv (445 MB unpacked; sizes not stored) streamed into the scratchpad and deleted.

## Next
Unchanged for the owner: a credentialed TNS cone search. A host spectroscopic redshift for id 5954 would decide the
member reading; no public NIRSpec spectrum within 1″ of g85 exists in DJA v4.6.
