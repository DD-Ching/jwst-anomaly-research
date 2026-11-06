# SMACS 0723 jwst 3.0.0: galaxy-stratum re-run without SIMBAD — plan (written before the run)

Follow-up to `PLAN.md` / `RESULT.md` (#742), whose ranks were a single mixed stratum because CDS was unreachable.

- **Deviation (stated before the run):** CDS (SIMBAD TAP, XMatch) still resets connections from the cloud session;
  the ESA Gaia archive TAP answers. Star separation (D-012) therefore runs with `classify.services: [gaia]` only:
  astrometric Gaia stars plus any Gaia source within 0.3″. The SIMBAD rule "a nearer SIMBAD galaxy wins" cannot
  apply, so some compact cluster members may be labelled stars. Bound: the count of `gaia_position`-only stars.
- **Hypothesis:** the six surviving SMACS records keep galaxy-stratum ranks within ±3 of their recorded ones.
- **New entries:** galaxy top-20 sources absent from the 2.0.1 records are listed with the contact sheet seen; none
  is called a candidate before `/vet-candidate`.
- **What changes the plan:** a recorded source falling out of the galaxy top 50, or > 5 `gaia_position` stars among
  the top 20 neighbours (ASSUMPTIONs), makes the result undetermined until SIMBAD is reachable.
