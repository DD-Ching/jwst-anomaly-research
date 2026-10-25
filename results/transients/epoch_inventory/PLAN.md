# Two-epoch inventory of public JWST NIRCam imaging on the D-027 cluster fields (plan, frozen before any query)

**Question.** D-027's two-epoch line says "next only a new epoch pair" for every field it has tested. Does MAST now
hold public NIRCam imaging epochs on those fields that no tested pair uses?

**Method (frozen).** One batched `astroquery.mast` `Observations.query_criteria` per field centre (radius 3′,
`obs_collection=JWST`, `instrument_name=NIRCAM/IMAGE`, `calib_level=3`, `dataRights=PUBLIC`), metadata only, no
product download. An epoch = one (proposal, observation) with `t_min`; filters listed per epoch.

**Rule (frozen).** A pair is "new" when (a) both epochs are public today, (b) they share at least one filter,
(c) Δt ≥ 30 d, and (d) at least one of them is not an epoch already used by a tested pair in CHANGELOG / TASKS
(D-027 addenda). Output: a table of new pairs per field, ranked by the latest epoch's `t_min`.

**What changes the plan.** ≥ 1 new pair on a field with a validated pipeline (`transient_grid.py`, `arc_diff.py`)
→ queue that pair as the next transient run. 0 new pairs → the line stays waiting, with the next public release
dates (`t_obs_release`) recorded.

**Ordinary explanations are not tested here** (metadata only; no science reading).
