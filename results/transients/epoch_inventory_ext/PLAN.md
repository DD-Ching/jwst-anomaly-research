# Two-epoch NIRCam inventory on 15 further lensing clusters (plan, frozen before any query)

**Question.** The 19-field inventory (`../epoch_inventory/`, #539) and the Abell 2744 footprint check left no
public two-epoch pair through a strong-lensing core that the D-027 runs have not used; the next known input is
Abell 370 from 2026-12-29. Do other well-modelled lensing clusters hold public NIRCam epoch pairs?

**Hypothesis / ordinary expectation.** Most clusters were imaged once (single GTO/GO/CANUCS-type visits); repeat
epochs exist mainly where a lensed transient or a time-delay quasar was followed (MACS J0138, SDSS J1004). Pairs
far from the core (parallels, flanking pointings) are expected and are not lens-model tests.

**Fields (frozen; MAST-resolvable names).** Lens-consistency fields of this repo not in the first inventory:
MACS J0138.0-2155, RX J0152.7-1357, MACS J0358.8-2955, RXC J2211.7-0350, SPT-CL J0546-5345, MACS J2129.4-0741.
Further well-studied lensing clusters: MACS J1206.2-0847, Abell 1689, MACS J0329.6-0211, Abell 383, Abell 209,
MACS J0308.9+2645, SDSS J1004+4112, SDSS J1029+2623, SDSS J2222+2745.

**Method (frozen).** `scripts/epoch_inventory.py --field-set ext`: the same batched query (3′ cone, JWST NIRCam
level-3 images), the same epoch definition and the same new-pair rule (a)–(d) as `../epoch_inventory/PLAN.md`.
Used epochs are read from `results/transients` excluding both inventory directories (they list every epoch, not
used ones). A field whose name MAST cannot resolve is recorded as `unresolved`, not dropped silently.
Core test (frozen, as in the Abell 2744 footprint check): for each new pair, the overlap of both epochs' `s_region`
footprints with a 1′ disc (ASSUMPTION) around the resolved field centre, per shared filter; `core` when > 1
arcmin² in some shared filter, `overlap` when the footprints overlap elsewhere, `none` otherwise.

**What changes the plan.** ≥ 1 `core` pair → queue it as the next D-027 grid (`transient_grid.py`, plan first).
0 `core` pairs → the two-epoch line waits on `upcoming.json` release dates; record them.

Metadata only (`observed`); no product download, no science reading.
