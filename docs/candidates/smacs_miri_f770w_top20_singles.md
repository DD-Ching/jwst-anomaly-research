# SMACS 0723 MIRI: the five new top-20 F770W singles (vetting record)

## 2026-11-07: vetting (`results/miri_confirm/PLAN_vet.md`)

- **Run** `20261107T042021Z-8f4189cd`, config `configs/smacs_jwst3_miri_confirm.yaml`, commit `ad18cd6` (main after
  #756). It re-runs the #754 run `20261107T022631Z-036f249e` (outputs lost with that cloud session); the top 20 and
  the five ranks are identical. Sample `smacs0723_miri`, galaxy stratum. MIRI L3 catalogs are jwst 3.0.0.
- **Scores** (model_prediction): `_221` rank 2 (0.977), `_239` 4 (0.965), `_176` 13 (0.902), `_118` 16 (0.896),
  `_240` 18 (0.879). Four are led by F770W concentration (`ref_log_ci_50_30` +4.0, at the flat-profile ceiling
  2.409 per #755); `_176` by low sharpness (−3.2).
- **Cross-match** (observed): SIMBAD TAP (CfA mirror), NED TAP `objdir` and Gaia DR3 (ESA TAP), queried
  2026-11-07T04:30Z, radius 3″ (CDS XMatch still fails, #751). NED `z`, `zflag`, `z_bibcode` from a positional
  `objdir` query (1″) at 04:35Z. No Gaia DR3 source within 3″ of any of the five.
- **Cross-band association** (derived, `results/miri_confirm/vet_xband.json`): nearest source of each other MIRI
  band's L3 catalog to the F770W centroid. The run merges bands at 0.1″ (`stages.catalog.merge_radius_arcsec`).

| uid | F1000W | F1500W | F1800W | footprint (F770W px) |
|---|---|---|---|---|
| `_221` | 0.53″ | 0.66″ | 0.63″ | 897 |
| `_239` | 0.26″ | 0.40″ | 16.1″ | 1035 |
| `_240` | 0.26″ | 1.45″ | 16.7″ | 250 |
| `_118` | 0.18″ | 0.27″ | 15.6″ | 49 |
| `_176` | 5.2″ | 7.3″ | 9.2″ | 10 |

Contact sheets (`bands.png`, all four bands) were viewed: `_221`, `_239`, `_240` and `_118` are the same object in
F770W and F1000W (and F1500W for `_221`, `_239`, `_118`); `_176` is visible by eye in F1000W but not in its catalog.

### Tests

| Test | Tier | Outcome | Evidence | Label |
|---|---|---|---|---|
| Edge / low weight | A | does not explain (all five) | gate `quality_ok`, rel_weight 1.02–1.03, edge 18–39″; cutouts `ok` | derived |
| Bright star / spikes | A | does not explain (all five) | nearest classified star 11.9–24.2″; no spikes in the cutouts; no Gaia source within 3″ | observed |
| Single-band artifact | A | does not explain (all five) | each is a galaxy seen in ≥ 2 bands on the cutouts | observed |
| Cross-band association | B | **explains "single band"** for `_221`, `_239`, `_240`, `_118` | counterparts at 0.18–0.66″ (table), all beyond the 0.1″ merge radius; centroids of extended or lensed light differ between bands | derived |
| Background annulus inside footprint | B | likely contributes for `_221`, `_239` (hypothesis, #755) | equivalent radii 16.9 / 18.2 px > 8.9 px annulus inner radius; CI at the ceiling | hypothesis |
| Bad deblending | A | contributes for `_221`, `_240` | `_221` marks the diffuse band of a large galaxy 0.5″ from its knot (where the other bands centroid); `_240` is a knot in the BCG halo (#755 bright-adjacent) | observed |
| Known object `_221` | C | **known** | NED `SMACS J0723:[NDA2023] 01216` 0.10″, z 1.1616 ± 0.0003, `SLS`, 2023MNRAS.525.2087B; SIMBAD `[CAR2023] Beret` (LeG) 0.73″ | observed |
| Known object `_239` | C | **known: the BCG** | 0.53″ from the Mahler+2022 model centre; NED `[MJR2023] J110.82688-73.45463` 0.36″, z 0.391 ± 0.001 `SLS` (2025A&A...699A.105A); `WISEA J072318.38-732717.0` z 0.3914 `SLS`; cluster MCXC J0723.3-7327 z 0.404 | observed |
| Known object `_240` | C | **known: cluster member** | 1.8″ from the centre; NED `[MJR2023] J110.82639-73.45499` 0.18″, z 0.3915 `SLS` (2023MNRAS.525.1867N) | observed |
| Known object `_118` | C | **known: lensed multiple image** | Mahler+2022 ICLv0 image 3.1 at 0.19″ (system 3, z 1.9914; `f36a41c`, sha256 `35d568c4…`); NED `[NDA2023] 01149` 0.06″, z 1.992 ± 0.0003 `SLS` (2023MNRAS.525.1867N); SIMBAD `[CLL2024] … 3.3` (LeG\|LeI) 0.40″ | observed |
| Lens-model consistency `_118` | C | consistent | it is a constraint image of the model at its catalogued position; major axis 28–29° from tangential (catalog and moments), r 18.4″ | derived |
| Known object `_176` | C | probable (photometric) | NED `[ICR2022] 103` 0.32″, z 0.35 `PSE` (2022ApJ...940L..24I); DJA eazy z_phot 0.338 (95 %: 0.31–0.35), 0.10″ | observed / model_prediction |
| Why `_176` ranks | C | explains (hypothesis) | low sharpness = a resolved faint galaxy (10 px footprint, aper50 S/N ≈ 65); F1000W light visible but below its catalog threshold | hypothesis |
| Tangential arcs (other four) | C | no lensing anomaly | `_221` 11–26° from tangential at 36″ (a z 1.16 galaxy, a single image here); `_239`/`_240` are cluster galaxies at the centre | derived |

### Verdicts

- `_221`: `catalog effect: cross-band association (0.1″ merge radius vs 0.53–0.66″ offsets) on a deblended large
  galaxy at spectroscopic z = 1.1616`. Store: `explained`.
- `_239`: `known population: the SMACS 0723 BCG (z = 0.391, spectroscopic)`; "single band" is the same catalog
  effect. Store: `known_object`.
- `_240`: `known population: cluster member (z = 0.3915, spectroscopic) in the BCG halo`. Store: `known_object`.
- `_118`: `known population: lensed multiple image 3.1 of system 3 (z = 1.99, spectroscopic)`; consistent with
  the published model. Store: `known_object`.
- `_176`: `known population: faint galaxy at z_phot ≈ 0.34` (hypothesis; photometric only). Store: `triaged`.

None is unexplained; none is an anomaly or exotic-lensing candidate.

**Plan trigger met** (`PLAN_vet.md`): 4 / 5 owe their "single band" status to the same catalog effect, so the next
step is a gate change, not more vetting: associate MIRI bands by segment overlap or a size-scaled radius instead
of a fixed 0.1″ (TASKS Then 7).

**Regenerate** (cutouts via S3 byte ranges; photo-z table per SOURCES.md "Vetting (cycle 7)"):
```
jwst-anomaly run --config configs/smacs_jwst3_miri_confirm.yaml
python scripts/vet_evidence.py --run-dir $JWST_ANOMALY_OUTPUTS/runs/<run_id> \
  --sample smacs0723_miri --config configs/smacs_jwst3_miri_confirm.yaml --center 110.826750 -73.454628 \
  --lens-images https://raw.githubusercontent.com/guillaumemahler/SMACS0723-mahler2022/f36a41c365865df252a6da08cd79955f9cd68f11/ICLv0/arcs.dat \
  --photoz $JWST_ANOMALY_DATA/cache/dja/photoz-v7.4/smacs0723-grizli-v7.4-fix.eazypy.zout.fits \
  --uid jw02736-o002_t001_miri_f770w_221 --uid jw02736-o002_t001_miri_f770w_239 \
  --uid jw02736-o002_t001_miri_f770w_176 --uid jw02736-o002_t001_miri_f770w_118 --uid jw02736-o002_t001_miri_f770w_240
```
`vet_xband.json`: nearest-neighbour search in the four cached `jw02736-o002_t001_miri_<band>_cat.ecsv` files
(`sky_centroid`, astropy `separation`).

**Open questions:** does a size-aware band association change the MIRI top 20 beyond these four; is `_176`'s
F1000W light measurable with forced photometry.
