# Data inventory: distance self-consistency and causal event networks

Research note, 2026-10-09. Everything below was opened and checked on 2026-10-09 unless it is marked
**unverified**. Counts marked "(counted)" come from files downloaded in this session. All other counts
are quoted from the source. Provenance labels: *observed* (catalogue measurement), *derived*
(posterior or model output), *assumption*.

## 1. Distance self-consistency

The aim is per-object consistency residuals: two independent distance or path-length measures of the
same object, each with its own error, compared with each other. This is not a cosmology fit.
Nearly every published analysis below collapses its objects into one shared parameter (H0, Ωm, η).
None publishes a leave-one-out residual per object with the full covariance. That is what an
outlier search would add.

### 1.1 Inventory

| Family | Dataset | N | Two measures per object | Access | Per-object values and errors | Mainstream analysis already done |
|---|---|---|---|---|---|---|
| Time-delay lenses (galaxy scale) | H0LiCOW public release | 6 lensed quasars | D_Δt (time delay + lens model) and D_d (stellar kinematics): joint chains for RXJ1131, PG1115, J1206; a parametric (D_d, D_Δt) fit for B1608; D_Δt only for HE0435 and WFI2033 | `github.com/shsuyu/H0LiCOW-public`, folder `h0licow_distance_chains/` (~120 MB; PG1115 alone is 84 MB). Last commit 2025-05-14 | Yes. *Derived* posterior samples | Joint H0 from all 6 lenses (Wong+2020, arXiv:1907.04869) |
| Time-delay lenses | TDCOSMO 2025 | 8 quasars, plus 11 SLACS and 4 SL2S kinematics-only lenses | D_Δt per lens; D_d given explicitly only for SDSS1206 (`final_D_d.npy`, `final_D_dt.npy`). For the others, kinematics enter as pre-processed `*_const_processed.pkl` likelihoods inside a hierarchical mass-sheet (λ_MST) model | `github.com/TDCOSMO/TDCOSMO2025_public` (last commit 2026-01-21). Cosmology chains are in `chains_export/*.h5` | Partly. Per-lens D_Δt chains are files; the kinematic distance is only implicit | H0 = 71.6 (+3.9/−3.3) with Pantheon+ Ωm. TDCOSMO 2025, A&A 704, A63, arXiv:2506.03023 |
| Lensed SN, cluster | SN Refsdal (MACS J1149) | 1 SN, 5 images | Time delay SX–S1 = 376.0 (+5.6/−5.5) d and magnification ratios; compared with lens-model predictions | Tables in Kelly+2023 ApJ (arXiv:2305.06377) | Yes, in the paper tables | H0 = 64.8 (+4.4/−4.3) from 8 models, or 66.6 (+4.1/−3.3) from the 2 best models (Kelly+2023 Science, arXiv:2305.06367) |
| Lensed SN Ia, cluster | SN H0pe (PLCK G165) | 1 SN, 3 images | Time delays (≈ −116.6 d and −48.6 d relative to image b) **and** absolute magnifications from the standard candle (μa = 4.3, μb = 7.6, μc = 6.4, each with an error) | Pierel+2024 ApJ 967, 50 (arXiv:2403.18954) | Yes | H0 = 71.8 (+9.2/−8.1) from 7 lens models weighted by agreement with the SN μ (Pascale+2025, arXiv:2403.18902) |
| Lensed SN Ia, cluster | SN Encore and SN Requiem (MACS J0138) | 2 SNe, same host | Encore: Δt(1b,1a) = −39.8 (+3.9/−3.3) d. Requiem: the 4th image is predicted for 2037 ± 2 (Rodney+2021, arXiv:2106.08935) | Pierel+2026 ApJ 998 (arXiv:2509.12301); 7 blind lens models (arXiv:2509.12319) | Yes | H0 = 66.9 (+11.2/−8.1) |
| Lensed SN Ia, galaxy | SN Zwicky (2022qmx) | 1 SN, 4 images | Standard-candle μ per image vs lens-model μ: two images are ≈1.7 and ≈0.9 mag brighter than predicted. Delays ≲ 2 d | LensWatch I (arXiv:2211.03772) and II (ApJ 980, 172; arXiv:2409.17239) | Yes | Interpreted as microlensing or substructure |
| Lensed SN Ia, galaxy | iPTF16geu | 1 SN, 4 images | Total μ = 67.8 (+2.6/−2.9); delays consistent with ≈0 (≈1 d errors); per-image flux anomaly vs smooth models | Dhawan+2020 MNRAS (arXiv:1907.06756); Goobar+2017 Science (arXiv:1611.00014) | Yes | Only a lower limit on H0 |
| GW bright siren | GW170817 + NGC 4993 | 1 | d_L(GW) = 40 Mpc (GWOSC default) vs SBF distance 40.7 ± 1.4 ± 1.9 Mpc (Cantiello+2018, arXiv:1801.06080). GW–GRB arrival difference +1.74 ± 0.05 s (arXiv:1710.05834, ApJL 848, L13) | GWOSC event API; papers | Yes | H0 (arXiv:1710.05835); the speed of gravity bounded to −3e−15 … +7e−16 |
| GW catalogue | GWTC cumulative (GWOSC) | 391 events (counted), 282 with a d_L summary. Release: GWTC-5.0 (O4b), 2026-05-26 | d_L posterior per event; a redshift only from a counterpart (GW170817) or statistically from galaxy catalogues or the mass spectrum | CSV `gwosc.org/api/v2/catalogs/GWTC/events?include-default-parameters=true&format=csv` (110 KB). PE samples and sky maps on Zenodo (GWTC-5.0: 10.5281/zenodo.20348005 and .20348006; GWTC-4.0 O4a: zenodo 17014085) | d_L median with 90% bounds in the CSV; full samples in HDF5 (large) | Dark- or spectral-siren H0 (e.g. GWTC-4.0 cosmology, arXiv:2509.04348). Per-event residuals have little power without a redshift |
| Localized FRBs | FRBs/FRB repo `frb/data/Galaxies/public_hosts.csv` + `frb/data/FRBs/*.json` | 94 hosts with z (counted; z_max = 1.016; projects DSA 30, CRAFT 28, CHIME 21) | Dispersion measure (path-integrated n_e) vs host spectroscopic z | `github.com/FRBs/FRB` (BSD-3, last commit 2026-05-06). Per-FRB JSON has DM, DM_err, DM_ISM, error ellipse | Yes for DM; the cosmic-DM scatter and host-DM are model terms (*assumption*) | Macquart relation (Macquart+2020 Nature, arXiv:2005.13161). The CHIME Outriggers host sample (ApJS 280, 6) adds low z |
| Clusters, SZ + X-ray D_A | Bonamente+2006 | 38 clusters | D_A from SZ+X-ray vs z | Paper tables (ApJ 647, 25; astro-ph/0512349). Not in VizieR (checked) | Yes, in the paper | H0. Distance-duality tests combine these with SNe (e.g. Uzan+2004; Holanda+2011) |
| Clusters, SZ + X-ray | Kozmanyan+2019 (Planck + XMM) | 61 clusters, z < 0.5 | D_A per cluster from resolved y and X-ray profiles | A&A 621, A34 (arXiv:1809.09560). Per-cluster table **unverified** | Unknown | H0 = 67 ± 3 |
| Clusters, WL vs hydrostatic mass | CoMaLit LC² (VizieR J/MNRAS/450/3665) + CoMaLit-I comparison | LC²-single 505 rows; Sigma catalogue 564 unique (counted via VizieR) | WL mass vs X-ray HE mass for the overlapping clusters | VizieR TAP | Yes, with literature heterogeneity | ~15% (WL) and ~25% (HE) intrinsic scatter; group-to-group bias up to ~40% (Sereno & Ettori 2015, arXiv:1407.7868) |
| Clusters, SZ mass | Planck PSZ2 (VizieR J/A+A/594/A27) | 1653 (counted) | M_SZ (hydrostatic-calibrated) to cross-match with WL masses | VizieR | Yes | "Hydrostatic mass bias" (1−b) |

### 1.2 What a per-object search adds, and what limits it

- **Lenses:** D_Δt/D_d = (1+z_d)·D_s/D_ds does not depend on H0 and depends only weakly on Ωm. Each lens
  therefore gives a nearly cosmology-free internal consistency number. TDCOSMO lets λ_MST absorb
  lens-to-lens differences as a population scatter, so a lens that needs an unusual λ is not reported as
  an outlier. The leave-one-out predictive residual per lens is the new product.
- **Lensed SNe:** the number of objects is tiny (≈5), and microlensing (SN Zwicky, iPTF16geu) is the
  ordinary explanation for magnification residuals. Each case is a vetted single-object test, not a
  population.
- **FRBs:** the per-object residual DM_obs − DM_ISM − DM_halo − ⟨DM_cosmic⟩(z) has a heavy-tailed,
  z-dependent intrinsic scatter (the F parameter) plus a log-normal host term. "Honest errors" must
  therefore be the predictive distribution, not a Gaussian.
- **GW:** without a redshift, the residual is a check of d_L against a galaxy-catalogue prior. Only
  GW170817 has a clean two-measure test.
- **Clusters:** triaxiality and clumping bias D_A at the 10–20% level per object, the ordinary
  systematic. Outliers are expected to be mostly merging clusters.

## 2. Causal event network

### 2.1 Catalogues

| Events | Dataset | N | Time precision | Localization | Access |
|---|---|---|---|---|---|
| GRB (GBM) | HEASARC `fermigbrst` | 4,390 bursts (TAP count), MJD 54661.1–61320.6 (2008-07 to 2026-10) | Trigger time, sub-second | `Error_Radius` = 1σ statistical, in degrees; 0 means localized by another instrument. Typically degrees. Whether a systematic term is included is not stated | HEASARC TAP `heasarc.gsfc.nasa.gov/xamin/vo/tap`, table updated 2026-10-07 |
| GRB (Swift) | Swift GRB table; HEASARC `swiftgrb` (872 rows) | Hundreds to ~1,000+; the page gives no total | Trigger time | BAT and XRT RA/Dec with 90% error radii (arcmin and arcsec) | `swift.gsfc.nasa.gov/archive/grb_table/` (tab-delimited; updated 2025-07-07) |
| FRB | CHIME/FRB Catalog 2 | 4,539 bursts from 3,641 sources (83 repeaters), 2018-07-25 to 2023-09-15 | Cat 1 schema has `mjd_400`/`mjd_inf` with `_err` (ms-level). Cat 2 columns are **unverified** | O(10′) per the paper; `ra_err`, `dec_err` | CANFAR DOI 10.11570/25.0066 → `data/table/chimefrbcat2.csv` (3.87 MB) plus `.fits`, `.json`. `data/exposure/chimefrbcat2_exposure.h5` (206 MB). The download was reset by the proxy in this sandbox |
| FRB, localized | FRBs/FRB repo | ~190 FRB JSONs, 94 hosts with z | Varies | Arcsecond (error ellipse in the JSON) | git |
| GW | GWTC (GWOSC) | 391 confident (GWTC-1 to GWTC-5.0) | GPS time, ms | Sky maps typically 10s to 1000s deg², only in the Zenodo PE and sky-map releases | GWOSC API (above) |
| ν (alert tracks) | ICECAT-1 v4 | 348 tracks (counted), 2011-05-14 to 2023-10-14 | Microsecond timestamp (`START`, `EVENTMJD`) | Asymmetric 90% CL RA/Dec rectangle (~1–7°); HEALPix ΔLLH maps in the per-year tars (20–38 MB each) | Harvard Dataverse doi:10.7910/DVN/SCRUCD; summary `IceCube_Gold_Bronze_Tracks.tab` (67 KB) |
| ν (all tracks) | IceCube PS 2008–2018 | 1,134,450 events (TAP count) | Event MJD | Per-event angular error (~1° and worse) | HEASARC `icecubepsc`; Dataverse doi:10.7910/DVN/VKL316. Superseded by IceTracks-DR2, 2008–2022 (doi:10.7910/DVN/MMIIZA, released 2026-10-05; smearing files ~1 GB each) |
| SN / TDE / transients | TNS public objects | ~160 K (Jan 2025, per the TNS page) | Discovery date; the true onset is uncertain by days (last non-detection) | Arcsec | Daily `tns_public_objects.csv.zip`. Requires a TNS login or bot `tns_marker` header with API key. No explicit licence |
| Fast X-ray transients | Einstein Probe | ~128 transients in year 1 (Wu+2025, quoted via secondary sources; **unverified**) | Seconds | WXT ~arcmin; FXT ~arcsec | No public WXT transient catalogue was found. FXT level 1–3 data are public after a 1-year proprietary period (first batch: 1,615 obs IDs, via ESA EP archive/NADC/NSSDC). Practical route: GCN Circulars bulk archive `gcn.nasa.gov/circulars/archive.json.tar.gz` (31.4 MB; ~45,870 circulars, which also carry EP, SVOM, GRB and GW follow-ups) |

### 2.2 Ordinary sources of spurious event–event correlation, and standard controls

| Effect | How it creates fake links | Usual control |
|---|---|---|
| Instrument duty cycle (SAA passages for GBM; GW detector livetime; CHIME up/downtime) | Event times cluster where the instruments are on; two catalogues with correlated uptime (e.g. both down for upgrades) give apparent temporal clustering | Scramble times within each instrument's own good-time intervals; compare only overlapping livetime |
| Sky exposure / declination dependence (CHIME transit, IceCube zenith acceptance, GBM Earth occultation ~30% of sky) | Spatial overdensities at declinations favored by both instruments | Declination-preserving scrambling (keep Dec, randomize RA or time, i.e. sidereal scrambling); exposure maps (CHIME `exposure.h5`, IceCube effective area per declination) |
| Sun/Moon and Galactic-plane constraints (optical TNS, X-ray pointing limits) | Optical transients avoid the Sun, giving an annual RA modulation shared by any optical catalogue | Scramble within the same day of year, or model the visibility window |
| Alert follow-up chains (GW or ν alert → EP/Swift/ZTF targeted search → TNS report) | A downstream detection caused by the alert itself: a causal link created by the observing strategy, not the sky | Exclude events found in targeted follow-up (check `reporting_group`/circular references), or flag "found because of" edges from the GCN circular graph |
| Shared or duplicate triggers (GBM bursts localized by Swift: `Error_Radius = 0`; one GRB in several catalogues; IceCat events in several streams, `OTHER_I3TYPES`; repeater bursts from one FRB source) | Self-matches and the same physical event counted as many nodes | Deduplicate by source ID before building edges; treat repeaters as one node with many times |
| Background events (IceCat `CR_VETO`; signalness ~0.3–0.5; GBM non-GRB triggers; GW P_astro) | Raise the noise floor and give false "hubs" | Weight edges by signalness or P_astro; drop CR-veto events |
| Look-elsewhere / trials over time windows and radii | Inflated significance | Pre-registered windows; post-trial p-values from many scrambles. Curtin+2023 (ApJ 954, 154) Monte-Carlo-simulated GRBs and FRBs with the catalogue's Dec and time distributions; Masaoka+2026 (arXiv:2603.24983; CHIME Cat 2 × ICECAT-1) accounted for Dec-dependent exposure and trials |
| Localization error mis-modelling (GBM systematics, non-Gaussian IceCube contours) | Too many or too few spatial matches | Use catalogue likelihood maps where available (IceCat HEALPix, GW sky maps), or inflate radii and test the sensitivity |

## 3. Smallest first tests

Both tests use less than 200 MB, run on a laptop in under an hour, need no credentials, and have a
built-in positive and null control.

### Test A (distance consistency): per-lens D_Δt vs D_d leave-one-out residuals (~125 MB)

1. Download `h0licow_distance_chains/` and `MontePython_cosmo_sampling/data/timedelay_6lenses/B1608_Dd_Ddt_params.dat`
   from `shsuyu/H0LiCOW-public` (~120 MB), plus `TDCOSMO_sample/TDCOSMO_data/SDSS1206+4332/final_D_*.npy`
   from TDCOSMO2025_public (a few MB).
2. For the 4 lenses with joint (D_d, D_Δt), compute the H0-free ratio
   R_i = D_Δt/[(1+z_d)·D_d] = D_s/D_ds. The redshifts come from the papers or the likelihood code.
   Compare R_i with the flat-ΛCDM prediction over a broad Ωm prior (0.1–0.5). The residual is the
   predictive p-value of each lens's (D_d, D_Δt) samples.
3. For all 6 lenses: fit H0 on 5 lenses (D_Δt only), predict the 6th, and report z = Δ/σ_total from KDE
   densities. Expected outcome: all |z| < 2. H0LiCOW reported mutual consistency, so this test calibrates
   the machinery, not a discovery.
   - Null control: shuffle the lens redshifts across lenses; residuals should get worse.
   - Provenance: *derived* posteriors only.

The cheaper fallback (<5 MB) is the FRB DM–z per-object predictive residual on the 94 `public_hosts.csv`
hosts, using the `frb` package's ⟨DM_cosmic⟩(z) and a published F-parameter scatter model.

### Test B (event network): GRB × ν × GW × FRB time–space pair excess with scrambled nulls (~5 MB)

1. Fetch the GBM `fermigbrst` table (TAP; trigger time, RA, Dec, Error_Radius; <2 MB), ICECAT-1
   `IceCube_Gold_Bronze_Tracks.tab` (67 KB), the GWTC CSV (110 KB; times only), and CHIME Cat 2
   `chimefrbcat2.csv` (3.9 MB; if CANFAR stays unreachable, use the Cat 1 CSV).
2. Build edges for pairs with |Δt| ≤ {10 s, 1 d, 7 d} and angular separation ≤ the quadrature sum of
   the 90% radii (convert GBM 1σ to 90%, ×2.15 for a 2-D Gaussian). Restrict to overlapping livetime.
3. Null: 10⁴ scrambles that keep each event's declination and permute times within the same catalogue
   and the same calendar year. This preserves exposure and duty-cycle structure without the 206 MB
   CHIME exposure file. Report the observed edge count vs the null per window and catalogue pair.
   - Positive control: GW170817 × GBM 170817A (Δt = 1.74 s) must appear in the 10 s window.
   - Negative control: GBM bursts with `Error_Radius = 0` vs Swift-localized bursts expose shared-trigger
     duplicates.
   - Compare the CHIME × ICECAT result with Masaoka+2026 (no significant pair; best post-trial p = 0.076).

## Sources (accessed 2026-10-09)

- GWOSC GWTC event list and API, https://gwosc.org/eventapi/html/GWTC/ ; GWTC-5.0 docs https://gwosc.org/GWTC-5.0/ (release 2026-05-26; DOI 10.7935/bk00-6a89); GWTC-4.0 PE https://zenodo.org/records/17014085
- LVK+Fermi+INTEGRAL 2017, ApJL 848, L13, doi:10.3847/2041-8213/aa920c (arXiv:1710.05834); LVK 2017 Nature standard siren, doi:10.1038/nature24471 (arXiv:1710.05835); Cantiello+2018 ApJL, doi:10.3847/2041-8213/aaad64 (arXiv:1801.06080); GWTC-4.0 cosmology arXiv:2509.04348
- H0LiCOW public, https://github.com/shsuyu/H0LiCOW-public (HEAD 2025-05-14); Wong+2020 MNRAS, doi:10.1093/mnras/stz3094 (arXiv:1907.04869)
- TDCOSMO 2025, A&A 704, A63, doi:10.1051/0004-6361/202555801 (arXiv:2506.03023v4); https://github.com/TDCOSMO/TDCOSMO2025_public (HEAD 2026-01-21)
- Kelly+2023 ApJ, doi:10.3847/1538-4357/ac4ccb (arXiv:2305.06377); Kelly+2023 Science 380, abh1322, doi:10.1126/science.abh1322 (arXiv:2305.06367)
- Pierel+2024 ApJ 967, 50, doi:10.3847/1538-4357/ad3c43 (arXiv:2403.18954); Pascale+ (arXiv:2403.18902)
- Pierel+ SN Encore (arXiv:2509.12301; ApJ 998, 2026); Suyu+ lens models (arXiv:2509.12319); Rodney+2021 (arXiv:2106.08935)
- Pierel+2023 LensWatch I, doi:10.3847/1538-4357/acc7a6 (arXiv:2211.03772); Larison+2025 ApJ 980, 172, doi:10.3847/1538-4357/ada776 (arXiv:2409.17239)
- Goobar+2017 Science 356, 291, doi:10.1126/science.aal2729 (arXiv:1611.00014); Dhawan+2020 MNRAS, doi:10.1093/mnras/stz2965 (arXiv:1907.06756)
- FRBs/FRB repo, https://github.com/FRBs/FRB (HEAD 996fcda, 2026-05-06, BSD-3); Macquart+2020 Nature, doi:10.1038/s41586-020-2300-2 (arXiv:2005.13161)
- CHIME/FRB Catalog 2, arXiv:2601.09399 (2026-01-14); data https://www.canfar.net/storage/list/AstroDataCitationDOI/CISTI.CANFAR/25.0066/data ; Cat 1 schema https://chime-frb-open-data.github.io/catalog/
- Bonamente+2006 ApJ 647, 25, doi:10.1086/505291 (astro-ph/0512349); Kozmanyan+2019 A&A 621, A34 (arXiv:1809.09560); Uzan+2004 PRD 70, 083533 (astro-ph/0405620); Holanda+2011 A&A 528, L14 (arXiv:1003.5906)
- Sereno & Ettori 2015 MNRAS, doi:10.1093/mnras/stv810 (arXiv:1407.7868); VizieR J/MNRAS/450/3665, J/MNRAS/450/3675, J/A+A/594/A27
- HEASARC fermigbrst https://heasarc.gsfc.nasa.gov/W3Browse/fermi/fermigbrst.html ; Swift GRB table https://swift.gsfc.nasa.gov/archive/grb_table/ ; HEASARC TAP https://heasarc.gsfc.nasa.gov/xamin/vo/tap
- ICECAT-1, Abbasi+2023 ApJS 269, 25, doi:10.3847/1538-4365/acfa95 (arXiv:2304.01174); data doi:10.7910/DVN/SCRUCD (v4); IceCube PS 2008–2018 doi:10.7910/DVN/VKL316 (arXiv:2101.09836); IceTracks-DR2 doi:10.7910/DVN/MMIIZA
- TNS getting started, https://www.wis-tns.org/content/tns-getting-started
- Einstein Probe data policy, https://www.cosmos.esa.int/web/einstein-probe/data ; GCN Circulars archive https://gcn.nasa.gov/circulars
- Curtin+2023 ApJ 954, 154 (arXiv:2208.00803); Masaoka+2026 (arXiv:2603.24983)
