# Neutrino Frontier: public data audit (2026-10-10)

Scope: public data for testing neutrino hypotheses (multimessenger timing, flavor/propagation, nonlocal
correlations, anomaly detection). Already-pinned E1 inputs (ICECAT-1, GBM, GWTC, CHIME Cat 2) are in SOURCES.md
"E1 causal event network" and D-074; they are only re-checked here.

Legend: **[V]** verified by a request from the cloud container on 2026-10-10 (HTTP code given); **[B]** believed
from memory or from secondary text, not checked. Sizes are Dataverse/HTTP `Content-Length` values. Nothing > 20 MB
was downloaded. Fetched files sit in the session scratchpad only.

## Summary

| Source | Version / DOI | Size | Key fields | Flavor info | Host (code) | Recommended use |
|---|---|---|---|---|---|---|
| ICECAT-1 alert tracks | v4, doi:10.7910/DVN/SCRUCD | CSV 64 KB; FITS tars 5–38 MB/yr | UTC time µs; RA/Dec + asym. 90 % box; E (TeV); FAR; signalness | none (track; CNN topology scores) | dataverse.harvard.edu 200 [V] | keep (E1); alert-level timing |
| IceTracks-DR2 | v3.1, doi:10.7910/DVN/MMIIZA | 3.65 GB total; events 4.8–19.3 MB/season | MJD 1e-8 d; RA/Dec/AngErr (≥ 0.2°); log10 E_µ; az/zen; run/event | none (track) | dataverse 200 [V] | **primary** event-level sample; good-run lists |
| IceTracks-DR1 (10-yr PS) | v2.0, doi:10.7910/DVN/VKL316 (also HEASARC `icecubepsc`) | 409 MB; events 3.8–14.9 MB | as DR2 without run | none | dataverse 200, HEASARC TAP 200 [V] | superseded; TAP for server-side cone/time cuts |
| IceCat-2 | preliminary only, doi:10.7910/DVN/RX28YT | 3.6 MB (7 events) | maps + table | none | dataverse 200 [V] | not a catalogue yet; watch |
| HESE 12-yr (DirectFit) | v2.0, doi:10.7910/DVN/PZNO2T | data.tab 32 KB; FITS tars 1.9–2.2 GB | MJD; RA/Dec (rad); deposited E; asym. PDF params | topology Track/Shower | dataverse 200 [V] | cascades/starting events with times |
| HESE 7.5-yr | doi:10.21234/4EQJ-BB17 | zip 79 MB | deposited E, zenith, length | cascade/track/double-cascade | icecube.wisc.edu 200 [V] | flavor-morphology only; **no time, no RA** |
| Flavor composition 11.4 yr | v2.0, doi:10.7910/DVN/CBNMEB | 100 MB (posterior 97 MB) | binned histograms, Aeff, flavor contours | binned 3-morphology; no events | dataverse 200 [V] | flavor-ratio benchmark/priors |
| Galactic plane release | doi:10.21234/8v5d-rn16 | zip 27 MB | per-event LLH terms, Aeff | cascades (no flavor split) | icecube.wisc.edu 200 [V] | later; > 20 MB not opened |
| TXS 0506+056 2008–17 | doi:10.21234/B4QG92 | zip 31.5 KB | MJD; RA/Dec; σ; log10 E_reco (≤ 3° of TXS) | none | icecube.wisc.edu 200 [V] | **benchmark warm-up** |
| GW O3 × IceCube | doi:10.7910/DVN/34B5AP | 0.7 MB | per-GW p, ν E, RA/Dec, r90, Δt | none | dataverse 200 [V] | GW–ν null benchmark |
| GCN AMON gold/bronze table | live HTML | 399 KB | rev, date, time 0.01 s, RA/Dec, r90/r50, E, signalness, FAR | none | gcn.gsfc.nasa.gov 200 [V] | post-ICECAT-1 alerts (2023-10 →) |
| GCN AMON cascade table | live HTML | 42 KB | as above + sky-map URL | cascade (topology) | gcn.gsfc.nasa.gov 200 [V] | cascade alerts 2020 → |
| GCN Circulars | archive.json.tar.gz, daily | 31.4 MB (not fetched) | JSON per circular | n/a | gcn.nasa.gov 200 [V] | retractions/revisions; per-ID JSON |
| ANTARES 2007–17 PS tracks | v1.1, doi:10.5072/FK2/HZQTC5 (**test prefix**) | 316 KB | MJD 1e-4 d; RA/Dec 0.1°; Nhit; Beta (ang. err) | none (track) | opendata.km3net.de 200 [V] | southern-sky cross-check |
| KM3-230213A | v1.0, doi:10.5072/FK2/JW72C9 (test prefix) | ~53 MB (50 MB sim file) | hits, VOEvent, flux models | single muon track | opendata.km3net.de 200 [V] | single-event context only |
| KM3NeT/ORCA6 oscillation | v2.3, doi:10.5072/FK2/Y0UXVW | 14 MB ROOT | binned response, χ² landscapes | flavor via oscillation fit | opendata.km3net.de 200 [V] | propagation (GeV oscillation) only |
| SN1987A events | papers only (no official MR file found) | ~25 rows | time, E, angle per event | ν̄e (IBD-dominated) [B] | inspirehep 200 [V] | hand-transcribe, cite tables |
| SNEWS | alerts only | — | — | — | snews2.org 200 [V] | no event archive; nothing to test |
| Swift GRBs (BA compilation) | HEASARC `swiftgrbba` | 2,043 rows | trigger MJD, RA/Dec, error radius | n/a | HEASARC TAP 200 [V] | add to E1 GRB channel |
| Fermi GBM, GWTC, CHIME Cat 2 | as SOURCES.md E1 | — | — | n/a | 200 [V] | as D-074 |

## Per-source notes

**ICECAT-1** (Abbasi+ 2023, ApJS 269, 25, doi:10.3847/1538-4365/acfa95, arXiv:2304.01174) [V]. Dataverse v4.0,
released 2023-11-09, CC0. `IceCube_Gold_Bronze_Tracks.tab` 348 rows, last event IC231014A; md5 matches the repo pin.
Columns: NAME, RUNID, EVENTID, START (UTC to µs), EVENTMJD, I3TYPE (gfu-gold/bronze, ehe-gold, hese-gold/bronze),
RA/DEC with asymmetric 90 % error box, ENERGY (TeV, most probable, E^-2.19 assumed), FAR, SIGNAL,
CASCADE/SKIMMING/START/STOP/THRGOING_SCR (CNN topology scores), CR_VETO, OTHER_I3TYPES. Per-year FITS tars carry
HEALPix ΔLLH maps. Selection: 2011–2019 entries are archival re-selections and 2019 → are real-time. Signalness is
~30–50 %, so most alerts are atmospheric. There is no uptime list. Errors are a rectangle, not a PDF. Most events
are northern (Earth-filtered); the energy threshold depends on declination. Flavor: none. A track is νµ-CC dominated
but not flavor-tagged [B].

**IceTracks-DR2** [V]: this is the real name. It is the successor of the 10-year point-source track release
("IceTracks-DR1"), not an alert catalogue. Paper: Abbasi+ arXiv:2605.19040 (2026-05-18), "IceCube Second Track Data
Release IceTracks-DR2: Data from 2008-2022 for Neutrino Source Searches". Dataverse doi:10.7910/DVN/MMIIZA, first
published 2026-05-21, current v3.1 (2026-10-05), README last updated 2026-03-13, CC0. Coverage 2008-04-06 to
2022-05-23: 14 seasons, 1,643,355 events, 4,963 d livetime (paper Table 1).
- Files: `events/<season>_exp.tab` (whitespace table; run, event, subevent, MJD (1e-8 d ≈ 1 ms), log10(E/GeV) muon
  energy proxy, AngErr (deg, floor 0.2°), RA, Dec, Azimuth, Zenith); `uptime/<season>_exp.tab` (good-run
  MJD_start/stop); `irfs/` effective areas (~148 KB each) and smearing matrices (0.6–1.04 GB each, the bulk of the
  3.65 GB).
- Uses the same event numbering as ICECAT-1, so alerts can be matched to sample events.
- Selection effects for timing:
  - The seasons overlap by weeks (test processing), so dedupe on (run, event, subevent).
  - The azimuth response varies ~10 % below 1-day timescales. The README recommends scrambling time and recomputing
    RA from the local coordinates.
  - The southern sky (Dec < −5°) is dominated by atmospheric muons with a high energy cut.
  - One duplicated event is known in IC86-2015.
  - Do not combine DR2 with other releases.
- Benchmark tooling: SkyLLH (PyPI `skyllh` 26.1.0, 2026-09-09 [V]).

**IceTracks-DR1 / 10-year PS** [V]: Dataverse doi:10.7910/DVN/VKL316 v2.0 (a repost of the 2021 release; the
paper's own data DOI is 10.21234/CPKQ-K003), arXiv:2101.09836. 2008-04-06 to 2018-07-08, the same fields minus `run`.
HEASARC `icecubepsc` holds the same 1,134,450 events (MJD 54562.38–58307.97). Server-side ADQL cone/time cuts avoid
any download. Superseded by DR2 (DR2 has < 50 % event overlap in IC79/IC86-2011); use it only to reproduce older
papers.

**IceCat-2** [V]: no full release exists as of 2026-10-10. There is only "IceCat-2 Preliminary Data Release: Seven
Selected Events" (doi:10.7910/DVN/RX28YT, 2025-12-22; IC-170922A, NGC 7469 pair, three TDE-associated, IC-230724A),
plus the ICRC2025 proceedings arXiv:2507.06176. ICECAT-1 ends 2023-10-14. Later alerts come only from GCN.

**HESE.**
- 12-year DirectFit (doi:10.7910/DVN/PZNO2T, Yuan & Chirkin PoS ICRC2023 1030) [V]. `data.tab` has 164 events
  (109 Shower, 55 Track). Columns: id, mjd, ra, dec (radians), f0–f7 (asymmetric directional PDF parameters),
  reconstruction, energy (deposited), drlogl. Per-event FITS maps come only as 2.2 GB tars.
- 7.5-year (Abbasi+ PRD 104, 022002, arXiv:2011.03545; doi:10.21234/4EQJ-BB17) [V]. 102 events with
  `recoMorphology` 0/1/2 = cascade/track/double cascade (the double cascade is the ντ signature). There is **no event
  time or RA** (zenith only), so it cannot be used for timing. The zip is 79 MB.

Starting events see both hemispheres with a self-veto. Cascades have angular errors of ~10° or more [B].

**Flavor composition** (arXiv:2510.24957, 11.4 yr, f_e:f_µ:f_τ = 0.30:0.37:0.33) [V abstract], data
doi:10.7910/DVN/CBNMEB v2.0. Binned only: MESE track / cascade / double-cascade histograms, effective areas,
smearing JSON, flavor contours, posterior samples (97 MB). No per-event list. It is the only public, quantitative
flavor-ratio product.

**Benchmark-linked small releases** [V]:
- TXS 0506+056 (doi:10.21234/B4QG92, 31.5 KB). Six samples IC40 to IC86c (2008-04-05 to 2017-10-31). It holds events
  within 3° of TXS (MJD, RA, Dec, σ, log10 E_reco), per-sample Aeff at that Dec, and tabulated Fig. S4/S5.
- Alert catalogue to IceCube-170922A (doi:10.21234/B4KS6S, 1.9 KB).
- GW O3 × IceCube (doi:10.7910/DVN/34B5AP, arXiv:2601.07595). It holds the upper-limit table, per-GW best-neutrino
  table (p-value, E, RA/Dec, r90, Δt within ±500 s) and two ±500 s neutrino lists with their sky maps.

**GCN.**
- AMON gold/bronze table (https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html) [V]: 821 notice rows.
  - Dates: 2019-02-05 to 2026-09-30.
  - 633 distinct event times, 77 after the end of ICECAT-1.
  - Rows are per revision (`Rev`), not deduplicated, and RunNum_EventNum differs between duplicates.
  - 2019 alone has 447 distinct times against 19–34 per year later, so it holds test or non-alert entries. Filter
    against the circulars.
  - The energy column changed units: 2019 rows look like GeV (8.7e4 for IC190221A, which ICECAT gives as 56 TeV).
    This is inferred, not documented.
  - Newer rows state that the position error is statistical only.
  - No retraction flag in the table.
- Cascade table: 49 rows, 2020-07-01 to 2026-08-24, with duplicates (e.g. 134259/134262 at the same time), 90 %/50 %
  radii (stat+sys), FITS sky-map links.
- Retractions appear only as Circulars (e.g. GCN 43866, "IceCube cascade alert 142200_1517143 retraction").
- Circulars: per-ID JSON `https://gcn.nasa.gov/circulars/<id>.json` (fields subject, eventId, body, createdOn,
  bibcode) [V]. The daily archive `archive.json.tar.gz` is 31.4 MB, so it was not fetched here (it is fine under the
  owner-machine 200 MB rule).
- GCN also lists a HEASARC IceCube archive and Dataverse as IceCube data archives.

**KM3NeT / ANTARES** (Dataverse 6.3 at https://opendata.km3net.de, 3 datasets) [V].
- All DOIs use **10.5072/FK2**, the DataCite *test* prefix. doi.org redirects them to datacite.org/testprefix, so they
  are not persistent. Cite them by the portal URL plus version plus access date.
- `open-data.km3net.org` gives a proxy 502 (CONNECT failed, no x-deny-reason; probably not a live host).
  antares.in2p3.fr reset once and then returned 200.
- ANTARES 2007–2017 PS (CC BY 4.0; Illuminati+ PoS ICRC2019 920; Albert+ ApJL 863, L30 2018): `events_tabulated.tab`
  has 8,753 rows (the description says 8,754), columns ID, Decl, RA, Nhit, Beta, MJD. MJD is at 1e-4 d (≈ 9 s).
  Coordinates are at 0.1°. Beta is the angular-error estimate (median 0.4°). There is no uptime list (3,125 d
  livetime only), so a timing null needs time permutation rather than an exposure model.
- KM3-230213A (Nature 638, 376, 2025, doi:10.1038/s41586-024-08543-1 [V]): one ~100+ PeV muon event (energy
  [B]).
- ORCA6 oscillation: binned response for propagation fits, not events.

**SN1987A** [V for references, B for content]. No official machine-readable event list was found (Zenodo and arXiv
searches returned only talks and posters).
- Primary tables:
  - Kamiokande-II: Hirata+ 1987 PRL 58, 1490 (10.1103/PhysRevLett.58.1490) and Hirata+ 1988 PRD 38, 448
    (10.1103/PhysRevD.38.448);
  - IMB: Bionta+ 1987 PRL 58, 1494 (10.1103/PhysRevLett.58.1494) and Bratton+ 1988 PRD 37, 3361 (angles);
  - Baksan: Alexeyev+ 1988 PLB 205, 209 (10.1016/0370-2693(88)91651-6).
- Compiled table: Loredo & Lamb 2002 PRD 65, 063002 (astro-ph/0107260).
- Expect ~11–12 (KII) + 8 (IMB) + 5 (Baksan) events [B]. Absolute clocks: IMB is good; KII has a ±1 min absolute
  offset; Baksan has −54 s/+2 s [B].
- Transcribe once, with a checksum and a per-row source citation, and label it "observed (transcribed)".

**SNEWS** [V]: it issues alerts only (GCN notice type SNEWS; ~52 test alerts per year, 0.03 real per year). There is
no public event archive and no Galactic SN since 1987, so there is no data product to test.

**Swift / GBM / GW / FRB** [V]:
- Swift: HEASARC `swiftgrbba` (Burst Advocate compilation) has 2,043 rows, MJD 53356.3 to 61320.24 (GRB 261007A).
  It includes non-BAT triggers (e.g. SVOM), so filter `trigger_obs`. `swiftgrb` is the old BAT catalogue: 872 rows,
  2004-12 to 2012-12.
- `fermigbrst` now has 4,391 rows (the repo pin has 4,390; the table grows).
- gwosc.org, CANFAR (CHIME Cat 2) and Zenodo all returned 200.

## Benchmarks a small project can reproduce

| # | Benchmark | Reference | Number to reproduce | Data |
|---|---|---|---|---|
| 1 | **TXS 0506+056 2014–15 flare, fixed box window** | IceTracks-DR2 paper §V and Table 6, arXiv:2605.19040 | box T0 = MJD 57020, ΔT = 185 d: n̂s = 12.7, γ̂ = 2.3 (SkyLLH). Internal tools give 12.56 / 2.26, pre-trial p = 4.18e-3 | DR2 IC86 season events + uptime + IC86 Aeff/smearing |
| 2 | NGC 1068 time-integrated | same, Table 8 | n̂s = 80.1, γ̂ = 3.2, pre-trial p = 1.3e-7 (5.1σ) with DR2+SkyLLH; published 71.1 / 3.1 / 2.1e-6 | all DR2 seasons |
| 3 | TXS flare, original | Aartsen+ 2018 Science 361, 147 (doi:10.1126/science.aat2890) | 3.5σ time-dependent excess, Sep 2014 – Mar 2015 [V abstract]; 13 ± 5 signal events, Gaussian σ_T ≈ 110 d [B] | 31.5 KB TXS release |
| 4 | GRB stacking null | Aartsen+ 2017 ApJ 843, 112 (arXiv:1702.06868) | 1,172 GRBs, no significant correlation [V abstract] | needs IceCube's own event selection; not reproducible exactly from public tracks |
| 5 | GW–ν null | Abbasi+ 2023 ApJ 944, 80 (arXiv:2208.09532); O3 deep search arXiv:2601.07595 | per-GW p-values and Δt in `table1.tab`, e.g. GPS 1262142545.615: p = 3.8e-4, Δt = −222 s [V] | DVN/34B5AP |
| 6 | FRB–ν null | Aartsen+ 2018 ApJ 857, 117; 2020 ApJ 890, 111; Masaoka+ arXiv:2603.24983 (best post-trial p = 0.076, FRB20190630C–IC190629A) | null / p = 0.076 | ICECAT-1 × CHIME Cat 2 (already in repo) |

**Best first benchmark: #1.**
- It is a public-tool reproduction that the collaboration itself published, with exact window and fit values.
- The events file for IC86_IV (2014 season) is 17 MB.
- Cost: the IC86 smearing matrix is 817 MB. Stream it in the cloud, or state the reason on the owner machine
  (> 200 MB rule).
- Warm-up: #3's 31.5 KB release, which gives the in-window event list near TXS with no IRF download.
- #6 is the cheapest null to reproduce with the existing E1 inputs.

## Ready-to-paste SOURCES.md lines (accessed 2026-10-10)

```
## Neutrino Frontier data audit (accessed 2026-10-10; docs/neutrino_frontier/data_audit.md)
- IceCube IceTracks-DR2, Harvard Dataverse doi:10.7910/DVN/MMIIZA (v3.1, 2026-10-05; CC0), 2008-04-06 to
  2022-05-23, 1,643,355 track events; paper arXiv:2605.19040.
- IceCube IceTracks-DR1 (10-yr PS), doi:10.7910/DVN/VKL316 (v2.0) / data DOI 10.21234/CPKQ-K003, arXiv:2101.09836;
  same events as HEASARC TAP table `icecubepsc` (1,134,450 rows).
- IceCube HESE 12-yr DirectFit, doi:10.7910/DVN/PZNO2T (v2.0; `data.tab`, 164 events), PoS(ICRC2023)1030.
- IceCube HESE 7.5-yr, doi:10.21234/4EQJ-BB17 (zip 78,998,573 B), PRD 104, 022002 (arXiv:2011.03545).
- IceCube flavor composition 11.4 yr, doi:10.7910/DVN/CBNMEB (v2.0), arXiv:2510.24957.
- IceCube TXS 0506+056 2008-2017 events, doi:10.21234/B4QG92 (zip 31,520 B), Science 361, 147
  (doi:10.1126/science.aat2890).
- IceCube GW O3 joint search replication data, doi:10.7910/DVN/34B5AP (v1.0), arXiv:2601.07595.
- IceCat-2 preliminary (7 events), doi:10.7910/DVN/RX28YT (v1.0); proceedings arXiv:2507.06176.
- GCN AMON IceCube gold/bronze and cascade tables, https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html,
  https://gcn.gsfc.nasa.gov/amon_icecube_cascade_events.html (live pages; snapshot sha256 to be pinned on use).
- GCN Circulars, https://gcn.nasa.gov/circulars/<id>.json and archive.json.tar.gz (31.4 MB, daily).
- ANTARES 2007-2017 point-source tracks, https://opendata.km3net.de dataset "ANTARES 2007-2017 Point Source
  Analysis" v1.1 (DOI 10.5072/FK2/HZQTC5 is a DataCite test prefix, not persistent; CC BY 4.0).
- KM3-230213A event data, https://opendata.km3net.de v1.0 (test-prefix DOI 10.5072/FK2/JW72C9); Nature 638, 376
  (doi:10.1038/s41586-024-08543-1).
- Swift GRB Burst Advocate compilation, HEASARC TAP table `swiftgrbba` (2,043 rows on 2026-10-10).
- SkyLLH, PyPI `skyllh` 26.1.0 (2026-09-09), https://github.com/icecube/skyllh.
- SN1987A: Hirata+ 1988 PRD 38, 448; Bionta+ 1987 PRL 58, 1494; Bratton+ 1988 PRD 37, 3361; Alexeyev+ 1988
  PLB 205, 209; compilation Loredo & Lamb 2002 PRD 65, 063002 (astro-ph/0107260).
```
