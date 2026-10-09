# S1 burst twins (Fermi GBM)

This is a wide-separation, any-delay light-curve twin search over the Fermi GBM burst catalogue. It tests
hypothesis S1 from docs/hypotheses/round-1/summary.md. Decision: DECISIONS.md "S1 burst twins" (D-071). Sources:
SOURCES.md "S1 burst twins". A match is an anomaly, never evidence of new physics.

Reproduce (about 3 min of streaming, then about 7 min on 4 cores):

    python scripts/s1_ingest.py            # stream bcat files -> lc_<YYYY>.ecsv.gz (resumable; raw FITS never stored)
    python scripts/s1_twins.py --null 100 --inject 5000 --plots <scratch dir>

| File | Provenance | Content |
|---|---|---|
| `catalogue.ecsv.gz` | observed | HEASARC `fermigbrst` subset (TAP response sha256 in meta) |
| `lc_<YYYY>.ecsv.gz` | derived | per-burst two-band light curves (bcat HDU 2), bcat file name, size and sha256 |
| `bursts.ecsv.gz` | derived | per-burst pulse count, peak S/N, eligibility |
| `pairs_top.ecsv` | model_prediction | top 300 real pairs by ρ (s = 1) plus flagged pairs, with every vetting column |
| `null.json` | model_prediction | pulse-shuffled surrogate null: per-catalogue maxima, ρ histograms (real and null) |
| `injections.ecsv.gz` | simulated | synthetic twins through the whole chain |
| `summary.json` | derived | headline numbers, efficiencies, limits, positive controls |

Light-curve encoding: band b's flux and error are integer tenths of the unit `u<b>` (ph cm⁻² s⁻¹), and `nan`
marks an invalid bin. Bin i covers [t_lo + i·dt, t_lo + (i+1)·dt) s from the trigger. Band 1 is 50–300 keV.
Band 2 is 10–50 plus 300–1000 keV. `jwst_anomaly.burst_twins.decode` reads the strings back.
