# Exotic-lens signatures: what the literature predicts and what our data can test

Survey of 2026-10-08 (D-023). Every reference was opened on arXiv; titles are exact. This page only says what
would distinguish an exotic lens from ordinary lensing. Any detection would be a hypothesis to test, never a
conclusion (charter).

## Signatures

| Lens | Predicted signature | Testable here? |
|---|---|---|
| Ellis wormhole (the "1/r^n" lens family, n = 2) | **Demagnification**: sources beyond about 2/(n+1) Einstein radii appear fainter than unlensed. Images are still stretched tangentially. Light curves dip by ~4% just outside each Einstein-ring crossing. The image centroid shift is a few µas and traces an oval path. | Static demagnification against a cluster model: in principle, from one epoch. Light curves and µas shifts: no. |
| Negative mass / negative convergence | **Radially stretched images**; an umbra (eclipsed region); light curves with a central dip flanked by spikes; an early, not delayed, arrival time | Radial stretching plus demagnification around one dark centre that no positive-mass perturber can fit: yes, in principle |
| Alcubierre-type warp bubble | No published imaging or light-curve prediction for a distant observer, only views from or near the bubble and gravitational waves from bubble collapse | No; a model would have to come first |

Key references:
- Abe, arXiv:1009.6084, "Gravitational Microlensing by the Ellis Wormhole".
- Kitamura, Nakajima & Asada, arXiv:1211.0379, "Demagnifying gravitational lenses toward hunting a clue of exotic
  matter and energy".
- Izumi et al., arXiv:1305.5037, "Gravitational lensing shear by an exotic lens object with negative convergence
  or negative mass".
- Cramer et al., arXiv:astro-ph/9409051, "Natural Wormholes as Gravitational Lenses".
- Safonova, Torres & Romero, arXiv:gr-qc/0105070, "Microlensing by natural wormholes: theory and simulations".
- Asada, arXiv:1711.01730, "Gravitational lensing by exotic objects" (review).
- Warp drives:
  - Clark, Hiscock & Larson, arXiv:gr-qc/9907019, "Null geodesics in the Alcubierre warp drive spacetime: the view
    from the bridge";
  - Clough, Dietrich & Khan, arXiv:2406.02466 (gravitational waves from warp-drive collapse).

## Where the signatures overlap with ordinary effects

The literature itself names these degeneracies, and each must be excluded first.
- **Radial images:** they also form near a cluster's radial critical curve, and voids give negative convergence
  (Izumi et al.).
- **Demagnified images:** also from saddle-point images and from microlensing by intracluster stars, e.g. Kelly
  et al., arXiv:1706.10279, "Extreme magnification of a star at redshift 1.5 by a galaxy-cluster lens".
- **Flux-ratio anomalies:** also from dark-matter substructure (Mao & Schneider, arXiv:astro-ph/9707187) or dust.
- **Two-epoch brightness changes:** not diagnostic of exotic lenses. They are degenerate with microlensing by
  stars, intrinsic variability with time delays, and supernovae.

## Existing observational limits

- Takahashi & Asada, arXiv:1303.1301 (ApJL 768, L16): no exotic lensing among about 50,000 SDSS quasars.
  - Negative-mass objects: n < 10⁻⁸ h³ Mpc⁻³ above 10¹⁵ M☉, and n < 10⁻⁴ h³ Mpc⁻³ above 10¹² M☉.
  - Ellis wormholes: n < 10⁻⁴ h³ Mpc⁻³ for throat radii 10–10⁴ pc.
- Anchordoqui et al., arXiv:astro-ph/9904399, "In search for natural wormholes": 631 BATSE gamma-ray-burst
  profiles.

**Consequence for this project.** The realistic outcome is a null result with better limits, which is still
publishable science. The search therefore tests lensing consistency in cluster fields against published models:
- First, ordinary checks: image positions, parities, flux ratios, and arc orientation against the predicted
  shear.
- Only then the exotic-specific patterns from the table, after ordinary explanations fail (TASKS "Now (M3)").

## Screens and results (D-031)

`scripts/exotic_screens.py fluxratio` (demagnification) and `radial` (radial arcs around a dark centre). Every
hit goes to `/vet-candidate`.

| Field | fluxratio (compact images) | radial (convergence centres) | Verdict |
|---|---|---|---|
| WHL0137 / Sunrise (RELICS Lenstool maps, D-033) | no image list | 3 against a null mean of 2.2 (p95 5); max 3 lines, p 0.885 | null |
| SMACS 0723 (ICLv2) | 6 compact images consistent, 0 flags (49 resolved) | 12 against a null mean of 9.0 (p95 15); max 4 lines, p 0.945 | null |
| MACS1149 (HFF CATS, D-037) | 19 consistent, 0 flags (106 resolved; no photo-z) | with CANUCS photo-z: 12 against a null mean of 12.8 (p95 20); max 3 lines, p 1.0 (without: 32 against 30.8, p 0.91) | null |
| MACS0416 (HFF CATS, D-038) | not run (image list gated) | with CANUCS photo-z: 77 against a null mean of 85.0 (p95 103); max 7 lines, p 0.225 (without: 98 against 99.6, p 0.29) | null |
| MACS0717 (HFF CATS, D-039) | 26 consistent; system 65 flags explained (catalogue offset; μ model-dependent) | 11 against a null mean of 16.7; max 5 lines, p ≥ 0.70 (no photo-z) | null |
