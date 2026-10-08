"""Galaxy count maps on HEALPix pixels for the W5 count-deficit screen (D-063).

Work in progress: a thin :class:`signatures.CatalogueSurvey`-style adapter that aggregates the
Legacy Surveys DR10 Tractor catalogue per ``nest4096`` pixel server-side (NOIRLab Astro Data Lab
TAP), with per-pixel mask and depth from the same catalogue. Exotic physics is a hypothesis; a flag
is an anomaly to vet, and a null becomes an injection-calibrated limit.
"""
