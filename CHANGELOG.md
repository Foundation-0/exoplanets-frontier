# Changelog

All notable campaigns, methodological innovations, falsifications, and repository milestones are logged chronologically below.

## 2026-11-13: Phase 3 Unified Sectors 1-108 Run and Pre-Gaia-DR4 Catalog (`9266751`)

- **Live TESS Sector 108 capture (`results/convergence_certificates.json`)**: Automatically indexed newly released STScI MAST Sector 108 (`s_max = 108`, `1,109.0 s` runtime in Run `#43343464806`) after replacing `Observations.get_product_list` with direct deterministic MAST product URIs (`converge_io.py`, `5652d80`).
- **Pre-Gaia-DR4 catalog (`catalog_665.py`, `results/eumemic_gaia_dr4_readiness.csv`, `results/platinum_followup_targets.csv`)**: Cross-matched 3,543 candidate entries from `eumemic/exoplanets` against MAST TIC v8.2 and CDS VizieR `I/355/gaiadr3` ahead of the 25 Nov 2026 ESA archive maintenance window, isolating 22 Platinum and Gold low-blend targets (including 6 bright `180-303 K` temperate candidates led by `TIC 181015118` and `TIC 246897668` queued for full SPOC transit-shape vetting).
- **Full Cycle 1-8 ephemeris lock and VizieR WCS pixel audit (`results/CONVERGENCE_REPORT.md`)**:
  - Rescued `Tier C` candidate `TIC 166848769` (`Gaia DR3 4858465154105701120`, `RUWE = 1.00`) across 8 available sectors (`11` transits): `P_ref = 7.6568169 d` (`dP = -9.86 s`), raising stacked SNR to `15.73` (`dBIC_lin = +190.9`) and localizing the transit to `0.16 +/- 1.16 arcsec` (`0.008 px`, `ON_TARGET_CLEAN`, nearest plausible neighbor `40.62 arcsec` away).
  - Resolved `TIC 31281821` as an off-target eclipsing binary blend (`50.90 +/- 1.56 arcsec` shift onto a neighbor at `40.62 arcsec`, `NEIGHBOR_BEB_LOCALIZED`), and uncovered close unresolved Gaia DR3 neighbors on `TIC 332644900` (`2.89 arcsec`) and `TIC 253923533` (`3.25 arcsec`).
- **Two-pass cool-dwarf survey (`despot.py`, `PLAIN_GUIDE.md`)**: Audited 33 cool dwarfs through Sector 108 (`0` new signals above `MES 9.5`; all residuals bounded at `MES <= 6.78`), and added `PLAIN_GUIDE.md` explaining all methods and results in everyday language.

## 2026-11-12: Phase 2 Ethos Capture and Convergent Technical Resolution (`412043e`)

- Added `METHODS.md`, `CLAIMS.md`, `HANDOFF.md`, `despot.py`, `ttv_track.py`, and `blend_gaia.py`; verified blind recovery of sub-Earth TOI-540 b (`P = 1.239126 d`, `MES = 14.33`) with `0/10` scrambled null false alarms.

## 2026-11-11: Phase 1 TESS Sectors 1-107 Blind Verification (`5d473fe`)

- Tested 39 held-out Sector 106 and 107 pairs (`22 CONFIRMED_ON_SCHEDULE`, `0 NON_DETECTION_4SIGMA`, reported in `eumemic/exoplanets#20`) and vetted 156 multi-sector cool dwarfs.
