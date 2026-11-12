# Changelog

All notable campaigns, methodological innovations, falsifications, and repository milestones are logged chronologically below.

## 2026-11-12: Phase 2 Ethos Capture and Convergent Technical Resolution (`Foundation-0/exoplanets-frontier`)

- **Organization migration**: Repository transferred into `Foundation-0/exoplanets-frontier` for shared multi-collaborator astronomy work alongside `Foundation-0/astr0`, `Foundation-0/gaia-dr4`, `Foundation-0/lost-eclipses`, `Foundation-0/cosmic-denoising`, and `Foundation-0/certified-orbits`.
- **Ethos and epistemic documentation**: Added `METHODS.md`, `CLAIMS.md`, `HANDOFF.md`, and `CHANGELOG.md` capturing all algorithms, five verification gates (`astr0` protocol), stalled technical bottlenecks, negative/falsified results, and cross-repo handoff entry points.
- **Convergent engine (`despot.py`, `ttv_track.py`, `blend_gaia.py`, `converge_io.py`, `converge.py`, `.github/workflows/converge.yml`)**:
  - Implemented `harmonic_despot`: per-segment periodogram + asymmetric Huber IRLS 4-harmonic Fourier starspot pre-whitener validated on known active M-dwarf planet host TOI-540 (`TIC 200322593`) and paired with time-reversed/sector-permuted null tests (`results/benchmark_and_null.json`).
  - Implemented wide-window (`+-240 min`) 8-year S1-S107 multi-sector transit timing extractor and 4-model $\Delta\text{BIC}$ ephemeris/TTV solver (`extract_transit_times`, `fit_ephemeris_models`), resolving the Phase 1 `+-60 min` grid-clipping bottleneck and enabling de-drifted multi-planet sibling searches (`results/ttv_ephemeris_audit.csv`, `results/sibling_residual_search.csv`).
  - Implemented proper-motion-propagated Gaia DR3 scene registration and multi-sector FITS WCS difference-image blend localizer (`multisector_wcs_blend_audit`), eliminating static neighbor dilution bias on `likely nearby FP` and `Tier C` repeaters (`results/blend_localization.csv`).
  - Added SHA-256 verification certificate output (`results/convergence_certificates.json`) and finite-guarded diagnostic plots in `results/figures/`.

## 2026-11-11: Phase 1 TESS Sectors 1-107 Blind Verification and Multi-Sector Dwarf Campaign

- **Modular MAST pipeline (`bls_search.py`, `vetter.py`, `mast_data.py`, `campaigns.py`, `deep_search.py`, `ledger.py`, `frontier.py`, `mdwarf_hunt.py`)**:
  - Indexed NASA STScI MAST SPOC and `Tesscut` FFI observations through Sector 107.
  - Executed out-of-sample Sector 106 and 107 blind verification across 39 (`eumemic/exoplanets` v3.0, S106/S107) pairs (`results/blind_verification_s106_s107.csv`): 22 `CONFIRMED_ON_SCHEDULE` ($\text{SNR} \ge 3.0$), 1 `TIMING_OFFSET_TTV` (`TIC 70557792`), 5 `MARGINAL_RECOVERY`, 9 `LOW_SNR_INCONCLUSIVE`, and 0 `NON_DETECTION_4SIGMA`. Reported upstream in `eumemic/exoplanets#20`.
  - Stacked 156 uncatalogued cool dwarfs and USP sibling hosts across up to 18 SPOC sectors (`results/discovery_triggers_vetted.csv`); vetted and falsified 19 $\text{MES} \ge 9.5$ triggers, including resolving cross-sector `SAP_FLUX` normalization (`norm_sap`) and starspot modulation (`spot_mod`) on `TIC 31281821`.
