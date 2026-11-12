# Team Handoff and Cross-Repository Integration Guide

Per our team collaboration ethos, every computation, method, verified result, and negative/falsified path in `Foundation-0/exoplanets-frontier` is structured for immediate continuation by any researcher or automated workflow across `Foundation-0` and `eumemic/exoplanets`.

## 1. The Three Ethos Questions

### What did we learn, decide, or verify?
1. **Out-of-sample confirmation works at scale**: Blind testing on newly released 2026 TESS Sectors 106 and 107 confirmed 22 candidate-sector transit recoveries on schedule (`0` high-sensitivity $4\sigma$ non-detections) across `eumemic/exoplanets` v3.0 candidates (`results/blind_verification_s106_s107.csv`).
2. **Why promising runs stalled on technicalities (and how they converge)**:
   - *Ephemeris window clipping (`+-60 min`)* hid multi-year period drift ($\delta P \sim 2\text{ to }15\text{ s}$ accumulated over $>1,000$ orbits between 2018 and 2026), artificially depressing fixed-ephemeris SNR on real repeaters (`TIC 166848769`, `TIC 332644900`, `TIC 206666344`, `TIC 321011164`) while `TIC 70557792` exhibited a sign-reversing timing shift. Wide-window (`+-240 min`) multi-sector timing (`ttv_track.py`) separates linear $\Delta P$ rescue from quadratic nodal precession and sinusoidal dynamical TTVs.
   - *Dilution photocenter bias* in 21-arcsecond TESS pixels caused single-sector centroid checks to penalize real on-target transits in crowded fields (`TIC 115692561` and `TIC 253923533`, both recovered on schedule in S106/S107 despite upstream `likely nearby FP` tags). Registering pixel coordinates to proper-motion-propagated Gaia DR3 scenes (`blend_gaia.py`) eliminates static neighbor pull.
   - *Active M-dwarf starspot waves* (`P_rot = 0.27 to 1.1 d`) blinded standard median-detrended BLS across all 19 initial M-dwarf triggers (`TIC 31281821` et al.). Asymmetric Huber 4-harmonic Fourier pre-whitening (`harmonic_despot` in `despot.py`) strips 1 to 3 percent starspot modulations without attenuating narrow planetary transits, verified on benchmark sub-Earth TOI-540 b (`TIC 200322593`).

### Where does it live now?
- **Core baseline modules**: `bls_search.py`, `vetter.py`, `mast_data.py`, `campaigns.py`, `deep_search.py`, `ledger.py`, `frontier.py`, `mdwarf_hunt.py`
- **Phase 2 convergent resolution engine**: `despot.py`, `ttv_track.py`, `blend_gaia.py`, `converge_io.py`, `converge.py`, `.github/workflows/converge.yml`
- **Documentation and epistemic state**: `README.md`, `METHODS.md`, `CLAIMS.md`, `HANDOFF.md`, `CHANGELOG.md`
- **Empirical data products and SHA-256 certificates**:
  - `results/blind_verification_s106_s107.csv` and `results/DISCOVERIES.md`
  - `results/benchmark_and_null.json` (`astr0` Rules 2 and 3 benchmark + scrambled null test)
  - `results/ttv_ephemeris_audit.csv` and `results/ttv_transit_epochs.csv` (8-year S1-S107 O-C timing & $\Delta\text{BIC}$ models)
  - `results/sibling_residual_search.csv` (de-drifted multi-planet sibling and single-transit scan)
  - `results/blend_localization.csv` (WCS + epoch-propagated Gaia DR3 sub-pixel blend audit)
  - `results/despotted_mdwarf_audit.csv` (harmonic pre-whitened M-dwarf planet rescue survey)
  - `results/CONVERGENCE_REPORT.md` and `results/convergence_certificates.json`
  - `results/figures/*.png` (multi-panel diagnostic phase folds, O-C curves, de-spotted light curves, and WCS difference images)

### Who can build on it next? (Direct Cross-Repo Entry Points)

1. **`Foundation-0/gaia-dr4` (Immediate Priority before 25 Nov 2026 ESA Archive Cutoff)**:
   - Import `gaia_epoch_scene` and `multisector_wcs_blend_audit` from `blend_gaia.py` to batch-audit all 665 `eumemic/exoplanets` candidates before the ESA Gaia archive goes offline on 25 November 2026, caching `source_id`, epoch-propagated pixel coordinates, `ruwe`, and `req_eb_depth` for instant join against Gaia DR4 epoch photometry on 2 December 2026.
2. **`Foundation-0/lost-eclipses` & `Foundation-0/certified-orbits`**:
   - Reuse `extract_transit_times` and `fit_ephemeris_models` in `ttv_track.py` to measure secular eclipse depth decay and quadratic orbital/nodal precession ($dP/dE$) across TESS Cycles 1 to 8 on GCVS disappearing-eclipse candidates, and feed verified non-linear TTV systems (`dBIC_ttv >= 10`) into N-body resonant orbit certification.
3. **`eumemic/exoplanets` Upstream Catalog (`v3.1` / Sector 108+)**:
   - Incorporate the 8-year refined linear ephemerides (`P_ref_d`, `T0_ref_BJD`) from `results/ttv_ephemeris_audit.csv` to rescue phase-drifted `Tier C` candidates (`TIC 166848769`, `TIC 115241103`, `TIC 332644900`), update dispositions on `likely nearby FP` targets audited in `results/blend_localization.csv`, and insert `harmonic_despot` upstream of FFI BLS on cool dwarfs ($T_{\text{eff}} < 4000\text{ K}$).
4. **Automated Sector 108+ Watch**:
   - Trigger `.github/workflows/campaign.yml` and `.github/workflows/converge.yml` whenever STScI MAST releases a new TESS SPOC/FFI sector to test predicted transit windows in `results/ttv_ephemeris_audit.csv`.
