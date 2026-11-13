# Epistemic Claims Ledger

Every empirical or methodological claim produced in `Foundation-0/exoplanets-frontier` is tracked below with an explicit status label (`VERIFIED`, `FALSIFIED`, `CONJECTURE`, or `KNOWN`), quantitative evidence, falsification criteria, and artifact paths pinned by SHA-256 hashes in `results/convergence_certificates.json` (commit `412043e`).

## Status Definitions

- **`VERIFIED`**: Empirically confirmed on live archival observations (NASA STScI MAST TESS SPOC or FFI data) with reproducible code committed on `main` and passing benchmark and null checks.
- **`FALSIFIED`**: Directly tested on live archival data and rejected by quantitative vetting gates or out-of-sample observations. Recording falsified hypotheses prevents duplicate compute across the team.
- **`CONJECTURE`**: Physically motivated hypothesis supported by partial or low-sector-count evidence requiring additional epochs, radial velocities, or Gaia DR4 astrometry to confirm.
- **`KNOWN`**: Published benchmark result reproduced independently by our pipeline to validate sensitivity before blind discovery runs.

## Claims Registry

### `C-01`: Out-of-Sample Sector 106 and 107 Blind Confirmation of 22 `eumemic/exoplanets` v3.0 Candidates
- **Status**: `VERIFIED`
- **Evidence**: Across 39 (`v3.0` candidate, Sector 106/107) pairs (`results/blind_verification_s106_s107.csv`, 37 with in-transit coverage), 22 recovered transits on schedule at $\text{SNR} \ge 3.0$, led by `TIC 224279805` ($\text{SNR} = 21.72$), `TIC 321011164` ($\text{SNR} = 9.83$), `TIC 70557792` ($\text{SNR} = 8.28$ in S107), dual-sector `Tier C` repeaters `TIC 115241103` ($\text{SNR} = 6.10, 3.34$) and `TIC 166848769` ($\text{SNR} = 6.05, 6.90$), `TIC 170788495` ($P = 23.91\text{ d}$), and `TIC 170638116` ($P = 28.25\text{ d}$), with zero $4\sigma$ non-detections.
- **Falsifier**: Re-extraction of Sector 106 and 107 photometry showing $\text{SNR} < 2.0$ or off-target eclipse localization in Gaia DR4.

### `C-02`: Falsification of 19 Naive-Detrended M-Dwarf Triggers and 24 Harmonic-Prewhitened M-Dwarf Residuals
- **Status**: `FALSIFIED` (no transiting planet with $\text{MES}_{\text{corr}} \ge 9.5$ in this 24-star active/cool M-dwarf sample) / `VERIFIED` (spot-suppression and vetter accuracy)
- **Evidence**: All 19 initial $\text{MES} \ge 9.5$ triggers in `results/discovery_triggers_vetted.csv` were stellar rotators, contact binaries, or high-RUWE multiples. Re-analysing 10 low-RUWE rotators plus 14 fresh S106-S107 cool M-dwarfs with `harmonic_despot` (`results/despotted_mdwarf_audit.csv`) suppressed spot waves (for example `TIC 627018` RMS `10,662 -> 3,786 ppm`, `TIC 1439071` `8,662 -> 3,630 ppm`, and `TIC 31281821` trigger dropping from `MES = 369.8` to `5.57`), leaving all 24 residual peaks at `MES = 4.66 to 6.78`, consistent with the scrambled Gaussian floor.
- **Falsifier**: Detection of a coherent box transit at $\text{MES}_{\text{corr}} \ge 9.5$ on any of these 24 TICs after two-pass unmasked+masked spot subtraction.

### `C-03`: Blind Benchmark Recovery of Sub-Earth TOI-540 b (`TIC 200322593`) and Zero-False-Alarm Scrambled Null Control
- **Status**: `KNOWN` / `VERIFIED` (`results/benchmark_and_null.json`, SHA-256 `6c2a1823...`)
- **Evidence**: On SPOC 2-minute light curves of active M4V star `TIC 200322593` ($P_{\text{rot}} = 0.726\text{ d}$), standard running-median BLS locked onto the starspot wave (`naive_P_d = 0.726307 d`), whereas `harmonic_despot` blindly recovered rocky planet TOI-540 b at `P = 1.239126 d` (`1.9 s` from published `1.239149 d`) with `MES = 14.33` and depth `1,374.5 ppm`. Ten time-reversed, sector-permuted null trials produced `0/10` false alarms above `MES 9.5` (`max MES = 5.77`).
- **Falsifier**: Non-zero false-alarm rate above `MES 9.5` on expanded 100-trial scrambled controls.

### `C-04`: Rescue of Four Phase-Clipped Candidates via Wide-Window (`+-240 min`) Multi-Epoch Ephemeris Refinement
- **Status**: `VERIFIED` (`results/ttv_ephemeris_audit.csv`, SHA-256 `a81f7b1c...`)
- **Evidence**: Widening transit midpoint extraction from `+-60 min` to `+-240 min` resolved cumulative period drift that had depressed fixed-ephemeris SNR in Phase 1:
  - `TIC 166848769` (`Tier C`, 8 transits): `P_ref = 7.6574415 d` (`dP = +44.11 s`), reducing timing RMS from `59.8 min` to `10.5 min` (`dBIC_lin = +212.5`, `dBIC_ttv = -3.4`) and doubling stacked SNR from `6.29` to `13.49`.
  - `TIC 332644900` (`Tier C`, 8 transits): `P_ref = 1.3290094 d` (`dP = -5.24 s`), reducing RMS from `79.7 min` to `30.9 min` (`dBIC_lin = +136.0`) and doubling stacked SNR from `2.67` (previously `LOW_SNR_INCONCLUSIVE`) to `5.39`.
  - `TIC 70557792` (`Tier S`, 8 transits including SPOC S30): `P_ref = 7.0488122 d` (`dP = -1.45 s`) improved fit by `dBIC_lin = +21.0` (`SNR 9.46 -> 10.22`), while sinusoidal TTV added only `dBIC_ttv = +4.4`, reclassifying our earlier 2-sector TTV flag primarily as linear period refinement.
  - `TIC 321011164` (`SNR 7.07 -> 10.24`, `dBIC_lin = +223.0`) and `TIC 7780510` (`SNR 0.84 -> 9.13`, `dBIC_lin = +179.6`) recovered strongly, though 2-sector FFI sampling leaves discrete cycle-count alias degeneracy requiring intermediate Cycle 1-7 FFI cutouts.
- **Falsifier**: Sector 108+ or intermediate Cycle 2-6 FFI transits missing the interpolated linear ephemeris by $> 3\sigma$.

### `C-05`: Candidate Non-Linear Timing Variation and Residual 14.35-Day Signal on `TIC 332644900`
- **Status**: `CONJECTURE` (`results/ttv_ephemeris_audit.csv`, `results/sibling_residual_search.csv`)
- **Evidence**: After removing `-5.24 s/orbit` linear drift on `TIC 332644900`, the 8 measured transit epochs retain `30.9 min` scatter with `dBIC_ttv = +18.7`, and the primary-masked residual light curve shows a `P = 14.3538 d` (`MES = 11.41`, single-event peak `SNR = 13.83` at `BTJD 4262.8802`) trigger. Low per-transit SNR on this `455 ppm` FFI host could mimic timing scatter; intermediate FFI sectors are required before claiming a resonant multi-planet system.
- **Falsifier**: Full-sector FFI extraction showing purely Gaussian timing residuals (`dBIC_ttv < 6`) or spacecraft scattered-light coincidence at `BTJD 4262.88`.

### `C-06`: WCS Sub-Pixel On-Target Localization of `likely nearby FP` Repeater `TIC 253923533` and Cloud Gaia TAP Finding
- **Status**: `VERIFIED` (`results/blend_localization.csv`, SHA-256 `7f735d0b...`)
- **Evidence**: FITS `WCS` difference-image centroiding localized `TIC 224279805` to `0.55 +/- 2.20 arcsec` (`0.026 px`, `ON_TARGET_CLEAN`) and upstream `likely nearby FP` repeater `TIC 253923533` ($P = 1.156013\text{ d}$) to `4.90 +/- 2.20 arcsec` (`0.23 px`, `ON_TARGET_CLEAN`), supporting promotion of `TIC 253923533` pending high-resolution follow-up, while `TIC 31281821` shifted by `24.06 +/- 1.56 arcsec` (`1.15 px`). Live execution also demonstrated that ESA's synchronous Gaia TAP server drops/throttles automated cloud-runner queries (`gaia_sid = none`), motivating a CDS VizieR (`I/355/gaiadr3`) + MAST TIC mirror query layer for `Foundation-0/gaia-dr4` ahead of the 25 Nov 2026 ESA maintenance cutoff.
- **Falsifier**: Sub-arcsecond ground-based transit imaging or Gaia DR4 epoch photometry resolving an eclipsing binary within `< 5 arcsec` of `TIC 253923533`.

### C-07: Multi-Sector Screen Rejection of Six Issue #5 Signals and Inverted-Control Table Matches

- **Status**: Recorded screen result: 0/6 signal rows on four stars survived. Their physical classification is not established by this result alone.
- **Evidence**: `results/hz_issue5_audit.csv` at analysis revision `883e28a7cecde05fdb6679bf574dbe78668332e7` records six rejected Issue #5 rows, with three to seven loaded sectors per row. Its SHA-256 is `28e36fbd280ddefc69a3adb9e4089e99a034c395230a16c7c39a9331781dd14c`. The seeded search refits period within 1.8 percent of each proposed value, plus epoch and duration. This is reanalysis, without a demonstrated held-out or blinded design.
- **Control-table provenance**: Every row records `validation/inverted_lightcurves_vetted.csv` in `src_files`. All six TIC/period pairs match that upstream table at revision `e2eee928fd621f75c1621f6b4206d34bb222ae59`. Upstream code explicitly presents the TIC 37904598 row as an inverted-light-curve false alarm. These are control-table matches, not established candidate discoveries.
- **Limits**: The scanner records matching paths but does not pin the historical upstream checkout or retain each row's control and vetting context. It therefore does not prove exclusive historical origin or a calibrated 100 percent rejection rate. The present run analyzes ordinary light curves, rather than replaying the original inverted controls.
- **Comparison signals**: The periodic dimming recovered on TIC 166848769 and TIC 224279805 also fails the configured screen, including radius and other vetoes. Those fits do not independently establish genuine planets.
- **Next check**: Preserve the upstream revision, exact signal rows and control labels, then calibrate recovery and rejection using separately identified planet and nonplanet controls through the intended analysis. No planet-absence or survey-completeness claim follows from the current zero-survivor result.
