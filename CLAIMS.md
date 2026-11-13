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
### C-07: Falsification of All 6 Issue #5 Bright Temperate Leads (0/6 Survived) and Isolation of Inverted-Flux Control Leak
- **Status**: FALSIFIED (as planet candidates; 0/6 survived) / VERIFIED (blind multi-sector physics rejection and provenance isolation)
- **Evidence**: Live S1-S108 multi-sector coherent BLS, two-pass starspot pre-whitening, 8-gate transit symmetry vetting, and Gaia DR3 + FITS WCS pixel localization (hz_core.py, hz_vet.py, commit 883e28a, results/hz_issue5_audit.csv SHA-256 28e36fbd280ddefc69a3adb9e4089e99a034c395230a16c7c39a9331781dd14c) rejected 6 out of 6 bright temperate leads (TIC 181015118 at 25.262 d and 18.318 d; active ~150 Myr AB Doradus M-dwarf TIC 246897668 at 24.752 d; TIC 37904598 at 17.049 d; TIC 146636926 at 26.056 d and 23.920 d). All four hosts are rapidly rotating active M-dwarfs (P_rot = 0.939-3.830 d) whose apparent signal collapsed below threshold after starspot removal (despot_MES = 3.77-6.58), failed odd/even and inverted/out-of-transit symmetry gates, and localized 17.7 arcsec off-target onto a neighboring background blend on TIC 37904598 (NEIGHBOR_BEB_LOCALIZED). Full-tree CSV provenance scanning (scan_eumemic_origins) independently proved all 6 rows originated solely in validation/inverted_lightcurves_vetted.csv (an upside-down light-curve negative-control table), confirming 100% blind rejection of leaked negative controls while recovering genuine on-target multi-sector warm giant transits on TIC 166848769 (despot_MES = 11.83, 0.37 arcsec ON_TARGET_CLEAN) and TIC 224279805 (despot_MES = 29.16, 2.56 arcsec ON_TARGET_CLEAN).
- **Falsifier**: Multi-sector SPOC/FFI photometry showing flat out-of-transit baseline and MES_corr >= 9.5 symmetric U-shaped transits at any of the 6 periods.
### C-08: Purification of Pre-Gaia-DR4 Candidate Catalog (2,421 Genuine Candidates; 348 Platinum and Gold Clean Targets)
- **Status**: VERIFIED (live VizieR I/355/gaiadr3 + MAST TIC rebuild in results/GAIA_DR4_CATALOG_AUDIT.md, commit 4059fc1)
- **Evidence**: Replacing alphabetical multi-CSV deduplication with priority candidate-only ingestion (dr4_core.py, refresh_dr4.py) eliminated all 1,122 leaked validation and negative-control rows (0 rows remaining from validation/inverted_lightcurves_vetted.csv; 0 overlap with Issue #5 stars) and restored true upstream tiers across all 2,421 genuine candidates in results/eumemic_gaia_dr4_readiness.csv (SHA-256 f332cf6c2140de2cd9e3c9eaf443d1d1d9eada981cbdbfbc57df3b9471361746: 291 Tier S, 551 Tier C, 19 reported, 1,109 likely nearby FP, and 451 Tier N/L1/L/M/legacy). Scaling 8-worker VizieR Gaia DR3 cone checks across 2,210 unique candidate fields (91.3% live Gaia DR3 coverage) while propagating multi-sector binary/blend flags expanded results/platinum_followup_targets.csv (SHA-256 e2bca4f9e22160d1676bf2678f20bfbd1acc29d022e5fff522401c0ad49161b5) from 22 to 348 genuine low-blend planet candidates (137 PLATINUM_SOLO_VERIFIED + 211 GOLD_LOW_BLEND), while confirming 0 small habitable-zone (180-325 K, <= 3.2 Re) candidates exist in the upstream FFI catalog (coolest clean candidates: TIC 399695544 at 426 K and TIC 409696823 at 453 K).
- **Falsifier**: Gaia DR4 epoch astrometry or photometry (2 Dec 2026) revealing unresolved sub-arcsecond binary companions or background eclipsing blends among the 137 PLATINUM_SOLO_VERIFIED hosts.
