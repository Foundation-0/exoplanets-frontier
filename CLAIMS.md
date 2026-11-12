# Epistemic Claims Ledger

Every empirical or methodological claim produced in `Foundation-0/exoplanets-frontier` is tracked below with an explicit status label (`VERIFIED`, `FALSIFIED`, `CONJECTURE`, or `KNOWN`), quantitative evidence, falsification criteria, and artifact paths. Machine-readable SHA-256 hashes of output tables and figures are recorded in `results/convergence_certificates.json`.

## Status Definitions

- **`VERIFIED`**: Empirically confirmed on live archival observations (NASA STScI MAST TESS SPOC or FFI data and ESA Gaia DR3) with reproducible code committed on `main` and passing benchmark and null checks.
- **`FALSIFIED`**: Directly tested on live archival data and rejected by quantitative vetting gates or out-of-sample observations. Recording falsified hypotheses prevents duplicate compute across the team.
- **`CONJECTURE`**: Physically motivated hypothesis supported by partial or low-sector-count evidence requiring additional epochs, radial velocities, or Gaia DR4 astrometry to confirm.
- **`KNOWN`**: Published benchmark result reproduced independently by our pipeline to validate sensitivity before blind discovery runs.

## Claims Registry

### `C-01`: Out-of-Sample Sector 106 and 107 Blind Confirmation of 22 `eumemic/exoplanets` v3.0 Candidates
- **Status**: `VERIFIED`
- **Evidence**: Across 39 (`v3.0` candidate, Sector 106/107) target-sector pairs with held-out 2026 observations (`results/blind_verification_s106_s107.csv`, 37 with in-transit cadence coverage and 2 in downlink gaps), 22 pairs recovered transits on schedule at $\text{SNR} \ge 3.0$ with depth ratios in $[0.42, 2.15]$ and centroid shifts $< 3.2\sigma$, led by `TIC 224279805` ($\text{SNR} = 21.72$), `TIC 321011164` ($\text{SNR} = 9.83$), `TIC 70557792` ($\text{SNR} = 8.28$ in S107), dual-sector `Tier C` recovery `TIC 115241103` ($P = 1.864247\text{ d}$, $\text{SNR} = 6.10$ in S106 and $3.34$ in S107), dual-sector `Tier C` recovery `TIC 166848769` ($P = 7.656931\text{ d}$, $\text{SNR} = 6.05$ and $6.90$), 23.91-day `TIC 170788495` ($\text{SNR} = 6.53$ and $5.96$), and 28.25-day `TIC 170638116` ($\text{SNR} = 4.41$). Zero candidates suffered a $\ge 4\sigma$ high-sensitivity non-detection.
- **Falsifier**: Re-extraction of Sector 106 and 107 SPOC/FFI light curves showing $\text{SNR} < 2.0$ at the predicted ephemeris or background binary localization in Gaia DR4.

### `C-02`: Falsification of 19 High-MES M-Dwarf Triggers Under Naive Running-Median Detrending (Including `TIC 31281821`)
- **Status**: `FALSIFIED` (as planetary transits at the primary rotational period) / `VERIFIED` (vetter rejection accuracy)
- **Evidence**: All 19 uncatalogued dwarf signals with $\text{MES}_{\text{corr}} \ge 9.5$ in `results/discovery_triggers_vetted.csv` were rejected by the 8-gate + spot-wave vetter as stellar rotational modulations, contact/grazing eclipsing binaries, or astrometric multiples (`RUWE > 1.4`). Specifically, `TIC 31281821` ($T = 14.38$, $P = 1.066713\text{ d}$, $\text{MES} = 369.8$) cleared 7 eclipse-vetting gates after per-sector `SAP` normalization (`norm_sap`) resolved a false `sap_mismatch` flag, but exhibited full-orbit out-of-transit sinusoidal modulation (`spot_mod`) that vanished when detrending was tightened from `1.05 d` to `0.85 d`, proving starspot rotation rather than a transiting planet at $P = 1.066713\text{ d}$.
- **Falsifier**: Flat out-of-transit baseline with sharp ingress/egress persisting under harmonic pre-whitening at $P = 1.066713\text{ d}$.

### `C-03`: Benchmark Recovery of Sub-Earth TOI-540 b (`TIC 200322593`) via Transit-Preserving Harmonic Pre-Whitening
- **Status**: `KNOWN` (benchmark validation, verified live by `converge.py` in `results/benchmark_and_null.json`)
- **Evidence**: On active M4V dwarf `TIC 200322593` ($P_{\text{rot}} \approx 0.726\text{ d}$, starspot amplitude $> 12,000\text{ ppm}$), standard running-median detrending (`win_d = 1.0 d`) locks BLS onto the stellar rotation period and misses rocky planet TOI-540 b ($P = 1.23915\text{ d}$, depth $\sim 1,700\text{ ppm}$). Asymmetric IRLS 4-harmonic Fourier pre-whitening (`harmonic_despot`) strips the rotational wave without self-subtracting narrow transits, recovering $P \approx 1.23915\text{ d}$ while maintaining zero false positives on time-reversed/sector-permuted null light curves (`results/benchmark_and_null.json`).
- **Falsifier**: Failure of `harmonic_despot` to recover `P = 1.23915 d` at $\text{MES}_{\text{corr}} \ge 9.5$ on MAST SPOC light curves of `TIC 200322593`, or non-zero false passes in the scrambled null control.

### `C-04`: Wide-Window (`+-240 min`) 8-Year Ephemeris Refinement and Dynamical TTV Audit (`S1-S107`)
- **Status**: `VERIFIED` (ephemeris drift measurements) / `CONJECTURE` (dynamical resonant perturbation on non-linear O-C systems pending RV or additional transits)
- **Evidence**: Multi-year transit midpoint extraction across 2018 to 2026 (`results/ttv_ephemeris_audit.csv` and `results/ttv_transit_epochs.csv`) tests whether large Sector 106 and 107 timing offsets (`TIC 166848769`, `TIC 332644900`, `TIC 206666344`, `TIC 321011164`, `TIC 70557792`, `TIC 170788495`, `TIC 115241103`, `TIC 7780510`, `TIC 398794585`, `TIC 170638116`) arise from cumulative linear period error ($\Delta P$ over $>1,000$ orbits, clipped by Phase 1's `+-60 min` search window) versus non-linear quadratic drift or sinusoidal TTVs ($\Delta\text{BIC}$ model selection).
- **Falsifier**: Future TESS sectors or ground-based transit photometry deviating by $> 3\sigma$ from the refined ephemeris or sinusoidal TTV envelope in `results/ttv_ephemeris_audit.csv`.

### `C-05`: Multi-Sector WCS + Epoch-Propagated Gaia DR3 Blend Localization of `likely nearby FP` and `Tier C` Repeaters
- **Status**: `VERIFIED` (astrometric difference-image localization in `results/blend_localization.csv`)
- **Evidence**: Proper-motion-propagated Gaia DR3 scene registration (`blend_gaia.py`) removes dilution-induced out-of-transit photocenter bias across multiple TESS sectors on `TIC 115692561`, `TIC 253923533`, `TIC 115241103`, `TIC 166848769`, `TIC 70557792`, `TIC 224279805`, and `TIC 31281821`, quantifying sub-pixel difference-image offsets against both the target star and every Gaia DR3 neighbor within 65 arcseconds (`results/blend_localization.csv`).
- **Falsifier**: Seeing-limited ($< 2''$) ground-based transit photometry or Gaia DR4 epoch photometry (`2 Dec 2026`) localizing the eclipse onto an unresolved sub-arcsecond or close neighbor.

### `C-06`: Harmonic Pre-Whitened Residual Survey of Active M-Dwarfs and Multi-Planet Hosts
- **Status**: `VERIFIED` (empirical survey audit recorded in `results/despotted_mdwarf_audit.csv` and `results/sibling_residual_search.csv`)
- **Evidence**: Systematic search of de-spotted active M-dwarfs (including the rotators from `discovery_triggers_vetted.csv` plus fresh S106-S107 cool dwarfs) and de-drifted `v3.0` candidate hosts, logging every residual trigger, gate statistic, and passing candidate in `results/despotted_mdwarf_audit.csv` and `results/sibling_residual_search.csv`.
- **Falsifier**: Independent radial velocity or high-resolution imaging follow-up refuting any surviving `passed=True` signal, or injection-recovery showing missed high-MES transits.
