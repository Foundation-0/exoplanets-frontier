# Issue 5 Habitable-Zone Multi-Sector Vetting Report (1396.0 s)

- 6 Bright Temperate Leads: 0/6 survived all S1-S108 transit-shape, starspot, and WCS gates.
- Verified/Cool Tier S/C Benchmarks: 0/5 survived.

## 1. Issue 5 Bright Temperate Leads (TIC 181015118, 246897668, 37904598, 146636926)

| tic | lane | src_files | known_nasa | Tmag | P_cand_d | n_sec_used | n_transits | p_rot_d | std_MES | despot_MES | Rp_Re | Teq_K | wcs_off_as | plausible_bebs | survived_hz_planet | vetoes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 181015118 | ISSUE5_HZ_PLATINUM | validation/inverted_lightcurves_vetted.csv | False | 8.99 | 25.26197 | 6 | 6 | 0.939 | 8.43 | 6.08 | 4.19 | 258.0 | 3.26 | 0 | False | low_mes(8.4);low_sde(3.3);odd_even(9.53s);inv_lc(9.9);sec_ecl(39.9);sap_mismatch |
| 181015118 | ISSUE5_HZ_GOLD | validation/inverted_lightcurves_vetted.csv | False | 8.99 | 18.31799 | 6 | 9 | 0.939 | 8.97 | 5.57 | 4.01 | 286.0 | 7.6 | 0 | False | low_mes(9.0);low_sde(3.6);odd_even(7.68s);inv_lc(10.2);sec_ecl(34.2);sap_mismatch |
| 246897668 | ISSUE5_HZ_PLATINUM | validation/inverted_lightcurves_vetted.csv | False | 9.24 | 24.75243 | 6 | 6 | 1.3774 | 7.57 | 3.77 | 3.85 | 282.0 | 5.87 | 0 | False | low_mes(7.6);low_sde(3.1);inv_lc(8.3);sec_ecl(31.5);sap_mismatch |
| 37904598 | ISSUE5_HZ_GOLD | validation/inverted_lightcurves_vetted.csv | False | 11.21 | 17.04918 | 3 | 5 | 3.8302 | 4.08 | 4.46 | 1.54 | 303.0 | 17.7 | 3 | False | low_mes(4.5);low_sde(4.4);inv_lc(25.8);tpf(17.7as);wcs(NEIGHBOR_BEB_LOCALIZED) |
| 146636926 | ISSUE5_HZ_GOLD | validation/inverted_lightcurves_vetted.csv | False | 11.47 | 26.05645 | 7 | 5 | 2.3721 | 14.39 | 6.58 | 5.3 | 249.0 | 1.17 | 0 | False | low_sde(2.9);odd_even(3.33s);sec_ecl(43.1);sap_mismatch |
| 146636926 | ISSUE5_HZ_GOLD | validation/inverted_lightcurves_vetted.csv | False | 11.47 | 23.92026 | 7 | 7 | 2.3721 | 14.96 | 6.23 | 4.98 | 258.0 | 1.26 | 0 | False | low_sde(4.3);odd_even(11.70s);sec_ecl(62.1);sap_mismatch |

## 2. Cool Tier S/C & Rescued Comparison Worlds

| tic | lane | src_files | known_nasa | Tmag | P_cand_d | n_sec_used | n_transits | p_rot_d | std_MES | despot_MES | Rp_Re | Teq_K | wcs_off_as | plausible_bebs | survived_hz_planet | vetoes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 166848769 | VERIFIED_COOL_TIER_C | candidates/ephemerides.csv;candidates/predictions.csv;results_public/candidates_allsky.csv;results_public/products/fits.csv;results_public/products/measurements_allsky.csv;results_public/products/triceratops.csv;results_public/products/triceratops_plan.csv;results_public/products/vetting.csv;validation/gaia_companion_limits.csv | True | 13.47 | 7.65682 | 4 | 13 | 5.4441 | 11.58 | 11.83 | 8.94 | 854.0 | 0.37 | 1 | False | rp(8.94Re);inv_lc(11.2);sec_ecl(4.3) |
| 224279805 | VERIFIED_COOL_REPORTED | candidates/ephemerides.csv;candidates/predictions.csv;results_public/candidates_allsky.csv;results_public/overrides.csv;results_public/products/fits.csv;results_public/products/measurements_allsky.csv;results_public/products/triceratops.csv;results_public/products/triceratops_plan.csv;results_public/products/vetting.csv | True | 12.15 | 22.98417 | 6 | 5 | 5.3663 | 33.02 | 29.16 | 11.05 | 661.0 | 2.56 | 1 | False | rp(11.05Re);odd_even(5.05s) |
| 170788495 | VERIFIED_COOL_TIER_S | candidates/ephemerides.csv;candidates/predictions.csv;results_public/candidates_allsky.csv;results_public/products/fits.csv;results_public/products/measurements_allsky.csv;results_public/products/triceratops.csv;results_public/products/triceratops_plan.csv;results_public/products/vetting.csv;validation/findability.csv;validation/gaia_companion_limits.csv | True | 11.59 | 23.90688 | 5 | 2 | 6.7332 | 19.4 | 16.21 | 6.71 | 482.0 | 4.92 | 1 | False | low_sde(3.6);odd_even(41.03s);sec_ecl(17.3);single_ev |
| 140575951 | V3_COOL_C | results_public/candidates_allsky.csv;results_public/products/measurements_allsky.csv;results_public/products/triceratops.csv;results_public/products/triceratops_plan.csv;results_public/products/vetting.csv | True | 12.55 | 21.565 | 6 | 9 | 6.336 | 7.48 | 7.45 | 4.1 | 410.0 | 9.98 | 6 | False | low_mes(7.4);low_sde(3.4);odd_even(7.47s);inv_lc(10.6);single_ev;cent(35.9as) |
| 31796864 | V3_COOL_S | candidates/ephemerides.csv;candidates/predictions.csv;results_public/candidates_allsky.csv;results_public/products/fits.csv;results_public/products/measurements_allsky.csv;results_public/products/spoc_2min.csv;results_public/products/triceratops.csv;results_public/products/triceratops_plan.csv;results_public/products/vetting.csv;validation/findability.csv;validation/spoc_2min_recheck.csv | True | 11.96 | 25.99096 | 7 | 7 | 6.1423 | 10.63 | 9.61 | 3.87 | 428.0 | 17.19 | 0 | False | low_sde(4.2);inv_lc(15.4);sec_ecl(5.4);single_ev;tpf(17.2as);wcs(OFFSET_MARGINAL) |
