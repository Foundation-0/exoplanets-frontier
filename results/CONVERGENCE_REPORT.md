# Phase 3 Synthesis (Purified): Genuine-Only Gaia DR4 Catalog, Cycle 1-8 Ephemerides & Two-Pass Cool-Dwarf Survey

Purified results_public/candidates_allsky.csv + VizieR I/355/gaiadr3 catalog rebuild (425.6 s; validation/inverted_lightcurves_vetted.csv purged).

- Genuine Candidate Catalog (results/eumemic_gaia_dr4_readiness.csv): 2421 candidates indexed ({'MULTI_SECTOR_VETTED_WARNING': 978, 'SILVER_CROWDED_SG1': 654, 'BINARY_OR_EVOLVED_WARNING': 230, 'GOLD_LOW_BLEND': 211, 'CROWDED_INDEXED': 191, 'PLATINUM_SOLO_VERIFIED': 137, 'LOW_CONTAM_INDEXED': 20}); 348 Platinum/Gold Clean Targets (0 temperate 180-325 K).

## 1. Top Purified Platinum Solo-Verified & Gold Low-Blend Candidates

| tic | upstream_tier | gaia_dr3_id | P_d | Rp_Re | Teq_K | is_hz_temperate | TSM | Tmag | ruwe | plausible_beb_65as | frontier_note | dr4_class |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 219329329 | S | 6507925182103208960 | 3.48114 | 11.2 | 1355.0 | False | 262.6 | 12.41 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 220162618 | S | 6230934043874288384 | 3.804428 | 12.09 | 1175.0 | False | 228.0 | 12.81 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 357504233 | S | 1264085364286632832 | 3.937343 | 11.83 | 1326.0 | False | 196.4 | 12.41 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 3886680 | S | 2524040649365516032 | 10.691905 | 10.49 | 899.0 | False | 120.7 | 12.68 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 224279805 | reported | 6534130392483208320 | 22.984167 | 11.21 | 661.0 | False | 147.0 | 12.15 | 1.0 | 1 | verified_on_target(4.78as) | GOLD_LOW_BLEND |
| 335955084 | C | 6305546460554805248 | 0.664381 | 5.68 | 1994.0 | False | 130.3 | 11.92 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 269746198 | C | 2851297159952973568 | 0.714152 | 8.16 | 2193.0 | False | 178.0 | 12.79 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 441301617 | C | 6793097197943841536 | 0.754367 | 8.21 | 1753.0 | False | 262.0 | 12.98 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 310261847 | C | 1300932163598098176 | 1.039797 | 12.25 | 1793.0 | False | 438.8 | 12.93 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 119968345 | C | 1907225724382651008 | 1.101738 | 8.21 | 1449.0 | False | 228.7 | 12.82 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 306484531 | C | 100461144774946944 | 1.772903 | 9.17 | 1813.0 | False | 146.0 | 12.9 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 293197619 | C | 2705701131576855040 | 1.799106 | 6.14 | 1864.0 | False | 146.9 | 10.98 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 305720469 | C | 6175362389744459520 | 1.806433 | 9.63 | 1739.0 | False | 146.1 | 12.8 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 437430304 | C | 6248561178072606592 | 3.418897 | 9.55 | 1470.0 | False | 266.5 | 11.14 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 241163081 | C | 927101882348867584 | 3.629423 | 11.1 | 1383.0 | False | 216.9 | 12.92 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 91051499 | C | 2844181326857165440 | 3.866608 | 9.31 | 1304.0 | False | 275.8 | 11.59 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 399695544 | C | 2562807711373323392 | 11.537094 | 8.66 | 426.0 | False | 124.7 | 13.37 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 355072530 | C | 272635972674562176 | 13.698274 | 9.09 | 780.0 | False | 154.8 | 12.66 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 409696823 | reported | 6456309089775016192 | 24.949846 | 6.47 | 453.0 | False | 89.4 | 11.7 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 79278415 | S | 6320381307658461184 | 3.348204 | 6.55 | 1328.0 | False | 88.0 | 12.03 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 156380270 | C | 3884882276685154944 | 1.521754 | 6.4 | 899.0 | False | 117.3 | 13.45 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 241537628 | S | 2481055723313271168 | 14.787554 | 7.76 | 825.0 | False | 84.1 | 12.13 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 187272399 | C | 3620212506276610432 | 1.762564 | 3.66 | 901.0 | False | 113.2 | 11.33 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 372946271 | C | 1759254067988060288 | 2.956309 | 8.47 | 1478.0 | False | 111.1 | 12.84 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 376665619 | C | 1321149880329990144 | 4.102 | 5.85 | 1188.0 | False | 107.6 | 11.21 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 46801759 | C | 5744530307567928448 | 4.690051 | 7.19 | 1232.0 | False | 106.9 | 12.53 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 437709086 | S | 2822231810591184896 | 2.925309 | 6.76 | 1504.0 | False | 76.2 | 12.63 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 270498731 | C | 3619863209471389696 | 0.898221 | 5.58 | 1916.0 | False | 101.8 | 12.57 | 1.0 | 1 | clean | GOLD_LOW_BLEND |
| 365061655 | C | 2566022064897681152 | 0.997125 | 6.3 | 1834.0 | False | 101.4 | 12.95 | 1.0 | 0 | clean | PLATINUM_SOLO_VERIFIED |
| 343942488 | S | 2823623208196929920 | 0.653612 | 8.84 | 2173.0 | False | 205.4 | 13.3 | 1.0 | 0 | clean | GOLD_LOW_BLEND |

## 2. Full 2018-2026 (Cycles 1-8) Ephemeris & TTV Lock

| tic | upstream_tier | n_sec_avail | n_sec_used | n_epochs_fit | P_orig_d | P_ref_d | dP_sec | rms0_min | rms_lin_min | dBIC_lin | dBIC_ttv | snr_stale | snr_refined | class |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 70557792 | S | 9 | 4 | 8 | 7.048829 | 7.0488256 | -0.292 | 38.42 | 24.4 | 29.5 | 4.0 | 9.34 | 9.81 | LINEAR_DRIFT_RESCUED |
| 166848769 | C | 8 | 4 | 11 | 7.656931 | 7.6568169 | -9.856 | 50.83 | 11.99 | 190.9 | -2.5 | 9.71 | 15.73 | LINEAR_DRIFT_RESCUED |
| 332644900 | C | 7 | 4 | 7 | 1.32907 | 1.3290876 | 1.522 | 64.26 | 38.09 | 66.3 | 30.4 | 2.63 | 4.88 | DYNAMICAL_TTV_CANDIDATE |
| 206666344 | S | 7 | 5 | 7 | 10.510106 | 10.5101753 | 5.99 | 27.06 | 22.31 | 1.7 | 1.0 | 8.19 | 8.11 | ON_SCHEDULE_STABLE |
| 321011164 | S | 10 | 3 | 3 | 22.578361 | 22.5781906 | -14.723 | 44.27 | 3.1 | 81.5 | 0.0 | 9.99 | 11.36 | LINEAR_DRIFT_RESCUED |
| 170788495 | S | 9 | 5 | 6 | 23.907109 | 23.9068806 | -19.731 | 22.92 | 11.4 | 7.7 | -0.2 | 10.16 | 10.02 | ON_SCHEDULE_STABLE |
| 115241103 | C | 10 | 5 | 4 | 1.864247 | 1.8642585 | 0.992 | 26.43 | 30.83 | 2.8 | 0.0 | 9.94 | 9.25 | ON_SCHEDULE_STABLE |
| 7780510 | S | 8 | 5 | 5 | 15.321453 | 15.3211834 | -23.291 | 73.95 | 92.29 | 135.7 | 0.0 | 2.21 | 10.91 | LINEAR_DRIFT_RESCUED |
| 253923533 | likely nearby FP | 6 | 4 | 9 | 1.156013 | 1.1560261 | 1.131 | 40.11 | 18.39 | 511.1 | 33.0 | 6.28 | 6.54 | DYNAMICAL_TTV_CANDIDATE |
| 115692561 | likely nearby FP | 9 | 4 | 13 | 0.850322 | 0.8503703 | 4.173 | 141.76 | 125.45 | 1036.2 | 446.0 | 7.92 | 4.93 | DYNAMICAL_TTV_CANDIDATE |

## 3. VizieR Gaia DR3 + WCS Sub-Pixel Blend Audit

| tic | upstream_tier | P_d | gaia_sid | ruwe | pm_masyr | n_sec | off_target_as | err_as | plausible_nbs | closest_nb_sep_as | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 115692561 | likely nearby FP | 0.85037 | 5004730780766383104 | 1.0 | 133.7 | 2 | 13.83 | 1.56 | 4 | 3.41 | OFFSET_MARGINAL |
| 253923533 | likely nearby FP | 1.156026 | 6600240625047103872 | 1.0 | 17.2 | 1 | 5.08 | 2.2 | 3 | 3.25 | ON_TARGET_UNRESOLVED_NEIGHBOR |
| 115241103 | C | 1.864259 | 4997818391680626048 | 1.0 | 61.9 | 3 | 15.21 | 1.27 | 2 | 4.52 | OFFSET_MARGINAL |
| 166848769 | C | 7.656817 | 4858465154105701120 | 1.0 | 10.0 | 4 | 0.16 | 1.16 | 1 | 40.62 | ON_TARGET_CLEAN |
| 70557792 | S | 7.048826 | 5069952436136567296 | 1.0 | 57.2 | 2 | 22.58 | 1.95 | 0 | -1.0 | OFFSET_MARGINAL |
| 224279805 | reported | 22.984167 | 6534130392483208320 | 1.0 | 31.6 | 1 | 4.78 | 2.2 | 1 | 45.53 | ON_TARGET_CLEAN |
| 332644900 | C | 1.329088 | 5096477810560211072 | 1.0 | 57.5 | 3 | 1.42 | 1.55 | 3 | 2.89 | ON_TARGET_UNRESOLVED_NEIGHBOR |
| 31281821 | MDWARF_ROT | 1.066713 | 2906668221931697024 | 1.0 | 22.3 | 2 | 50.9 | 1.56 | 2 | 40.62 | NEIGHBOR_BEB_LOCALIZED |

## 4. Two-Pass De-Spotted Cool-Dwarf Survey

| tic | lane | p_rot_d | raw_ppm | clean_ppm | P_d | MES | passed | vetoes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1274935 | TWO_PASS_COOL_DWARF | 2.97803 | 6327.8 | 2901.7 | 21.163593 | 6.78 | False | sub_threshold |
| 168675553 | TWO_PASS_ROTATOR_CHECK | 0.47139 | 45596.8 | 43965.2 | 0.235799 | 6.74 | False | sub_threshold |
| 172081181 | TWO_PASS_COOL_DWARF | 3.08786 | 21390.2 | 21849.0 | 3.55371 | 6.55 | False | sub_threshold |
| 209462764 | TWO_PASS_ROTATOR_CHECK | 0.27009 | 51019.1 | 51642.8 | 2.654723 | 6.23 | False | sub_threshold |
| 31281821 | TWO_PASS_ROTATOR_CHECK | 0.83561 | 26014.1 | 24242.9 | 0.936417 | 6.15 | False | sub_threshold |
| 65345801 | TWO_PASS_COOL_DWARF | 0.89946 | 13976.5 | 13580.6 | 2.473478 | 6.01 | False | sub_threshold |
| 178241984 | TWO_PASS_COOL_DWARF | 0.96532 | 12695.4 | 12360.6 | 0.6627 | 5.96 | False | sub_threshold |
| 121026331 | TWO_PASS_COOL_DWARF | 5.47396 | 5147.8 | 5219.0 | 2.768738 | 5.9 | False | sub_threshold |
| 168750720 | TWO_PASS_COOL_DWARF | 0.77817 | 11976.8 | 10816.1 | 0.915935 | 5.9 | False | sub_threshold |
| 94988963 | TWO_PASS_COOL_DWARF | 3.13176 | 6418.8 | 6471.9 | 1.217594 | 5.87 | False | sub_threshold |
| 92303437 | TWO_PASS_COOL_DWARF | 0.56719 | 18097.9 | 17088.5 | 1.70527 | 5.84 | False | sub_threshold |
| 117626744 | TWO_PASS_COOL_DWARF | 0.25857 | 14537.8 | 14910.9 | 3.406102 | 5.81 | False | sub_threshold |
| 13220060 | TWO_PASS_COOL_DWARF | 3.23037 | 5585.5 | 5722.2 | 2.782322 | 5.79 | False | sub_threshold |
| 97507947 | TWO_PASS_COOL_DWARF | 0.83997 | 20611.4 | 21134.5 | 0.690268 | 5.78 | False | sub_threshold |
| 268639026 | TWO_PASS_COOL_DWARF | 0.8966 | 23750.8 | 23609.5 | 2.183344 | 5.71 | False | sub_threshold |
| 178349868 | TWO_PASS_COOL_DWARF | 1.43579 | 9176.7 | 6996.5 | 1.011643 | 5.56 | False | sub_threshold |
| 72667734 | TWO_PASS_COOL_DWARF | 0.71133 | 13932.9 | 12024.0 | 1.225051 | 5.52 | False | sub_threshold |
| 92250728 | TWO_PASS_COOL_DWARF | 4.00011 | 4638.5 | 4621.2 | 1.535724 | 5.5 | False | sub_threshold |
| 121710826 | TWO_PASS_COOL_DWARF | 2.59396 | 9197.8 | 7698.1 | 1.831433 | 5.44 | False | sub_threshold |
| 172525833 | TWO_PASS_COOL_DWARF | 1.38176 | 17135.1 | 16770.8 | 0.933069 | 5.41 | False | sub_threshold |
| 89255318 | TWO_PASS_COOL_DWARF | 2.26331 | 9327.8 | 8941.3 | 1.046739 | 5.38 | False | sub_threshold |
| 77995349 | TWO_PASS_COOL_DWARF | 1.22251 | 6674.9 | 5975.8 | 2.062462 | 5.38 | False | sub_threshold |
