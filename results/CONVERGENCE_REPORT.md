# Phase 3 Synthesis: 665-Candidate Gaia DR4 Catalog & Cycle 1-8 Ephemerides (1109.0 s)

- Catalog: 3543 indexed ({'CROWDED_INDEXED': 1927, 'BINARY_OR_EVOLVED_WARNING': 844, 'LOW_CONTAM_INDEXED': 724, 'SILVER_CROWDED_SG1': 26, 'GOLD_LOW_BLEND': 19, 'PLATINUM_SOLO_VERIFIED': 3}), 22 Platinum/Gold (6 temperate HZ 180-325 K).
- TOI-540b Benchmark: P=1.239126 d (MES=14.33), Null=0/10 | Cool-Dwarfs (33): 0 passed.

## 1. Top Platinum Solo-Verified & Temperate Candidates

| tic | gaia_dr3_id | upstream_tier | P_d | Rp_Re | Teq_K | is_hz_temperate | TSM | Tmag | ruwe | plausible_beb_65as | dr4_class |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 181015118 | 5382868669400908928 | nan | 25.26197 | 2.28 | 258.0 | True | 97.4 | 8.99 | 1.0 | 0 | PLATINUM_SOLO_VERIFIED |
| 181015118 | 5382868669400908928 | nan | 18.317993 | 1.78 | 288.0 | True | 64.6 | 8.99 | 1.0 | 1 | GOLD_LOW_BLEND |
| 246897668 | 2980822034888930176 | nan | 24.752432 | 2.1 | 282.0 | True | 58.1 | 9.24 | 1.0 | 0 | PLATINUM_SOLO_VERIFIED |
| 37904598 | 3199132439405197696 | nan | 17.049176 | 2.77 | 303.0 | True | 44.3 | 11.21 | 1.0 | 1 | GOLD_LOW_BLEND |
| 146636926 | 4947259307983420928 | nan | 26.056445 | 2.55 | 251.0 | True | 42.7 | 11.47 | 1.0 | 1 | GOLD_LOW_BLEND |
| 146636926 | 4947259307983420928 | nan | 23.920264 | 2.11 | 258.0 | True | 29.8 | 11.47 | 1.0 | 1 | GOLD_LOW_BLEND |
| 79278415 | 6320381307658461184 | S | 3.348204 | 6.55 | 1328.0 | False | 88.0 | 12.03 | 1.0 | 1 | GOLD_LOW_BLEND |
| 241537628 | 2481055723313271168 | S | 14.787554 | 7.76 | 825.0 | False | 84.1 | 12.13 | 1.0 | 0 | GOLD_LOW_BLEND |
| 211405734 | 6783075699211325696 | nan | 7.297485 | 7.53 | 444.0 | False | 157.0 | 13.48 | 1.0 | 0 | GOLD_LOW_BLEND |
| 84191887 | 1953355906525258752 | nan | 2.12798 | 5.57 | 598.0 | False | 147.2 | 13.5 | 1.0 | 1 | GOLD_LOW_BLEND |
| 283669507 | 1780614933214780800 | nan | 2.180834 | 4.19 | 1772.0 | False | 152.6 | 9.0 | 1.0 | 0 | PLATINUM_SOLO_VERIFIED |
| 105655329 | 1283573314535218560 | nan | 1.800672 | 4.58 | 1356.0 | False | 154.2 | 10.55 | 1.0 | 0 | GOLD_LOW_BLEND |
| 355072530 | 272635972674562176 | C | 13.698274 | 9.09 | 780.0 | False | 154.8 | 12.66 | 1.0 | 1 | GOLD_LOW_BLEND |
| 399948046 | 12893530478007552 | likely nearby FP | 15.913365 | 8.27 | 607.0 | False | 122.6 | 12.82 | 1.0 | 1 | GOLD_LOW_BLEND |
| 244252435 | 3212738719703024384 | nan | 0.566432 | 3.43 | 1235.0 | False | 133.0 | 11.05 | 1.0 | 1 | GOLD_LOW_BLEND |
| 156380270 | 3884882276685154944 | C | 1.521754 | 6.4 | 899.0 | False | 117.3 | 13.45 | 1.0 | 0 | GOLD_LOW_BLEND |
| 79390722 | 6321177835113922816 | likely nearby FP | 4.423039 | 5.49 | 1117.0 | False | 117.2 | 10.98 | 1.0 | 1 | GOLD_LOW_BLEND |
| 16983459 | 1370728607676424192 | likely nearby FP | 0.811675 | 6.15 | 1915.0 | False | 116.9 | 12.53 | 1.0 | 1 | GOLD_LOW_BLEND |
| 437709086 | 2822231810591184896 | S | 2.925309 | 6.76 | 1504.0 | False | 76.2 | 12.63 | 1.0 | 1 | GOLD_LOW_BLEND |
| 1400125 | 5683477411879049216 | nan | 0.948647 | 2.37 | 1341.0 | False | 113.3 | 9.47 | 1.0 | 1 | GOLD_LOW_BLEND |
| 372946271 | 1759254067988060288 | C | 2.956309 | 8.47 | 1478.0 | False | 111.1 | 12.84 | 1.0 | 1 | GOLD_LOW_BLEND |
| 343942488 | 2823623208196929920 | S | 0.653612 | 8.84 | 2173.0 | False | 205.4 | 13.3 | 1.0 | 0 | GOLD_LOW_BLEND |

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
