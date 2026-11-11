# TESS Sectors 1-106 Verification and Discovery Ledger

Empirical measurements from NASA STScI MAST 2-minute SPOC light curves and Target Pixel Files (`2026-11-11T23:51:07Z`, `2764.0 s`).

- **Highest MAST SPOC sector indexed**: `S106` (new sectors tested: `[106]`)
- **Out-of-sample (`v3.0`, S106/S107) pairs tested**: `39` -> `{'CONFIRMED_ON_SCHEDULE': 22, 'LOW_SNR_INCONCLUSIVE': 9, 'MARGINAL_RECOVERY': 5, 'DATA_GAP': 2, 'TIMING_OFFSET_TTV': 1}`
- **Uncatalogued S106-S107 dwarf signals found / stacked across S1-S106**: `128` found, `128` stacked
- **Hosts searched for 0.155-0.50 d USP siblings**: `28` | **Passed all 8 gates + TPF**: `0`

## 1. Out-of-Sample Blind Verification (Sectors 106-107)

| TIC | Tier | Sec | Period (d) | Pred (ppm) | Obs (ppm) | SNR | Offset (m) | Cent | Verdict |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| `224279805` | `reported(FFI)` | S106 | 22.98417 | 7247 | 6760+/-311 | 21.72 | +0.0 | 0.56 | `CONFIRMED_ON_SCHEDULE` |
| `321011164` | `S(FFI)` | S106 | 22.57836 | 2330 | 2336+/-238 | 9.83 | -14.0 | 2.69 | `CONFIRMED_ON_SCHEDULE` |
| `70557792` | `S(FFI)` | S107 | 7.04883 | 467 | 666+/-80 | 8.28 | +10.0 | 1.21 | `CONFIRMED_ON_SCHEDULE` |
| `166848769` | `C(FFI)` | S107 | 7.65693 | 6813 | 4178+/-605 | 6.90 | -48.0 | 0.99 | `CONFIRMED_ON_SCHEDULE` |
| `321011164` | `S(FFI)` | S107 | 22.57836 | 2330 | 5192+/-796 | 6.52 | -60.0 | 14.34 | `MARGINAL_RECOVERY` |
| `469942354` | `S(FFI)` | S106 | 3.59347 | 619 | 741+/-118 | 6.28 | -4.0 | 0.14 | `CONFIRMED_ON_SCHEDULE` |
| `115241103` | `C(FFI)` | S106 | 1.86425 | 364 | 389+/-64 | 6.10 | -6.0 | 1.63 | `CONFIRMED_ON_SCHEDULE` |
| `166848769` | `C(FFI)` | S106 | 7.65693 | 6813 | 3330+/-551 | 6.05 | -58.0 | 0.56 | `CONFIRMED_ON_SCHEDULE` |
| `170788495` | `S(FFI)` | S107 | 23.90711 | 1326 | 2168+/-364 | 5.96 | -28.0 | 1.39 | `CONFIRMED_ON_SCHEDULE` |
| `146654196` | `S(FFI)` | S106 | 4.89736 | 3173 | 2831+/-478 | 5.93 | +2.0 | 1.55 | `CONFIRMED_ON_SCHEDULE` |
| `170788495` | `S(FFI)` | S106 | 23.90711 | 1326 | 1664+/-290 | 5.75 | -8.0 | 2.23 | `CONFIRMED_ON_SCHEDULE` |
| `398794585` | `S(FFI)` | S107 | 9.54315 | 1337 | 1116+/-205 | 5.43 | +18.0 | 1.08 | `CONFIRMED_ON_SCHEDULE` |
| `170638116` | `S(FFI)` | S106 | 28.24869 | 1224 | 1277+/-290 | 4.41 | +16.0 | 1.99 | `CONFIRMED_ON_SCHEDULE` |
| `115692561` | `likely nearby FP(FFI)` | S106 | 0.85032 | 213 | 149+/-35 | 4.26 | -2.0 | 1.98 | `CONFIRMED_ON_SCHEDULE` |
| `115273010` | `S(FFI)` | S106 | 4.31417 | 3138 | 1856+/-456 | 4.08 | +4.0 | 2.42 | `CONFIRMED_ON_SCHEDULE` |
| `270312820` | `S(FFI)` | S106 | 6.54573 | 651 | 709+/-179 | 3.95 | -18.0 | 1.50 | `CONFIRMED_ON_SCHEDULE` |
| `382532095` | `S(FFI)` | S106 | 9.00052 | 486 | 386+/-103 | 3.75 | -20.0 | 1.48 | `CONFIRMED_ON_SCHEDULE` |
| `115692561` | `likely nearby FP(FFI)` | S107 | 0.85032 | 213 | 163+/-47 | 3.46 | +0.0 | 1.30 | `CONFIRMED_ON_SCHEDULE` |
| `206666344` | `S(FFI)` | S106 | 10.51011 | 1118 | 764+/-223 | 3.43 | +60.0 | 1.74 | `CONFIRMED_ON_SCHEDULE` |
| `115241103` | `C(FFI)` | S107 | 1.86425 | 364 | 249+/-75 | 3.34 | +4.0 | 0.56 | `CONFIRMED_ON_SCHEDULE` |
| `229115129` | `S(FFI)` | S106 | 1.73794 | 923 | 660+/-203 | 3.26 | -18.0 | 1.04 | `CONFIRMED_ON_SCHEDULE` |
| `253923533` | `likely nearby FP(FFI)` | S106 | 1.15601 | 1039 | 788+/-246 | 3.20 | +2.0 | 0.22 | `CONFIRMED_ON_SCHEDULE` |
| `218962245` | `S(FFI)` | S106 | 12.10051 | 1221 | 870+/-286 | 3.05 | -14.0 | 1.07 | `CONFIRMED_ON_SCHEDULE` |
| `11650366` | `reported(FFI)` | S107 | 7.09264 | 1584 | 913+/-302 | 3.02 | +6.0 | 4.15 | `MARGINAL_RECOVERY` |
| `70557792` | `S(FFI)` | S106 | 7.04883 | 467 | 265+/-94 | 2.81 | -20.0 | 1.18 | `TIMING_OFFSET_TTV` |
| `12377308` | `S(FFI)` | S107 | 9.13668 | 921 | 451+/-169 | 2.67 | -2.0 | 0.37 | `MARGINAL_RECOVERY` |
| `7780510` | `S(FFI)` | S107 | 15.32145 | 1108 | 470+/-198 | 2.37 | +52.0 | 0.68 | `MARGINAL_RECOVERY` |
| `11650366` | `reported(FFI)` | S106 | 7.09264 | 1584 | 877+/-397 | 2.21 | -2.0 | 1.77 | `MARGINAL_RECOVERY` |

## 2. New Uncatalogued Planet Candidates (All 8 Gates + TPF Passed)

Zero triggers in this batch passed all 8 gates simultaneously; see audit table below.

### Vetted Trigger Audit (Top Signals & Veto Flags)

| TIC | Lane | Period (d) | Depth (ppm) | Rp (Re) | MES | Passed | Vetoes |
| :--- | :--- | ---: | ---: | ---: | ---: | :--- | :--- |
| `77842052` | `MDWARF_T>13.5` | 0.38191 | 12249 | 3.26 | 27.7 | `False` | `dur_ratio(3.18);odd_even(3.57s);inv_lc(20.0);sec_ecl(61.3);sap_mismatch` |
| `31353324` | `MDWARF_T>13.5` | 0.54063 | 15684 | 4.43 | 25.7 | `False` | `dur_ratio(2.69);odd_even(3.30s);inv_lc(23.0);sec_ecl(63.3);rv(5.4)` |
| `166833387` | `MDWARF_T>13.5` | 0.55063 | 21088 | 4.70 | 25.7 | `False` | `dur_ratio(2.70);sec_ecl(22.9)` |
| `209462764` | `MDWARF_T>13.5` | 0.26992 | 14234 | 3.01 | 23.8 | `False` | `dur_ratio(3.97);odd_even(3.00s);sec_ecl(25.0)` |
| `168675553` | `MDWARF_T>13.5` | 0.47155 | 28420 | 4.19 | 22.9 | `False` | `dur_ratio(3.47);sec_ecl(31.2);sap_mismatch` |
| `304096841` | `MDWARF_T>13.5` | 0.46931 | 9240 | 3.81 | 22.8 | `False` | `odd_even(2.85s);inv_lc(28.7);sec_ecl(34.8);sap_mismatch` |
| `321087099` | `MDWARF_T>13.5` | 1.09091 | 8391 | 3.17 | 17.4 | `False` | `inv_lc(13.6);sec_ecl(23.2);sap_mismatch` |
| `7600145` | `MDWARF_T>13.5` | 0.32503 | 4284 | 1.70 | 14.6 | `False` | `dur_ratio(2.71);inv_lc(12.6);sec_ecl(25.2);sap_mismatch;ruwe(1.79)` |
| `11614485` | `S106_S107_STACK` | 6.00147 | 14795 | 4.87 | 14.1 | `False` | `odd_even(3.69s);inv_lc(11.0);sec_ecl(69.4)` |
| `142938655` | `MDWARF_T>13.5` | 1.00406 | 12304 | 6.60 | 14.1 | `False` | `alias;inv_lc(12.2);sec_ecl(22.4);sap_mismatch;rv(6.2)` |
| `129739941` | `S106_S107_STACK` | 4.35175 | 807 | 2.96 | 13.2 | `False` | `odd_even(8.37s);inv_lc(14.7);ruwe(1.41)` |
| `9166136` | `S106_S107_STACK` | 6.08063 | 1216 | 3.38 | 13.0 | `False` | `inv_lc(10.3);cent(85.1as)` |
| `77959225` | `MDWARF_T>13.5` | 0.27326 | 4109 | 2.96 | 11.8 | `False` | `dur_ratio(2.78);odd_even(2.97s);inv_lc(17.2);sec_ecl(25.1);sap_mismatch;ruwe(1.45)` |
| `31281821` | `MDWARF_T>13.5` | 1.06671 | 12948 | 6.16 | 11.6 | `False` | `sap_mismatch` |
| `31118029` | `MDWARF_T>13.5` | 0.31726 | 3135 | 2.51 | 11.5 | `False` | `odd_even(5.21s);inv_lc(11.6);sec_ecl(12.4);sap_mismatch;ruwe(17.42)` |
| `44794000` | `MDWARF_T>13.5` | 0.53412 | 22092 | 9.50 | 10.8 | `False` | `rp(9.50Re);dur_ratio(3.19);inv_lc(11.2);sec_ecl(22.1);sap_mismatch;ruwe(7.72)` |
| `1439071` | `S106_S107_STACK` | 7.51590 | 5379 | 3.86 | 10.7 | `False` | `odd_even(7.40s);inv_lc(9.4);sec_ecl(46.8);rv(39.9)` |
| `139637517` | `MDWARF_T>13.5` | 0.51885 | 10056 | 3.52 | 10.7 | `False` | `inv_lc(9.7);sec_ecl(11.9)` |
| `44892011` | `MDWARF_T>13.5` | 0.83075 | 4007 | 1.43 | 10.1 | `False` | `dur_ratio(4.08);odd_even(2.93s);inv_lc(8.3);sec_ecl(5.5)` |
