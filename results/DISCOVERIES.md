# TESS Sectors 1-107 Verification and Discovery Ledger

Empirical measurements from NASA STScI MAST 2-minute SPOC light curves and Target Pixel Files (`2026-11-11T19:39:13Z`, `182.1 s`).

- **Highest MAST SPOC sector indexed**: `S107` (new sectors tested: `[106, 107]`)
- **Out-of-sample (`v3.0`, S106/S107) pairs tested**: `3` -> `{'LOW_SNR_INCONCLUSIVE': 2, 'MARGINAL_RECOVERY': 1}`
- **Uncatalogued S106-S107 dwarf signals found / stacked across S1-S107**: `15` found, `15` stacked
- **Hosts searched for 0.155-0.50 d USP siblings**: `3` | **Passed all 8 gates + TPF**: `0`

## 1. Out-of-Sample Blind Verification (Sectors 106-107)

| TIC | Tier | Sec | Period (d) | Pred (ppm) | Obs (ppm) | SNR | Offset (m) | Cent | Verdict |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| `32816844` | `C` | S107 | 14.45236 | 80 | 92+/-24 | 3.91 | -14.0 | 5.16 | `MARGINAL_RECOVERY` |

## 2. New Uncatalogued Planet Candidates (All 8 Gates + TPF Passed)

Zero triggers in this batch passed all 8 gates simultaneously; see audit table below.

### Vetted Trigger Audit (Top Signals & Veto Flags)

| TIC | Lane | Period (d) | Depth (ppm) | Rp (Re) | MES | Passed | Vetoes |
| :--- | :--- | ---: | ---: | ---: | ---: | :--- | :--- |
| `627018` | `S106_S107_STACK` | 0.84487 | 23335 | 7.13 | 33.9 | `False` | `inv_lc(29.8);sec_ecl(58.6);sap_mismatch;rv(27.8)` |
| `671393` | `S106_S107_STACK` | 0.44206 | 27504 | 13.38 | 26.4 | `False` | `rp(13.38Re);odd_even(3.66s);inv_lc(24.4);sec_ecl(95.5)` |
| `1219543` | `S106_S107_STACK` | 0.57789 | 26713 | 6.06 | 26.0 | `False` | `sec_ecl(14.9)` |
| `1273249` | `S106_S107_STACK` | 1.80401 | 4301 | 6.55 | 21.3 | `False` | `inv_lc(30.6);sec_ecl(38.5);sap_mismatch` |
| `672614` | `S106_S107_STACK` | 0.78185 | 21106 | 5.15 | 17.3 | `False` | `sec_ecl(24.5)` |
| `1308855` | `S106_S107_STACK` | 0.47882 | 125 | nan | 16.0 | `False` | `rp(nanRe);dur_ratio(nan);inv_lc(12.5);sec_ecl(26.6);sap_mismatch;cent(15.7as)` |
| `11614485` | `S106_S107_STACK` | 6.00147 | 14795 | 4.87 | 14.1 | `False` | `odd_even(3.69s);inv_lc(11.0);sec_ecl(69.4)` |
| `9166136` | `S106_S107_STACK` | 6.08063 | 1216 | 3.38 | 13.0 | `False` | `inv_lc(10.3);cent(85.1as)` |
| `684071` | `S106_S107_STACK` | 14.38900 | 14986 | 18.58 | 12.2 | `False` | `rp(18.58Re);inv_lc(17.1);sec_ecl(5.4)` |
| `1439071` | `S106_S107_STACK` | 7.51590 | 5379 | 3.86 | 10.7 | `False` | `odd_even(7.40s);inv_lc(9.4);sec_ecl(46.8);rv(39.9)` |
| `709290` | `S106_S107_STACK` | 1.26209 | 22194 | 7.95 | 10.4 | `False` | `odd_even(3.60s);inv_lc(12.5);sec_ecl(32.6)` |
