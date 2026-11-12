# Phase 2 Convergent Discovery and Technical Resolution Report

Live STScI MAST Sectors 1-107 and ESA Gaia DR3 run (`947.3 s`).

## 1. Benchmark & Null Control (`astr0` Rules 2-3)

- **TOI-540 (`TIC 200322593`, $P_b=1.239149\text{ d}$)**: Naive P=`0.726307 d` -> De-Spotted P=`1.239126 d` (`MES=14.33`, recovered=`True`, RMS `2828.7->2484.1 ppm`). Null false alarms: `0/10` (max MES `5.77`).

## 2. 8-Year S1-S107 Wide-Window (`+-240 min`) Ephemeris & TTV Audit

| TIC | Tier | P_orig | P_ref | dP (s) | Ep | RMS0 (m) | RMS_lin (m) | dBIC_lin | dBIC_ttv | SNR_0 -> SNR_ref | Class |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- | :--- |
| `70557792` | `S` | 7.048829 | 7.0488122 | -1.45 | 8 | 32.6 | 23.5 | 21.0 | 4.4 | 9.46 -> **10.22** | `LINEAR_DRIFT_RESCUED` |
| `166848769` | `C` | 7.656931 | 7.6574415 | +44.11 | 8 | 59.8 | 10.5 | 212.5 | -3.4 | 6.29 -> **13.49** | `LINEAR_DRIFT_RESCUED` |
| `332644900` | `C` | 1.329070 | 1.3290094 | -5.24 | 8 | 79.7 | 30.9 | 136.0 | 18.7 | 2.67 -> **5.39** | `DYNAMICAL_TTV_CANDIDATE` |
| `206666344` | `S` | 10.510106 | 10.5101683 | +5.39 | 2 | 29.4 | 0.0 | 2.7 | 0.0 | 5.08 -> **5.86** | `ON_SCHEDULE_STABLE` |
| `321011164` | `S` | 22.578361 | 22.5548367 | -2032.50 | 2 | 90.6 | 0.0 | 223.0 | 0.0 | 7.07 -> **10.24** | `LINEAR_DRIFT_RESCUED` |
| `170788495` | `S` | 23.907109 | 23.9261646 | +1646.40 | 2 | 18.3 | 0.0 | 2.8 | 0.0 | 7.92 -> **8.32** | `ON_SCHEDULE_STABLE` |
| `115241103` | `C` | 1.864247 | 1.8642470 | +0.00 | 1 | 0.0 | 0.0 | 0.0 | 0.0 | 4.15 -> **4.15** | `INSUFFICIENT_EPOCHS` |
| `7780510` | `S` | 15.321453 | 15.3920190 | +6096.90 | 2 | 105.0 | 0.0 | 179.6 | 0.0 | 0.84 -> **9.13** | `LINEAR_DRIFT_RESCUED` |
| `398794585` | `S` | 9.543154 | 9.5431540 | +0.00 | 1 | 0.0 | 0.0 | 0.0 | 0.0 | 46.64 -> **46.64** | `INSUFFICIENT_EPOCHS` |
| `170638116` | `S` | 28.248693 | 28.2486930 | +0.00 | 1 | 0.0 | 0.0 | 0.0 | 0.0 | 2.51 -> **2.51** | `INSUFFICIENT_EPOCHS` |

## 3. Multi-Sector WCS + Gaia DR3 Sub-Pixel Blend Audit

| TIC | Tier | P (d) | Gaia DR3 ID | RUWE | PM | Sec | Off (") | Plaus NBs | NB Sep (") | Verdict |
| :--- | :--- | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| `115692561` | `likely nearby FP` | 0.850322 | `none` | 1.00 | 0.0 | 2 | 11.51+/-1.56 | 0 | -1.0 | `OFFSET_MARGINAL` |
| `253923533` | `likely nearby FP` | 1.156013 | `none` | 1.00 | 0.0 | 1 | 4.90+/-2.20 | 0 | -1.0 | `ON_TARGET_CLEAN` |
| `115241103` | `C` | 1.864247 | `none` | 1.00 | 0.0 | 2 | 11.09+/-1.56 | 0 | -1.0 | `OFFSET_MARGINAL` |
| `166848769` | `C` | 7.657442 | `none` | 1.00 | 0.0 | 0 | -1.00+/--1.00 | 0 | -1.0 | `NO_DIFF_SNR` |
| `70557792` | `S` | 7.048812 | `none` | 1.00 | 0.0 | 2 | 13.03+/-1.56 | 0 | -1.0 | `OFFSET_MARGINAL` |
| `224279805` | `reported` | 22.984167 | `none` | 1.00 | 0.0 | 1 | 0.55+/-2.20 | 0 | -1.0 | `ON_TARGET_CLEAN` |
| `31281821` | `MDWARF_ROT` | 1.066713 | `none` | 1.00 | 0.0 | 2 | 24.06+/-1.56 | 0 | -1.0 | `OFFSET_MARGINAL` |

## 4. De-Spotted Active M-Dwarf Residual Survey

| TIC | Lane | P_rot (d) | Raw->Clean (ppm) | Resid P (d) | MES | Passed | Vetoes |
| :--- | :--- | ---: | :--- | ---: | ---: | :--- | :--- |
| `1439071` | `DESPOT_FRESH_MDWARF` | 2.5211 | 8662->3630 | 17.77262 | 6.78 | `False` | `sub_threshold` |
| `209462764` | `DESPOT_ROTATOR_RESCUE` | 0.2700 | 51019->51643 | 2.65472 | 6.47 | `False` | `sub_threshold` |
| `44892011` | `DESPOT_ROTATOR_RESCUE` | 1.4058 | 21775->21847 | 2.07420 | 6.27 | `False` | `sub_threshold` |
| `1441539` | `DESPOT_FRESH_MDWARF` | 0.9204 | 8653->6791 | 2.53850 | 6.15 | `False` | `sub_threshold` |
| `663735` | `DESPOT_FRESH_MDWARF` | 9.8039 | 5593->5634 | 1.63632 | 6.07 | `False` | `sub_threshold` |
| `612870` | `DESPOT_FRESH_MDWARF` | 3.3676 | 5937->5747 | 1.10468 | 5.94 | `False` | `sub_threshold` |
| `304096841` | `DESPOT_ROTATOR_RESCUE` | 0.6029 | 11178->9885 | 1.13239 | 5.89 | `False` | `sub_threshold` |
| `674490` | `DESPOT_FRESH_MDWARF` | 4.5931 | 20990->21210 | 1.03640 | 5.88 | `False` | `sub_threshold` |
| `1396655` | `DESPOT_FRESH_MDWARF` | 4.4840 | 14128->14310 | 0.70130 | 5.85 | `False` | `sub_threshold` |
| `1449937` | `DESPOT_FRESH_MDWARF` | 0.5563 | 12240->12395 | 4.68506 | 5.80 | `False` | `sub_threshold` |
| `672614` | `DESPOT_FRESH_MDWARF` | 0.7830 | 22853->20571 | 1.03631 | 5.71 | `False` | `sub_threshold` |
| `1447223` | `DESPOT_FRESH_MDWARF` | 5.6762 | 3521->3557 | 1.67230 | 5.69 | `False` | `sub_threshold` |
| `31281821` | `DESPOT_ROTATOR_RESCUE` | 0.8364 | 26577->24588 | 0.93642 | 5.57 | `False` | `sub_threshold` |
| `7598799` | `DESPOT_FRESH_MDWARF` | 1.9209 | 22148->22151 | 10.94327 | 5.57 | `False` | `sub_threshold` |
| `31353324` | `DESPOT_ROTATOR_RESCUE` | 0.5408 | 18523->16216 | 1.89601 | 5.54 | `False` | `sub_threshold` |
| `77842052` | `DESPOT_ROTATOR_RESCUE` | 0.4507 | 16922->14659 | 1.45810 | 5.44 | `False` | `sub_threshold` |
| `166833387` | `DESPOT_ROTATOR_RESCUE` | 0.5503 | 28933->27973 | 0.74347 | 5.31 | `False` | `sub_threshold` |
| `168675553` | `DESPOT_ROTATOR_RESCUE` | 0.4713 | 45597->43886 | 1.57319 | 5.25 | `False` | `sub_threshold` |
| `321087099` | `DESPOT_ROTATOR_RESCUE` | 1.0865 | 11899->10996 | 1.19805 | 5.18 | `False` | `sub_threshold` |
| `7222629` | `DESPOT_FRESH_MDWARF` | 6.2646 | 2984->3046 | 0.70907 | 5.10 | `False` | `sub_threshold` |
| `1525835` | `DESPOT_FRESH_MDWARF` | 4.2179 | 22328->22862 | 1.67907 | 5.06 | `False` | `sub_threshold` |
| `1534458` | `DESPOT_FRESH_MDWARF` | 4.4465 | 12882->13176 | 1.45505 | 4.87 | `False` | `sub_threshold` |
| `627018` | `DESPOT_FRESH_MDWARF` | 0.8449 | 10662->3786 | 2.88173 | 4.78 | `False` | `sub_threshold` |
| `139637517` | `DESPOT_ROTATOR_RESCUE` | 0.5167 | 33339->33796 | 1.80983 | 4.66 | `False` | `sub_threshold` |
