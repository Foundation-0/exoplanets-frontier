# exoplanets-frontier

This repository tests transit predictions and searches for uncatalogued exoplanets in NASA TESS starlight through Sector 107. It extends [`eumemic/exoplanets`](https://github.com/eumemic/exoplanets) (`v3.0`, which stops at Sector 105 and periods of 0.5 to 30 days) into newly released sectors and shorter orbits. Researchers and collaborators can run the code to verify existing candidates and inspect new multi-sector transit signals.

## Campaigns

1. **Out-of-sample Sector 106 and 107 verification**: Tests every `eumemic/exoplanets` (`v3.0`) candidate observed at 2-minute SPOC cadence in TESS Sectors 106 and 107 against its pre-committed ephemeris (`T0_BJD`, `P_d`, `depth_ppm`, `T14_h`).
2. **Multi-year stacking of unpromoted Sector 106 and 107 signals**: Filters NASA Ames SPOC single-sector Threshold Crossing Events (TCEs) from Sectors 106 and 107 to remove all known planets (`PSCompPars`), TESS Objects of Interest (`TOI`), community candidates (`CTOI`), and `eumemic/exoplanets` entries, then stacks each dwarf star across all available SPOC sectors from 2018 to 2026 (Sectors 1 to 107).
3. **Ultra-short-period (0.155 to 0.50 day) sibling search**: Masks known candidate transits on multi-sector host stars and searches unbinned 2-minute SPOC light curves for inner rocky companions below the 0.5-day floor used in `v3.0`.

## Eight-Gate Automated Vetting

Every candidate must clear eight false-positive checks before entering `results/vetted_new_candidates.csv`:

1. Multi-sector coverage (`>= 3` transits across `>= 2` TESS sectors) and Pont red-noise corrected `MES_corr >= 9.5`, `SDE >= 7.5`.
2. Inverted light-curve symmetry (`MES_inv <= 0.70 * MES_corr` or `< 7.5`).
3. Odd-even transit depth consistency (`<= 2.8 sigma`).
4. Secondary eclipse veto across orbital phases `0.18` to `0.82` (`<= 4.2 sigma`).
5. Limb-darkened transit geometry (`0.45 <= Rp <= 8.0` Earth radii, non-grazing shape, Keplerian duration match).
6. Single-event dominance veto (`<= 55%` of `MES^2` from any single transit).
7. Flux-weighted centroid stability (`MOM_CENTR - POS_CORR`) and 2D `_tp.fits` Target Pixel File difference-image localization (`<= 15 arcsec`).
8. Gaia DR3 astrometry (`RUWE <= 1.4`, radial-velocity scatter `<= 5 km/s`) and neighbor false-positive bounds (`FPP < 0.50`, `NFPP < 0.10`; `Tier S-Val` requires `FPP < 0.015`, `NFPP < 0.001`).

## Usage

```bash
pip install numpy scipy pandas matplotlib astropy astroquery requests
python frontier.py --self-test
python frontier.py --out-dir results --max-tce 48 --max-usp 36
```
