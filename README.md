# Exoplanets Frontier (`Foundation-0/exoplanets-frontier`)

Reproducible transit verification, dynamical timing analysis, sub-pixel Gaia blend localization, and harmonic-prewhitened planet discovery on live NASA STScI MAST TESS Sectors 1-107 photometry and ESA Gaia DR3 astrometry, building on [eumemic/exoplanets](https://github.com/eumemic/exoplanets) within the `Foundation-0` astronomy program (`astr0`, `gaia-dr4`, `lost-eclipses`, `certified-orbits`, `cosmic-denoising`).

## Repository Navigation and Team Capture Surfaces

Every method, verified claim, falsified hypothesis, and reusable module is captured on shared repository surfaces for immediate continuation:

- **[`METHODS.md`](METHODS.md)**: Mathematical and algorithmic reference covering two-stage alias-resolving BLS, the 8-gate astrophysical vetter, the 4 technical bottlenecks encountered in Phase 1, and the convergent techniques built to resolve them.
- **[`CLAIMS.md`](CLAIMS.md)**: Numbered epistemic ledger (`C-01` through `C-06+`) tagged `VERIFIED`, `FALSIFIED`, `CONJECTURE`, or `KNOWN` with explicit falsifiers and SHA-256 certificate pins.
- **[`HANDOFF.md`](HANDOFF.md)**: Cross-repository handoff guide connecting our data products and code directly to `Foundation-0/gaia-dr4` (ahead of the 25 Nov 2026 ESA archive maintenance window), `Foundation-0/lost-eclipses`, `Foundation-0/certified-orbits`, and `eumemic/exoplanets`.
- **[`CHANGELOG.md`](CHANGELOG.md)**: Chronological record of campaigns, bug fixes, and empirical milestones.
- **Empirical Ledgers & Certificates (`results/`)**:
  - [`results/CONVERGENCE_REPORT.md`](results/CONVERGENCE_REPORT.md) & [`results/convergence_certificates.json`](results/convergence_certificates.json): Phase 2 benchmark recovery, scrambled null test, 8-year S1-S107 wide-window ephemeris/TTV audit, multi-sector WCS + Gaia DR3 blend localization, and de-spotted M-dwarf residual planet search.
  - [`results/DISCOVERIES.md`](results/DISCOVERIES.md) & [`results/blind_verification_s106_s107.csv`](results/blind_verification_s106_s107.csv): Out-of-sample Sector 106 and 107 blind test of 39 `v3.0` candidate-sector pairs (`22 CONFIRMED_ON_SCHEDULE`, `0 NON_DETECTION_4SIGMA`; reported in `eumemic/exoplanets#20`).
  - [`results/figures/`](results/figures/): Diagnostic phase-fold cards, multi-year O-C timing curves, harmonic de-spotting comparisons, and WCS-registered sub-pixel difference images.

## Module Map

| Module | Role |
| :--- | :--- |
| `bls_search.py` | Orbit-aware detrending, Pont red-noise $\beta$, and two-stage multi-year alias-resolving AstroPy `BoxLeastSquares` |
| `vetter.py` | Eight-gate false-positive vetter, Sector 106-107 blind verifier, and 4-panel phase-fold diagnostic plotter |
| `mast_data.py` | Live STScI MAST SPOC FITS/TPF/DVT indexer, multi-sector light-curve stacker, and catalog cross-matcher |
| `deep_search.py` | `Tesscut` 9x9 FFI aperture extractor, high-pass centroid detrender, and held-out sector blind evaluator |
| `mdwarf_hunt.py` | Per-sector `SAP_FLUX` normalizer (`norm_sap`), out-of-transit spot-wave auditor (`spot_amp`), and cool-dwarf survey |
| `despot.py` | Asymmetric Huber 4-harmonic starspot pre-whitener (`harmonic_despot`), scrambled null test, and de-spot plots |
| `ttv_track.py` | Wide-window (`+-240 min`) S1-S107 O-C timing extractor, 4-model $\Delta\text{BIC}$ ephemeris/TTV solver, and sibling scan |
| `blend_gaia.py` | Proper-motion-propagated Gaia DR3 scene registration, multi-sector FITS WCS sub-pixel difference centroiding, and O-C/WCS plots |
| `converge_io.py` / `converge.py` | End-to-end runner executing Benchmark + Null (`astr0` Rules 2-3), 8-Year TTV/Drift Rescue, Gaia WCS Blend Audit, and De-Spotted M-Dwarf Survey with SHA-256 certificates |

## Running Locally or via GitHub Actions

```bash
uv venv .venv
uv pip install --python .venv/bin/python numpy scipy pandas matplotlib astropy astroquery
.venv/bin/python -m unittest ledger.py
.venv/bin/python converge.py
```

Live archival campaigns execute reproducibly via `.github/workflows/campaign.yml`, `.github/workflows/mdwarf_hunt.yml`, and `.github/workflows/converge.yml`.
