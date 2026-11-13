# Seeded transit-search control

The unchanged seeded search recovers a fixed 25-day box transit at the correct period with corrected MES 49.76, but its SDE is 7.23. The configured SDE threshold is 7.5, so this signal would fail that gate. A fixed 17-day box transit passes the same gate.

MES measures the fitted transit strength relative to noise. Here SDE measures how far the largest search score rises above the middle of the scores in the searched period window. The two statistics are not interchangeable.

| Injected period | Recovered period | Corrected MES | SDE | SDE gate |
| --- | --- | --- | --- | --- |
| 17 days | 17.0015229138 days | 58.63565 | 8.58803 | Pass |
| 25 days | 25.0010843180 days | 49.76145 | 7.23137 | Fail |

Both cases were fixed before execution: 80 days of continuous coverage at 30-minute cadence, epoch day 5, stellar radius 0.55 solar radii, box depth 1,000 ppm, and independent Gaussian noise with standard deviation 100 ppm. Both use the same noise vector from seed 314159. The injected duration is the unchanged `exp_dur` function's output. No parameters were adjusted after seeing the results.

The script extracts the original `coherent_seed_bls` function from `hz_core.py` and `exp_dur` and `pont_beta` from `bls_search.py`. It executes those functions with NumPy, math and Astropy. This avoids importing the archive pipeline.

The script verifies both source hashes against the reviewed revision before execution. A deliberate source change is rejected before the control runs. Source hashes and software versions are saved in the result JSON. A second execution produces identical JSON bytes.

This demonstrates that the current threshold can reject a strong, correctly recovered synthetic transit. It does not establish which part of the search design causes the rejection. In particular, comparison with a different period window or a calibrated threshold has not been performed.

No detrending, starspot removal, real observation gaps, archive retrieval or full vetter was run. The result is not a survey completeness estimate or an assessment of any real target. All six Issue 5 lead rows have additional vetoes, so this control does not rescue a candidate. It supports calibrating the seeded SDE gate before using this screen to infer sensitivity.

The environment is isolated under `observations/.venv`. Direct dependencies are NumPy and Astropy. The resolved transitive versions are frozen in `synthetic_seed_requirements.txt`.

```sh
uv venv .venv
uv pip install --python .venv/bin/python -r synthetic_seed_requirements.txt
.venv/bin/python synthetic_seed_control.py --repo /path/to/exoplanets-frontier
```

The required `--repo` option selects the source checkout, independently of the script location.

## Artifact hashes

| File | SHA-256 |
| --- | --- |
| `synthetic_seed_control.py` | `4f3671c95efd5f9d1175dca5e493c39e7276a1ac68232ba2bde066f63395cc26` |
| `synthetic_seed_control.json` | `68b8ffd939e5c76cbe9c5668438b91905fab61f7a8d61427ec66e9652fc448bc` |
| `synthetic_seed_requirements.txt` | `d28050faeb756f09cf94733efc808d9a861ba330f4f0ae05d1378053a27b8c19` |
