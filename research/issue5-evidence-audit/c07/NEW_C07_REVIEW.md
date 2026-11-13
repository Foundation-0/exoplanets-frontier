# Review of the revised C-07 claim

C-07 at `3636257c9852fcd805a84734a18516526cd1d746` needs a narrower interpretation. The negative-control connection has primary-source support. The claims of exclusive historical origin, blinded rejection and recovery of genuine planets exceed the evidence. The proposed replacement preserves the zero-survivor result and the upstream table matches.

## What is established

All six Issue #5 rows record only `validation/inverted_lightcurves_vetted.csv` in `src_files`. Each TIC/period pair matches a row in that upstream file at revision `e2eee928fd621f75c1621f6b4206d34bb222ae59`, within 0.00001 day of the rounded frontier period. The table has Git blob SHA `3cbc987a0b700d94db756769fd23ff6e3085aaaf`. Five of the six matched upstream rows already have `n_fail > 0`. TIC 37904598 has `n_fail = 0`.

| TIC | Frontier period (days) | Upstream period (days) | Upstream failed tests |
| --- | --- | --- | --- |
| 181015118 | 25.26197 | 25.26196970935829 | 7 |
| 181015118 | 18.31799 | 18.3179925909046 | 6 |
| 246897668 | 24.75243 | 24.7524323734834 | 7 |
| 37904598 | 17.04918 | 17.0491756357338 | 0 |
| 146636926 | 26.05645 | 26.056444564621607 | 6 |
| 146636926 | 23.92026 | 23.920263616157246 | 4 |


The upstream [film preparation code](https://github.com/eumemic/exoplanets/blob/e2eee928fd621f75c1621f6b4206d34bb222ae59/video/tools/prep_lc.py#L121-L148) reads this exact table, selects TIC 37904598, inverts its light curve and describes the signal as a false alarm.

The upstream [README](https://github.com/eumemic/exoplanets/blob/e2eee928fd621f75c1621f6b4206d34bb222ae59/README.md#L897-L902) explains the inverted-light-curve control and its limits. This supports the table's control interpretation beyond its filename alone.

## What needs correction

| C-07 interpretation | Evidence boundary |
| --- | --- |
| All six originated solely in the control file | The saved ledger reports one matching path. `scan_eumemic_origins` searches `.csv` files, excludes names containing `funnel` or `sector`, silently ignores read failures and drops upstream `n_fail`, `fails`, `known_status` and `known_detail`. It has no recorded upstream revision. It cannot establish exhaustive historical origin. |
| Independent provenance proof | `audit_target` matches TIC and periods within 0.20 day, then falls back to any row for the TIC. That operation finds associations. It does not reconstruct a signal's generation or selection history. |
| 100 percent blind rejection of negative controls | The six periods are explicit seeds. The search refits epoch and duration on ordinary archive flux. It does not replay the original inverted control curves or document blinding. Six failures among selected signals do not calibrate general rejection performance. |
| Genuine warm giant transits recovered | Both comparison rows fail `survived_hz_planet`. TIC 166848769 has radius, inverted-flux and secondary-eclipse vetoes. TIC 224279805 has radius and odd/even vetoes. Strong periodic dimming and small reported image offsets do not independently establish planets. |
| Every signal failed odd/even and inverted symmetry | The vetoes differ by row. TIC 246897668 lacks an odd/even veto. Both TIC 146636926 rows lack an inverted-flux veto. |
| Every host's apparent signal collapsed after spot removal | For TIC 37904598 the reported MES rises from 4.08 to 4.46. It stays below threshold. The other rows lose MES, but this alone is not a physical classification. |

No new transit analysis or full upstream-tree scan was performed for this review. The six matching rows support a specific table association. They do not recover the unrecorded upstream checkout used by the original run. Existing audit and synthetic-control reviews remain unchanged.

## Pinned sources

- [C-07 under review](https://github.com/Foundation-0/exoplanets-frontier/blob/3636257c9852fcd805a84734a18516526cd1d746/CLAIMS.md)
- [Saved frontier result rows](https://github.com/Foundation-0/exoplanets-frontier/blob/883e28a7cecde05fdb6679bf574dbe78668332e7/results/hz_issue5_audit.csv)
- [Provenance scanner](https://github.com/Foundation-0/exoplanets-frontier/blob/883e28a7cecde05fdb6679bf574dbe78668332e7/hz_vet.py#L28-L46)
- [Signal association and seeded analysis](https://github.com/Foundation-0/exoplanets-frontier/blob/883e28a7cecde05fdb6679bf574dbe78668332e7/hz_core.py#L72-L108)
- [Upstream inverted-light-curve table](https://github.com/eumemic/exoplanets/blob/e2eee928fd621f75c1621f6b4206d34bb222ae59/validation/inverted_lightcurves_vetted.csv)

`source_manifest.json` records the fetched URLs, revisions, Git blob hashes and SHA-256 values. Every saved source file matches its GitHub blob SHA. `c07_matched_rows.json` preserves the six joins and the upstream statuses that the frontier scanner omits.
