# Independent review of the Issue 5 evidence audit

Verdict: **PUBLIC_READY** for the exact files below. The audit supports a correction to how the saved screen is interpreted. It does not validate the transit pipeline or establish any planet's physical status. No blocking scientific or implementation finding remains in the proposed increment.

The source base is `883e28a7cecde05fdb6679bf574dbe78668332e7`. The tracked source and result files reviewed here match that base. The three proposed audit files were reviewed as uncommitted additions.

The bounded command `python3 research/issue5-evidence-audit/audit_snapshot.py` exits successfully. Its parsed output matches the saved `evidence.json`. All 40 historical result-manifest entries match their SHA-256 values. The extracted original loader accepts a synthetic upstream candidate with empty external catalog stubs. This demonstrates the field's possible origin, not the origin of any real row.

| Claim | Independent check |
| --- | --- |
| Zero surviving proposed signals | Six rows, four unique TIC identifiers, all rejected in the saved CSV. |
| Zero surviving comparison signals | Five rows, all rejected. Their labels do not establish independent confirmed-planet controls. |
| Mixed catalog membership | `mast_data.py:49-71` combines TOIs, the composite planet table, eclipsing binaries and upstream candidate CSVs without preserving origin or disposition. |
| First stored period | `hz_core.py:102` uses the first period for the target, without a signal-specific source match. |
| No explicit habitable-zone gate | `hz_core.py:108` copies acceptance after generic, spot and position tests. `vetter.py` computes temperature but does not gate on a habitable-zone boundary. |
| Comparison temperatures and radii | Reported temperatures are 410-854 K. Two radii exceed eight Earth radii. Both rows have additional vetoes. |
| Sector counts | The revised note correctly calls 3-7 the reported loaded-sector count. `stack_dual` can discard sectors before analysis. |
| Independent prediction | The seeded fit searches within 1.8 percent of the proposed period and refits epoch and duration. The note correctly avoids calling this a held-out ephemeris prediction. |

The mock uses a deliberately narrow pandas substitute and empty remote tables. It is sufficient for the demonstrated upstream-only path. It does not reproduce remote catalog behavior. Hash equality verifies bytes, not scientific validity. No light curves were downloaded, and no transit analysis was rerun.

A cheap next diagnostic is available before a broader injection campaign. All six lead rows contain a low-SDE veto. `coherent_seed_bls` computes this contrast statistic within a narrow period window, while `search_bls` normally computes it across a broader first-stage search. Both feed the same 7.5 threshold in the generic vetter. 

A small deterministic synthetic transit test through the seeded function could test whether valid long-period signals clear that threshold. This is a calibration question, not a demonstrated defect or a rescue of these six rows. Every lead has other vetoes too. The follow-up numerical test is separate from this publication verdict and will use an isolated environment.

The cross-project connection in the note remains limited to preserving evidence behind labels. Nothing reviewed here supports a shared physical mechanism with gravity or geometry.

## Exact reviewed bytes

Paths are relative to the source repository. A change to any proposed audit file requires a fresh check of the affected claim.

| File | SHA-256 |
| --- | --- |
| `research/issue5-evidence-audit/REVIEW.md` | `28e3cbe544a4c69d53c5dba13a30d2a41b1039d3b2c9180b9c990f3c84eb0e71` |
| `research/issue5-evidence-audit/audit_snapshot.py` | `cffff78a1289088ba221033b57cbb500f65a6a2b4cb41141b9f84d6bf9666001` |
| `research/issue5-evidence-audit/evidence.json` | `6b56ef99ffed9befff30f8f66ef84186e69c72dee6a3d657dcfefaf427200f94` |
| `hz_core.py` | `7977d141a473281c328c944cd5772670c205a059032e16af3e37a64e2e6c9a7e` |
| `hz_vet.py` | `32b809fd34a4d0d126624ad05ef786f6e70bf1e08d33b49a0fc87f27c76977f2` |
| `mast_data.py` | `9cf8770d0b6bc550b87ca67c0eca537edf79b8b37ce6fb43b5eca5dc382ec240` |
| `vetter.py` | `cf910ed563e14831a6cfaff2d8309cdc05848fbf84928e0b2a73f535124b982f` |
| `bls_search.py` | `a5bb48b81db07d1d633b39c72cfea340e82e4f8c6f46b8ba8a755f64f3f40977` |
| `converge_io.py` | `36be87209c978115540ab88f7d96de35c4f95ea584cc8d4b484338c15f72109b` |
| `results/hz_issue5_audit.csv` | `28e36fbd280ddefc69a3adb9e4089e99a034c395230a16c7c39a9331781dd14c` |
| `results/HZ_ISSUE5_REPORT.md` | `506d1577c2e4ec08c905b34a79e9f50479b1e3643c09b7b8698828cf5697b7c4` |
| `results/convergence_certificates.json` | `dc1c7b9f5a4f92a418548b34c4939f3782e31ce4bdc308ac5ffc226c50e15a38` |

Reviewed by gpt-forge, independent subagent.
