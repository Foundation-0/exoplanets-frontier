# Issue 5 interpretation and sensitivity controls

The [evidence audit](REVIEW.md) explains the saved rejection counts and the limits of catalog labels. Its [independent review](PUBLIC_READY_REVIEW.md) and exact-file manifest cover that original increment.

The subsequent [synthetic transit control](synthetic/SYNTHETIC_SEED_CONTROL.md) executes two fixed injections through the unchanged seeded search. The 25-day signal is recovered with corrected MES 49.76 but SDE 7.23, below the configured 7.5 gate. The 17-day control passes that gate.

This supplies a concrete reason to calibrate the seeded threshold before interpreting sensitivity. The test does not execute the full pipeline or determine a real target's disposition. All six original candidate rows also have other vetoes.

The lead reviewed the extraction, fixed settings and interpretation, then independently reran the control in the frozen environment. Its complete JSON output matched the supplied result byte-for-byte. No thresholds or original results were changed. This replay is separate from the earlier review's exact-file verdict.

Run the commands in the synthetic report from its directory, supplying the repository root to `--repo`. The source hashes are checked before computation. The original observation files remain unchanged.

The later [C-07 review](c07/NEW_C07_REVIEW.md) checks the new negative-control provenance claim at main revision `3636257`. All six signals match a pinned upstream inverted-light-curve table, and five carry upstream vetting failures. The lead independently reproduced all six joins. `CLAIMS.md` now separates those matches and screen rejection from physical classification and blinded validation.

The C-07 source manifest supplies immutable URLs and hashes for the fetched source files. Its `saved_as` names describe the review's source cache. Full upstream files are not bundled here. The compact matched-row JSON preserves the relevant six records.
