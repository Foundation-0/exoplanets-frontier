# What the Issue 5 screen establishes

The committed multi-sector screen rejects all six proposed temperate signals, which belong to four distinct stars. It also rejects five comparison signals. This is useful candidate triage. The screen does not establish that every target lacks a planet, measure survey completeness, or independently validate a planetary system.

This review binds to main revision `883e28a7cecde05fdb6679bf574dbe78668332e7`. The original results and analysis code are unchanged. The accompanying offline audit checks the ledger and result hashes and executes a mocked provenance control. It does not download light curves or repeat the full transit analysis.

## Counts and scope

The [committed report](../../results/HZ_ISSUE5_REPORT.md) records zero survivors among six Issue 5 rows and five comparison rows. The first group contains four unique TIC identifiers. Each row reports three to seven loaded sectors. The processing can discard sectors. The ledger does not give each branch's retained sector count, and S1–S108 names the search era.

The [runner](../../hz_core.py) searches within 1.8% of each proposed period, refits the epoch and duration, and compares standard and starspot-filtered data. A separate blind search is recorded. This is a multi-sector reanalysis of proposed signals. It is not a new held-out prediction at a frozen ephemeris.

The earlier TOI-540 period-recovery benchmark and ten scrambled-null trials remain useful checks of their particular procedure. They do not calibrate completeness or false-rejection rates for these longer-period signals, different stars, and all gates in the Issue 5 analysis.

## Two labels need careful interpretation

`known_nasa` is computed as membership in `kt`. The [catalog loader](../../mast_data.py) populates that set from NASA TOIs, the composite planet table, an eclipsing-binary catalog, and upstream candidate CSVs. It does not retain a source or disposition in the set.

The accompanying control supplies empty external catalogs and one synthetic upstream candidate. The unchanged loader still adds that identifier to `kt` and its period to `kp`. Therefore `known_nasa=True` cannot establish NASA membership, independent confirmation, or even independence from the candidate list being tested.

This is a counterexample to the field's interpretation. It does not determine which source added each real target. Likewise, `known_P_d` selects the first stored period without source provenance or a required planet match. Agreement with that value need not be an independent ephemeris check.

The [NASA TOI definitions](https://exoplanetarchive.ipac.caltech.edu/docs/API_TOI_columns.html) distinguish confirmed and known planets from candidates, ambiguous candidates, false positives and false alarms. A future crossmatch must retain that distinction and associate a disposition with the particular signal.

`survived_hz_planet` copies the generic vetter's acceptance flag after additional spot and position gates. Neither it nor the [generic vetter](../../vetter.py) applies an explicit habitable-zone boundary. Equilibrium temperature is computed for reporting. Thus the flag describes screen survival, rather than verified habitability or planetary status.

The five comparison rows have reported equilibrium temperatures from 410 to 854 K. Two fail the configured upper radius cut of eight Earth radii, and both also fail other gates. These are selected comparison signals, not an established positive-control sample from which one can infer a zero recovery rate for known planets.

## The next useful test

Before expanding this search, freeze a source-qualified control set with independent dispositions and ephemerides. Include true planets within the screen's intended period, radius, signal strength and stellar-variability range, plus characterized nonplanetary signals. Separate planet recovery from selection for a temperate follow-up program.

Inject signals into actual light curves before detrending, then measure recovery through the full standard and starspot-filtered pipeline. Preserve the exact periods, epochs, sectors, quality masks, catalog revisions and veto reasons. This can test whether aggressive spot removal or multiple gates erase valid transits.

These controls can distinguish useful rejection of contaminated signals from an overly restrictive screen. The present zero-survivor result alone cannot. No new planet discovery, statistical planet validation, or universal false-positive verdict follows from this audit.

## Reproduce the bounded audit

From the repository root:

```sh
python3 research/issue5-evidence-audit/audit_snapshot.py
```

The output records source hashes, row counts, all existing result-manifest checks and the mocked loader counterexample. Compare it with [the saved output](evidence.json). A matching hash establishes byte identity, not scientific validity. The mock controls only catalog-membership semantics, not the behavior of remote catalogs.

The broader connection to the gravity and geometry work is methodological: a summary label must retain the information needed for its intended claim. This shared issue does not identify a common physical mechanism.
