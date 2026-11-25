#!/usr/bin/env python3
"""Fail-closed scientific & provenance regression tests for Foundation-0/exoplanets-frontier.

Covers all six hardened non-conflicting pipeline modules:
1. hz_core.py    - FFI SAP/crowd provenance (no flux cloning / hardcoded crowd) & NaN/domain BLS guards.
2. mdwarf_hunt.py - norm_sap preserves invalid/sparse cadences as NaN; spot_amp NaN/domain guards.
3. cool_io.py    - FFI SAP/crowd provenance & calibrated coarse-vs-fine SDE + matched duration grids.
4. ttv_track.py  - 140-trial look-elsewhere BIC penalty, bounded TTV amplitude, & NaN/epoch guards.
5. report_v3.py  - deterministic CONVERGENCE_REPORT.md & canonical sorted JSON provenance manifest.
6. vetter.py     - alias veto without MES>=15 amnesty, NaN fail-closed, close-BEB nfpp floor,
                   multi-epoch blind confirmation, and Pont-beta red-noise blind SNR.
"""

from __future__ import annotations

import contextlib
import csv
import hashlib
import importlib
import io
import json
import math
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _ensure_optional_stubs() -> None:
    """Install lightweight stdlib stubs ONLY for optional plotting/network/table deps absent in minimal envs."""
    if "matplotlib" not in sys.modules:
        try:
            importlib.import_module("matplotlib")
        except ImportError:
            mpl = types.ModuleType("matplotlib")
            mpl.use = lambda *_a, **_k: None
            plt = types.ModuleType("matplotlib.pyplot")
            mpl.pyplot = plt
            sys.modules["matplotlib"] = mpl
            sys.modules["matplotlib.pyplot"] = plt

    if "pandas" not in sys.modules:
        try:
            importlib.import_module("pandas")
        except ImportError:
            pd = types.ModuleType("pandas")

            class _MiniSeries:
                def __init__(self, vals):
                    self._vals = list(vals)

                def value_counts(self):
                    counts: dict = {}
                    for v in self._vals:
                        counts[v] = counts.get(v, 0) + 1
                    return _MiniCounts(counts)

            class _MiniCounts:
                def __init__(self, d):
                    self._d = dict(d)

                def to_dict(self):
                    return dict(self._d)

            class _MiniDataFrame:
                def __init__(self, rows=None):
                    self._rows = [dict(r) for r in (rows or [])]

                def __getitem__(self, key):
                    return _MiniSeries([r.get(key) for r in self._rows])

                def to_csv(self, path, index=False):
                    p = Path(path)
                    if not self._rows:
                        p.write_text("\n", encoding="utf-8")
                        return
                    cols = list(self._rows[0].keys())
                    with p.open("w", encoding="utf-8", newline="") as fh:
                        w = csv.DictWriter(fh, fieldnames=cols)
                        w.writeheader()
                        w.writerows(self._rows)

            pd.DataFrame = _MiniDataFrame
            pd.notna = lambda v: v is not None and not (isinstance(v, float) and math.isnan(v))
            sys.modules["pandas"] = pd

    if "astropy" not in sys.modules:
        try:
            importlib.import_module("astropy")
        except ImportError:
            ap = types.ModuleType("astropy")
            coords = types.ModuleType("astropy.coordinates")
            coords.SkyCoord = type("SkyCoord", (), {})
            io_mod = types.ModuleType("astropy.io")
            fits_mod = types.ModuleType("astropy.io.fits")
            table_mod = types.ModuleType("astropy.table")
            table_mod.vstack = lambda xs: xs
            units_mod = types.ModuleType("astropy.units")
            units_mod.deg = 1.0
            units_mod.arcsec = 1.0
            ts_mod = types.ModuleType("astropy.timeseries")
            ts_mod.BoxLeastSquares = type("BoxLeastSquares", (), {})
            sys.modules["astropy"] = ap
            sys.modules["astropy.coordinates"] = coords
            sys.modules["astropy.io"] = io_mod
            sys.modules["astropy.io.fits"] = fits_mod
            sys.modules["astropy.table"] = table_mod
            sys.modules["astropy.units"] = units_mod
            sys.modules["astropy.timeseries"] = ts_mod

    if "astroquery" not in sys.modules:
        try:
            importlib.import_module("astroquery")
        except ImportError:
            aq = types.ModuleType("astroquery")
            for sub in (
                "astroquery.mast",
                "astroquery.gaia",
                "astroquery.vizier",
                "astroquery.ipac",
                "astroquery.ipac.nexsci",
                "astroquery.ipac.nexsci.nasa_exoplanet_archive",
            ):
                m = types.ModuleType(sub)
                sys.modules[sub] = m
            sys.modules["astroquery"] = aq
            sys.modules["astroquery.mast"].Catalogs = type("Catalogs", (), {"TIMEOUT": 45})
            sys.modules["astroquery.mast"].Observations = type("Observations", (), {"TIMEOUT": 45})
            sys.modules["astroquery.mast"].Tesscut = type("Tesscut", (), {})
            sys.modules["astroquery.gaia"].Gaia = type("Gaia", (), {})
            sys.modules["astroquery.vizier"].Vizier = type("Vizier", (), {})
            sys.modules["astroquery.ipac.nexsci.nasa_exoplanet_archive"].NasaExoplanetArchive = type(
                "NasaExoplanetArchive", (), {}
            )


try:
    importlib.import_module("astropy.timeseries")
    _HAS_REAL_ASTROPY = True
except ImportError:
    _HAS_REAL_ASTROPY = False

_ensure_optional_stubs()

import cool_io  # noqa: E402
import hz_core  # noqa: E402
import mdwarf_hunt  # noqa: E402
import report_v3  # noqa: E402
import ttv_track  # noqa: E402
import vetter  # noqa: E402


def _make_synthetic_lc(period=2.45, t0=1500.4, depth=0.0020, dur=0.08, sectors=(12, 106), noise=2e-4, seed=11):
    rng = np.random.default_rng(seed)
    ts, fs, fes, secs = [], [], [], []
    for idx, sec in enumerate(sectors):
        t_start = 1500.0 + idx * 120.0
        t = np.arange(t_start, t_start + 26.0, 10.0 / 1440.0)
        f = 1.0 + rng.normal(0.0, noise, len(t))
        ph = ((t - t0 + 0.5 * period) % period) - 0.5 * period
        f[np.abs(ph) <= 0.5 * dur] -= depth
        ts.append(t)
        fs.append(f)
        fes.append(np.full(len(t), noise))
        secs.append(np.full(len(t), sec, dtype=int))
    return np.concatenate(ts), np.concatenate(fs), np.concatenate(fes), np.concatenate(secs)


class FakeBLSResult:
    def __init__(self, period, power, depth, transit_time, duration):
        self.period = np.asarray(period, dtype=float)
        self.power = np.asarray(power, dtype=float)
        self.depth = np.asarray(depth, dtype=float)
        self.transit_time = np.asarray(transit_time, dtype=float)
        self.duration = np.asarray(duration, dtype=float)


class MdwarfHuntHardeningTest(unittest.TestCase):
    def test_norm_sap_preserves_nan_nonpositive_and_sparse_sectors_as_nan(self):
        rng = np.random.default_rng(17)
        sap_s1 = 1000.0 + rng.normal(0.0, 2.0, 120)
        sap_s1[5] = np.nan
        sap_s1[12] = 0.0
        sap_s1[20] = -42.0
        sap_s1[33] = np.inf
        # Sparse sector (< 30 cadences) must remain NaN, never fabricated 1.0.
        sap_s2 = np.full(15, 1000.0)
        sap = np.concatenate([sap_s1, sap_s2])
        sec = np.concatenate([np.full(120, 10, dtype=int), np.full(15, 20, dtype=int)])

        out = mdwarf_hunt.norm_sap(sap, sec)
        self.assertEqual(out.shape, sap.shape)
        for bad_idx in (5, 12, 20, 33):
            self.assertTrue(np.isnan(out[bad_idx]), f"index {bad_idx} must remain NaN, got {out[bad_idx]}")
        self.assertTrue(np.all(np.isnan(out[120:])), "sparse sector must remain NaN, not fabricated 1.0")
        good_mask = np.ones(120, dtype=bool)
        good_mask[[5, 12, 20, 33]] = False
        self.assertTrue(np.all(np.isfinite(out[:120][good_mask])))
        self.assertAlmostEqual(float(np.median(out[:120][good_mask])), 1.0, delta=0.05)

    def test_norm_sap_rejects_mismatched_lengths(self):
        with self.assertRaises(ValueError):
            mdwarf_hunt.norm_sap(np.ones(10), np.ones(9, dtype=int))

    def test_spot_amp_rejects_invalid_period_and_filters_nan_flux(self):
        t = np.linspace(0.0, 20.0, 500)
        f = 1.0 + 0.004 * np.sin(2.0 * np.pi * t / 3.5)
        f[10] = np.nan
        f[50] = np.inf
        amp = mdwarf_hunt.spot_amp(t, f, 3.5, 0.0, 0.08)
        self.assertTrue(np.isfinite(amp))
        self.assertAlmostEqual(amp, 0.008, delta=5e-4)

        for bad_p, bad_t0, bad_dur in ((0.0, 0.0, 0.05), (-2.0, 0.0, 0.05), (math.nan, 0.0, 0.05), (3.5, math.nan, 0.05), (3.5, 0.0, -0.1)):
            with self.subTest(p=bad_p, t0=bad_t0, dur=bad_dur):
                # Batch-safe fail-closed: returns math.inf (rejected by amp > 0.48*dep without aborting batch).
                self.assertEqual(mdwarf_hunt.spot_amp(t, f, bad_p, bad_t0, bad_dur), math.inf)
        self.assertEqual(mdwarf_hunt.spot_amp(t, f[:-1], 3.5, 0.0, 0.08), math.inf)


class HzCoreAndCoolIoHardeningTest(unittest.TestCase):
    def test_hz_load_multi_ffi_marks_sap_and_crowd_nan_not_cloned(self):
        t = np.linspace(1500.0, 1525.0, 200)
        ffi_rec = {
            "sec": 12,
            "t": t,
            "f": np.full(200, 1000.0),
            "fe": np.full(200, 0.5),
            "q": np.zeros(200, dtype=int),
            "mdx": np.zeros(200),
            "mdy": np.zeros(200),
        }
        with mock.patch.object(hz_core, "fetch_raw_spoc", return_value=[]), \
             mock.patch.object(hz_core, "discover_all_sectors", return_value=[12]), \
             mock.patch.object(hz_core, "ffi_cutout_lc", return_value=ffi_rec):
            raw, _ = hz_core.load_multi(12345, 10.0, -20.0)
        self.assertEqual(len(raw), 1)
        self.assertTrue(np.all(np.isnan(raw[0]["sap"])))
        self.assertTrue(np.isnan(raw[0]["crowd"]))

        std, _ = hz_core.stack_dual(raw, 1.5, 1.2)
        self.assertIsNotNone(std)
        self.assertTrue(np.all(np.isnan(std["sap"])))
        self.assertTrue(np.isnan(std["crowd"]))

    def test_cool_io_load_cycles_ffi_marks_sap_and_crowd_nan(self):
        t = np.linspace(1500.0, 1525.0, 150)
        ffi_rec = {
            "sec": 12,
            "t": t,
            "f": np.full(150, 500.0),
            "fe": np.full(150, 0.5),
            "q": np.zeros(150, dtype=int),
            "mdx": np.zeros(150),
            "mdy": np.zeros(150),
        }
        with mock.patch.object(cool_io, "fetch_raw_spoc", return_value=[]), \
             mock.patch.object(cool_io, "discover_all_sectors", return_value=[12]), \
             mock.patch.object(cool_io, "ffi_cutout_lc", return_value=ffi_rec):
            raw, _ = cool_io.load_cycles_1_to_8(12345, 10.0, -20.0)
        self.assertEqual(len(raw), 1)
        self.assertTrue(np.all(np.isnan(raw[0]["sap"])))
        self.assertTrue(np.isnan(raw[0]["crowd"]))

    def test_calibrated_bls_sde_uses_broadband_reference_and_recovers_real_bls_transits(self):
        # 1. Synthetic unit check: a strong transit peak lobe inside the +-1.5% fine window
        #    raises pw_fn median/MAD, which MUST NOT suppress SDE when off-window pw_bg is present.
        pw_bg = np.array([1.0 - 0.5 / 1.4826, 1.0, 1.0 + 0.5 / 1.4826] * 100)
        pw_fn_strong = np.concatenate([
            np.array([8.0 - 3.0 / 1.4826, 8.0, 8.0 + 3.0 / 1.4826] * 100),
            np.array([25.0]),
        ])
        sde_bg = cool_io.calibrated_bls_sde(pw_bg, pw_fn_strong)
        self.assertAlmostEqual(sde_bg, (25.0 - 1.0) / 0.5, delta=0.1)

        # Fallback without broadband background masks top-10% peak lobe instead of letting the peak crush SDE.
        pw_fn_only = np.concatenate([np.full(90, 1.0) + np.tile([-0.2, 0.0, 0.2], 30), np.linspace(5.0, 20.0, 10)])
        self.assertGreater(cool_io.calibrated_bls_sde([], pw_fn_only), 15.0)

        # All-NaN or non-positive fine power fails closed to 0.0.
        self.assertEqual(cool_io.calibrated_bls_sde(pw_bg, np.full(20, np.nan)), 0.0)
        self.assertEqual(cool_io.calibrated_bls_sde(pw_bg, np.zeros(20)), 0.0)

        # 2. Real astropy BoxLeastSquares null-vs-injection recovery check (runs whenever astropy is installed).
        if _HAS_REAL_ASTROPY:
            rng_null = np.random.default_rng(101)
            t = np.linspace(1500.0, 1554.0, 4000)
            fe = np.full(len(t), 3e-4)
            f_null = 1.0 + rng_null.normal(0.0, 3e-4, len(t))
            res_null = cool_io.coherent_planet_bls(t, f_null, fe, 3.3, 0.55)
            if res_null is not None:
                self.assertLess(res_null["sde"], 6.0, f"null noise falsely cleared SDE gate: {res_null['sde']}")

            for depth_ppm in (300.0, 1000.0, 2000.0, 5000.0):
                with self.subTest(depth_ppm=depth_ppm):
                    rng = np.random.default_rng(int(depth_ppm) + 7)
                    f_inj = 1.0 + rng.normal(0.0, 3e-4, len(t))
                    ph = ((t - 1501.1 + 0.5 * 3.3) % 3.3) - 0.5 * 3.3
                    f_inj[np.abs(ph) <= 0.04] -= depth_ppm * 1e-6
                    res_inj = cool_io.coherent_planet_bls(t, f_inj, fe, 3.3, 0.55)
                    self.assertIsNotNone(res_inj)
                    self.assertGreater(res_inj["sde"], 8.0, f"depth={depth_ppm} ppm crushed to SDE={res_inj['sde']}")
                    self.assertAlmostEqual(res_inj["p"], 3.3, delta=0.01)

    def test_coherent_bls_guards_invalid_inputs_and_matches_duration_grids(self):
        t, f, fe, _ = _make_synthetic_lc(period=4.0, t0=1501.0, depth=0.002, dur=0.08, sectors=(12, 13))
        self.assertIsNone(hz_core.coherent_seed_bls(t, f, fe, math.nan, 0.5))
        self.assertIsNone(hz_core.coherent_seed_bls(t, f, fe, -1.0, 0.5))
        self.assertIsNone(cool_io.coherent_planet_bls(t, f, fe, math.nan, 0.5))
        self.assertIsNone(cool_io.coherent_planet_bls(t, f, fe, 0.0, 0.5))

        seen_calls = []

        class RecordingBLS:
            def __init__(self, tt, ff, dy=None):
                self.tt = tt

            def autoperiod(self, durs, minimum_period, maximum_period, frequency_factor=2.5):
                return np.linspace(minimum_period, maximum_period, 50)

            def power(self, periods, durations, objective="snr"):
                p_arr = np.asarray(periods, dtype=float)
                d_arr = np.asarray(durations, dtype=float)
                seen_calls.append((p_arr.copy(), d_arr.copy()))
                # Return zero power to verify clean None fail-closed on zero/all-NaN power.
                return FakeBLSResult(
                    p_arr,
                    np.full(len(p_arr), np.nan),
                    np.zeros(len(p_arr)),
                    np.full(len(p_arr), 1501.0),
                    np.full(len(p_arr), float(d_arr[0])),
                )

        ts_mod = sys.modules["astropy.timeseries"]
        with mock.patch.object(ts_mod, "BoxLeastSquares", RecordingBLS):
            self.assertIsNone(hz_core.coherent_seed_bls(t, f, fe, 4.0, 0.5))
            seen_calls.clear()
            self.assertIsNone(cool_io.coherent_planet_bls(t, f, fe, 4.0, 0.5))

        # cool_io must call power twice (broadband bg excluding [p_lo, p_hi] + fine) with IDENTICAL durations.
        self.assertEqual(len(seen_calls), 2)
        bg_p, bg_d = seen_calls[0]
        fn_p, fn_d = seen_calls[1]
        self.assertTrue(np.array_equal(bg_d, fn_d), "broadband and fine BLS must use identical duration grids")
        p_lo, p_hi = 0.985 * 4.0, 1.015 * 4.0
        self.assertFalse(np.any((bg_p >= p_lo) & (bg_p <= p_hi)), "broadband bg must exclude fine signal window")


class TtvTrackHardeningTest(unittest.TestCase):
    def test_fit_ephemeris_applies_look_elsewhere_penalty_and_bounds_collinear_amp(self):
        # 1. Exact look-elsewhere penalty verification on a clean 12-epoch sinusoid.
        epochs = np.arange(0, 24, 2, dtype=float)
        p_sup = 14.0
        true_amp = 14.0
        dts = true_amp * np.sin(2.0 * np.pi * epochs / p_sup)
        rows = [{"epoch": int(e), "dt_min": float(y), "err_min": 2.0, "snr_peak": 6.0} for e, y in zip(epochs, dts)]
        fit = ttv_track.fit_ephemeris_models(rows, rows, 3.2, 1500.0)
        self.assertEqual(fit["class"], "DYNAMICAL_TTV_CANDIDATE")
        self.assertAlmostEqual(fit["ttv_amp_m"], true_amp, delta=1.5)

        # Recompute expected dBIC_ttv with 2*ln(140) look-elsewhere penalty vs unpenalized.
        n = len(rows)
        w = np.full(n, 1.0 / 4.0)
        Xl = np.column_stack([np.ones(n), epochs])
        cl, *_ = np.linalg.lstsq(Xl * np.sqrt(w)[:, None], dts * np.sqrt(w), rcond=None)
        bicl = float(np.sum(w * (dts - Xl @ cl) ** 2)) + 2.0 * math.log(n)
        unpenalized_upper = bicl - 5.0 * math.log(n)
        penalized_upper = bicl - (5.0 * math.log(n) + 2.0 * math.log(140))
        self.assertLessEqual(fit["dBIC_ttv"], round(penalized_upper, 1) + 0.2)
        self.assertLess(fit["dBIC_ttv"], unpenalized_upper - 9.0)

        # 2. Clustered-epoch collinearity attack: 2 pairs of close epochs separated by a huge gap
        # with tiny 2-min linear residuals must NOT explode ttv_amp_m to hundreds/thousands of minutes.
        clustered_e = [0, 1, 2, 300, 301, 302]
        clustered_y = [1.5, -1.5, 0.5, -1.5, 1.5, -0.5]
        c_rows = [{"epoch": e, "dt_min": y, "err_min": 2.0, "snr_peak": 5.0} for e, y in zip(clustered_e, clustered_y)]
        c_fit = ttv_track.fit_ephemeris_models(c_rows, c_rows, 2.5, 1500.0)
        self.assertLessEqual(c_fit["ttv_amp_m"], 25.0)
        self.assertNotEqual(c_fit["class"], "DYNAMICAL_TTV_CANDIDATE")

    def test_fit_ephemeris_and_extract_reject_nan_and_degenerate_inputs(self):
        self.assertEqual(ttv_track.extract_transit_times([], -1.0, 1500.0, 0.08), ([], []))
        self.assertEqual(ttv_track.extract_transit_times([], math.nan, 1500.0, 0.08), ([], []))
        for bad_p, bad_t0 in ((0.0, 1500.0), (math.nan, 1500.0), (2.5, math.nan)):
            inv = ttv_track.fit_ephemeris_models([], [], bad_p, bad_t0)
            self.assertEqual(inv["class"], "INSUFFICIENT_EPOCHS")
            self.assertEqual(inv["source"], "invalid")
            self.assertTrue(np.isfinite(inv["P_ref"]) and np.isfinite(inv["T0_ref_BJD"]))

        # NaN rows and identical epochs must return INSUFFICIENT_EPOCHS cleanly without NaN propagation.
        bad_rows = [
            {"epoch": 5, "dt_min": 1.0, "err_min": 2.0, "snr_peak": 5.0},
            {"epoch": 5, "dt_min": 1.5, "err_min": 2.0, "snr_peak": 5.0},
            {"epoch": 10, "dt_min": math.nan, "err_min": 2.0, "snr_peak": 5.0},
            {"epoch": 12, "dt_min": 1.0, "err_min": 0.0, "snr_peak": 5.0},
        ]
        res = ttv_track.fit_ephemeris_models(bad_rows, [], 2.5, 1500.0)
        self.assertEqual(res["class"], "INSUFFICIENT_EPOCHS")
        self.assertTrue(np.isfinite(res["P_ref"]))


class ReportV3DeterminismTest(unittest.TestCase):
    def test_save_phase3_outputs_is_bit_deterministic_and_rejects_nan_json(self):
        bench = {
            "despotted_P_d": 1.239126,
            "despotted_MES": 14.33,
            "null_false_alarms_ge_9p5": 0,
            "null_trials": 10,
        }
        cat_a = [{"dr4_class": "GOLD_LOW_BLEND"}, {"dr4_class": "PLATINUM_SOLO_VERIFIED"}]
        cat_b = [{"dr4_class": "PLATINUM_SOLO_VERIFIED"}, {"dr4_class": "GOLD_LOW_BLEND"}]
        plat = [{
            "tic": 100, "gaia_dr3_id": 1, "upstream_tier": "S", "P_d": 3.5, "Rp_Re": 2.1,
            "Teq_K": 310.0, "is_hz_temperate": True, "TSM": 42.0, "Tmag": 11.2, "ruwe": 1.01,
            "plausible_beb_65as": 0, "dr4_class": "PLATINUM_SOLO_VERIFIED",
        }]
        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            p1, p2 = Path(d1), Path(d2)
            with contextlib.redirect_stdout(io.StringIO()):
                report_v3.save_phase3_outputs(p1, p1 / "figures", 1000.0, 108, bench, cat_a, plat, [], [], [], [], [], [])
                report_v3.save_phase3_outputs(p2, p2 / "figures", 500.0, 108, bench, cat_b, plat, [], [], [], [], [], [])

            md1 = (p1 / "CONVERGENCE_REPORT.md").read_bytes()
            md2 = (p2 / "CONVERGENCE_REPORT.md").read_bytes()
            self.assertEqual(md1, md2, "CONVERGENCE_REPORT.md must be bit-identical regardless of elapsed_s or class order")
            self.assertNotIn(b" s)", md1)

            c1 = json.loads((p1 / "convergence_certificates.json").read_text(encoding="utf-8"))
            c2 = json.loads((p2 / "convergence_certificates.json").read_text(encoding="utf-8"))
            self.assertEqual(c1["files_sha256"], c2["files_sha256"])
            self.assertEqual(c1["manifest_sha256"], c2["manifest_sha256"])
            self.assertEqual(c1["code_sha256"], c2["code_sha256"])
            self.assertEqual(c1["manifest_sha256"], report_v3.canonical_files_digest(c1["files_sha256"]))
            self.assertEqual(
                c1["files_sha256"]["CONVERGENCE_REPORT.md"],
                hashlib.sha256(md1).hexdigest(),
            )

            # NaN in benchmark must fail closed BEFORE writing CSV/MD/JSON (leaving existing files intact).
            prev_md = (p1 / "CONVERGENCE_REPORT.md").read_bytes()
            prev_cert = (p1 / "convergence_certificates.json").read_bytes()
            bad_bench = dict(bench, despotted_MES=float("nan"))
            with self.assertRaises(ValueError):
                with contextlib.redirect_stdout(io.StringIO()):
                    report_v3.save_phase3_outputs(p1, p1 / "figures", 1000.0, 108, bad_bench, cat_a, plat, [], [], [], [], [], [])
            self.assertEqual((p1 / "CONVERGENCE_REPORT.md").read_bytes(), prev_md)
            self.assertEqual((p1 / "convergence_certificates.json").read_bytes(), prev_cert)


class VetterHardeningTest(unittest.TestCase):
    def _clean_pk(self, period=2.45, t0=1500.4, dur=0.08, dep=0.0020, mes=22.0, mes_c=20.0, sde=14.0):
        return {"p": period, "t0": t0, "dur": dur, "dep": dep, "sc": 2e-4, "mes": mes, "mes_c": mes_c, "sde": sde, "ntr": 12}

    def test_alias_veto_fires_even_when_mes_ge_15(self):
        t, f, fe, sec = _make_synthetic_lc(period=1.001, t0=1500.2, depth=0.0020, dur=0.06, sectors=(12, 106))
        pk = self._clean_pk(period=1.001, t0=1500.2, dur=0.06, dep=0.0020, mes=28.0, mes_c=25.0, sde=16.0)
        with mock.patch.object(vetter, "search_bls", return_value=None):
            rep = vetter.vet_signal(101, pk, t, f, fe, sec, "TEST", meta={"rstar": 0.55, "teff": 4000, "tmag": 11.0, "ruwe": 1.0})
        self.assertFalse(rep["passed"])
        self.assertIn("alias", rep["vetoes"])

    def test_vet_signal_fails_closed_on_nan_and_nonpositive_pk(self):
        t, f, fe, sec = _make_synthetic_lc()
        with mock.patch.object(vetter, "search_bls", return_value=None):
            for field, val in (("mes_c", math.nan), ("sde", math.nan), ("dep", math.nan), ("dep", -0.001), ("p", 0.0), ("mes", math.nan)):
                with self.subTest(field=field, val=val):
                    pk = self._clean_pk()
                    pk[field] = val
                    rep = vetter.vet_signal(102, pk, t, f, fe, sec, "TEST")
                    self.assertFalse(rep["passed"])
                    self.assertEqual(rep["tier"], "rejected")

    def test_close_unresolved_beb_not_discounted_to_018_and_blocks_tier_s_val(self):
        t, f, fe, sec = _make_synthetic_lc(period=2.45, t0=1500.4, depth=0.0012, dur=0.07, sectors=(12, 106), noise=1.5e-4)
        pk = self._clean_pk(period=2.45, t0=1500.4, dur=0.07, dep=0.0012, mes=22.0, mes_c=18.0, sde=15.0)
        close_nb = [{"sep_arcsec": 1.5, "tmag": 12.0}]
        with mock.patch.object(vetter, "search_bls", return_value=None):
            rep = vetter.vet_signal(
                103, pk, t, f, fe, sec, "TEST",
                meta={"rstar": 0.55, "teff": 4000, "tmag": 11.0, "ruwe": 1.0},
                nbs=close_nb,
                tpf_off=-1.0,
            )
        # Without 0.18 collapse, loc_w = exp(-0.5*(1.5/13)^2) ~ 0.993 -> nfpp ~ 0.00262 (vs buggy 0.00047).
        self.assertGreater(rep["NFPP"], 0.0020)
        self.assertNotEqual(rep["tier"], "Tier S-Val")
        self.assertEqual(rep["tier"], "Tier N")

    def test_nan_neighbor_and_missing_sap_crowdsap_fail_closed(self):
        t, f, fe, sec = _make_synthetic_lc(period=2.45, t0=1500.4, depth=0.0012, dur=0.07, sectors=(12, 106), noise=1.5e-4)
        pk = self._clean_pk(period=2.45, t0=1500.4, dur=0.07, dep=0.0012, mes=22.0, mes_c=18.0, sde=15.0)
        with mock.patch.object(vetter, "search_bls", return_value=None):
            # 1. Missing CROWDSAP (NaN) and all-NaN SAP must emit no_crowdsap + no_sap notes and cap below Tier S-Val.
            rep_ffi = vetter.vet_signal(
                104, pk, t, f, fe, sec, "TEST_FFI",
                sap=np.full(len(f), np.nan),
                crowd=math.nan,
                meta={"rstar": 0.55, "teff": 4000, "tmag": 11.0, "ruwe": 1.0},
            )
            self.assertTrue(rep_ffi["passed"])
            self.assertIn("no_crowdsap", rep_ffi["notes"])
            self.assertIn("no_sap", rep_ffi["notes"])
            self.assertEqual(rep_ffi["tier"], "Tier S")

            # 2. Neighbor with NaN separation or NaN tmag must NOT be silently skipped (fail-closed -> Tier N, NFPP > 0.002).
            for bad_nb in ([{"sep_arcsec": math.nan, "tmag": 12.0}], [{"sep_arcsec": 5.0, "tmag": math.nan}]):
                with self.subTest(bad_nb=bad_nb):
                    rep_nb = vetter.vet_signal(
                        105, pk, t, f, fe, sec, "TEST_NB",
                        sap=f.copy(),
                        crowd=0.98,
                        meta={"rstar": 0.55, "teff": 4000, "tmag": 11.0, "ruwe": 1.0},
                        nbs=bad_nb,
                    )
                    self.assertGreater(rep_nb["NFPP"], 0.0020)
                    self.assertEqual(rep_nb["tier"], "Tier N")
                    self.assertIn("unknown_neighbor_meta", rep_nb["notes"])

    def test_blind_check_sector_requires_two_epochs_applies_pont_beta_and_guards_nan(self):
        # 1. Single-epoch dip (P = 22.98 d in a 25-day sector -> n_cov == 1) must NOT be CONFIRMED_ON_SCHEDULE.
        t1, f1, _, _ = _make_synthetic_lc(period=22.98, t0=1512.0, depth=0.0020, dur=0.10, sectors=(106,), noise=1.5e-4)
        b1 = vetter.blind_check_sector(201, "S", 22.98, 1512.0 + vetter.BTJD0, 2000.0, 2.4, 106, t1, f1)
        self.assertEqual(b1["n_cov"], 1)
        self.assertNotEqual(b1["verdict"], "CONFIRMED_ON_SCHEDULE")
        self.assertEqual(b1["verdict"], "MARGINAL_RECOVERY")

        # 2. Multi-epoch dip (P = 2.5 d in 25-day sector -> n_cov >= 2) confirms cleanly.
        t2, f2, _, _ = _make_synthetic_lc(period=2.5, t0=1501.0, depth=0.0020, dur=0.08, sectors=(106,), noise=1.5e-4)
        b2 = vetter.blind_check_sector(202, "S", 2.5, 1501.0 + vetter.BTJD0, 2000.0, 1.92, 106, t2, f2)
        self.assertGreaterEqual(b2["n_cov"], 2)
        self.assertEqual(b2["verdict"], "CONFIRMED_ON_SCHEDULE")

        # 3. Red-noise correlated light curve must inflate err_ppm via Pont beta and reduce snr.
        red_wave = 0.0008 * np.sin(2.0 * np.pi * t2 / 1.37)
        b2_red = vetter.blind_check_sector(202, "S", 2.5, 1501.0 + vetter.BTJD0, 2000.0, 1.92, 106, t2, f2 + red_wave)
        self.assertGreater(b2_red["err_ppm"], b2["err_ppm"] * 1.5)
        self.assertLess(b2_red["snr"], b2["snr"])

        # 4. Non-finite parameters and NaN-contaminated flux fail closed without NaN output.
        f2_nan = f2.copy()
        f2_nan[::25] = np.nan
        b2_nan = vetter.blind_check_sector(202, "S", 2.5, 1501.0 + vetter.BTJD0, 2000.0, 1.92, 106, t2, f2_nan)
        self.assertTrue(np.isfinite(b2_nan["snr"]))
        self.assertEqual(b2_nan["verdict"], "CONFIRMED_ON_SCHEDULE")

        b_bad = vetter.blind_check_sector(203, "S", math.nan, 1501.0, 2000.0, 1.92, 106, t2, f2)
        self.assertEqual(b_bad["verdict"], "NO_DATA")


if __name__ == "__main__":
    unittest.main()
