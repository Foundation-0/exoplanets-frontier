"""Markdown discovery ledger writer and end-to-end synthetic self-tests."""
from __future__ import annotations
import unittest
from pathlib import Path
import numpy as np
from bls_search import BTJD0, RSUN_RE, detrend_lc, exp_dur, search_bls
from vetter import blind_check_sector, vet_signal

def write_ledger(out_dir: Path, summ: dict, bres: list, vpass: list, vall: list) -> None:
    s_max, new_s, vc = summ["s_max"], summ["new_sectors"], summ["blind_counts"]
    md = [f"# TESS Sectors 1-{s_max} Verification and Discovery Ledger", "",
          f"Empirical measurements from NASA STScI MAST 2-minute SPOC light curves and Target Pixel Files (`{summ['utc']}`, `{summ['sec_s']} s`).", "",
          f"- **Highest MAST SPOC sector indexed**: `S{s_max}` (new sectors tested: `{new_s}`)",
          f"- **Out-of-sample (`v3.0`, S106/S107) pairs tested**: `{summ['blind_pairs']}` -> `{vc}`",
          f"- **Uncatalogued S106-S107 dwarf signals found / stacked across S1-S{s_max}**: `{summ['unpromoted_tce_pool_total']}` found, `{summ['tce_stacked']}` stacked",
          f"- **Hosts searched for 0.155-0.50 d USP siblings**: `{summ['usp_hosts']}` | **Passed all 8 gates + TPF**: `{summ['new_passed']}`", "",
          "## 1. Out-of-Sample Blind Verification (Sectors 106-107)", "",
          "| TIC | Tier | Sec | Period (d) | Pred (ppm) | Obs (ppm) | SNR | Offset (m) | Cent | Verdict |",
          "| :--- | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |"]
    for r in [x for x in bres if x["verdict"] in ("CONFIRMED_ON_SCHEDULE", "TIMING_OFFSET_TTV", "MARGINAL_RECOVERY")][:35]:
        md.append(f"| `{r['tic']}` | `{r['tier']}` | S{r['sector']} | {r['P_d']:.5f} | {r['pred_ppm']:.0f} | {r['obs_ppm']:.0f}+/-{r['err_ppm']:.0f} | {r['snr']:.2f} | {r['dt_min']:+.1f} | {r['cent_sig']:.2f} | `{r['verdict']}` |")
    nd = [x for x in bres if x["verdict"] == "NON_DETECTION_4SIGMA"]
    if nd:
        md += ["", "### High-Sensitivity Non-Detections (>= 4 Sigma Deficit)", "",
               "| TIC | Tier | Sec | Period (d) | Pred (ppm) | Obs (ppm) | N | Verdict |", "| :--- | :--- | :--- | ---: | ---: | ---: | ---: | :--- |"]
        for r in nd[:20]:
            md.append(f"| `{r['tic']}` | `{r['tier']}` | S{r['sector']} | {r['P_d']:.5f} | {r['pred_ppm']:.0f} | {r['obs_ppm']:.0f}+/-{r['err_ppm']:.0f} | {r['n_cov']} | `{r['verdict']}` |")
    md += ["", "## 2. New Uncatalogued Planet Candidates (All 8 Gates + TPF Passed)", ""]
    if vpass:
        md += ["| TIC | Lane | Period (d) | T0 (BJD) | Depth (ppm) | Dur (h) | Rp (Re) | MES | Sec | OE | Off (\") | RUWE | Teq (K) | TSM | FPP | Tier |",
               "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- |"]
        for c in vpass:
            md.append(f"| `{c['tic']}` | `{c['lane']}` | {c['P_d']:.6f} | {c['T0_BJD']:.5f} | {c['depth_ppm']:.0f}+/-{c['err_ppm']:.0f} | {c['dur_h']:.2f} | {c['Rp_Re']:.2f} | {c['MES']:.1f} | {c['n_sec']} | {c['oe_sig']:.2f} | {c['cent_as']:.1f} | {c['ruwe']:.2f} | {c['Teq_K']:.0f} | {c['TSM']:.1f} | {c['FPP']:.4f} | {c['tier']} |")
    else:
        md.append("Zero triggers in this batch passed all 8 gates simultaneously; see audit table below.")
    if vall:
        md += ["", "### Vetted Trigger Audit (Top Signals & Veto Flags)", "",
               "| TIC | Lane | Period (d) | Depth (ppm) | Rp (Re) | MES | Passed | Vetoes |", "| :--- | :--- | ---: | ---: | ---: | ---: | :--- | :--- |"]
        for c in vall[:25]:
            md.append(f"| `{c['tic']}` | `{c['lane']}` | {c['P_d']:.5f} | {c['depth_ppm']:.0f} | {c['Rp_Re']:.2f} | {c['MES']:.1f} | `{c['passed']}` | `{c['vetoes']}` |")
    (out_dir / "DISCOVERIES.md").write_text("\n".join(md) + "\n")

class SelfTests(unittest.TestCase):
    def test_synthetic_usp_and_blind_check(self):
        rng = np.random.default_rng(7)
        t = np.concatenate([np.arange(1500.0, 1525.0, 2.0 / 1440.0), np.arange(3800.0, 3825.0, 2.0 / 1440.0)])
        sec = np.concatenate([np.full(len(t) // 2, 12, int), np.full(len(t) - len(t) // 2, 106, int)])
        f = 1.0 + 0.0015 * np.sin(2.0 * np.pi * t / 6.2) + rng.normal(0.0, 3e-4, len(t))
        p0, t0 = 0.31842, 1500.112; dur = exp_dur(p0, 0.42)
        f[np.abs(((t - t0 + 0.5 * p0) % p0) - 0.5 * p0) <= 0.5 * dur] -= (1.15 / (0.42 * RSUN_RE)) ** 2
        tc, fc, fec, sc = [], [], [], []
        for s in (12, 106):
            m = sec == s; a, b, c = detrend_lc(t[m], f[m], np.full(np.sum(m), 3e-4), np.zeros(np.sum(m), int), 0.9)
            tc.append(a); fc.append(b); fec.append(c); sc.append(np.full(len(a), s, int))
        t2, f2, fe2, s2 = np.concatenate(tc), np.concatenate(fc), np.concatenate(fec), np.concatenate(sc)
        pk = search_bls(t2, f2, fe2, 0.16, 0.50, 0.42)
        self.assertAlmostEqual(pk["p"], p0, delta=2e-5)
        rep = vet_signal(999, pk, t2, f2, fe2, s2, "USP", meta={"rstar": 0.42, "teff": 3550, "tmag": 11.2, "ruwe": 1.04}, tpf_off=2.1)
        self.assertTrue(rep["passed"], rep["vetoes"])
        b_ok = blind_check_sector(123, "S", p0, t0 + BTJD0, 620.0, dur * 24.0, 106, t2[s2 == 106], f2[s2 == 106])
        self.assertEqual(b_ok["verdict"], "CONFIRMED_ON_SCHEDULE")
