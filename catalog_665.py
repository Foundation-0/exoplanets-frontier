"""TOI-540b benchmark runner and re-export of purified Pre-Gaia-DR4 catalog builder."""
from __future__ import annotations
from pathlib import Path
import numpy as np
from bls_search import detrend_lc, search_bls
from converge_io import fetch_raw_spoc
from despot import harmonic_despot, plot_despot_card, scrambled_null_trial
from dr4_core import build_665_catalog


def run_benchmark(fig_d: Path):
    raw = fetch_raw_spoc(200322593, max_sec=4); tn, fn, fen, td, fd, fed, sd = [], [], [], [], [], [], []
    for l in sorted(raw, key=lambda x: x["sec"]):
        a, b, c = detrend_lc(l["t"], l["f"], l["fe"], l["q"], 1.0)
        if len(a) >= 80: tn.append(a); fn.append(b); fen.append(c)
        a2, b2, c2, _ = harmonic_despot(l["t"], l["f"], l["fe"], l["q"])
        if len(a2) >= 80: td.append(a2); fd.append(b2); fed.append(c2); sd.append(np.full(len(a2), l["sec"], int))
    pkn = search_bls(np.concatenate(tn), np.concatenate(fn), np.concatenate(fen), 0.25, 4.0, 0.19) if tn else None
    tc, fc, fec, sc = (np.concatenate(td), np.concatenate(fd), np.concatenate(fed), np.concatenate(sd)) if td else (np.array([]),) * 4
    pkd = search_bls(tc, fc, fec, 0.25, 4.0, 0.19) if len(tc) else None
    if pkd and tn:
        plot_despot_card(fig_d / "benchmark_TOI540_despot.png", 200322593, f"TOI-540b (Naive P={(pkn['p'] if pkn else -1.0):.4f}d -> Clean P={pkd['p']:.5f}d, MES={pkd['mes_c']:.1f})", tn[0], fn[0], tc, fc, pkd["p"], pkd["t0"], pkd["dur"])
    nm = [scrambled_null_trial(tc, fc, fec, sc, 0.19, 100 + k) for k in range(10)] if len(tc) else []
    return {"benchmark_tic": 200322593, "published_P_d": 1.239149, "naive_P_d": round(pkn["p"], 6) if pkn else -1.0,
            "naive_locked_on_spot": bool(pkn and abs(pkn["p"] - 1.239149) > 0.02), "despotted_P_d": round(pkd["p"], 6) if pkd else -1.0,
            "despotted_MES": round(pkd["mes_c"], 2) if pkd else 0.0, "recovered_within_0p1pct": bool(pkd and abs(pkd["p"] - 1.239149) / 1.239149 < 0.001),
            "null_trials": len(nm), "null_max_MES": max(nm) if nm else 0.0, "null_false_alarms_ge_9p5": sum(1 for m in nm if m >= 9.5)}


__all__ = ["build_665_catalog", "run_benchmark"]
