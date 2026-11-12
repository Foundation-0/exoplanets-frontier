#!/usr/bin/env python3
"""Unified Phase 3 Runner: 665-Candidate Gaia DR4 Catalog, Cycle 1-8 Ephemerides & Two-Pass Cool-Dwarf Hunt."""
from __future__ import annotations
import concurrent.futures as cf, json, time
from pathlib import Path
import numpy as np, pandas as pd
from blend_gaia import multisector_wcs_blend_audit
from bls_search import BTJD0
from campaigns import gv
from catalog_665 import build_665_catalog, run_benchmark
from converge_io import despot_one_mdwarf, fetch_spanning_lcs, plot_ttv_and_blend, q_tic
from deep_search import load_all_eumemic
from despot import plot_despot_card
from mast_data import index_mast_sectors, load_catalogs
from report_v3 import save_phase3_outputs
from ttv_track import extract_transit_times, fit_ephemeris_models, restack_and_siblings

TTV_TICS = [70557792, 166848769, 332644900, 206666344, 321011164, 170788495, 115241103, 7780510, 253923533, 115692561]
BLEND_TICS = [115692561, 253923533, 115241103, 166848769, 70557792, 224279805, 332644900, 31281821]
ROT_CHECK = [209462764, 168675553, 31281821]

def main(out_dir: Path = Path("results")):
    t0w = time.time(); out_dir.mkdir(parents=True, exist_ok=True); fig_d = out_dir / "figures"; fig_d.mkdir(exist_ok=True)
    s_max, new_s, smaps, _, _ = index_mast_sectors(); s_max = max(107, s_max); kt, _, _ = load_catalogs(); dfe = load_all_eumemic(); e_map = {}
    for _, r in (dfe.iterrows() if not dfe.empty else []): tid = int(r["_tid"]); kt.add(tid); e_map.setdefault(tid, r)
    bench = run_benchmark(fig_d); (out_dir / "benchmark_and_null.json").write_text(json.dumps(bench, indent=2))
    cat_rows, tmap = build_665_catalog(dfe)
    plat = [r for r in cat_rows if r["dr4_class"] in ("PLATINUM_SOLO_VERIFIED", "GOLD_LOW_BLEND")]
    pd.DataFrame(cat_rows).to_csv(out_dir / "eumemic_gaia_dr4_readiness.csv", index=False)
    pd.DataFrame(plat[:50]).to_csv(out_dir / "platinum_followup_targets.csv", index=False)
    for r in q_tic([t for t in set(TTV_TICS + BLEND_TICS) if t not in tmap]): tmap[int(r["ID"])] = r
    coords = {tid: (float(r["ra"]), float(r["dec"]), float(r["rad"]) if r["rad"] is not None and np.isfinite(r["rad"]) else 0.85) for tid, r in tmap.items() if r["ra"] is not None and np.isfinite(r["ra"])}
    ttv_r, ep_r, sib_r, snaps, blend_r, md_r, vpass = [], [], [], {}, [], [], []
    for tid in TTV_TICS:
        if tid not in e_map or tid not in coords: continue
        er, (ra, dec, rs) = e_map[tid], coords[tid]; p0, t0b = gv(er, "P_d", "period_d", "period"), gv(er, "T0_BJD", "t0_bjd", "T0_BTJD", "t0")
        dur_d, t0_bt = max(gv(er, "T14_h", "dur_h", d=2.2) / 24.0, 0.03), (t0b - BTJD0 if t0b > 2400000 else t0b)
        lcs, all_s = fetch_spanning_lcs(tid, ra, dec, win_d=max(1.3, min(2.4, dur_d * 8.0)))
        sr, ir = extract_transit_times(lcs, p0, t0_bt, dur_d, 240.0); ep_r.extend([{"tic": tid, **x} for x in sr])
        fit = fit_ephemeris_models(sr, ir, p0, t0_bt); rsib = restack_and_siblings(lcs, p0, t0_bt, fit["P_ref"], fit["T0_ref_BJD"] - BTJD0, dur_d, rs)
        ttv_r.append({"tic": tid, "upstream_tier": str(er.get("tier", er.get("Tier", "v3"))), "P_orig_d": round(p0, 6), "P_ref_d": fit["P_ref"], "T0_ref_BJD": fit["T0_ref_BJD"], "dP_sec": fit["dP_sec"], "n_sec_avail": len(all_s), "n_sec_used": len(lcs), "n_epochs_fit": fit["n_ep"], "rms0_min": fit["rms0_m"], "rms_lin_min": fit["rms_lin_m"], "dBIC_lin": fit["dBIC_lin"], "dBIC_quad": fit["dBIC_quad"], "dBIC_ttv": fit["dBIC_ttv"], "ttv_amp_min": fit["ttv_amp_m"], "snr_stale": rsib["snr_stale"], "snr_refined": rsib["snr_ref"], "class": fit["class"]})
        sb = rsib["sib"]; sib_r.append({"tic": tid, "host_P_ref_d": fit["P_ref"], "sib_P_d": round(sb["p"], 6) if sb else -1.0, "sib_MES": round(sb["mes_c"], 2) if sb else 0.0, "single_dip_snr": rsib["single_snr"], "single_dip_btjd": rsib["single_btjd"]})
        snaps[tid] = (sr, fit, lcs, dur_d, all_s)
    for tid in BLEND_TICS:
        if tid not in coords: continue
        ra, dec, _ = coords[tid]
        if tid in e_map:
            er = e_map[tid]; pu, t0b = gv(er, "P_d", "period_d", "period"), gv(er, "T0_BJD", "t0_bjd", "T0_BTJD", "t0")
            tu, dd, dp, ut = (t0b - BTJD0 if t0b > 2400000 else t0b), max(gv(er, "T14_h", "dur_h", d=2.0) / 24.0, 0.03), gv(er, "depth_ppm", "depth", d=500.0), str(er.get("tier", er.get("Tier", "v3")))
            dp = dp * 1e6 if dp < 1.0 else dp; shint = snaps[tid][4] if tid in snaps else None
            if tid in snaps: pu, tu = snaps[tid][1]["P_ref"], snaps[tid][1]["T0_ref_BJD"] - BTJD0
        else: pu, tu, dd, dp, ut, shint = 1.0667132, 2458439.143382 - BTJD0, 1.67 / 24.0, 12948.0, "MDWARF_ROT", [106, 107]
        br, bsnap = multisector_wcs_blend_audit(tid, ra, dec, pu, tu, dd, dp, shint); br["upstream_tier"], br["P_d"] = ut, round(pu, 6); blend_r.append(br)
        sr, fit, lcs, ddx = (snaps[tid][:4] if tid in snaps else ([], {"class": br["verdict"], "dP_sec": 0.0}, fetch_spanning_lcs(tid, ra, dec, 1.2)[0], dd))
        if lcs or bsnap: plot_ttv_and_blend(fig_d / f"converge_TIC{tid}.png", tid, sr, fit, np.concatenate([x["t"] for x in lcs]) if lcs else np.array([]), np.concatenate([x["f"] for x in lcs]) if lcs else np.array([]), pu, tu, ddx, bsnap)
    ov = [t for t in sorted(set(smaps[new_s[-2]]) & set(smaps[new_s[-1]])) if t not in kt and t not in set(ROT_CHECK)] if len(new_s) >= 2 else []
    fresh = [{"tic": int(r["ID"]), "rstar": float(r["rad"]), "lane": "TWO_PASS_COOL_DWARF"} for r in q_tic(ov[::max(1, len(ov) // 900)][:900]) if all(r[k] is not None and np.isfinite(r[k]) for k in ("Teff", "rad", "Tmag")) and 2950 <= float(r["Teff"]) <= 4200 and 0.14 <= float(r["rad"]) <= 0.58 and 9.0 <= float(r["Tmag"]) <= 14.6]
    fresh.sort(key=lambda x: x["rstar"])
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        for out in ex.map(despot_one_mdwarf, [{"tic": t, "rstar": 0.32, "lane": "TWO_PASS_ROTATOR_CHECK"} for t in ROT_CHECK] + fresh[:30]):
            if not out: continue
            rep, pl = out; md_r.append(rep)
            if rep.get("passed"):
                vpass.append(rep)
                if pl: plot_despot_card(fig_d / f"discovery_TIC{rep['tic']}.png", rep["tic"], f"NEW CANDIDATE ({rep['tier']}, P={rep['P_d']:.5f}d, Rp={rep['Rp_Re']:.2f}Re, MES={rep['MES']:.1f})", pl["raw_t"], pl["raw_f"], pl["t"], pl["f"], rep["P_d"], pl["t0"], pl["dur"])
            elif pl and rep["tic"] in set(ROT_CHECK):
                plot_despot_card(fig_d / f"despot_TIC{rep['tic']}.png", rep["tic"], f"Two-Pass ({rep['raw_ppm']:.0f}->{rep['clean_ppm']:.0f} ppm, MES={rep['MES']:.1f})", pl["raw_t"], pl["raw_f"], pl["t"], pl["f"], max(rep["P_d"], 0.25), pl["t0"], pl["dur"])
    save_phase3_outputs(out_dir, fig_d, t0w, s_max, bench, cat_rows, plat, ttv_r, ep_r, sib_r, blend_r, md_r, vpass)

if __name__ == "__main__":
    main()
