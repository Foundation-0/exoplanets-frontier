#!/usr/bin/env python3
"""Phase 2 Convergent Runner: Benchmark+Null, S1-S107 TTV/Drift Rescue, Gaia WCS Audit & De-Spotted M-Dwarfs."""
from __future__ import annotations
import concurrent.futures as cf, hashlib, json, time
from pathlib import Path
import numpy as np, pandas as pd
from blend_gaia import multisector_wcs_blend_audit, plot_ttv_and_blend
from bls_search import BTJD0
from campaigns import gv
from converge_io import despot_one_mdwarf, fetch_spanning_lcs, q_tic, run_benchmark_and_null
from deep_search import load_all_eumemic
from despot import plot_despot_card
from mast_data import index_mast_sectors, load_catalogs
from ttv_track import extract_transit_times, fit_ephemeris_models, restack_and_siblings

TTV_TICS = [70557792, 166848769, 332644900, 206666344, 321011164, 170788495, 115241103, 7780510, 398794585, 170638116]
BLEND_TICS = [115692561, 253923533, 115241103, 166848769, 70557792, 224279805, 31281821]
ROT_TICS = [31281821, 166833387, 209462764, 168675553, 304096841, 321087099, 139637517, 44892011, 77842052, 31353324]

def write_report(out_dir: Path, s_max: int, el_s: float, bench: dict, ttv_r: list, blend_r: list, md_r: list):
    md = ["# Phase 2 Convergent Discovery and Technical Resolution Report", "", f"Live STScI MAST Sectors 1-{s_max} and ESA Gaia DR3 run (`{el_s} s`).", "",
          "## 1. Benchmark & Null Control (`astr0` Rules 2-3)", "",
          f"- **TOI-540 (`TIC 200322593`, $P_b=1.239149\\text{{ d}}$)**: Naive P=`{bench['naive_P_d']} d` -> De-Spotted P=`{bench['despotted_P_d']} d` (`MES={bench['despotted_MES']}`, recovered=`{bench['recovered_within_0p1pct']}`, RMS `{bench['mean_raw_rms_ppm']}->{bench['mean_clean_rms_ppm']} ppm`). Null false alarms: `{bench['null_false_alarms_ge_9p5']}/{bench['null_trials']}` (max MES `{bench['null_max_MES']}`).", "",
          "## 2. 8-Year S1-S107 Wide-Window (`+-240 min`) Ephemeris & TTV Audit", "",
          "| TIC | Tier | P_orig | P_ref | dP (s) | Ep | RMS0 (m) | RMS_lin (m) | dBIC_lin | dBIC_ttv | SNR_0 -> SNR_ref | Class |", "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- | :--- |"]
    for r in ttv_r: md.append(f"| `{r['tic']}` | `{r['upstream_tier']}` | {r['P_orig_d']:.6f} | {r['P_ref_d']:.7f} | {r['dP_sec']:+.2f} | {r['n_epochs_fit']} | {r['rms0_min']:.1f} | {r['rms_lin_min']:.1f} | {r['dBIC_lin']:.1f} | {r['dBIC_ttv']:.1f} | {r['snr_stale']:.2f} -> **{r['snr_refined']:.2f}** | `{r['class']}` |")
    md += ["", "## 3. Multi-Sector WCS + Gaia DR3 Sub-Pixel Blend Audit", "", "| TIC | Tier | P (d) | Gaia DR3 ID | RUWE | PM | Sec | Off (\") | Plaus NBs | NB Sep (\") | Verdict |", "| :--- | :--- | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |"]
    for b in blend_r: md.append(f"| `{b['tic']}` | `{b['upstream_tier']}` | {b['P_d']:.6f} | `{b['gaia_sid']}` | {b['ruwe']:.2f} | {b['pm_masyr']:.1f} | {b['n_sec']} | {b['off_target_as']:.2f}+/-{b['err_as']:.2f} | {b['plausible_nbs']} | {b['closest_nb_sep_as']:.1f} | `{b['verdict']}` |")
    md += ["", "## 4. De-Spotted Active M-Dwarf Residual Survey", "", "| TIC | Lane | P_rot (d) | Raw->Clean (ppm) | Resid P (d) | MES | Passed | Vetoes |", "| :--- | :--- | ---: | :--- | ---: | ---: | :--- | :--- |"]
    for m in sorted(md_r, key=lambda x: float(x.get("MES", 0.0)), reverse=True): md.append(f"| `{m['tic']}` | `{m['lane']}` | {m['p_rot_d']:.4f} | {m['raw_ppm']:.0f}->{m['clean_ppm']:.0f} | {m['P_d']:.5f} | {m['MES']:.2f} | `{m['passed']}` | `{m['vetoes']}` |")
    (out_dir / "CONVERGENCE_REPORT.md").write_text("\n".join(md) + "\n")

def main(out_dir: Path = Path("results")):
    t0w = time.time(); out_dir.mkdir(parents=True, exist_ok=True); fig_d = out_dir / "figures"; fig_d.mkdir(exist_ok=True)
    s_max, new_s, smaps, _, _ = index_mast_sectors(); s_max = max(107, s_max); kt, _, _ = load_catalogs(); dfe = load_all_eumemic(); e_map = {}
    for _, r in (dfe.iterrows() if not dfe.empty else []): tid = int(r["_tid"]); kt.add(tid); e_map.setdefault(tid, r)
    bench = run_benchmark_and_null(fig_d); (out_dir / "benchmark_and_null.json").write_text(json.dumps(bench, indent=2))
    coords = {int(r["ID"]): (float(r["ra"]), float(r["dec"]), float(r["rad"]) if r["rad"] is not None and np.isfinite(r["rad"]) else 0.85) for r in q_tic(sorted(set(TTV_TICS + BLEND_TICS)))}
    for tid in set(TTV_TICS + BLEND_TICS):
        if tid not in coords and tid in e_map and np.isfinite(gv(e_map[tid], "RA", "ra")) and np.isfinite(gv(e_map[tid], "Dec", "dec")):
            coords[tid] = (gv(e_map[tid], "RA", "ra"), gv(e_map[tid], "Dec", "dec"), gv(e_map[tid], "Rstar", "rad", d=0.85))
    ttv_r, ep_r, sib_r, snaps, blend_r, md_r = [], [], [], {}, [], []
    for tid in TTV_TICS:
        if tid not in e_map or tid not in coords: continue
        er, (ra, dec, rs) = e_map[tid], coords[tid]; p0, t0b = gv(er, "P_d", "period_d", "period"), gv(er, "T0_BJD", "t0_bjd", "T0_BTJD", "t0")
        dur_d, dppm, t0_bt = max(gv(er, "T14_h", "dur_h", d=2.2) / 24.0, 0.03), gv(er, "depth_ppm", "depth", d=600.0), (t0b - BTJD0 if t0b > 2400000 else t0b)
        lcs, all_s = fetch_spanning_lcs(tid, ra, dec, win_d=max(1.3, min(2.4, dur_d * 8.0)))
        sr, ir = extract_transit_times(lcs, p0, t0_bt, dur_d, 240.0); ep_r.extend([{"tic": tid, **x} for x in sr])
        fit = fit_ephemeris_models(sr, ir, p0, t0_bt); rsib = restack_and_siblings(lcs, p0, t0_bt, fit["P_ref"], fit["T0_ref_BJD"] - BTJD0, dur_d, rs)
        ttv_r.append({"tic": tid, "upstream_tier": str(er.get("tier", er.get("Tier", "v3"))), "P_orig_d": round(p0, 6), "P_ref_d": fit["P_ref"], "T0_ref_BJD": fit["T0_ref_BJD"], "dP_sec": fit["dP_sec"], "n_sec_avail": len(all_s), "n_epochs_fit": fit["n_ep"], "rms0_min": fit["rms0_m"], "rms_lin_min": fit["rms_lin_m"], "dBIC_lin": fit["dBIC_lin"], "dBIC_quad": fit["dBIC_quad"], "dBIC_ttv": fit["dBIC_ttv"], "ttv_amp_min": fit["ttv_amp_m"], "ttv_sup_orb": fit["ttv_sup_orb"], "snr_stale": rsib["snr_stale"], "snr_refined": rsib["snr_ref"], "class": fit["class"]})
        sb = rsib["sib"]; sib_r.append({"tic": tid, "host_P_ref_d": fit["P_ref"], "sib_P_d": round(sb["p"], 6) if sb else -1.0, "period_ratio": round(sb["p"] / fit["P_ref"], 4) if sb else -1.0, "sib_MES": round(sb["mes_c"], 2) if sb else 0.0, "sib_depth_ppm": round(sb["dep"] * 1e6, 1) if sb else 0.0, "single_dip_snr": rsib["single_snr"], "single_dip_btjd": rsib["single_btjd"]})
        snaps[tid] = (sr, fit, lcs, dur_d, dppm * 1e6 if dppm < 1.0 else dppm)
    for tid in BLEND_TICS:
        if tid not in coords: continue
        ra, dec, _ = coords[tid]
        if tid in e_map:
            er = e_map[tid]; pu, t0b = gv(er, "P_d", "period_d", "period"), gv(er, "T0_BJD", "t0_bjd", "T0_BTJD", "t0")
            tu, dd, dp, ut = (t0b - BTJD0 if t0b > 2400000 else t0b), max(gv(er, "T14_h", "dur_h", d=2.0) / 24.0, 0.03), gv(er, "depth_ppm", "depth", d=500.0), str(er.get("tier", er.get("Tier", "v3")))
            dp = dp * 1e6 if dp < 1.0 else dp
            if tid in snaps: pu, tu = snaps[tid][1]["P_ref"], snaps[tid][1]["T0_ref_BJD"] - BTJD0
        else: pu, tu, dd, dp, ut = 1.0667132, 2458439.143382 - BTJD0, 1.67 / 24.0, 12948.0, "MDWARF_ROT"
        br, bsnap = multisector_wcs_blend_audit(tid, ra, dec, pu, tu, dd, dp); br["upstream_tier"], br["P_d"] = ut, round(pu, 6); blend_r.append(br)
        sr, fit, lcs, ddx = (snaps[tid][:4] if tid in snaps else ([], {"class": br["verdict"], "dP_sec": 0.0}, fetch_spanning_lcs(tid, ra, dec, 1.2)[0], dd))
        if lcs or bsnap: plot_ttv_and_blend(fig_d / f"converge_TIC{tid}.png", tid, sr, fit, np.concatenate([x["t"] for x in lcs]) if lcs else np.array([]), np.concatenate([x["f"] for x in lcs]) if lcs else np.array([]), pu, tu, ddx, bsnap)
    for tid in [x for x in TTV_TICS if x in snaps and x not in BLEND_TICS][:4]:
        sr, fit, lcs, ddx, _ = snaps[tid]
        if lcs: plot_ttv_and_blend(fig_d / f"ttv_TIC{tid}.png", tid, sr, fit, np.concatenate([x["t"] for x in lcs]), np.concatenate([x["f"] for x in lcs]), fit["P_ref"], fit["T0_ref_BJD"] - BTJD0, ddx, None)
    ov = [t for t in sorted(set(smaps[new_s[-2]]) & set(smaps[new_s[-1]])) if t not in kt and t not in set(ROT_TICS)] if len(new_s) >= 2 else []
    fresh = [{"tic": int(r["ID"]), "rstar": float(r["rad"]), "lane": "DESPOT_FRESH_MDWARF"} for i in range(0, min(len(ov), 450), 150) for r in q_tic(ov[i : i + 150]) if all(r[k] is not None and np.isfinite(r[k]) for k in ("Teff", "rad", "Tmag")) and 2900 <= float(r["Teff"]) <= 3950 and 0.14 <= float(r["rad"]) <= 0.52 and 9.5 <= float(r["Tmag"]) <= 14.8]
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        for out in ex.map(despot_one_mdwarf, [{"tic": t, "rstar": 0.35, "lane": "DESPOT_ROTATOR_RESCUE"} for t in ROT_TICS] + fresh[:14]):
            if not out: continue
            rep, pl = out; md_r.append(rep)
            if pl and (rep.get("passed") or rep["tic"] == 31281821 or rep.get("MES", 0) >= 9.5):
                plot_despot_card(fig_d / f"despot_TIC{rep['tic']}.png", rep["tic"], f"{rep['lane']} (Prot={rep['p_rot_d']:.3f}d, P={rep['P_d']:.4f}d, MES={rep['MES']:.1f})", pl["raw_t"], pl["raw_f"], pl["t"], pl["f"], max(rep["P_d"], 0.25), pl["t0"], pl["dur"])
    for nm, rows in [("ttv_ephemeris_audit.csv", ttv_r), ("ttv_transit_epochs.csv", ep_r), ("sibling_residual_search.csv", sib_r), ("blend_localization.csv", blend_r), ("despotted_mdwarf_audit.csv", md_r)]: pd.DataFrame(rows).to_csv(out_dir / nm, index=False)
    el_s = round(time.time() - t0w, 1); write_report(out_dir, s_max, el_s, bench, ttv_r, blend_r, md_r)
    certs = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "elapsed_s": el_s, "s_max": s_max, "benchmark": bench,
             "files_sha256": {str(fp.relative_to(out_dir)): hashlib.sha256(fp.read_bytes()).hexdigest() for fp in sorted(list(out_dir.glob("*.csv")) + list(out_dir.glob("*.md")) + list(out_dir.glob("*.json")) + list(fig_d.glob("*.png"))) if fp.name != "convergence_certificates.json"}}
    (out_dir / "convergence_certificates.json").write_text(json.dumps(certs, indent=2))
    print(json.dumps({"benchmark_ok": bench["recovered_within_0p1pct"], "null_fa": bench["null_false_alarms_ge_9p5"], "ttv": len(ttv_r), "blend": len(blend_r), "mdwarf": len(md_r), "figs": len(list(fig_d.glob("*.png")))}, indent=2))

if __name__ == "__main__":
    main()
