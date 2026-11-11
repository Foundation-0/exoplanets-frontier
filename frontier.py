#!/usr/bin/env python3
"""Campaign runner for TESS Sectors 1-107 verification and discovery."""
from __future__ import annotations
import argparse, concurrent.futures as cf, json, time, unittest
from pathlib import Path
import numpy as np, pandas as pd
from bls_search import BTJD0, detrend_lc, search_bls
from ledger import SelfTests, write_ledger
from vetter import plot_card

def run_live(out_dir: Path, max_tce: int = 36, max_usp: int = 28):
    from campaigns import run_blind, vet_one
    from mast_data import dl_lc, find_dvt_seeds, index_mast_sectors, load_catalogs
    out_dir.mkdir(parents=True, exist_ok=True); fig_d = out_dir / "figures"; fig_d.mkdir(exist_ok=True); t0w = time.time()
    s_max, new_s, smaps, obs_by_s, dv_multi = index_mast_sectors()
    kt, kp, dfe = load_catalogs()
    print(f"[info] s_max={s_max} new_s={new_s} e_rows={len(dfe)}")

    bres, bhosts, e_by_tic = run_blind(dfe, new_s, smaps, kt, kp, fig_d)
    pd.DataFrame(bres).to_csv(out_dir / "blind_verification_s106_s107.csv", index=False)

    tce_by_tic = find_dvt_seeds(new_s, obs_by_s, dv_multi, kt, max_dvt=max(44, max_tce + 10))
    if len(tce_by_tic) < 20 and len(new_s) >= 2:
        ov = [t for t in sorted(set(smaps[new_s[-2]]) & set(smaps[new_s[-1]])) if t not in kt and t not in tce_by_tic][:80]
        def _scan_ov(tid):
            l1, l2 = dl_lc(smaps[new_s[-2]][tid], new_s[-2]), dl_lc(smaps[new_s[-1]][tid], new_s[-1])
            if not (l1 and l2) or (np.isfinite(l2["teff"]) and l2["teff"] > 5400) or (np.isfinite(l2["rstar"]) and l2["rstar"] > 1.05): return None
            t1, f1, e1 = detrend_lc(l1["t"], l1["f"], l1["fe"], l1["q"], 1.1); t2, f2, e2 = detrend_lc(l2["t"], l2["f"], l2["fe"], l2["q"], 1.1)
            if len(t1) < 80 or len(t2) < 80: return None
            pk = search_bls(np.concatenate([t1, t2]), np.concatenate([f1, f2]), np.concatenate([e1, e2]), 0.18, 16.0, l2["rstar"] if np.isfinite(l2["rstar"]) else 0.65, fast_inv=True)
            return {"tic": tid, "p": pk["p"], "snr": pk["mes_c"]} if (pk and pk["mes_c"] >= 7.5) else None
        with cf.ThreadPoolExecutor(max_workers=4) as ex:
            for x in ex.map(_scan_ov, ov):
                if x and x["tic"] not in tce_by_tic: tce_by_tic[x["tic"]] = x

    vall, vpass = [], []
    tce_list = sorted(tce_by_tic.values(), key=lambda x: x["snr"], reverse=True)[:max_tce]
    usp_hosts = sorted(bhosts.keys())[:max_usp]
    for pool, lane, pfx in [(tce_list, "S106_S107_STACK", "discovery"), (usp_hosts, "USP_SIBLING", "usp")]:
        with cf.ThreadPoolExecutor(max_workers=4) as ex:
            for res in ex.map(lambda it: vet_one(it, lane, e_by_tic, kp), pool):
                if not res: continue
                rep, ms, tpfi = res; vall.append(rep)
                if rep["passed"] and ms:
                    vpass.append(rep)
                    plot_card(fig_d / f"{pfx}_TIC{rep['tic']}.png", rep["tic"], f"{rep['tier']} (MES={rep['MES']:.1f}, Rp={rep['Rp_Re']:.2f} Re, FPP={rep['FPP']:.4f})", rep["P_d"], rep["T0_BJD"] - BTJD0, rep["dur_h"] / 24.0, rep["depth_ppm"], ms["t"], ms["f"], ms["sec"], tpfi)

    vall.sort(key=lambda x: x["MES"], reverse=True); vpass.sort(key=lambda x: x["MES"], reverse=True)
    pd.DataFrame(vall).to_csv(out_dir / "discovery_triggers_vetted.csv", index=False)
    pd.DataFrame(vpass).to_csv(out_dir / "vetted_new_candidates.csv", index=False)
    vc = pd.DataFrame(bres)["verdict"].value_counts().to_dict() if bres else {}
    summ = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "sec_s": round(time.time() - t0w, 1), "s_max": s_max, "new_sectors": new_s,
            "blind_pairs": len(bres), "blind_counts": vc, "unpromoted_tce_pool_total": len(tce_by_tic), "tce_stacked": len(tce_list),
            "usp_hosts": len(usp_hosts), "triggers_vetted": len(vall), "new_passed": len(vpass)}
    (out_dir / "campaign_summary.json").write_text(json.dumps(summ, indent=2))
    write_ledger(out_dir, summ, bres, vpass, vall)
    return summ

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out-dir", type=Path, default=Path("results"))
    ap.add_argument("--max-tce", "--max-tce-targets", dest="max_tce", type=int, default=36)
    ap.add_argument("--max-usp", "--max-usp-targets", dest="max_usp", type=int, default=28)
    args = ap.parse_args()
    if args.self_test: unittest.main(argv=[__file__], verbosity=2)
    else: print(json.dumps(run_live(args.out_dir, args.max_tce, args.max_usp), indent=2))
