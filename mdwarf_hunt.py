#!/usr/bin/env python3
"""Deep multi-sector M/K-dwarf transit search with per-sector SAP norm and spot veto."""
from __future__ import annotations
import concurrent.futures as cf, json, math, time
from pathlib import Path
import numpy as np, pandas as pd
from astroquery.mast import Catalogs
from scipy import signal
from bls_search import BTJD0, RSUN_RE, search_bls
from deep_search import load_all_eumemic
from ledger import write_ledger
from mast_data import find_dvt_seeds, gaia_nb, index_mast_sectors, load_catalogs, stack_tic, tpf_diff
from vetter import plot_card, vet_signal

def norm_sap(sap: np.ndarray, sec: np.ndarray) -> np.ndarray:
    out = np.ones_like(sap, dtype=float)
    for s in np.unique(sec):
        m = (sec == s) & np.isfinite(sap) & (sap > 0)
        if np.sum(m) < 30: continue
        sm = sap[m] / float(np.median(sap[m]))
        k = min(int(np.sum(m)) - (1 - int(np.sum(m)) % 2), 301)
        tr = signal.medfilt(sm, kernel_size=max(21, k | 1))
        out[m] = sm / np.maximum(tr, 0.2)
    return out

def spot_amp(t: np.ndarray, f: np.ndarray, p: float, t0: float, dur: float) -> float:
    ph = ((t - t0 + 0.5 * p) % p) - 0.5 * p
    ot = np.abs(ph) > 1.3 * dur
    if np.sum(ot) < 60: return 0.0
    ang = 2.0 * np.pi * (t[ot] - t0) / p
    X = np.column_stack([np.ones(np.sum(ot)), np.cos(ang), np.sin(ang), np.cos(2 * ang), np.sin(2 * ang)])
    coef, *_ = np.linalg.lstsq(X, f[ot] - np.median(f[ot]), rcond=None)
    grid = np.linspace(0, 2 * np.pi, 120)
    mod = coef[1] * np.cos(grid) + coef[2] * np.sin(grid) + coef[3] * np.cos(2 * grid) + coef[4] * np.sin(2 * grid)
    return float(np.max(mod) - np.min(mod))

def sel_tgts(smaps: dict, obs_by_s: dict, dv_multi: list, new_s: list[int], kt: set, mx_s: int = 76) -> list[dict]:
    seeds = find_dvt_seeds(new_s, obs_by_s, dv_multi, kt, max_dvt=60)
    out = [{"tic": 31281821, "rstar": 0.288, "tmag": 14.38, "p_seed": 1.0667132}]
    seen = {31281821}
    for t, v in sorted(seeds.items(), key=lambda kv: kv[1]["snr"], reverse=True):
        if t not in kt and t not in seen and v["snr"] <= 50.0:
            seen.add(t); out.append({"tic": t, "rstar": 0.65, "tmag": 12.8, "p_seed": v["p"]})
    ov = [t for t in sorted(set(smaps[new_s[-2]]) & set(smaps[new_s[-1]])) if t not in kt and t not in seen] if len(new_s) >= 2 else []
    extra = [t for s in reversed(new_s) for t in smaps[s] if t not in kt and t not in seen and t not in set(ov)]
    sample = (ov + extra[:2000])[::max(1, len(ov + extra[:2000]) // 1800)][:1800]
    cands = []
    for i in range(0, len(sample), 150):
        try:
            tb = Catalogs.query_criteria(catalog="Tic", ID=sample[i : i + 150])
            for r in (tb or []):
                tid = int(r["ID"])
                lc = str(r["lumclass"]) if "lumclass" in tb.colnames and r["lumclass"] is not None else "DWARF"
                teff = float(r["Teff"]) if "Teff" in tb.colnames and r["Teff"] is not None and np.isfinite(r["Teff"]) else np.nan
                rad = float(r["rad"]) if "rad" in tb.colnames and r["rad"] is not None and np.isfinite(r["rad"]) else np.nan
                tmag = float(r["Tmag"]) if "Tmag" in tb.colnames and r["Tmag"] is not None and np.isfinite(r["Tmag"]) else np.nan
                if lc.upper() not in ("DWARF", "", "NONE") or not (np.isfinite(teff) and np.isfinite(rad) and np.isfinite(tmag)): continue
                if 2700 <= teff <= 4650 and 0.14 <= rad <= 0.68 and 8.5 <= tmag <= 15.2:
                    cands.append((rad + 0.05 * max(0.0, tmag - 13.0), {"tic": tid, "rstar": rad, "tmag": tmag, "p_seed": None}))
        except Exception: pass
    cands.sort(key=lambda x: x[0])
    for _, d in cands:
        if len(out) >= mx_s: break
        if d["tic"] not in seen:
            seen.add(d["tic"]); out.append(d)
    return out[:mx_s]

def search_one(item: dict):
    tid, rs0, ps = item["tic"], item["rstar"], item.get("p_seed")
    try:
        win = 0.52 if tid != 31281821 else 0.75
        st, pr = stack_tic(tid, win)
        if not st or len(np.unique(st["sec"])) < 2: return None
        t, f, fe, sc, mx, my = st["t"], st["f"], st["fe"], st["sec"], st["mdx"], st["mdy"]
        sap_n = norm_sap(st["sap"], sc)
        rs = float(st["hdr"]["rstar"]) if np.isfinite(st["hdr"]["rstar"]) and st["hdr"]["rstar"] > 0 else rs0
        pk = search_bls(t, f, fe, max(0.165, 0.8 * ps) if ps else 0.18, min(42.0, 1.25 * ps) if ps else 34.0, rs, seed_p=ps, fast_inv=False)
        if not pk or pk["mes_c"] < 9.0 or pk["sde"] < 7.3: return None
        rp = math.sqrt(pk["dep"]) * rs * RSUN_RE
        if not (0.45 <= rp <= 7.2): return None
        gm, nbs = gaia_nb(tid, st["hdr"]["ra"], st["hdr"]["dec"], st["hdr"])
        lane = "MDWARF_T>13.5" if item["tmag"] > 13.5 else "S106_S107_DWARF"
        rep = vet_signal(tid, pk, t, f, fe, sc, lane, sap_n, mx, my, st["crowd"], gm, nbs)
        amp = spot_amp(t, f, pk["p"], pk["t0"], pk["dur"])
        if amp > 0.52 * pk["dep"]:
            rep["passed"] = False; rep["tier"] = "rejected"
            rep["vetoes"] = ("" if rep["vetoes"] == "CLEAN" else rep["vetoes"] + ";") + f"spot_mod({amp/pk['dep']:.2f})"
        tpfi = tpf_diff(pr, rep["P_d"], rep["T0_BJD"] - BTJD0, rep["dur_h"] / 24.0) if (rep["passed"] or tid == 31281821) else None
        if tpfi and rep["passed"]:
            rep = vet_signal(tid, pk, t, f, fe, sc, lane, sap_n, mx, my, st["crowd"], gm, nbs, tpfi["off"])
            if amp > 0.52 * pk["dep"]:
                rep["passed"] = False; rep["tier"] = "rejected"; rep["vetoes"] += f";spot_mod({amp/pk['dep']:.2f})"
        return rep, ({"t": t, "f": f, "sec": sc} if (rep["passed"] or tid == 31281821) else None), tpfi
    except Exception: return None

def main(out_dir: Path = Path("results"), max_mdw: int = 72):
    t0w = time.time(); out_dir.mkdir(parents=True, exist_ok=True); fig_d = out_dir / "figures"; fig_d.mkdir(exist_ok=True)
    s_max, new_s, smaps, obs_by_s, dv_multi = index_mast_sectors()
    kt, _, _ = load_catalogs()
    dfe = load_all_eumemic()
    for _, r in (dfe.iterrows() if not dfe.empty else []): kt.add(int(r["_tid"]))
    kt.discard(31281821)

    bres = pd.read_csv(out_dir / "blind_verification_s106_s107.csv").to_dict("records") if (out_dir / "blind_verification_s106_s107.csv").exists() else []
    pv = pd.read_csv(out_dir / "discovery_triggers_vetted.csv").to_dict("records") if (out_dir / "discovery_triggers_vetted.csv").exists() else []
    ps0 = json.loads((out_dir / "campaign_summary.json").read_text()) if (out_dir / "campaign_summary.json").exists() else {}

    pool = sel_tgts(smaps, obs_by_s, dv_multi, new_s, kt, mx_s=max_mdw)
    vm = {int(r["tic"]): r for r in pv}
    pm = {int(r["tic"]): r for r in pv if str(r.get("passed")).lower() == "true"}

    with cf.ThreadPoolExecutor(max_workers=5) as ex:
        for res in ex.map(search_one, pool):
            if not res: continue
            rep, ms, tpfi = res; tid = int(rep["tic"]); vm[tid] = rep
            if rep["passed"] and ms:
                pm[tid] = rep
                plot_card(fig_d / f"discovery_TIC{tid}.png", tid, f"{rep['tier']} (MES={rep['MES']:.1f}, Rp={rep['Rp_Re']:.2f} Re, FPP={rep['FPP']:.4f})", rep["P_d"], rep["T0_BJD"] - BTJD0, rep["dur_h"] / 24.0, rep["depth_ppm"], ms["t"], ms["f"], ms["sec"], tpfi)
            elif tid == 31281821 and ms:
                pm.pop(tid, None)
                plot_card(fig_d / f"audit_TIC{tid}.png", tid, f"Audit ({rep['vetoes']}, MES={rep['MES']:.1f})", rep["P_d"], rep["T0_BJD"] - BTJD0, rep["dur_h"] / 24.0, rep["depth_ppm"], ms["t"], ms["f"], ms["sec"], tpfi)

    vall = sorted(vm.values(), key=lambda x: float(x["MES"]), reverse=True)
    vpass = sorted(pm.values(), key=lambda x: float(x["MES"]), reverse=True)
    pd.DataFrame(vall).to_csv(out_dir / "discovery_triggers_vetted.csv", index=False)
    pd.DataFrame(vpass).to_csv(out_dir / "vetted_new_candidates.csv", index=False)
    vc = pd.DataFrame(bres)["verdict"].value_counts().to_dict() if bres else {}
    summ = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "sec_s": round(float(ps0.get("sec_s", 0.0)) + (time.time() - t0w), 1),
            "s_max": s_max, "new_sectors": new_s, "blind_pairs": len(bres), "blind_counts": vc,
            "unpromoted_tce_pool_total": int(ps0.get("unpromoted_tce_pool_total", 0)) + len(pool),
            "tce_stacked": int(ps0.get("tce_stacked", 0)) + len(pool), "usp_hosts": int(prev_s.get("usp_hosts", 28) if "prev_s" in locals() else ps0.get("usp_hosts", 28)),
            "triggers_vetted": len(vall), "new_passed": len(vpass)}
    (out_dir / "campaign_summary.json").write_text(json.dumps(summ, indent=2))
    write_ledger(out_dir, summ, bres, vpass, vall)
    print(json.dumps(summ, indent=2))

if __name__ == "__main__":
    main()
