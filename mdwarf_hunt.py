#!/usr/bin/env python3
"""Fast 24-target search of quiet K/M dwarfs + TIC 31281821 TPF audit & S106-107 plots."""
from __future__ import annotations
import concurrent.futures as cf, json, math, time
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from astroquery.mast import Catalogs
from scipy import signal
from bls_search import BTJD0, detrend_lc, search_bls
from deep_search import ffi_cutout_lc, load_all_eumemic
from ledger import write_ledger
from mast_data import gaia_nb, index_mast_sectors, load_catalogs, stack_tic, tpf_diff
from vetter import vet_signal

def norm_sap(sap, sec):
    out = np.ones_like(sap, dtype=float)
    for s in np.unique(sec):
        m = (sec == s) & np.isfinite(sap) & (sap > 0)
        if np.sum(m) < 30: continue
        sm = sap[m] / float(np.median(sap[m])); k = min(int(np.sum(m)) - (1 - int(np.sum(m)) % 2), 301)
        out[m] = sm / np.maximum(signal.medfilt(sm, kernel_size=max(21, k | 1)), 0.2)
    return out

def spot_amp(t, f, p, t0, dur):
    ph = ((t - t0 + 0.5 * p) % p) - 0.5 * p; ot = np.abs(ph) > 1.3 * dur
    if np.sum(ot) < 60: return 0.0
    a = 2.0 * np.pi * (t[ot] - t0) / p
    X = np.column_stack([np.ones(np.sum(ot)), np.cos(a), np.sin(a), np.cos(2 * a), np.sin(2 * a)])
    c, *_ = np.linalg.lstsq(X, f[ot] - np.median(f[ot]), rcond=None); g = np.linspace(0, 2 * np.pi, 120)
    m = c[1] * np.cos(g) + c[2] * np.sin(g) + c[3] * np.cos(2 * g) + c[4] * np.sin(2 * g)
    return float(np.max(m) - np.min(m))

def save_card(out_p, tic, tag, p_d, t0_bt, dur_d, dep_ppm, t, f, sec, tpf=None):
    out_p.parent.mkdir(parents=True, exist_ok=True)
    ph = (((t - t0_bt + 0.5 * p_d) % p_d) - 0.5 * p_d) * 24.0; orb = np.round((t - t0_bt) / p_d).astype(int); fp = (f - np.median(f)) * 1e6
    nc = 4 if tpf else 3; fig, ax = plt.subplots(1, nc, figsize=(4.4 * nc, 3.8), dpi=135); fig.patch.set_facecolor("#0f172a")
    for a in ax: a.set_facecolor("#1e293b"); a.tick_params(colors="#cbd5e1", labelsize=8)
    wh = min(max(3.5 * dur_d * 24.0, 2.2), 0.45 * p_d * 24.0); zm = np.abs(ph) <= wh; bins = np.linspace(-wh, wh, 25); bi = np.digitize(ph[zm], bins) - 1
    ax[0].scatter(ph[zm], fp[zm], s=4, alpha=0.22, color="#94a3b8")
    pts = [(0.5 * (bins[b] + bins[b + 1]), float(np.mean(fp[zm][bi == b]))) for b in range(24) if np.sum(bi == b) >= 2]
    if pts: ax[0].plot([x for x, _ in pts], [y for _, y in pts], "o-", color="#38bdf8", ms=4, lw=1.2)
    ax[0].axhline(-dep_ppm, color="#f43f5e", ls="--", lw=1.1); ax[0].set_xlim(-wh, wh); ax[0].set_title(f"Phase Fold P={p_d:.5f}d", color="#f8fafc", fontsize=9)
    for par, col in [(1, "#34d399"), (0, "#fbbf24")]:
        pm = zm & ((orb % 2) == par); pbi = np.digitize(ph[pm], bins) - 1
        op = [(0.5 * (bins[b] + bins[b + 1]), float(np.mean(fp[pm][pbi == b]))) for b in range(24) if np.sum(pbi == b) >= 2]
        if op: ax[1].plot([x for x, _ in op], [y for _, y in op], "o-", color=col, ms=3.5, lw=1)
    ax[1].axhline(-dep_ppm, color="#f43f5e", ls="--", lw=1); ax[1].set_xlim(-wh, wh); ax[1].set_title("Odd vs Even", color="#f8fafc", fontsize=9)
    oph = ((t - t0_bt + 0.5 * p_d) % p_d) / p_d - 0.5; fb = np.linspace(-0.5, 0.5, 49); fbi = np.digitize(oph, fb) - 1
    fpt = [(0.5 * (fb[b] + fb[b + 1]), float(np.mean(fp[fbi == b]))) for b in range(48) if np.sum(fbi == b) >= 3]
    if fpt: ax[2].plot([x for x, _ in fpt], [y for _, y in fpt], "o-", color="#a78bfa", ms=3, lw=1)
    ax[2].axvline(0.5, color="#fb7185", ls=":"); ax[2].set_xlim(-0.5, 0.5); ax[2].set_title(f"Orbit ({len(np.unique(sec))} sec)", color="#f8fafc", fontsize=9)
    if tpf:
        ax[3].imshow(tpf["diff"], origin="lower", cmap="inferno")
        ax[3].plot(*tpf["sxy"], "+", color="#38bdf8", ms=10, mew=2); ax[3].plot(*tpf["dxy"], "x", color="#34d399", ms=8, mew=2)
        ax[3].set_title(f"S{tpf['sec']} TPF ({tpf['off']:.1f}\")", color="#f8fafc", fontsize=9)
    fig.suptitle(f"TIC {tic}: {tag}", color="#f8fafc", fontsize=10, fontweight="bold", y=1.01)
    fig.tight_layout(); fig.savefig(out_p, bbox_inches="tight", facecolor=fig.get_facecolor()); plt.close(fig)

def search_one(item):
    tid, rs0, ps = item["tic"], item["rstar"], item.get("p_seed")
    try:
        st, pr = stack_tic(tid, 0.65 if tid != 31281821 else 0.85)
        if not st or len(np.unique(st["sec"])) < 2: return None
        t, f, fe, sc, mx, my = st["t"], st["f"], st["fe"], st["sec"], st["mdx"], st["mdy"]; sap_n = norm_sap(st["sap"], sc)
        rs = float(st["hdr"]["rstar"]) if np.isfinite(st["hdr"]["rstar"]) and st["hdr"]["rstar"] > 0 else rs0
        pk = search_bls(t, f, fe, max(0.18, 0.85 * ps) if ps else 0.22, min(35.0, 1.2 * ps) if ps else 28.0, rs, seed_p=ps, fast_inv=False)
        if not pk or pk["mes_c"] < 8.8: return None
        gm, nbs = gaia_nb(tid, st["hdr"]["ra"], st["hdr"]["dec"], st["hdr"])
        lane = "MDWARF_T>13.5" if item["tmag"] > 13.5 else "S106_S107_KM_DWARF"
        rep = vet_signal(tid, pk, t, f, fe, sc, lane, sap_n, mx, my, st["crowd"], gm, nbs)
        amp = spot_amp(t, f, pk["p"], pk["t0"], pk["dur"])
        tpfi = tpf_diff(pr, rep["P_d"], rep["T0_BJD"] - BTJD0, rep["dur_h"] / 24.0) if (rep["passed"] or tid == 31281821) else None
        if tpfi: rep = vet_signal(tid, pk, t, f, fe, sc, lane, sap_n, mx, my, st["crowd"], gm, nbs, tpfi["off"])
        if amp > 0.52 * pk["dep"]:
            rep["passed"] = False; rep["tier"] = "rejected"
            rep["vetoes"] = ("" if rep["vetoes"] == "CLEAN" else rep["vetoes"] + ";") + f"spot_mod({amp/pk['dep']:.2f})"
        return rep, ({"t": t, "f": f, "sec": sc} if (rep["passed"] or tid == 31281821) else None), tpfi
    except Exception: return None

def main(out_dir: Path = Path("results"), max_new: int = 22):
    t0w = time.time(); out_dir.mkdir(parents=True, exist_ok=True); fig_d = out_dir / "figures"; fig_d.mkdir(exist_ok=True)
    s_max, new_s, smaps, _, _ = index_mast_sectors(); kt, _, _ = load_catalogs(); dfe = load_all_eumemic()
    for _, r in (dfe.iterrows() if not dfe.empty else []): kt.add(int(r["_tid"]))
    bres = pd.read_csv(out_dir / "blind_verification_s106_s107.csv").to_dict("records") if (out_dir / "blind_verification_s106_s107.csv").exists() else []
    pv = pd.read_csv(out_dir / "discovery_triggers_vetted.csv").to_dict("records") if (out_dir / "discovery_triggers_vetted.csv").exists() else []
    ps0 = json.loads((out_dir / "campaign_summary.json").read_text()) if (out_dir / "campaign_summary.json").exists() else {}
    already = {int(r["tic"]) for r in pv}

    for br in [x for x in bres if x.get("verdict") == "CONFIRMED_ON_SCHEDULE"][:4]:
        try:
            btid, bs = int(br["tic"]), int(br["sector"])
            tb0 = Catalogs.query_criteria(catalog="Tic", ID=[btid])
            er0 = dfe[dfe["_tid"] == btid].iloc[0]
            lc0 = ffi_cutout_lc(btid, float(tb0[0]["ra"]), float(tb0[0]["dec"]), bs)
            if lc0:
                tc, fc, _ = detrend_lc(lc0["t"], lc0["f"], lc0["fe"], lc0["q"], 1.5)
                t0v = float(er0.get("T0_BJD", er0.get("t0_bjd", 2459000.0))); t0_bt = t0v - BTJD0 if t0v > 2400000 else t0v
                save_card(fig_d / f"blind_S{bs}_TIC{btid}.png", btid, f"S{bs} Blind Confirmed (SNR={br['snr']:.1f}, Obs={br['obs_ppm']:.0f} ppm)", float(br["P_d"]), t0_bt, float(er0.get("T14_h", 2.4)) / 24.0, float(br["obs_ppm"]), tc, fc, np.full(len(tc), bs))
        except Exception: pass

    ov = [t for t in sorted(set(smaps[new_s[-2]]) & set(smaps[new_s[-1]])) if t not in kt and t not in already] if len(new_s) >= 2 else []
    sample = ov[::max(1, len(ov) // 750)][:750]; cands = []
    for i in range(0, len(sample), 150):
        try:
            tb = Catalogs.query_criteria(catalog="Tic", ID=sample[i : i + 150])
            for r in (tb or []):
                lc = str(r["lumclass"]) if "lumclass" in tb.colnames and r["lumclass"] is not None else "DWARF"
                teff, rad, tmag = (float(r[k]) if k in tb.colnames and r[k] is not None and np.isfinite(r[k]) else np.nan for k in ("Teff", "rad", "Tmag"))
                if lc.upper() in ("DWARF", "", "NONE") and 3450 <= teff <= 4700 and 0.35 <= rad <= 0.72 and 9.0 <= tmag <= 14.8:
                    cands.append((rad + 0.05 * max(0.0, tmag - 12.8), {"tic": int(r["ID"]), "rstar": rad, "tmag": tmag, "p_seed": None}))
        except Exception: pass
    cands.sort(key=lambda x: x[0])
    pool = [{"tic": 31281821, "rstar": 0.288, "tmag": 14.38, "p_seed": 1.0667132}] + [d for _, d in cands[:max_new]]

    vm = {int(r["tic"]): r for r in pv}; pm = {int(r["tic"]): r for r in pv if str(r.get("passed")).lower() == "true"}
    with cf.ThreadPoolExecutor(max_workers=5) as ex:
        for res in ex.map(search_one, pool):
            if not res: continue
            rep, ms, tpfi = res; tid = int(rep["tic"]); vm[tid] = rep
            if rep["passed"] and ms:
                pm[tid] = rep
                save_card(fig_d / f"discovery_TIC{tid}.png", tid, f"{rep['tier']} (MES={rep['MES']:.1f}, Rp={rep['Rp_Re']:.2f} Re, FPP={rep['FPP']:.4f})", rep["P_d"], rep["T0_BJD"] - BTJD0, rep["dur_h"] / 24.0, rep["depth_ppm"], ms["t"], ms["f"], ms["sec"], tpfi)
            elif tid == 31281821 and ms:
                pm.pop(tid, None)
                save_card(fig_d / f"audit_TIC{tid}.png", tid, f"Audit ({rep['vetoes']}, MES={rep['MES']:.1f})", rep["P_d"], rep["T0_BJD"] - BTJD0, rep["dur_h"] / 24.0, rep["depth_ppm"], ms["t"], ms["f"], ms["sec"], tpfi)

    vall = sorted(vm.values(), key=lambda x: float(x["MES"]), reverse=True); vpass = sorted(pm.values(), key=lambda x: float(x["MES"]), reverse=True)
    pd.DataFrame(vall).to_csv(out_dir / "discovery_triggers_vetted.csv", index=False)
    pd.DataFrame(vpass).to_csv(out_dir / "vetted_new_candidates.csv", index=False)
    vc = pd.DataFrame(bres)["verdict"].value_counts().to_dict() if bres else {}
    summ = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "sec_s": round(float(ps0.get("sec_s", 0.0)) + (time.time() - t0w), 1),
            "s_max": s_max, "new_sectors": new_s, "blind_pairs": len(bres), "blind_counts": vc,
            "unpromoted_tce_pool_total": int(ps0.get("unpromoted_tce_pool_total", 0)) + len(pool),
            "tce_stacked": int(ps0.get("tce_stacked", 0)) + len(pool), "usp_hosts": int(ps0.get("usp_hosts", 28)),
            "triggers_vetted": len(vall), "new_passed": len(vpass)}
    (out_dir / "campaign_summary.json").write_text(json.dumps(summ, indent=2))
    write_ledger(out_dir, summ, bres, vpass, vall); print(json.dumps(summ, indent=2))

if __name__ == "__main__":
    main()
