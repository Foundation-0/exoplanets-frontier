"""Blind Sector 106-107 check and multi-year stacking worker routines."""
from __future__ import annotations
import concurrent.futures as cf, re
import numpy as np
from bls_search import BTJD0, detrend_lc, search_bls
from mast_data import dl_lc, gaia_nb, stack_tic, tpf_diff
from vetter import blind_check_sector, plot_card, vet_signal

def gv(row, *ks, d=np.nan):
    cl = {str(k).strip().lower(): k for k in row.index}
    for k in ks:
        c = cl.get(k.lower())
        if c is not None and row[c] is not None and str(row[c]) != "nan":
            s = re.sub(r"[^\d.+-eE]", "", str(row[c])) if k.lower() == "tic" else row[c]
            try: return float(s)
            except Exception: pass
    return float(d)

def run_blind(dfe, new_s, smaps, kt, kp, fig_d):
    btasks, e_by_tic = [], {}
    for _, row in (dfe.iterrows() if not dfe.empty else []):
        tv, pv = gv(row, "TIC", "tic", "tic_id"), gv(row, "P_d", "period_d", "period")
        if not np.isfinite(tv): continue
        tid = int(tv); kt.add(tid)
        if np.isfinite(pv) and pv > 0: kp.setdefault(tid, []).append(pv)
        tr = str(row.get("tier", row.get("Tier", row.get("src", "earlier"))) or "earlier").strip()
        if any(x in tr.lower() for x in ("rej", "false", "eb", "l2")): continue
        e_by_tic.setdefault(tid, []).append(row)
        for s in new_s:
            if tid in smaps[s]: btasks.append((tid, tr, row, s, smaps[s][tid]))
    bres, bhosts, bplots = [], {}, []
    with cf.ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(dl_lc, u, s): (tid, tr, row, s) for tid, tr, row, s, u in btasks}
        for fut in cf.as_completed(futs):
            tid, tr, row, s = futs[fut]; lc = fut.result()
            if not lc: continue
            bhosts[tid] = row
            pd_, t0b = gv(row, "P_d", "period_d", "period"), gv(row, "T0_BJD", "t0_bjd", "T0_BTJD", "t0_btjd", "t0")
            dppm, durh = gv(row, "depth_ppm", "depth", d=500.0), gv(row, "T14_h", "dur_h", "duration_h", d=2.2)
            if not (np.isfinite(pd_) and np.isfinite(t0b)): continue
            if dppm < 1.0: dppm *= 1e6
            tc, fc, _ = detrend_lc(lc["t"], lc["f"], lc["fe"], lc["q"], win_d=max(1.2, min(2.4, (durh / 24.0) * 8.0)))
            im = np.clip(np.searchsorted(lc["t"], tc), 0, len(lc["t"]) - 1) if len(tc) else np.array([], int)
            rb = blind_check_sector(tid, tr, pd_, t0b, dppm, durh, s, tc, fc, lc["mdx"][im] if len(tc) else None, lc["mdy"][im] if len(tc) else None)
            bres.append(rb)
            if rb["verdict"] in ("CONFIRMED_ON_SCHEDULE", "TIMING_OFFSET_TTV") and len(bplots) < 12:
                t0_bt = t0b - BTJD0 if t0b > 2400000 else t0b
                bplots.append((fig_d / f"blind_S{s}_TIC{tid}.png", tid, f"S{s} Blind Confirmed (SNR={rb['snr']:.1f}, Obs={rb['obs_ppm']:.0f} ppm)", pd_, t0_bt, durh / 24.0, rb["obs_ppm"], tc, fc, np.full(len(tc), s)))
    for ap in bplots: plot_card(*ap)
    bres.sort(key=lambda r: r["snr"], reverse=True)
    return bres, bhosts, e_by_tic

def vet_one(item, lane, e_by_tic, kp):
    try:
        tid = item["tic"] if isinstance(item, dict) else item
        st, pr = stack_tic(tid, 1.3 if lane == "S106_S107_STACK" else 1.0)
        if not st: return None
        t, f, fe, sc, sap, mx, my = st["t"], st["f"], st["fe"], st["sec"], st["sap"], st["mdx"], st["mdy"]
        if lane == "USP_SIBLING":
            km = np.ones(len(t), bool)
            for sr in e_by_tic.get(tid, []):
                sp, st0, sd = gv(sr, "P_d", "period_d", "period"), gv(sr, "T0_BJD", "t0_bjd", "T0_BTJD", "t0_btjd", "t0"), gv(sr, "T14_h", "dur_h", "duration_h", d=2.2)
                if np.isfinite(sp) and np.isfinite(st0):
                    st0_bt = st0 - BTJD0 if st0 > 2400000 else st0
                    km &= np.abs(((t - st0_bt + 0.5 * sp) % sp) - 0.5 * sp) > (0.9 * sd / 24.0)
            if np.sum(km) < 200: return None
            t, f, fe, sc, sap, mx, my = t[km], f[km], fe[km], sc[km], sap[km], mx[km], my[km]
        rs = float(st["hdr"]["rstar"]) if np.isfinite(st["hdr"]["rstar"]) else 0.7
        pk = search_bls(t, f, fe, max(0.155, 0.75 * item["p"]), min(55.0, 1.35 * item["p"]), rs, seed_p=item["p"]) if lane == "S106_S107_STACK" else search_bls(t, f, fe, 0.155, 0.498, rs)
        if not pk or pk["mes_c"] < 8.6 or (lane == "USP_SIBLING" and any(abs(pk["p"] - p0 / kh) / pk["p"] < 0.01 for p0 in kp.get(tid, []) for kh in range(1, 9))): return None
        gm, nbs = gaia_nb(tid, st["hdr"]["ra"], st["hdr"]["dec"], st["hdr"])
        if str(gm.get("lumclass", "DWARF")).upper() not in ("DWARF", "", "NONE"): return None
        rep = vet_signal(tid, pk, t, f, fe, sc, lane, sap, mx, my, st["crowd"], gm, nbs)
        tpfi = tpf_diff(pr, rep["P_d"], rep["T0_BJD"] - BTJD0, rep["dur_h"] / 24.0) if rep["passed"] else None
        if tpfi: rep = vet_signal(tid, pk, t, f, fe, sc, lane, sap, mx, my, st["crowd"], gm, nbs, tpfi["off"])
        return rep, ({"t": t, "f": f, "sec": sc} if rep["passed"] else None), tpfi
    except Exception: return None
