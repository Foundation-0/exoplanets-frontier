"""Wide-window S1-S107 transit timing extractor, 4-model dBIC ephemeris/TTV solver, and sibling scan."""
from __future__ import annotations
import math
import numpy as np
from bls_search import BTJD0, search_bls

def extract_transit_times(lcs: list[dict], p0: float, t0_btjd: float, dur_d: float, max_shift_min: float = 240.0):
    lim = min(max_shift_min, 0.30 * p0 * 1440.0); grid = np.linspace(-lim, lim, max(61, int(round(lim)) + 1)); sec_r, ind_r = [], []
    for lc in sorted(lcs, key=lambda x: x["sec"]):
        t, f, s = lc["t"], lc["f"], int(lc["sec"])
        if len(t) < 60: continue
        ks = range(math.ceil((float(t[0]) - t0_btjd) / p0) - 1, math.floor((float(t[-1]) - t0_btjd) / p0) + 2)
        cov = [k for k in ks if np.sum(np.abs(t - (t0_btjd + k * p0)) <= (lim / 1440.0 + 0.5 * dur_d)) >= 2]
        if not cov: continue
        ph = ((t - t0_btjd + 0.5 * p0) % p0) - 0.5 * p0; loc = np.abs(ph) > (lim / 1440.0 + 0.8 * dur_d)
        if np.sum(loc) < 25: loc = np.abs(ph) > 1.3 * dur_d
        if np.sum(loc) < 15: continue
        mo = float(np.median(f[loc])); sc = max(1.4826 * float(np.median(np.abs(f[loc] - mo))), 4e-5) * 1e6
        it0 = np.abs(ph) <= 0.5 * dur_d; snr0 = ((mo - float(np.mean(f[it0]))) * 1e6 / (sc / math.sqrt(int(np.sum(it0))))) if np.sum(it0) >= 2 else 0.0
        snrs, deps = np.full(len(grid), -99.0), np.zeros(len(grid))
        for i, dm in enumerate(grid):
            mt = np.abs(ph - dm / 1440.0) <= 0.5 * dur_d; ni = int(np.sum(mt))
            if ni >= 2: deps[i] = (mo - float(np.mean(f[mt]))) * 1e6; snrs[i] = deps[i] / (sc / math.sqrt(ni))
        bi = int(np.argmax(snrs)); bsn, bdt = float(snrs[bi]), float(grid[bi])
        if 0 < bi < len(grid) - 1 and snrs[bi - 1] > -50 and snrs[bi + 1] > -50:
            den = snrs[bi - 1] - 2.0 * bsn + snrs[bi + 1]
            if den < -1e-4: bdt += 0.5 * (snrs[bi - 1] - snrs[bi + 1]) / den * float(grid[1] - grid[0])
        if bsn >= 2.2:
            sec_r.append({"sec": s, "mode": lc.get("mode", "FFI"), "epoch": int(round(float(np.median(cov)))), "n_cov": len(cov),
                          "dt_min": round(bdt, 2), "err_min": round(max(1.5, 0.52 * dur_d * 1440.0 / max(bsn, 1.0)), 2),
                          "snr_peak": round(bsn, 2), "snr_stale": round(snr0, 2), "depth_ppm": round(float(deps[bi]), 1)})
        for k in cov:
            tc0 = t0_btjd + k * p0; sub = np.abs(t - tc0) <= min(0.42 * p0, 4.5 * dur_d + lim / 1440.0)
            if np.sum(sub) < 18: continue
            ts, fs = t[sub], f[sub]; otk = np.abs(ts - (tc0 + bdt / 1440.0)) > 1.25 * dur_d
            if np.sum(otk) < 10: continue
            mk = float(np.median(fs[otk])); sk = max(1.4826 * float(np.median(np.abs(fs[otk] - mk))), 4e-5) * 1e6; bs2, bd2, bdp = -99.0, bdt, 0.0
            for dm in np.linspace(bdt - 55.0, bdt + 55.0, 56):
                ik = np.abs(ts - (tc0 + dm / 1440.0)) <= 0.5 * dur_d
                if np.sum(ik) >= 3:
                    dpk = (mk - float(np.mean(fs[ik]))) * 1e6; snk = dpk / (sk / math.sqrt(int(np.sum(ik))))
                    if snk > bs2: bs2, bd2, bdp = snk, float(dm), dpk
            if bs2 >= 3.0:
                ind_r.append({"sec": s, "mode": lc.get("mode", "FFI"), "epoch": int(k), "n_cov": 1, "dt_min": round(bd2, 2),
                              "err_min": round(max(2.0, 0.55 * dur_d * 1440.0 / bs2), 2), "snr_peak": round(bs2, 2), "snr_stale": round(snr0, 2), "depth_ppm": round(bdp, 1)})
    return sec_r, ind_r

def fit_ephemeris_models(sec_r: list[dict], ind_r: list[dict], p0: float, t0_btjd: float):
    use_i = len(ind_r) >= 6 and float(np.median([r["snr_peak"] for r in ind_r])) >= 3.6
    rows = [r for r in (ind_r if use_i else sec_r) if r["snr_peak"] >= 2.4]; n = len(rows)
    if n < 2:
        return {"n_ep": n, "source": "low", "P_ref": p0, "T0_ref_BJD": t0_btjd + BTJD0, "dP_sec": 0.0, "rms0_m": 0.0, "rms_lin_m": 0.0,
                "ttv_amp_m": 0.0, "ttv_sup_orb": 0.0, "dBIC_lin": 0.0, "dBIC_quad": 0.0, "dBIC_ttv": 0.0, "class": "INSUFFICIENT_EPOCHS"}
    E = np.array([float(r["epoch"]) for r in rows]); y = np.array([float(r["dt_min"]) for r in rows])
    w = 1.0 / (np.maximum(np.array([float(r["err_min"]) for r in rows]), 1.5) ** 2); bic0 = float(np.sum(w * y ** 2))
    Xl = np.column_stack([np.ones(n), E]); cl, *_ = np.linalg.lstsq(Xl * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
    rl = y - Xl @ cl; chil = float(np.sum(w * rl ** 2)); bicl = chil + 2.0 * math.log(n)
    p_ref, t0_ref, dbic_q, dbic_t, amp_t, sup_t = p0 + float(cl[1]) / 1440.0, t0_btjd + float(cl[0]) / 1440.0, 0.0, 0.0, 0.0, 0.0
    if n >= 4 and len(np.unique(E)) >= 4:
        Xq = np.column_stack([np.ones(n), E, 0.5 * E ** 2]); cq, *_ = np.linalg.lstsq(Xq * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
        dbic_q = bicl - (float(np.sum(w * (y - Xq @ cq) ** 2)) + 3.0 * math.log(n))
    if n >= 6 and len(np.unique(E)) >= 6:
        bchi = chil
        for ps in np.geomspace(3.5, max(16.0, 1.3 * float(np.max(E) - np.min(E))), 140):
            Xt = np.column_stack([np.ones(n), E, np.cos(2.0 * np.pi * E / ps), np.sin(2.0 * np.pi * E / ps)])
            ct, *_ = np.linalg.lstsq(Xt * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None); ct_chi = float(np.sum(w * (y - Xt @ ct) ** 2))
            if ct_chi < bchi: bchi, amp_t, sup_t = ct_chi, float(math.hypot(ct[2], ct[3])), float(ps)
        dbic_t = bicl - (bchi + 5.0 * math.log(n))
    dbic_l = bic0 - bicl
    cls = ("DYNAMICAL_TTV_CANDIDATE" if (n >= 6 and dbic_t >= 10.0 and amp_t >= 4.0) else "QUADRATIC_PRECESSION_CANDIDATE" if (n >= 5 and dbic_q >= 10.0) else
           "LINEAR_DRIFT_RESCUED" if (dbic_l >= 10.0 and float(np.max(np.abs(y))) >= 16.0) else "ON_SCHEDULE_STABLE")
    return {"n_ep": n, "source": "per_transit" if use_i else "sector_stacked", "P_ref": round(p_ref, 7), "T0_ref_BJD": round(t0_ref + BTJD0, 6),
            "dP_sec": round(float(cl[1]) * 60.0, 3), "rms0_m": round(float(np.sqrt(np.mean(y ** 2))), 2), "rms_lin_m": round(float(np.sqrt(np.mean(rl ** 2))), 2),
            "ttv_amp_m": round(amp_t, 2), "ttv_sup_orb": round(sup_t, 1), "dBIC_lin": round(dbic_l, 1), "dBIC_quad": round(dbic_q, 1), "dBIC_ttv": round(dbic_t, 1), "class": cls}

def restack_and_siblings(lcs: list[dict], p0: float, t0_bt: float, p_ref: float, t0_ref_bt: float, dur_d: float, rs: float):
    valid = [x for x in lcs if len(x["t"]) >= 60]
    if not valid: return {"snr_stale": 0.0, "snr_ref": 0.0, "sib": None, "single_snr": 0.0, "single_btjd": 0.0}
    t = np.concatenate([x["t"] for x in valid]); f = np.concatenate([x["f"] for x in valid]); fe = np.concatenate([x["fe"] for x in valid])
    o = np.argsort(t); t, f, fe = t[o], f[o], fe[o]
    def _gs(per, t0x):
        ph = ((t - t0x + 0.5 * per) % per) - 0.5 * per; it, ot = np.abs(ph) <= 0.5 * dur_d, np.abs(ph) > 1.3 * dur_d
        if np.sum(it) < 4 or np.sum(ot) < 30: return 0.0
        return float((np.median(f[ot]) - np.mean(f[it])) / (max(1.4826 * float(np.median(np.abs(f[ot] - np.median(f[ot])))), 4e-5) / math.sqrt(int(np.sum(it)))))
    msk = np.abs(((t - t0_ref_bt + 0.5 * p_ref) % p_ref) - 0.5 * p_ref) > 0.85 * dur_d; tm, fm, fem = t[msk], f[msk], fe[msk]
    sib = search_bls(tm, fm, fem, 0.22, min(26.0, float(tm[-1] - tm[0]) / 2.5), rs) if len(tm) >= 200 else None
    b1, bt = 0.0, 0.0
    if len(tm) >= 120:
        med = float(np.median(fm)); sc = max(1.4826 * float(np.median(np.abs(fm - med))), 5e-5)
        for tc in tm[::max(1, len(tm) // 350)]:
            sl = tm[max(0, int(np.searchsorted(tm, tc)) - 4) : int(np.searchsorted(tm, tc)) + 4]
            if len(sl) < 5 or float(np.max(np.diff(sl))) > 0.2: continue
            wi = np.abs(tm - tc) <= 0.06
            if 4 <= np.sum(wi) <= 90:
                s1 = float((med - np.mean(fm[wi])) / (sc / math.sqrt(int(np.sum(wi)))))
                if s1 > b1: b1, bt = s1, float(tc)
    return {"snr_stale": round(_gs(p0, t0_bt), 2), "snr_ref": round(_gs(p_ref, t0_ref_bt), 2), "sib": sib, "single_snr": round(b1, 2), "single_btjd": round(bt, 4)}
