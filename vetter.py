"""Eight-gate false-positive vetter, Sector 106-107 blind verifier, and 4-panel plots."""
from __future__ import annotations
import math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from bls_search import BTJD0, PIX_AS, RSUN_RE, exp_dur, search_bls

def vet_signal(tic, pk, t, f, fe, sec, lane, sap=None, mdx=None, mdy=None, crowd=1.0, meta=None, nbs=None, tpf_off=-1.0):
    meta = meta or {}; rs = float(meta.get("rstar", 0.65) or 0.65)
    teff, tmag = float(meta.get("teff", 4200) or 4200), float(meta.get("tmag", 11.5) or 11.5)
    ruwe, rve = float(meta.get("ruwe", np.nan)), float(meta.get("rv_err_kms", np.nan))
    p, t0, dur, dep, sc = pk["p"], pk["t0"], pk["dur"], pk["dep"], pk["sc"]
    cs = float(np.clip(crowd if np.isfinite(crowd) and crowd > 0 else 1.0, 0.15, 1.0))
    rp = math.sqrt(dep / cs) * rs * RSUN_RE; drt = dur / max(exp_dur(p, rs), 1e-3); vet = []
    if pk["mes_c"] < 9.5: vet.append(f"low_mes({pk['mes_c']:.1f})")
    if pk["sde"] < 7.5: vet.append(f"low_sde({pk['sde']:.1f})")
    if not (0.45 <= rp <= 8.0): vet.append(f"rp({rp:.2f}Re)")
    if not (0.20 <= drt <= 2.6): vet.append(f"dur_ratio({drt:.2f})")
    if any(abs(p - a) / a < 0.008 for a in (0.5, 1.0, 6.85, 13.7)) and pk["mes_c"] < 15.0: vet.append("alias")

    ph = ((t - t0 + 0.5 * p) % p) - 0.5 * p; orb = np.round((t - t0) / p).astype(int)
    it, ot = np.abs(ph) <= 0.5 * dur, np.abs(ph) > 1.25 * dur
    core = np.abs(ph) <= 0.22 * dur; wing = (np.abs(ph) > 0.28 * dur) & (np.abs(ph) <= 0.5 * dur)
    if np.sum(core) >= 4 and np.sum(wing) >= 4 and dep > 0.0025:
        d_core = float(np.median(f[ot]) - np.mean(f[core])); d_wing = float(np.median(f[ot]) - np.mean(f[wing]))
        if d_core > 0 and (d_wing / d_core) < 0.28: vet.append("grazing_v")

    io_, ie_ = it & ((orb % 2) != 0), it & ((orb % 2) == 0); oes = 0.0
    if np.sum(io_) >= 3 and np.sum(ie_) >= 3:
        do, de = float(np.median(f[ot]) - np.mean(f[io_])), float(np.median(f[ot]) - np.mean(f[ie_]))
        oes = abs(do - de) / math.hypot(sc / math.sqrt(int(np.sum(io_))), sc / math.sqrt(int(np.sum(ie_))))
        if do <= 0 or de <= 0 or oes > 2.8: vet.append(f"odd_even({oes:.2f}s)")

    ipk = search_bls(t, 2.0 * float(np.median(f)) - f, fe, max(0.15, 0.72 * p), min(65.0, 1.38 * p), rs, fast_inv=True)
    if ipk and ipk["mes_c"] > 0.70 * pk["mes_c"] and ipk["mes_c"] > 7.5: vet.append(f"inv_lc({ipk['mes_c']:.1f})")

    ssnr = max([float(np.median(f[ot]) - np.mean(f[sm])) / (sc / math.sqrt(int(np.sum(sm)))) for c in np.linspace(0.18, 0.82, 31) for sm in [np.abs(((t - t0) / p) % 1.0 - c) <= (0.5 * dur / p)] if np.sum(sm) >= 4] or [0.0])
    if ssnr > 4.2 and ssnr > 0.25 * pk["mes"]: vet.append(f"sec_ecl({ssnr:.1f})")

    evq = [(max(0.0, float(np.median(f[ot]) - np.mean(f[it & (orb == o)]))) / (sc / math.sqrt(int(np.sum(it & (orb == o)))))) ** 2 for o in np.unique(orb[it]) if np.sum(it & (orb == o)) >= 1]
    if (max(evq) / max(sum(evq), 1e-9) if evq else 1.0) > 0.55: vet.append("single_ev")
    usec = np.unique(sec[it])
    if len(usec) < 2: vet.append("1_sector")

    if sap is not None and len(sap) == len(f):
        vs = np.isfinite(sap) & (sap > 0)
        if np.sum(vs & it) >= 5 and np.sum(vs & ot) >= 30:
            los = vs & (np.abs(ph) > 1.2 * dur) & (np.abs(ph) <= min(3.5 * dur, 0.45 * p))
            if np.sum(los) < 15: los = vs & ot
            ms = float(np.median(sap[los])); ds = (ms - float(np.mean(sap[vs & it]))) / max(ms, 1e-6)
            if dep > 0 and ds > 0 and not (0.25 <= (ds / cs) / dep <= 2.6): vet.append("sap_mismatch")

    ca = csig = 0.0
    if mdx is not None and mdy is not None and len(mdx) == len(f):
        wx = wy = wsx = wsy = 0.0
        for s in usec:
            si = it & (sec == s) & np.isfinite(mdx) & np.isfinite(mdy)
            sl = (sec == s) & (np.abs(ph) > 1.2 * dur) & (np.abs(ph) <= min(4.0 * dur, 0.45 * p)) & np.isfinite(mdx) & np.isfinite(mdy)
            if np.sum(si) >= 3 and np.sum(sl) >= 15:
                sxs = max(1.4826 * float(np.median(np.abs(mdx[sl] - np.median(mdx[sl])))), 1e-4) / math.sqrt(int(np.sum(si)))
                sys_ = max(1.4826 * float(np.median(np.abs(mdy[sl] - np.median(mdy[sl])))), 1e-4) / math.sqrt(int(np.sum(si)))
                wx += float(np.mean(mdx[si]) - np.median(mdx[sl])) / sxs ** 2; wy += float(np.mean(mdy[si]) - np.median(mdy[sl])) / sys_ ** 2
                wsx += 1.0 / sxs ** 2; wsy += 1.0 / sys_ ** 2
        if wsx > 0 and wsy > 0:
            mx, my = wx / wsx, wy / wsy; sh = math.hypot(mx, my) * PIX_AS
            csig = math.hypot(mx * math.sqrt(wsx), my * math.sqrt(wsy)); ca = (sh / max(dep, 2e-4)) if csig > 2.2 else sh
            if csig > 3.5 and ca > 12.0: vet.append(f"cent({ca:.1f}as)")
    if tpf_off >= 0 and tpf_off > 15.0: vet.append(f"tpf({tpf_off:.1f}as)")
    if np.isfinite(ruwe) and ruwe > 1.4: vet.append(f"ruwe({ruwe:.2f})")
    if np.isfinite(rve) and rve > 5.0: vet.append(f"rv({rve:.1f})")

    cnb, nfpp, effc = 0, 0.0, (tpf_off if tpf_off >= 0 else ca)
    for nb in (nbs or []):
        sp = float(nb.get("sep_arcsec", 99.0)); ff = 10.0 ** (-0.4 * (float(nb.get("tmag", 18.0)) - tmag)); pc = math.exp(-0.5 * (sp / 24.0) ** 2)
        if 0.5 * (ff * pc) / (1.0 + ff * pc) >= 0.7 * dep:
            cnb += 1; nfpp += 0.012 * pc * (math.exp(-0.5 * ((sp - effc) / 13.0) ** 2) if (csig > 2.0 or tpf_off > 8.0) else 0.18) * 0.22
    heb = 0.0025 * math.exp(0.5 * min(oes, 3.0) ** 2 - 0.5) * (2.2 if (np.isfinite(ruwe) and ruwe > 1.25) else 1.0) * math.exp(-0.22 * max(0.0, pk["mes_c"] - 9.0))
    nfpp = float(np.clip(nfpp, 0.0, 0.99)); fpp = float(np.clip(heb + nfpp, 1e-5, 0.99))
    if fpp >= 0.50 or nfpp >= 0.10: vet.append(f"fpp({fpp:.3f})")

    insol = ((rs ** 2) * ((teff / 5772.0) ** 4)) / max(((max(rs ** 0.9, 0.08) * (p / 365.25) ** 2) ** (1.0 / 3.0)) ** 2, 1e-6)
    teq = 278.5 * (insol ** 0.25) * (0.7 ** 0.25)
    mp, sf = (max(0.15, rp ** 3.5), 0.190) if rp < 1.5 else (2.69 * (rp ** 0.93), 1.26)
    tsm = sf * ((rp ** 3) * teq / (mp * rs ** 2)) * (10.0 ** (-0.2 * (tmag - (0.75 if teff < 4200 else 0.45))))
    ok = len(vet) == 0
    tier = "rejected" if not ok else ("Tier N" if (cnb > 0 and nfpp >= 0.001) else ("Tier S-Val" if (fpp < 0.015 and nfpp < 0.001 and pk["mes_c"] >= 10.0) else "Tier S"))
    return {"tic": tic, "lane": lane, "P_d": round(p, 7), "T0_BJD": round(t0 + BTJD0, 6), "depth_ppm": round(dep * 1e6, 1),
            "err_ppm": round((sc / math.sqrt(max(int(np.sum(it)), 1))) * 1e6, 1), "dur_h": round(dur * 24.0, 2), "Rp_Re": round(rp, 2),
            "MES": round(pk["mes_c"], 2), "SDE": round(pk["sde"], 1), "n_tr": pk["ntr"], "n_sec": int(len(usec)), "oe_sig": round(oes, 2),
            "cent_as": round(effc, 1), "ruwe": round(ruwe, 3) if np.isfinite(ruwe) else -1.0, "Teff": round(teff, 0), "Rstar": round(rs, 3),
            "Tmag": round(tmag, 2), "Teq_K": round(teq, 0), "TSM": round(float(tsm), 1), "FPP": round(fpp, 5), "NFPP": round(nfpp, 5),
            "passed": ok, "tier": tier, "vetoes": ";".join(vet) if vet else "CLEAN"}

defde blind_check_sector(tic, tier, p_d, t0_bjd, pred_ppm, dur_h, sec, t, f, mdx=None, mdy=None):
    t0 = t0_bjd - BTJD0 if t0_bjd > 2400000 else t0_bjd; dur = max(dur_h / 24.0, 0.025)
    if len(t) < 50: return {"tic": tic, "tier": tier, "sector": sec, "P_d": round(p_d, 6), "pred_ppm": round(pred_ppm, 1), "obs_ppm": 0.0, "err_ppm": 0.0, "snr": 0.0, "dt_min": 0.0, "cent_sig": 0.0, "n_cov": 0, "verdict": "NO_DATA"}
    preds = [t0 + k * p_d for k in range(math.ceil((float(t[0]) - t0) / p_d), math.floor((float(t[-1]) - t0) / p_d) + 1)]
    cov = [tc for tc in preds if np.sum(np.abs(t - tc) <= 0.5 * dur) >= 2]
    if not cov: return {"tic": tic, "tier": tier, "sector": sec, "P_d": round(p_d, 6), "pred_ppm": round(pred_ppm, 1), "obs_ppm": 0.0, "err_ppm": 0.0, "snr": 0.0, "dt_min": 0.0, "cent_sig": 0.0, "n_cov": 0, "verdict": "DATA_GAP"}
    ph = ((t - t0 + 0.5 * p_d) % p_d) - 0.5 * p_d; it = np.abs(ph) <= 0.5 * dur
    loc = (np.abs(ph) > 1.2 * dur) & (np.abs(ph) <= min(4.0 * dur, 0.45 * p_d))
    if np.sum(loc) < 25: loc = np.abs(ph) > 1.2 * dur
    sc = max(1.4826 * float(np.median(np.abs(f[loc] - np.median(f[loc])))), 3e-5) * 1e6
    ob = float(np.median(f[loc]) - np.mean(f[it])) * 1e6; er = sc * math.sqrt(1.0 / max(int(np.sum(it)), 1) + 1.0 / max(int(np.sum(loc)), 1))
    snr = ob / max(er, 1.0); bdt, bsn = 0.0, -99.0
    for dm in np.linspace(-60.0, 60.0, 61):
        mt = np.abs(ph - dm / 1440.0) <= 0.5 * dur
        if np.sum(mt) >= 2:
            st = (float(np.median(f[loc]) - np.mean(f[mt])) * 1e6) / (sc / math.sqrt(int(np.sum(mt))))
            if st > bsn: bsn, bdt = st, float(dm)
    csig = 0.0
    if mdx is not None and mdy is not None and len(mdx) == len(f):
        vc = np.isfinite(mdx) & np.isfinite(mdy)
        if np.sum(vc & it) >= 3 and np.sum(vc & loc) >= 20:
            sx = max(1.4826 * float(np.median(np.abs(mdx[vc & loc] - np.median(mdx[vc & loc])))), 1e-4)
            sy = max(1.4826 * float(np.median(np.abs(mdy[vc & loc] - np.median(mdy[vc & loc])))), 1e-4)
            csig = math.hypot(float(np.mean(mdx[vc & it] - np.median(mdx[vc & loc]))) / (sx / math.sqrt(int(np.sum(vc & it)))),
                              float(np.mean(mdy[vc & it] - np.median(mdy[vc & loc]))) / (sy / math.sqrt(int(np.sum(vc & it)))))
    dr, defs = ob / max(pred_ppm, 1.0), (pred_ppm - ob) / max(er, 1.0)
    vd = ("CONFIRMED_ON_SCHEDULE" if (snr >= 3.0 and 0.42 <= dr <= 2.15 and csig < 3.2) else
          "TIMING_OFFSET_TTV" if (bsn >= 3.5 and abs(bdt) >= 18.0 and csig < 3.2) else
          "NON_DETECTION_4SIGMA" if (pred_ppm / max(er, 1.0) >= 5.0 and snr < 2.0 and defs >= 4.0 and bsn < 3.0) else
          "MARGINAL_RECOVERY" if (snr >= 2.0 and 0.38 <= dr <= 2.25) else "LOW_SNR_INCONCLUSIVE")
    return {"tic": tic, "tier": tier, "sector": sec, "P_d": round(p_d, 6), "pred_ppm": round(pred_ppm, 1), "obs_ppm": round(ob, 1),
            "err_ppm": round(er, 1), "snr": round(snr, 2), "dt_min": round(bdt, 1), "cent_sig": round(csig, 2), "n_cov": len(cov), "verdict": vd}

def plot_card(out_p, tic, tag, p_d, t0_btjd, dur_d, dep_ppm, t, f, sec, tpf=None):
    out_p.parent.mkdir(parents=True, exist_ok=True)
    ph = (((t - t0_btjd + 0.5 * p_d) % p_d) - 0.5 * p_d) * 24.0; orb = np.round((t - t0_btjd) / p_d).astype(int); fp = (f - np.median(f)) * 1e6
    nc = 4 if tpf else 3; fig, ax = plt.subplots(1, nc, figsize=(4.6 * nc, 4.0), dpi=140); fig.patch.set_facecolor("#0f172a")
    for a in ax: a.set_facecolor("#1e293b"); a.tick_params(colors="#cbd5e1", labelsize=8)
    wh = min(max(3.5 * dur_d * 24.0, 2.2), 0.45 * p_d * 24.0); zm = np.abs(ph) <= wh; bins = np.linspace(-wh, wh, 27); bi = np.digitize(ph[zm], bins) - 1
    ax[0].scatter(ph[zm], fp[zm], s=4, alpha=0.22, color="#94a3b8")
    pts = [(0.5 * (bins[b] + bins[b + 1]), float(np.mean(fp[zm][bi == b])), float(np.std(fp[zm][bi == b]) / math.sqrt(int(np.sum(bi == b))))) for b in range(26) if np.sum(bi == b) >= 3]
    if pts: ax[0].errorbar([x for x, _, _ in pts], [y for _, y, _ in pts], yerr=[e for _, _, e in pts], fmt="o", color="#38bdf8", ms=4, lw=1)
    ax[0].axhline(-dep_ppm, color="#f43f5e", ls="--", lw=1.2); ax[0].set_xlim(-wh, wh); ax[0].set_title(f"Phase Fold P={p_d:.5f}d", color="#f8fafc", fontsize=9)
    for par, col in [(1, "#34d399"), (0, "#fbbf24")]:
        pm = zm & ((orb % 2) == par); pbi = np.digitize(ph[pm], bins) - 1
        op = [(0.5 * (bins[b] + bins[b + 1]), float(np.mean(fp[pm][pbi == b]))) for b in range(26) if np.sum(pbi == b) >= 2]
        if op: ax[1].plot([x for x, _ in op], [y for _, y in op], "o-", color=col, ms=3.8, lw=1)
    ax[1].axhline(-dep_ppm, color="#f43f5e", ls="--", lw=1); ax[1].set_xlim(-wh, wh); ax[1].set_title("Odd (green) vs Even (gold)", color="#f8fafc", fontsize=9)
    oph = ((t - t0_btjd + 0.5 * p_d) % p_d) / p_d - 0.5; fb = np.linspace(-0.5, 0.5, 55); fbi = np.digitize(oph, fb) - 1
    fpt = [(0.5 * (fb[b] + fb[b + 1]), float(np.mean(fp[fbi == b]))) for b in range(54) if np.sum(fbi == b) >= 4]
    if fpt: ax[2].plot([x for x, _ in fpt], [y for _, y in fpt], "o-", color="#a78bfa", ms=3, lw=1)
    ax[2].axvline(0.5, color="#fb7185", ls=":"); ax[2].set_xlim(-0.5, 0.5); ax[2].set_title(f"Full Orbit ({len(np.unique(sec))} sec)", color="#f8fafc", fontsize=9)
    if tpf:
        ax[3].imshow(tpf["diff"], origin="lower", cmap="inferno")
        ax[3].plot(*tpf["sxy"], "+", color="#38bdf8", ms=10, mew=2); ax[3].plot(*tpf["dxy"], "x", color="#34d399", ms=8, mew=2)
        ax[3].set_title(f"S{tpf['sec']} Diff ({tpf['off']:.1f}\")", color="#f8fafc", fontsize=9)
    
