"""Coherent S1-S108 multi-sector BLS, starspot pre-whitening & WCS blend core for Issue 5."""
from __future__ import annotations
import concurrent.futures as cf, math
from pathlib import Path
import numpy as np
from blend_gaia import multisector_wcs_blend_audit
from bls_search import detrend_lc, exp_dur, pont_beta, search_bls
from converge_io import discover_all_sectors, fetch_raw_spoc
from deep_search import ffi_cutout_lc
from despot import harmonic_despot, plot_despot_card
from mdwarf_hunt import norm_sap, spot_amp
from vetter import vet_signal


def load_multi(tid: int, ra: float, dec: float):
    raw = fetch_raw_spoc(tid, max_sec=14); ok = np.isfinite(ra) and np.isfinite(dec)
    alls = discover_all_sectors(ra, dec) if ok else [x["sec"] for x in raw]
    if len(raw) < 3 and ok and alls:
        have, n = {x["sec"] for x in raw}, len(alls)
        want = [s for s in sorted({alls[i] for i in [0, min(1, n - 1), n // 3, (2 * n) // 3, max(0, n - 2), n - 1]}) if s not in have][:6]
        with cf.ThreadPoolExecutor(max_workers=3) as ex:
            for r in ex.map(lambda s: ffi_cutout_lc(tid, ra, dec, s), want):
                if r: raw.append({**r, "sap": r["f"].copy(), "crowd": 0.92})
    return sorted(raw, key=lambda x: x["sec"]), alls


def stack_dual(raw: list[dict], win_s: float, pol_d: float):
    def _cat(dsp_mode: bool):
        ta, fa, fea, sa, sapa, mxa, mya, crs, infs = [], [], [], [], [], [], [], [], []
        for l in raw:
            a, b, c, *rest = harmonic_despot(l["t"], l["f"], l["fe"], l["q"], polish_d=pol_d) if dsp_mode else (*detrend_lc(l["t"], l["f"], l["fe"], l["q"], win_s), None)
            if len(a) >= 60:
                im = np.clip(np.searchsorted(l["t"], a), 0, len(l["t"]) - 1)
                ta.append(a); fa.append(b); fea.append(c); sa.append(np.full(len(a), l["sec"], int))
                sapa.append(l["sap"][im]); mxa.append(l["mdx"][im]); mya.append(l["mdy"][im]); crs.append(l["crowd"])
                if rest[0]: infs.append(rest[0])
        if not ta: return None
        sc = np.concatenate(sa)
        return {"t": np.concatenate(ta), "f": np.concatenate(fa), "fe": np.concatenate(fea), "sec": sc,
                "sap": norm_sap(np.concatenate(sapa), sc), "mdx": np.concatenate(mxa), "mdy": np.concatenate(mya),
                "crowd": float(np.median(crs)) if crs else 0.95, "infs": infs}
    return _cat(False), _cat(True)


def coherent_seed_bls(t, f, fe, pc: float, rs: float):
    from astropy.timeseries import BoxLeastSquares
    if len(t) < 120 or (t[-1] - t[0]) < 1.35 * pc: return None
    p_lo, p_hi = 0.982 * pc, min(1.018 * pc, float(t[-1] - t[0]) / 1.8)
    if p_hi <= p_lo * 1.002: return None
    d0 = max(0.018, 0.38 * exp_dur(pc, rs))
    durs = np.geomspace(d0, min(0.35 * p_lo, max(d0 * 1.2, 1.7 * exp_dur(pc, rs))), 7)
    bls = BoxLeastSquares(t, f, dy=fe)
    ps = bls.autoperiod(durs, minimum_period=p_lo, maximum_period=p_hi, frequency_factor=2.5)
    if len(ps) > 48000: ps = 1.0 / np.linspace(1.0 / p_hi, 1.0 / p_lo, 48000)[::-1]
    res = bls.power(ps, durs[durs < 0.85 * p_lo], objective="snr")
    pw = np.where(np.asarray(res.depth) > 0, np.asarray(res.power, float), 0.0)
    bi = int(np.nanargmax(pw)); mp = float(np.nanmedian(pw))
    sde = float((np.nanmax(pw) - mp) / max(1.4826 * float(np.nanmedian(np.abs(pw - mp))), 1e-6))
    bp, bt0, bdur = float(res.period[bi]), float(res.transit_time[bi]), float(res.duration[bi])
    st = bls.compute_stats(bp, bdur, bt0)
    ph = ((t - bt0 + 0.5 * bp) % bp) - 0.5 * bp; it, ot = np.abs(ph) <= 0.5 * bdur, np.abs(ph) > 1.2 * bdur
    ntr = int(np.sum(np.asarray(st["per_transit_count"]) >= 2)) or int(np.sum(np.asarray(st["per_transit_count"]) >= 1))
    if np.sum(it) < 4 or np.sum(ot) < 30 or ntr < 2: return None
    dep = float(np.median(f[ot]) - np.mean(f[it]))
    if dep <= 0: return None
    sc = max(1.4826 * float(np.median(np.abs(f[ot] - np.median(f[ot])))), 1e-6)
    mes = dep / (sc * math.sqrt(1.0 / int(np.sum(it)) + 1.0 / int(np.sum(ot))))
    mc = mes / pont_beta(t[ot], f[ot] - np.median(f[ot]), bdur)
    return {"p": bp, "t0": bt0, "dur": bdur, "dep": dep, "sc": sc, "mes": mes, "mes_c": mc, "sde": sde, "ntr": ntr, "stats": st}


def audit_target(tid: int, pc: float, lane: str, tr, hits: list[dict], kt: set, kp: dict, fig_d: Path) -> dict:
    mh = [h for h in hits if h["tic"] == tid and abs(h["p"] - pc) < 0.20] or [h for h in hits if h["tic"] == tid]
    sfs = ";".join(sorted({h["rel"] for h in mh})) if mh else "none"
    ut = ";".join(sorted({h["tier"] for h in mh if h["tier"] and h["tier"] != "nan"})) or "un-tiered"
    tf = lambda k, d: float(tr[k]) if (tr is not None and k in tr.colnames and tr[k] is not None and np.isfinite(tr[k])) else d
    ra, dec = tf("ra", mh[0]["ra"] if mh else np.nan), tf("dec", mh[0]["dec"] if mh else np.nan)
    rs, teff, tmag = max(tf("rad", 0.55), 0.08), max(tf("Teff", 4200.0), 2500.0), tf("Tmag", 10.5)
    raw, alls = load_multi(tid, ra, dec)
    if len(raw) < 2:
        return {"tic": tid, "lane": lane, "src_files": sfs, "upstream_tier": ut, "known_nasa": tid in kt, "Tmag": round(tmag, 2),
                "P_cand_d": round(pc, 5), "n_sec_used": len(raw), "n_transits": 0, "p_rot_d": 0.0, "std_MES": 0.0, "despot_MES": 0.0,
                "Rp_Re": -1.0, "Teq_K": -1.0, "wcs_off_as": -1.0, "plausible_bebs": -1, "survived_hz_planet": False, "vetoes": "insufficient_sectors"}
    std, dsp = stack_dual(raw, max(1.8, min(3.4, 9.0 * exp_dur(pc, rs))), max(1.6, 7.5 * exp_dur(pc, rs)))
    pks = coherent_seed_bls(std["t"], std["f"], std["fe"], pc, rs) if std else None
    pkd = coherent_seed_bls(dsp["t"], dsp["f"], dsp["fe"], pc, rs) if dsp else None
    prot = float(np.median([x["p_rot_d"] for x in (dsp["infs"] if dsp else []) if x["p_rot_d"] > 0] or [0.0]))
    pk, st = (pkd, dsp) if (pkd and (not pks or pkd["mes_c"] >= 0.85 * pks["mes_c"])) else (pks, std)
    wr, _ = multisector_wcs_blend_audit(tid, ra, dec, pk["p"], pk["t0"], pk["dur"], pk["dep"] * 1e6, alls[-4:]) if (pk and np.isfinite(ra)) else ({"off_target_as": -1.0, "verdict": "NONE", "ruwe": 1.0, "plausible_nbs": -1}, None)
    meta = {"rstar": rs, "teff": teff, "tmag": tmag, "ruwe": float(wr.get("ruwe", 1.0))}
    rep = vet_signal(tid, pk, st["t"], st["f"], st["fe"], st["sec"], lane, st["sap"], st["mdx"], st["mdy"], st["crowd"], meta, tpf_off=float(wr.get("off_target_as", -1.0))) if (pk and st) else None
    if rep and st:
        amp = spot_amp(st["t"], st["f"], pk["p"], pk["t0"], pk["dur"])
        if amp > 0.48 * pk["dep"]: rep["passed"] = False; rep["vetoes"] = ("" if rep["vetoes"] == "CLEAN" else rep["vetoes"] + ";") + f"spot_wave({amp/max(pk['dep'],1e-6):.2f})"
        if wr.get("verdict") in ("NEIGHBOR_BEB_LOCALIZED", "OFFSET_MARGINAL") and float(wr.get("off_target_as", -1.0)) > 12.0: rep["passed"] = False; rep["vetoes"] = ("" if rep["vetoes"] == "CLEAN" else rep["vetoes"] + ";") + f"wcs({wr['verdict']})"
    pka = search_bls(dsp["t"], dsp["f"], dsp["fe"], 0.35, 27.0, rs) if dsp else None
    repa = vet_signal(tid, pka, dsp["t"], dsp["f"], dsp["fe"], dsp["sec"], lane + "_BLIND", dsp["sap"], dsp["mdx"], dsp["mdy"], dsp["crowd"], meta) if pka else None
    if repa and dsp and (spot_amp(dsp["t"], dsp["f"], pka["p"], pka["t0"], pka["dur"]) > 0.48 * pka["dep"] or (prot > 0 and any(abs(pka["p"] - k * prot) / (k * prot) < 0.025 for k in (0.5, 1.0, 2.0)))): repa["passed"] = False
    if dsp and (pk or pka):
        px, vt = (pk if pk else pka), (rep["vetoes"] if rep else "no_multi_sector_dip")
        plot_despot_card(fig_d / f"hz_TIC{tid}_P{int(round(pc))}d.png", tid, f"{lane} Seed={pc:.2f}d MES={(pk['mes_c'] if pk else 0.0):.1f} ({vt})", raw[0]["t"], raw[0]["f"], dsp["t"], dsp["f"], px["p"], px["t0"], px["dur"])
    return {"tic": tid, "lane": lane, "src_files": sfs, "upstream_tier": ut, "known_nasa": bool(tid in kt), "known_P_d": round(kp[tid][0], 4) if tid in kp and kp[tid] else -1.0,
            "Tmag": round(tmag, 2), "Teff": round(teff, 0), "Rstar": round(rs, 3), "P_cand_d": round(pc, 5), "n_sec_used": len(raw), "n_transits": int(pk["ntr"]) if pk else 0,
            "p_rot_d": round(prot, 4), "std_MES": round(pks["mes_c"], 2) if pks else 0.0, "despot_MES": round(pkd["mes_c"], 2) if pkd else 0.0,
            "P_fit_d": rep["P_d"] if rep else -1.0, "depth_ppm": rep["depth_ppm"] if rep else 0.0, "Rp_Re": rep["Rp_Re"] if rep else -1.0, "Teq_K": rep["Teq_K"] if rep else -1.0,
            "wcs_off_as": wr.get("off_target_as", -1.0), "plausible_bebs": wr.get("plausible_nbs", -1), "wcs_verdict": wr.get("verdict", "NONE"),
            "blind_alt_P_d": repa["P_d"] if repa else -1.0, "blind_alt_MES": repa["MES"] if repa else 0.0, "blind_alt_passed": bool(repa and repa["passed"]),
            "survived_hz_planet": bool(rep and rep["passed"]), "vetoes": rep["vetoes"] if rep else "no_multi_sector_dip"}
