"""Cycle 1-8 SPOC/FFI sector loader and broadband-calibrated coherent BLS for cool worlds."""
from __future__ import annotations
import concurrent.futures as cf, math
import numpy as np
from bls_search import detrend_lc, exp_dur, pont_beta
from converge_io import discover_all_sectors, fetch_raw_spoc
from deep_search import ffi_cutout_lc
from despot import harmonic_despot


def load_cycles_1_to_8(tid: int, ra: float, dec: float):
    raw = fetch_raw_spoc(tid, max_sec=8)
    for x in raw:
        x["mode"] = "SPOC"
    ok = np.isfinite(ra) and np.isfinite(dec)
    alls = discover_all_sectors(ra, dec) if ok else sorted({x["sec"] for x in raw})
    have = {x["sec"] for x in raw}
    if ok and alls and (len(raw) < 4 or alls[0] not in have or alls[-1] not in have):
        n = len(alls)
        ix = sorted({0, min(1, n - 1), n // 4, n // 2, (3 * n) // 4, max(0, n - 2), n - 1})
        want = [alls[i] for i in ix if alls[i] not in have][:6]
        with cf.ThreadPoolExecutor(max_workers=4) as ex:
            for r in ex.map(lambda s: ffi_cutout_lc(tid, ra, dec, s), want):
                if r:
                    raw.append({**r, "sap": np.full(len(r["f"]), np.nan, dtype=float), "crowd": np.nan, "mode": "FFI"})
    return sorted(raw, key=lambda x: x["sec"]), alls


def calibrated_bls_sde(pw_bg, pw_fn) -> float:
    """Compute conservative SDE calibrated across both broadband and fine grids."""
    fn = np.asarray(pw_fn, dtype=float)
    fn = fn[np.isfinite(fn)]
    if len(fn) == 0:
        return 0.0
    peak = float(np.max(fn))
    if not (np.isfinite(peak) and peak > 0.0):
        return 0.0
    mp_f = float(np.median(fn))
    mad_f = max(1.4826 * float(np.median(np.abs(fn - mp_f))), 1e-6)
    bg = np.asarray(pw_bg, dtype=float)
    bg = bg[np.isfinite(bg)]
    if len(bg) >= 20:
        mp_b = float(np.median(bg))
        mad_b = max(1.4826 * float(np.median(np.abs(bg - mp_b))), 1e-6)
        mp_ref = max(mp_b, mp_f)
        mad_ref = max(mad_b, mad_f, 1e-6)
    else:
        mp_ref, mad_ref = mp_f, mad_f
    sde = float((peak - mp_ref) / mad_ref)
    return sde if np.isfinite(sde) else 0.0


def stack_clean(raw: list[dict], win_s: float, pol_d: float):
    lcs_s, ts, fs, fes, scs, td, fd, fed, scd, infs = [], [], [], [], [], [], [], [], [], []
    for l in raw:
        a, b, c = detrend_lc(l["t"], l["f"], l["fe"], l["q"], win_s)
        if len(a) >= 55:
            ts.append(a); fs.append(b); fes.append(c); scs.append(np.full(len(a), l["sec"], int))
            lcs_s.append({"sec": l["sec"], "t": a, "f": b, "fe": c, "mode": l.get("mode", "FFI")})
        a2, b2, c2, inf = harmonic_despot(l["t"], l["f"], l["fe"], l["q"], polish_d=pol_d)
        if len(a2) >= 55:
            td.append(a2); fd.append(b2); fed.append(c2); scd.append(np.full(len(a2), l["sec"], int)); infs.append(inf)
    std = {"t": np.concatenate(ts), "f": np.concatenate(fs), "fe": np.concatenate(fes), "sec": np.concatenate(scs)} if ts else None
    dsp = {"t": np.concatenate(td), "f": np.concatenate(fd), "fe": np.concatenate(fed), "sec": np.concatenate(scd), "infs": infs} if td else None
    return lcs_s, std, dsp


def coherent_planet_bls(t, f, fe, pc: float, rs: float):
    from astropy.timeseries import BoxLeastSquares
    if not (np.isfinite(pc) and pc > 0 and np.isfinite(rs) and rs > 0):
        return None
    t, f, fe = np.asarray(t, float), np.asarray(f, float), np.asarray(fe, float)
    valid = np.isfinite(t) & np.isfinite(f) & np.isfinite(fe) & (fe > 0)
    t, f, fe = t[valid], f[valid], fe[valid]
    if len(t) < 100:
        return None
    o = np.argsort(t); t, f, fe = t[o], f[o], fe[o]
    fb = float(t[-1] - t[0])
    p_lo, p_hi = 0.985 * pc, min(1.015 * pc, fb / 1.65)
    if fb < 1.32 * pc or p_hi <= p_lo * 1.001:
        return None
    d0 = max(0.020, 0.40 * exp_dur(pc, rs))
    durs = np.geomspace(d0, min(0.35 * p_lo, max(d0 * 1.25, 1.75 * exp_dur(pc, rs))), 7)
    vdf = durs[durs < 0.85 * p_lo]
    durs_use = vdf if len(vdf) else np.array([0.12 * p_lo], dtype=float)
    bls = BoxLeastSquares(t, f, dy=fe)
    bg_lo = max(0.75, 0.42 * pc, float(np.max(durs_use)) / 0.84)
    bg_hi = min(34.0, fb / 1.65)
    pw_bg = np.array([], dtype=float)
    if bg_hi > bg_lo * 1.05:
        p_bg = 1.0 / np.linspace(1.0 / bg_hi, 1.0 / bg_lo, 2800)[::-1]
        p_bg = p_bg[(p_bg < p_lo) | (p_bg > p_hi)]
        if len(p_bg) >= 20:
            r_bg = bls.power(p_bg, durs_use, objective="snr")
            pw_bg = np.where(np.isfinite(r_bg.depth) & (np.asarray(r_bg.depth) > 0) & np.isfinite(r_bg.power), np.asarray(r_bg.power, float), 0.0)
    p_fn = 1.0 / np.linspace(1.0 / p_hi, 1.0 / p_lo, 12000)[::-1]
    r_fn = bls.power(p_fn, durs_use, objective="snr")
    pw_f = np.where(np.isfinite(r_fn.depth) & (np.asarray(r_fn.depth) > 0) & np.isfinite(r_fn.power), np.asarray(r_fn.power, float), 0.0)
    if len(pw_f) == 0 or not np.any(pw_f > 0):
        return None
    bi = int(np.nanargmax(pw_f)); bp, bt0, bdur = float(r_fn.period[bi]), float(r_fn.transit_time[bi]), float(r_fn.duration[bi])
    if not (np.isfinite(bp) and bp > 0 and np.isfinite(bt0) and np.isfinite(bdur) and bdur > 0):
        return None
    st = bls.compute_stats(bp, bdur, bt0)
    ph = ((t - bt0 + 0.5 * bp) % bp) - 0.5 * bp; it, ot = np.abs(ph) <= 0.5 * bdur, np.abs(ph) > 1.25 * bdur
    ntr = int(np.sum(np.asarray(st["per_transit_count"]) >= 2)) or int(np.sum(np.asarray(st["per_transit_count"]) >= 1))
    if np.sum(it) < 4 or np.sum(ot) < 30 or ntr < 2:
        return None
    dep = float(np.median(f[ot]) - np.mean(f[it]))
    if not (np.isfinite(dep) and dep > 0):
        return None
    sc = max(1.4826 * float(np.median(np.abs(f[ot] - np.median(f[ot])))), 1e-6)
    mes = dep / (sc * math.sqrt(1.0 / int(np.sum(it)) + 1.0 / int(np.sum(ot))))
    beta = max(1.0, float(pont_beta(t[ot], f[ot] - np.median(f[ot]), bdur)))
    mes_c = mes / beta
    sde = calibrated_bls_sde(pw_bg, pw_f)
    if not (np.isfinite(mes) and np.isfinite(mes_c) and np.isfinite(sde)):
        return None
    return {"p": bp, "t0": bt0, "dur": bdur, "dep": dep, "sc": sc, "mes": mes,
            "mes_c": mes_c, "sde": sde, "ntr": ntr}
