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
                    raw.append({**r, "sap": r["f"].copy(), "crowd": 0.94, "mode": "FFI"})
    return sorted(raw, key=lambda x: x["sec"]), alls


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
    fb = float(t[-1] - t[0]) if len(t) >= 100 else 0.0
    p_lo, p_hi = 0.985 * pc, min(1.015 * pc, fb / 1.65)
    if fb < 1.32 * pc or p_hi <= p_lo * 1.001:
        return None
    d0 = max(0.020, 0.40 * exp_dur(pc, rs))
    durs = np.geomspace(d0, min(0.35 * p_lo, max(d0 * 1.25, 1.75 * exp_dur(pc, rs))), 7)
    bls = BoxLeastSquares(t, f, dy=fe)
    bg_lo, bg_hi = max(0.75, 0.42 * pc), min(34.0, fb / 1.65)
    p_bg = 1.0 / np.linspace(1.0 / max(bg_hi, bg_lo * 1.2), 1.0 / bg_lo, 2800)[::-1]
    vd = durs[durs < 0.85 * float(np.min(p_bg))]
    r_bg = bls.power(p_bg, vd if len(vd) else [0.12 * bg_lo], objective="snr")
    pw_bg = np.where(np.asarray(r_bg.depth) > 0, np.asarray(r_bg.power, float), 0.0)
    mp_b = float(np.nanmedian(pw_bg)); mad_b = max(1.4826 * float(np.nanmedian(np.abs(pw_bg - mp_b))), 1e-6)
    p_fn = 1.0 / np.linspace(1.0 / p_hi, 1.0 / p_lo, 12000)[::-1]
    vdf = durs[durs < 0.85 * p_lo]
    r_fn = bls.power(p_fn, vdf if len(vdf) else [0.12 * p_lo], objective="snr")
    pw_f = np.where(np.asarray(r_fn.depth) > 0, np.asarray(r_fn.power, float), 0.0)
    bi = int(np.nanargmax(pw_f)); bp, bt0, bdur = float(r_fn.period[bi]), float(r_fn.transit_time[bi]), float(r_fn.duration[bi])
    st = bls.compute_stats(bp, bdur, bt0)
    ph = ((t - bt0 + 0.5 * bp) % bp) - 0.5 * bp; it, ot = np.abs(ph) <= 0.5 * bdur, np.abs(ph) > 1.25 * bdur
    ntr = int(np.sum(np.asarray(st["per_transit_count"]) >= 2)) or int(np.sum(np.asarray(st["per_transit_count"]) >= 1))
    if np.sum(it) < 4 or np.sum(ot) < 30 or ntr < 2:
        return None
    dep = float(np.median(f[ot]) - np.mean(f[it]))
    if dep <= 0:
        return None
    sc = max(1.4826 * float(np.median(np.abs(f[ot] - np.median(f[ot])))), 1e-6)
    mes = dep / (sc * math.sqrt(1.0 / int(np.sum(it)) + 1.0 / int(np.sum(ot))))
    return {"p": bp, "t0": bt0, "dur": bdur, "dep": dep, "sc": sc, "mes": mes,
            "mes_c": mes / pont_beta(t[ot], f[ot] - np.median(f[ot]), bdur), "sde": float((np.nanmax(pw_f) - mp_b) / mad_b), "ntr": ntr}
