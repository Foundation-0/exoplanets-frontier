"""Detrending and two-stage multi-year AstroPy BoxLeastSquares transit search."""
from __future__ import annotations
import math
import numpy as np
from scipy import signal

BTJD0, PIX_AS, RSUN_RE = 2457000.0, 21.0, 109.076

def detrend_lc(t, f, fe, q=None, win_d=1.2):
    t, f, fe = np.asarray(t, float), np.asarray(f, float), np.asarray(fe, float)
    m = np.isfinite(t) & np.isfinite(f) & (f > 0) & (True if q is None else (np.asarray(q, int) == 0))
    t, f, fe = t[m], f[m], fe[m]
    if len(t) < 80: return np.array([]), np.array([]), np.array([])
    o = np.argsort(t); t, f, fe = t[o], f[o], fe[o]
    ek = ((t - t[0]) >= 0.12) & ((t[-1] - t) >= 0.12)
    for g in np.where(np.diff(t) > 0.5)[0]:
        ek &= ~(((t <= t[g]) & ((t[g] - t) < 0.12)) | ((t >= t[g + 1]) & ((t - t[g + 1]) < 0.12)))
    if np.sum(ek) >= 80: t, f, fe = t[ek], f[ek], fe[ek]
    fn = f / float(np.median(f))
    sp = max(1.4826 * float(np.median(np.abs(np.diff(fn)))), 5e-5)
    fen = np.where(np.isfinite(fe) & (fe > 0), fe / float(np.median(f)), sp)
    to, fo, feo = [], [], []
    for s in np.split(np.arange(len(t)), np.where(np.diff(t) > 0.4)[0] + 1):
        if len(s) < 40: continue
        ts, fs, fes = t[s], fn[s], fen[s]
        wp = max(21, min(int(round(win_d / max(float(np.median(np.diff(ts))), 1e-4))) | 1, len(ts) - (1 - len(ts) % 2)))
        inl, tr = np.ones(len(fs), bool), np.ones_like(fs)
        for _ in range(3):
            if np.sum(inl) < 15: break
            tr = signal.medfilt(np.interp(ts, ts[inl], fs[inl]), kernel_size=wp)
            r = fs / tr - 1.0; mad = max(1.4826 * float(np.median(np.abs(r[inl] - np.median(r[inl])))), 5e-5)
            inl = (r < 3.5 * mad) & (r > -2.6 * mad)
        det = fs / tr; rf = det - 1.0
        kp = rf < 4.5 * max(1.4826 * float(np.median(np.abs(rf - np.median(rf)))), 5e-5)
        to.append(ts[kp]); fo.append(det[kp]); feo.append(fes[kp] / tr[kp])
    return (np.concatenate(to), np.concatenate(fo), np.concatenate(feo)) if to else (np.array([]), np.array([]), np.array([]))

def exp_dur(p, rs=0.7):
    a = 215.032 * ((max(rs ** 0.9, 0.08) * (p / 365.25) ** 2) ** (1.0 / 3.0))
    return (p / math.pi) * math.asin(min(max(rs, 0.08) / max(a, 1e-3), 0.65))

def pont_beta(t, res, dur):
    s1 = 1.4826 * float(np.median(np.abs(res - np.median(res))))
    ed = np.arange(float(t[0]), float(t[-1]), max(dur, 0.01))
    if s1 <= 0 or len(ed) < 12: return 1.0
    ix = np.digitize(t, ed) - 1; v = (ix >= 0) & (ix < len(ed) - 1)
    c = np.bincount(ix[v], minlength=len(ed) - 1); sm = np.bincount(ix[v], weights=res[v], minlength=len(ed) - 1)
    g = c >= 3
    if np.sum(g) < 10: return 1.0
    bm = sm[g] / c[g]
    return float(np.clip(1.4826 * float(np.median(np.abs(bm - np.median(bm)))) / max(s1 / math.sqrt(float(np.median(c[g]))), 1e-9), 1.0, 5.0))

def search_bls(t, f, fe, p_min=0.16, p_max=25.0, rs=0.65, seed_p=None, fast_inv=False):
    from astropy.timeseries import BoxLeastSquares
    fb = float(t[-1] - t[0]) if len(t) >= 80 else 0.0
    p_hi = min(p_max, fb / 2.4)
    if fb <= p_min * 1.4 or p_hi <= p_min * 1.02: return None
    d0 = max(0.012, 0.35 * exp_dur(p_min, rs))
    durs = np.geomspace(d0, max(d0 * 1.15, min(0.42 * p_min, 1.65 * exp_dur(p_hi, rs))), 7)
    bw = (t - t[0]) <= 60.0
    for at in t[::max(1, len(t) // 20)]:
        w = (t >= at) & (t <= at + 60.0)
        if np.sum(w) > np.sum(bw): bw = w
    ts1, fs1, fes1 = (t[bw], f[bw], fe[bw]) if (np.sum(bw) >= 120 and fb > 70.0) else (t, f, fe)
    bs1 = max(float(ts1[-1] - ts1[0]), 18.0); seeds, sde = [], 10.0
    if seed_p and p_min <= seed_p <= p_hi:
        seeds.append((float(seed_p), min(fb, 27.0)))
    else:
        phi1 = min(p_hi, max(p_min * 1.15, bs1 * 0.90)); mx = 1500 if fast_inv else 22000
        b1 = BoxLeastSquares(ts1, fs1, dy=fes1)
        ps1 = b1.autoperiod(durs, minimum_period=p_min, maximum_period=phi1, frequency_factor=3.0)
        if len(ps1) > mx: ps1 = 1.0 / np.linspace(1.0 / phi1, 1.0 / p_min, mx)[::-1]
        vd = durs[durs < 0.85 * float(np.min(ps1))]
        r1 = b1.power(ps1, vd if len(vd) else [0.15 * float(np.min(ps1))], objective="snr")
        pw = np.where(np.asarray(r1.depth) > 0, np.asarray(r1.power, float), 0.0)
        mp = float(np.nanmedian(pw)); sde = float((np.nanmax(pw) - mp) / max(1.4826 * float(np.nanmedian(np.abs(pw - mp))), 1e-6))
        for ix in np.argsort(pw)[::-1]:
            cp = float(r1.period[ix])
            if all(abs(cp - c0) / c0 > 0.004 for c0, _ in seeds):
                seeds.append((cp, bs1))
                if len(seeds) >= (3 if fast_inv else 9): break
    if not seeds: return None
    comb = [np.linspace(max(1.0 / p_hi, 1.0 / cp - max(2.0 / bsd, 0.008 / cp if seed_p else 0.0)),
                        min(1.0 / p_min, 1.0 / cp + max(2.0 / bsd, 0.008 / cp if seed_p else 0.0)),
                        400 if fast_inv else 2000) for cp, bsd in seeds]
    allp = 1.0 / np.unique(np.concatenate(comb))[::-1]
    bf = BoxLeastSquares(t, f, dy=fe); vd = durs[durs < 0.85 * float(np.min(allp))]
    rf = bf.power(allp, vd if len(vd) else [0.15 * float(np.min(allp))], objective="snr")
    pwf = np.where(np.asarray(rf.depth) > 0, np.asarray(rf.power, float), 0.0); bi = int(np.nanargmax(pwf))
    bp, bt0, bdur = float(rf.period[bi]), float(rf.transit_time[bi]), float(rf.duration[bi])
    if seed_p:
        mpf = float(np.nanmedian(pwf)); sde = max(sde, float((np.nanmax(pwf) - mpf) / max(1.4826 * float(np.nanmedian(np.abs(pwf - mpf))), 1e-6)))
    st = bf.compute_stats(bp, bdur, bt0)
    ph = ((t - bt0 + 0.5 * bp) % bp) - 0.5 * bp; it, ot = np.abs(ph) <= 0.5 * bdur, np.abs(ph) > 1.2 * bdur
    ntr = int(np.sum(np.asarray(st["per_transit_count"]) >= 2)) or int(np.sum(np.asarray(st["per_transit_count"]) >= 1))
    if np.sum(it) < 5 or np.sum(ot) < 30 or ntr < 3: return None
    dep = float(np.median(f[ot]) - np.mean(f[it]))
    if dep <= 0: return None
    sc = max(1.4826 * float(np.median(np.abs(f[ot] - np.median(f[ot])))), 1e-6)
    mes = dep / (sc * math.sqrt(1.0 / int(np.sum(it)) + 1.0 / int(np.sum(ot))))
    beta = pont_beta(t[ot], f[ot] - np.median(f[ot]), bdur)
    return {"p": bp, "t0": bt0, "dur": bdur, "dep": dep, "sc": sc, "mes": mes, "mes_c": mes / beta, "sde": sde, "ntr": ntr, "stats": st}
