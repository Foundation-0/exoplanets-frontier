"""Asymmetric Huber 4-harmonic starspot pre-whitener, null control, and de-spotting diagnostic plot."""
from __future__ import annotations
import math
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from bls_search import search_bls

def harmonic_despot(t, f, fe, q=None, max_modes: int = 3, n_harm: int = 4, polish_d: float = 1.15):
    t, f, fe = np.asarray(t, float), np.asarray(f, float), np.asarray(fe, float)
    m = np.isfinite(t) & np.isfinite(f) & (f > 0) & (True if q is None else (np.asarray(q, int) == 0))
    t, f, fe = t[m], f[m], fe[m]
    if len(t) < 90: return np.array([]), np.array([]), np.array([]), {"p_rot_d": 0.0, "raw_ppm": 0.0, "clean_ppm": 0.0}
    o = np.argsort(t); t, f, fe = t[o], f[o], fe[o]
    ek = ((t - t[0]) >= 0.10) & ((t[-1] - t) >= 0.10)
    for g in np.where(np.diff(t) > 0.5)[0]:
        ek &= ~(((t <= t[g]) & ((t[g] - t) < 0.10)) | ((t >= t[g + 1]) & ((t - t[g + 1]) < 0.10)))
    if np.sum(ek) >= 90: t, f, fe = t[ek], f[ek], fe[ek]
    med0 = float(np.median(f)); fn = f / med0; raw_ppm = 1.4826 * float(np.median(np.abs(fn - 1.0))) * 1e6
    sp = max(1.4826 * float(np.median(np.abs(np.diff(fn)))), 5e-5)
    fen = np.where(np.isfinite(fe) & (fe > 0), fe / med0, sp)
    to, fo, feo, prots, fgrid = [], [], [], [], np.linspace(0.12, 5.6, 850)
    for s in np.split(np.arange(len(t)), np.where(np.diff(t) > 0.45)[0] + 1):
        if len(s) < 55: continue
        ts, fs, fes = t[s], fn[s].copy(), fen[s].copy()
        mad0 = max(1.4826 * float(np.median(np.abs(fs - np.median(fs)))), 5e-5)
        kp = (fs - np.median(fs)) < 4.2 * mad0; ts, fs, fes = ts[kp], fs[kp], fes[kp]
        if len(ts) < 50: continue
        tau = ts - float(np.mean(ts))
        for _ in range(max_modes):
            y = fs - float(np.median(fs)); mad = max(1.4826 * float(np.median(np.abs(y))), 4e-5); inl = y > -2.2 * mad
            if np.sum(inl) < 40: break
            ang = 2.0 * np.pi * tau[inl, None] * fgrid[None, :]; yc = y[inl] - float(np.mean(y[inl]))
            amp = 2.0 * np.hypot(np.dot(yc, np.cos(ang)), np.dot(yc, np.sin(ang))) / max(np.sum(inl), 1); bi = int(np.argmax(amp))
            if float(amp[bi]) < 5.2 * mad * math.sqrt(2.0 / np.sum(inl)): break
            fine = np.linspace(max(0.10, float(fgrid[bi]) - 0.018), float(fgrid[bi]) + 0.018, 41)
            af = 2.0 * np.pi * tau[inl, None] * fine[None, :]
            f0 = float(fine[int(np.argmax(np.hypot(np.dot(yc, np.cos(af)), np.dot(yc, np.sin(af)))))]); prots.append(1.0 / f0)
            cols = [np.ones(len(tau)), tau] + [fn_ for h in range(1, n_harm + 1) if (h * f0 <= 9.0 or h == 1) for fn_ in (np.cos(2.0 * np.pi * h * f0 * tau), np.sin(2.0 * np.pi * h * f0 * tau))]
            X = np.column_stack(cols); w = np.ones(len(tau), bool); mod = np.ones_like(fs)
            for _it in range(4):
                if np.sum(w) < X.shape[1] + 10: break
                coef, *_ = np.linalg.lstsq(X[w], fs[w], rcond=None); mod = np.clip(X @ coef, 0.35, 1.85); res = fs / mod - 1.0
                rm = max(1.4826 * float(np.median(np.abs(res[w] - np.median(res[w])))), 4e-5); w = (res > -2.2 * rm) & (res < 3.6 * rm)
            fs /= mod; fes /= mod
        wp = max(21, min(int(round(polish_d / max(float(np.median(np.diff(ts))), 1e-4))) | 1, len(ts) - (1 - len(ts) % 2)))
        inl, tr = np.ones(len(fs), bool), np.ones_like(fs)
        for _ in range(2):
            if np.sum(inl) < 20: break
            tr = signal.medfilt(np.interp(ts, ts[inl], fs[inl]), kernel_size=wp); r = fs / tr - 1.0
            rm = max(1.4826 * float(np.median(np.abs(r[inl] - np.median(r[inl])))), 4e-5); inl = (r > -2.3 * rm) & (r < 3.5 * rm)
        to.append(ts); fo.append(fs / tr); feo.append(fes / tr)
    if not to: return np.array([]), np.array([]), np.array([]), {"p_rot_d": 0.0, "raw_ppm": raw_ppm, "clean_ppm": raw_ppm}
    tc, fc, fec = np.concatenate(to), np.concatenate(fo), np.concatenate(feo)
    cln = 1.4826 * float(np.median(np.abs(fc - np.median(fc)))) * 1e6
    return tc, fc, fec, {"p_rot_d": round(float(np.median(prots)), 5) if prots else 0.0, "raw_ppm": round(raw_ppm, 1), "clean_ppm": round(cln, 1)}

def scrambled_null_trial(t, f, fe, sec, rs: float = 0.42, seed: int = 42) -> float:
    rng = np.random.default_rng(seed); fn = np.empty_like(f)
    for s in np.unique(sec):
        m = sec == s; sub = f[m][::-1].copy()
        fn[m] = np.roll(sub, int(rng.integers(5, len(sub) - 5))) if len(sub) > 20 else sub
    pk = search_bls(t, fn, fe, 0.22, 18.0, rs, fast_inv=True)
    return round(float(pk["mes_c"]), 2) if pk else 0.0

def plot_despot_card(out_p: Path, tic: int, tag: str, t_raw, f_raw, t_cln, f_cln, p_d: float, t0_bt: float, dur_d: float):
    out_p.parent.mkdir(parents=True, exist_ok=True)
    m1, m2 = np.isfinite(t_raw) & np.isfinite(f_raw), np.isfinite(t_cln) & np.isfinite(f_cln)
    t_raw, f_raw, t_cln, f_cln = t_raw[m1], f_raw[m1], t_cln[m2], f_cln[m2]
    fig, ax = plt.subplots(1, 3, figsize=(13.0, 3.6), dpi=135); fig.patch.set_facecolor("#0f172a")
    for a in ax: a.set_facecolor("#1e293b"); a.tick_params(colors="#cbd5e1", labelsize=8)
    if len(t_raw) > 20:
        w0 = (t_raw - t_raw[0]) <= 12.0; ax[0].scatter(t_raw[w0] - t_raw[0], (f_raw[w0] / np.median(f_raw[w0]) - 1.0) * 1e6, s=3, alpha=0.45, color="#fb7185")
    if len(t_cln) > 20:
        w1 = (t_cln - t_cln[0]) <= 12.0; ax[1].scatter(t_cln[w1] - t_cln[0], (f_cln[w1] - np.median(f_cln[w1])) * 1e6, s=3, alpha=0.45, color="#38bdf8")
        ph = (((t_cln - t0_bt + 0.5 * p_d) % p_d) - 0.5 * p_d) * 24.0; fp = (f_cln - np.median(f_cln)) * 1e6
        wh = min(max(4.0 * dur_d * 24.0, 2.0), 0.45 * p_d * 24.0); zm = np.abs(ph) <= wh; bins = np.linspace(-wh, wh, 25); bi = np.digitize(ph[zm], bins) - 1
        ax[2].scatter(ph[zm], fp[zm], s=4, alpha=0.22, color="#94a3b8")
        pts = [(0.5 * (bins[b] + bins[b + 1]), float(np.mean(fp[zm][bi == b]))) for b in range(24) if np.sum(bi == b) >= 2]
        if pts: ax[2].plot([x for x, _ in pts], [y for _, y in pts], "o-", color="#34d399", ms=4, lw=1.3)
        ax[2].set_xlim(-wh, wh)
    ax[0].set_title("Naive Median (Spot Wave)", color="#f8fafc", fontsize=8.5)
    ax[1].set_title("Harmonic Pre-Whitened", color="#f8fafc", fontsize=8.5)
    ax[2].set_title(f"Fold P={p_d:.5f}d", color="#f8fafc", fontsize=8.5)
    fig.suptitle(f"TIC {tic}: {tag}", color="#f8fafc", fontsize=9.5, fontweight="bold", y=1.01)
    fig.tight_layout(); fig.savefig(out_p, bbox_inches="tight", facecolor=fig.get_facecolor()); plt.close(fig)
