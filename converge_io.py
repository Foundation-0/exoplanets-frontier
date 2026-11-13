"""MAST Cycle 1-8 sector resolver, fast SPOC/FFI light-curve fetcher & M-dwarf worker."""
from __future__ import annotations
import concurrent.futures as cf, math, re
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Catalogs, Observations, Tesscut
from bls_search import detrend_lc, search_bls
from deep_search import ffi_cutout_lc
from despot import harmonic_despot, plot_despot_card, scrambled_null_trial
from mast_data import dl_lc, gaia_nb
from mdwarf_hunt import norm_sap, spot_amp
from vetter import vet_signal

Observations.TIMEOUT = 45
Catalogs.TIMEOUT = 45

def q_tic(ids):
    out, ids_l = [], list(ids)
    for i in range(0, len(ids_l), 150):
        try: out.extend(list(Catalogs.query_criteria(catalog="Tic", ID=ids_l[i : i + 150]) or []))
        except Exception: pass
    return out

def discover_all_sectors(ra: float, dec: float) -> list[int]:
    secs = set()
    try:
        st = Tesscut.get_sectors(coordinates=SkyCoord(ra=ra * u.deg, dec=dec * u.deg), radius=0.02 * u.deg)
        c = next((k for k in st.colnames if "sec" in k.lower()), None) if st is not None and len(st) else None
        if c: secs.update(int(r[c]) for r in st if int(r[c]) < 500)
    except Exception: pass
    if len(secs) < 3:
        dra = 1.15 / max(math.cos(math.radians(dec)), 0.08)
        try:
            ob = Observations.query_criteria(obs_collection="TESS", provenance_name="SPOC", dataproduct_type="timeseries",
                                             s_ra=[max(0.0, ra - dra), min(360.0, ra + dra)], s_dec=[max(-90.0, dec - 1.15), min(90.0, dec + 1.15)])
            for r in (ob or []):
                m = re.search(r"-s(\d{4})-", str(r["obs_id"]))
                if m and int(m.group(1)) < 500: secs.add(int(m.group(1)))
        except Exception: pass
    return sorted(secs) or [106, 107]

def fetch_raw_spoc(tid: int, max_sec: int = 10):
    try:
        obs = Observations.query_criteria(obs_collection="TESS", dataproduct_type="timeseries", provenance_name="SPOC", target_name=str(tid))
        if obs is None or len(obs) == 0: return []
        seen, u0 = set(), []
        for r in sorted(obs, key=lambda x: str(x["obs_id"])):
            oid = str(r["obs_id"])
            if oid.endswith("_fast") or not re.search(r"-s\d{4}-\d+-\d+-s$", oid): continue
            m = re.search(r"-s(\d{4})-", oid); s = int(m.group(1)) if m else 0
            if s and s != 1751 and s not in seen: seen.add(s); u0.append((s, f"mast:TESS/product/{oid}_lc.fits"))
        if not u0: return []
        ix = sorted({0, 1, len(u0) // 3, (2 * len(u0)) // 3, len(u0) - 2, len(u0) - 1} | set(range(max(0, len(u0) - 4), len(u0)))) if len(u0) > max_sec else list(range(len(u0)))
        with cf.ThreadPoolExecutor(max_workers=4) as ex: return [x for x in ex.map(lambda su: dl_lc(su[1], su[0]), [u0[i] for i in ix[:max_sec]]) if x]
    except Exception: return []

def fetch_spanning_lcs(tid: int, ra: float, dec: float, win_d: float = 1.4):
    raw = fetch_raw_spoc(tid, max_sec=8); lcs, have = [], {x["sec"] for x in raw}
    for l in raw:
        tc, fc, fec = detrend_lc(l["t"], l["f"], l["fe"], l["q"], win_d)
        if len(tc) >= 60: lcs.append({"sec": l["sec"], "t": tc, "f": fc, "fe": fec, "mode": "SPOC"})
    alls = discover_all_sectors(ra, dec); n = len(alls)
    picks = sorted({alls[i] for i in [0, min(1, n - 1), n // 4, n // 2, (3 * n) // 4, max(0, n - 2), n - 1]}) if n else []
    want = [s for s in picks if s not in have][:6]
    def _f(s):
        r = ffi_cutout_lc(tid, ra, dec, s)
        if not r: return None
        tc, fc, fec = detrend_lc(r["t"], r["f"], r["fe"], r["q"], win_d)
        return {"sec": s, "t": tc, "f": fc, "fe": fec, "mode": "FFI"} if len(tc) >= 60 else None
    if want:
        with cf.ThreadPoolExecutor(max_workers=3) as ex: lcs.extend([x for x in ex.map(_f, want) if x])
    return sorted(lcs, key=lambda x: x["sec"]), alls

def despot_one_mdwarf(item):
    try:
        tid, rs0 = item["tic"], item["rstar"]; raw = fetch_raw_spoc(tid, max_sec=8)
        if len(raw) < 2: return None
        ta, fa, fea, sa, sapa, mxa, mya, cr, infs = [], [], [], [], [], [], [], [], []
        for l in sorted(raw, key=lambda x: x["sec"]):
            tc, fc, fec, inf = harmonic_despot(l["t"], l["f"], l["fe"], l["q"])
            if len(tc) < 60: continue
            im = np.clip(np.searchsorted(l["t"], tc), 0, len(l["t"]) - 1)
            ta.append(tc); fa.append(fc); fea.append(fec); sa.append(np.full(len(tc), l["sec"], int)); sapa.append(l["sap"][im]); mxa.append(l["mdx"][im]); mya.append(l["mdy"][im]); cr.append(l["crowd"]); infs.append(inf)
        if len(ta) < 2: return None
        t, f, fe, sc = np.concatenate(ta), np.concatenate(fa), np.concatenate(fea), np.concatenate(sa); hdr = raw[-1]
        rs = float(hdr["rstar"]) if np.isfinite(hdr["rstar"]) and hdr["rstar"] > 0 else rs0; pk = search_bls(t, f, fe, 0.20, 22.0, rs)
        prot = float(np.median([x["p_rot_d"] for x in infs if x["p_rot_d"] > 0] or [0.0]))
        r0, c0 = round(float(np.mean([x["raw_ppm"] for x in infs])), 1), round(float(np.mean([x["clean_ppm"] for x in infs])), 1)
        if not pk or pk["mes_c"] < 8.5:
            return {"tic": tid, "lane": item["lane"], "p_rot_d": round(prot, 5), "raw_ppm": r0, "clean_ppm": c0, "P_d": round(pk["p"], 6) if pk else -1.0, "MES": round(pk["mes_c"], 2) if pk else 0.0, "Rp_Re": -1.0, "Teq_K": -1.0, "passed": False, "vetoes": "sub_threshold"}, None
        gm, nbs = gaia_nb(tid, hdr["ra"], hdr["dec"], hdr)
        rep = vet_signal(tid, pk, t, f, fe, sc, item["lane"], norm_sap(np.concatenate(sapa), sc), np.concatenate(mxa), np.concatenate(mya), float(np.median(cr)), gm, nbs)
        amp = spot_amp(t, f, pk["p"], pk["t0"], pk["dur"])
        if amp > 0.48 * pk["dep"] or (prot > 0 and any(abs(pk["p"] - k * prot) / (k * prot) < 0.02 for k in (0.5, 1.0, 2.0))):
            rep["passed"] = False; rep["tier"] = "rejected"; rep["vetoes"] = ("" if rep["vetoes"] == "CLEAN" else rep["vetoes"] + ";") + f"resid_rot({amp/max(pk['dep'],1e-6):.2f})"
        rep["p_rot_d"], rep["raw_ppm"], rep["clean_ppm"] = round(prot, 5), r0, c0
        return rep, {"raw_t": raw[0]["t"], "raw_f": raw[0]["f"], "t": t, "f": f, "dur": pk["dur"], "t0": pk["t0"]}
    except Exception: return None

def plot_ttv_and_blend(out_p: Path, tic: int, sec_r: list[dict], fit: dict, t, f, p_ref: float, t0_bt: float, dur_d: float, snap=None):
    out_p.parent.mkdir(parents=True, exist_ok=True); m = np.isfinite(t) & np.isfinite(f); t, f = t[m], f[m]; nc = 3 if snap else 2
    fig, ax = plt.subplots(1, nc, figsize=(4.5 * nc, 3.7), dpi=135); fig.patch.set_facecolor("#0f172a")
    for a in ax: a.set_facecolor("#1e293b"); a.tick_params(colors="#cbd5e1", labelsize=8)
    if sec_r:
        E, dt, er = np.array([r["epoch"] for r in sec_r], float), np.array([r["dt_min"] for r in sec_r], float), np.array([r["err_min"] for r in sec_r], float)
        ax[0].errorbar(E, dt, yerr=er, fmt="o", color="#38bdf8", ms=5, lw=1.2)
        eg = np.linspace(float(np.min(E)), float(np.max(E)), 100); sl = fit["dP_sec"] / 60.0
        ax[0].plot(eg, float(np.mean(dt - sl * E)) + sl * eg, "--", color="#34d399", lw=1.3, label=f"dP={fit['dP_sec']:+.2f}s")
        ax[0].axhline(0.0, color="#94a3b8", ls=":", lw=0.9); ax[0].legend(fontsize=7, facecolor="#0f172a", labelcolor="#f8fafc")
    ax[0].set_title(f"S1-S108 O-C ({fit['class']})", color="#f8fafc", fontsize=8.5)
    if len(t) >= 40:
        ph = (((t - t0_bt + 0.5 * p_ref) % p_ref) - 0.5 * p_ref) * 24.0; fp = (f - np.median(f)) * 1e6
        wh = min(max(3.5 * dur_d * 24.0, 2.4), 0.45 * p_ref * 24.0); zm = np.abs(ph) <= wh; bins = np.linspace(-wh, wh, 25); bi = np.digitize(ph[zm], bins) - 1
        ax[1].scatter(ph[zm], fp[zm], s=4, alpha=0.22, color="#94a3b8")
        pts = [(0.5 * (bins[b] + bins[b + 1]), float(np.mean(fp[zm][bi == b]))) for b in range(24) if np.sum(bi == b) >= 2]
        if pts: ax[1].plot([x for x, _ in pts], [y for _, y in pts], "o-", color="#fbbf24", ms=4, lw=1.3)
        ax[1].set_xlim(-wh, wh)
    ax[1].set_title(f"Refined Fold P={p_ref:.6f}d", color="#f8fafc", fontsize=8.5)
    if snap:
        ax[2].imshow(snap["diff"], origin="lower", cmap="inferno"); ax[2].plot(*snap["txy"], "+", color="#38bdf8", ms=11, mew=2); ax[2].plot(*snap["dxy"], "x", color="#34d399", ms=9, mew=2)
        for nx, ny, okb in snap["nbs"]: ax[2].plot(nx, ny, "o" if okb else ".", color="#f43f5e" if okb else "#cbd5e1", ms=5 if okb else 3)
        ax[2].set_title(f"S{snap['sec']} WCS Diff (SNR={snap['snr']:.1f})", color="#f8fafc", fontsize=8.5)
    fig.suptitle(f"TIC {tic}: Ephemeris & Gaia WCS Audit", color="#f8fafc", fontsize=9.5, fontweight="bold", y=1.01)
    fig.tight_layout(); fig.savefig(out_p, bbox_inches="tight", facecolor=fig.get_facecolor()); plt.close(fig)
