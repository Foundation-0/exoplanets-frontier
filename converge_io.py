"""MAST spanning light-curve fetchers, TOI-540 benchmark runner, M-dwarf despot worker & report writer."""
from __future__ import annotations
import concurrent.futures as cf, re
from pathlib import Path
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

def q_tic(ids):
    try: return list(Catalogs.query_criteria(catalog="Tic", ID=list(ids)) or [])
    except Exception: return []

def fetch_raw_spoc(tid: int, max_sec: int = 10):
    obs = Observations.query_criteria(obs_collection="TESS", dataproduct_type="timeseries", provenance_name="SPOC", target_name=str(tid))
    if obs is None or len(obs) == 0: obs = Observations.query_criteria(obs_collection="TESS", provenance_name="SPOC", obs_id=f"*{tid:016d}*")
    if obs is None or len(obs) == 0: return []
    pr = Observations.get_product_list(obs); seen, u0 = set(), []
    for r in sorted([x for x in pr if str(x["productFilename"]).endswith(("_s_lc.fits", "-s_lc.fits"))], key=lambda x: str(x["productFilename"])):
        m = re.search(r"-s(\d{4})-", str(r["productFilename"])); s = int(m.group(1)) if m else 0
        if s and s != 1751 and s not in seen: seen.add(s); u0.append((s, str(r["dataURI"])))
    ix = sorted({0, 1, len(u0) // 3, (2 * len(u0)) // 3, len(u0) - 2, len(u0) - 1} | set(range(max(0, len(u0) - 4), len(u0)))) if len(u0) > max_sec else list(range(len(u0)))
    with cf.ThreadPoolExecutor(max_workers=4) as ex: return [x for x in ex.map(lambda su: dl_lc(su[1], su[0]), [u0[i] for i in ix[:max_sec]]) if x]

def fetch_spanning_lcs(tid: int, ra: float, dec: float, win_d: float = 1.4):
    raw = fetch_raw_spoc(tid, max_sec=8); lcs, have = [], {x["sec"] for x in raw}
    for l in raw:
        tc, fc, fec = detrend_lc(l["t"], l["f"], l["fe"], l["q"], win_d)
        if len(tc) >= 60: lcs.append({"sec": l["sec"], "t": tc, "f": fc, "fe": fec, "mode": "SPOC"})
    try:
        st = Tesscut.get_sectors(coordinates=SkyCoord(ra=ra * u.deg, dec=dec * u.deg))
        c = next((k for k in st.colnames if "sec" in k.lower()), None) if st is not None and len(st) else None
        alls = sorted({int(r[c]) for r in st}) if c else []
    except Exception: alls = [106, 107]
    want = [s for s in sorted(set(alls[:2] + alls[len(alls) // 2 - 1 : len(alls) // 2 + 1] + alls[-2:])) if s not in have][:6] if alls else []
    def _f(s):
        r = ffi_cutout_lc(tid, ra, dec, s)
        if not r: return None
        tc, fc, fec = detrend_lc(r["t"], r["f"], r["fe"], r["q"], win_d)
        return {"sec": s, "t": tc, "f": fc, "fe": fec, "mode": "FFI"} if len(tc) >= 60 else None
    if want:
        with cf.ThreadPoolExecutor(max_workers=3) as ex: lcs.extend([x for x in ex.map(_f, want) if x])
    return sorted(lcs, key=lambda x: x["sec"]), alls

def run_benchmark_and_null(fig_d: Path):
    raw = fetch_raw_spoc(200322593, max_sec=4); tn, fn, fen, td, fd, fed, sd, infs = [], [], [], [], [], [], [], []
    for l in sorted(raw, key=lambda x: x["sec"]):
        a, b, c = detrend_lc(l["t"], l["f"], l["fe"], l["q"], 1.0)
        if len(a) >= 80: tn.append(a); fn.append(b); fen.append(c)
        a2, b2, c2, inf = harmonic_despot(l["t"], l["f"], l["fe"], l["q"])
        if len(a2) >= 80: td.append(a2); fd.append(b2); fed.append(c2); sd.append(np.full(len(a2), l["sec"], int)); infs.append(inf)
    pkn = search_bls(np.concatenate(tn), np.concatenate(fn), np.concatenate(fen), 0.25, 4.0, 0.19) if tn else None
    tc, fc, fec, sc = (np.concatenate(td), np.concatenate(fd), np.concatenate(fed), np.concatenate(sd)) if td else (np.array([]),) * 4
    pkd = search_bls(tc, fc, fec, 0.25, 4.0, 0.19) if len(tc) else None
    if pkd and tn:
        plot_despot_card(fig_d / "benchmark_TOI540_despot.png", 200322593, f"TOI-540b (Naive P={(pkn['p'] if pkn else -1.0):.4f}d -> Clean P={pkd['p']:.5f}d, MES={pkd['mes_c']:.1f})", tn[0], fn[0], tc, fc, pkd["p"], pkd["t0"], pkd["dur"])
    nm = [scrambled_null_trial(tc, fc, fec, sc, 0.19, 100 + k) for k in range(10)] if len(tc) else []
    return {"benchmark_tic": 200322593, "published_P_d": 1.239149, "naive_P_d": round(pkn["p"], 6) if pkn else -1.0,
            "naive_locked_on_spot": bool(pkn and abs(pkn["p"] - 1.239149) > 0.02), "despotted_P_d": round(pkd["p"], 6) if pkd else -1.0,
            "despotted_MES": round(pkd["mes_c"], 2) if pkd else 0.0, "despotted_depth_ppm": round(pkd["dep"] * 1e6, 1) if pkd else 0.0,
            "recovered_within_0p1pct": bool(pkd and abs(pkd["p"] - 1.239149) / 1.239149 < 0.001),
            "mean_raw_rms_ppm": round(float(np.mean([x["raw_ppm"] for x in infs])), 1) if infs else 0.0,
            "mean_clean_rms_ppm": round(float(np.mean([x["clean_ppm"] for x in infs])), 1) if infs else 0.0,
            "null_trials": len(nm), "null_max_MES": max(nm) if nm else 0.0, "null_false_alarms_ge_9p5": sum(1 for m in nm if m >= 9.5)}

def despot_one_mdwarf(item):
    tid, rs0 = item["tic"], item["rstar"]; raw = fetch_raw_spoc(tid, max_sec=10)
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
    if not pk or pk["mes_c"] < 7.5:
        return {"tic": tid, "lane": item["lane"], "p_rot_d": round(prot, 5), "raw_ppm": r0, "clean_ppm": c0, "P_d": round(pk["p"], 6) if pk else -1.0, "MES": round(pk["mes_c"], 2) if pk else 0.0, "passed": False, "vetoes": "sub_threshold"}, None
    gm, nbs = gaia_nb(tid, hdr["ra"], hdr["dec"], hdr)
    rep = vet_signal(tid, pk, t, f, fe, sc, item["lane"], norm_sap(np.concatenate(sapa), sc), np.concatenate(mxa), np.concatenate(mya), float(np.median(cr)), gm, nbs)
    amp = spot_amp(t, f, pk["p"], pk["t0"], pk["dur"])
    if amp > 0.48 * pk["dep"] or (prot > 0 and any(abs(pk["p"] - k * prot) / (k * prot) < 0.02 for k in (0.5, 1.0, 2.0))):
        rep["passed"] = False; rep["tier"] = "rejected"; rep["vetoes"] = ("" if rep["vetoes"] == "CLEAN" else rep["vetoes"] + ";") + f"resid_rot({amp/max(pk['dep'],1e-6):.2f})"
    rep["p_rot_d"], rep["raw_ppm"], rep["clean_ppm"] = round(prot, 5), r0, c0
    return rep, {"raw_t": raw[0]["t"], "raw_f": raw[0]["f"], "t": t, "f": f, "dur": pk["dur"], "t0": pk["t0"]}
