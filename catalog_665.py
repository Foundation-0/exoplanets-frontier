"""665-candidate Pre-Gaia-DR4 readiness, JWST TSM, temperate HZ & VizieR blend cataloger."""
from __future__ import annotations
import concurrent.futures as cf, math
from pathlib import Path
import numpy as np, pandas as pd
from blend_gaia import propagate_scene, query_gaia_cone
from bls_search import RSUN_RE, detrend_lc, search_bls
from campaigns import gv
from converge_io import fetch_raw_spoc, q_tic
from despot import harmonic_despot, plot_despot_card, scrambled_null_trial

def run_benchmark(fig_d: Path):
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
            "despotted_MES": round(pkd["mes_c"], 2) if pkd else 0.0, "recovered_within_0p1pct": bool(pkd and abs(pkd["p"] - 1.239149) / 1.239149 < 0.001),
            "null_trials": len(nm), "null_max_MES": max(nm) if nm else 0.0, "null_false_alarms_ge_9p5": sum(1 for m in nm if m >= 9.5)}

def build_665_catalog(dfe: pd.DataFrame):
    cands = dfe[dfe["src_file"].str.contains("candidates", case=False, na=False)].copy() if "src_file" in dfe.columns else dfe.copy()
    if len(cands) < 100: cands = dfe.copy()
    cands = cands.drop_duplicates(subset=["_tid", "_p"]).reset_index(drop=True)
    tmap = {int(r["ID"]): r for r in q_tic(sorted(cands["_tid"].unique()))}
    rows = []
    for _, er in cands.iterrows():
        tid = int(er["_tid"]); tr = tmap.get(tid); p_d = gv(er, "P_d", "period_d", "period", d=-1.0)
        if p_d <= 0: continue
        dppm = gv(er, "depth_ppm", "depth", d=600.0); dppm = dppm * 1e6 if dppm < 1.0 else dppm
        utier = str(er.get("tier", er.get("Tier", er.get("status", "v3")))).strip()
        if any(w in utier.lower() for w in ("rej", "false", "eb")): continue
        ra = float(tr["ra"]) if tr is not None and tr["ra"] is not None and np.isfinite(tr["ra"]) else gv(er, "RA", "ra", d=np.nan)
        dec = float(tr["dec"]) if tr is not None and tr["dec"] is not None and np.isfinite(tr["dec"]) else gv(er, "Dec", "dec", d=np.nan)
        pma = float(tr["pmRA"]) if tr is not None and "pmRA" in tr.colnames and tr["pmRA"] is not None and np.isfinite(tr["pmRA"]) else 0.0
        pmd = float(tr["pmDEC"]) if tr is not None and "pmDEC" in tr.colnames and tr["pmDEC"] is not None and np.isfinite(tr["pmDEC"]) else 0.0
        rs = float(tr["rad"]) if tr is not None and tr["rad"] is not None and np.isfinite(tr["rad"]) and float(tr["rad"]) > 0 else gv(er, "Rstar", "rad", d=0.9)
        teff = float(tr["Teff"]) if tr is not None and tr["Teff"] is not None and np.isfinite(tr["Teff"]) and float(tr["Teff"]) > 0 else gv(er, "Teff", "teff", d=5300.0)
        tmag = float(tr["Tmag"]) if tr is not None and tr["Tmag"] is not None and np.isfinite(tr["Tmag"]) else gv(er, "Tmag", "tmag", d=12.5)
        crat = float(tr["contratio"]) if tr is not None and "contratio" in tr.colnames and tr["contratio"] is not None and np.isfinite(tr["contratio"]) else 0.05
        gid = str(tr["GAIA"]) if tr is not None and "GAIA" in tr.colnames and tr["GAIA"] is not None and str(tr["GAIA"]).strip() not in ("", "--", "None") else "unmatched"
        lum = str(tr["lumclass"]) if tr is not None and "lumclass" in tr.colnames and tr["lumclass"] is not None else "DWARF"
        rp = math.sqrt(max(dppm * 1e-6, 1e-6) * (1.0 + max(crat, 0.0))) * rs * RSUN_RE
        a_au = ((max(rs ** 0.9, 0.08) * (p_d / 365.25) ** 2) ** (1.0 / 3.0))
        insol = ((rs ** 2) * ((teff / 5772.0) ** 4)) / max(a_au ** 2, 1e-6)
        teq = 278.5 * (insol ** 0.25) * (0.7 ** 0.25)
        mp, sf = (max(0.15, rp ** 3.5), 0.190) if rp < 1.5 else (2.69 * (rp ** 0.93), 1.26)
        tsm = sf * ((rp ** 3) * teq / max(mp * rs ** 2, 1e-4)) * (10.0 ** (-0.2 * (tmag - (0.75 if teff < 4200 else 0.45))))
        dt = 2026.92 - 2015.5
        ra26 = ra + (pma / (3.6e6 * max(math.cos(math.radians(dec)), 0.05))) * dt if np.isfinite(ra) and np.isfinite(dec) else ra
        dc26 = dec + (pmd / 3.6e6) * dt if np.isfinite(dec) else dec
        score = (50.0 if 180 <= teq <= 325 and 0.5 <= rp <= 2.8 else 0.0) + min(tsm, 120.0) * 0.5 + (25.0 if crat < 0.05 else 0.0) + (20.0 if "s" in utier.lower() or "rep" in utier.lower() else 0.0)
        rows.append({"tic": tid, "gaia_dr3_id": gid, "upstream_tier": utier, "P_d": round(p_d, 6), "depth_ppm": round(dppm, 1),
                     "Rp_Re": round(rp, 2), "Teq_K": round(teq, 0), "is_hz_temperate": bool(180 <= teq <= 325 and rp <= 3.2),
                     "TSM": round(float(tsm), 1), "Tmag": round(tmag, 2), "Teff": round(teff, 0), "Rstar": round(rs, 3),
                     "lumclass": lum, "contratio": round(crat, 4), "ra_2026": round(ra26, 6) if np.isfinite(ra26) else -1.0,
                     "dec_2026": round(dc26, 6) if np.isfinite(dc26) else -1.0, "pm_masyr": round(math.hypot(pma, pmd), 1),
                     "ruwe": -1.0, "nss": -1, "plausible_beb_65as": -1, "priority_score": round(score, 1), "dr4_class": "CATALOG_INDEXED"})
    rows.sort(key=lambda x: x["priority_score"], reverse=True)
    def _vz(r):
        if r["ra_2026"] < 0: return r
        raw = query_gaia_cone(r["ra_2026"], r["dec_2026"]); tg, nbs = propagate_scene(raw, r["ra_2026"], r["dec_2026"], 2026.92, r["depth_ppm"])
        if tg["sid"] != "none" and not tg["sid"].startswith("TIC:"): r["gaia_dr3_id"] = tg["sid"]
        r["ruwe"], r["nss"] = tg["ruwe"], int(tg.get("nss", 0))
        r["plausible_beb_65as"] = sum(1 for nb in nbs if nb["req_eb_depth"] <= 1.0 and nb["sep_as"] <= 55.0)
        bad = (r["ruwe"] > 1.4) or (r["nss"] > 0) or (r["Rstar"] > 2.2) or ("giant" in r["lumclass"].lower()) or (r["Rp_Re"] > 9.5)
        r["dr4_class"] = ("BINARY_OR_EVOLVED_WARNING" if bad else
                          "PLATINUM_SOLO_VERIFIED" if (r["plausible_beb_65as"] == 0 and r["contratio"] < 0.05 and 0.5 <= r["Rp_Re"] <= 4.5) else
                          "GOLD_LOW_BLEND" if (r["plausible_beb_65as"] <= 1 and r["contratio"] < 0.15) else "SILVER_CROWDED_SG1")
        return r
    with cf.ThreadPoolExecutor(max_workers=5) as ex: top_v = list(ex.map(_vz, rows[:95]))
    for r in rows[95:]:
        bad = (r["Rstar"] > 2.2) or ("giant" in r["lumclass"].lower()) or (r["Rp_Re"] > 9.5)
        r["dr4_class"] = "BINARY_OR_EVOLVED_WARNING" if bad else ("LOW_CONTAM_INDEXED" if r["contratio"] < 0.05 else "CROWDED_INDEXED")
    return top_v + rows[95:], tmap
