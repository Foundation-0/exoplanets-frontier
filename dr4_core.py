"""Purified Pre-Gaia-DR4 candidate loader, VizieR I/355/gaiadr3 scene audit & TSM ranker."""
from __future__ import annotations
import concurrent.futures as cf, math, re, time
from pathlib import Path
import numpy as np, pandas as pd
from blend_gaia import propagate_scene, query_gaia_cone
from bls_search import RSUN_RE
from campaigns import gv
from converge_io import q_tic


def load_genuine_eumemic_candidates(root: Path = Path("_eumemic"), fallback_df: pd.DataFrame | None = None) -> pd.DataFrame:
    frames = []
    if root.exists():
        cfiles = []
        for fp in sorted(root.rglob("*.csv")):
            rel = str(fp.relative_to(root)).replace("\\", "/")
            rl = rel.lower()
            if any(b in rl for b in ("validation/", "inverted", "false", "reject", "funnel", "sector", "overrides", "products/", "plan")):
                continue
            if "candidates" in fp.name.lower():
                cfiles.append((0 if fp.name.lower() == "candidates_allsky.csv" else 1, rel, fp))
        for _, rel, fp in sorted(cfiles):
            try:
                df = pd.read_csv(fp)
                cl = {c.strip().lower() for c in df.columns}
                if ("tic" in cl or "tic_id" in cl) and ("p_d" in cl or "period_d" in cl or "period" in cl):
                    df["src_file"] = rel; frames.append(df)
            except Exception:
                pass
    if not frames and fallback_df is not None and not fallback_df.empty:
        df0 = fallback_df.copy()
        if "src_file" in df0.columns:
            m = df0["src_file"].astype(str).str.contains("candidates", case=False, na=False) & ~df0["src_file"].astype(str).str.contains("validation|inverted|false|reject", case=False, na=False)
            df0 = df0[m].copy()
        frames.append(df0)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    out["_tid"] = [int(gv(r, "TIC", "tic", "tic_id")) if np.isfinite(gv(r, "TIC", "tic", "tic_id")) else -1 for _, r in out.iterrows()]
    out["_p"] = [round(gv(r, "P_d", "period_d", "period"), 3) for _, r in out.iterrows()]
    return out[(out["_tid"] > 0) & (out["_p"] > 0)].drop_duplicates(subset=["_tid", "_p"], keep="first").reset_index(drop=True)


def load_frontier_flags(out_dir: Path = Path("results")) -> dict[int, str]:
    flags: dict[int, str] = {}
    if (out_dir / "blend_localization.csv").exists():
        for _, r in pd.read_csv(out_dir / "blend_localization.csv").iterrows():
            vd, off = str(r.get("verdict", "")), float(r.get("off_target_as", 0.0))
            if vd in ("NEIGHBOR_BEB_LOCALIZED", "OFFSET_MARGINAL") and off > 12.0:
                flags[int(r["tic"])] = f"wcs_{vd}({off}as)"
            elif vd == "ON_TARGET_CLEAN":
                flags[int(r["tic"])] = f"verified_on_target({off}as)"
    if (out_dir / "hz_issue5_audit.csv").exists():
        for _, r in pd.read_csv(out_dir / "hz_issue5_audit.csv").iterrows():
            vt = str(r.get("vetoes", ""))
            m_oe = re.search(r"odd_even\(([0-9.]+)s\)", vt)
            oe_hi = bool(m_oe and float(m_oe.group(1)) >= 6.5)
            if oe_hi or any(k in vt for k in ("NEIGHBOR_BEB", "OFFSET_MARGINAL", "cent(")):
                flags[int(r["tic"])] = f"multi_sec_flag:{vt[:48]}"
    return flags


def build_665_catalog(dfe: pd.DataFrame | None = None, out_dir: Path = Path("results")):
    cands = load_genuine_eumemic_candidates(Path("_eumemic"), dfe)
    tmap = {int(r["ID"]): r for r in q_tic(sorted(cands["_tid"].unique()))} if not cands.empty else {}
    fflags = load_frontier_flags(out_dir); rows = []
    for _, er in cands.iterrows():
        tid = int(er["_tid"]); tr = tmap.get(tid); p_d = gv(er, "P_d", "period_d", "period", d=-1.0)
        if p_d <= 0: continue
        sfile = str(er.get("src_file", "results_public/candidates_allsky.csv"))
        if any(b in sfile.lower() for b in ("validation", "inverted", "false", "reject")): continue
        utier = str(er.get("tier", er.get("Tier", er.get("status", "candidate")))).strip()
        if any(w in utier.lower() for w in ("rej", "false", "inverted", "eb")): continue
        dppm = gv(er, "depth_ppm", "depth", d=600.0); dppm = dppm * 1e6 if dppm < 1.0 else dppm
        tf = lambda k, d: float(tr[k]) if (tr is not None and k in tr.colnames and tr[k] is not None and np.isfinite(tr[k])) else d
        ra, dec = tf("ra", gv(er, "RA", "ra", d=np.nan)), tf("dec", gv(er, "Dec", "dec", d=np.nan))
        pma, pmd = tf("pmRA", 0.0), tf("pmDEC", 0.0)
        rs, teff, tmag = max(tf("rad", gv(er, "Rstar", "rad", d=0.9)), 0.08), max(tf("Teff", gv(er, "Teff", "teff", d=5300.0)), 2500.0), tf("Tmag", gv(er, "Tmag", "tmag", d=12.5))
        crat = max(tf("contratio", 0.05), 0.0)
        gid = str(tr["GAIA"]) if (tr is not None and "GAIA" in tr.colnames and tr["GAIA"] is not None and str(tr["GAIA"]).strip() not in ("", "--", "None")) else "unmatched"
        lum = str(tr["lumclass"]) if (tr is not None and "lumclass" in tr.colnames and tr["lumclass"] is not None) else "DWARF"
        rp = math.sqrt(max(dppm * 1e-6, 1e-6) * (1.0 + crat)) * rs * RSUN_RE
        a_au = ((max(rs ** 0.9, 0.08) * (p_d / 365.25) ** 2) ** (1.0 / 3.0))
        teq = 278.5 * ((((rs ** 2) * ((teff / 5772.0) ** 4)) / max(a_au ** 2, 1e-6)) ** 0.25) * (0.7 ** 0.25)
        mp, sf = (max(0.15, rp ** 3.5), 0.190) if rp < 1.5 else (2.69 * (rp ** 0.93), 1.26)
        tsm = sf * ((rp ** 3) * teq / max(mp * rs ** 2, 1e-4)) * (10.0 ** (-0.2 * (tmag - (0.75 if teff < 4200 else 0.45))))
        dt = 2026.92 - 2015.5
        ra26 = ra + (pma / (3.6e6 * max(math.cos(math.radians(dec)), 0.05))) * dt if (np.isfinite(ra) and np.isfinite(dec)) else ra
        dc26 = dec + (pmd / 3.6e6) * dt if np.isfinite(dec) else dec
        is_nfp = "fp" in utier.lower() or "nearby" in utier.lower()
        score = (50.0 if 180 <= teq <= 325 and 0.5 <= rp <= 2.8 else 0.0) + min(tsm, 120.0) * 0.5 + (25.0 if crat < 0.05 else 0.0) + (25.0 if utier in ("S", "reported") else (10.0 if utier == "C" else -30.0 if is_nfp else 0.0))
        rows.append({"tic": tid, "src_file": sfile, "gaia_dr3_id": gid, "upstream_tier": utier, "P_d": round(p_d, 6), "depth_ppm": round(dppm, 1),
                     "Rp_Re": round(rp, 2), "Teq_K": round(teq, 0), "is_hz_temperate": bool(180 <= teq <= 325 and rp <= 3.2),
                     "TSM": round(float(tsm), 1), "Tmag": round(tmag, 2), "Teff": round(teff, 0), "Rstar": round(rs, 3),
                     "lumclass": lum, "contratio": round(crat, 4), "ra_2026": round(ra26, 6) if np.isfinite(ra26) else -1.0,
                     "dec_2026": round(dc26, 6) if np.isfinite(dc26) else -1.0, "pm_masyr": round(math.hypot(pma, pmd), 1),
                     "ruwe": -1.0, "nss": -1, "plausible_beb_65as": -1, "frontier_note": fflags.get(tid, "clean"),
                     "priority_score": round(score, 1), "dr4_class": "CATALOG_INDEXED"})
    rows.sort(key=lambda x: x["priority_score"], reverse=True)
    cone_cache: dict[int, list[dict]] = {}; seen_t: set[int] = set()
    u_q = [(r["tic"], r["ra_2026"], r["dec_2026"]) for r in rows if r["ra_2026"] >= 0 and not (r["tic"] in seen_t or seen_t.add(r["tic"]))]
    t0_vz = time.time()
    with cf.ThreadPoolExecutor(max_workers=8) as ex:
        for tk, rc in ex.map(lambda x: (x[0], None if time.time() - t0_vz > 210.0 else query_gaia_cone(x[1], x[2])), u_q):
            if rc is not None: cone_cache[tk] = rc
    for r in rows:
        fn, ut_l = r["frontier_note"], r["upstream_tier"].lower()
        flagged = fn.startswith("multi_sec_flag:") or fn.startswith("wcs_") or (("fp" in ut_l or "nearby" in ut_l) and not fn.startswith("verified_on_target"))
        raw = cone_cache.get(r["tic"])
        if raw is not None and r["ra_2026"] >= 0:
            tg, nbs = propagate_scene(raw, r["ra_2026"], r["dec_2026"], 2026.92, r["depth_ppm"])
            if tg["sid"] != "none" and not tg["sid"].startswith("TIC:"): r["gaia_dr3_id"] = tg["sid"]
            r["ruwe"], r["nss"] = tg["ruwe"], int(tg.get("nss", 0))
            r["plausible_beb_65as"] = sum(1 for nb in nbs if nb["req_eb_depth"] <= 1.0 and nb["sep_as"] <= 55.0)
            bad = (r["ruwe"] > 1.4) or (r["nss"] > 0) or (r["Rstar"] > 2.2) or ("giant" in r["lumclass"].lower()) or (r["Rp_Re"] > 13.5)
            r["dr4_class"] = ("BINARY_OR_EVOLVED_WARNING" if bad else "MULTI_SECTOR_VETTED_WARNING" if flagged else
                              "PLATINUM_SOLO_VERIFIED" if (r["plausible_beb_65as"] == 0 and r["contratio"] < 0.05 and 0.5 <= r["Rp_Re"] <= 13.0) else
                              "GOLD_LOW_BLEND" if (r["plausible_beb_65as"] <= 1 and r["contratio"] < 0.15 and 0.5 <= r["Rp_Re"] <= 13.0) else "SILVER_CROWDED_SG1")
        else:
            bad = (r["Rstar"] > 2.2) or ("giant" in r["lumclass"].lower()) or (r["Rp_Re"] > 13.5)
            r["dr4_class"] = "BINARY_OR_EVOLVED_WARNING" if bad else ("MULTI_SECTOR_VETTED_WARNING" if flagged else ("LOW_CONTAM_INDEXED" if r["contratio"] < 0.05 else "CROWDED_INDEXED"))
    return rows, tmap
