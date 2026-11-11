"""MAST TESS FITS/TPF/DVT retrieval, archive cross-matching, and Gaia DR3 queries."""
from __future__ import annotations
import concurrent.futures as cf, math, os, re, tempfile, threading
from pathlib import Path
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import vstack
import astropy.units as u
from astroquery.gaia import Gaia
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astroquery.mast import Catalogs, Observations
from astroquery.vizier import Vizier
from bls_search import PIX_AS, RSUN_RE, detrend_lc

_TMP, _LK, _CNT = Path(tempfile.mkdtemp(prefix="tess_")), threading.Lock(), 0

def dl_uri(uri: str) -> Path | None:
    global _CNT
    with _LK: _CNT += 1; p = _TMP / f"f_{os.getpid()}_{_CNT}.fits"
    for _ in range(2):
        try:
            Observations.download_file(uri, local_path=str(p))
            if p.exists() and p.stat().st_size > 1000: return p
        except Exception: pass
    p.unlink(missing_ok=True); return None

def index_mast_sectors():
    sm, obs_s, dvm = {}, {}, []
    for s in range(105, 113):
        try: ob = Observations.query_criteria(obs_collection="TESS", provenance_name="SPOC", dataproduct_type="timeseries", sequence_number=s)
        except Exception: ob = None
        if ob is None or len(ob) == 0:
            if s >= 106 and sm: break
            continue
        mp, rows = {}, []
        for r in ob:
            oid, tn = str(r["obs_id"]), re.sub(r"\D", "", str(r["target_name"]))
            if not tn or oid.endswith("_fast"): continue
            tid = int(tn)
            if re.search(r"-s\d{4}-s\d{4}-", oid):
                if s >= 106: dvm.append((tid, r))
            elif re.search(r"-s\d{4}-\d+-\d+-s$", oid):
                mp[tid] = f"mast:TESS/product/{oid}_lc.fits"; rows.append((tid, r))
        if mp: sm[s] = mp; obs_s[s] = rows
    alls = sorted(sm); ns = [s for s in alls if s >= 106] or alls[-2:]
    return (max(alls) if alls else 105), ns, {s: sm[s] for s in ns}, {s: obs_s[s] for s in ns}, dvm

def load_catalogs(eu_dir: Path = Path("_eumemic/results_public")):
    kt, kp, er = set(), {}, []
    def _a(tv, pv=None):
        s = re.sub(r"\D", "", str(tv))
        if s:
            kt.add(int(s))
            if pv is not None and pd.notna(pv) and float(pv) > 0: kp.setdefault(int(s), []).append(float(pv))
    for tb, ct, cp, wh in [("toi", "tid", "pl_orbper", None), ("pscomppars", "tic_id", "pl_orbper", "tic_id is not null")]:
        try:
            for r in NasaExoplanetArchive.query_criteria(table=tb, select=f"{ct},{cp}", **({"where": wh} if wh else {})):
                _a(r[ct], float(r[cp]) if r[cp] is not None and np.isfinite(r[cp]) else None)
        except Exception as e: print(f"[warn] NEA {tb}: {e}")
    try:
        vt = Vizier(columns=["TIC"], row_limit=10000).get_catalogs("J/ApJS/258/16/tess-ebs")
        if vt:
            for r in vt[0]: _a(r["TIC"])
    except Exception: pass
    for fn, tg in [("candidates_allsky.csv", "allsky"), ("candidates.csv", "earlier")]:
        fp = eu_dir / fn
        if fp.exists():
            de = pd.read_csv(fp); de["src"] = tg; er.append(de)
            for _, r in de.iterrows(): _a(r.get("TIC"), r.get("P_d"))
    return kt, kp, (pd.concat(er, ignore_index=True) if er else pd.DataFrame())

def find_dvt_seeds(ns, obs_s, dvm, kt, max_dvt: int = 60):
    du, seen, cobs = [], set(), []
    for tid, r in dvm + [x for s in reversed(ns) for x in obs_s.get(s, [])]:
        if tid not in kt and tid not in seen:
            seen.add(tid); cobs.append(r)
            if len(cobs) >= 500: break
    for i in range(0, len(cobs), 125):
        try:
            dv = Observations.filter_products(Observations.get_product_list(vstack(cobs[i : i + 125])), productSubGroupDescription=["DVT"])
            for r in dv:
                m = re.search(r"-0*(\d+)-\d+_dvt\.fits$", str(r["productFilename"]))
                if int(r["size"]) <= 22_000_000 and m and int(m.group(1)) not in kt:
                    du.append((int(m.group(1)), str(r["dataURI"])))
                    if len(du) >= max_dvt: break
        except Exception as e: print(f"[warn] DVT: {e}")
        if len(du) >= max_dvt: break
    out = {}
    def _rd(tu):
        tid, uri = tu; p = dl_uri(uri)
        if not p: return None
        try:
            with fits.open(p, memmap=False) as h:
                h0, h1 = h[0].header, h[1].header
                per, dep = float(h1.get("TPERIOD", 0.0) or 0.0), float(h1.get("TDEPTH", 0.0) or 0.0)
                teff, rs = float(h0.get("TEFF", 5000.0) or 5000.0), float(h0.get("RADIUS", 0.85) or 0.85)
                snr = float(h1.get("MODELSNR", h1.get("MAXMES", 9.0)) or 9.0)
                rp = math.sqrt(max(dep, 0.0) * 1e-6) * rs * RSUN_RE
                if not (0.16 <= per <= 45.0) or any(abs(per - a) / a < 0.012 for a in (0.5, 1.0, 6.85, 13.7)): return None
                return {"tic": tid, "p": per, "snr": snr} if (teff <= 6500 and rs <= 1.55 and 0.5 <= rp <= 6.0) else None
        except Exception: return None
        finally: p.unlink(missing_ok=True)
    with cf.ThreadPoolExecutor(max_workers=4) as ex:
        for x in ex.map(_rd, du[:max_dvt]):
            if x and (x["tic"] not in out or x["snr"] > out[x["tic"]]["snr"]): out[x["tic"]] = x
    return out

def dl_lc(uri: str, sh: int = 0):
    p = dl_uri(uri)
    if not p: return None
    try:
        with fits.open(p, memmap=False) as h:
            h0, h1, d = h[0].header, h[1].header, h[1].data; t = np.array(d["TIME"], float)
            p1 = np.where(np.isfinite(d["POS_CORR1"]), d["POS_CORR1"], 0.0) if "POS_CORR1" in d.names else np.zeros_like(t)
            p2 = np.where(np.isfinite(d["POS_CORR2"]), d["POS_CORR2"], 0.0) if "POS_CORR2" in d.names else np.zeros_like(t)
            return {"sec": int(h0.get("SECTOR", sh) or sh), "t": t, "f": np.array(d["PDCSAP_FLUX"], float), "fe": np.array(d["PDCSAP_FLUX_ERR"], float),
                    "sap": np.array(d["SAP_FLUX"], float), "q": np.array(d["QUALITY"], int), "mdx": np.array(d["MOM_CENTR1"], float) - p1,
                    "mdy": np.array(d["MOM_CENTR2"], float) - p2, "crowd": float(h1.get("CROWDSAP", 1.0) or 1.0),
                    "teff": float(h0.get("TEFF", np.nan) or np.nan), "rstar": float(h0.get("RADIUS", np.nan) or np.nan),
                    "tmag": float(h0.get("TESSMAG", np.nan) or np.nan), "ra": float(h0.get("RA_OBJ", np.nan) or np.nan), "dec": float(h0.get("DEC_OBJ", np.nan) or np.nan)}
    except Exception: return None
    finally: p.unlink(missing_ok=True)

def stack_tic(tid: int, win_d: float = 1.2):
    obs = Observations.query_criteria(obs_collection="TESS", dataproduct_type="timeseries", provenance_name="SPOC", target_name=str(tid))
    if obs is None or len(obs) == 0: obs = Observations.query_criteria(obs_collection="TESS", provenance_name="SPOC", obs_id=f"*{tid:016d}*")
    if obs is None or len(obs) == 0: return None, []
    pr = Observations.get_product_list(obs); seen, uris = set(), []
    for r in sorted([x for x in pr if str(x["productFilename"]).endswith(("_s_lc.fits", "-s_lc.fits"))], key=lambda x: str(x["productFilename"]), reverse=True):
        m = re.search(r"-s(\d{4})-", str(r["productFilename"])); s = int(m.group(1)) if m else 0
        if s in seen or s == 1751: continue
        seen.add(s); uris.append((s, str(r["dataURI"])))
        if len(uris) >= 18: break
    with cf.ThreadPoolExecutor(max_workers=4) as ex: lcs = [x for x in ex.map(lambda su: dl_lc(su[1], su[0]), uris) if x]
    ta, fa, fea, sa, sapa, mxa, mya, cr = [], [], [], [], [], [], [], []
    for l in sorted(lcs, key=lambda x: x["sec"]):
        tc, fc, fec = detrend_lc(l["t"], l["f"], l["fe"], l["q"], win_d)
        if len(tc) < 60: continue
        im = np.clip(np.searchsorted(l["t"], tc), 0, len(l["t"]) - 1)
        ta.append(tc); fa.append(fc); fea.append(fec); sa.append(np.full(len(tc), l["sec"], int))
        sapa.append(l["sap"][im]); mxa.append(l["mdx"][im]); mya.append(l["mdy"][im]); cr.append(l["crowd"])
    if len(ta) < 2: return None, pr
    tc = np.concatenate(ta); oi = np.argsort(tc)
    return {"t": tc[oi], "f": np.concatenate(fa)[oi], "fe": np.concatenate(fea)[oi], "sec": np.concatenate(sa)[oi],
            "sap": np.concatenate(sapa)[oi], "mdx": np.concatenate(mxa)[oi], "mdy": np.concatenate(mya)[oi], "crowd": float(np.median(cr)), "hdr": lcs[-1]}, pr

def gaia_nb(tid: int, ra: float, dec: float, hdr: dict):
    gm, nbs = {"teff": hdr["teff"], "rstar": hdr["rstar"], "tmag": hdr["tmag"]}, []
    try:
        tb = Catalogs.query_region(SkyCoord(ra=ra * u.deg, dec=dec * u.deg), radius=60.0 * u.arcsec, catalog="Tic")
        for r in (tb or []):
            rid, sep = int(r["ID"]), float(r["dstArcSec"]) if "dstArcSec" in tb.colnames else 0.0
            if rid == tid:
                for ko, ki in [("teff", "Teff"), ("rstar", "rad"), ("mstar", "mass"), ("tmag", "Tmag"), ("dist_pc", "d")]:
                    if ki in tb.colnames and r[ki] is not None and np.isfinite(r[ki]): gm[ko] = float(r[ki])
                gm["lumclass"] = str(r["lumclass"]) if "lumclass" in tb.colnames and r["lumclass"] is not None else "DWARF"
            elif sep <= 60.0 and r["Tmag"] is not None and np.isfinite(r["Tmag"]): nbs.append({"tic": rid, "sep_arcsec": sep, "tmag": float(r["Tmag"])})
    except Exception: pass
    try:
        dg = Gaia.launch_job(f"SELECT TOP 3 ruwe,radial_velocity_error FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.6f},{dec:.6f},0.0035)) ORDER BY phot_g_mean_mag ASC").get_results()
        if dg is not None and len(dg) > 0:
            if dg["ruwe"][0] is not None and np.isfinite(dg["ruwe"][0]): gm["ruwe"] = float(dg["ruwe"][0])
            if dg["radial_velocity_error"][0] is not None and np.isfinite(dg["radial_velocity_error"][0]): gm["rv_err_kms"] = float(dg["radial_velocity_error"][0])
    except Exception: pass
    return gm, nbs

def tpf_diff(pr, p_d: float, t0_btjd: float, dur_d: float):
    tpr = sorted([x for x in pr if str(x["productFilename"]).endswith(("_s_tp.fits", "-s_tp.fits"))], key=lambda x: str(x["productFilename"]), reverse=True)
    p = dl_uri(str(tpr[0]["dataURI"])) if tpr else None
    if not p: return None
    try:
        with fits.open(p, memmap=False) as h:
            sec = int(h[0].header.get("SECTOR", 0) or 0); t = np.array(h[1].data["TIME"], float); q = np.array(h[1].data["QUALITY"], int); fc = np.array(h[1].data["FLUX"], float)
        g = np.isfinite(t) & (q == 0); t, fc = t[g], fc[g]; ph = ((t - t0_btjd + 0.5 * p_d) % p_d) - 0.5 * p_d
        it, ol = np.abs(ph) <= 0.45 * dur_d, (np.abs(ph) > 1.25 * dur_d) & (np.abs(ph) <= min(3.5 * dur_d, 0.45 * p_d))
        if np.sum(it) < 3 or np.sum(ol) < 15: return None
        mo, mi = np.nanmedian(fc[ol], axis=0), np.nanmedian(fc[it], axis=0); df = mo - mi; yy, xx = np.indices(mo.shape)
        py, px = np.unravel_index(int(np.nanargmax(mo)), mo.shape)
        ws = np.clip(np.where((np.abs(yy - py) <= 2) & (np.abs(xx - px) <= 2) & np.isfinite(mo), mo - np.nanpercentile(mo, 20), 0.0), 0.0, None)
        nd = max(1.4826 * float(np.nanmedian(np.abs(df - np.nanmedian(df)))), 1e-4); wd = np.where(np.isfinite(df) & (df > 1.8 * nd), df, 0.0)
        if np.sum(ws) <= 0 or np.sum(wd) <= 0: return None
        sx, sy = float(np.sum(xx * ws) / np.sum(ws)), float(np.sum(yy * ws) / np.sum(ws)); dx, dy = float(np.sum(xx * wd) / np.sum(wd)), float(np.sum(yy * wd) / np.sum(wd))
        return {"sec": sec, "diff": df, "sxy": (sx, sy), "dxy": (dx, dy), "off": math.hypot(dx - sx, dy - sy) * PIX_AS}
    except Exception: return None
    finally: p.unlink(missing_ok=True)
