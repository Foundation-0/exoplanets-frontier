"""FFI Tesscut Sector 106-107 blind verifier and planet-depth seed discovery."""
from __future__ import annotations
import concurrent.futures as cf, math, re
from pathlib import Path
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import vstack
import astropy.units as u
from astroquery.mast import Catalogs, Observations, Tesscut
from scipy import signal
from bls_search import BTJD0, RSUN_RE, detrend_lc, search_bls
from campaigns import gv
from mast_data import dl_lc, dl_uri
from vetter import blind_check_sector, plot_card

def hp_cent(c):
    c = np.asarray(c, float); v = np.isfinite(c)
    if np.sum(v) < 35: return c
    tr = signal.medfilt(np.interp(np.arange(len(c)), np.where(v)[0], c[v]), kernel_size=min(len(c) | 1, 301))
    return c - tr

def load_all_eumemic(root: Path = Path("_eumemic")) -> pd.DataFrame:
    frames = []
    for fp in sorted(root.rglob("*.csv")):
        if "funnel" in fp.name.lower() or "sector" in fp.name.lower(): continue
        try:
            df = pd.read_csv(fp)
            cl = {c.strip().lower() for c in df.columns}
            if ("tic" in cl or "tic_id" in cl) and ("p_d" in cl or "period_d" in cl or "period" in cl):
                df["src_file"] = fp.name; frames.append(df)
                print(f"[info] loaded {fp}: {len(df)} rows")
        except Exception: pass
    if not frames: return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    out["_tid"] = [int(gv(r, "TIC", "tic", "tic_id")) if np.isfinite(gv(r, "TIC", "tic", "tic_id")) else -1 for _, r in out.iterrows()]
    out["_p"] = [round(gv(r, "P_d", "period_d", "period"), 3) for _, r in out.iterrows()]
    return out[out["_tid"] > 0].drop_duplicates(subset=["_tid", "_p"]).reset_index(drop=True)

def ffi_cutout_lc(tid: int, ra: float, dec: float, sec: int):
    try:
        hdus = Tesscut.get_cutouts(coordinates=SkyCoord(ra=ra * u.deg, dec=dec * u.deg), size=9, sector=sec)
        if not hdus: return None
        d = hdus[0][1].data; t = np.array(d["TIME"], float); q = np.array(d["QUALITY"], int); fc = np.array(d["FLUX"], float)
        g = np.isfinite(t) & (q == 0)
        if np.sum(g) < 80: return None
        t, fc = t[g], fc[g]; med = np.nanmedian(fc, axis=0); yy, xx = np.indices(med.shape)
        cy, cx = med.shape[0] // 2, med.shape[1] // 2
        sub = np.where((np.abs(yy - cy) <= 1) & (np.abs(xx - cx) <= 1) & np.isfinite(med), med, -np.inf)
        py, px = np.unravel_index(int(np.argmax(sub)), med.shape)
        ap = (np.abs(yy - py) <= 1) & (np.abs(xx - px) <= 1); bg = (np.abs(yy - py) >= 3) | (np.abs(xx - px) >= 3)
        bkg = np.nanmedian(fc[:, bg], axis=1)
        fl = np.nansum(fc[:, ap] - bkg[:, None], axis=1)
        w = np.clip(fc[:, ap] - bkg[:, None], 0.0, None); ws = np.maximum(np.sum(w, axis=1), 1e-6)
        mx = np.sum(w * xx[ap][None, :], axis=1) / ws; my = np.sum(w * yy[ap][None, :], axis=1) / ws
        fe = np.full_like(fl, max(1.4826 * float(np.nanmedian(np.abs(np.diff(fl)))), 1e-3))
        for h in hdus: h.close()
        return {"sec": sec, "t": t, "f": fl, "fe": fe, "q": np.zeros(len(t), int), "mdx": mx, "mdy": my}
    except Exception: return None

def run_full_blind(dfe: pd.DataFrame, new_s: list[int], smaps: dict, obs_by_s: dict, kt: set, kp: dict, fig_d: Path):
    e_by_tic, cand_rows = {}, []
    for _, row in dfe.iterrows():
        tid, pv = int(row["_tid"]), gv(row, "P_d", "period_d", "period")
        kt.add(tid)
        if np.isfinite(pv) and pv > 0: kp.setdefault(tid, []).append(pv)
        tr = str(row.get("tier", row.get("Tier", row.get("src_file", "v3"))) or "v3").strip()
        if any(x in tr.lower() for x in ("rej", "false", "eb", "l2")): continue
        e_by_tic.setdefault(tid, []).append(row); cand_rows.append((tid, tr, row))

    coords = {}
    need_c = [t for t, _, r in cand_rows if not (np.isfinite(gv(r, "RA", "ra", "ra_deg")) and np.isfinite(gv(r, "Dec", "dec", "dec_deg")))]
    for i in range(0, len(need_c), 100):
        try:
            tb = Catalogs.query_criteria(catalog="Tic", ID=need_c[i : i + 100])
            for r in (tb or []): coords[int(r["ID"])] = (float(r["ra"]), float(r["dec"]))
        except Exception: pass
    for t, _, r in cand_rows:
        ra, dc = gv(r, "RA", "ra", "ra_deg"), gv(r, "Dec", "dec", "dec_deg")
        if np.isfinite(ra) and np.isfinite(dc): coords[t] = (ra, dc)

    sec_sky = {}
    for s in new_s:
        ras = [float(r["s_ra"]) for _, r in obs_by_s.get(s, []) if r["s_ra"] is not None and np.isfinite(r["s_ra"])]
        dcs = [float(r["s_dec"]) for _, r in obs_by_s.get(s, []) if r["s_dec"] is not None and np.isfinite(r["s_dec"])]
        if ras: sec_sky[s] = SkyCoord(ra=np.array(ras) * u.deg, dec=np.array(dcs) * u.deg)

    tasks = []
    for tid, tr, row in cand_rows:
        for s in new_s:
            if tid in smaps[s]: tasks.append(("SPOC2m", tid, tr, row, s, smaps[s][tid]))
            elif tid in coords and s in sec_sky:
                c0 = SkyCoord(ra=coords[tid][0] * u.deg, dec=coords[tid][1] * u.deg)
                if float(np.min(c0.separation(sec_sky[s]).deg)) <= 2.4:
                    tasks.append(("FFI", tid, tr, row, s, coords[tid]))
    print(f"[info] blind candidate-sector tasks queued: {len(tasks)} (from {len(cand_rows)} v3 candidates)")

    bres, bhosts, bplots = [], {}, []
    def _eval(tk):
        mode, tid, tr, row, s, arg = tk
        lc = dl_lc(arg, s) if mode == "SPOC2m" else ffi_cutout_lc(tid, arg[0], arg[1], s)
        if not lc: return None
        pd_, t0b = gv(row, "P_d", "period_d", "period"), gv(row, "T0_BJD", "t0_bjd", "T0_BTJD", "t0_btjd", "t0")
        dppm, durh = gv(row, "depth_ppm", "depth", d=500.0), gv(row, "T14_h", "dur_h", "duration_h", d=2.2)
        if not (np.isfinite(pd_) and np.isfinite(t0b)): return None
        if dppm < 1.0: dppm *= 1e6
        tc, fc, _ = detrend_lc(lc["t"], lc["f"], lc["fe"], lc["q"], win_d=max(1.2, min(2.4, (durh / 24.0) * 8.0)))
        im = np.clip(np.searchsorted(lc["t"], tc), 0, len(lc["t"]) - 1) if len(tc) else np.array([], int)
        mx, my = (hp_cent(lc["mdx"])[im], hp_cent(lc["mdy"])[im]) if len(tc) else (None, None)
        rb = blind_check_sector(tid, f"{tr}({mode})", pd_, t0b, dppm, durh, s, tc, fc, mx, my)
        return tid, row, pd_, t0b, durh, tc, fc, s, rb

    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        for out in ex.map(_eval, tasks[:75]):
            if not out: continue
            tid, row, pd_, t0b, durh, tc, fc, s, rb = out
            if rb["verdict"] == "NO_DATA": continue
            bhosts[tid] = row; bres.append(rb)
            if rb["verdict"] in ("CONFIRMED_ON_SCHEDULE", "TIMING_OFFSET_TTV") and len(bplots) < 14:
                t0_bt = t0b - BTJD0 if t0b > 2400000 else t0b
                bplots.append((fig_d / f"blind_S{s}_TIC{tid}.png", tid, f"S{s} Blind {rb['verdict']} (SNR={rb['snr']:.1f}, Obs={rb['obs_ppm']:.0f} ppm)", pd_, t0_bt, durh / 24.0, rb["obs_ppm"], tc, fc, np.full(len(tc), s)))
    for ap in bplots: plot_card(*ap)
    bres.sort(key=lambda r: r["snr"], reverse=True)
    return bres, bhosts, e_by_tic
