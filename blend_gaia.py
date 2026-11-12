"""Proper-motion-propagated Gaia DR3 WCS sub-pixel blend localizer and O-C/WCS figure plotter."""
from __future__ import annotations
import math
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from bls_search import BTJD0

def query_gaia_cone(ra: float, dec: float) -> list[dict]:
    from astroquery.gaia import Gaia
    q = (f"SELECT TOP 40 source_id,ra,dec,pmra,pmdec,phot_g_mean_mag,phot_rp_mean_mag,ruwe "
         f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.6f},{dec:.6f},0.0181)) ORDER BY phot_g_mean_mag ASC")
    out = []
    try:
        for r in (Gaia.launch_job(q).get_results() or []):
            g = float(r["phot_g_mean_mag"]) if r["phot_g_mean_mag"] is not None and np.isfinite(r["phot_g_mean_mag"]) else 17.0
            rp = float(r["phot_rp_mean_mag"]) if r["phot_rp_mean_mag"] is not None and np.isfinite(r["phot_rp_mean_mag"]) else (g - 0.5)
            out.append({"sid": str(r["source_id"]), "ra0": float(r["ra"]), "dec0": float(r["dec"]),
                        "pma": float(r["pmra"]) if r["pmra"] is not None and np.isfinite(r["pmra"]) else 0.0,
                        "pmd": float(r["pmdec"]) if r["pmdec"] is not None and np.isfinite(r["pmdec"]) else 0.0,
                        "tmag": round(g - 0.0052 - 0.58 * float(np.clip(g - rp, -0.5, 2.5)), 2),
                        "ruwe": round(float(r["ruwe"]) if r["ruwe"] is not None and np.isfinite(r["ruwe"]) else 1.0, 3)})
    except Exception: pass
    return out

def propagate_scene(raw_stars: list[dict], ra: float, dec: float, epoch_yr: float, dep_ppm: float):
    if not raw_stars: return {"sid": "none", "ra": ra, "dec": dec, "tmag": 12.0, "ruwe": 1.0, "pm_masyr": 0.0}, []
    dt = epoch_yr - 2016.0; stars = []
    for s in raw_stars:
        rae = s["ra0"] + (s["pma"] / (3.6e6 * max(math.cos(math.radians(s["dec0"])), 0.05))) * dt
        dce = s["dec0"] + (s["pmd"] / 3.6e6) * dt
        sep = math.hypot((rae - ra) * math.cos(math.radians(dec)), dce - dec) * 3600.0
        stars.append({**s, "ra": rae, "dec": dce, "pm_masyr": round(math.hypot(s["pma"], s["pmd"]), 1), "sep0_as": round(sep, 2)})
    c0 = [s for s in stars if s["sep0_as"] <= 5.5] or [min(stars, key=lambda s: s["sep0_as"])]
    targ = min(c0, key=lambda s: s["tmag"] + 0.35 * s["sep0_as"]); nbs, tap = [], 1.0
    for s in stars:
        if s["sid"] == targ["sid"]: continue
        sp = math.hypot((s["ra"] - targ["ra"]) * math.cos(math.radians(targ["dec"])), s["dec"] - targ["dec"]) * 3600.0
        fr = 10.0 ** (-0.4 * (s["tmag"] - targ["tmag"])); aw = fr * math.exp(-0.5 * (sp / 23.0) ** 2); tap += aw
        nbs.append({**s, "sep_as": round(sp, 2), "flux_ratio": fr, "ap_w": aw})
    for nb in nbs: nb["req_eb_depth"] = round((dep_ppm * 1e-6) * tap / max(nb["ap_w"], 1e-9), 4)
    return targ, sorted(nbs, key=lambda x: x["sep_as"])

def multisector_wcs_blend_audit(tic: int, ra: float, dec: float, p_d: float, t0_btjd: float, dur_d: float, dep_ppm: float, secs_hint=None):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astropy.wcs import WCS
    from astroquery.mast import Tesscut
    raw_g = query_gaia_cone(ra, dec); targ, nbs = propagate_scene(raw_g, ra, dec, 2026.7, dep_ppm)
    c0 = SkyCoord(ra=targ["ra"] * u.deg, dec=targ["dec"] * u.deg); secs = list(secs_hint or [])
    if not secs:
        try:
            st = Tesscut.get_sectors(coordinates=c0)
            col = next((c for c in st.colnames if "sec" in c.lower()), None) if st is not None and len(st) else None
            all_s = sorted({int(r[col]) for r in st}) if col else [106, 107]
            secs = sorted(set(all_s[:1] + all_s[len(all_s) // 2 : len(all_s) // 2 + 1] + all_s[-2:]))
        except Exception: secs = [106, 107]
    dra_l, ddec_l, w_l, snap = [], [], [], None
    for s in secs[:4]:
        try:
            hdus = Tesscut.get_cutouts(coordinates=c0, size=11, sector=int(s))
            if not hdus: continue
            h = hdus[0]; wcs = WCS(h[2].header).celestial; d = h[1].data
            t, q, fc = np.array(d["TIME"], float), np.array(d["QUALITY"], int), np.array(d["FLUX"], float)
            for hh in hdus: hh.close()
            g = np.isfinite(t) & (q == 0)
            if np.sum(g) < 60: continue
            t, fc = t[g], fc[g]; ep_yr = 2000.0 + (BTJD0 + float(np.nanmedian(t)) - 2451545.0) / 365.25
            t_ep, nbs_ep = propagate_scene(raw_g, ra, dec, ep_yr, dep_ppm)
            ph = ((t - t0_btjd + 0.5 * p_d) % p_d) - 0.5 * p_d
            it, ol = np.abs(ph) <= 0.48 * dur_d, (np.abs(ph) > 1.2 * dur_d) & (np.abs(ph) <= min(4.0 * dur_d, 0.45 * p_d))
            if np.sum(it) < 3 or np.sum(ol) < 15: continue
            mo, mi = np.nanmedian(fc[ol], axis=0), np.nanmedian(fc[it], axis=0); df = mo - mi; yy, xx = np.indices(mo.shape)
            tx0, ty0 = [float(v) for v in wcs.world_to_pixel_values(t_ep["ra"], t_ep["dec"])]
            wo = np.where((np.hypot(xx - tx0, yy - ty0) <= 2.3) & np.isfinite(mo), np.clip(mo - float(np.nanpercentile(mo[np.isfinite(mo)], 20)), 0.0, None), 0.0)
            if np.sum(wo) <= 0: continue
            ox, oy = float(np.sum(xx * wo) / np.sum(wo)), float(np.sum(yy * wo) / np.sum(wo))
            mx, my, mw = tx0, ty0, 1.0
            for nb in nbs_ep[:12]:
                if nb["sep_as"] <= 48.0:
                    nx, ny = [float(v) for v in wcs.world_to_pixel_values(nb["ra"], nb["dec"])]
                    if math.hypot(nx - tx0, ny - ty0) <= 2.3: mx += nx * nb["flux_ratio"]; my += ny * nb["flux_ratio"]; mw += nb["flux_ratio"]
            rx, ry = float(np.clip(ox - mx / mw, -1.3, 1.3)), float(np.clip(oy - my / mw, -1.3, 1.3)); tx, ty = tx0 + rx, ty0 + ry
            nd = max(1.4826 * float(np.nanmedian(np.abs(df[np.isfinite(df)] - np.nanmedian(df[np.isfinite(df)])))), 1e-4)
            wd = np.where((np.hypot(xx - tx, yy - ty) <= 3.3) & np.isfinite(df) & (df > 1.6 * nd), df, 0.0)
            snr_d = float(np.nanmax(np.where(np.hypot(xx - tx, yy - ty) <= 3.0, df, 0.0))) / nd
            if np.sum(wd) <= 0 or snr_d < 2.1: continue
            dx, dy = float(np.sum(xx * wd) / np.sum(wd)), float(np.sum(yy * wd) / np.sum(wd))
            rad_, dcd_ = [float(v) for v in wcs.pixel_to_world_values(dx - rx, dy - ry)]
            dra_l.append((rad_ - t_ep["ra"]) * math.cos(math.radians(t_ep["dec"])) * 3600.0)
            ddec_l.append((dcd_ - t_ep["dec"]) * 3600.0); w_l.append(1.0 / (max(2.2, 16.5 / max(snr_d, 1.5)) ** 2))
            if snap is None or snr_d > snap["snr"]:
                nbp = [(nx + rx, ny + ry, nb["req_eb_depth"] <= 1.0) for nb in nbs_ep[:10] for nx, ny in [[float(v) for v in wcs.world_to_pixel_values(nb["ra"], nb["dec"])]] if 0 <= nx + rx <= df.shape[1] - 1 and 0 <= ny + ry <= df.shape[0] - 1]
                snap = {"sec": int(s), "diff": np.nan_to_num(df, nan=0.0), "txy": (tx, ty), "dxy": (dx, dy), "nbs": nbp, "snr": snr_d}
        except Exception: continue
    plaus = [nb for nb in nbs if nb["req_eb_depth"] <= 1.0 and nb["sep_as"] <= 55.0]
    if not w_l:
        return {"tic": tic, "gaia_sid": targ["sid"], "ruwe": targ["ruwe"], "pm_masyr": targ["pm_masyr"], "n_sec": 0, "off_target_as": -1.0,
                "err_as": -1.0, "plausible_nbs": len(plaus), "closest_nb_sid": plaus[0]["sid"] if plaus else "none", "closest_nb_sep_as": plaus[0]["sep_as"] if plaus else -1.0, "verdict": "NO_DIFF_SNR"}, None
    wt = float(np.sum(w_l)); mdra, mddec = float(np.sum(np.array(dra_l) * np.array(w_l)) / wt), float(np.sum(np.array(ddec_l) * np.array(w_l)) / wt)
    off_t, err_t, bnb, boff = math.hypot(mdra, mddec), 1.0 / math.sqrt(wt), None, 999.0
    for nb in plaus:
        d_nb = math.hypot(mdra - (nb["ra"] - targ["ra"]) * math.cos(math.radians(targ["dec"])) * 3600.0, mddec - (nb["dec"] - targ["dec"]) * 3600.0)
        if d_nb < boff: bnb, boff = nb, d_nb
    vd = ("NEIGHBOR_BEB_LOCALIZED" if (bnb and off_t > max(9.5, 2.4 * err_t) and boff < min(0.65 * off_t, 12.0)) else
          "ON_TARGET_UNRESOLVED_NEIGHBOR" if (off_t <= max(8.5, 2.5 * err_t) and any(nb["sep_as"] < 7.5 for nb in plaus)) else
          "ON_TARGET_CLEAN" if off_t <= max(8.5, 2.5 * err_t) else "OFFSET_MARGINAL")
    return {"tic": tic, "gaia_sid": targ["sid"], "ruwe": targ["ruwe"], "pm_masyr": targ["pm_masyr"], "n_sec": len(w_l), "off_target_as": round(off_t, 2),
            "err_as": round(err_t, 2), "plausible_nbs": len(plaus), "closest_nb_sid": bnb["sid"] if bnb else "none", "closest_nb_sep_as": bnb["sep_as"] if bnb else -1.0, "verdict": vd}, snap

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
    ax[0].set_title(f"S1-S107 O-C ({fit['class']})", color="#f8fafc", fontsize=8.5)
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
