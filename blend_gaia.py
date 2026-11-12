"""Cloud-proof VizieR I/355/gaiadr3 + MAST TIC neighbor localizer and WCS difference-image audit."""
from __future__ import annotations
import math
import numpy as np
from bls_search import BTJD0

def query_gaia_cone(ra: float, dec: float) -> list[dict]:
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astroquery.mast import Catalogs
    from astroquery.vizier import Vizier
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg); out = []
    try:
        vz = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "pmRA", "pmDE", "Gmag", "RPmag", "RUWE", "sepsi", "NSS", "e_RV"], row_limit=45)
        res = vz.query_region(c, radius=65.0 * u.arcsec, catalog="I/355/gaiadr3")
        if res and len(res) > 0:
            for r in res[0]:
                g = float(r["Gmag"]) if "Gmag" in r.colnames and r["Gmag"] is not None and np.isfinite(r["Gmag"]) else 17.0
                rp = float(r["RPmag"]) if "RPmag" in r.colnames and r["RPmag"] is not None and np.isfinite(r["RPmag"]) else (g - 0.5)
                out.append({"sid": str(r["Source"]), "ra0": float(r["RA_ICRS"]), "dec0": float(r["DE_ICRS"]),
                            "pma": float(r["pmRA"]) if "pmRA" in r.colnames and r["pmRA"] is not None and np.isfinite(r["pmRA"]) else 0.0,
                            "pmd": float(r["pmDE"]) if "pmDE" in r.colnames and r["pmDE"] is not None and np.isfinite(r["pmDE"]) else 0.0,
                            "tmag": round(g - 0.0052 - 0.58 * float(np.clip(g - rp, -0.5, 2.5)), 2),
                            "ruwe": round(float(r["RUWE"]) if "RUWE" in r.colnames and r["RUWE"] is not None and np.isfinite(r["RUWE"]) else 1.0, 3),
                            "sepsi": round(float(r["sepsi"]) if "sepsi" in r.colnames and r["sepsi"] is not None and np.isfinite(r["sepsi"]) else 0.0, 2),
                            "nss": int(r["NSS"]) if "NSS" in r.colnames and r["NSS"] is not None and np.isfinite(r["NSS"]) else 0,
                            "rv_err": round(float(r["e_RV"]) if "e_RV" in r.colnames and r["e_RV"] is not None and np.isfinite(r["e_RV"]) else -1.0, 2), "ep0": 2016.0})
    except Exception: pass
    if out: return out
    try:
        for r in (Catalogs.query_region(c, radius=65.0 * u.arcsec, catalog="Tic") or []):
            tm = float(r["Tmag"]) if r["Tmag"] is not None and np.isfinite(r["Tmag"]) else 17.0
            gid = str(r["GAIA"]) if "GAIA" in r.colnames and r["GAIA"] is not None and str(r["GAIA"]).strip() not in ("", "--", "None") else f"TIC:{r['ID']}"
            out.append({"sid": gid, "ra0": float(r["ra"]), "dec0": float(r["dec"]),
                        "pma": float(r["pmRA"]) if "pmRA" in r.colnames and r["pmRA"] is not None and np.isfinite(r["pmRA"]) else 0.0,
                        "pmd": float(r["pmDEC"]) if "pmDEC" in r.colnames and r["pmDEC"] is not None and np.isfinite(r["pmDEC"]) else 0.0,
                        "tmag": round(tm, 2), "ruwe": 1.0, "sepsi": 0.0, "nss": 0, "rv_err": -1.0, "ep0": 2015.5})
    except Exception: pass
    return out

def propagate_scene(raw_stars: list[dict], ra: float, dec: float, epoch_yr: float, dep_ppm: float):
    if not raw_stars: return {"sid": "none", "ra": ra, "dec": dec, "tmag": 12.0, "ruwe": 1.0, "pm_masyr": 0.0, "sepsi": 0.0, "nss": 0, "rv_err": -1.0}, []
    stars = []
    for s in raw_stars:
        dt = epoch_yr - float(s.get("ep0", 2016.0))
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
    raw_g = query_gaia_cone(ra, dec); targ, nbs = propagate_scene(raw_g, ra, dec, 2026.8, dep_ppm)
    c0 = SkyCoord(ra=targ["ra"] * u.deg, dec=targ["dec"] * u.deg); secs = list(secs_hint or []) or [106, 107]
    dra_l, ddec_l, w_l, snap = [], [], [], None
    for s in secs[-4:]:
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
