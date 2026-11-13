"""Eight-gate S1-S108 transit-shape, ephemeris & WCS centroid evaluator for cool Platinum/Gold worlds."""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
from blend_gaia import multisector_wcs_blend_audit
from bls_search import BTJD0, RSUN_RE, exp_dur
from converge_io import plot_ttv_and_blend
from cool_io import coherent_planet_bls, load_cycles_1_to_8, stack_clean
from mdwarf_hunt import spot_amp
from ttv_track import extract_transit_times, fit_ephemeris_models


def vet_cool_candidate(rp_row: dict, fig_d: Path) -> dict:
    tid, pc = int(rp_row["tic"]), float(rp_row["P_d"])
    ra, dec = float(rp_row["ra_2026"]), float(rp_row["dec_2026"])
    rs, teff, tmag = max(float(rp_row["Rstar"]), 0.10), max(float(rp_row["Teff"]), 2600.0), float(rp_row["Tmag"])
    crat, ut, dcls = max(float(rp_row.get("contratio", 0.02)), 0.0), str(rp_row["upstream_tier"]), str(rp_row["dr4_class"])
    raw, alls = load_cycles_1_to_8(tid, ra, dec)
    miss = lambda ec, vt: {"tic": tid, "upstream_tier": ut, "dr4_class": dcls, "Tmag": round(tmag, 2), "Teff": round(teff, 0),
                           "Rstar": round(rs, 3), "P_cat_d": round(pc, 6), "P_ref_d": -1.0, "T0_ref_BJD": -1.0, "dP_sec": 0.0,
                           "n_sec_avail": len(alls), "n_sec_used": len(raw), "n_sec_hit": 0, "n_transits": 0, "std_MES": 0.0,
                           "despot_MES": 0.0, "SDE": 0.0, "depth_ppm": 0.0, "Rp_fit_Re": -1.0, "Teq_fit_K": -1.0, "odd_even_s": -1.0,
                           "sec_ecl_snr": 0.0, "wcs_off_as": -1.0, "wcs_verdict": "NONE", "ephem_class": ec, "survived": False, "vetoes": vt}
    if len(raw) < 2:
        return miss("LOW_DATA", "insufficient_sectors")
    ed = exp_dur(pc, rs)
    lcs_s, std, dsp = stack_clean(raw, max(1.6, min(3.4, 9.0 * ed)), max(1.5, 7.5 * ed))
    pks = coherent_planet_bls(std["t"], std["f"], std["fe"], pc, rs) if std else None
    pkd = coherent_planet_bls(dsp["t"], dsp["f"], dsp["fe"], pc, rs) if dsp else None
    prot = float(np.median([x["p_rot_d"] for x in (dsp["infs"] if dsp else []) if x["p_rot_d"] > 0] or [0.0]))
    pk, st = (pkd, dsp) if (pkd and (not pks or pkd["mes_c"] >= 0.90 * pks["mes_c"])) else (pks, std)
    if not pk or not st:
        return miss("NO_DIP", "no_coherent_dip")
    sr, ir = extract_transit_times(lcs_s, pk["p"], pk["t0"], pk["dur"], 210.0)
    fit = fit_ephemeris_models(sr, ir, pk["p"], pk["t0"])
    p_u, t0_u = (float(fit["P_ref"]), float(fit["T0_ref_BJD"]) - BTJD0) if fit["n_ep"] >= 3 else (pk["p"], pk["t0"])
    wr, bsnap = multisector_wcs_blend_audit(tid, ra, dec, p_u, t0_u, pk["dur"], pk["dep"] * 1e6, alls[-4:])
    t, f, sec = st["t"], st["f"], st["sec"]; dur, dep, sc = pk["dur"], pk["dep"], pk["sc"]
    ph = ((t - t0_u + 0.5 * p_u) % p_u) - 0.5 * p_u; orb = np.round((t - t0_u) / p_u).astype(int)
    it, ot = np.abs(ph) <= 0.5 * dur, np.abs(ph) > 1.25 * dur; med_o = float(np.median(f[ot]))
    rp = math.sqrt(max(dep, 1e-6) * (1.0 + crat)) * rs * RSUN_RE
    a_au = ((max(rs ** 0.9, 0.08) * (p_u / 365.25) ** 2) ** (1.0 / 3.0))
    teq = 278.5 * ((((rs ** 2) * ((teff / 5772.0) ** 4)) / max(a_au ** 2, 1e-6)) ** 0.25) * (0.7 ** 0.25)
    vet = []
    if pk["mes_c"] < 8.5: vet.append(f"low_mes({pk['mes_c']:.1f})")
    if pk["sde"] < 6.8: vet.append(f"low_sde({pk['sde']:.1f})")
    if not (0.50 <= rp <= 13.2): vet.append(f"rp({rp:.2f}Re)")
    if not (0.22 <= dur / max(exp_dur(p_u, rs), 1e-3) <= 2.5): vet.append("dur_ratio")
    core, wing = np.abs(ph) <= 0.22 * dur, (np.abs(ph) > 0.28 * dur) & (np.abs(ph) <= 0.5 * dur)
    if np.sum(core) >= 4 and np.sum(wing) >= 4 and dep > 0.0025:
        dc, dw = med_o - float(np.mean(f[core])), med_o - float(np.mean(f[wing]))
        if dc > 0 and (dw / dc) < 0.26: vet.append("grazing_v")
    btm = np.abs(ph) <= 0.32 * dur; io_, ie_ = btm & ((orb % 2) != 0), btm & ((orb % 2) == 0); oes = 0.0
    if np.sum(io_) >= 3 and np.sum(ie_) >= 3:
        do, de = med_o - float(np.median(f[io_])), med_o - float(np.median(f[ie_]))
        so = max(1.4826 * float(np.median(np.abs(f[io_] - np.median(f[io_])))), sc) / math.sqrt(int(np.sum(io_)))
        se = max(1.4826 * float(np.median(np.abs(f[ie_] - np.median(f[ie_])))), sc) / math.sqrt(int(np.sum(ie_)))
        oes = abs(do - de) / math.hypot(so, se)
        if do <= 0 or de <= 0 or oes > 3.2: vet.append(f"odd_even({oes:.2f}s)")
    ssnr = up_snr = 0.0
    for c in np.linspace(0.16, 0.84, 33):
        sm = np.abs(((t - t0_u) / p_u) % 1.0 - c) <= (0.5 * dur / p_u)
        if np.sum(sm) >= 4:
            ssnr = max(ssnr, (med_o - float(np.mean(f[sm]))) / (sc / math.sqrt(int(np.sum(sm)))))
            if len(np.unique(orb[sm])) >= 2: up_snr = max(up_snr, (float(np.mean(f[sm])) - med_o) / (sc / math.sqrt(int(np.sum(sm)))))
    if ssnr > 4.2 and ssnr > 0.25 * pk["mes"]: vet.append(f"sec_ecl({ssnr:.1f})")
    if up_snr > 0.75 * pk["mes"] and up_snr > 7.5: vet.append(f"inv_bump({up_snr:.1f})")
    evq = [(max(0.0, med_o - float(np.mean(f[it & (orb == o)]))) / (sc / math.sqrt(int(np.sum(it & (orb == o)))))) ** 2 for o in np.unique(orb[it]) if np.sum(it & (orb == o)) >= 1]
    if (max(evq) / max(sum(evq), 1e-9) if evq else 1.0) > 0.58: vet.append("single_ev")
    n_sec_hit = int(len(np.unique(sec[it])))
    if n_sec_hit < 2: vet.append("1_sector")
    amp = spot_amp(t, f, p_u, t0_u, dur)
    if amp > 0.48 * dep: vet.append(f"spot_wave({amp/max(dep,1e-6):.2f})")
    if prot > 0 and any(abs(p_u - k * prot) / (k * prot) < 0.02 for k in (0.5, 1.0, 2.0)): vet.append("rot_alias")
    wvd, woff = str(wr.get("verdict", "NONE")), float(wr.get("off_target_as", -1.0))
    if wvd == "NEIGHBOR_BEB_LOCALIZED" or (wvd == "OFFSET_MARGINAL" and woff > 11.5): vet.append(f"wcs({wvd}:{woff:.1f}as)")
    if float(wr.get("ruwe", 1.0)) > 1.4: vet.append(f"ruwe({wr['ruwe']:.2f})")
    plot_ttv_and_blend(fig_d / f"cool_TIC{tid}.png", tid, sr, fit, t, f, p_u, t0_u, dur, bsnap)
    ok_s = len(vet) == 0
    return {**miss(fit["class"], "CLEAN" if ok_s else ";".join(vet)), "P_ref_d": round(p_u, 7), "T0_ref_BJD": round(t0_u + BTJD0, 5),
            "dP_sec": fit["dP_sec"], "n_sec_hit": n_sec_hit, "n_transits": max(pk["ntr"], len(evq)),
            "std_MES": round(pks["mes_c"], 2) if pks else 0.0, "despot_MES": round(pkd["mes_c"], 2) if pkd else 0.0, "SDE": round(pk["sde"], 2),
            "depth_ppm": round(dep * 1e6, 1), "Rp_fit_Re": round(rp, 2), "Teq_fit_K": round(teq, 0), "odd_even_s": round(oes, 2),
            "sec_ecl_snr": round(ssnr, 2), "wcs_off_as": woff, "wcs_verdict": wvd, "survived": ok_s}
