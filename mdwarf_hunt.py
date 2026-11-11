#!/usr/bin/env python3
"""Deep multi-sector search of uncatalogued S106-S107 M-dwarfs (including T > 13.5)."""
from __future__ import annotations
import concurrent.futures as cf, json, math, time
from pathlib import Path
import numpy as np, pandas as pd
from astroquery.mast import Catalogs
from bls_search import BTJD0, RSUN_RE, search_bls
from deep_search import load_all_eumemic
from ledger import write_ledger
from mast_data import gaia_nb, index_mast_sectors, load_catalogs, stack_tic, tpf_diff
from vetter import plot_card, vet_signal

def select_uncatalogued_mdwarfs(smaps: dict, new_s: list[int], kt: set, max_stars: int = 60) -> list[dict]:
    ov = [t for t in sorted(set(smaps[new_s[-2]]) & set(smaps[new_s[-1]])) if t not in kt] if len(new_s) >= 2 else []
    if len(ov) < 400:
        extra = [t for s in reversed(new_s) for t in smaps[s] if t not in kt and t not in set(ov)]
        ov += extra[: max(0, 1500 - len(ov))]
    sample = ov[::max(1, len(ov) // 1500)][:1500]
    mdw = []
    for i in range(0, len(sample), 150):
        try:
            tb = Catalogs.query_criteria(catalog="Tic", ID=sample[i : i + 150])
            for r in (tb or []):
                tid = int(r["ID"])
                lc = str(r["lumclass"]) if "lumclass" in tb.colnames and r["lumclass"] is not None else "DWARF"
                teff = float(r["Teff"]) if "Teff" in tb.colnames and r["Teff"] is not None and np.isfinite(r["Teff"]) else np.nan
                rad = float(r["rad"]) if "rad" in tb.colnames and r["rad"] is not None and np.isfinite(r["rad"]) else np.nan
                tmag = float(r["Tmag"]) if "Tmag" in tb.colnames and r["Tmag"] is not None and np.isfinite(r["Tmag"]) else np.nan
                if lc.upper() not in ("DWARF", "", "NONE") or not (np.isfinite(teff) and np.isfinite(rad) and np.isfinite(tmag)):
                    continue
                if 2600 <= teff <= 4150 and 0.12 <= rad <= 0.54 and 9.0 <= tmag <= 15.6:
                    score = (0 if tmag > 13.5 else 1, rad + 0.04 * max(0.0, tmag - 13.5))
                    mdw.append({"tic": tid, "teff": teff, "rstar": rad, "tmag": tmag, "score": score})
        except Exception as e:
            print(f"[warn] TIC chunk: {e}")
    mdw.sort(key=lambda d: d["score"])
    print(f"[info] selected {min(len(mdw), max_stars)} / {len(mdw)} uncatalogued multi-sector M-dwarfs from {len(sample)} S106/S107 stars")
    return mdw[:max_stars]

def search_one_mdwarf(item: dict):
    tid, rs0 = item["tic"], item["rstar"]
    try:
        st, pr = stack_tic(tid, 1.05)
        if not st or len(np.unique(st["sec"])) < 2:
            return None
        t, f, fe, sc, sap, mx, my = st["t"], st["f"], st["fe"], st["sec"], st["sap"], st["mdx"], st["mdy"]
        rs = float(st["hdr"]["rstar"]) if np.isfinite(st["hdr"]["rstar"]) and st["hdr"]["rstar"] > 0 else rs0
        pk = search_bls(t, f, fe, 0.165, 38.0, rs, fast_inv=False)
        if not pk or pk["mes_c"] < 8.8 or pk["sde"] < 7.2:
            return None
        rp = math.sqrt(pk["dep"]) * rs * RSUN_RE
        if not (0.45 <= rp <= 6.5):
            return None
        gm, nbs = gaia_nb(tid, st["hdr"]["ra"], st["hdr"]["dec"], st["hdr"])
        lane = "MDWARF_T>13.5" if item["tmag"] > 13.5 else "MDWARF_STACK"
        rep = vet_signal(tid, pk, t, f, fe, sc, lane, sap, mx, my, st["crowd"], gm, nbs)
        tpfi = tpf_diff(pr, rep["P_d"], rep["T0_BJD"] - BTJD0, rep["dur_h"] / 24.0) if rep["passed"] else None
        if tpfi:
            rep = vet_signal(tid, pk, t, f, fe, sc, lane, sap, mx, my, st["crowd"], gm, nbs, tpfi["off"])
        return rep, ({"t": t, "f": f, "sec": sc} if rep["passed"] else None), tpfi
    except Exception:
        return None

def main(out_dir: Path = Path("results"), max_mdw: int = 60):
    t0w = time.time(); out_dir.mkdir(parents=True, exist_ok=True); fig_d = out_dir / "figures"; fig_d.mkdir(exist_ok=True)
    s_max, new_s, smaps, _, _ = index_mast_sectors()
    kt, _, _ = load_catalogs()
    dfe = load_all_eumemic()
    for _, r in (dfe.iterrows() if not dfe.empty else []):
        kt.add(int(r["_tid"]))

    bres = pd.read_csv(out_dir / "blind_verification_s106_s107.csv").to_dict("records") if (out_dir / "blind_verification_s106_s107.csv").exists() else []
    prev_v = pd.read_csv(out_dir / "discovery_triggers_vetted.csv").to_dict("records") if (out_dir / "discovery_triggers_vetted.csv").exists() else []
    prev_s = json.loads((out_dir / "campaign_summary.json").read_text()) if (out_dir / "campaign_summary.json").exists() else {}

    mdw_pool = select_uncatalogued_mdwarfs(smaps, new_s, kt, max_stars=max_mdw)
    vall_map = {int(r["tic"]): r for r in prev_v}
    vpass = [r for r in prev_v if str(r.get("passed")).lower() == "true"]

    with cf.ThreadPoolExecutor(max_workers=5) as ex:
        for res in ex.map(search_one_mdwarf, mdw_pool):
            if not res:
                continue
            rep, ms, tpfi = res
            vall_map[int(rep["tic"])] = rep
            if rep["passed"] and ms:
                vpass.append(rep)
                plot_card(fig_d / f"discovery_TIC{rep['tic']}.png", rep["tic"], f"{rep['tier']} (MES={rep['MES']:.1f}, Rp={rep['Rp_Re']:.2f} Re, FPP={rep['FPP']:.4f})", rep["P_d"], rep["T0_BJD"] - BTJD0, rep["dur_h"] / 24.0, rep["depth_ppm"], ms["t"], ms["f"], ms["sec"], tpfi)

    vall = sorted(vall_map.values(), key=lambda x: float(x["MES"]), reverse=True)
    vpass = sorted({int(x["tic"]): x for x in vpass}.values(), key=lambda x: float(x["MES"]), reverse=True)
    pd.DataFrame(vall).to_csv(out_dir / "discovery_triggers_vetted.csv", index=False)
    pd.DataFrame(vpass).to_csv(out_dir / "vetted_new_candidates.csv", index=False)

    vc = pd.DataFrame(bres)["verdict"].value_counts().to_dict() if bres else {}
    summ = {
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "sec_s": round(float(prev_s.get("sec_s", 0.0)) + (time.time() - t0w), 1),
        "s_max": s_max,
        "new_sectors": new_s,
        "blind_pairs": len(bres),
        "blind_counts": vc,
        "unpromoted_tce_pool_total": int(prev_s.get("unpromoted_tce_pool_total", 0)) + len(mdw_pool),
        "tce_stacked": int(prev_s.get("tce_stacked", 0)) + len(mdw_pool),
        "mdwarfs_searched": len(mdw_pool),
        "usp_hosts": int(prev_s.get("usp_hosts", 28)),
        "triggers_vetted": len(vall),
        "new_passed": len(vpass),
    }
    (out_dir / "campaign_summary.json").write_text(json.dumps(summ, indent=2))
    write_ledger(out_dir, summ, bres, vpass, vall)
    print(json.dumps(summ, indent=2))

if __name__ == "__main__":
    main()
    main()
