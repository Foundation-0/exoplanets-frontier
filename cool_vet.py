#!/usr/bin/env python3
# Vet the coolest purified Platinum and Gold candidates across TESS Sectors 1-108
from __future__ import annotations
import hashlib, json, time
from pathlib import Path
import pandas as pd
from cool_core import vet_cool_candidate

MANDATORY_TICS = [399695544, 409696823, 224279805]


def select_coolest_targets(df_plat: pd.DataFrame, max_n: int = 10) -> list[dict]:
    df_s = df_plat.sort_values(["Teq_K", "priority_score"], ascending=[True, False]).reset_index(drop=True)
    picked, seen = [], set()
    for tid in MANDATORY_TICS:
        sub = df_s[df_s["tic"] == tid]
        if not sub.empty:
            r = sub.iloc[0].to_dict()
            key = (int(r["tic"]), round(float(r["P_d"]), 2))
            if key not in seen:
                seen.add(key); picked.append(r)
    for _, r0 in df_s.iterrows():
        if len(picked) >= max_n:
            break
        r = r0.to_dict()
        key = (int(r["tic"]), round(float(r["P_d"]), 2))
        if key not in seen:
            seen.add(key); picked.append(r)
    picked.sort(key=lambda x: float(x["Teq_K"]))
    return picked


def main(out_dir: Path = Path("results")):
    t0 = time.time()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_d = out_dir / "figures"
    fig_d.mkdir(exist_ok=True)
    plat_p = out_dir / "platinum_followup_targets.csv"
    assert plat_p.exists(), "results/platinum_followup_targets.csv missing"
    df_plat = pd.read_csv(plat_p)
    targets = select_coolest_targets(df_plat, max_n=10)
    rows = [vet_cool_candidate(r, fig_d) for r in targets]
    pd.DataFrame(rows).to_csv(out_dir / "cool_platinum_vetted.csv", index=False)
    surv = [r for r in rows if r["survived"]]
    pd.DataFrame(surv).to_csv(out_dir / "verified_cool_survivors.csv", index=False)

    note_map, bad_tics = {}, set()
    for r in rows:
        tid = int(r["tic"])
        if r["survived"]:
            note_map[tid] = f"verified_S1_S108_survivor(P={r['P_ref_d']}d,off={r['wcs_off_as']}as)"
        elif any(k in r["vetoes"] for k in ("odd_even", "wcs(", "sec_ecl", "grazing_v", "spot_wave", "ruwe", "inv_bump")):
            note_map[tid] = f"multi_sec_flag:{r['vetoes'][:48]}"
            bad_tics.add(tid)

    for csv_nm, keep_only_clean in [("eumemic_gaia_dr4_readiness.csv", False), ("platinum_followup_targets.csv", True)]:
        fp = out_dir / csv_nm
        if not fp.exists():
            continue
        dfc = pd.read_csv(fp)
        for tid, nt in note_map.items():
            m = dfc["tic"] == tid
            dfc.loc[m, "frontier_note"] = nt
            if tid in bad_tics:
                dfc.loc[m, "dr4_class"] = "MULTI_SECTOR_VETTED_WARNING"
        if keep_only_clean:
            dfc = dfc[dfc["dr4_class"].isin(["PLATINUM_SOLO_VERIFIED", "GOLD_LOW_BLEND"])].reset_index(drop=True)
        dfc.to_csv(fp, index=False)

    el = round(time.time() - t0, 1)
    vt_k = [float(r["Teq_fit_K"]) for r in rows if float(r["Teq_fit_K"]) > 0] or [float(r["Teq_K"]) for r in targets] or [0.0]
    cols = ["tic", "upstream_tier", "dr4_class", "Tmag", "P_ref_d", "dP_sec", "n_sec_used", "n_sec_hit", "n_transits",
            "std_MES", "despot_MES", "SDE", "Rp_fit_Re", "Teq_fit_K", "odd_even_s", "wcs_off_as", "wcs_verdict", "survived", "vetoes"]
    def _tb(rs):
        if not rs:
            return ["(none)"]
        s = " | "
        return ["| " + s.join(cols) + " |", "| " + s.join(["---"] * len(cols)) + " |"] + [
            "| " + s.join(str(r.get(c, "")) for c in cols) + " |" for r in rs
        ]
    md = [
        f"# S1-S108 Coolest Purified Platinum & Gold Candidate Audit ({el} s)",
        "",
        f"- Coolest Purified Platinum and Gold Targets Audited: {len(rows)} (Teq span {min(vt_k):.0f}-{max(vt_k):.0f} K).",
        f"- Verified S1-S108 Multi-Sector Planet Survivors (results/verified_cool_survivors.csv): {len(surv)} / {len(rows)}.",
        f"- Newly Quarantined Binary/Blend Systems Removed from Platinum/Gold Queue: {sorted(bad_tics)}.",
        "",
        "## 1. Verified Multi-Sector Cool Planet Survivors (Sectors 1-108)",
        "",
    ] + _tb(surv) + [
        "",
      "## 2. Full Audit Ledger of Top 10 Coolest Platinum & Gold Candidates",  
        "",
    ] + _tb(rows)
    (out_dir / "COOL_SURVIVORS_REPORT.md").write_text("\n".join(md) + "\n")
    cp = out_dir / "convergence_certificates.json"
    certs = json.loads(cp.read_text()) if cp.exists() else {}
    certs.update({
        "cool_vet_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "cool_vet_elapsed_s": el,
        "cool_targets_audited": len(rows),
        "cool_verified_survivors": len(surv),
        "cool_survivor_tics": [int(r["tic"]) for r in surv],
    })
    certs["files_sha256"] = {
        str(fp.relative_to(out_dir)): hashlib.sha256(fp.read_bytes()).hexdigest()
        for fp in sorted(list(out_dir.glob("*.csv")) + list(out_dir.glob("*.md")) + list(out_dir.glob("*.json")) + list(fig_d.glob("*.png")))
        if fp.name != "convergence_certificates.json"
    }
    cp.write_text(json.dumps(certs, indent=2))
    print(json.dumps({"el": el, "audited": len(rows), "survived": len(surv), "surv_tics": [r["tic"] for r in surv]}, indent=2))

if __name__ == "__main__":
    main()
    main()
