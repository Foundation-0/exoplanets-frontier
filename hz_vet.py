#!/usr/bin/env python3
"""Issue 5 Runner: S1-S108 habitable-zone audit of the 6 bright leads + cool Tier S/C benchmarks."""
from __future__ import annotations
import hashlib, json, time
from pathlib import Path
import numpy as np, pandas as pd
from campaigns import gv
from converge_io import q_tic
from hz_core import audit_target
from mast_data import load_catalogs

ISSUE5_LEADS = [
    (181015118, 25.261970, "ISSUE5_HZ_PLATINUM"),
    (181015118, 18.317993, "ISSUE5_HZ_GOLD"),
    (246897668, 24.752432, "ISSUE5_HZ_PLATINUM"),
    (37904598, 17.049176, "ISSUE5_HZ_GOLD"),
    (146636926, 26.056445, "ISSUE5_HZ_GOLD"),
    (146636926, 23.920264, "ISSUE5_HZ_GOLD"),
]

COOL_BENCHMARKS = [
    (166848769, 7.656817, "VERIFIED_COOL_TIER_C"),
    (224279805, 22.984167, "VERIFIED_COOL_REPORTED"),
    (170788495, 23.906881, "VERIFIED_COOL_TIER_S"),
]


def scan_eumemic_origins(root: Path = Path("_eumemic")) -> list[dict]:
    hits = []
    for fp in sorted(root.rglob("*.csv")):
        if "funnel" in fp.name.lower() or "sector" in fp.name.lower():
            continue
        try:
            df = pd.read_csv(fp)
            cl = {c.strip().lower() for c in df.columns}
            if not (("tic" in cl or "tic_id" in cl) and ("p_d" in cl or "period_d" in cl or "period" in cl)):
                continue
            rel = str(fp.relative_to(root))
            for _, r in df.iterrows():
                tid, pd_ = gv(r, "TIC", "tic", "tic_id", d=-1.0), gv(r, "P_d", "period_d", "period", d=-1.0)
                if np.isfinite(tid) and tid > 0 and np.isfinite(pd_) and pd_ > 0:
                    tr = str(r.get("tier", r.get("Tier", r.get("status", "un-tiered")))).strip()
                    hits.append({"tic": int(tid), "p": float(pd_), "rel": rel, "tier": tr, "ra": gv(r, "RA", "ra"), "dec": gv(r, "Dec", "dec")})
        except Exception:
            pass
    return hits


def main(out_dir: Path = Path("results")):
    t0w = time.time()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_d = out_dir / "figures"
    fig_d.mkdir(exist_ok=True)
    kt, kp, _ = load_catalogs()
    hits = scan_eumemic_origins()
    ex = list(COOL_BENCHMARKS)
    seen = {(t, round(p, 2)) for t, p, _ in ISSUE5_LEADS + ex}
    if (out_dir / "eumemic_gaia_dr4_readiness.csv").exists():
        dfc = pd.read_csv(out_dir / "eumemic_gaia_dr4_readiness.csv")
        v3c = dfc[dfc["upstream_tier"].astype(str).isin(["S", "C", "reported"]) & (dfc["Teq_K"] <= 430) & (dfc["Rp_Re"] <= 4.5)].sort_values("Teq_K")
        for _, r in v3c.head(4).iterrows():
            key = (int(r["tic"]), round(float(r["P_d"]), 2))
            if key not in seen:
                seen.add(key)
                ex.append((int(r["tic"]), float(r["P_d"]), f"V3_COOL_{r['upstream_tier']}"))
    tgs = ISSUE5_LEADS + ex
    tmap = {int(r["ID"]): r for r in q_tic(sorted({t for t, _, _ in tgs}))}
    rows = [audit_target(tid, pc, ln, tmap.get(tid), hits, kt, kp, fig_d) for tid, pc, ln in tgs]
    pd.DataFrame(rows).to_csv(out_dir / "hz_issue5_audit.csv", index=False)
    i5 = [r for r in rows if r["lane"].startswith("ISSUE5")]
    v3r = [r for r in rows if not r["lane"].startswith("ISSUE5")]
    n1 = sum(1 for r in i5 if r["survived_hz_planet"])
    n2 = sum(1 for r in v3r if r["survived_hz_planet"])
    el = round(time.time() - t0w, 1)
    cols = ["tic", "lane", "src_files", "known_nasa", "Tmag", "P_cand_d", "n_sec_used", "n_transits", "p_rot_d",
            "std_MES", "despot_MES", "Rp_Re", "Teq_K", "wcs_off_as", "plausible_bebs", "survived_hz_planet", "vetoes"]
    def _tb(rs):
        if not rs:
            return ["(none)"]
        s = " | "
        return ["| " + s.join(cols) + " |", "| " + s.join(["---"] * len(cols)) + " |"] + [
            "| " + s.join(str(r.get(c, "")) for c in cols) + " |" for r in rs
        ]
    md = [
        f"# Issue 5 Habitable-Zone Multi-Sector Vetting Report ({el} s)",
        "",
        f"- 6 Bright Temperate Leads: {n1}/6 survived all S1-S108 transit-shape, starspot, and WCS gates.",
        f"- Verified/Cool Tier S/C Benchmarks: {n2}/{len(v3r)} survived.",
        "",
        "## 1. Issue 5 Bright Temperate Leads (TIC 181015118, 246897668, 37904598, 146636926)",
        "",
    ] + _tb(i5) + [
        "",
        "## 2. Cool Tier S/C & Rescued Comparison Worlds",
        "",
    ] + _tb(v3r)
    (out_dir / "HZ_ISSUE5_REPORT.md").write_text("\n".join(md) + "\n")
    cp = out_dir / "convergence_certificates.json"
    certs = json.loads(cp.read_text()) if cp.exists() else {}
    certs.update({
        "issue5_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "issue5_elapsed_s": el,
        "issue5_elapsed_targets": len(rows),
        "issue5_leads_survived": n1,
        "v3_cool_survived": n2,
    })
    certs["files_sha256"] = {
        str(fp.relative_to(out_dir)): hashlib.sha256(fp.read_bytes()).hexdigest()
        for fp in sorted(list(out_dir.glob("*.csv")) + list(out_dir.glob("*.md")) + list(out_dir.glob("*.json")) + list(fig_d.glob("*.png")))
        if fp.name != "convergence_certificates.json"
    }
    cp.write_text(json.dumps(certs, indent=2))
    print(json.dumps({"el": el, "i5_ok": n1, "v3_ok": n2, "tot": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
