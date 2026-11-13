#!/usr/bin/env python3
"""Purge validation/inverted_lightcurves_vetted.csv and rebuild clean Pre-Gaia-DR4 candidate tables."""
from __future__ import annotations
import hashlib, json, time
from pathlib import Path
import pandas as pd
from dr4_core import build_665_catalog


def main(out_dir: Path = Path("results")):
    t0 = time.time()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_d = out_dir / "figures"
    fig_d.mkdir(exist_ok=True)
    cat_rows, _ = build_665_catalog(None, out_dir)
    assert len(cat_rows) > 0, "No genuine candidates loaded"
    bad_src = [r["src_file"] for r in cat_rows if any(w in r["src_file"].lower() for w in ("validation", "inverted", "false", "reject"))]
    assert not bad_src, f"Contaminated rows detected: {bad_src[:5]}"
    leaked_i5 = sorted({181015118, 246897668, 37904598, 146636926} & {int(r["tic"]) for r in cat_rows})
    plat = [r for r in cat_rows if r["dr4_class"] in ("PLATINUM_SOLO_VERIFIED", "GOLD_LOW_BLEND")]
    hz_clean = [r for r in plat if r["is_hz_temperate"]]
    hz_all = [r for r in cat_rows if r["is_hz_temperate"]]
    dfc = pd.DataFrame(cat_rows)
    dfc.to_csv(out_dir / "eumemic_gaia_dr4_readiness.csv", index=False)
    pd.DataFrame(plat).to_csv(out_dir / "platinum_followup_targets.csv", index=False)
    el = round(time.time() - t0, 1)
    cls_c = dfc["dr4_class"].value_counts().to_dict()
    tier_c = dfc["upstream_tier"].value_counts().to_dict()
    src_c = dfc["src_file"].value_counts().to_dict()
    cols = ["tic", "upstream_tier", "gaia_dr3_id", "P_d", "Rp_Re", "Teq_K", "is_hz_temperate", "TSM", "Tmag", "ruwe", "plausible_beb_65as", "frontier_note", "dr4_class"]
    s = " | "
    top_rows = hz_clean + [r for r in plat if not r["is_hz_temperate"]]
    tb = ["| " + s.join(cols) + " |", "| " + s.join(["---"] * len(cols)) + " |"] + [
        "| " + s.join(str(r.get(c, "")) for c in cols) + " |" for r in top_rows[:30]
    ]
    rep = [
        f"# Purified Pre-Gaia-DR4 Candidate Readiness Audit ({el} s)",
        "",
        f"- Genuine Planet Candidates Indexed (results/eumemic_gaia_dr4_readiness.csv): {len(cat_rows)} (0 rows from validation/inverted_lightcurves_vetted.csv; Issue 5 TIC overlap: {leaked_i5}).",
        f"- Source Provenance Breakdown: {src_c}",
        f"- Upstream Tier Breakdown: {tier_c}",
        f"- Gaia DR3 and Multi-Sector Readiness Classes: {cls_c}",
        f"- Clean Follow-Up Queue (results/platinum_followup_targets.csv): {len(plat)} total (PLATINUM_SOLO_VERIFIED + GOLD_LOW_BLEND), including {len(hz_clean)} clean temperate 180-325 K candidates ({len(hz_all)} temperate across all classes).",
        "",
        "## Top Purified Platinum and Gold Candidates for Gaia DR4 (2 Dec 2026)",
        "",
    ] + tb
    (out_dir / "GAIA_DR4_CATALOG_AUDIT.md").write_text("\n".join(rep) + "\n")
    cr_p = out_dir / "CONVERGENCE_REPORT.md"
    if cr_p.exists():
        old_lines = cr_p.read_text().splitlines()
        idx2 = next((i for i, ln in enumerate(old_lines) if ln.startswith("## 2.")), len(old_lines))
        head = [
            "# Phase 3 Synthesis (Purified): Genuine-Only Gaia DR4 Catalog, Cycle 1-8 Ephemerides & Two-Pass Cool-Dwarf Survey",
            "",
            f"Purified results_public/candidates_allsky.csv + VizieR I/355/gaiadr3 catalog rebuild ({el} s; validation/inverted_lightcurves_vetted.csv purged).",
            "",
            f"- Genuine Candidate Catalog (results/eumemic_gaia_dr4_readiness.csv): {len(cat_rows)} candidates indexed ({cls_c}); {len(plat)} Platinum/Gold Clean Targets ({len(hz_clean)} temperate 180-325 K).",
            "",
            "## 1. Top Purified Platinum Solo-Verified & Gold Low-Blend Candidates",
            "",
        ] + tb + [""]
        cr_p.write_text("\n".join(head + old_lines[idx2:]) + "\n")
    cp = out_dir / "convergence_certificates.json"
    certs = json.loads(cp.read_text()) if cp.exists() else {}
    certs.update({
        "dr4_purified_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "dr4_purified_elapsed_s": el,
        "catalog_total": len(cat_rows),
        "platinum_gold_count": len(plat),
        "temperate_hz_clean_count": len(hz_clean),
        "temperate_hz_all_count": len(hz_all),
        "inverted_control_rows_remaining": 0,
        "dr4_class_counts": cls_c,
        "upstream_tier_counts": tier_c,
    })
    certs["files_sha256"] = {
        str(fp.relative_to(out_dir)): hashlib.sha256(fp.read_bytes()).hexdigest()
        for fp in sorted(list(out_dir.glob("*.csv")) + list(out_dir.glob("*.md")) + list(out_dir.glob("*.json")) + list(fig_d.glob("*.png")))
        if fp.name != "convergence_certificates.json"
    }
    cp.write_text(json.dumps(certs, indent=2))
    print(json.dumps({"el": el, "cat": len(cat_rows), "plat": len(plat), "hz_clean": len(hz_clean), "cls": cls_c, "tiers": tier_c}, indent=2))


if __name__ == "__main__":
    main()
