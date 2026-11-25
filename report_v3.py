from __future__ import annotations
import hashlib, json, time
from pathlib import Path
import pandas as pd

def _tbl(rows: list, cols: list) -> list[str]:
    if not rows: return ["(none)"]
    s = " | "
    out = ["| " + s.join(cols) + " |", "| " + s.join(["---"] * len(cols)) + " |"]
    for r in rows:
        out.append("| " + s.join(str(r.get(c, "")) for c in cols) + " |")
    return out

def canonical_files_digest(files_sha256: dict[str, str]) -> str:
    """Deterministic SHA-256 over canonical sorted filename->sha256 mapping."""
    canon = json.dumps(
        {str(k): str(v) for k, v in sorted(files_sha256.items())},
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )
    return hashlib.sha256(canon.encode("utf-8")).hexdigest()


def save_phase3_outputs(out_dir: Path, fig_d: Path, t0w: float, s_max: int, bench: dict, cat_rows: list, plat: list, ttv_r: list, ep_r: list, sib_r: list, blend_r: list, md_r: list, vpass: list):
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_d.mkdir(parents=True, exist_ok=True)
    for nm, rows in [("ttv_ephemeris_audit.csv", ttv_r), ("ttv_transit_epochs.csv", ep_r), ("sibling_residual_search.csv", sib_r), ("blend_localization.csv", blend_r), ("despotted_mdwarf_audit.csv", md_r), ("vetted_new_candidates.csv", vpass)]:
        pd.DataFrame(rows).to_csv(out_dir / nm, index=False)
    el_s = round(time.time() - t0w, 1)
    cls_raw = pd.DataFrame(cat_rows)["dr4_class"].value_counts().to_dict() if cat_rows else {}
    cls_c = {str(k): int(v) for k, v in sorted(cls_raw.items())}
    hz_p = [r for r in plat if r["is_hz_temperate"]]
    top_p = (hz_p + [x for x in plat if not x["is_hz_temperate"]])[:25]
    top_m = sorted(md_r, key=lambda x: (bool(x.get("passed")), float(x.get("MES", 0.0))), reverse=True)[:22]
    md = ["# Phase 3 Synthesis: 665-Candidate Gaia DR4 Catalog & Cycle 1-8 Ephemerides", "",
          f"- Catalog: {len(cat_rows)} indexed ({cls_c}), {len(plat)} Platinum/Gold ({len(hz_p)} temperate HZ 180-325 K).",
          f"- TOI-540b Benchmark: P={bench['despotted_P_d']} d (MES={bench['despotted_MES']}), Null={bench['null_false_alarms_ge_9p5']}/{bench['null_trials']} | Cool-Dwarfs ({len(md_r)}): {len(vpass)} passed.", "",
          "## 1. Top Platinum Solo-Verified & Temperate Candidates", ""] + _tbl(top_p, ["tic", "gaia_dr3_id", "upstream_tier", "P_d", "Rp_Re", "Teq_K", "is_hz_temperate", "TSM", "Tmag", "ruwe", "plausible_beb_65as", "dr4_class"]) + [
          "", "## 2. Full 2018-2026 (Cycles 1-8) Ephemeris & TTV Lock", ""] + _tbl(ttv_r, ["tic", "upstream_tier", "n_sec_avail", "n_sec_used", "n_epochs_fit", "P_orig_d", "P_ref_d", "dP_sec", "rms0_min", "rms_lin_min", "dBIC_lin", "dBIC_ttv", "snr_stale", "snr_refined", "class"]) + [
          "", "## 3. VizieR Gaia DR3 + WCS Sub-Pixel Blend Audit", ""] + _tbl(blend_r, ["tic", "upstream_tier", "P_d", "gaia_sid", "ruwe", "pm_masyr", "n_sec", "off_target_as", "err_as", "plausible_nbs", "closest_nb_sep_as", "verdict"]) + [
          "", "## 4. Two-Pass De-Spotted Cool-Dwarf Survey", ""] + _tbl(top_m, ["tic", "lane", "p_rot_d", "raw_ppm", "clean_ppm", "P_d", "MES", "passed", "vetoes"])
    (out_dir / "CONVERGENCE_REPORT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    artifact_paths = sorted(
        (
            fp for fp in list(out_dir.glob("*.csv")) + list(out_dir.glob("*.md")) + list(out_dir.glob("*.json")) + list(fig_d.glob("*.png"))
            if fp.name != "convergence_certificates.json" and fp.is_file()
        ),
        key=lambda fp: fp.relative_to(out_dir).as_posix(),
    )
    files_sha256 = {
        fp.relative_to(out_dir).as_posix(): hashlib.sha256(fp.read_bytes()).hexdigest()
        for fp in artifact_paths
    }
    certs = {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "elapsed_s": el_s, "s_max": s_max,
             "catalog_total": len(cat_rows), "platinum_gold_count": len(plat), "temperate_hz_clean_count": len(hz_p),
             "new_mdwarf_passed": len(vpass), "benchmark": bench,
             "manifest_sha256": canonical_files_digest(files_sha256),
             "files_sha256": files_sha256}
    (out_dir / "convergence_certificates.json").write_text(
        json.dumps(certs, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"cat": len(cat_rows), "plat": len(plat), "hz": len(hz_p), "ttv": len(ttv_r), "blend": len(blend_r), "mdwarf": len(md_r), "new_pass": len(vpass)}, indent=2, sort_keys=True))
