#!/usr/bin/env python3
"""Offline ledger checks and a mocked provenance counterexample, not a transit replay."""
from __future__ import annotations

import ast
import csv
import hashlib
import json
from pathlib import Path
import re
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SOURCE_HEAD = "883e28a7cecde05fdb6679bf574dbe78668332e7"


def catalog_membership_control():
    """Execute the original loader with empty external catalogs and one local row."""
    tree = ast.parse((ROOT / "mast_data.py").read_text())
    loader = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "load_catalogs")

    class Frame:
        def __init__(self, rows):
            self.rows = rows

        def __setitem__(self, key, value):
            for row in self.rows:
                row[key] = value

        def iterrows(self):
            return enumerate(self.rows)

    class PandasStub:
        @staticmethod
        def notna(value):
            return value is not None

        @staticmethod
        def read_csv(path):
            with path.open() as handle:
                return Frame(list(csv.DictReader(handle)))

        @staticmethod
        def concat(frames, **kwargs):
            return Frame([row for frame in frames for row in frame.rows])

        DataFrame = staticmethod(lambda: Frame([]))

    class EmptyArchive:
        @staticmethod
        def query_criteria(**kwargs):
            return []

    class EmptyVizier:
        def __init__(self, **kwargs):
            pass

        def get_catalogs(self, catalog):
            return []

    env = {"Path": Path, "re": re, "pd": PandasStub,
           "NasaExoplanetArchive": EmptyArchive, "Vizier": EmptyVizier}
    exec(compile(ast.Module(body=[loader], type_ignores=[]), "mast_data.py:load_catalogs", "exec"), env)
    with tempfile.TemporaryDirectory() as tmp:
        local = Path(tmp)
        (local / "candidates.csv").write_text("TIC,P_d\n999999999,25.0\n")
        known, periods, _ = env["load_catalogs"](local)
    assert 999999999 in known and periods[999999999] == [25.0]
    return {"external_catalog_rows": 0, "upstream_only_tic": 999999999,
            "member_of_kt": True, "stored_period_days": 25.0,
            "meaning": "kt membership can arise solely from the upstream candidate CSV"}


def main():
    rows = list(csv.DictReader((ROOT / "results/hz_issue5_audit.csv").open()))
    leads = [r for r in rows if r["lane"].startswith("ISSUE5")]
    comparisons = [r for r in rows if not r["lane"].startswith("ISSUE5")]
    certificate = json.loads((ROOT / "results/convergence_certificates.json").read_text())
    checks = {name: hashlib.sha256((ROOT / "results" / name).read_bytes()).hexdigest() == digest
              for name, digest in certificate["files_sha256"].items()}
    assert all(checks.values()), "Historical result bytes do not match their manifest"
    assert len(leads) == 6 and len({r["tic"] for r in leads}) == 4
    assert len(comparisons) == 5 and all(r["survived_hz_planet"] == "False" for r in rows)
    assert certificate["issue5_leads_survived"] == 0
    assert certificate["v3_cool_survived"] == 0
    files = ["hz_core.py", "hz_vet.py", "mast_data.py", "vetter.py", "converge_io.py",
             "results/hz_issue5_audit.csv", "results/convergence_certificates.json"]
    report = {
        "reviewed_source_head": SOURCE_HEAD,
        "evidence_type": "offline aggregate and byte checks plus mocked original-loader control",
        "transit_pipeline_rerun": False,
        "result_hashes_checked": len(checks), "all_result_hashes_match": all(checks.values()),
        "lead_signal_rows": len(leads), "lead_unique_stars": len({r["tic"] for r in leads}),
        "lead_survivors": 0, "comparison_signal_rows": len(comparisons), "comparison_survivors": 0,
        "reported_loaded_sectors_range": [min(int(r["n_sec_used"]) for r in rows), max(int(r["n_sec_used"]) for r in rows)],
        "comparison_teq_range_kelvin": [min(float(r["Teq_K"]) for r in comparisons), max(float(r["Teq_K"]) for r in comparisons)],
        "catalog_membership_control": catalog_membership_control(),
        "inputs_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in files},
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
