from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


EXPECTED = {
    "SPARC-MIF-01": {"n": 11, "icc": 0.9519634455573124, "lo": 0.8634736002632606, "hi": 0.9924102622364553, "loo_min": 0.942, "loo_max": 0.975},
    "SPARC-MIF-02": {"n": 11, "icc": 0.8949273737591393, "lo": 0.710416975382867, "hi": 0.9720043828901788, "loo_min": 0.878, "loo_max": 0.954},
    "SPARC-MIF-03": {"n": 9, "icc": 0.7003021676873524, "lo": 0.0, "hi": 0.9902665804577716, "loo_min": 0.594, "loo_max": 0.976},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def close(a: float, b: float, tolerance: float = 1e-12) -> bool:
    return abs(a - b) <= tolerance


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-dir", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.output_dir.resolve()
    expected_dir = (args.expected_dir or (root / "expected")).resolve()
    errors = []

    with (output / "primary_icc.tsv").open("r", encoding="utf-8", newline="") as handle:
        rows = {row["feature_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if set(rows) != set(EXPECTED):
        errors.append("primary feature set mismatch")
    for feature, expected in EXPECTED.items():
        row = rows.get(feature, {})
        if not row:
            continue
        checks = {
            "n": int(row["n_patients"]) == expected["n"],
            "icc": close(float(row["icc"]), expected["icc"]),
            "ci_lower": close(float(row["ci_lower"]), expected["lo"]),
            "ci_upper": close(float(row["ci_upper"]), expected["hi"]),
            "loo_min_rounded": round(float(row["loo_min_icc"]), 3) == expected["loo_min"],
            "loo_max_rounded": round(float(row["loo_max_icc"]), 3) == expected["loo_max"],
        }
        errors.extend(f"{feature}:{name}" for name, passed in checks.items() if not passed)

    with (output / "bootstrap_icc_draws.tsv").open("r", encoding="utf-8", newline="") as handle:
        draws = list(csv.DictReader(handle, delimiter="\t"))
    if len(draws) != 6000:
        errors.append(f"bootstrap row count {len(draws)} != 6000")
    for feature in EXPECTED:
        feature_rows = [row for row in draws if row["feature_id"] == feature]
        if len(feature_rows) != 2000 or [int(row["replicate"]) for row in feature_rows] != list(range(1, 2001)):
            errors.append(f"{feature}:bootstrap sequence")

    preflight = json.loads((output / "source_preflight.json").read_text(encoding="utf-8"))
    if preflight.get("status") != "PASS" or preflight.get("expression_files_read") != 0:
        errors.append("source preflight boundary")
    manifest = json.loads((output / "analysis_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("analysis_status") != "SOURCE_TO_RESULT_RECONSTRUCTION_PASS":
        errors.append("analysis manifest status")
    if manifest.get("patient_or_region_level_outputs_written") != 0:
        errors.append("output boundary")

    expected_files = sorted(path for path in expected_dir.iterdir() if path.is_file())
    for expected_path in expected_files:
        actual = output / expected_path.name
        if not actual.is_file() or sha256(actual) != sha256(expected_path):
            errors.append(f"hash mismatch: {expected_path.name}")

    report = {
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "numerical_claims_checked": 18,
        "bootstrap_rows_checked": len(draws),
        "expected_files_checked": len(expected_files),
    }
    print(json.dumps(report, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
