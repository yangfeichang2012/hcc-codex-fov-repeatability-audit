from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import pickle
import shutil
import tempfile
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path, PurePosixPath

import numpy as np


SEED = 20260729
BOOTSTRAPS = 2000
ME_MEMBER = "Supporting_data/microenvironment_annotations.pkl"
FEATURES = {
    "SPARC-MIF-01": "Intratumoral lymphocyte presence",
    "SPARC-MIF-02": "Tumor-endothelium interaction",
    "SPARC-MIF-03": "LA exposure level",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def safe_member(name: str) -> bool:
    member = PurePosixPath(name)
    return not member.is_absolute() and ".." not in member.parts


def verify_source(path: Path, spec: dict) -> None:
    if not path.is_file():
        raise RuntimeError(f"Missing source file: {path}")
    if path.stat().st_size != int(spec["bytes"]):
        raise RuntimeError(f"Byte-count mismatch for {path.name}")
    actual = sha256(path)
    if actual != spec["sha256"]:
        raise RuntimeError(f"SHA-256 mismatch for {path.name}: {actual}")


def download_file(spec: dict, target: Path) -> None:
    if target.exists():
        try:
            verify_source(target, spec)
            return
        except RuntimeError:
            target.unlink()
    partial = target.with_suffix(target.suffix + ".part")
    if partial.exists():
        partial.unlink()
    request = urllib.request.Request(spec["url"], headers={"User-Agent": "PRJ-2026-009-reproducibility/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    partial.replace(target)
    verify_source(target, spec)


def acquire_sources(source_dir: Path, specs: list[dict], download: bool) -> dict[str, Path]:
    source_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for spec in specs:
        path = source_dir / spec["name"]
        if download:
            download_file(spec, path)
        verify_source(path, spec)
        paths[spec["name"]] = path
    return paths


def read_metadata(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = ["region_id", "region_label", "patient_id", "state", "response", "species", "treatment"]
        if reader.fieldnames != expected:
            raise RuntimeError(f"Unexpected metadata columns: {reader.fieldnames}")
        rows = list(reader)
    if len(rows) != 48:
        raise RuntimeError(f"Expected 48 metadata rows, found {len(rows)}")
    if any(row["species"] != "Human" for row in rows):
        raise RuntimeError("Non-human row in the frozen CODEX cohort")
    if len({row["region_id"] for row in rows}) != 48:
        raise RuntimeError("Metadata region IDs are not unique")
    return rows


def extract_minimal_inputs(codex_zip: Path, resources_zip: Path, root: Path) -> tuple[list[Path], Path, list[dict]]:
    cell_dir = root / "cell_types"
    cell_dir.mkdir(parents=True)
    inventory = []
    with zipfile.ZipFile(codex_zip) as archive:
        members = archive.infolist()
        if any(not safe_member(item.filename) for item in members):
            raise RuntimeError("Unsafe path in CODEX archive")
        if any(item.flag_bits & 0x1 for item in members):
            raise RuntimeError("Encrypted member in CODEX archive")
        selected = sorted(
            (item for item in members if item.filename.endswith(".cell_types.csv")),
            key=lambda item: item.filename,
        )
        if len(selected) != 48:
            raise RuntimeError(f"Expected 48 cell-type members, found {len(selected)}")
        for item in selected:
            target = cell_dir / PurePosixPath(item.filename).name
            with archive.open(item) as source, target.open("wb") as output:
                shutil.copyfileobj(source, output)
            inventory.append({
                "region_id": target.name.removesuffix(".cell_types.csv"),
                "member": item.filename,
                "bytes": target.stat().st_size,
                "sha256": sha256(target),
            })
    me_file = root / "microenvironment_annotations.pkl"
    with zipfile.ZipFile(resources_zip) as archive:
        if ME_MEMBER not in archive.namelist():
            raise RuntimeError(f"Missing {ME_MEMBER}")
        if not safe_member(ME_MEMBER):
            raise RuntimeError("Unsafe microenvironment member path")
        with archive.open(ME_MEMBER) as source, me_file.open("wb") as output:
            shutil.copyfileobj(source, output)
    return sorted(cell_dir.glob("*.cell_types.csv")), me_file, inventory


def read_cell_types(path: Path) -> list[str]:
    labels = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["CELL_ID", "CELL_TYPE"]:
            raise RuntimeError(f"Unexpected cell-type columns in {path.name}")
        for expected_id, row in enumerate(reader, start=1):
            if int(row["CELL_ID"]) != expected_id:
                raise RuntimeError(f"Non-sequential CELL_ID in {path.name}")
            labels.append(row["CELL_TYPE"])
    return labels


def ratio(numerator: int, denominator: int) -> float:
    return float(numerator / denominator) if denominator else float("nan")


def region_features(cell_types: list[str], me_values: object) -> tuple[dict[str, float], dict[str, tuple[int, int]]]:
    mes = np.asarray(me_values, dtype=float)
    cells = np.asarray(cell_types, dtype=object)
    if len(cells) != len(mes):
        raise RuntimeError("Cell-type/ME length mismatch")
    endothelial = cells == "Endothelial cells"
    counts = {
        "SPARC-MIF-01": (int(np.isin(mes, [3.0]).sum()), int(np.isin(mes, [3.0, 4.0, 4.1]).sum())),
        "SPARC-MIF-02": (
            int((endothelial & np.isin(mes, [2.0, 4.0, 4.1])).sum()),
            int((endothelial & np.isin(mes, [2.0, 2.1, 2.2, 2.4, 2.5, 3.0, 4.0, 4.1])).sum()),
        ),
        "SPARC-MIF-03": (int(np.isin(mes, [7.0, 7.1]).sum()), int(np.isin(mes, [0.0, 7.0, 7.1]).sum())),
    }
    return {feature: ratio(*value) for feature, value in counts.items()}, counts


def icc_balanced_two(values: np.ndarray) -> tuple[float, float, float, float, bool]:
    if values.ndim != 2 or values.shape[1] != 2:
        raise RuntimeError(f"Frozen ICC requires an n x 2 matrix, got {values.shape}")
    n, k = values.shape
    patient_means = values.mean(axis=1)
    grand_mean = float(values.mean())
    ms_between = float(k * np.square(patient_means - grand_mean).sum() / (n - 1))
    ms_within = float(np.square(values - patient_means[:, None]).sum() / (n * (k - 1)))
    raw = float((ms_between - ms_within) / (ms_between + (k - 1) * ms_within))
    return max(0.0, raw), raw, ms_between, ms_within, raw < 0


def percentile(values: np.ndarray, q: float) -> float:
    return float(np.percentile(values[np.isfinite(values)], q))


def build_primary(metadata: list[dict], region_values: dict[str, dict[str, float]]) -> tuple[list[dict], list[dict], list[dict]]:
    pre = defaultdict(list)
    for row in metadata:
        if row["state"] == "Pre":
            pre[row["patient_id"]].append(row)
    candidate_ids = sorted(patient for patient, rows in pre.items() if len(rows) >= 2)
    if len(candidate_ids) != 11:
        raise RuntimeError(f"Expected 11 pretreatment paired-FOV patients, found {len(candidate_ids)}")

    primary_rows = []
    bootstrap_rows = []
    loo_rows = []
    rng = np.random.default_rng(SEED)

    for feature_id, feature_name in FEATURES.items():
        eligible = []
        for patient in candidate_ids:
            rows = sorted(pre[patient], key=lambda row: row["region_id"])
            finite_rows = [row for row in rows if math.isfinite(region_values[row["region_id"]][feature_id])]
            if len(finite_rows) >= 2:
                eligible.append((patient, finite_rows[:2]))
        if len(eligible) < 8:
            raise RuntimeError(f"Frozen N gate failed for {feature_id}")
        values = np.asarray(
            [[region_values[rows[0]["region_id"]][feature_id], region_values[rows[1]["region_id"]][feature_id]] for _, rows in eligible],
            dtype=float,
        )
        estimate, raw, ms_between, ms_within, truncated = icc_balanced_two(values)
        draws = []
        n = len(eligible)
        for replicate in range(1, BOOTSTRAPS + 1):
            indices = rng.integers(0, n, size=n)
            draw, _, _, _, _ = icc_balanced_two(values[indices, :])
            draws.append(draw)
            bootstrap_rows.append({"feature_id": feature_id, "replicate": replicate, "icc": draw})
        draw_array = np.asarray(draws, dtype=float)
        slots = {patient: f"P{index:02d}" for index, (patient, _) in enumerate(eligible, start=1)}
        loo_values = []
        for omit_index, (patient, _) in enumerate(eligible):
            keep = np.arange(n) != omit_index
            loo, _, _, _, _ = icc_balanced_two(values[keep, :])
            loo_values.append(loo)
            loo_rows.append({
                "feature_id": feature_id,
                "omitted_patient_slot": slots[patient],
                "n_patients": int(keep.sum()),
                "icc": loo,
            })
        primary_rows.append({
            "feature_id": feature_id,
            "feature_name": feature_name,
            "n_patients": n,
            "n_regions": n * 2,
            "patient_coverage": n / len(candidate_ids),
            "icc": estimate,
            "icc_raw": raw,
            "ms_between": ms_between,
            "ms_within": ms_within,
            "boundary_truncated": truncated,
            "bootstrap_replicates": BOOTSTRAPS,
            "seed": SEED,
            "bootstrap_median": percentile(draw_array, 50),
            "ci_lower": percentile(draw_array, 2.5),
            "ci_upper": percentile(draw_array, 97.5),
            "boundary_zero_fraction": float((draw_array == 0).mean()),
            "loo_min_icc": min(loo_values),
            "loo_max_icc": max(loo_values),
        })
    return primary_rows, bootstrap_rows, loo_rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not args.download and args.source_dir is None:
        parser.error("Use --download or provide --source-dir")
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "source_manifest.json").read_text(encoding="utf-8"))
    source_dir = args.source_dir or (root / ".source_cache")
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise RuntimeError(f"Output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    sources = acquire_sources(source_dir.resolve(), manifest["files"], args.download)
    metadata = read_metadata(sources["metadata.csv"])

    with tempfile.TemporaryDirectory(prefix="prj009_source_rebuild_") as temp:
        cell_files, me_file, inventory = extract_minimal_inputs(
            sources["CODEX-MIF.zip"], sources["notebook_and_resources_v2.zip"], Path(temp)
        )
        # The pinned public resource is trusted only after its archive hash has passed.
        with me_file.open("rb") as handle:
            annotations = pickle.load(handle)
        if sorted(int(key) for key in annotations) != list(range(48)):
            raise RuntimeError("ME annotation keys are not exactly 0..47")
        region_ids = [path.name.removesuffix(".cell_types.csv") for path in cell_files]
        if region_ids != sorted(row["region_id"] for row in metadata):
            raise RuntimeError("Metadata-to-cell-table region mapping mismatch")
        region_values = {}
        zero_denominators = {feature: 0 for feature in FEATURES}
        for index, path in enumerate(cell_files):
            values, counts = region_features(read_cell_types(path), annotations[index])
            region_values[region_ids[index]] = values
            for feature, (_, denominator) in counts.items():
                zero_denominators[feature] += int(denominator == 0)
        primary, bootstrap, loo = build_primary(metadata, region_values)
        preflight = {
            "schema_version": "1.0.0",
            "status": "PASS",
            "source_record": manifest["source_record"],
            "source_files": [
                {"name": spec["name"], "bytes": sources[spec["name"]].stat().st_size, "sha256": sha256(sources[spec["name"]])}
                for spec in manifest["files"]
            ],
            "cell_type_files": len(cell_files),
            "microenvironment_sha256": sha256(me_file),
            "metadata_rows": len(metadata),
            "all_region_lengths_match": all(item["bytes"] > 0 for item in inventory),
            "expression_files_extracted": 0,
            "expression_files_read": 0,
            "zero_denominator_regions": zero_denominators,
        }

    write_tsv(output_dir / "primary_icc.tsv", primary, list(primary[0]))
    write_tsv(output_dir / "bootstrap_icc_draws.tsv", bootstrap, ["feature_id", "replicate", "icc"])
    write_tsv(output_dir / "leave_one_patient_out.tsv", loo, ["feature_id", "omitted_patient_slot", "n_patients", "icc"])
    write_json(output_dir / "source_preflight.json", preflight)
    analysis_manifest = {
        "schema_version": "1.0.0",
        "project_id": "PRJ-2026-009",
        "analysis_status": "SOURCE_TO_RESULT_RECONSTRUCTION_PASS",
        "analysis_unit": "patient",
        "scope": "pretreatment within-biopsy paired-CODEX-FOV ICC-only audit",
        "authorized_features": list(FEATURES),
        "seed": SEED,
        "bootstrap_replicates_per_feature": BOOTSTRAPS,
        "primary_candidate_patients_with_two_regions": 11,
        "expression_files_read": 0,
        "patient_or_region_level_outputs_written": 0,
        "primary_results": primary,
        "prohibited_layers_executed": 0,
    }
    write_json(output_dir / "analysis_manifest.json", analysis_manifest)
    files = sorted(path for path in output_dir.iterdir() if path.is_file())
    (output_dir / "SHA256SUMS.txt").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in files), encoding="utf-8"
    )
    print(json.dumps({"status": "PASS", "output_dir": str(output_dir), "files": len(files) + 1}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
