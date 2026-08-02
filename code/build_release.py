from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path


FIXED_TIME = (2026, 8, 2, 0, 0, 0)
EXCLUDED_PARTS = {".git", ".source_cache", "reproduced_outputs", "__pycache__"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def included(root: Path, path: Path) -> bool:
    relative = path.relative_to(root)
    return not any(part in EXCLUDED_PARTS for part in relative.parts) and path.suffix != ".pyc" and path.name != "SHA256SUMS.json"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    files = sorted(path for path in root.rglob("*") if path.is_file() and included(root, path))
    manifest = {
        "schema_version": "1.0.0",
        "release": "v1.0.0",
        "files": [
            {"path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in files
        ],
    }
    manifest_path = root / "SHA256SUMS.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    files.append(manifest_path)
    files.sort()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in files:
            name = path.relative_to(root).as_posix()
            info = zipfile.ZipInfo(name, FIXED_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            archive.writestr(info, path.read_bytes(), compresslevel=9)
    report = {"status": "PASS", "members": len(files), "bytes": output.stat().st_size, "sha256": sha256(output), "output": str(output)}
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
