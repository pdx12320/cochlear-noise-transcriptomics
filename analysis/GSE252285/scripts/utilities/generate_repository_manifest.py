#!/usr/bin/env python3
"""Generate deterministic file inventory and SHA256 manifest for this repository."""

from __future__ import annotations

import csv
import hashlib
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_DIR = REPO_ROOT / "manifests"
SHA_PATH = MANIFEST_DIR / "SHA256SUMS.tsv"
INCLUDED_PATH = MANIFEST_DIR / "INCLUDED_FILES.tsv"
SKIP_PARTS = {".git", ".local", "__pycache__", ".pytest_cache"}


def included_files() -> list[Path]:
    files: list[Path] = []
    for path in REPO_ROOT.rglob("*"):
        if not path.is_file() or any(part in SKIP_PARTS for part in path.parts):
            continue
        if path == SHA_PATH:
            continue
        files.append(path)
    return sorted(files, key=lambda path: path.relative_to(REPO_ROOT).as_posix())


def suffix(path: Path) -> str:
    name = path.name.lower()
    for compound in (".tsv.gz", ".tar.gz", ".fastq.gz"):
        if name.endswith(compound):
            return compound
    return path.suffix.lower() or "[none]"


def write_inventory(files: list[Path]) -> None:
    summary: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
    for path in files:
        relative = path.relative_to(REPO_ROOT)
        top_level = relative.parts[0] if len(relative.parts) > 1 else "[root]"
        key = (top_level, suffix(path))
        summary[key][0] += 1
        summary[key][1] += path.stat().st_size

    with INCLUDED_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("top_level", "extension", "file_count", "total_bytes"))
        for (top_level, extension), (count, total_bytes) in sorted(summary.items()):
            writer.writerow((top_level, extension, count, total_bytes))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_sha_manifest(files: list[Path]) -> None:
    with SHA_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sha256", "bytes", "relative_path"))
        for path in files:
            writer.writerow(
                (sha256(path), path.stat().st_size, path.relative_to(REPO_ROOT).as_posix())
            )


def main() -> None:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    files = included_files()
    write_inventory(files)
    # Re-scan because INCLUDED_FILES.tsv was just created or updated.
    write_sha_manifest(included_files())
    print(f"MANIFEST_FILES={sum(1 for _ in SHA_PATH.open(encoding='utf-8')) - 1}")


if __name__ == "__main__":
    main()
