#!/usr/bin/env python3
"""Maintainer utility: hash the staged repository files, excluding the hash manifest."""
import csv
import hashlib
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parents[1]
manifest = root/'manifests/REPOSITORY_FILES.tsv'
files = subprocess.check_output(['git','ls-files','-z'], cwd=root).decode().split('\0')
rows = []
for name in sorted(x for x in files if x and x != 'manifests/REPOSITORY_FILES.tsv'):
    path = root/name
    if not path.is_file():
        raise SystemExit('Stage the completed changes first; missing tracked file: '+name)
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    rows.append({'path': name, 'bytes': path.stat().st_size, 'sha256': h.hexdigest()})
manifest.parent.mkdir(exist_ok=True)
with manifest.open('w',newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=['path','bytes','sha256'], delimiter='\t',lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
print(f'Hashed {len(rows)} staged file paths; stage {manifest.relative_to(root)} before committing.')
