#!/usr/bin/env python3
"""Verify and reassemble versioned release objects without replacing existing files."""
import argparse
import csv
import hashlib
import os
from pathlib import Path
import shutil
import tempfile

BLOCK = 8 * 1024 * 1024


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(BLOCK), b''):
            h.update(block)
    return h.hexdigest()


def safe_path(base, relative):
    value = Path(relative)
    if value.is_absolute() or '..' in value.parts or '\\' in relative:
        raise ValueError('Unsafe relative path in release manifest')
    base = Path(base).resolve()
    result = (base / value).resolve()
    if result == base or base not in result.parents:
        raise ValueError('Manifest destination escapes output directory')
    return result


def read_manifest(path):
    with Path(path).open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    required = {'logical_object', 'asset_name', 'part_index', 'part_count',
                'asset_bytes', 'asset_sha256', 'object_bytes', 'object_sha256',
                'relative_destination', 'release_url'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('Missing release-manifest columns or rows')
    seen = set()
    for row in rows:
        if Path(row['asset_name']).name != row['asset_name'] or '\\' in row['asset_name']:
            raise ValueError('Asset names must be basenames')
        if row['asset_name'] in seen:
            raise ValueError('Duplicate asset name')
        seen.add(row['asset_name'])
        for key in ['part_index', 'part_count', 'asset_bytes', 'object_bytes']:
            row[key] = int(row[key])
        for key in ['asset_sha256', 'object_sha256']:
            if len(row[key]) != 64 or any(c not in '0123456789abcdef' for c in row[key]):
                raise ValueError('Invalid SHA256 digest')
    return rows


def verify(path, size, sha256):
    path = Path(path)
    if not path.is_file() or path.stat().st_size != size or digest(path) != sha256:
        raise ValueError('Size/SHA256 mismatch: ' + path.name)


def reassemble(rows, assets_dir, output_dir):
    groups = {}
    for row in rows:
        groups.setdefault(row['logical_object'], []).append(row)
    outputs = []
    for name, pieces in groups.items():
        pieces.sort(key=lambda row: row['part_index'])
        first = pieces[0]
        if [row['part_index'] for row in pieces] != list(range(1, first['part_count'] + 1)):
            raise ValueError('Incomplete or duplicate parts: ' + name)
        for row in pieces:
            if any(row[k] != first[k] for k in ['part_count', 'object_bytes', 'object_sha256', 'relative_destination']):
                raise ValueError('Inconsistent object metadata: ' + name)
            verify(safe_path(assets_dir, row['asset_name']), row['asset_bytes'], row['asset_sha256'])
        target = safe_path(output_dir, first['relative_destination'])
        if target.exists():
            verify(target, first['object_bytes'], first['object_sha256'])
            print('Already verified: ' + first['relative_destination'])
            outputs.append(target)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp_name = tempfile.mkstemp(prefix='.' + target.name + '.', suffix='.partial', dir=target.parent)
        tmp = Path(tmp_name)
        try:
            with os.fdopen(fd, 'wb') as output:
                for row in pieces:
                    with safe_path(assets_dir, row['asset_name']).open('rb') as source:
                        shutil.copyfileobj(source, output, BLOCK)
                output.flush()
                os.fsync(output.fileno())
            verify(tmp, first['object_bytes'], first['object_sha256'])
            # A concurrent writer must not cause an unchecked replacement.
            if target.exists():
                verify(target, first['object_bytes'], first['object_sha256'])
            else:
                os.link(tmp, target)
            print('Verified: ' + first['relative_destination'])
            outputs.append(target)
        finally:
            tmp.unlink(missing_ok=True)
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=Path(__file__).resolve().parents[1] / 'manifests' / 'RELEASE_OBJECTS.tsv')
    parser.add_argument('--assets-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    reassemble(read_manifest(args.manifest), args.assets_dir, args.output_dir)


if __name__ == '__main__':
    main()
