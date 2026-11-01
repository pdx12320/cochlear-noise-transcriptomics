#!/usr/bin/env python3
"""Download SHA256-verified v1.0.0 assets and reassemble the archived analysis objects."""
import argparse
import os
from pathlib import Path
import shutil
import tempfile
import urllib.request
from reassemble_release_objects import BLOCK, read_manifest, reassemble, safe_path, verify


def fetch(row, assets_dir):
    target = safe_path(assets_dir, row['asset_name'])
    if target.exists():
        verify(target, row['asset_bytes'], row['asset_sha256'])
        print('Already downloaded: ' + row['asset_name'])
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix='.' + target.name + '.', suffix='.partial', dir=target.parent)
    tmp = Path(tmp_name)
    try:
        request = urllib.request.Request(row['release_url'], headers={'User-Agent': 'cochlear-noise-transcriptomics-archive/1.0.0'})
        with urllib.request.urlopen(request, timeout=60) as response, os.fdopen(fd, 'wb') as output:
            fd = None
            shutil.copyfileobj(response, output, BLOCK)
            output.flush()
            os.fsync(output.fileno())
        verify(tmp, row['asset_bytes'], row['asset_sha256'])
        if target.exists():
            verify(target, row['asset_bytes'], row['asset_sha256'])
        else:
            os.link(tmp, target)
        print('Downloaded and verified: ' + row['asset_name'])
    finally:
        if fd is not None:
            os.close(fd)
        tmp.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=Path(__file__).resolve().parents[1] / 'manifests' / 'RELEASE_OBJECTS.tsv')
    parser.add_argument('--assets-dir', type=Path, default=Path('data/release_downloads'))
    parser.add_argument('--output-dir', type=Path, default=Path('data/objects'))
    parser.add_argument('--download-only', action='store_true')
    args = parser.parse_args()
    rows = read_manifest(args.manifest)
    for row in rows:
        fetch(row, args.assets_dir)
    if not args.download_only:
        reassemble(rows, args.assets_dir, args.output_dir)
    print('The supporting .tar.gz is verified but not automatically extracted; extract it into a new work directory.')


if __name__ == '__main__':
    main()
