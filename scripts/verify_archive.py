#!/usr/bin/env python3
"""Verify the versioned repository payload using only Python's standard library."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

with (ROOT/'manifests/REPOSITORY_FILES.tsv').open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
assert rows, 'Empty repository manifest'
recorded = {row['path'] for row in rows}
assert len(recorded) == len(rows), 'Duplicate repository manifest path'
actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
          if p.is_file() and p != ROOT/'manifests/REPOSITORY_FILES.tsv'
          and '.git' not in p.relative_to(ROOT).parts
          and '__pycache__' not in p.parts and p.suffix != '.pyc'}
assert recorded == actual, f'Incomplete repository manifest: missing={sorted(actual-recorded)}, extra={sorted(recorded-actual)}'
for row in rows:
    rel = Path(row['path'])
    assert not rel.is_absolute() and '..' not in rel.parts, 'Unsafe manifest path'
    path = ROOT/rel
    assert path.is_file(), f'Missing file: {rel}'
    assert path.stat().st_size == int(row['bytes']), f'Size mismatch: {rel}'
    assert digest(path) == row['sha256'], f'Checksum mismatch: {rel}'

counts_path = ROOT/'analysis/GSE252285/tables/condition_descriptive/CONDITION_TOTAL_CELLS.tsv'
with counts_path.open() as handle:
    counts = {r['condition']: int(r['total_cells']) for r in csv.DictReader(handle, delimiter='\t')}
assert counts == {'Control': 10273, 'HHL': 16631, 'NIHL': 9207}, counts

required = ['README.md','CITATION.cff','LICENSE','THIRD_PARTY_NOTICES.md',
            'manuscript/main_scientific_reports.tex','manuscript/main_scientific_reports.pdf',
            'manuscript/supplementary.pdf','data/README.md','docs/REPRODUCING.md',
            'docs/VALIDATION.md','analysis/GSE160639/README.md',
            'analysis/GSE160639/deseq2/scripts/authorized_bulk_deseq2_reestimate.R',
            'analysis/GSE160639/deseq2/environment/CONDA_EXPLICIT.txt',
            'analysis/GSE252285/inputs/annotation_reference.tsv']
for name in required:
    assert (ROOT/name).is_file(), f'Missing entry point: {name}'

bulk = ROOT/'analysis/GSE160639'
with (bulk/'deseq2/results/NEW_DESEQ2_SUMMARY.tsv').open() as handle:
    de_rows = list(csv.DictReader(handle, delimiter='\t'))
expected_de = {'94dB_vs_70dB': 5777, '105dB_vs_70dB': 889,
               '105dB_vs_94dB': 4177, 'overall_dose_LRT': 6398}
assert len(de_rows) == 4
assert {r['comparison']: int(r['significant']) for r in de_rows} == expected_de
assert all(int(r['genes']) == 19673 for r in de_rows)
with (bulk/'enrichment/results/TARGETED_27_RESULTS.tsv').open() as handle:
    pathways = list(csv.DictReader(handle, delimiter='\t'))
assert len(pathways) == 27
assert len({(r['contrast'], r['pathway']) for r in pathways}) == 27
assert len({r['pathway'] for r in pathways}) == 9
positive = {contrast: sum(r['contrast'] == contrast and
            float(r['global_27_test_FDR']) < 0.05 for r in pathways)
            for contrast in expected_de if contrast != 'overall_dose_LRT'}
assert positive == {'94dB_vs_70dB': 8, '105dB_vs_70dB': 9, '105dB_vs_94dB': 0}

with (ROOT/'analysis/GSE168041/sensitivity/summary/LOO_FOCUSED_ALL_56.tsv').open() as handle:
    focused = list(csv.DictReader(handle, delimiter='\t'))
assert len(focused) == 56
assert all(r['same_direction_as_full'].lower() == 'true' for r in focused)
assert sum(float(r['LOO_FDR']) < .05 for r in focused) == 53
with (ROOT/'analysis/GSE168041/sensitivity/summary/ALL_48_FIT_RECORDS.tsv').open() as handle:
    assert len(list(csv.DictReader(handle, delimiter='\t'))) == 48
with (ROOT/'analysis/GSE252285/annotation_sensitivity/null_final/C19_VS_ORIGINAL_MATCHED_NULL.tsv').open() as handle:
    annotation = list(csv.DictReader(handle, delimiter='\t'))
assert len(annotation) == 12
assert all(r['significant_FDR05_original'] == r['significant_FDR05_c19_revised'] for r in annotation)
assert sum(float(r['empirical_FDR_c19_revised']) < .05 for r in annotation) == 10
assert len(list((ROOT/'manuscript/figures').glob('*.pdf'))) == 14
print(json.dumps({'repository_files_verified': len(rows), 'GSE252285_cells': sum(counts.values()),
                  'GSE160639_genes': 19673, 'GSE168041_LOO_fits': 48, 'GSE168041_focused_estimates': 56, 'GSE252285_annotation_tests': 12, 'GSE160639_pathway_tests': len(pathways),
                  'GSE160639_global_FDR05_by_contrast': positive,
                  'status': 'PASS', 'scope': 'Archive checksums and declared structure; no upstream analysis is executed.'}, indent=2))
