# Manifests

- `CURRENT_ANALYSIS_PARAMETERS.tsv`: fixed parameters for the selected 36,111-cell branch.
- `CURRENT_RELEASE_STATUS.tsv`: current cell, cluster, annotation and condition-mapping status.
- `SHA256SUMS.tsv`: SHA256 and byte size for every Git candidate except itself.
- `INCLUDED_FILES.tsv`: module/extension inventory.
- `EXCLUDED_FILES.tsv`: raw, large, private and superseded-result exclusions.
- `SOURCE_PROVENANCE.tsv`: provenance of the selected modules.
- `REPRODUCIBILITY_CHECK.tsv`: repository validation results.
- `source/`: executed resolved-condition parameters, input-object hashes, output manifests and consistency checks. Machine paths use `<GSE252285_PROJECT_ROOT>`.

The SHA256 file cannot include its own digest because that would be recursive.
