# GSE252285 cluster 19 annotation sensitivity

This review-requested derived analysis changes only cluster 19 (319 cells) from UNRESOLVED to Fibrocyte_mesenchymal. The original selected 36,111-cell object and primary labels are preserved. The original eleven mean-expression localization assignments all remain; Cebpb becomes additionally recurrent only in the derived sensitivity. Apold1 retains the highest tested broad-category mean and detection fraction in all libraries, while cluster 19 itself has strong Apold1 expression. The data do not support endothelial-specific expression.

All six originally tested program-cell-type pairs retain their significance classifications in both matching modes: five significant and MAPK/AP1-to-endothelial nonsignificant. The fibrocyte TNFA/NFKB delta changes from 0.0379826 to 0.0393930. Primary MAPK q changes from 0.0723928 to 0.0724928; its gene-length-matched q remains 0.0816918. Complete statistics are retained, including nonsignificant and changed values.

## Contents

- `scripts/`: executed parametric summary/sensitivity scripts and the original portable matched-null implementation for ranking-cache preparation.
- `results/`: compact expression/UCell comparisons, cluster 19 evidence and independent validations.
- `null_final/`: compact matched-null results, matching diagnostics and runtime records.
- `validation/`: input hashes, source/public hashes, code-adaptation and resource records, independent checks, figure QA and the release-asset digest.
- `runtime/`: publication copies of executed logs; only machine-local storage paths were sanitized in logs and session records. Original records are retained separately. Each source/public digest and transformation count is recorded.
- `EXPRESSION_CONCLUSIONS.md`, `FULL_NULL_CONCLUSIONS.md`: full interpretation and limits.

Download the [v1.1.0 companion asset](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/download/v1.1.0/GSE252285_cluster19_annotation_sensitivity_full_records_20261101.tar.gz). Its size is 12,166,547 bytes and its SHA256 is `3fc12921550f0e7adc70588baa464eca4b019fa44eb775153e07847545117b6f`. Verify the downloaded file before extraction:

```sh
python3 - <<'PYVERIFY'
from pathlib import Path
import hashlib
p = Path("GSE252285_cluster19_annotation_sensitivity_full_records_20261101.tar.gz")
assert p.stat().st_size == 12166547
assert hashlib.sha256(p.read_bytes()).hexdigest() == "3fc12921550f0e7adc70588baa464eca4b019fa44eb775153e07847545117b6f"
print("Size and SHA256 verified")
PYVERIFY
```

The unchanged base selected Seurat object is provided by the [v1.0.0 release](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/tag/v1.0.0), in numbered parts with a reassembly/checksum script. Its complete-object SHA256 is `24909e74b584dcb33a113d815234cda0861ed6fea07f78cf66fd99d5324d76b8`.

The companion `GSE252285_cluster19_annotation_sensitivity_full_records_20261101.tar.gz` supplies both full null distributions, all eight program/matching-mode checkpoints, the original cached six-program scores, cell-level metadata, compact t-SNE/expression/UCell plotting data and detectable-gene covariates. Its internal `ASSET_CONTENTS.tsv` records every member hash. Extract it into this directory to restore the documented relative layout. Scientific data are byte-identical copies; the asset contains no private invocation logs.

## Reproduce

Use the original recorded R 4.4.3 / UCell 2.10.1 environment and the other package versions in `null_final/sessionInfo.txt`. No dependency was installed for this sensitivity. The original UCell user library was used; a pre-existing library containing the data-only GenomeInfoDbData dependency was appended after the original runtime libraries. If reproducing in that arrangement, supply `R_LIBS_USER` and optional `EXISTING_DATA_LIBRARY` explicitly. These environment paths are local configuration, not scientific parameters.

Required source inputs are identified by basename, size and SHA256 in `validation/INPUT_HASHES.json`:

1. The selected object from the root repository release (SHA256 `24909e74b584dcb33a113d815234cda0861ed6fea07f78cf66fd99d5324d76b8`).
2. The cached scores included in the companion asset (`inputs/UCELL_SCORES_PER_CELL.tsv.gz`).
3. The same frozen Hallmark/Reactome GMT files under the root repository's `analysis/GSE252285/complete_results/`, and GENCODE vM23 GTF documented in `analysis/GSE252285/portable_scripts/README.md`. Their original licenses and notices continue to apply; this bundle does not relicense third-party references.
4. The original UCell ranking cache, or rankings rebuilt with the original method as below.
5. The original matched-null outputs under the root repository's `analysis/GSE252285/complete_results/` for the numerical reproduction check.

Use a new output directory for each invocation. The scripts reject an existing destination.

```sh
Rscript --vanilla scripts/summarize_annotation_sensitivity.R \
  ORIGINAL_OBJECT_RDS inputs/UCELL_SCORES_PER_CELL.tsv.gz NEW_SUMMARY_DIRECTORY

# Only if the original ranking cache is unavailable: the original 100-iteration
# method prepares work/ucell_ranks_maxRank1500.rds in a new run directory.
N_ITER=100 Rscript --vanilla scripts/original_method_run_matched_null.R \
  ORIGINAL_OBJECT_RDS FROZEN_GMT_DIRECTORY GENCODE_VM23_GTF_GZ NEW_CACHE_PREPARATION_DIRECTORY

N_ITER=100 Rscript --vanilla scripts/run_matched_null_annotation_sensitivity.R \
  ORIGINAL_OBJECT_RDS FROZEN_GMT_DIRECTORY GENCODE_VM23_GTF_GZ \
  NEW_BENCHMARK_DIRECTORY RANK_CACHE_RDS ORIGINAL_MATCHED_NULL_OUTPUT_DIRECTORY

N_ITER=10000 Rscript --vanilla scripts/run_matched_null_annotation_sensitivity.R \
  ORIGINAL_OBJECT_RDS FROZEN_GMT_DIRECTORY GENCODE_VM23_GTF_GZ \
  NEW_FINAL_DIRECTORY RANK_CACHE_RDS ORIGINAL_MATCHED_NULL_OUTPUT_DIRECTORY
```

The executed full sensitivity used one OpenMP/BLAS thread, a 24-GiB address-space cap and a 2,400-second timeout; it took 1,069.85 seconds and peaked at 7.24 GiB. The separately supplied original cache-preparation implementation retains its original ranking parallelism. The existing cache was reused in this run; cache rebuilding was not rerun as part of this sensitivity.

The benchmark is a resource and equivalence check, not a scientific result. The final run retains the original 10,000 draws per program/matching mode, seed 20260823, expression/detection strata, gene-length sensitivity, UCell maxRank 1500, average ties, six original test pairs and per-mode six-test BH correction. MYC/E2F are not added to the null family. Paired original-label evaluation reproduced all 120,000 archived null statistics to within 5.55e-17; independent empirical P/BH calculations agreed within 1.1e-19.
