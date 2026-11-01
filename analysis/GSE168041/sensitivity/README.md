# Animal-level robustness of the GSE168041 results

This compact repository module contains the executed scripts, exact runtime records, all seven focused pathway/model summaries, per-animal expression profiles and validation records. The primary four-versus-four-animal results remain unchanged. All 56 focused deletion estimates retained their direction; fibrocyte TNF–NF-kB remained below FDR 0.05 in 7/8 SGN and 6/8 lateral-wall deletions. The other five focused models remained below FDR 0.05 in all eight deletions. Full results and exceptions are in `summary/`.

The full frozen inputs and all 48 complete DE/GSEA outputs are distributed in the Release asset **`GSE168041_LOO_Sensitivity_20261101.tar.gz`** (67,555,809 bytes), SHA-256:

`e5aee91475f6932e23bf66eae16f5842998d18aef6d951f454fa0b1b42f6f811`

Download and extract the versioned release asset:

```sh
curl -fL -o GSE168041_LOO_Sensitivity_20261101.tar.gz \
  https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/download/v1.1.0/GSE168041_LOO_Sensitivity_20261101.tar.gz
sha256sum -c RELEASE_ASSET_SHA256SUMS.txt
mkdir -p extracted
tar -xzf GSE168041_LOO_Sensitivity_20261101.tar.gz -C extracted
```

The archive already contains the top-level directory `GSE168041_LOO_Sensitivity_20261101/`. Its README describes the full layout. To rerun, create a different directory and copy only `inputs/`, `scripts/` and `FROZEN_INPUT_SHA256.tsv` into it; do not overwrite delivered fits. Restore or supply the exact Python and R environments using `environment/python/README.md` and `environment/fgsea/README.md`, then run the parameterized controller as documented.

`FROZEN_INPUT_SHA256.tsv` lists all 58 staged inputs. The original source hashes are in `metadata/ORIGINAL_INPUT_SHA256.tsv`. `metadata/SOURCE_PUBLIC_SHA256.tsv` records copied-file source/public hashes and path-only text sanitization; its paths refer to the full Release tree. `metadata/FULL_RELEASE_TREE_SHA256SUMS.txt` verifies that extracted full tree. This compact copy omits large `inputs/`, complete `fits/` and raw execution logs by design; its own file hashes are in `REPOSITORY_SHA256SUMS.txt`.

The methods were planned before these sensitivity results. All 48 fits were executed and independently checked. A second clean environment restoration was not executed. LOO ranges describe overlapping subset sensitivity, not independent replication, confidence intervals or statistical power. Tied-rank warnings are preserved in `summary/WARNINGS_BY_FIT.tsv`; no PyDESeq2 warnings or unavailable pathway families occurred.
