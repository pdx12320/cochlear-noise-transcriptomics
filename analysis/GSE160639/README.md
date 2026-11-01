# GSE160639: current raw-count analysis

Executed 1 November 2026. The workflow starts from the public whole-cochlea count matrix, retains the fixed 12 samples (four per exposure group) and 19,673-gene universe, estimates median-ratio normalization, and fits DESeq2 models. Three Wald-statistic contrasts feed nine targeted pathways, with BH adjustment across all 27 pathway comparisons. The manuscript and Figures 3 and S1 use these results.

## Files

- `deseq2/inputs/`: deposited 35,496-gene count matrix, 12-sample manifest and fixed gene identifiers/annotations.
- `deseq2/scripts/`: executed DESeq2 script and publication-object metadata cleaner.
- `deseq2/results/`: four full gene-result tables, normalized counts, size factors, diagnostics and parameters.
- `deseq2/environment/`: exact Linux runtime lock, installed R package versions, source hashes and restoration commands.
- `deseq2/docs/`: execution logs, independent checks and model-object equivalence records.
- `enrichment/`: executed enrichment script, frozen MSigDB 2026.1.Mm inputs, full ranks, memberships, leading edges and 27 results.

## Run

Use a separate output directory for each execution. The DESeq2 stage requires R 4.6.1 / DESeq2 1.52.0; the enrichment stage requires R 4.5.3 / fgsea 1.36.2. Select the matching runtime before each stage. Restore instructions for the DESeq2 environment are in [`deseq2/environment/README.md`](deseq2/environment/README.md); enrichment versions are recorded with its results and environment files.

From the repository root, in the DESeq2 runtime:

```sh
Rscript --vanilla analysis/GSE160639/deseq2/scripts/authorized_bulk_deseq2_reestimate.R \
  --preflight-only analysis/GSE160639/deseq2/inputs/GSE160639_RawReadCount.txt.gz \
  analysis/GSE160639/deseq2/inputs reproduced/GSE160639_DE
Rscript --vanilla analysis/GSE160639/deseq2/scripts/authorized_bulk_deseq2_reestimate.R \
  --execute-approved analysis/GSE160639/deseq2/inputs/GSE160639_RawReadCount.txt.gz \
  analysis/GSE160639/deseq2/inputs reproduced/GSE160639_DE
```

In the enrichment runtime, the following command uses the adopted current DE tables. The script verifies all input hashes and the fixed pathway membership before computing:

```sh
Rscript --vanilla analysis/GSE160639/enrichment/authorized_targeted_bulk_from_rawcounts.R \
  --preflight-only analysis/GSE160639/deseq2/results \
  analysis/GSE160639/deseq2/inputs/03b_GSE160639_sample_manifest_locked.csv \
  analysis/GSE160639/enrichment/frozen_gene_sets reproduced/GSE160639_GSEA
Rscript --vanilla analysis/GSE160639/enrichment/authorized_targeted_bulk_from_rawcounts.R \
  --execute-approved analysis/GSE160639/deseq2/results \
  analysis/GSE160639/deseq2/inputs/03b_GSE160639_sample_manifest_locked.csv \
  analysis/GSE160639/enrichment/frozen_gene_sets reproduced/GSE160639_GSEA
```

The new DE output can replace the archived input directory when its three result-table hashes match the recorded execution. Both scripts require a nonexistent output directory; preflight mode creates no result directory. Use one OpenMP/BLAS thread. Recorded execution limits were one statistical process, 8 GiB and 1,200 seconds for DESeq2, and one process, 4 GiB and 900 seconds for enrichment.

## Results and checks

The current gene-level BH-positive counts are 5,777 (94 vs 70 dB), 889 (105 vs 70 dB), 4,177 (105 vs 94 dB), and 6,398 for the overall-dose LRT. Targeted enrichment uses 17,441 ranked symbols per contrast. Its global-FDR-positive counts are 8/9, 9/9 and 0/9, respectively.

Independent checks reproduced size factors, normalized counts, Wald/LRT P values and gene-level BH adjustment, and verified all 27 pathway identities, ranking choices, overlaps, leading-edge membership and enrichment BH adjustment. Detailed records are in `deseq2/docs/` and `enrichment/INDEPENDENT_STATISTICAL_AUDIT.json`.

## Release assets

- [Fixed R package source archives](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/download/v1.0.0/Bulk_R_Source_Packages_20261101.tar)
- [DESeq2 normalized and fitted objects](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/download/v1.0.0/GSE160639_DESeq2_Objects_20261101.tar.gz)

The model bundle contains publication copies with machine-local formula/function environment references removed. All assays, row/column metadata, dispersion predictions and re-extracted Wald/LRT results were verified identical. The bundle records each original and public SHA256 and the exact transformation; original analysis objects remain preserved separately.
