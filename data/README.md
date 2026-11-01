# Data access

## Public source data

| Dataset | GEO record | Role |
|---|---|---|
| GSE168041 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168041 | Replicated single-cell counts, two cochlear compartments |
| GSE160639 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160639 | Whole-cochlea count matrix and exposure contrasts |
| GSE252285 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252285 | Three pooled-library expression-localization profiles |

Public FASTQ files and full software runtimes are not duplicated in Git. The source accessions, sample records and archived input hashes identify the inputs used by each workflow.

## Processed analysis objects

The data-only [v1.0.0 collection](https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/tag/v1.0.0) contains the selected 36,111-cell GSE252285 annotated Seurat object, both GSE168041 annotated compartment objects, and GSE168041 pseudobulk/count/model records. The object manifest records byte counts and SHA256 values. The Seurat object is split into numbered binary parts because a single GitHub release asset must be smaller than 2 GiB; reassembly reproduces the original file byte for byte.

See `manifests/RELEASE_OBJECTS.tsv` and the retrieval script under `scripts/`. Downloaded objects belong under `data/objects/` and are ignored by Git.

## Bulk model objects and runtime sources

The release additionally provides `GSE160639_DESeq2_Objects_20261101.tar.gz`, containing the current normalized object and Wald/LRT fits, plus `Bulk_R_Source_Packages_20261101.tar`, containing the exact R package source distributions used for the DESeq2 environment. Their source/public object hashes, equivalence checks and restoration instructions are under [`analysis/GSE160639/deseq2/`](../analysis/GSE160639/deseq2/). Machine-local environment references were removed from the public model copies; all numerical data and re-extracted results are identical.

The complete release-asset inventory is recorded in `manifests/RELEASE_ASSETS.tsv`. The public bulk count matrix and fixed metadata are small enough to be included directly in Git under `analysis/GSE160639/deseq2/inputs/`.

## Manuscript data

`manuscript/source_inputs/` contains the inputs used to generate the current figures. `manuscript/source_data/` contains the exact plotted-panel tables. Current bulk estimates and their full executed workflow are separately identified under `analysis/GSE160639/`. Superseded pathway tables and internal version comparisons are excluded.

Three GSE252285 libraries contribute descriptive condition profiles. Cells are not used as independent biological replicates for condition inference. The sample-level GSE168041 and GSE160639 designs are documented in the manuscript.
