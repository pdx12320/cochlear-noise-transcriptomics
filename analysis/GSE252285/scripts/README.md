# Scripts for the selected 36,111-cell result

## Current modules

| Directory | Purpose | Main entry point |
|---|---|---|
| `stage1/` | Public input audit and upstream raw/QC provenance | `01_input_qc.R`, then `02_seurat_qc.R`; see `../inputs/README.md` |
| `no_doublet_sensitivity/` | Exact 36,111-cell robust-QC subset, PCA, Louvain, t-SNE, marker annotation, Top20 and resolution audit | `run_no_doublet_robust_qc_tsne.sh` |
| `four_cluster_audit/` | Marker/reference audit of clusters 6/9/13/19 and targeted erythroid/co-expression checks | numbered R scripts |
| `condition_descriptive/` | Verified suffix mapping, condition composition and Top20 descriptive summaries | `17_condition_descriptive_top20.R`, then `18_finalize_condition_release.R` |
| `utilities/` | Output validation, format conversion and Git repository manifests | individual scripts |

Machine-specific paths were replaced by `GSE252285_PROJECT_ROOT`, `PROJECT_DIR`, `SCRIPT_DIR` or `PAPER_PDF`. Numerical thresholds and executed scientific parameters were not changed.

## Reproduce the current branch from its audited inputs

```bash
export GSE252285_PROJECT_ROOT=/path/to/new/GSE252285_run
export PROJECT_DIR="$GSE252285_PROJECT_ROOT"
export SCRIPT_DIR="$PWD/scripts/no_doublet_sensitivity"
export R_BIN=/path/to/Rscript
export SC_THREADS=4

bash scripts/no_doublet_sensitivity/run_no_doublet_robust_qc_tsne.sh
```

Required external inputs are:

```text
06_objects/01_raw_seurat.rds
04_tables/QC_CELL_METADATA.tsv.gz
04_tables/QC_THRESHOLDS.tsv
```

The wrapper also reads the supplied `inputs/annotation_reference.tsv`. Override its location with `ANNOTATION_REFERENCE_TSV` if needed.

The script checks 39,092 QC metadata rows, exactly 36,111 passing barcodes, 32,285 features, absence of scDblFinder fields and the expected library totals before analysis. The supplied barcode-level reference labels are required annotation priors. Canonical marker support and reference-label purity jointly determine the unchanged adopted labels. The full earlier object is not required; no earlier cell-inclusion or doublet-class field is imported. Optional legacy comparison outputs are disabled by default (`RUN_LEGACY_COMPARISON=false`).

## Resolved-condition description

```bash
export SCRIPT_DIR="$PWD/scripts/no_doublet_sensitivity"
"$R_BIN" scripts/condition_descriptive/17_condition_descriptive_top20.R
"$R_BIN" scripts/condition_descriptive/18_finalize_condition_release.R
```

The condition script verifies `GEO_-3/LIB3=Control`, `GEO_-1/LIB1=HHL` and `GEO_-2/LIB2=NIHL` from cell-name suffixes. It never guesses from library numbering and performs no inferential condition test.

## Upstream reconstruction from Matrix Market files

`stage1/` retains the input-audit and upstream QC code required to create the raw object and QC metadata from:

```text
GSE252285_matrix.mtx.gz
GSE252285_features.tsv.gz
GSE252285_barcodes.tsv.gz
```

The first two upstream stages produce the raw object, library-map template and QC records needed by the selected workflow. The second stage also computes a doublet audit; those classes are ignored by the 36,111-cell branch. The full stage1 wrapper additionally runs historical downstream stages that are unnecessary for the current branch. The complete current-input inventory and lightweight adapter checks are documented in `../inputs/README.md`.

## Historical layout note

Some executed module scripts expect active files under `${PROJECT_DIR}/03_scripts`. For a rerun, use a new work directory and copy the selected scripts there or adjust `SCRIPT_DIR`. Do not overwrite a frozen analysis project.
