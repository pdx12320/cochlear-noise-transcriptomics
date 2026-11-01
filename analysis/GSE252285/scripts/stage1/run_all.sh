#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C
export LANG=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${GSE252285_PROJECT_ROOT:-${PROJECT_DIR:-$(pwd)}}"
ROOT_DIR="$(cd "$ROOT_DIR" && pwd)"
cd "$ROOT_DIR"

for d in 00_manifest 01_input 02_processed 03_scripts 04_tables 05_figures 06_objects 07_logs 08_methods 09_qc 10_release; do
  mkdir -p "$d"
done

LOG_FILE="$ROOT_DIR/07_logs/pipeline.log"
: > "$LOG_FILE"
exec > >(tee -a "$LOG_FILE") 2>&1

export PROJECT_DIR="$ROOT_DIR"
export GSE252285_SCRIPT_DIR="$SCRIPT_DIR"
export GSE_MATRIX="${GSE_MATRIX:-$ROOT_DIR/01_input/GSE252285_matrix.mtx.gz}"
export GSE_FEATURES="${GSE_FEATURES:-$ROOT_DIR/01_input/GSE252285_features.tsv.gz}"
export GSE_BARCODES="${GSE_BARCODES:-$ROOT_DIR/01_input/GSE252285_barcodes.tsv.gz}"
export R_LIBS_USER="${R_LIBS_USER:-$ROOT_DIR/.Rlib}"
mkdir -p "$R_LIBS_USER"

THREADS="${SC_THREADS:-4}"
export SC_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

if [[ -x "$ROOT_DIR/work/r_env/bin/Rscript" ]]; then
  R_BIN="$ROOT_DIR/work/r_env/bin/Rscript"
elif command -v Rscript >/dev/null 2>&1; then
  R_BIN="$(command -v Rscript)"
else
  echo "ERROR: Rscript not found. Expected project-local work/r_env/bin/Rscript."
  exit 2
fi

for input in "$GSE_MATRIX" "$GSE_FEATURES" "$GSE_BARCODES"; do
  if [[ ! -r "$input" ]]; then
    echo "STOP_INPUT_AUDIT: missing or unreadable input: $input"
    printf '%s\n' STOP_INPUT_AUDIT > "$ROOT_DIR/00_manifest/STOP_INPUT_AUDIT"
    exit 3
  fi
done

echo "GSE252285 pipeline start: $(date -Is)"
echo "Project: $ROOT_DIR"
echo "Rscript: $R_BIN"
echo "Threads: $THREADS"
sha256sum "$GSE_MATRIX" "$GSE_FEATURES" "$GSE_BARCODES" > "$ROOT_DIR/00_manifest/INPUT_SHA256.txt"

for script in \
  01_input_qc.R \
  02_seurat_qc.R \
  03_clustering.R \
  04_markers_annotation.R \
  05_top20.R \
  06_composition.R; do
  echo "[$(date -Is)] START $script"
  "$R_BIN" "$SCRIPT_DIR/$script"
  echo "[$(date -Is)] END $script"
done

"$R_BIN" -e 'sink(file.path(Sys.getenv("PROJECT_DIR"),"07_logs","sessionInfo.txt")); sessionInfo(); sink()'
"$R_BIN" -e 'pkgs <- c("R","Seurat","SeuratObject","Matrix","data.table","dplyr","ggplot2","patchwork","scDblFinder","SingleCellExperiment","pheatmap","BiocParallel"); versions <- vapply(pkgs, function(x) if (x=="R") R.version.string else if (requireNamespace(x,quietly=TRUE)) as.character(packageVersion(x)) else "NOT_INSTALLED", character(1)); write.table(data.frame(software=names(versions),version=unname(versions)), file=file.path(Sys.getenv("PROJECT_DIR"),"07_logs","software_versions.txt"), sep="\t", quote=FALSE, row.names=FALSE)'

required_outputs=(
  00_manifest/INPUT_AUDIT.tsv
  00_manifest/BARCODE_SUFFIX_COUNTS.tsv
  00_manifest/library_map_template.tsv
  04_tables/QC_THRESHOLDS.tsv
  04_tables/QC_FLOW.tsv
  04_tables/QC_CELL_METADATA.tsv.gz
  04_tables/DOUBLET_RESULTS.tsv.gz
  04_tables/ALL_CLUSTER_MARKERS.tsv.gz
  04_tables/TOP20_MARKERS_PER_CLUSTER.tsv
  04_tables/TOP50_MARKERS_PER_CLUSTER.tsv
  04_tables/BROAD_ANNOTATION.tsv
  04_tables/TOP20_FEATURE_AUDIT.tsv
  04_tables/TOP20_BY_BROAD_CELLTYPE.tsv
  04_tables/CELLTYPE_COMPOSITION_DESCRIPTIVE.tsv
  05_figures/DOUBLET_RATE_BY_LIBRARY.pdf
  05_figures/DOUBLET_RATE_BY_LIBRARY.png
  05_figures/UMAP_BY_CLUSTER_RES0.6.pdf
  05_figures/UMAP_BY_CLUSTER_RES0.6.png
  05_figures/TOP_MARKER_HEATMAP.pdf
  05_figures/TOP_MARKER_HEATMAP.png
  05_figures/UMAP_BY_BROAD_ANNOTATION.pdf
  05_figures/UMAP_BY_BROAD_ANNOTATION.png
  05_figures/TOP20_BROAD_CELLTYPE_DOTPLOT.pdf
  05_figures/TOP20_BROAD_CELLTYPE_DOTPLOT.png
  05_figures/TOP20_FEATUREPLOT_ATLAS.pdf
  05_figures/TOP20_FEATUREPLOT_ATLAS.png
  05_figures/TOP20_VLNPLOT.pdf
  05_figures/TOP20_VLNPLOT.png
  05_figures/TOP20_Apold1_FeaturePlot_HIGHRES.pdf
  05_figures/TOP20_Errfi1_FeaturePlot_HIGHRES.pdf
  05_figures/CELLTYPE_COMPOSITION_STACKED_BAR.pdf
  05_figures/CELLTYPE_COMPOSITION_STACKED_BAR.png
  06_objects/01_raw_seurat.rds
  06_objects/02_qc_seurat.rds
  06_objects/03_clustered_seurat.rds
  06_objects/04_annotated_seurat.rds
  07_logs/sessionInfo.txt
  07_logs/software_versions.txt
  08_methods/METHODS_SINGLE_CELL.md
  10_release/MODULE_SUMMARY.md
)
for output in "${required_outputs[@]}"; do
  if [[ ! -s "$ROOT_DIR/$output" ]]; then
    echo "ERROR: required output missing or empty: $output"
    exit 4
  fi
done

printf '%s\n' SUCCESS > "$ROOT_DIR/10_release/SUCCESS"
echo "GSE252285 pipeline SUCCESS: $(date -Is)"
