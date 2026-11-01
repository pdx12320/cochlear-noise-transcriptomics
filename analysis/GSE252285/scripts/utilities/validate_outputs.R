suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

root <- normalizePath(commandArgs(trailingOnly = TRUE)[1], mustWork = TRUE)
p <- function(...) file.path(root, ...)

required <- c(
  "run_all.sh",
  paste0("03_scripts/0", 1:6, c("_input_qc.R", "_seurat_qc.R", "_clustering.R",
                                "_markers_annotation.R", "_top20.R", "_composition.R")),
  "00_manifest/INPUT_AUDIT.tsv", "00_manifest/BARCODE_SUFFIX_COUNTS.tsv",
  "00_manifest/library_map_template.tsv", "04_tables/QC_THRESHOLDS.tsv",
  "04_tables/QC_FLOW.tsv", "04_tables/QC_CELL_METADATA.tsv.gz",
  "04_tables/DOUBLET_RESULTS.tsv.gz", "04_tables/ALL_CLUSTER_MARKERS.tsv.gz",
  "04_tables/TOP20_MARKERS_PER_CLUSTER.tsv", "04_tables/TOP50_MARKERS_PER_CLUSTER.tsv",
  "04_tables/BROAD_ANNOTATION.tsv", "04_tables/TOP20_FEATURE_AUDIT.tsv",
  "04_tables/TOP20_BY_BROAD_CELLTYPE.tsv", "04_tables/CELLTYPE_COMPOSITION_DESCRIPTIVE.tsv",
  paste0("06_objects/0", 1:4, c("_raw_seurat.rds", "_qc_seurat.rds",
                                "_clustered_seurat.rds", "_annotated_seurat.rds")),
  "07_logs/pipeline.log", "07_logs/sessionInfo.txt", "07_logs/software_versions.txt",
  "08_methods/METHODS_SINGLE_CELL.md", "10_release/MODULE_SUMMARY.md", "10_release/SUCCESS"
)

figure_stems <- c(
  "QC_nFeature_VIOLIN", "QC_nCount_VIOLIN", "QC_percent_mt_VIOLIN",
  "QC_nFeature_vs_nCount_SCATTER", "QC_nCount_vs_percent_mt_SCATTER",
  "QC_DENSITY_HISTOGRAM", "QC_LIBRARY_CELL_COUNTS", "DOUBLET_RATE_BY_LIBRARY",
  "PCA_ELBOW", "UMAP_BY_CLUSTER_RES0.6", "UMAP_BY_LIBRARY", "PCA_BY_LIBRARY",
  "CLUSTER_SIZE_BARPLOT", "CLUSTER_BY_LIBRARY_COMPOSITION",
  "TOP_MARKER_HEATMAP", "SELECTED_MARKER_DOTPLOT_BY_CLUSTER",
  "UMAP_BY_BROAD_ANNOTATION", "UMAP_BROAD_ANNOTATION_SPLIT_LIBRARY",
  "CANONICAL_MARKER_DOTPLOT", "CANONICAL_MARKER_HEATMAP",
  "TOP20_BROAD_CELLTYPE_DOTPLOT", "TOP20_MEAN_EXPRESSION_HEATMAP",
  "TOP20_FEATUREPLOT_ATLAS", "TOP20_VLNPLOT",
  "TOP20_Apold1_FeaturePlot_HIGHRES", "TOP20_Apold1_VlnPlot_HIGHRES",
  "TOP20_Apold1_DotPlot_HIGHRES", "TOP20_Errfi1_FeaturePlot_HIGHRES",
  "TOP20_Errfi1_VlnPlot_HIGHRES", "TOP20_Errfi1_DotPlot_HIGHRES",
  "CELLTYPE_COMPOSITION_STACKED_BAR"
)
required <- c(required, as.vector(outer(file.path("05_figures", figure_stems),
                                        c(".pdf", ".png"), paste0)))
paths <- p(required)
missing <- required[!file.exists(paths) | file.info(paths)$size <= 0]
if (length(missing)) stop("Missing/empty required outputs: ", paste(missing, collapse = ", "))

flow <- fread(p("04_tables", "QC_FLOW.tsv"))
libs <- fread(p("04_tables", "CELLS_PER_LIBRARY_FINAL.tsv"))
ann <- fread(p("04_tables", "BROAD_ANNOTATION.tsv"))
top_audit <- fread(p("04_tables", "TOP20_FEATURE_AUDIT.tsv"))
comp <- fread(p("04_tables", "CELLTYPE_COMPOSITION_DESCRIPTIVE.tsv"))
markers <- fread(p("04_tables", "ALL_CLUSTER_MARKERS.tsv.gz"))
doublets <- fread(p("04_tables", "DOUBLET_RESULTS.tsv.gz"))
qc_meta <- fread(p("04_tables", "QC_CELL_METADATA.tsv.gz"))
suffix <- fread(p("00_manifest", "BARCODE_SUFFIX_COUNTS.tsv"))
map <- fread(p("00_manifest", "library_map_template.tsv"))
obj <- readRDS(p("06_objects", "04_annotated_seurat.rds"))

stopifnot(
  identical(readLines(p("10_release", "SUCCESS"), warn = FALSE)[1], "SUCCESS"),
  nrow(qc_meta) == 39092L,
  sum(suffix$cell_number) == 39092L,
  nrow(doublets) == flow[step == "QC_RETAINED", n_cells],
  ncol(obj) == flow[step == "FINAL_CELLS", n_cells],
  nrow(obj) == 32285L,
  length(unique(obj$seurat_clusters)) == 29L,
  all(c("barcode", "barcode_suffix", "library_id", "condition", "final_broad_label",
        "annotation_confidence") %in% colnames(obj[[]])),
  !anyNA(obj$final_broad_label),
  !anyNA(obj$annotation_confidence),
  nrow(ann) == 29L,
  uniqueN(markers$cluster) == 29L,
  nrow(top_audit) == 20L,
  identical(top_audit$gene, c("Dusp1", "Fos", "Nfkbia", "Arrdc2", "Apold1", "Sik1",
                               "Per1", "Ccn1", "Cebpb", "Fosb", "Rasd1", "Gadd45g", "Depp1",
                               "Serpine1", "Errfi1", "Nr4a1", "Sgk1", "Rgs1", "Egr1", "Jun")),
  all(top_audit$status %in% c("FOUND", "ALIAS_REQUIRED", "NOT_FOUND")),
  sum(top_audit$status == "NOT_FOUND") == 0L,
  all(map$condition == "UNKNOWN")
)

feature_md <- obj[["RNA"]][[]]
stopifnot(nrow(feature_md) == nrow(obj), !anyNA(feature_md$feature_id),
          !anyNA(feature_md$original_symbol), !anyDuplicated(rownames(obj)))

comp_check <- comp[, .(n = sum(cell_count)), by = library_id]
setnames(comp_check, "n", "cell_number")
setkey(comp_check, library_id)
setkey(libs, library_id)
stopifnot(identical(comp_check$cell_number, libs$cell_number))

cat("VALIDATION_SUCCESS\n")
cat("required_files=", length(required), "\n", sep = "")
cat("object_dimension=", nrow(obj), "x", ncol(obj), "\n", sep = "")
cat("clusters=", uniqueN(obj$seurat_clusters), "\n", sep = "")
cat("markers=", nrow(markers), "\n", sep = "")
cat("unresolved_clusters=", paste(ann[final_broad_label == "UNRESOLVED", cluster], collapse = ","), "\n", sep = "")
cat("broad_counts\n")
print(as.data.table(table(obj$final_broad_label))[order(-N)])
