source(file.path(Sys.getenv("SCRIPT_DIR", "."), "branch_common.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(png)
})

stage_message("Validation: no-doublet robust-QC t-SNE branch")
object_path <- branch_path("06_objects", "02_no_doublet_annotated_tsne.rds")
obj <- readRDS(object_path)
checks <- list(
  cell_count_36111 = ncol(obj) == 36111L,
  gene_count_32285 = nrow(obj) == 32285L,
  no_scDblFinder_metadata = !any(grepl("scDblFinder", colnames(obj[[]]), ignore.case = TRUE)),
  doublet_filter_flag_NO = all(obj$doublet_filter_applied == "NO"),
  tsne_present = "tsne" %in% Reductions(obj),
  umap_absent = !"umap" %in% Reductions(obj),
  main_cluster_field_present = "RNA_snn_res.0.6" %in% colnames(obj[[]]),
  annotation_complete = !anyNA(obj$broad_annotation_no_doublet_sensitivity),
  annotation_count_conserved = length(obj$broad_annotation_no_doublet_sensitivity) == 36111L,
  exact_HVG_count = length(VariableFeatures(obj)) == 2000L
)

top20 <- fread(branch_path("04_tables", "TOP20_FEATURE_AUDIT.tsv"))
checks$top20_20_rows <- nrow(top20) == 20L
checks$top20_all_found <- all(top20$status %in% c("FOUND", "ALIAS_REQUIRED"))
flow <- fread(branch_path("04_tables", "CELL_FLOW.tsv"))
checks$flow_final_36111 <- flow[step == "FINAL_CELLS", n_cells] == 36111L
checks$flow_doublets_removed_zero <- flow[step == "DOUBLETS_REMOVED", n_cells] == 0L
composition <- fread(branch_path("04_tables", "CELLTYPE_COMPOSITION_DESCRIPTIVE.tsv"))
checks$composition_conserved <- sum(composition$cell_count) == 36111L
annotation <- fread(branch_path("04_tables", "BROAD_ANNOTATION_SENSITIVITY.tsv"))
checks$cluster_annotation_one_to_one <- nrow(annotation) == uniqueN(obj$seurat_clusters_no_doublet)
res <- fread(branch_path("04_tables", "CLUSTER_RESOLUTION_AUDIT.tsv"))
checks$one_main_resolution <- sum(res$selected_for_main) == 1L

required_stems <- c(
  "PCA_ELBOW", "TSNE_BY_CLUSTER_RES0.6", "TSNE_BY_LIBRARY",
  "TSNE_TECHNICAL_COVARIATE_AUDIT", "CLUSTER_SIZE_BARPLOT",
  "CLUSTER_BY_LIBRARY_COMPOSITION", "TSNE_BY_BROAD_ANNOTATION_SENSITIVITY",
  "CANONICAL_MARKER_DOTPLOT_BY_CLUSTER", "LINEAGE_SCORE_HEATMAP",
  "TOP20_DOTPLOT",
  "TOP20_MEAN_EXPRESSION_HEATMAP", "TOP20_FEATUREPLOT_ATLAS", "TOP20_VLNPLOT",
  "APOLD1_FEATUREPLOT", "APOLD1_VLNPLOT", "APOLD1_DOTPLOT",
  "ERRFI1_FEATUREPLOT", "ERRFI1_VLNPLOT", "ERRFI1_DOTPLOT",
  "CELLTYPE_COMPOSITION_STACKED_BAR"
)
if (run_legacy_comparison) {
  required_stems <- c(required_stems, "ADDITIONAL_CELL_FRACTION_BY_CLUSTER")
}
figure_paths <- unlist(lapply(required_stems, function(stem) {
  c(branch_path("05_figures", paste0(stem, ".pdf")),
    branch_path("05_figures", paste0(stem, ".png")))
}))
checks$all_expected_figures_nonempty <- all(file.exists(figure_paths) & file.info(figure_paths)$size > 1000)

png_paths <- figure_paths[grepl("\\.png$", figure_paths)]
png_dims <- rbindlist(lapply(png_paths, function(path) {
  img <- png::readPNG(path, native = TRUE, info = TRUE)
  data.table(file = basename(path), width_px = dim(img)[2], height_px = dim(img)[1],
             nonempty = length(unique(as.vector(img))) > 1)
}))
write_tsv(png_dims, branch_path("09_qc", "FIGURE_PNG_TECHNICAL_QA.tsv"))
checks$all_pngs_nonempty_pixels <- all(png_dims$nonempty)

check_dt <- data.table(check = names(checks), passed = unlist(checks, use.names = FALSE))
write_tsv(check_dt, branch_path("09_qc", "VALIDATION_CHECKS.tsv"))
if (!all(check_dt$passed)) {
  print(check_dt[passed == FALSE])
  stop("Validation failed")
}

summary_path <- branch_path("10_release", "MODULE_SUMMARY.md")
summary_lines <- readLines(summary_path)
summary_lines <- sub("STATUS: SUCCESS_COMPUTATION_PENDING_FINAL_VALIDATION",
                     "STATUS: SUCCESS", summary_lines, fixed = TRUE)
writeLines(summary_lines, summary_path)

stage_message("Validation SUCCESS: ", nrow(check_dt), " checks passed; ",
              nrow(png_dims), " PNGs passed technical pixel/dimension QA")
