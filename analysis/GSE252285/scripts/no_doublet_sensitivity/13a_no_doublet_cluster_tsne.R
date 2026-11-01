source(file.path(Sys.getenv("SCRIPT_DIR", "."), "branch_common.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(patchwork)
})

stage_message("Stage 13a: exact robust-QC subset, no doublet filtering, PCA/Louvain/t-SNE")
raw_path <- file.path(project_dir, "06_objects", "01_raw_seurat.rds")
qc_meta_path <- file.path(project_dir, "04_tables", "QC_CELL_METADATA.tsv.gz")
qc_threshold_path <- file.path(project_dir, "04_tables", "QC_THRESHOLDS.tsv")

qc_cols <- c("barcode", "library_id", "nCount_RNA", "nFeature_RNA", "percent.mt",
             "min_nFeature_RNA", "max_nFeature_RNA", "max_nCount_RNA",
             "max_percent_mt", "pass_min_features", "pass_max_features",
             "pass_max_counts", "pass_mito", "qc_pass")
qc <- fread(qc_meta_path, select = qc_cols, showProgress = FALSE)
if (nrow(qc) != 39092L) stop("Expected 39,092 QC metadata rows; observed ", nrow(qc))
if (any(grepl("scDblFinder", names(qc), ignore.case = TRUE))) {
  stop("Doublet fields must not be imported into this branch")
}
keep <- qc[qc_pass %in% TRUE, barcode]
if (length(keep) != 36111L || uniqueN(keep) != 36111L) {
  stop("Expected 36,111 unique robust-QC barcodes; observed ", length(keep))
}

raw <- readRDS(raw_path)
if (nrow(raw) != 32285L || ncol(raw) != 39092L) {
  stop("Raw object invariant failed: ", nrow(raw), " genes x ", ncol(raw), " cells")
}
if (!all(keep %in% colnames(raw))) stop("Some robust-QC barcodes are absent from raw object")

raw$percent.mt <- qc$percent.mt[match(colnames(raw), qc$barcode)]
raw$qc_pass <- qc$qc_pass[match(colnames(raw), qc$barcode)]
obj <- subset(raw, cells = keep)
rm(raw)
gc(verbose = FALSE)

if (ncol(obj) != 36111L) stop("Robust-QC subset did not conserve 36,111 cells")
if (any(grepl("scDblFinder", colnames(obj[[]]), ignore.case = TRUE))) {
  stop("New object unexpectedly contains scDblFinder metadata")
}
obj$doublet_filter_applied <- "NO"
obj$analysis_branch <- "NO_DOUBLETFILTER_SENSITIVITY"
obj$embedding_for_display <- "tSNE_only_no_UMAP"
DefaultAssay(obj) <- "RNA"
set.seed(seed)

stage_message("NormalizeData / 2,000 HVGs / ScaleData / PCA")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,
                     verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000,
                            verbose = FALSE)
if (length(VariableFeatures(obj)) != 2000L) {
  stop("Expected exactly 2,000 variable features; observed ", length(VariableFeatures(obj)))
}
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, seed.use = seed,
              verbose = FALSE)

pc_stdev <- Stdev(obj[["pca"]])
pc_var <- pc_stdev^2
pca_variance <- data.table(
  PC = seq_along(pc_var),
  variance = pc_var,
  variance_fraction = pc_var / sum(pc_var),
  cumulative_variance_fraction = cumsum(pc_var / sum(pc_var)),
  selected_for_main = seq_along(pc_var) <= 30
)
write_tsv(pca_variance, branch_path("04_tables", "PCA_VARIANCE.tsv"))
p_elbow <- ElbowPlot(obj, ndims = 50) +
  geom_vline(xintercept = 30, linetype = 2, colour = "#D55E00") +
  labs(title = "PCA elbow: no-doublet-filter sensitivity",
       subtitle = "Thirty PCs were pre-specified; not selected to match a cluster count") +
  theme_branch()
save_branch_dual(p_elbow, "PCA_ELBOW", 8, 6)

main_dims <- 1:30
stage_message("FindNeighbors / fixed Louvain resolution audit")
set.seed(seed)
obj <- FindNeighbors(obj, dims = main_dims, k.param = 20, verbose = FALSE)
set.seed(seed)
obj <- FindClusters(obj, resolution = c(0.4, 0.6, 0.8, 1.0), algorithm = 1,
                    random.seed = seed, n.start = 10, n.iter = 10,
                    verbose = FALSE)
main_cluster_col <- "RNA_snn_res.0.6"
if (!main_cluster_col %in% colnames(obj[[]])) stop("Missing ", main_cluster_col)
Idents(obj) <- main_cluster_col
obj$seurat_clusters_no_doublet <- as.character(Idents(obj))

stage_message("RunTSNE from PCs 1:30; perplexity=30; UMAP is not run")
set.seed(seed)
obj <- RunTSNE(obj, reduction = "pca", dims = main_dims, seed.use = seed,
               perplexity = 30, check_duplicates = FALSE, verbose = FALSE)
if ("umap" %in% Reductions(obj)) stop("UMAP must not exist in the new branch object")

res_cols <- paste0("RNA_snn_res.", c("0.4", "0.6", "0.8", "1"))
res_cols <- res_cols[res_cols %in% colnames(obj[[]])]
resolution_audit <- rbindlist(lapply(res_cols, function(rc) {
  data.table(
    resolution = sub("^RNA_snn_res\\.", "", rc),
    number_of_clusters = uniqueN(obj[[rc, drop = TRUE]]),
    cell_number = ncol(obj),
    selected_for_main = rc == main_cluster_col
  )
}))
write_tsv(resolution_audit, branch_path("04_tables", "CLUSTER_RESOLUTION_AUDIT.tsv"))

qc_flow <- data.table(
  step = c("INITIAL_CELLS", "ROBUST_QC_RETAINED", "DOUBLETS_REMOVED", "FINAL_CELLS"),
  n_cells = c(39092L, 36111L, 0L, 36111L),
  definition = c("official processed matrix barcodes",
                 "exact existing per-library robust-QC pass set",
                 "not performed in this sensitivity branch",
                 "all robust-QC cells; no doublet filter")
)
write_tsv(qc_flow, branch_path("04_tables", "CELL_FLOW.tsv"))
file.copy(qc_threshold_path, branch_path("09_qc", "QC_THRESHOLDS_USED.tsv"), overwrite = TRUE)

lib_counts <- as.data.table(obj[[]])[, .(cell_number = .N), by = library_id][order(library_id)]
write_tsv(lib_counts, branch_path("04_tables", "CELLS_PER_LIBRARY.tsv"))
expected_lib <- data.table(library_id = c("LIB1", "LIB2", "LIB3"),
                           expected = c(16631L, 9207L, 10273L))
lib_check <- merge(lib_counts, expected_lib, by = "library_id", all = TRUE)
if (any(lib_check$cell_number != lib_check$expected)) stop("Library QC cell-count invariant failed")

cluster_sizes <- as.data.table(table(obj$seurat_clusters_no_doublet))
setnames(cluster_sizes, c("cluster", "cell_number"))
cluster_sizes[, cluster := as.character(cluster)]
write_tsv(cluster_sizes, branch_path("04_tables", "CLUSTER_SIZES_RES0.6.tsv"))

cluster_lib <- as.data.table(obj[[]])[, .(cell_number = .N),
  by = .(cluster = seurat_clusters_no_doublet, library_id)]
cluster_lib[, fraction_within_cluster := cell_number / sum(cell_number), by = cluster]
cluster_lib[, fraction_within_library := cell_number / sum(cell_number), by = library_id]
write_tsv(cluster_lib, branch_path("04_tables", "CLUSTER_BY_LIBRARY_COMPOSITION.tsv"))

tab <- table(obj$seurat_clusters_no_doublet, obj$library_id)
chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
cramers_v <- sqrt(as.numeric(chi$statistic) /
                    (sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1)))
confound_status <- if (is.finite(cramers_v) && cramers_v >= 0.5) {
  "POTENTIAL_LIBRARY_CONDITION_CONFOUNDING"
} else {
  "NO_EXTREME_CLUSTER_LIBRARY_ASSOCIATION_BY_PRESET_CRAMERS_V_0.5"
}
write_tsv(data.table(
  metric = c("cluster_library_cramers_v", "preset_warning_threshold", "status",
             "condition_mapping_status"),
  value = c(sprintf("%.6f", cramers_v), "0.5", confound_status,
            ifelse(any(obj$condition == "UNKNOWN"), "CONDITION_MAPPING_UNRESOLVED", "MAPPED_BY_USER"))
), branch_path("09_qc", "LIBRARY_CONFOUNDING_AUDIT.tsv"))

tsne_df <- as.data.table(Embeddings(obj, "tsne"), keep.rownames = "barcode")
tsne_df[, `:=`(
  cluster = obj$seurat_clusters_no_doublet[match(barcode, colnames(obj))],
  library_id = obj$library_id[match(barcode, colnames(obj))],
  nCount_RNA = obj$nCount_RNA[match(barcode, colnames(obj))],
  nFeature_RNA = obj$nFeature_RNA[match(barcode, colnames(obj))],
  percent.mt = obj$percent.mt[match(barcode, colnames(obj))]
)]
write_tsv(tsne_df, branch_path("04_tables", "TSNE_EMBEDDINGS.tsv.gz"))

p_cluster <- DimPlot(obj, reduction = "tsne", group.by = "seurat_clusters_no_doublet",
                     label = TRUE, repel = TRUE, raster = TRUE) +
  labs(title = "Robust-QC atlas without doublet filtering",
       subtitle = "Louvain resolution 0.6; t-SNE is for visualization only") +
  theme_branch()
p_library <- DimPlot(obj, reduction = "tsne", group.by = "library_id",
                     raster = TRUE, shuffle = TRUE, seed = seed) +
  scale_colour_manual(values = c(LIB1 = "#0072B2", LIB2 = "#D55E00", LIB3 = "#009E73")) +
  labs(title = "t-SNE by library",
       subtitle = "Condition mapping unresolved; no integration was applied") +
  theme_branch()
save_branch_dual(p_cluster, "TSNE_BY_CLUSTER_RES0.6", 10, 8)
save_branch_dual(p_library, "TSNE_BY_LIBRARY", 9, 7)

p_count <- FeaturePlot(obj, reduction = "tsne", features = "nCount_RNA",
                       raster = TRUE, order = TRUE) + theme_branch()
p_feature <- FeaturePlot(obj, reduction = "tsne", features = "nFeature_RNA",
                         raster = TRUE, order = TRUE) + theme_branch()
p_mito <- FeaturePlot(obj, reduction = "tsne", features = "percent.mt",
                      raster = TRUE, order = TRUE) + theme_branch()
save_branch_dual(p_count + p_feature + p_mito,
                 "TSNE_TECHNICAL_COVARIATE_AUDIT", 16, 5.5)

p_size <- ggplot(cluster_sizes, aes(reorder(cluster, -cell_number), cell_number)) +
  geom_col(fill = "#0072B2") +
  geom_text(aes(label = cell_number), vjust = -0.2, size = 2.7) +
  labs(x = "Cluster (resolution 0.6)", y = "Cells",
       title = "Cluster sizes: no-doublet-filter sensitivity") +
  theme_branch()
save_branch_dual(p_size, "CLUSTER_SIZE_BARPLOT", 12, 6)

p_comp <- ggplot(cluster_lib, aes(cluster, fraction_within_cluster, fill = library_id)) +
  geom_col() +
  scale_fill_manual(values = c(LIB1 = "#0072B2", LIB2 = "#D55E00", LIB3 = "#009E73")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Library fraction within cluster",
       title = "Cluster composition by library") +
  theme_branch()
save_branch_dual(p_comp, "CLUSTER_BY_LIBRARY_COMPOSITION", 12, 6)

write_tsv(data.table(
  parameter = c("analysis_status", "seed", "input_cells", "doublet_filter",
                "normalization", "variable_features", "PCA_computed", "PCs_used",
                "neighbors_k", "clustering_algorithm", "main_resolution",
                "resolution_audit", "display_embedding", "tSNE_perplexity",
                "integration", "condition_mapping"),
  value = c("NO_DOUBLETFILTER_SENSITIVITY_NOT_FROZEN", "42", "36111",
            "NONE", "LogNormalize_scale.factor_10000", "vst_2000", "50", "1:30",
            "20", "Louvain_algorithm_1", "0.6", "0.4,0.6,0.8,1.0",
            "tSNE_only_no_UMAP", "30", "NONE", "CONDITION_MAPPING_UNRESOLVED")
), branch_path("00_manifest", "ANALYSIS_PARAMETERS.tsv"))

object_path <- branch_path("06_objects", "01_no_doublet_clustered_tsne.rds")
saveRDS(obj, object_path, compress = FALSE)
stage_message("Stage 13a complete: cells=", ncol(obj),
              "; clusters(res0.6)=", uniqueN(obj$seurat_clusters_no_doublet),
              "; Cramer's V=", sprintf("%.3f", cramers_v),
              "; object=", object_path)
