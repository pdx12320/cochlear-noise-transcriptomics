source(file.path(Sys.getenv("GSE252285_SCRIPT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "03_scripts")), "common.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
})

stage_message("Stage 03: unintegrated normalization, PCA, UMAP, and clustering")
obj <- readRDS(p("06_objects", "02_qc_seurat.rds"))
DefaultAssay(obj) <- "RNA"
set.seed(seed)

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,
                     verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000,
                            verbose = FALSE)
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
write_tsv(pca_variance, p("04_tables", "PCA_VARIANCE.tsv"))

p_elbow <- ElbowPlot(obj, ndims = 50) +
  geom_vline(xintercept = 30, linetype = 2, colour = "red") +
  labs(subtitle = "Main analysis pre-specified at 30 PCs") + theme_sc()
save_dual(p_elbow, "PCA_ELBOW", 8, 6)

main_dims <- 1:30
set.seed(seed)
obj <- FindNeighbors(obj, dims = main_dims, k.param = 20, verbose = FALSE)
set.seed(seed)
obj <- RunUMAP(obj, dims = main_dims, seed.use = seed, n.neighbors = 30,
               min.dist = 0.3, verbose = FALSE)
set.seed(seed)
obj <- FindClusters(obj, resolution = c(0.4, 0.6, 0.8, 1.0), algorithm = 1,
                    random.seed = seed, n.start = 10, n.iter = 10,
                    verbose = FALSE)

main_cluster_col <- "RNA_snn_res.0.6"
if (!main_cluster_col %in% colnames(obj[[]])) {
  stop("Expected clustering column not found: ", main_cluster_col)
}
Idents(obj) <- main_cluster_col
obj$seurat_clusters <- as.character(Idents(obj))

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
write_tsv(resolution_audit, p("04_tables", "CLUSTER_RESOLUTION_AUDIT.tsv"))
write_tsv(data.table(
  parameter = c("seed", "normalization", "variable_features", "PCA_computed",
                "PCs_used", "neighbors_k", "UMAP_n_neighbors", "UMAP_min_dist",
                "clustering_algorithm", "main_resolution", "integration"),
  value = c("42", "LogNormalize_scale.factor_10000", "vst_2000", "50", "1:30",
            "20", "30", "0.3", "Louvain_algorithm_1", "0.6", "NONE_unintegrated_atlas")
), p("00_manifest", "CLUSTERING_PARAMETERS.tsv"))

cluster_lib <- as.data.table(obj[[]])[,
  .(cell_number = .N), by = .(cluster = seurat_clusters, library_id)]
cluster_lib[, fraction_within_cluster := cell_number / sum(cell_number), by = cluster]
cluster_lib[, fraction_within_library := cell_number / sum(cell_number), by = library_id]
write_tsv(cluster_lib, p("04_tables", "CLUSTER_BY_LIBRARY_COMPOSITION.tsv"))

tab <- table(obj$seurat_clusters, obj$library_id)
chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
cramers_v <- sqrt(as.numeric(chi$statistic) /
                    (sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1)))
confounding_status <- if (is.finite(cramers_v) && cramers_v >= 0.5) {
  "POTENTIAL_LIBRARY_CONDITION_CONFOUNDING"
} else {
  "NO_EXTREME_CLUSTER_LIBRARY_ASSOCIATION_BY_PRESET_CRAMERS_V_0.5"
}
write_tsv(data.table(
  metric = c("cluster_library_cramers_v", "preset_warning_threshold", "status",
             "condition_mapping_status"),
  value = c(sprintf("%.6f", cramers_v), "0.5", confounding_status,
            ifelse(any(obj$condition == "UNKNOWN"), "CONDITION_MAPPING_UNRESOLVED", "MAPPED_BY_USER"))
), p("09_qc", "LIBRARY_CONFOUNDING_AUDIT.tsv"))

p_umap_cluster <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                          label = TRUE, repel = TRUE, raster = TRUE) +
  labs(title = "Unintegrated atlas: resolution 0.6") + theme_sc()
p_umap_lib <- DimPlot(obj, reduction = "umap", group.by = "library_id",
                      raster = TRUE) + theme_sc()
p_pca_lib <- DimPlot(obj, reduction = "pca", dims = c(1, 2), group.by = "library_id",
                     raster = TRUE) + theme_sc()
save_dual(p_umap_cluster, "UMAP_BY_CLUSTER_RES0.6", 9, 7)
save_dual(p_umap_lib, "UMAP_BY_LIBRARY", 9, 7)
save_dual(p_pca_lib, "PCA_BY_LIBRARY", 9, 7)

cluster_sizes <- as.data.table(table(obj$seurat_clusters))
setnames(cluster_sizes, c("cluster", "cell_number"))
cluster_sizes[, cluster := as.character(cluster)]
write_tsv(cluster_sizes, p("04_tables", "CLUSTER_SIZES_RES0.6.tsv"))
p_size <- ggplot(cluster_sizes, aes(reorder(cluster, -cell_number), cell_number)) +
  geom_col(fill = "#2C7FB8") + geom_text(aes(label = cell_number), vjust = -0.15, size = 3) +
  labs(x = "Cluster (resolution 0.6)", y = "Cells") + theme_sc()
save_dual(p_size, "CLUSTER_SIZE_BARPLOT", 11, 6)

p_comp <- ggplot(cluster_lib, aes(cluster, fraction_within_cluster, fill = library_id)) +
  geom_col() + scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Cluster", y = "Library fraction within cluster") + theme_sc()
save_dual(p_comp, "CLUSTER_BY_LIBRARY_COMPOSITION", 11, 6)

saveRDS(obj, p("06_objects", "03_clustered_seurat.rds"), compress = FALSE)
stage_message("Stage 03 complete: resolution 0.6 clusters=", uniqueN(obj$seurat_clusters),
              "; Cramer's V=", sprintf("%.3f", cramers_v))
