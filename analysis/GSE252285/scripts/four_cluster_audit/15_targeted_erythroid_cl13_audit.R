source(file.path(Sys.getenv("SCRIPT_DIR", "."), "branch_common.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

set.seed(42L)
audit_name <- "cluster_6_9_13_targeted_final_audit"
audit_path <- function(section, ...) file.path(branch_path(section), audit_name, ...)
for (section in c("00_manifest", "04_tables", "05_figures", "07_logs",
                  "08_methods", "09_qc", "10_release")) {
  dir.create(audit_path(section), recursive = TRUE, showWarnings = FALSE)
}
save_audit_dual <- function(plot, stem, width = 10, height = 7, dpi = 300) {
  ggsave(audit_path("05_figures", paste0(stem, ".pdf")), plot = plot,
         width = width, height = height, units = "in", limitsize = FALSE, bg = "white")
  ggsave(audit_path("05_figures", paste0(stem, ".png")), plot = plot,
         width = width, height = height, units = "in", limitsize = FALSE,
         dpi = dpi, device = ragg::agg_png, bg = "white")
}

stage_message("Targeted erythroid and cluster-13 co-expression audit; read-only upstream object")
target_path <- branch_path("06_objects", "02_no_doublet_annotated_tsne.rds")
doublet_path <- file.path(project_dir, "04_tables", "DOUBLET_RESULTS.tsv.gz")
qc_threshold_path <- file.path(project_dir, "04_tables", "QC_THRESHOLDS.tsv")
required <- c(target_path, doublet_path, qc_threshold_path)
if (any(!file.exists(required))) stop("Missing input: ", paste(required[!file.exists(required)], collapse = ", "))

obj <- readRDS(target_path)
DefaultAssay(obj) <- "RNA"
cluster_field <- "RNA_snn_res.0.6"
if (nrow(obj) != 32285L || ncol(obj) != 36111L) stop("Target object dimension invariant failed")
if (!all(c("counts", "data") %in% Layers(obj[["RNA"]]))) stop("Required sparse RNA layers missing")
if (!all(c("pca", "tsne") %in% Reductions(obj)) || "umap" %in% Reductions(obj)) stop("Reduction invariant failed")
cluster_vec <- as.character(obj[[cluster_field, drop = TRUE]])
if (uniqueN(cluster_vec) != 28L) stop("Fixed resolution cluster invariant failed")
target_clusters <- c("6", "9", "13")
target_counts <- data.table(cluster = cluster_vec)[cluster %in% target_clusters, .N, by = cluster]
target_counts <- target_counts[order(as.integer(cluster))]
if (!identical(target_counts$N, c(2359L, 1210L, 838L))) stop("Target cluster cell counts changed")

counts <- LayerData(obj[["RNA"]], layer = "counts")
data_mat <- LayerData(obj[["RNA"]], layer = "data")
meta <- as.data.table(obj[[]], keep.rownames = "barcode")
meta[, cluster := cluster_vec]
meta[, audit_group := fifelse(cluster %in% target_clusters, paste0("cluster_", cluster), "REST")]

historical_doublet <- fread(doublet_path, showProgress = FALSE)
setnames(historical_doublet,
         c("scDblFinder_score", "scDblFinder_class"),
         c("historical_scDblFinder_score", "historical_scDblFinder_class"),
         skip_absent = TRUE)
if (!all(c("barcode", "historical_scDblFinder_score", "historical_scDblFinder_class") %in% names(historical_doublet))) {
  stop("Historical doublet table schema mismatch")
}
meta <- merge(meta, historical_doublet[, .(barcode, historical_scDblFinder_score,
                                            historical_scDblFinder_class)],
              by = "barcode", all.x = TRUE, sort = FALSE)
meta[, object_order__ := match(barcode, colnames(obj))]
setorder(meta, object_order__)
meta[, object_order__ := NULL]
if (!identical(meta$barcode, colnames(obj))) stop("Metadata order invariant failed")

# Existing QC thresholds remain unchanged; descriptive within-library q99 is used only as a doublet-like audit.
q99 <- meta[, .(
  descriptive_q99_nCount = as.numeric(quantile(nCount_RNA, 0.99, names = FALSE)),
  descriptive_q99_nFeature = as.numeric(quantile(nFeature_RNA, 0.99, names = FALSE))
), by = library_id]
meta <- merge(meta, q99, by = "library_id", all.x = TRUE, sort = FALSE)
meta[, object_order__ := match(barcode, colnames(obj))]
setorder(meta, object_order__)
meta[, object_order__ := NULL]
meta[, descriptive_high_count_feature := nCount_RNA > descriptive_q99_nCount &
       nFeature_RNA > descriptive_q99_nFeature]
write_tsv(q99, audit_path("04_tables", "DESCRIPTIVE_LIBRARY_Q99_FOR_DOUBLETLIKE_AUDIT.tsv"))

requested_erythroid <- c("Klf1", "Gata1", "Alas2", "Gypa", "Slc4a1", "Hba-a1", "Hba-a2")
hbb_genes <- grep("^Hbb-", rownames(counts), value = TRUE)
hb_genes <- grep("^Hb[ab]-", rownames(counts), value = TRUE)
erythroid_genes <- unique(c(requested_erythroid, hbb_genes))
feature_audit <- data.table(
  requested_pattern = c(requested_erythroid, "Hbb-*", "hemoglobin_fraction:^Hb[ab]-"),
  matched_features = c(requested_erythroid, paste(hbb_genes, collapse = ","), paste(hb_genes, collapse = ",")),
  status = c(ifelse(requested_erythroid %in% rownames(counts), "FOUND", "NOT_FOUND"),
             ifelse(length(hbb_genes), "FOUND", "NOT_FOUND"),
             ifelse(length(hb_genes), "FOUND", "NOT_FOUND"))
)
write_tsv(feature_audit, audit_path("04_tables", "ERYTHROID_FEATURE_AUDIT.tsv"))
erythroid_present <- intersect(erythroid_genes, rownames(counts))
if (length(intersect(requested_erythroid, erythroid_present)) < 7L || length(hbb_genes) < 1L) {
  stop("Required erythroid genes unavailable")
}

gene_rows <- list()
group_indices <- list(
  `6` = which(cluster_vec == "6"),
  `9` = which(cluster_vec == "9"),
  REST = which(!cluster_vec %in% c("6", "9"))
)
for (group_name in names(group_indices)) {
  idx <- group_indices[[group_name]]
  cluster_total_counts <- sum(meta$nCount_RNA[idx])
  hb_total <- sum(counts[hb_genes, idx, drop = FALSE])
  for (gene in erythroid_present) {
    x_count <- counts[gene, idx, drop = TRUE]
    x_data <- data_mat[gene, idx, drop = TRUE]
    gene_total <- sum(x_count)
    gene_rows[[length(gene_rows) + 1L]] <- data.table(
      group = group_name,
      gene = gene,
      cell_number = length(idx),
      n_expressing = sum(x_count > 0),
      pct_expression = mean(x_count > 0),
      mean_counts_per_cell = mean(x_count),
      mean_log_normalized_expression = mean(x_data),
      cluster_total_UMI_fraction = gene_total / cluster_total_counts,
      hemoglobin_UMI_contribution = if (gene %in% hb_genes && hb_total > 0) gene_total / hb_total else NA_real_
    )
  }
}
erythroid_gene_stats <- rbindlist(gene_rows)
write_tsv(erythroid_gene_stats, audit_path("04_tables", "CL6_CL9_ERYTHROID_GENE_EXPRESSION.tsv"))

# Cell-level hemoglobin fraction and coherent non-hemoglobin erythroid program.
ery_cells <- which(cluster_vec %in% c("6", "9"))
ery_meta <- copy(meta[ery_cells])
ery_meta[, hemoglobin_counts := as.numeric(Matrix::colSums(counts[hb_genes, ery_cells, drop = FALSE]))]
ery_meta[, hemoglobin_fraction := hemoglobin_counts / pmax(nCount_RNA, 1)]
core_non_hb <- c("Klf1", "Gata1", "Alas2", "Gypa", "Slc4a1")
mature_core <- c("Alas2", "Gypa", "Slc4a1")
core_detect <- t(as.matrix(counts[core_non_hb, ery_cells, drop = FALSE] > 0))
colnames(core_detect) <- core_non_hb
ery_meta[, core_non_hb_detected := rowSums(core_detect)]
ery_meta[, mature_core_detected := rowSums(core_detect[, mature_core, drop = FALSE])]
ery_meta[, erythroid_TF_detected := rowSums(core_detect[, c("Klf1", "Gata1"), drop = FALSE])]
ery_meta[, hemoglobin_detected := hemoglobin_counts > 0]
ery_meta[, ambient_only_pattern := hemoglobin_detected & core_non_hb_detected == 0]
ery_meta[, coherent_mature_erythroid := mature_core_detected == length(mature_core) & hemoglobin_detected]

ery_cell_output <- ery_meta[, .(
  barcode, cluster, library_id, nCount_RNA, nFeature_RNA, percent.mt,
  historical_scDblFinder_score, historical_scDblFinder_class,
  descriptive_high_count_feature, hemoglobin_counts, hemoglobin_fraction,
  core_non_hb_detected, mature_core_detected, erythroid_TF_detected,
  hemoglobin_detected, ambient_only_pattern, coherent_mature_erythroid
)]
fwrite(ery_cell_output, audit_path("04_tables", "CL6_CL9_ERYTHROID_CELL_AUDIT.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

ery_summary <- ery_meta[, .(
  cell_number = .N,
  mean_hemoglobin_fraction = mean(hemoglobin_fraction),
  median_hemoglobin_fraction = median(hemoglobin_fraction),
  q25_hemoglobin_fraction = quantile(hemoglobin_fraction, 0.25),
  q75_hemoglobin_fraction = quantile(hemoglobin_fraction, 0.75),
  pct_hemoglobin_fraction_ge_10pct = mean(hemoglobin_fraction >= 0.10),
  pct_hemoglobin_fraction_ge_25pct = mean(hemoglobin_fraction >= 0.25),
  pct_hemoglobin_fraction_ge_50pct = mean(hemoglobin_fraction >= 0.50),
  pct_ambient_only_pattern = mean(ambient_only_pattern),
  pct_coherent_mature_erythroid = mean(coherent_mature_erythroid),
  pct_core_non_hb_ge3 = mean(core_non_hb_detected >= 3),
  pct_any_erythroid_TF = mean(erythroid_TF_detected >= 1),
  median_nCount_RNA = median(nCount_RNA),
  median_nFeature_RNA = median(nFeature_RNA),
  median_percent_mt = median(percent.mt),
  pct_nFeature_lt300 = mean(nFeature_RNA < 300),
  pct_nFeature_lt500 = mean(nFeature_RNA < 500),
  historical_doublet_rate = mean(historical_scDblFinder_class == "doublet", na.rm = TRUE),
  descriptive_high_count_feature_rate = mean(descriptive_high_count_feature)
), by = cluster]
ery_summary <- ery_summary[order(as.integer(cluster))]
write_tsv(ery_summary, audit_path("04_tables", "CL6_CL9_ERYTHROID_PROGRAM_SUMMARY.tsv"))

# Cluster 13 same-cell hematopoietic/endothelial co-expression.
cl13_idx <- which(cluster_vec == "13")
cl13 <- copy(meta[cl13_idx])
coexpr_genes <- c("Ptprc", "Pecam1", "Egfl7", "Cdh5", "Kdr")
if (any(!coexpr_genes %in% rownames(counts))) stop("Co-expression feature missing")
co_count <- t(as.matrix(counts[coexpr_genes, cl13_idx, drop = FALSE] > 0))
colnames(co_count) <- coexpr_genes
co_log <- t(as.matrix(data_mat[coexpr_genes, cl13_idx, drop = FALSE]))
colnames(co_log) <- paste0(coexpr_genes, "_lognorm")
cl13 <- cbind(cl13, as.data.table(co_count), as.data.table(co_log))
endo_genes <- c("Pecam1", "Egfl7", "Cdh5", "Kdr")
cl13[, endothelial_marker_count := rowSums(.SD), .SDcols = endo_genes]
cl13[, endothelial_module_mean_lognorm := rowMeans(.SD),
     .SDcols = paste0(endo_genes, "_lognorm")]
cl13[, coexpression_group := fcase(
  Ptprc & endothelial_marker_count >= 2, "Ptprc+_Endo>=2",
  Ptprc & endothelial_marker_count == 1, "Ptprc+_Endo=1",
  Ptprc & endothelial_marker_count == 0, "Ptprc+_Endo=0",
  !Ptprc & endothelial_marker_count >= 2, "Ptprc-_Endo>=2",
  default = "Ptprc-_Endo<=1"
)]

event_summary <- function(label, flag) data.table(
  event = label, cell_count = sum(flag), percentage = mean(flag)
)
coexpr_summary <- rbindlist(list(
  event_summary("Ptprc+", cl13$Ptprc),
  event_summary("Pecam1+", cl13$Pecam1),
  event_summary("Egfl7+", cl13$Egfl7),
  event_summary("Cdh5+", cl13$Cdh5),
  event_summary("Kdr+", cl13$Kdr),
  event_summary("Ptprc+ & Pecam1+", cl13$Ptprc & cl13$Pecam1),
  event_summary("Ptprc+ & Egfl7+", cl13$Ptprc & cl13$Egfl7),
  event_summary("Ptprc+ & Cdh5+", cl13$Ptprc & cl13$Cdh5),
  event_summary("Ptprc+ & Kdr+", cl13$Ptprc & cl13$Kdr),
  event_summary("Ptprc+ & >=1 endothelial marker", cl13$Ptprc & cl13$endothelial_marker_count >= 1),
  event_summary("Ptprc+ & >=2 endothelial markers", cl13$Ptprc & cl13$endothelial_marker_count >= 2),
  event_summary("Ptprc+ & >=3 endothelial markers", cl13$Ptprc & cl13$endothelial_marker_count >= 3),
  event_summary("Ptprc+ & all 4 endothelial markers", cl13$Ptprc & cl13$endothelial_marker_count == 4),
  event_summary("Ptprc- & >=2 endothelial markers", !cl13$Ptprc & cl13$endothelial_marker_count >= 2)
))
write_tsv(coexpr_summary, audit_path("04_tables", "CL13_COEXPRESSION_SUMMARY.tsv"))

pair_rows <- rbindlist(lapply(endo_genes, function(gene) {
  a <- cl13$Ptprc
  b <- cl13[[gene]]
  expected <- mean(a) * mean(b)
  observed <- mean(a & b)
  data.table(
    marker_pair = paste0("Ptprc+", gene),
    ptprc_pct = mean(a), endothelial_marker_pct = mean(b),
    observed_coexpression_pct = observed,
    independence_expected_pct = expected,
    observed_to_expected_ratio = observed / pmax(expected, 1e-12)
  )
}))
write_tsv(pair_rows, audit_path("04_tables", "CL13_PAIRWISE_COEXPRESSION_ENRICHMENT.tsv"))

coexpr_group_audit <- cl13[, .(
  cell_number = .N,
  fraction_of_cluster13 = .N / nrow(cl13),
  median_nCount_RNA = median(nCount_RNA),
  median_nFeature_RNA = median(nFeature_RNA),
  median_percent_mt = median(percent.mt),
  historical_scDblFinder_available = sum(!is.na(historical_scDblFinder_class)),
  historical_doublet_number = sum(historical_scDblFinder_class == "doublet", na.rm = TRUE),
  historical_doublet_rate = mean(historical_scDblFinder_class == "doublet", na.rm = TRUE),
  median_historical_scDblFinder_score = median(historical_scDblFinder_score, na.rm = TRUE),
  descriptive_high_count_feature_number = sum(descriptive_high_count_feature),
  descriptive_high_count_feature_rate = mean(descriptive_high_count_feature),
  LIB1_fraction = mean(library_id == "LIB1"),
  LIB2_fraction = mean(library_id == "LIB2"),
  LIB3_fraction = mean(library_id == "LIB3")
), by = coexpression_group]
setorder(coexpr_group_audit, -cell_number)
write_tsv(coexpr_group_audit, audit_path("04_tables", "CL13_COEXPRESSION_QC_DOUBLET_AUDIT.tsv"))

cl13_output <- cl13[, c(
  "barcode", "cluster", "library_id", "nCount_RNA", "nFeature_RNA", "percent.mt",
  "historical_scDblFinder_score", "historical_scDblFinder_class",
  "descriptive_high_count_feature", coexpr_genes, paste0(coexpr_genes, "_lognorm"),
  "endothelial_marker_count", "endothelial_module_mean_lognorm", "coexpression_group"
), with = FALSE]
fwrite(cl13_output, audit_path("04_tables", "CL13_COEXPRESSION_CELL_AUDIT.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

# Common target-cluster QC, library and historical-doublet context.
target_meta <- meta[cluster %in% target_clusters]
target_qc <- target_meta[, .(
  cell_number = .N,
  LIB1_fraction = mean(library_id == "LIB1"),
  LIB2_fraction = mean(library_id == "LIB2"),
  LIB3_fraction = mean(library_id == "LIB3"),
  median_nCount_RNA = median(nCount_RNA),
  q25_nCount_RNA = quantile(nCount_RNA, 0.25),
  q75_nCount_RNA = quantile(nCount_RNA, 0.75),
  median_nFeature_RNA = median(nFeature_RNA),
  q25_nFeature_RNA = quantile(nFeature_RNA, 0.25),
  q75_nFeature_RNA = quantile(nFeature_RNA, 0.75),
  median_percent_mt = median(percent.mt),
  q25_percent_mt = quantile(percent.mt, 0.25),
  q75_percent_mt = quantile(percent.mt, 0.75),
  historical_scDblFinder_match_rate = mean(!is.na(historical_scDblFinder_class)),
  historical_doublet_rate = mean(historical_scDblFinder_class == "doublet", na.rm = TRUE),
  median_historical_scDblFinder_score = median(historical_scDblFinder_score, na.rm = TRUE),
  descriptive_high_count_feature_rate = mean(descriptive_high_count_feature)
), by = cluster]
target_qc <- target_qc[order(as.integer(cluster))]
write_tsv(target_qc, audit_path("04_tables", "TARGET_CLUSTERS_QC_LIBRARY_DOUBLET_AUDIT.tsv"))

# Figures.
plot_gene_stats <- erythroid_gene_stats[group %in% c("6", "9", "REST")]
plot_gene_stats[, group := factor(group, levels = c("6", "9", "REST"), labels = c("cluster 6", "cluster 9", "rest"))]
plot_gene_stats[, gene := factor(gene, levels = rev(erythroid_present))]
p_ery_dot <- ggplot(plot_gene_stats, aes(group, gene, size = 100 * pct_expression,
                                         colour = mean_log_normalized_expression)) +
  geom_point() +
  scale_size_area(max_size = 8, limits = c(0, 100), breaks = c(25, 50, 75, 100)) +
  scale_colour_viridis_c(option = "C") +
  labs(x = NULL, y = NULL, size = "% expressing", colour = "Mean log-normalized",
       title = "Erythroid program audit",
       subtitle = "Non-hemoglobin regulators and membrane/heme genes distinguish a coherent program from ambient globin") +
  theme_branch()
save_audit_dual(p_ery_dot, "CL6_CL9_ERYTHROID_GENE_DOTPLOT", 9, 8)

hb_plot_dt <- copy(meta)
hb_plot_dt[, hemoglobin_counts := as.numeric(Matrix::colSums(counts[hb_genes, , drop = FALSE]))]
hb_plot_dt[, hemoglobin_fraction := hemoglobin_counts / pmax(nCount_RNA, 1)]
hb_plot_dt[, plot_group := fcase(cluster == "6", "cluster 6", cluster == "9", "cluster 9",
                                cluster == "13", "cluster 13", default = "rest")]
hb_plot_dt[, plot_group := factor(plot_group, levels = c("cluster 6", "cluster 9", "cluster 13", "rest"))]
p_hb <- ggplot(hb_plot_dt, aes(hemoglobin_fraction, colour = plot_group, fill = plot_group)) +
  geom_density(alpha = 0.12, linewidth = 0.8, adjust = 1.2) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_colour_manual(values = c("cluster 6" = "#D55E00", "cluster 9" = "#CC79A7",
                                 "cluster 13" = "#0072B2", "rest" = "grey45")) +
  scale_fill_manual(values = c("cluster 6" = "#D55E00", "cluster 9" = "#CC79A7",
                               "cluster 13" = "#0072B2", "rest" = "grey70")) +
  labs(x = "Hemoglobin UMI fraction", y = "Density", colour = NULL, fill = NULL,
       title = "Per-cell hemoglobin fraction") + theme_branch()
save_audit_dual(p_hb, "CL6_CL9_HEMOGLOBIN_FRACTION_DENSITY", 10, 7)

group_levels <- c("Ptprc+_Endo>=2", "Ptprc+_Endo=1", "Ptprc+_Endo=0",
                  "Ptprc-_Endo>=2", "Ptprc-_Endo<=1")
coexpr_group_audit[, coexpression_group := factor(coexpression_group, levels = group_levels)]
p_co_count <- ggplot(coexpr_group_audit, aes(coexpression_group, cell_number, fill = coexpression_group)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = paste0(cell_number, "\n", sprintf("%.1f%%", 100 * fraction_of_cluster13))),
            vjust = -0.2, size = 3.5) +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(x = NULL, y = "Cells", title = "Cluster 13 same-cell Ptprc/endothelial co-expression") +
  theme_branch() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
p_co_doublet <- ggplot(coexpr_group_audit, aes(coexpression_group, 100 * historical_doublet_rate,
                                               fill = coexpression_group)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * historical_doublet_rate)), vjust = -0.2, size = 3.5) +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(x = NULL, y = "Historical scDblFinder doublets (%)",
       title = "Historical doublet-call enrichment (audit only)") +
  theme_branch() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_audit_dual(p_co_count / p_co_doublet, "CL13_COEXPRESSION_DOUBLET_AUDIT", 11, 10)

cl13_long <- melt(cl13,
  id.vars = c("barcode", "coexpression_group"),
  measure.vars = paste0(coexpr_genes, "_lognorm"),
  variable.name = "gene", value.name = "lognorm")
cl13_long[, gene := sub("_lognorm$", "", gene)]
cl13_dot <- cl13_long[, .(pct_expression = mean(lognorm > 0), mean_lognorm = mean(lognorm)),
                        by = .(coexpression_group, gene)]
cl13_dot[, coexpression_group := factor(coexpression_group, levels = group_levels)]
cl13_dot[, gene := factor(gene, levels = rev(coexpr_genes))]
p_cl13_dot <- ggplot(cl13_dot, aes(coexpression_group, gene, size = 100 * pct_expression,
                                   colour = mean_lognorm)) +
  geom_point() + scale_size_area(max_size = 10, limits = c(0, 100)) +
  scale_colour_viridis_c(option = "C") +
  labs(x = NULL, y = NULL, size = "% expressing", colour = "Mean log-normalized",
       title = "Cluster 13 co-expression groups") +
  theme_branch() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_audit_dual(p_cl13_dot, "CL13_COEXPRESSION_GROUP_DOTPLOT", 10, 6)

lib_dt <- target_meta[, .N, by = .(cluster, library_id)]
lib_dt[, fraction := N / sum(N), by = cluster]
lib_dt[, cluster := factor(cluster, levels = target_clusters)]
p_lib <- ggplot(lib_dt, aes(cluster, fraction, fill = library_id)) +
  geom_col() + scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(LIB1 = "#0072B2", LIB2 = "#D55E00", LIB3 = "#009E73")) +
  labs(x = "Cluster", y = "Cell fraction", fill = NULL,
       title = "Target-cluster library composition") + theme_branch()
save_audit_dual(p_lib, "TARGET_CLUSTERS_LIBRARY_COMPOSITION", 8, 6)

qc_long <- melt(target_meta[, .(cluster, nCount_RNA, nFeature_RNA, percent.mt)],
                id.vars = "cluster", variable.name = "metric", value.name = "value")
qc_long[, cluster := factor(cluster, levels = target_clusters)]
p_qc <- ggplot(qc_long, aes(cluster, value, fill = cluster)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_viridis_d(option = "D", end = 0.85) +
  labs(x = "Cluster", y = NULL, title = "Target-cluster QC distributions") + theme_branch()
save_audit_dual(p_qc, "TARGET_CLUSTERS_QC_VIOLIN", 12, 5)

parameters <- data.table(
  parameter = c("target_object", "cluster_field", "target_clusters", "target_cell_counts",
                "RNA_counts_layer", "RNA_data_layer", "hemoglobin_gene_rule",
                "historical_scDblFinder_use", "descriptive_doubletlike_rule",
                "QC_rerun", "PCA_rerun", "clustering_rerun", "object_modified",
                "SGN_recovery_attempted"),
  value = c(target_path, cluster_field, paste(target_clusters, collapse = ","),
            paste(target_counts$N, collapse = ","), "read_only", "read_only", "^Hb[ab]-",
            "AUDIT_ONLY_NO_THRESHOLD_CHANGE", "within-library nCount and nFeature both above q99",
            "NO", "NO", "NO", "NO", "NO")
)
write_tsv(parameters, audit_path("00_manifest", "TARGETED_AUDIT_PARAMETERS.tsv"))
file.copy(qc_threshold_path, audit_path("00_manifest", "ORIGINAL_QC_THRESHOLDS_READ_ONLY.tsv"), overwrite = TRUE)
capture.output(sessionInfo(), file = audit_path("07_logs", "sessionInfo.txt"))
stage_message("Targeted evidence audit complete")
