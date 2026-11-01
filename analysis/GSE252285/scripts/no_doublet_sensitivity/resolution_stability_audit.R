source(file.path(Sys.getenv("SCRIPT_DIR", "."), "branch_common.R"))
suppressPackageStartupMessages({
  library(Seurat)
})

stage_message("Post hoc resolution stability/technical audit; no clustering is rerun")
obj <- readRDS(branch_path("06_objects", "02_no_doublet_annotated_tsne.rds"))
meta <- as.data.table(obj[[]], keep.rownames = "barcode")
res_map <- c(`0.4` = "RNA_snn_res.0.4", `0.6` = "RNA_snn_res.0.6",
             `0.8` = "RNA_snn_res.0.8", `1.0` = "RNA_snn_res.1")
if (!all(res_map %in% names(meta))) stop("One or more audited resolution fields are missing")

pairs <- list(c("0.4", "0.6"), c("0.6", "0.8"), c("0.8", "1.0"))
transition_rows <- list()
summary_rows <- list()
for (pair in pairs) {
  low_res <- pair[1]
  high_res <- pair[2]
  low <- as.character(meta[[res_map[[low_res]]]])
  high <- as.character(meta[[res_map[[high_res]]]])
  dt <- data.table(low_cluster = low, high_cluster = high)[, .N,
    by = .(low_cluster, high_cluster)]
  dt[, fraction_of_low_cluster := N / sum(N), by = low_cluster]
  dt[, fraction_of_high_cluster := N / sum(N), by = high_cluster]
  dt[, `:=`(low_resolution = low_res, high_resolution = high_res)]
  transition_rows[[length(transition_rows) + 1L]] <- dt

  high_parent <- dt[, .SD[which.max(N)], by = high_cluster]
  low_children <- dt[N >= 20 & fraction_of_low_cluster >= 0.05,
                     .(substantial_children = uniqueN(high_cluster)), by = low_cluster]
  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    low_resolution = low_res,
    high_resolution = high_res,
    low_cluster_number = uniqueN(low),
    high_cluster_number = uniqueN(high),
    additional_cluster_number = uniqueN(high) - uniqueN(low),
    adjusted_rand_index = adjusted_rand_index(low, high),
    median_dominant_parent_fraction_of_high_cluster = median(high_parent$fraction_of_high_cluster),
    minimum_dominant_parent_fraction_of_high_cluster = min(high_parent$fraction_of_high_cluster),
    lower_clusters_with_two_or_more_substantial_children = sum(low_children$substantial_children >= 2),
    interpretation = "Higher resolution is an audit only; marker coherence is required before considering a change"
  )
}
transitions <- rbindlist(transition_rows)
stability <- rbindlist(summary_rows)
write_tsv(transitions, branch_path("04_tables", "RESOLUTION_PARENT_CHILD_TRANSITIONS.tsv"))
write_tsv(stability, branch_path("04_tables", "RESOLUTION_STABILITY_AUDIT.tsv"))

technical_rows <- list()
global_medians <- c(nCount_RNA = median(meta$nCount_RNA),
                    nFeature_RNA = median(meta$nFeature_RNA),
                    percent.mt = median(meta$percent.mt))
global_mads <- c(nCount_RNA = mad(meta$nCount_RNA, constant = 1.4826),
                 nFeature_RNA = mad(meta$nFeature_RNA, constant = 1.4826),
                 percent.mt = mad(meta$percent.mt, constant = 1.4826))
for (res in names(res_map)) {
  field <- res_map[[res]]
  dt <- meta[, .(
    cell_number = .N,
    median_nCount_RNA = as.numeric(median(nCount_RNA)),
    median_nFeature_RNA = as.numeric(median(nFeature_RNA)),
    median_percent_mt = as.numeric(median(percent.mt)),
    max_library_fraction = max(table(library_id)) / .N,
    dominant_library = names(which.max(table(library_id)))
  ), by = .(cluster = as.character(get(field)))]
  dt[, `:=`(
    resolution = res,
    robust_z_median_nCount = (median_nCount_RNA - global_medians[["nCount_RNA"]]) /
      max(global_mads[["nCount_RNA"]], 1e-12),
    robust_z_median_nFeature = (median_nFeature_RNA - global_medians[["nFeature_RNA"]]) /
      max(global_mads[["nFeature_RNA"]], 1e-12),
    robust_z_median_percent_mt = (median_percent_mt - global_medians[["percent.mt"]]) /
      max(global_mads[["percent.mt"]], 1e-12)
  )]
  dt[, technical_warning := ifelse(
    abs(robust_z_median_nCount) >= 3 | abs(robust_z_median_nFeature) >= 3 |
      abs(robust_z_median_percent_mt) >= 3 | max_library_fraction >= 0.90,
    "POTENTIAL_TECHNICAL_OR_LIBRARY_DRIVER", "NONE"
  )]
  technical_rows[[length(technical_rows) + 1L]] <- dt
}
technical <- rbindlist(technical_rows)
setcolorder(technical, c("resolution", "cluster", setdiff(names(technical), c("resolution", "cluster"))))
write_tsv(technical, branch_path("04_tables", "TECHNICAL_COVARIATE_BY_RESOLUTION.tsv"))

plot_transition <- function(low_res, high_res, stem) {
  dt <- transitions[low_resolution == low_res & high_resolution == high_res]
  dt[, low_cluster := factor(low_cluster, levels = sort(unique(low_cluster)))]
  dt[, high_cluster := factor(high_cluster, levels = rev(sort(unique(high_cluster))))]
  p <- ggplot(dt, aes(low_cluster, high_cluster, fill = fraction_of_low_cluster)) +
    geom_tile(colour = "grey85", linewidth = 0.1) +
    scale_fill_viridis_c(option = "C", limits = c(0, 1)) +
    labs(x = paste0("Cluster at resolution ", low_res),
         y = paste0("Cluster at resolution ", high_res),
         fill = "fraction of\nlower cluster",
         title = paste0("Resolution transition: ", low_res, " to ", high_res)) +
    theme_branch()
  save_branch_dual(p, stem, 12, 9)
}
plot_transition("0.6", "0.8", "RESOLUTION_TRANSITION_0.6_TO_0.8")
plot_transition("0.8", "1.0", "RESOLUTION_TRANSITION_0.8_TO_1.0")

audit_lines <- c(
  "# Resolution selection audit",
  "",
  "The primary resolution remains the pre-specified 0.6. No PCA, graph, embedding, or clustering was rerun in this audit.",
  "",
  paste0("- Cluster counts: 0.4=", uniqueN(meta[[res_map[["0.4"]]]]),
         "; 0.6=", uniqueN(meta[[res_map[["0.6"]]]]),
         "; 0.8=", uniqueN(meta[[res_map[["0.8"]]]]),
         "; 1.0=", uniqueN(meta[[res_map[["1.0"]]]])),
  paste0("- Adjacent-resolution ARI: ", paste0(stability$low_resolution, "->",
         stability$high_resolution, "=", sprintf("%.3f", stability$adjusted_rand_index),
         collapse = "; ")),
  paste0("- Clusters with technical/library warning: ",
         paste(technical[technical_warning != "NONE", paste0(resolution, ":", cluster)],
               collapse = ", ")),
  "",
  "A higher resolution is not preferred merely because it yields more clusters. A change would require independent marker coherence and non-technical reproducibility for each split; those alternative-resolution marker analyses were not used to replace the fixed main result."
)
writeLines(audit_lines, branch_path("10_release", "RESOLUTION_SELECTION_AUDIT.md"))
stage_message("Resolution audit complete")
