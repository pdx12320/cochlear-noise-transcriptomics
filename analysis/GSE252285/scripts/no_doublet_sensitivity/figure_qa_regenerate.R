source(file.path(Sys.getenv("SCRIPT_DIR", "."), "branch_common.R"))
suppressPackageStartupMessages({
  library(Seurat)
})

stage_message("Figure-only QA regeneration; no analysis is recomputed")
obj <- readRDS(branch_path("06_objects", "02_no_doublet_annotated_tsne.rds"))
canonical_present <- intersect(unique(unlist(canonical_panels, use.names = FALSE)), rownames(obj))
cluster_levels <- as.character(sort(as.integer(unique(as.character(obj$seurat_clusters_no_doublet)))))
obj$cluster_numeric_plot <- factor(as.character(obj$seurat_clusters_no_doublet),
                                   levels = cluster_levels)
p_canonical <- DotPlot(obj, features = canonical_present, group.by = "cluster_numeric_plot") +
  labs(title = "Canonical lineage markers by no-doublet-filter cluster",
       subtitle = "Dot size: percent expressing; colour: scaled average expression",
       x = "Canonical marker", y = "Cluster") +
  theme_branch() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7))
save_branch_dual(p_canonical, "CANONICAL_MARKER_DOTPLOT_BY_CLUSTER", 18, 9)

transitions <- fread(branch_path("04_tables", "RESOLUTION_PARENT_CHILD_TRANSITIONS.tsv"))
plot_transition_numeric <- function(low_res, high_res, stem) {
  dt <- copy(transitions[low_resolution == low_res & high_resolution == high_res])
  low_levels <- as.character(sort(as.integer(unique(dt$low_cluster))))
  high_levels <- as.character(sort(as.integer(unique(dt$high_cluster))))
  dt[, low_cluster := factor(as.character(low_cluster), levels = low_levels)]
  dt[, high_cluster := factor(as.character(high_cluster), levels = rev(high_levels))]
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
plot_transition_numeric("0.6", "0.8", "RESOLUTION_TRANSITION_0.6_TO_0.8")
plot_transition_numeric("0.8", "1.0", "RESOLUTION_TRANSITION_0.8_TO_1.0")

write_tsv(data.table(
  figure = c("CANONICAL_MARKER_DOTPLOT_BY_CLUSTER",
             "RESOLUTION_TRANSITION_0.6_TO_0.8",
             "RESOLUTION_TRANSITION_0.8_TO_1.0"),
  issue_found = c("x-axis labels not rotated after theme application",
                  "cluster axis sorted lexicographically",
                  "cluster axis sorted lexicographically"),
  correction = c("explicit 60-degree label rotation and numeric cluster ordering",
                 "numeric cluster ordering", "numeric cluster ordering"),
  analysis_recomputed = "NO"
), branch_path("09_qc", "MANUAL_FIGURE_QA_CORRECTIONS.tsv"))
stage_message("Figure-only QA regeneration complete")
