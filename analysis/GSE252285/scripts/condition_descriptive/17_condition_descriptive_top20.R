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
analysis_name <- "condition_descriptive_top20_resolved_geo_mapping"
out_path <- function(section, ...) file.path(branch_path(section), analysis_name, ...)
for (section in c("00_manifest", "04_tables", "05_figures", "07_logs",
                  "08_methods", "09_qc", "10_release")) {
  dir.create(out_path(section), recursive = TRUE, showWarnings = FALSE)
}

save_dual <- function(plot, stem, width = 10, height = 7, dpi = 300) {
  ggsave(out_path("05_figures", paste0(stem, ".pdf")), plot = plot,
         width = width, height = height, units = "in", limitsize = FALSE, bg = "white")
  ggsave(out_path("05_figures", paste0(stem, ".png")), plot = plot,
         width = width, height = height, units = "in", limitsize = FALSE,
         dpi = dpi, device = ragg::agg_png, bg = "white")
}

stage_message("Resolved-condition descriptive Top20 analysis")
target_path <- branch_path("06_objects", "02_no_doublet_annotated_tsne.rds")
map_template_path <- file.path(project_dir, "00_manifest", "library_map_template.tsv")
if (!file.exists(target_path)) stop("Missing fixed object: ", target_path)
if (!file.exists(map_template_path)) stop("Missing original library map template")

obj <- readRDS(target_path)
DefaultAssay(obj) <- "RNA"
cluster_field <- "RNA_snn_res.0.6"
broad_field <- "broad_annotation_no_doublet_sensitivity"
if (nrow(obj) != 32285L || ncol(obj) != 36111L) stop("Fixed object dimension invariant failed")
if (!all(c(cluster_field, broad_field, "library_id", "barcode_suffix", "condition") %in%
         colnames(obj[[]]))) stop("Required metadata field missing")
if (!all(c("counts", "data") %in% Layers(obj[["RNA"]]))) stop("RNA counts/data layers missing")
if (!all(c("pca", "tsne") %in% Reductions(obj)) || "umap" %in% Reductions(obj)) {
  stop("Fixed PCA/t-SNE-only branch invariant failed")
}

meta <- as.data.table(obj[[]], keep.rownames = "cell_name")
meta[, parsed_suffix := sub("^.*(-[123])$", "\\1", cell_name)]
meta[!grepl("-[123]$", cell_name), parsed_suffix := NA_character_]
if (anyNA(meta$parsed_suffix) || !identical(meta$parsed_suffix, as.character(meta$barcode_suffix))) {
  stop("Cell-name suffix and barcode_suffix metadata do not agree")
}
meta[, cluster := as.character(get(cluster_field))]
meta[, broad_cell_type := as.character(get(broad_field))]
if (uniqueN(meta$cluster) != 28L || !setequal(unique(meta$cluster), as.character(0:27))) {
  stop("Expected fixed clusters 0-27")
}
if (any(as.character(meta$condition) != "UNKNOWN")) {
  stop("Legacy condition field was expected to remain UNKNOWN")
}

# User-resolved condition mapping is keyed to observed GEO barcode suffix, never inferred from LIB numbering.
user_condition_map <- data.table(
  geo_library_key = c("GEO_-3", "GEO_-1", "GEO_-2"),
  barcode_suffix = c("-3", "-1", "-2"),
  condition_resolved = c("Control", "HHL", "NIHL"),
  mapping_source = "USER_SPECIFIED_GEO_SUFFIX"
)
observed_crosswalk <- meta[, .(cell_count = .N),
                           by = .(barcode_suffix = as.character(barcode_suffix),
                                  library_id = as.character(library_id))]
if (observed_crosswalk[, uniqueN(library_id), by = barcode_suffix][, any(V1 != 1L)] ||
    observed_crosswalk[, uniqueN(barcode_suffix), by = library_id][, any(V1 != 1L)]) {
  stop("Object suffix-library crosswalk is not one-to-one")
}
template_map <- fread(map_template_path)[, .(barcode_suffix = as.character(suffix),
                                             library_id_template = library_id,
                                             legacy_condition_template = condition)]
mapping_audit <- merge(observed_crosswalk, template_map, by = "barcode_suffix", all = TRUE)
mapping_audit <- merge(mapping_audit, user_condition_map, by = "barcode_suffix", all = TRUE)
mapping_audit[, object_template_library_match := library_id == library_id_template]
mapping_audit[, mapping_complete := !is.na(condition_resolved) &
                !is.na(object_template_library_match) & object_template_library_match]
mapping_audit[, CONDITION_MAPPING_RESOLVED := "YES"]
mapping_audit[, condition_order := match(condition_resolved, c("Control", "HHL", "NIHL"))]
setorder(mapping_audit, condition_order)
mapping_audit[, condition_order := NULL]
if (nrow(mapping_audit) != 3L || any(!mapping_audit$mapping_complete) ||
    sum(mapping_audit$cell_count) != 36111L) stop("Resolved mapping audit failed")

meta[, geo_library_key := paste0("GEO_", as.character(barcode_suffix))]
meta[user_condition_map, on = "geo_library_key",
     condition_resolved := i.condition_resolved]
condition_levels <- c("Control", "HHL", "NIHL")
if (anyNA(meta$condition_resolved) || !setequal(unique(meta$condition_resolved), condition_levels)) {
  stop("Not all cells received a resolved condition")
}

broad_order_preferred <- c("Immune_myeloid", "Fibrocyte_mesenchymal", "Endothelial_vascular",
                           "Glial", "Supporting_epithelial", "Hair_sensory", "UNRESOLVED")
broad_levels <- c(intersect(broad_order_preferred, unique(meta$broad_cell_type)),
                  setdiff(sort(unique(meta$broad_cell_type)), broad_order_preferred))
expected_broad_counts <- data.table(
  broad_cell_type = broad_order_preferred,
  expected_n = c(27676L, 2997L, 275L, 264L, 149L, 24L, 4726L)
)
observed_broad_counts <- meta[, .(observed_n = .N), by = broad_cell_type]
broad_count_check <- merge(expected_broad_counts, observed_broad_counts,
                           by = "broad_cell_type", all = TRUE)
broad_count_check[, pass := expected_n == observed_n]
if (anyNA(broad_count_check$pass) || any(!broad_count_check$pass)) {
  stop("Existing broad annotation count invariant failed")
}

meta[, condition_resolved := factor(condition_resolved, levels = condition_levels)]
meta[, broad_cell_type := factor(broad_cell_type, levels = broad_levels)]
meta[, cluster := factor(cluster, levels = as.character(0:27))]
counts <- LayerData(obj[["RNA"]], layer = "counts")
data_mat <- LayerData(obj[["RNA"]], layer = "data")
if (!identical(meta$cell_name, colnames(counts)) || !identical(colnames(counts), colnames(data_mat))) {
  stop("Metadata and RNA matrix column order disagree")
}

top20 <- top20_fixed
top20_rank <- setNames(seq_along(top20), top20)
feature_audit <- data.table(
  top20_rank = seq_along(top20), gene = top20,
  status = fifelse(top20 %in% rownames(counts), "FOUND", "NOT_FOUND")
)
if (any(feature_audit$status != "FOUND")) stop("Frozen Top20 feature missing")

write_tsv(mapping_audit, out_path("00_manifest", "GEO_CONDITION_MAPPING_AUDIT.tsv"))
write_tsv(broad_count_check, out_path("04_tables", "EXISTING_BROAD_ANNOTATION_COUNT_CHECK.tsv"))
write_tsv(feature_audit, out_path("04_tables", "TOP20_FEATURE_AUDIT_CONDITION.tsv"))

# Condition x broad composition with explicit zero combinations.
composition_grid <- CJ(condition = condition_levels, broad_cell_type = broad_levels, unique = TRUE)
composition_observed <- meta[, .(cell_count = .N),
                             by = .(condition = as.character(condition_resolved),
                                    broad_cell_type = as.character(broad_cell_type))]
composition <- merge(composition_grid, composition_observed,
                     by = c("condition", "broad_cell_type"), all.x = TRUE)
composition[is.na(cell_count), cell_count := 0L]
condition_totals <- meta[, .(total_cells = .N),
                         by = .(condition = as.character(condition_resolved))]
condition_totals <- condition_totals[match(condition_levels, condition)]
composition <- merge(composition, condition_totals, by = "condition", all.x = TRUE)
composition[, `:=`(
  cell_type_fraction_within_condition = cell_count / total_cells,
  cell_type_percentage_within_condition = 100 * cell_count / total_cells,
  cell_count_flag = fifelse(cell_count < 20L, "LOW_CELL_COUNT", "ADEQUATE_CELL_COUNT")
)]
composition[, condition := factor(condition, levels = condition_levels)]
composition[, broad_cell_type := factor(broad_cell_type, levels = broad_levels)]
setorder(composition, condition, broad_cell_type)
write_tsv(condition_totals, out_path("04_tables", "CONDITION_TOTAL_CELLS.tsv"))
write_tsv(composition, out_path("04_tables", "CONDITION_BROAD_CELLTYPE_COMPOSITION.tsv"))

row_medians_sparse_subset <- function(x) {
  vapply(seq_len(nrow(x)), function(i) median(as.numeric(x[i, , drop = TRUE])), numeric(1))
}

calculate_expression_groups <- function(group_grid, index_function, group_columns,
                                        group_fraction_table = NULL) {
  rows <- vector("list", nrow(group_grid))
  for (i in seq_len(nrow(group_grid))) {
    group_row <- group_grid[i]
    idx <- index_function(group_row)
    if (length(idx)) {
      x_counts <- counts[top20, idx, drop = FALSE]
      x_data <- data_mat[top20, idx, drop = FALSE]
      n_expressing <- as.integer(Matrix::rowSums(x_counts > 0))
      result <- data.table(
        top20_rank = seq_along(top20),
        gene = top20,
        n_total_cells = length(idx),
        n_expressing = n_expressing,
        pct_expressing = 100 * n_expressing / length(idx),
        mean_expression = as.numeric(Matrix::rowMeans(x_data)),
        median_expression = row_medians_sparse_subset(x_data),
        cell_count_flag = fifelse(length(idx) < 20L, "LOW_CELL_COUNT", "ADEQUATE_CELL_COUNT")
      )
    } else {
      result <- data.table(
        top20_rank = seq_along(top20), gene = top20, n_total_cells = 0L,
        n_expressing = 0L, pct_expressing = NA_real_, mean_expression = NA_real_,
        median_expression = NA_real_, cell_count_flag = "LOW_CELL_COUNT"
      )
    }
    for (column in group_columns) result[, (column) := as.character(group_row[[column]])]
    rows[[i]] <- result
  }
  answer <- rbindlist(rows, use.names = TRUE)
  setcolorder(answer, c("top20_rank", "gene", group_columns, "n_total_cells",
                        "n_expressing", "pct_expressing", "mean_expression",
                        "median_expression", "cell_count_flag"))
  answer
}

broad_grid <- CJ(condition = condition_levels, broad_cell_type = broad_levels, unique = TRUE)
broad_stats <- calculate_expression_groups(
  broad_grid,
  function(group_row) which(as.character(meta$condition_resolved) == group_row$condition &
                              as.character(meta$broad_cell_type) == group_row$broad_cell_type),
  c("condition", "broad_cell_type")
)
broad_stats <- merge(
  broad_stats,
  composition[, .(condition = as.character(condition),
                  broad_cell_type = as.character(broad_cell_type),
                  cell_type_cell_fraction = cell_type_fraction_within_condition,
                  cell_type_cell_percentage = cell_type_percentage_within_condition)],
  by = c("condition", "broad_cell_type"), all.x = TRUE
)
broad_stats[, condition := factor(condition, levels = condition_levels)]
broad_stats[, broad_cell_type := factor(broad_cell_type, levels = broad_levels)]
broad_stats[, gene := factor(gene, levels = top20)]
setorder(broad_stats, top20_rank, broad_cell_type, condition)
setcolorder(broad_stats, c("top20_rank", "gene", "broad_cell_type", "condition",
                           "n_total_cells", "n_expressing", "pct_expressing",
                           "mean_expression", "median_expression",
                           "cell_type_cell_fraction", "cell_type_cell_percentage",
                           "cell_count_flag"))
write_tsv(broad_stats, out_path("04_tables", "TOP20_BY_BROAD_CELLTYPE_CONDITION.tsv"))

cluster_grid <- CJ(condition = condition_levels, cluster = as.character(0:27), unique = TRUE)
cluster_stats <- calculate_expression_groups(
  cluster_grid,
  function(group_row) which(as.character(meta$condition_resolved) == group_row$condition &
                              as.character(meta$cluster) == group_row$cluster),
  c("condition", "cluster")
)
cluster_stats[, condition := factor(condition, levels = condition_levels)]
cluster_stats[, cluster := factor(cluster, levels = as.character(0:27))]
cluster_stats[, gene := factor(gene, levels = top20)]
setorder(cluster_stats, top20_rank, cluster, condition)
write_tsv(cluster_stats, out_path("04_tables", "TOP20_BY_28CLUSTER_CONDITION.tsv"))

# Pooled Top1/Top2 broad cell types per fixed Top20, followed by three-condition detection proportions.
pooled_localization <- broad_stats[, .(
  pooled_n_total_cells = sum(n_total_cells),
  pooled_n_expressing = sum(n_expressing),
  pooled_pct_expressing = 100 * sum(n_expressing) / sum(n_total_cells),
  pooled_mean_expression = weighted.mean(mean_expression, n_total_cells)
), by = .(top20_rank, gene = as.character(gene), broad_cell_type = as.character(broad_cell_type))]
setorder(pooled_localization, top20_rank, -pooled_mean_expression, -pooled_pct_expressing)
top2_localization <- pooled_localization[, head(.SD, 2L), by = .(top20_rank, gene)]
top2_localization[, localization_rank := seq_len(.N), by = .(top20_rank, gene)]

condition_wide <- dcast(
  broad_stats[, .(gene = as.character(gene), broad_cell_type = as.character(broad_cell_type),
                  condition = as.character(condition), pct_expressing,
                  mean_expression, n_total_cells, cell_count_flag)],
  gene + broad_cell_type ~ condition,
  value.var = c("pct_expressing", "mean_expression", "n_total_cells", "cell_count_flag")
)
top2_localization <- merge(top2_localization, condition_wide,
                           by = c("gene", "broad_cell_type"), all.x = TRUE)
setorder(top2_localization, top20_rank, localization_rank)
write_tsv(top2_localization, out_path("04_tables", "TOP20_TOP1_TOP2_CELLTYPE_CONDITION.tsv"))

condition_top1 <- broad_stats[
  order(top20_rank, condition, -mean_expression, -pct_expressing),
  .(broad_cell_type = as.character(broad_cell_type[1]),
    n_total_cells = n_total_cells[1],
    pct_expressing = pct_expressing[1],
    mean_expression = mean_expression[1]),
  by = .(top20_rank, gene = as.character(gene), condition = as.character(condition))
]
stability <- rbindlist(lapply(seq_along(top20), function(i) {
  z <- condition_top1[gene == top20[i]]
  z <- z[match(condition_levels, condition)]
  same_top1 <- uniqueN(as.character(z$broad_cell_type)) == 1L
  supported <- all(z$n_total_cells >= 20L)
  data.table(
    top20_rank = i,
    gene = top20[i],
    Control_top1 = as.character(z$broad_cell_type[1]),
    Control_top1_n = z$n_total_cells[1],
    Control_top1_pct = z$pct_expressing[1],
    HHL_top1 = as.character(z$broad_cell_type[2]),
    HHL_top1_n = z$n_total_cells[2],
    HHL_top1_pct = z$pct_expressing[2],
    NIHL_top1 = as.character(z$broad_cell_type[3]),
    NIHL_top1_n = z$n_total_cells[3],
    NIHL_top1_pct = z$pct_expressing[3],
    top1_same_all_conditions = same_top1,
    all_top1_groups_n_ge_20 = supported,
    localization_stability = fifelse(same_top1 & supported, "STABLE_SUPPORTED",
                                     fifelse(same_top1, "STABLE_LOW_CELL_COUNT",
                                             "CONDITION_VARIABLE"))
  )
}))
write_tsv(stability, out_path("04_tables", "TOP20_LOCALIZATION_STABILITY.tsv"))

candidate_definitions <- data.table(
  candidate_order = 1:4,
  gene = c("Apold1", "Rasd1", "Depp1", "Errfi1"),
  broad_cell_type = c("Endothelial_vascular", "Endothelial_vascular",
                      "Endothelial_vascular", "Fibrocyte_mesenchymal"),
  candidate_context = c("Apold1 | Endothelial_vascular", "Rasd1 | Endothelial_vascular",
                        "Depp1 | Endothelial_vascular", "Errfi1 | Fibrocyte_mesenchymal")
)
candidate_stats <- merge(
  broad_stats[, .(gene = as.character(gene), broad_cell_type = as.character(broad_cell_type),
                  condition = as.character(condition), n_total_cells, n_expressing,
                  pct_expressing, mean_expression, median_expression,
                  cell_type_cell_fraction, cell_count_flag)],
  candidate_definitions, by = c("gene", "broad_cell_type")
)
candidate_stats[, condition := factor(condition, levels = condition_levels)]
candidate_stats[, candidate_context := factor(candidate_context,
                                               levels = candidate_definitions$candidate_context)]
setorder(candidate_stats, candidate_order, condition)
write_tsv(candidate_stats, out_path("04_tables", "CANDIDATE_CONTEXT_BY_CONDITION.tsv"))

low_broad <- unique(composition[cell_count < 20L,
                                .(grouping_level = "broad_cell_type",
                                  condition = as.character(condition),
                                  group = as.character(broad_cell_type), n_total_cells = cell_count,
                                  warning = "LOW_CELL_COUNT")])
low_cluster <- unique(cluster_stats[n_total_cells < 20L,
                                    .(grouping_level = "cluster",
                                      condition = as.character(condition),
                                      group = as.character(cluster), n_total_cells,
                                      warning = "LOW_CELL_COUNT")])
low_count_combinations <- rbindlist(list(low_broad, low_cluster), use.names = TRUE)
setorder(low_count_combinations, grouping_level, condition, group)
write_tsv(low_count_combinations, out_path("04_tables", "LOW_CELL_COUNT_COMBINATIONS.tsv"))

# Publication-style descriptive figures; no significance or inferential encodings.
condition_palette <- c(Control = "#0072B2", HHL = "#D55E00", NIHL = "#009E73")
composition_plot <- copy(composition)
composition_plot[, condition := factor(condition, levels = condition_levels)]
composition_plot[, broad_cell_type := factor(broad_cell_type, levels = broad_levels)]
p_stack <- ggplot(composition_plot, aes(condition, cell_type_fraction_within_condition,
                                         fill = broad_cell_type)) +
  geom_col(colour = "white", linewidth = 0.15) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_fill_manual(values = broad_palette[broad_levels], drop = FALSE) +
  labs(x = NULL, y = "Cells within condition (%)", fill = "Broad cell type",
       title = "Descriptive broad cell-type composition by resolved condition",
       subtitle = "GEO_-3=Control; GEO_-1=HHL; GEO_-2=NIHL; no abundance test") +
  theme_branch()
save_dual(p_stack, "CONDITION_BROAD_CELLTYPE_STACKED_BAR", 10, 7)

p_grouped <- ggplot(composition_plot,
                    aes(broad_cell_type, cell_type_percentage_within_condition, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  scale_y_continuous(limits = c(0, NA), labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = condition_palette) +
  labs(x = NULL, y = "Cells within condition (%)", fill = "Condition",
       title = "Broad cell-type fractions by resolved condition",
       subtitle = "Descriptive percentages only; no differential-abundance inference") +
  theme_branch() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_dual(p_grouped, "CONDITION_BROAD_CELLTYPE_GROUPED_BAR", 12, 7)

broad_plot <- copy(broad_stats[n_total_cells > 0L])
broad_plot[, condition := factor(condition, levels = condition_levels)]
broad_plot[, broad_cell_type := factor(broad_cell_type, levels = broad_levels)]
broad_plot[, gene := factor(gene, levels = rev(top20))]
p_broad_dot <- ggplot(broad_plot,
                      aes(broad_cell_type, gene, size = pct_expressing, colour = mean_expression)) +
  geom_point() +
  geom_point(data = broad_plot[cell_count_flag == "LOW_CELL_COUNT"],
             shape = 4, colour = "black", stroke = 0.55, size = 2.4) +
  facet_wrap(~condition, ncol = 1) +
  scale_size_area(max_size = 6.5, limits = c(0, 100), breaks = c(10, 25, 50, 75, 100)) +
  scale_colour_viridis_c(option = "C") +
  labs(x = NULL, y = NULL, size = "% expressing", colour = "Mean log-normalized",
       title = "Frozen Top20 by broad cell type and resolved condition",
       subtitle = "Dot area = detection; colour = mean expression; x overlay = n<20") +
  theme_branch() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_dual(p_broad_dot, "TOP20_BROAD_CELLTYPE_CONDITION_DOTPLOT", 13, 20)

p_broad_heat <- ggplot(broad_plot, aes(broad_cell_type, gene, fill = mean_expression)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  geom_point(data = broad_plot[cell_count_flag == "LOW_CELL_COUNT"],
             shape = 4, colour = "black", stroke = 0.6, size = 2.2) +
  facet_wrap(~condition, ncol = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(x = NULL, y = NULL, fill = "Mean log-normalized",
       title = "Frozen Top20 mean expression by broad cell type",
       subtitle = "Common colour scale across Control/HHL/NIHL; x overlay = n<20") +
  theme_branch() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_dual(p_broad_heat, "TOP20_BROAD_CELLTYPE_CONDITION_HEATMAP", 13, 20)

cluster_plot <- copy(cluster_stats[n_total_cells > 0L])
cluster_plot[, condition := factor(condition, levels = condition_levels)]
cluster_plot[, cluster := factor(cluster, levels = as.character(0:27))]
cluster_plot[, gene := factor(gene, levels = rev(top20))]
p_cluster_dot <- ggplot(cluster_plot,
                        aes(cluster, gene, size = pct_expressing, colour = mean_expression)) +
  geom_point() +
  geom_point(data = cluster_plot[cell_count_flag == "LOW_CELL_COUNT"],
             shape = 4, colour = "black", stroke = 0.45, size = 1.5) +
  facet_wrap(~condition, ncol = 1) +
  scale_size_area(max_size = 5.5, limits = c(0, 100), breaks = c(10, 25, 50, 75, 100)) +
  scale_colour_viridis_c(option = "C") +
  labs(x = "Fixed resolution-0.6 cluster", y = NULL, size = "% expressing",
       colour = "Mean log-normalized",
       title = "Frozen Top20 across all 28 clusters by resolved condition",
       subtitle = "Supplement-level descriptive plot; x overlay = n<20") +
  theme_branch()
save_dual(p_cluster_dot, "TOP20_28CLUSTER_CONDITION_DOTPLOT", 16, 20)

candidate_plot <- copy(candidate_stats)
candidate_plot[, condition := factor(condition, levels = condition_levels)]
p_candidate_pct <- ggplot(candidate_plot, aes(condition, pct_expressing, fill = condition)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%\nn=%d", pct_expressing, n_total_cells)),
            vjust = -0.18, size = 3.0) +
  facet_wrap(~candidate_context, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = condition_palette) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = "Expressing cells (%)",
       title = "Focused candidate detection by resolved condition",
       subtitle = "Each gene facet has its own zero-based y scale; descriptive only") +
  theme_branch()
save_dual(p_candidate_pct, "CANDIDATE_CONDITION_PCT_EXPRESSION_BAR", 11, 8)

p_candidate_mean <- ggplot(candidate_plot, aes(condition, mean_expression, fill = condition)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f\nn=%d", mean_expression, n_total_cells)),
            vjust = -0.18, size = 3.0) +
  facet_wrap(~candidate_context, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = condition_palette) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = "Mean log-normalized expression",
       title = "Focused candidate mean expression by resolved condition",
       subtitle = "Each gene facet has its own zero-based y scale; means include zero-expression cells") +
  theme_branch()
save_dual(p_candidate_mean, "CANDIDATE_CONDITION_MEAN_EXPRESSION_BAR", 11, 8)

for (candidate in candidate_definitions$candidate_context) {
  z <- candidate_plot[as.character(candidate_context) == candidate]
  p1 <- ggplot(z, aes(condition, pct_expressing, fill = condition)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%\n%d/%d", pct_expressing, n_expressing,
                                  n_total_cells)), vjust = -0.18) +
    scale_fill_manual(values = condition_palette) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.22))) +
    labs(x = NULL, y = "Expressing cells (%)", title = candidate) + theme_branch()
  p2 <- ggplot(z, aes(condition, mean_expression, fill = condition)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.3f", mean_expression)), vjust = -0.18) +
    scale_fill_manual(values = condition_palette) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Mean log-normalized expression",
         subtitle = "No p-values or cell-level condition inference") + theme_branch()
  safe_stem <- toupper(gsub("[^A-Za-z0-9]+", "_", candidate))
  save_dual(p1 / p2, paste0("FOCUSED_", safe_stem, "_BY_CONDITION"), 8, 9)
}

figure_descriptions <- data.table(
  figure = c(
    "CONDITION_BROAD_CELLTYPE_STACKED_BAR", "CONDITION_BROAD_CELLTYPE_GROUPED_BAR",
    "TOP20_BROAD_CELLTYPE_CONDITION_DOTPLOT", "TOP20_BROAD_CELLTYPE_CONDITION_HEATMAP",
    "TOP20_28CLUSTER_CONDITION_DOTPLOT", "CANDIDATE_CONDITION_PCT_EXPRESSION_BAR",
    "CANDIDATE_CONDITION_MEAN_EXPRESSION_BAR"
  ),
  alt_text = c(
    "Stacked bars show existing broad cell-type fractions for Control, HHL and NIHL.",
    "Grouped zero-baseline bars compare descriptive broad cell-type percentages across resolved conditions.",
    "Three condition facets show frozen Top20 detection by dot area and mean log-normalized expression by colour; x overlays mark groups with fewer than 20 cells.",
    "Three condition facets show frozen Top20 mean log-normalized expression using one common colour scale; x overlays mark groups with fewer than 20 cells.",
    "Supplement-level three-facet DotPlot shows frozen Top20 across fixed clusters 0-27; x overlays mark cluster-condition groups with fewer than 20 cells.",
    "Four focused candidate facets show descriptive percent expressing by resolved condition with group cell counts.",
    "Four focused candidate facets show descriptive mean log-normalized expression by resolved condition with group cell counts."
  )
)
write_tsv(figure_descriptions, out_path("00_manifest", "FIGURE_DESCRIPTIONS.tsv"))

parameters <- data.table(
  parameter = c("target_object", "fixed_cells", "fixed_cluster_field", "fixed_cluster_number",
                "existing_broad_field", "condition_mapping", "mapping_key", "Top20_order",
                "n_expressing_definition", "pct_expressing_definition",
                "mean_expression_definition", "median_expression_definition",
                "low_cell_count_rule", "condition_DE", "p_values",
                "differential_abundance", "QC_rerun", "PCA_rerun",
                "clustering_rerun", "annotation_rerun", "object_modified",
                "CONDITION_MAPPING_RESOLVED"),
  value = c(target_path, "36111", cluster_field, "28", broad_field,
            "GEO_-3=Control;GEO_-1=HHL;GEO_-2=NIHL", "observed barcode suffix",
            paste(top20, collapse = ","), "RNA counts > 0",
            "100 * n_expressing / n_total_cells",
            "arithmetic mean of Seurat RNA data layer including zeros",
            "median of Seurat RNA data layer including zeros",
            "n_total_cells < 20", "NOT_PERFORMED", "NOT_COMPUTED", "NOT_PERFORMED",
            "NO", "NO", "NO", "NO", "NO", "YES")
)
write_tsv(parameters, out_path("00_manifest", "CONDITION_DESCRIPTIVE_PARAMETERS.tsv"))

methods_lines <- c(
  "# Resolved-condition descriptive Top20 methods",
  "",
  "The existing 36,111-cell robust-QC Seurat object was read only. Fixed RNA_snn_res.0.6 clusters and broad_annotation_no_doublet_sensitivity labels were retained without rerunning or modifying QC, PCA, clustering or annotation.",
  "",
  "Conditions were assigned directly from observed barcode suffix using the user-resolved mapping GEO_-3=Control, GEO_-1=HHL and GEO_-2=NIHL. The object cell-name suffix, barcode_suffix metadata, library_id and original library_map_template.tsv were required to form a complete one-to-one audit. The legacy object condition field remained UNKNOWN and was not overwritten.",
  "",
  "For every fixed Top20 gene, expression was summarized for each resolved condition x existing broad cell type and condition x fixed cluster. n_expressing was defined as RNA counts > 0. pct_expressing was 100*n_expressing/n_total_cells. Mean and median expression used the Seurat RNA data layer and included zero-expression cells. Cell-type fraction was cell count divided by the total cells in that condition.",
  "",
  "Groups with fewer than 20 cells were explicitly marked LOW_CELL_COUNT. All comparisons are descriptive. No condition-level differential expression, p-values, significance tests, differential-abundance inference or cell-as-replicate claims were produced. Bars start at zero; DotPlot area encodes detection percentage and colour encodes mean expression using common scales across compared condition facets."
)
writeLines(methods_lines, out_path("08_methods", "METHODS_CONDITION_DESCRIPTIVE_TOP20.md"))

required_figure_stems <- c(
  "CONDITION_BROAD_CELLTYPE_STACKED_BAR", "CONDITION_BROAD_CELLTYPE_GROUPED_BAR",
  "TOP20_BROAD_CELLTYPE_CONDITION_DOTPLOT", "TOP20_BROAD_CELLTYPE_CONDITION_HEATMAP",
  "TOP20_28CLUSTER_CONDITION_DOTPLOT", "CANDIDATE_CONDITION_PCT_EXPRESSION_BAR",
  "CANDIDATE_CONDITION_MEAN_EXPRESSION_BAR",
  vapply(candidate_definitions$candidate_context, function(candidate) {
    paste0("FOCUSED_", toupper(gsub("[^A-Za-z0-9]+", "_", candidate)), "_BY_CONDITION")
  }, character(1))
)
required_figure_paths <- unlist(lapply(required_figure_stems, function(stem) {
  out_path("05_figures", paste0(stem, c(".pdf", ".png")))
}), use.names = FALSE)
valid_pdf <- function(path) identical(rawToChar(readBin(path, "raw", n = 4L)), "%PDF")
valid_png <- function(path) {
  identical(as.integer(readBin(path, "raw", n = 8L)),
            c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))
}
nonempty_broad <- broad_stats[n_total_cells > 0L]
validation <- rbindlist(list(
  data.table(check = "fixed_cell_count", observed = as.character(ncol(obj)), expected = "36111",
             pass = ncol(obj) == 36111L),
  data.table(check = "fixed_cluster_count", observed = as.character(uniqueN(meta$cluster)), expected = "28",
             pass = uniqueN(meta$cluster) == 28L),
  data.table(check = "mapping_rows_complete", observed = as.character(sum(mapping_audit$mapping_complete)),
             expected = "3", pass = nrow(mapping_audit) == 3L && all(mapping_audit$mapping_complete)),
  data.table(check = "condition_total_conservation", observed = as.character(sum(condition_totals$total_cells)),
             expected = "36111", pass = sum(condition_totals$total_cells) == 36111L),
  data.table(check = "composition_total_conservation", observed = as.character(sum(composition$cell_count)),
             expected = "36111", pass = sum(composition$cell_count) == 36111L),
  data.table(check = "composition_fraction_sums", observed = paste(round(composition[, sum(cell_type_fraction_within_condition), by = condition]$V1, 10), collapse = ","),
             expected = "1,1,1", pass = all(abs(composition[, sum(cell_type_fraction_within_condition), by = condition]$V1 - 1) < 1e-10)),
  data.table(check = "broad_annotation_count_invariants", observed = as.character(sum(broad_count_check$pass)),
             expected = "7", pass = all(broad_count_check$pass)),
  data.table(check = "top20_feature_presence", observed = as.character(sum(feature_audit$status == "FOUND")),
             expected = "20", pass = all(feature_audit$status == "FOUND")),
  data.table(check = "broad_expression_table_rows", observed = as.character(nrow(broad_stats)),
             expected = as.character(20L * length(broad_levels) * 3L),
             pass = nrow(broad_stats) == 20L * length(broad_levels) * 3L),
  data.table(check = "cluster_expression_table_rows", observed = as.character(nrow(cluster_stats)),
             expected = as.character(20L * 28L * 3L), pass = nrow(cluster_stats) == 20L * 28L * 3L),
  data.table(check = "top1_top2_rows", observed = as.character(nrow(top2_localization)), expected = "40",
             pass = nrow(top2_localization) == 40L),
  data.table(check = "candidate_context_rows", observed = as.character(nrow(candidate_stats)), expected = "12",
             pass = nrow(candidate_stats) == 12L),
  data.table(check = "expression_count_bounds", observed = as.character(sum(nonempty_broad$n_expressing <= nonempty_broad$n_total_cells)),
             expected = as.character(nrow(nonempty_broad)),
             pass = all(nonempty_broad$n_expressing <= nonempty_broad$n_total_cells)),
  data.table(check = "pct_expression_consistency", observed = "max_abs_error", expected = "<1e-10",
             pass = max(abs(nonempty_broad$pct_expressing -
                              100 * nonempty_broad$n_expressing / nonempty_broad$n_total_cells)) < 1e-10),
  data.table(check = "median_nonnegative", observed = as.character(min(nonempty_broad$median_expression)), expected = ">=0",
             pass = all(nonempty_broad$median_expression >= 0)),
  data.table(check = "low_count_flags_complete", observed = as.character(nrow(low_count_combinations)), expected = ">=3",
             pass = nrow(low_count_combinations) >= 3L && all(low_count_combinations$n_total_cells < 20L)),
  data.table(check = "figure_pair_count", observed = as.character(sum(file.exists(required_figure_paths))),
             expected = as.character(length(required_figure_paths)), pass = all(file.exists(required_figure_paths))),
  data.table(check = "figure_files_nonempty", observed = as.character(sum(file.info(required_figure_paths)$size > 0)),
             expected = as.character(length(required_figure_paths)), pass = all(file.info(required_figure_paths)$size > 0)),
  data.table(check = "pdf_signatures", observed = as.character(sum(vapply(required_figure_paths[grepl("\\.pdf$", required_figure_paths)], valid_pdf, logical(1)))),
             expected = as.character(sum(grepl("\\.pdf$", required_figure_paths))),
             pass = all(vapply(required_figure_paths[grepl("\\.pdf$", required_figure_paths)], valid_pdf, logical(1)))),
  data.table(check = "png_signatures", observed = as.character(sum(vapply(required_figure_paths[grepl("\\.png$", required_figure_paths)], valid_png, logical(1)))),
             expected = as.character(sum(grepl("\\.png$", required_figure_paths))),
             pass = all(vapply(required_figure_paths[grepl("\\.png$", required_figure_paths)], valid_png, logical(1))))
))
write_tsv(validation, out_path("09_qc", "CONDITION_DESCRIPTIVE_VALIDATION.tsv"))
if (any(!validation$pass)) stop("Condition descriptive validation failed")

stable_genes <- stability[localization_stability == "STABLE_SUPPORTED", gene]
summary_lines <- c(
  "# GSE252285 resolved-condition descriptive Top20 summary",
  "",
  "STATUS=SUCCESS",
  "FIXED_CELLS=36111",
  "FIXED_CLUSTERS_RES0.6=28",
  "CONDITION_MAPPING_RESOLVED=YES",
  "CONDITION_MAPPING=GEO_-3:Control;GEO_-1:HHL;GEO_-2:NIHL",
  paste0("CONDITION_TOTALS=Control:", condition_totals[condition == "Control", total_cells],
         ";HHL:", condition_totals[condition == "HHL", total_cells],
         ";NIHL:", condition_totals[condition == "NIHL", total_cells]),
  "TOP20_FOUND=20/20",
  paste0("STABLE_SUPPORTED_TOP1_GENES=", paste(stable_genes, collapse = ",")),
  paste0("LOW_CELL_COUNT_BROAD_COMBINATIONS=", nrow(low_broad)),
  paste0("LOW_CELL_COUNT_CLUSTER_COMBINATIONS=", nrow(low_cluster)),
  "CONDITION_LEVEL_DE=NOT_PERFORMED",
  "P_VALUES=NOT_COMPUTED",
  "DIFFERENTIAL_ABUNDANCE_INFERENCE=NOT_PERFORMED",
  "OBJECT_MODIFIED=NO",
  "",
  "Differences among Control, HHL and NIHL are descriptive cell-level summaries. Biological replication is unresolved, so no inferential condition claim is made."
)
writeLines(summary_lines, out_path("10_release", "CONDITION_DESCRIPTIVE_TOP20_SUMMARY.md"))

capture.output(sessionInfo(), file = out_path("07_logs", "sessionInfo.txt"))
manifest_roots <- vapply(c("00_manifest", "04_tables", "05_figures", "07_logs",
                           "08_methods", "09_qc", "10_release"),
                         function(section) out_path(section), character(1))
manifest_files <- sort(unlist(lapply(manifest_roots, function(root) {
  list.files(root, recursive = TRUE, full.names = TRUE)
}), use.names = FALSE))
manifest_files <- manifest_files[file.info(manifest_files)$isdir %in% FALSE]
manifest_files <- setdiff(manifest_files, out_path("00_manifest", "OUTPUT_MANIFEST.tsv"))
manifest <- data.table(
  relative_path = sub(paste0("^", project_dir, "/"), "", manifest_files),
  bytes = as.numeric(file.info(manifest_files)$size),
  md5 = unname(tools::md5sum(manifest_files))
)
write_tsv(manifest, out_path("00_manifest", "OUTPUT_MANIFEST.tsv"))
stage_message("Resolved-condition descriptive Top20 analysis complete")
