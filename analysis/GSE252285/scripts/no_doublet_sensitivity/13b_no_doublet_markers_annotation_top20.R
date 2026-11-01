source(file.path(Sys.getenv("SCRIPT_DIR", "."), "branch_common.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(patchwork)
  library(future)
})

stage_message("Stage 13b: markers, reference-supported conservative annotation, and Top20")
clustered_path <- branch_path("06_objects", "01_no_doublet_clustered_tsne.rds")
annotation_reference_path <- Sys.getenv(
  "ANNOTATION_REFERENCE_TSV",
  file.path(Sys.getenv("SCRIPT_DIR", "."), "..", "..", "inputs", "annotation_reference.tsv")
)
if (!file.exists(annotation_reference_path)) stop("Required annotation reference TSV is missing")
obj <- readRDS(clustered_path)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "RNA_snn_res.0.6"
if (ncol(obj) != 36111L) stop("Expected 36,111 no-doublet-filter cells")
if ("umap" %in% Reductions(obj)) stop("UMAP is forbidden in this branch")
if (!"tsne" %in% Reductions(obj)) stop("t-SNE reduction is missing")
if (any(grepl("scDblFinder", colnames(obj[[]]), ignore.case = TRUE))) {
  stop("No-doublet-filter object contains forbidden doublet metadata")
}
set.seed(seed)

marker_path <- branch_path("04_tables", "ALL_CLUSTER_MARKERS.tsv.gz")
if (file.exists(marker_path) && file.info(marker_path)$size > 1000) {
  stage_message("Reusing completed marker checkpoint: ", marker_path)
  markers <- fread(marker_path, showProgress = FALSE)
} else {
  stage_message("FindAllMarkers: Wilcoxon; positive markers; min.pct=0.25; log2FC threshold=0.25")
  marker_workers <- max(1L, min(4L, suppressWarnings(as.integer(Sys.getenv("SC_THREADS", "4")))))
  if (!is.finite(marker_workers)) marker_workers <- 4L
  if (.Platform$OS.type == "unix" && marker_workers > 1L) {
    future::plan(future::multicore, workers = marker_workers)
    options(future.globals.maxSize = 96 * 1024^3)
  } else {
    future::plan(future::sequential)
  }
  markers <- FindAllMarkers(
    obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
    test.use = "wilcox", return.thresh = 1, random.seed = seed,
    densify = FALSE, verbose = TRUE
  )
  future::plan(future::sequential)
  markers <- as.data.table(markers)
  if (!nrow(markers)) stop("FindAllMarkers returned no markers")
}
fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(markers))[1]
if (is.na(fc_col)) stop("No average log fold-change column found")
if (!all(c("pct.1", "pct.2") %in% names(markers))) stop("pct.1/pct.2 missing")
markers[, cluster := as.character(cluster)]
markers[, pct_diff := pct.1 - pct.2]
setorderv(markers, c("cluster", fc_col, "pct_diff", "p_val_adj"),
          c(1, -1, -1, 1), na.last = TRUE)
write_tsv(markers, marker_path)
top30 <- markers[, head(.SD, 30), by = cluster]
write_tsv(top30, branch_path("04_tables", "TOP30_MARKERS_PER_CLUSTER.tsv"))

canonical_all <- unique(unlist(canonical_panels, use.names = FALSE))
canonical_audit <- data.table(
  marker = canonical_all,
  status = ifelse(canonical_all %in% rownames(obj), "FOUND", "NOT_FOUND")
)
write_tsv(canonical_audit, branch_path("04_tables", "CANONICAL_MARKER_AUDIT.tsv"))
canonical_present <- intersect(canonical_all, rownames(obj))

data_mat <- LayerData(obj[["RNA"]], layer = "data")
cluster_ids <- sort(unique(as.character(Idents(obj))))
cluster_cells <- split(colnames(obj), as.character(Idents(obj)))
lineage_rows <- list()
for (cl in cluster_ids) {
  cells <- cluster_cells[[cl]]
  for (lineage in names(canonical_panels)) {
    genes <- intersect(canonical_panels[[lineage]], rownames(data_mat))
    if (!length(genes)) {
      lineage_rows[[length(lineage_rows) + 1L]] <- data.table(
        cluster = cl, lineage = lineage, panel_genes_present = 0L,
        mean_panel_expression = 0, mean_panel_pct = 0, module_score = 0,
        supported_marker_number = 0L, supported_markers = "NONE"
      )
      next
    }
    sm <- data_mat[genes, cells, drop = FALSE]
    avg <- Matrix::rowMeans(sm)
    pct <- Matrix::rowMeans(sm > 0)
    supported <- genes[avg >= 0.10 & pct >= 0.05]
    lineage_rows[[length(lineage_rows) + 1L]] <- data.table(
      cluster = cl,
      lineage = lineage,
      panel_genes_present = length(genes),
      mean_panel_expression = mean(avg),
      mean_panel_pct = mean(pct),
      module_score = mean(avg * pct),
      supported_marker_number = length(supported),
      supported_markers = if (length(supported)) paste(supported, collapse = ",") else "NONE"
    )
  }
}
lineage_scores <- rbindlist(lineage_rows)
write_tsv(lineage_scores, branch_path("04_tables", "LINEAGE_MODULE_SCORES.tsv"))

# The barcode-level reference preserves the prior labels used by the adopted annotation.
# No expression data or doublet scores/classes are imported from the earlier object.
reference <- fread(annotation_reference_path, colClasses = "character", showProgress = FALSE)
if (!all(c("barcode", "frozen_label") %in% names(reference))) {
  stop("Annotation reference requires barcode and frozen_label columns")
}
if (nrow(reference) != 31273L || uniqueN(reference$barcode) != 31273L ||
    anyNA(reference$barcode) || anyNA(reference$frozen_label)) {
  stop("Annotation reference must contain 31,273 unique, labeled barcodes")
}
if (run_legacy_comparison && !"frozen_cluster" %in% names(reference)) {
  stop("Optional legacy cluster comparison requires a frozen_cluster reference column")
}
reference_columns <- c("barcode", "frozen_label",
                       if (run_legacy_comparison) "frozen_cluster")
frozen_meta <- reference[, ..reference_columns]
rm(reference)
gc(verbose = FALSE)

new_meta <- data.table(
  barcode = colnames(obj),
  cluster = as.character(Idents(obj)),
  library_id = as.character(obj$library_id)
)
shared <- merge(new_meta, frozen_meta, by = "barcode", all = FALSE, sort = FALSE)
if (nrow(shared) != 31273L) stop("Expected 31,273 shared frozen cells; observed ", nrow(shared))
new_meta[, membership := ifelse(barcode %in% shared$barcode, "SHARED_WITH_FROZEN", "ADDITIONAL_ROBUST_QC_CELL")]
if (sum(new_meta$membership == "ADDITIONAL_ROBUST_QC_CELL") != 4838L) {
  stop("Expected 4,838 additional robust-QC cells")
}

prior_counts <- shared[, .N, by = .(cluster, frozen_label)]
prior_counts[, shared_cells := sum(N), by = cluster]
prior_counts[, frozen_label_fraction := N / shared_cells]
setorder(prior_counts, cluster, -N, frozen_label)
prior_top <- prior_counts[, .SD[1], by = cluster]
setnames(prior_top, "N", "majority_label_cells")
cluster_total <- new_meta[, .(cell_number = .N), by = cluster]
cluster_added <- new_meta[membership == "ADDITIONAL_ROBUST_QC_CELL",
                          .(additional_cells = .N), by = cluster]
prior_top <- merge(prior_top, cluster_total, by = "cluster", all = TRUE)
prior_top <- merge(prior_top, cluster_added, by = "cluster", all = TRUE)
prior_top[is.na(additional_cells), additional_cells := 0L]
prior_top[, additional_fraction := additional_cells / cell_number]
if (run_legacy_comparison) {
  write_tsv(prior_counts, branch_path("04_tables", "FROZEN_LABEL_COMPOSITION_BY_NEW_CLUSTER.tsv"))
}

top_marker_text <- top30[, .(
  top30_markers = paste(gene, collapse = ","),
  top10_specificity_markers = paste(head(gene[order(-pct_diff, -get(fc_col))], 10), collapse = ",")
), by = cluster]

annotation_rows <- lapply(cluster_ids, function(cl) {
  prior <- prior_top[cluster == cl]
  ev <- copy(lineage_scores[cluster == cl])
  setorder(ev, -module_score, -supported_marker_number, lineage)
  eligible <- ev[supported_marker_number >= 2]
  top_lineage <- if (nrow(eligible)) eligible$lineage[1] else "NONE"
  top_score <- if (nrow(eligible)) eligible$module_score[1] else 0
  top_support <- if (nrow(eligible)) eligible$supported_marker_number[1] else 0L
  strong_conflicts <- if (nrow(eligible) > 1L && top_score > 0) {
    eligible[lineage != top_lineage & module_score >= 0.75 * top_score]
  } else eligible[0]

  prior_label <- prior$frozen_label
  prior_ev <- ev[lineage == prior_label]
  prior_support <- if (nrow(prior_ev)) prior_ev$supported_marker_number[1] else 0L
  prior_markers <- if (nrow(prior_ev)) prior_ev$supported_markers[1] else "NONE"
  prior_score <- if (nrow(prior_ev)) prior_ev$module_score[1] else 0
  purity <- prior$frozen_label_fraction
  shared_n <- prior$shared_cells

  final_label <- "UNRESOLVED"
  confidence <- "UNRESOLVED"
  reason <- character()
  if (!length(shared_n) || is.na(shared_n) || shared_n < 20L) {
    reason <- c(reason, "fewer_than_20_shared_frozen_cells")
  } else if (!length(purity) || is.na(purity) || purity < 0.70) {
    reason <- c(reason, "frozen_label_purity_below_0.70")
  } else if (prior_label == "Other") {
    if (top_support >= 3L) {
      reason <- c(reason, paste0("canonical_lineage_signal_conflicts_with_Other:", top_lineage))
    } else {
      final_label <- "Other"
      confidence <- if (purity >= 0.90) "HIGH" else "MODERATE"
      reason <- c(reason, "high_frozen_Other_concordance_without_dominant_canonical_panel")
    }
  } else if (!prior_label %in% names(canonical_panels)) {
    reason <- c(reason, "historical_label_not_in_allowed_broad_panel")
  } else if (prior_support < 2L) {
    reason <- c(reason, "fewer_than_two_supporting_markers_for_historical_lineage")
  } else if (nrow(strong_conflicts) && any(strong_conflicts$lineage != prior_label)) {
    reason <- c(reason, paste0("strong_multi_lineage_conflict:",
                              paste(strong_conflicts$lineage, collapse = ",")))
  } else if (top_lineage != prior_label && top_score >= 1.5 * max(prior_score, 1e-12)) {
    reason <- c(reason, paste0("canonical_top_lineage_disagrees:", top_lineage))
  } else {
    final_label <- prior_label
    confidence <- if (purity >= 0.90 && prior_support >= 4L && top_lineage == prior_label) {
      "HIGH"
    } else if (purity >= 0.75 && prior_support >= 3L) {
      "MODERATE"
    } else {
      "LOW"
    }
    reason <- c(reason, "historical_v3_majority_and_multi_gene_canonical_support_agree")
  }

  conflict_text <- if (nrow(strong_conflicts)) {
    paste(paste0(strong_conflicts$lineage, "[", strong_conflicts$supported_markers, "]"),
          collapse = ";")
  } else "NONE"
  data.table(
    cluster = cl,
    cell_number = prior$cell_number,
    shared_frozen_cells = shared_n,
    additional_robust_qc_cells = prior$additional_cells,
    additional_fraction = prior$additional_fraction,
    frozen_v3_majority_label = prior_label,
    frozen_v3_majority_fraction = purity,
    canonical_top_lineage = top_lineage,
    canonical_top_supported_markers = if (nrow(eligible)) eligible$supported_markers[1] else "NONE",
    supporting_markers = prior_markers,
    conflicting_markers = conflict_text,
    broad_annotation_no_doublet_sensitivity = final_label,
    annotation_confidence_no_doublet_sensitivity = confidence,
    reason = paste(reason, collapse = ";"),
    warning = if (prior$additional_fraction >= 0.30 || final_label == "UNRESOLVED")
      "POTENTIAL_MIXED_PROFILE_WITHOUT_DOUBLET_FILTER" else "NONE"
  )
})
annotation <- rbindlist(annotation_rows)
annotation <- merge(annotation, top_marker_text, by = "cluster", all.x = TRUE)
annotation[, cluster_numeric__ := suppressWarnings(as.integer(cluster))]
setorder(annotation, cluster_numeric__)
annotation[, cluster_numeric__ := NULL]
write_tsv(annotation, branch_path("04_tables", "BROAD_ANNOTATION_SENSITIVITY.tsv"))

obj$broad_annotation_no_doublet_sensitivity <- annotation$broad_annotation_no_doublet_sensitivity[
  match(as.character(Idents(obj)), annotation$cluster)]
obj$annotation_confidence_no_doublet_sensitivity <- annotation$annotation_confidence_no_doublet_sensitivity[
  match(as.character(Idents(obj)), annotation$cluster)]
obj$annotation_version_no_doublet_sensitivity <- "marker_plus_frozen_concordance_no_doublet_sensitivity_v1"
if (anyNA(obj$broad_annotation_no_doublet_sensitivity)) stop("Annotation mapping created NA cells")

annotation_counts <- as.data.table(table(obj$broad_annotation_no_doublet_sensitivity))
setnames(annotation_counts, c("broad_cell_type", "cell_number"))
annotation_counts[, broad_cell_type := as.character(broad_cell_type)]
write_tsv(annotation_counts, branch_path("04_tables", "BROAD_CELLTYPE_COUNTS.tsv"))

cluster_lib <- fread(branch_path("04_tables", "CLUSTER_BY_LIBRARY_COMPOSITION.tsv"))
cluster_lib[, cluster := as.character(cluster)]
cluster_lib_wide <- dcast(cluster_lib, cluster ~ library_id,
                          value.var = "fraction_within_cluster", fill = 0)
annotation_all <- merge(annotation, cluster_lib_wide, by = "cluster", all.x = TRUE)
setcolorder(annotation_all, c("cluster", "cell_number", "LIB1", "LIB2", "LIB3",
                             setdiff(names(annotation_all), c("cluster", "cell_number", "LIB1", "LIB2", "LIB3"))))
write_tsv(annotation_all, branch_path("04_tables", "ANNOTATION_AUDIT_ALL_CLUSTERS.tsv"))

p_ann <- DimPlot(obj, reduction = "tsne", group.by = "broad_annotation_no_doublet_sensitivity",
                 label = TRUE, repel = TRUE, raster = TRUE) +
  scale_colour_manual(values = broad_palette, drop = FALSE) +
  labs(title = "Conservative broad annotation: no-doublet-filter sensitivity",
       subtitle = "Marker-supported sensitivity labels; not an annotation freeze") +
  theme_branch()
save_branch_dual(p_ann, "TSNE_BY_BROAD_ANNOTATION_SENSITIVITY", 11, 8)

p_canonical <- DotPlot(obj, features = canonical_present,
                       group.by = "seurat_clusters_no_doublet") +
  RotatedAxis() +
  labs(title = "Canonical lineage markers by no-doublet-filter cluster",
       subtitle = "Dot size: percent expressing; colour: scaled average expression") +
  theme_branch()
save_branch_dual(p_canonical, "CANONICAL_MARKER_DOTPLOT_BY_CLUSTER", 18, 9)

score_plot <- copy(lineage_scores)
score_plot[, module_z := safe_zscore(module_score), by = lineage]
score_plot[, cluster := factor(cluster, levels = rev(cluster_ids))]
p_score <- ggplot(score_plot, aes(lineage, cluster, fill = module_z)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(x = NULL, y = "Cluster", fill = "lineage-wise z",
       title = "Canonical lineage module scores") +
  theme_branch() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_branch_dual(p_score, "LINEAGE_SCORE_HEATMAP", 10, 9)

if (run_legacy_comparison) {
  added_plot <- annotation[, .(cluster, additional_fraction, warning)]
  p_added <- ggplot(added_plot, aes(reorder(cluster, -additional_fraction), additional_fraction,
                                    fill = warning)) +
    geom_col() +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = c(NONE = "#0072B2",
                                 POTENTIAL_MIXED_PROFILE_WITHOUT_DOUBLET_FILTER = "#D55E00")) +
    labs(x = "New cluster", y = "Cells absent from frozen 31,273-cell atlas",
         fill = NULL, title = "Additional robust-QC cell fraction by new cluster") +
    theme_branch()
  save_branch_dual(p_added, "ADDITIONAL_CELL_FRACTION_BY_CLUSTER", 12, 6)
}

stage_message("Fixed Top20 localization by no-doublet-filter sensitivity annotation")
feature_lut <- setNames(rownames(obj), toupper(rownames(obj)))
matched <- unname(feature_lut[toupper(top20_fixed)])
feature_audit <- data.table(
  fixed_rank = seq_along(top20_fixed),
  gene = top20_fixed,
  matched_feature = matched,
  status = ifelse(is.na(matched), "NOT_FOUND", ifelse(matched == top20_fixed, "FOUND", "ALIAS_REQUIRED"))
)
write_tsv(feature_audit, branch_path("04_tables", "TOP20_FEATURE_AUDIT.tsv"))
found_audit <- feature_audit[!is.na(matched_feature)]
if (nrow(found_audit) != 20L) stop("Expected 20/20 fixed Top20 genes; found ", nrow(found_audit))

counts_mat <- LayerData(obj[["RNA"]], layer = "counts")
label_field <- "broad_annotation_no_doublet_sensitivity"
broad_order <- c("Hair_sensory", "Supporting_epithelial", "Endothelial_vascular",
                 "Immune_myeloid", "Fibrocyte_mesenchymal", "Neuronal", "Glial",
                 "Other", "UNRESOLVED")
broad_present <- intersect(broad_order, unique(as.character(obj[[label_field, drop = TRUE]])))
top20_rows <- list()
for (i in seq_len(nrow(found_audit))) {
  gene <- found_audit$gene[i]
  feature <- found_audit$matched_feature[i]
  for (label in broad_present) {
    cells <- colnames(obj)[obj[[label_field, drop = TRUE]] == label]
    raw_vals <- counts_mat[feature, cells, drop = TRUE]
    log_vals <- data_mat[feature, cells, drop = TRUE]
    top20_rows[[length(top20_rows) + 1L]] <- data.table(
      fixed_rank = found_audit$fixed_rank[i],
      gene = gene,
      broad_cell_type = label,
      cell_number = length(cells),
      mean_expression = mean(raw_vals),
      average_log_normalized_expression = mean(log_vals),
      pct_expression = 100 * mean(raw_vals > 0)
    )
  }
}
top20_stats <- rbindlist(top20_rows)
top20_stats[, specificity_descriptive := {
  denom <- sum(average_log_normalized_expression)
  if (denom > 0) average_log_normalized_expression / denom else 0
}, by = gene]
top20_stats[, detectability := ifelse(pct_expression < 1, "LOW_DETECTABILITY", "DETECTABLE")]
top20_stats[, broad_order__ := match(broad_cell_type, broad_order)]
setorder(top20_stats, fixed_rank, broad_order__)
top20_stats[, broad_order__ := NULL]
write_tsv(top20_stats, branch_path("04_tables", "TOP20_BY_BROAD_CELLTYPE.tsv"))

candidate_rows <- lapply(c("Apold1", "Errfi1"), function(candidate_gene) {
  dt <- copy(top20_stats[gene == candidate_gene])
  setorder(dt, -mean_expression, -pct_expression)
  feature <- found_audit[gene == candidate_gene, matched_feature][1]
  vals <- as.numeric(counts_mat[feature, , drop = TRUE])
  positive <- sort(vals[vals > 0], decreasing = TRUE)
  top_n <- max(1L, ceiling(0.01 * length(vals)))
  top_share <- if (sum(positive) > 0) sum(head(positive, top_n)) / sum(positive) else NA_real_
  data.table(
    gene = candidate_gene,
    primary_cell_type = dt$broad_cell_type[1],
    primary_mean_expression = dt$mean_expression[1],
    primary_pct_expression = dt$pct_expression[1],
    second_cell_type = if (nrow(dt) >= 2) dt$broad_cell_type[2] else "NA",
    second_mean_expression = if (nrow(dt) >= 2) dt$mean_expression[2] else NA_real_,
    second_pct_expression = if (nrow(dt) >= 2) dt$pct_expression[2] else NA_real_,
    number_celltypes_pct_ge_10 = sum(dt$pct_expression >= 10),
    top_1pct_raw_count_share = top_share,
    few_cell_driver = ifelse(dt$pct_expression[1] < 1 || (!is.na(top_share) && top_share >= 0.5),
                             "YES", "NO"),
    widespread_across_broad_types = ifelse(sum(dt$pct_expression >= 10) >= 3, "YES", "NO")
  )
})
candidate_context <- rbindlist(candidate_rows)
write_tsv(candidate_context, branch_path("04_tables", "APOLD1_ERRFI1_CONTEXT.tsv"))

p_topdot <- DotPlot(obj, features = found_audit$matched_feature, group.by = label_field) +
  RotatedAxis() +
  labs(title = "Fixed Top20 by broad cell type",
       subtitle = "No-doublet-filter sensitivity labels; fixed gene order retained") +
  theme_branch()
save_branch_dual(p_topdot, "TOP20_DOTPLOT", 15, 7)

heat_dt <- copy(top20_stats)
heat_dt[, row_z := safe_zscore(average_log_normalized_expression), by = gene]
heat_dt[, gene := factor(gene, levels = rev(top20_fixed))]
heat_dt[, broad_cell_type := factor(broad_cell_type, levels = broad_present)]
p_topheat <- ggplot(heat_dt, aes(broad_cell_type, gene, fill = row_z)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "gene-wise z",
       title = "Fixed Top20 mean log-normalized expression") +
  theme_branch() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_branch_dual(p_topheat, "TOP20_MEAN_EXPRESSION_HEATMAP", 11, 9)

p_atlas <- FeaturePlot(obj, reduction = "tsne", features = found_audit$matched_feature,
                       ncol = 4, order = TRUE, raster = TRUE, keep.scale = "feature")
save_branch_dual(p_atlas, "TOP20_FEATUREPLOT_ATLAS", 18, 22, dpi = 250)
p_vln <- VlnPlot(obj, features = found_audit$matched_feature, group.by = label_field,
                 pt.size = 0, ncol = 4) + plot_layout(guides = "collect")
save_branch_dual(p_vln, "TOP20_VLNPLOT", 18, 22, dpi = 250)

for (candidate_gene in c("Apold1", "Errfi1")) {
  feature <- found_audit[gene == candidate_gene, matched_feature][1]
  p_fp <- FeaturePlot(obj, reduction = "tsne", features = feature, order = TRUE,
                      raster = TRUE) +
    labs(title = paste0(candidate_gene, " on t-SNE"), subtitle = "No-doublet-filter sensitivity") +
    theme_branch()
  p_vp <- VlnPlot(obj, features = feature, group.by = label_field, pt.size = 0) +
    labs(title = paste0(candidate_gene, " by broad cell type"), x = NULL) + theme_branch() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p_dp <- DotPlot(obj, features = feature, group.by = label_field) +
    labs(title = paste0(candidate_gene, " by broad cell type"), x = NULL) + theme_branch()
  save_branch_dual(p_fp, paste0(toupper(candidate_gene), "_FEATUREPLOT"), 9, 7)
  save_branch_dual(p_vp, paste0(toupper(candidate_gene), "_VLNPLOT"), 10, 7)
  save_branch_dual(p_dp, paste0(toupper(candidate_gene), "_DOTPLOT"), 8, 5)
}

composition <- as.data.table(obj[[]])[, .(cell_count = .N),
  by = .(library_id, broad_cell_type = get(label_field))]
composition[, fraction_within_library := cell_count / sum(cell_count), by = library_id]
setorder(composition, library_id, broad_cell_type)
write_tsv(composition, branch_path("04_tables", "CELLTYPE_COMPOSITION_DESCRIPTIVE.tsv"))
p_composition <- ggplot(composition,
                        aes(library_id, fraction_within_library, fill = broad_cell_type)) +
  geom_col() +
  scale_fill_manual(values = broad_palette, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Fraction within library", fill = "Broad cell type",
       title = "Descriptive composition only",
       subtitle = "No condition-level test; condition mapping unresolved") +
  theme_branch()
save_branch_dual(p_composition, "CELLTYPE_COMPOSITION_STACKED_BAR", 10, 7)

if (run_legacy_comparison) {
  ari <- adjusted_rand_index(shared$cluster, shared$frozen_cluster)
  cell_label_map <- annotation[, .(cluster, new_label = broad_annotation_no_doublet_sensitivity)]
  shared <- merge(shared, cell_label_map, by = "cluster", all.x = TRUE)
  label_concordance <- mean(shared$new_label == shared$frozen_label)
  comparison_summary <- data.table(
    metric = c("frozen_cells", "new_robust_qc_cells", "additional_cells_vs_frozen",
               "frozen_clusters_res0.6", "new_clusters_res0.6",
               "ARI_new_vs_frozen_clusters_on_shared_cells",
               "cell_label_concordance_on_shared_cells", "UMAP_run",
               "doublet_filter_applied", "annotation_freeze_complete_for_this_branch"),
    value = c("31273", "36111", "4838", "29", as.character(uniqueN(obj$seurat_clusters_no_doublet)),
              sprintf("%.6f", ari), sprintf("%.6f", label_concordance), "NO", "NO", "NO")
  )
  write_tsv(comparison_summary, branch_path("04_tables", "FROZEN_COMPARISON_SUMMARY.tsv"))

  cluster_overlap <- shared[, .N, by = .(new_cluster = cluster, frozen_cluster)]
  cluster_overlap[, fraction_within_new_cluster_shared := N / sum(N), by = new_cluster]
  write_tsv(cluster_overlap, branch_path("04_tables", "NEW_VS_FROZEN_CLUSTER_OVERLAP.tsv"))
  added_distribution <- merge(
    new_meta[membership == "ADDITIONAL_ROBUST_QC_CELL", .N, by = cluster],
    annotation[, .(cluster, broad_annotation_no_doublet_sensitivity,
                   annotation_confidence_no_doublet_sensitivity)],
    by = "cluster", all.x = TRUE
  )
  added_distribution[, fraction_of_additional_cells := N / sum(N)]
  write_tsv(added_distribution, branch_path("04_tables", "ADDITIONAL_CELL_DISTRIBUTION.tsv"))
}

final_object_path <- branch_path("06_objects", "02_no_doublet_annotated_tsne.rds")
saveRDS(obj, final_object_path, compress = FALSE)

capture.output(sessionInfo(), file = branch_path("07_logs", "sessionInfo.txt"))
version_lines <- c(
  paste("R", R.version.string),
  paste("Seurat", as.character(packageVersion("Seurat"))),
  paste("SeuratObject", as.character(packageVersion("SeuratObject"))),
  paste("Matrix", as.character(packageVersion("Matrix"))),
  paste("data.table", as.character(packageVersion("data.table"))),
  paste("ggplot2", as.character(packageVersion("ggplot2"))),
  paste("Rtsne", as.character(packageVersion("Rtsne")))
)
writeLines(version_lines, branch_path("07_logs", "software_versions.txt"))

cluster_n <- uniqueN(obj$seurat_clusters_no_doublet)
unresolved_clusters <- annotation[
  broad_annotation_no_doublet_sensitivity == "UNRESOLVED", cluster]
count_text <- paste0(annotation_counts$broad_cell_type, "=", annotation_counts$cell_number,
                     collapse = "; ")
summary_lines <- c(
  "# GSE252285 robust-QC / no-doublet-filter / t-SNE sensitivity",
  "",
  "- STATUS: SUCCESS_COMPUTATION_PENDING_FINAL_VALIDATION",
  "- ANALYSIS_STATUS: NO_DOUBLETFILTER_SENSITIVITY_NOT_FROZEN",
  "- INPUT_DIMENSION: 32285 genes x 39092 cells",
  "- ROBUST_QC_RETAINED: 36111",
  "- DOUBLETS_REMOVED: 0 (filter not run or used in this branch)",
  "- FINAL_CELLS: 36111",
  paste0("- NUMBER_OF_CLUSTERS_RES0.6: ", cluster_n),
  "- DISPLAY_EMBEDDING: t-SNE only; UMAP not run",
  paste0("- BROAD_CELLTYPE_COUNTS: ", count_text),
  paste0("- UNRESOLVED_CLUSTERS: ", if (length(unresolved_clusters))
         paste(unresolved_clusters, collapse = ",") else "NONE"),
  "- TOP20_FOUND: 20/20",
  paste0("- APOLD1_CONTEXT: ", candidate_context[gene == "Apold1", primary_cell_type],
         " (", sprintf("%.2f", candidate_context[gene == "Apold1", primary_pct_expression]),
         "% expressing); second=", candidate_context[gene == "Apold1", second_cell_type]),
  paste0("- ERRFI1_CONTEXT: ", candidate_context[gene == "Errfi1", primary_cell_type],
         " (", sprintf("%.2f", candidate_context[gene == "Errfi1", primary_pct_expression]),
         "% expressing); second=", candidate_context[gene == "Errfi1", second_cell_type]),
  if (run_legacy_comparison) c(
    paste0("- FROZEN_CLUSTER_ARI_ON_SHARED_CELLS: ", sprintf("%.6f", ari)),
    paste0("- FROZEN_LABEL_CONCORDANCE_ON_SHARED_CELLS: ", sprintf("%.2f%%", 100 * label_concordance))
  ),
  "- CONDITION_MAPPING_STATUS: CONDITION_MAPPING_UNRESOLVED",
  "- CONDITION_LEVEL_DE: NOT_PERFORMED",
  "- ANNOTATION_FREEZE_COMPLETE_FOR_THIS_BRANCH: NO",
  "- ANNOTATION_REFERENCE: required barcode-level prior labels; original source object unchanged",
  "",
  "## Interpretation warning",
  "",
  "The 4,838 additional robust-QC cells were retained without a doublet filter. The adopted annotation uses reference-label purity and canonical lineage evidence. A high additional-cell fraction contributes a warning; it does not independently change the label. t-SNE geometry was not used as annotation evidence.",
  "",
  "Reference labels enter the broad-annotation decision alongside canonical marker support. No scDblFinder score/class is imported for inclusion, clustering, markers, or label assignment."
)
writeLines(summary_lines, branch_path("10_release", "MODULE_SUMMARY.md"))

methods_lines <- c(
  "# Methods: robust-QC / no-doublet-filter / t-SNE sensitivity",
  "",
  "This is an explicitly non-frozen sensitivity branch. The immutable raw Seurat object and the existing robust-QC metadata were read without modification. The exact 36,111 cells passing the existing per-library robust rule were retained (LIB1 16,631; LIB2 9,207; LIB3 10,273). No doublet detection, doublet score, or doublet class was used in this branch.",
  "",
  "Counts were log-normalized with a scale factor of 10,000. Two thousand variable features were selected by the VST method, scaled, and used for 50-component PCA. PCs 1-30 were pre-specified for k-nearest-neighbor graph construction (k=20). Louvain clustering (algorithm 1) used a pre-specified primary resolution of 0.6; resolutions 0.4, 0.6, 0.8 and 1.0 were audited without selecting a resolution to reproduce any target cluster count. Random seed 42 was used.",
  "",
  "A t-SNE embedding was calculated from PCs 1-30 with perplexity 30 and seed 42. No UMAP or integration was run. t-SNE was used only to display the PCA-derived structure and not to define clusters, quantify distances, or support biological labels.",
  "",
  "Positive cluster markers were recomputed with Seurat FindAllMarkers using the Wilcoxon test, min.pct=0.25 and log2 fold-change threshold=0.25. Canonical lineage support required at least two genes with cluster-average log-normalized expression >=0.10 and expression in >=5% of cluster cells. A reference majority label among the 31,273 shared cells was used as an annotation prior. Clusters with <70% frozen-label purity, fewer than two lineage-supporting markers, or strong competing lineage modules were labeled UNRESOLVED. The resulting labels are sensitivity annotations and do not replace the frozen v3 annotation.",
  "",
  "The fixed Top20 gene order was retained exactly. Raw-count mean, mean log-normalized expression, percent expressing, cell count, descriptive specificity, and detectability were computed by sensitivity broad label. Library composition was descriptive; no condition-level differential expression or inferential test was performed because condition mapping remains unresolved."
)
writeLines(methods_lines, branch_path("08_methods", "METHODS_NO_DOUBLETFILTER_TSNE.md"))

stage_message("Stage 13b complete: markers=", nrow(markers),
              "; clusters=", cluster_n,
              "; unresolved clusters=", if (length(unresolved_clusters))
              paste(unresolved_clusters, collapse = ",") else "NONE",
              "; final object=", final_object_path)
