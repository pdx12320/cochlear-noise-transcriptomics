source(file.path(Sys.getenv("GSE252285_SCRIPT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "03_scripts")), "common.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
})

stage_message("Stage 05: fixed Top20 cellular localization")
obj <- readRDS(p("06_objects", "04_annotated_seurat.rds"))
DefaultAssay(obj) <- "RNA"

top20_genes <- c(
  "Dusp1", "Fos", "Nfkbia", "Arrdc2", "Apold1", "Sik1", "Per1", "Ccn1",
  "Cebpb", "Fosb", "Rasd1", "Gadd45g", "Depp1", "Serpine1", "Errfi1", "Nr4a1",
  "Sgk1", "Rgs1", "Egr1", "Jun"
)
feature_md <- as.data.table(obj[["RNA"]][[]], keep.rownames = "feature_name")

map_one <- function(gene) {
  if (gene %in% rownames(obj)) {
    return(data.table(gene = gene, mapped_feature = gene, status = "FOUND",
                      mapping_note = "exact feature name"))
  }
  exact_symbol <- feature_md[original_symbol == gene, feature_name]
  if (length(exact_symbol)) {
    return(data.table(gene = gene, mapped_feature = exact_symbol[1], status = "ALIAS_REQUIRED",
                      mapping_note = "mapped through original_symbol after make.unique"))
  }
  ci <- feature_md[tolower(original_symbol) == tolower(gene), feature_name]
  if (length(ci)) {
    return(data.table(gene = gene, mapped_feature = ci[1], status = "ALIAS_REQUIRED",
                      mapping_note = "case-insensitive original_symbol match"))
  }
  data.table(gene = gene, mapped_feature = NA_character_, status = "NOT_FOUND",
             mapping_note = "no exact or case-insensitive symbol match")
}
audit <- rbindlist(lapply(top20_genes, map_one))
audit[, fixed_rank := match(gene, top20_genes)]
setorder(audit, fixed_rank)

data_mat <- SeuratObject::LayerData(obj[["RNA"]], layer = "data")
count_mat <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")
broad_order <- c("Hair_sensory", "Supporting_epithelial", "Endothelial_vascular",
                 "Immune_myeloid", "Fibrocyte_mesenchymal", "Neuronal", "Glial",
                 "Other", "UNRESOLVED")
broad_present <- intersect(broad_order, unique(as.character(obj$final_broad_label)))

stats_rows <- list()
for (i in seq_along(top20_genes)) {
  gene_name <- top20_genes[i]
  feature <- audit[gene == gene_name, mapped_feature]
  feature_status <- audit[gene == gene_name, status]
  if (is.na(feature)) {
    stats_rows[[i]] <- data.table(
      gene = gene_name, fixed_rank = i, feature_status = feature_status,
      final_broad_label = broad_present, mean_expression = NA_real_,
      average_log_normalized_expression = NA_real_, pct_expression = NA_real_,
      cell_number = vapply(broad_present, function(lbl) sum(obj$final_broad_label == lbl), integer(1))
    )
    next
  }
  stats_rows[[i]] <- rbindlist(lapply(broad_present, function(lbl) {
    cells <- colnames(obj)[obj$final_broad_label == lbl]
    raw <- count_mat[feature, cells, drop = TRUE]
    norm <- data_mat[feature, cells, drop = TRUE]
    data.table(
      gene = gene_name, fixed_rank = i, feature_status = feature_status,
      final_broad_label = lbl,
      mean_expression = mean(raw),
      average_log_normalized_expression = mean(norm),
      pct_expression = mean(raw > 0),
      cell_number = length(cells)
    )
  }))
}
top20_stats <- rbindlist(stats_rows)
top20_stats[, specificity_score := {
  den <- sum(average_log_normalized_expression, na.rm = TRUE)
  if (den > 0) average_log_normalized_expression / den else 0
}, by = gene]
top20_stats[, broad_order__ := match(final_broad_label, broad_order)]
setorder(top20_stats, fixed_rank, broad_order__)
top20_stats[, broad_order__ := NULL]
write_tsv(top20_stats, p("04_tables", "TOP20_BY_BROAD_CELLTYPE.tsv"))

detectability <- rbindlist(lapply(seq_len(nrow(audit)), function(i) {
  feature <- audit$mapped_feature[i]
  detected <- if (is.na(feature)) 0L else sum(count_mat[feature, , drop = TRUE] > 0)
  data.table(
    gene = audit$gene[i],
    detected_cells = detected,
    overall_pct_expression = detected / ncol(obj),
    detectability = if (is.na(feature)) "NOT_FOUND" else
      if (detected < max(10, 0.01 * ncol(obj))) "LOW_DETECTABILITY" else "DETECTABLE"
  )
}))
audit <- merge(audit, detectability, by = "gene", sort = FALSE)
audit[, fixed_rank := match(gene, top20_genes)]
setorder(audit, fixed_rank)
write_tsv(audit, p("04_tables", "TOP20_FEATURE_AUDIT.tsv"))

found_audit <- audit[status != "NOT_FOUND"]
found_features <- found_audit$mapped_feature
names(found_features) <- found_audit$gene
if (length(found_features)) {
  p_dot <- DotPlot(obj, features = unname(found_features), group.by = "final_broad_label") +
    RotatedAxis() + labs(title = "Fixed Top20 by provisional broad cell type") + theme_sc()
  save_dual(p_dot, "TOP20_BROAD_CELLTYPE_DOTPLOT", 14, 7)

  heat_dt <- top20_stats[!is.na(average_log_normalized_expression)]
  heat_dt[, gene := factor(gene, levels = rev(top20_genes))]
  heat_dt[, final_broad_label := factor(final_broad_label, levels = broad_present)]
  p_heat <- ggplot(heat_dt, aes(final_broad_label, gene,
                               fill = average_log_normalized_expression)) +
    geom_tile() + scale_fill_viridis_c(option = "magma") +
    labs(x = NULL, y = NULL, fill = "Mean log-normalized") + theme_sc() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_dual(p_heat, "TOP20_MEAN_EXPRESSION_HEATMAP", 11, 9)

  p_feature <- FeaturePlot(obj, features = unname(found_features), reduction = "umap",
                           order = TRUE, min.cutoff = "q05", max.cutoff = "q95",
                           ncol = 4, raster = TRUE, combine = TRUE)
  save_dual(p_feature, "TOP20_FEATUREPLOT_ATLAS", 16, 24, dpi = 250)

  p_vln <- VlnPlot(obj, features = unname(found_features), group.by = "final_broad_label",
                   pt.size = 0, ncol = 4, combine = TRUE)
  save_dual(p_vln, "TOP20_VLNPLOT", 18, 24, dpi = 250)
}

focus_summary <- list()
for (gene in c("Apold1", "Errfi1")) {
  gene_name <- gene
  feature <- audit[gene == gene_name, mapped_feature]
  if (is.na(feature)) next
  p_fp <- FeaturePlot(obj, features = feature, reduction = "umap", order = TRUE,
                      min.cutoff = "q05", max.cutoff = "q95", raster = TRUE) +
    ggtitle(paste(gene, "expression")) + theme_sc()
  p_vp <- VlnPlot(obj, features = feature, group.by = "final_broad_label", pt.size = 0) +
    ggtitle(paste(gene, "by broad cell type")) + theme_sc() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p_dp <- DotPlot(obj, features = feature, group.by = "final_broad_label") +
    ggtitle(paste(gene, "by broad cell type")) + theme_sc()
  save_dual(p_fp, paste0("TOP20_", gene, "_FeaturePlot_HIGHRES"), 9, 8, dpi = 320)
  save_dual(p_vp, paste0("TOP20_", gene, "_VlnPlot_HIGHRES"), 11, 7, dpi = 320)
  save_dual(p_dp, paste0("TOP20_", gene, "_DotPlot_HIGHRES"), 8, 6, dpi = 320)

  gene_stats <- top20_stats[gene == gene_name]
  if (all(gene_stats$average_log_normalized_expression == 0, na.rm = TRUE)) {
    primary <- "NONE_DETECTED"
  } else {
    primary <- gene_stats[which.max(average_log_normalized_expression), final_broad_label]
  }
  focus_summary[[gene]] <- data.table(
    gene = gene,
    primary_expression_cell_type = primary,
    detectability = audit[gene == gene_name, detectability],
    basis = "maximum average log-normalized expression; descriptive only"
  )
}
write_tsv(rbindlist(focus_summary), p("04_tables", "APOLD1_ERRFI1_PRIMARY_CELLTYPE.tsv"))
stage_message("Stage 05 complete: found=", sum(audit$status != "NOT_FOUND"),
              "; not found=", sum(audit$status == "NOT_FOUND"))
