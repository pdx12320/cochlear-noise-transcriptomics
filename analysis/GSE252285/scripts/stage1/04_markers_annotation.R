source(file.path(Sys.getenv("GSE252285_SCRIPT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "03_scripts")), "common.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(future)
})

stage_message("Stage 04: cluster markers and provisional broad annotation")
obj <- readRDS(p("06_objects", "03_clustered_seurat.rds"))
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "RNA_snn_res.0.6"
set.seed(seed)

marker_workers <- max(1L, min(4L, suppressWarnings(as.integer(Sys.getenv("SC_THREADS", "4")))))
if (!is.finite(marker_workers)) marker_workers <- 4L
if (.Platform$OS.type == "unix" && marker_workers > 1L) {
  future::plan(future::multicore, workers = marker_workers)
  options(future.globals.maxSize = 64 * 1024^3)
} else {
  future::plan(future::sequential)
}

stage_message("Running FindAllMarkers (Wilcoxon; complete rows passing the declared pct/FC filters)")
markers <- FindAllMarkers(
  object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
  test.use = "wilcox", return.thresh = 1, random.seed = seed,
  densify = FALSE, verbose = TRUE
)
future::plan(future::sequential)
markers <- as.data.table(markers)
if (!nrow(markers)) stop("FindAllMarkers returned no markers")
write_tsv(markers, p("04_tables", "ALL_CLUSTER_MARKERS.tsv.gz"))

fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(markers))[1]
if (is.na(fc_col)) stop("Cannot identify average log-fold-change column")
setorderv(markers, c("cluster", "p_val_adj", fc_col), c(1, 1, -1), na.last = TRUE)
top20 <- markers[, head(.SD, 20), by = cluster]
top50 <- markers[, head(.SD, 50), by = cluster]
write_tsv(top20, p("04_tables", "TOP20_MARKERS_PER_CLUSTER.tsv"))
write_tsv(top50, p("04_tables", "TOP50_MARKERS_PER_CLUSTER.tsv"))

canonical <- list(
  Hair_sensory = c("Pou4f3", "Myo7a", "Otof", "Pvalb", "Slc26a5", "Ocm"),
  Supporting_epithelial = c("Sox2", "Gjb2", "Gjb6", "Krt8", "Krt18", "Epcam"),
  Endothelial_vascular = c("Pecam1", "Cdh5", "Kdr", "Emcn", "Vwf", "Esam"),
  Immune_myeloid = c("Ptprc", "Lyz2", "Tyrobp", "Aif1", "C1qa", "C1qb", "C1qc", "Cd68"),
  Fibrocyte_mesenchymal = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra"),
  Neuronal = c("Tubb3", "Snap25", "Syt1", "Rbfox3", "Nefh", "Prph"),
  Glial = c("Plp1", "Mpz", "S100b", "Sox10", "Fabp7")
)
canonical_all <- unique(unlist(canonical, use.names = FALSE))
canonical_present <- intersect(canonical_all, rownames(obj))
canonical_missing <- setdiff(canonical_all, rownames(obj))
write_tsv(data.table(
  marker = canonical_all,
  status = ifelse(canonical_all %in% canonical_present, "FOUND", "NOT_FOUND")
), p("04_tables", "CANONICAL_MARKER_AUDIT.tsv"))

data_mat <- SeuratObject::LayerData(obj[["RNA"]], layer = "data")
cluster_ids <- levels(Idents(obj))
cluster_cells <- split(colnames(obj), as.character(Idents(obj)))

evidence_rows <- list()
for (cl in cluster_ids) {
  cells <- cluster_cells[[cl]]
  for (lineage in names(canonical)) {
    genes <- intersect(canonical[[lineage]], rownames(data_mat))
    if (!length(genes)) {
      evidence_rows[[length(evidence_rows) + 1L]] <- data.table(
        cluster = cl, lineage = lineage, mean_panel_score = 0,
        supported_marker_number = 0L, supported_markers = "NONE"
      )
      next
    }
    sm <- data_mat[genes, cells, drop = FALSE]
    avg <- Matrix::rowMeans(sm)
    pct <- Matrix::rowMeans(sm > 0)
    supported <- names(avg)[avg >= 0.1 & pct >= 0.05]
    evidence_rows[[length(evidence_rows) + 1L]] <- data.table(
      cluster = cl,
      lineage = lineage,
      mean_panel_score = mean(avg * pct),
      supported_marker_number = length(supported),
      supported_markers = ifelse(length(supported), paste(supported, collapse = ","), "NONE")
    )
  }
}
evidence <- rbindlist(evidence_rows)
write_tsv(evidence, p("04_tables", "BROAD_ANNOTATION_EVIDENCE.tsv"))

cluster_sizes <- as.data.table(table(as.character(Idents(obj))))
setnames(cluster_sizes, c("cluster", "cell_number"))
cluster_sizes[, cluster := as.character(cluster)]
cluster_lib <- fread(p("04_tables", "CLUSTER_BY_LIBRARY_COMPOSITION.tsv"))

annotation_rows <- lapply(cluster_ids, function(cl) {
  ev <- copy(evidence[cluster == cl])
  setorder(ev, -mean_panel_score, -supported_marker_number)
  eligible <- ev[supported_marker_number >= 2]
  top_lineage <- if (nrow(eligible)) eligible$lineage[1] else NA_character_
  top_score <- if (nrow(eligible)) eligible$mean_panel_score[1] else 0
  conflicts <- if (nrow(eligible) > 1L) {
    eligible[lineage != top_lineage & mean_panel_score >= 0.60 * top_score]
  } else eligible[0]
  conflict_flag <- nrow(conflicts) > 0
  top_support_n <- if (nrow(eligible)) eligible$supported_marker_number[1] else 0L
  if (is.na(top_lineage) || conflict_flag) {
    final_label <- "UNRESOLVED"
    confidence <- "UNRESOLVED"
  } else {
    final_label <- top_lineage
    confidence <- if (top_support_n >= 4) "HIGH" else if (top_support_n >= 3) "MODERATE" else "LOW"
  }
  top_gene_text <- markers[as.character(cluster) == cl, paste(head(gene, 10), collapse = ",")]
  if (!length(top_gene_text) || !nzchar(top_gene_text)) top_gene_text <- "NONE"
  lib_text <- cluster_lib[cluster == cl,
    paste0(library_id, "=", cell_number, "(", sprintf("%.1f%%", 100 * fraction_within_cluster), ")",
           collapse = ";")]
  conflict_text <- if (nrow(conflicts)) {
    paste(paste0(conflicts$lineage, "[", conflicts$supported_markers, "]"), collapse = ";")
  } else "NONE"
  comments <- c(
    "PROVISIONAL_BROAD_ANNOTATION",
    paste0("library_distribution:", lib_text),
    if (conflict_flag) "POSSIBLE_MIXED_LINEAGE_SIGNAL_AFTER_DOUBLET_REMOVAL" else "no_strong_panel_conflict"
  )
  data.table(
    cluster = cl,
    cell_number = cluster_sizes[cluster == cl, cell_number],
    top_markers = top_gene_text,
    canonical_marker_support = if (nrow(eligible))
      paste0(eligible$lineage[1], "[", eligible$supported_markers[1], "]") else "NONE",
    conflicting_markers = conflict_text,
    final_broad_label = final_label,
    confidence = confidence,
    comments = paste(comments, collapse = "; ")
  )
})
annotation <- rbindlist(annotation_rows)
annotation[, cluster_numeric__ := suppressWarnings(as.integer(cluster))]
setorder(annotation, cluster_numeric__)
annotation[, cluster_numeric__ := NULL]
write_tsv(annotation, p("04_tables", "BROAD_ANNOTATION.tsv"))

obj$final_broad_label <- annotation$final_broad_label[
  match(as.character(Idents(obj)), annotation$cluster)]
obj$annotation_confidence <- annotation$confidence[
  match(as.character(Idents(obj)), annotation$cluster)]

# Marker heatmap uses a fixed display-only downsample and never alters the analysis object.
heat_features <- unique(markers[, head(gene, 5), by = cluster]$V1)
heat_features <- intersect(heat_features, rownames(obj))
if (length(heat_features)) {
  obj <- ScaleData(obj, features = heat_features, verbose = FALSE)
  set.seed(seed)
  heat_cells <- WhichCells(obj, downsample = 100, seed = seed)
  p_heat <- DoHeatmap(obj, features = heat_features, cells = heat_cells,
                      group.by = "seurat_clusters", raster = TRUE) + NoLegend()
  save_dual(p_heat, "TOP_MARKER_HEATMAP", 14, 12, dpi = 250)
}

selected_features <- canonical_present
p_selected <- DotPlot(obj, features = selected_features, group.by = "seurat_clusters") +
  RotatedAxis() + labs(title = "Selected canonical markers by cluster") + theme_sc()
save_dual(p_selected, "SELECTED_MARKER_DOTPLOT_BY_CLUSTER", 16, 8)

p_ann <- DimPlot(obj, reduction = "umap", group.by = "final_broad_label",
                 label = TRUE, repel = TRUE, raster = TRUE) + theme_sc()
p_ann_lib <- DimPlot(obj, reduction = "umap", group.by = "final_broad_label",
                     split.by = "library_id", label = FALSE, raster = TRUE,
                     ncol = length(unique(obj$library_id))) + theme_sc()
save_dual(p_ann, "UMAP_BY_BROAD_ANNOTATION", 10, 8)
save_dual(p_ann_lib, "UMAP_BROAD_ANNOTATION_SPLIT_LIBRARY", 18, 7)

p_cdot <- DotPlot(obj, features = canonical_present, group.by = "final_broad_label") +
  RotatedAxis() + labs(title = "Canonical markers by provisional broad annotation") + theme_sc()
save_dual(p_cdot, "CANONICAL_MARKER_DOTPLOT", 17, 8)

broad_levels <- intersect(c(names(canonical), "Other", "UNRESOLVED"), unique(obj$final_broad_label))
canon_avg_rows <- lapply(broad_levels, function(lbl) {
  cells <- colnames(obj)[obj$final_broad_label == lbl]
  vals <- Matrix::rowMeans(data_mat[canonical_present, cells, drop = FALSE])
  data.table(marker = names(vals), final_broad_label = lbl, mean_expression = as.numeric(vals))
})
canon_avg <- rbindlist(canon_avg_rows)
canon_avg[, z_expression := {
  z <- as.numeric(scale(mean_expression))
  ifelse(is.finite(z), z, 0)
}, by = marker]
canon_avg[, marker := factor(marker, levels = rev(canonical_present))]
canon_avg[, final_broad_label := factor(final_broad_label, levels = broad_levels)]
p_cheat <- ggplot(canon_avg, aes(final_broad_label, marker, fill = z_expression)) +
  geom_tile() + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                     midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "row z-score") + theme_sc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_dual(p_cheat, "CANONICAL_MARKER_HEATMAP", 12, 12)

saveRDS(obj, p("06_objects", "04_annotated_seurat.rds"), compress = FALSE)
stage_message("Stage 04 complete: markers=", nrow(markers),
              "; unresolved clusters=", sum(annotation$final_broad_label == "UNRESOLVED"))
