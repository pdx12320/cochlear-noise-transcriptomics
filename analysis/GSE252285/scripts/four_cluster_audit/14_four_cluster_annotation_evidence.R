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
audit_name <- "cluster_6_9_13_19_annotation_audit"
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

stage_message("Four-cluster annotation evidence audit; no QC/PCA/clustering/object modification")
target_path <- branch_path("06_objects", "02_no_doublet_annotated_tsne.rds")
reference_path <- file.path(project_dir, "02_processed", "reference_final_audit", "GSE217727",
                            "GSE217727_P7_reconstructed_reference.rds")
reference_provenance_path <- file.path(project_dir, "04_tables", "REFERENCE_GSE217727_PROVENANCE.tsv")
required <- c(target_path, reference_path, reference_provenance_path)
if (any(!file.exists(required))) stop("Missing input: ", paste(required[!file.exists(required)], collapse = ", "))

target_clusters <- c("6", "9", "13", "19")
cluster_field <- "RNA_snn_res.0.6"
obj <- readRDS(target_path)
DefaultAssay(obj) <- "RNA"
if (nrow(obj) != 32285L || ncol(obj) != 36111L) stop("Target object invariant failed")
if (!cluster_field %in% colnames(obj[[]])) stop("Fixed resolution field missing")
if (!all(c("pca", "tsne") %in% Reductions(obj)) || "umap" %in% Reductions(obj)) {
  stop("Expected unchanged PCA/t-SNE-only branch reductions")
}
cluster_vec <- as.character(obj[[cluster_field, drop = TRUE]])
if (uniqueN(cluster_vec) != 28L) stop("Expected fixed 28 clusters at resolution 0.6")
target_counts <- data.table(cluster = cluster_vec)[cluster %in% target_clusters, .N, by = cluster]
target_counts <- target_counts[order(as.integer(cluster))]
if (!identical(target_counts$N, c(2359L, 1210L, 838L, 319L))) {
  stop("Target cluster count invariant failed: ", paste(target_counts$N, collapse = ","))
}
Idents(obj) <- factor(cluster_vec, levels = as.character(0:27))

fixed_top20 <- top20_fixed
marker_cache <- audit_path("04_tables", "FOUR_CLUSTER_FINDMARKERS_ALL.tsv.gz")
reuse_markers <- file.exists(marker_cache) && file.info(marker_cache)$size > 1000
if (reuse_markers) {
  stage_message("Reusing completed four-cluster FindMarkers checkpoint")
  markers <- fread(marker_cache, showProgress = FALSE)
} else {
  marker_workers <- suppressWarnings(as.integer(Sys.getenv("FOUR_CLUSTER_MARKER_WORKERS", "4")))
  if (is.na(marker_workers)) marker_workers <- 4L
  marker_workers <- max(1L, min(4L, marker_workers))
  stage_message("Running FindMarkers(cluster vs rest) for clusters ",
                paste(target_clusters, collapse = ","), " with ", marker_workers, " fork workers")
  marker_list <- parallel::mclapply(target_clusters, function(cl) {
    ans <- FindMarkers(
      object = obj, ident.1 = cl, ident.2 = NULL, assay = "RNA",
      only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0,
      test.use = "wilcox", return.thresh = 1, densify = FALSE, verbose = FALSE
    )
    ans <- as.data.table(ans, keep.rownames = "gene")
    ans[, cluster := cl]
    ans
  }, mc.cores = marker_workers, mc.preschedule = FALSE)
  markers <- rbindlist(marker_list, use.names = TRUE, fill = TRUE)
}
fc_col <- intersect(c("avg_log2FC", "avg_logFC"), names(markers))[1]
if (is.na(fc_col)) stop("FindMarkers output lacks avg_log2FC")
if (fc_col != "avg_log2FC") setnames(markers, fc_col, "avg_log2FC")
markers[, cluster := as.character(cluster)]
if (!all(c("gene", "pct.1", "pct.2", "p_val", "p_val_adj") %in% names(markers))) {
  stop("FindMarkers output schema mismatch")
}
markers[, pct_diff := pct.1 - pct.2]
markers[, effect_specificity_score := avg_log2FC * pmax(pct_diff, 0)]
markers[, top20_candidate_excluded_from_annotation := gene %in% fixed_top20]
setcolorder(markers, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "pct_diff",
                      "effect_specificity_score", "p_val", "p_val_adj",
                      "top20_candidate_excluded_from_annotation"))
setorder(markers, cluster, -effect_specificity_score, -avg_log2FC, -pct_diff, gene)
fwrite(markers, marker_cache, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

top50 <- markers[, head(.SD, 50L), by = cluster]
top50_annotation <- markers[top20_candidate_excluded_from_annotation == FALSE,
                            head(.SD, 50L), by = cluster]
specific <- markers[
  top20_candidate_excluded_from_annotation == FALSE &
    avg_log2FC >= 0.50 & pct.1 >= 0.25 & pct_diff >= 0.20
]
write_tsv(top50, audit_path("04_tables", "FOUR_CLUSTER_TOP50_MARKERS.tsv"))
write_tsv(top50_annotation, audit_path("04_tables", "FOUR_CLUSTER_TOP50_ANNOTATION_MARKERS.tsv"))
write_tsv(specific, audit_path("04_tables", "FOUR_CLUSTER_SPECIFIC_MARKERS.tsv"))

# Expanded programs intentionally extend beyond the user's seed list and exclude fixed Top20 candidates.
programs <- list(
  SGN_canonical_neuronal = c("Tubb3", "Snap25", "Syt1", "Rbfox3", "Nefh", "Prph", "Nefl", "Nefm",
                             "Stmn2", "Elavl3", "Slc17a7", "Pou4f1", "Calb2", "Calb1", "Bhlhe22", "Epha4", "Isl1"),
  Schwann = c("Mpz", "Plp1", "Sox10", "S100b", "S100a1", "Gjb1", "Prx", "Pou3f1", "Mbp", "Mag", "Egr2"),
  Satellite_glia = c("Fabp7", "S100b", "Sox10", "Plp1", "Glul", "Slc1a3", "Kcnj10", "Gja1", "Apoe", "Sparc"),
  Chondrocyte = c("Coch", "Cald1", "Col2a1", "Acan", "Sox9", "Comp", "Col9a1", "Col9a2", "Col9a3", "Chad", "Chodl"),
  Fibrocyte_mesenchymal = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra", "Col6a1", "Pou3f4", "Tbx18", "Coch", "Otor", "Igfbp5", "Rgs2"),
  Supporting_epithelial = c("Sox2", "Gjb2", "Gjb6", "Krt8", "Krt18", "Epcam", "Krt19", "Fbxo2", "Cldn3", "Cldn9", "Sox9", "Cdh2", "Lbsp", "Bgn", "Tpm1", "Apod", "Cavin2", "Gas1"),
  Hair_sensory = c("Pou4f3", "Myo7a", "Otof", "Slc26a5", "Ocm", "Pvalb", "Lhfpl5", "Tmc1", "Grxcr2", "Xirp2", "Cnmd", "S100a16"),
  Endothelial = c("Pecam1", "Cdh5", "Kdr", "Emcn", "Vwf", "Esam", "Klf2", "Ptprb", "Cldn5", "Egfl7", "Flt1", "Klf4"),
  Pericyte_mural = c("Rgs5", "Pdgfrb", "Cspg4", "Kcnj8", "Abcc9", "Des", "Notch3", "Col4a1", "Col4a2", "Rbp1", "Mcam"),
  SMC = c("Acta2", "Tagln", "Myh11", "Cnn1", "Myocd", "Myl9", "Des", "Mustn1", "Smtn"),
  Immune_myeloid = c("Ptprc", "Lyz2", "Tyrobp", "Aif1", "C1qa", "C1qb", "C1qc", "Cd68", "Ctss", "Fcer1g", "Csf1r", "Fcgr3"),
  Hematopoietic_progenitor = c("Cd34", "Kit", "Flt3", "Myb", "Hoxa9", "Meis1", "Msi2",
                               "Runx1", "Gata2", "Procr", "Mpl", "Dntt", "Bcl11a", "Ikzf2",
                               "Mef2c", "Ptprc", "Ptprcap"),
  Stria_intermediate_melanocyte = c("Mitf", "Kit", "Dct", "Tyr", "Tyrp1", "Slc45a2", "Kcnj10", "Met"),
  Erythroid = c("Klf1", "Gata1", "Alas2", "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Gypa", "Slc4a1")
)
programs <- lapply(programs, function(x) unique(setdiff(x, fixed_top20)))

author_panels <- list(
  HC_author = c("Cnmd", "S100a16"),
  DC_PC_author = c("Col1a1", "Gas1"),
  IPhC_IBC_author = c("Cdh2", "Lbsp"),
  HeC_author = c("Bgn", "Tpm1"),
  OSC_ISC_author = c("Apod", "Cavin2"),
  SGN_author = c("Car1", "Ces2g"),
  SGC_author = c("Plk1", "Rad18"),
  Schwann_author = c("Mpz", "Plp1"),
  Chondrocyte_author = c("Coch", "Cald1"),
  CEC_author = c("S100a4", "Vcan"),
  Fibroblast_author = c("Igfbp5", "Rgs2"),
  SMC_author = c("Acta2", "Myl9"),
  SGN_type_IA_author = c("Ligp1"),
  SGN_type_IB_author = c("Mgst3"),
  SGN_type_IC_author = c("Hist1h1a"),
  SGN_type_II_author = c("Arhgdib")
)
author_source <- rbindlist(lapply(names(author_panels), function(panel) {
  data.table(
    author_cell_type = panel,
    gene = author_panels[[panel]],
    source = "Jiang et al., Cell Biology and Toxicology 2024, DOI 10.1007/s10565-024-09912-2",
    source_location = ifelse(grepl("SGN_type", panel), "Results page 12 / Supplementary Fig S1", "Results page 8 / Fig 2 text")
  )
}))
write_tsv(author_source, audit_path("00_manifest", "TARGET_PAPER_AUTHOR_MARKERS.tsv"))

data_mat <- LayerData(obj[["RNA"]], layer = "data")
all_panel_genes <- unique(c(unlist(programs, use.names = FALSE), unlist(author_panels, use.names = FALSE)))
present <- intersect(all_panel_genes, rownames(data_mat))
gene_stats_rows <- list()
for (cl in target_clusters) {
  inside <- which(cluster_vec == cl)
  outside <- which(cluster_vec != cl)
  for (gene in present) {
    xin <- data_mat[gene, inside, drop = TRUE]
    xout <- data_mat[gene, outside, drop = TRUE]
    gene_stats_rows[[length(gene_stats_rows) + 1L]] <- data.table(
      cluster = cl, gene = gene,
      mean_log_normalized_in = mean(xin),
      mean_log_normalized_out = mean(xout),
      mean_log_normalized_diff = mean(xin) - mean(xout),
      pct_in = mean(xin > 0), pct_out = mean(xout > 0),
      pct_diff = mean(xin > 0) - mean(xout > 0)
    )
  }
}
gene_stats <- rbindlist(gene_stats_rows)
program_map <- rbindlist(lapply(names(programs), function(program) {
  data.table(program = program, gene = intersect(programs[[program]], present))
}))
program_gene_stats <- merge(gene_stats, program_map, by = "gene", allow.cartesian = TRUE)
setcolorder(program_gene_stats, c("cluster", "program", "gene", setdiff(names(program_gene_stats), c("cluster", "program", "gene"))))
write_tsv(program_gene_stats, audit_path("04_tables", "FOUR_CLUSTER_PROGRAM_GENE_STATS.tsv"))

program_scores <- program_gene_stats[, .(
  genes_present = .N,
  supported_gene_number = sum(pct_in >= 0.10 & pct_diff >= 0.05 & mean_log_normalized_diff >= 0.05),
  supported_genes = {
    z <- .SD[pct_in >= 0.10 & pct_diff >= 0.05 & mean_log_normalized_diff >= 0.05]
    if (nrow(z)) paste(z[order(-pct_diff, -mean_log_normalized_diff), gene], collapse = ",") else "NONE"
  },
  program_effect_score = mean(pmax(pct_diff, 0) * pmax(mean_log_normalized_diff, 0)),
  mean_pct_in = mean(pct_in), mean_pct_diff = mean(pct_diff)
), by = .(cluster, program)]
program_scores[, program_rank := frank(-program_effect_score, ties.method = "min"), by = cluster]
program_scores <- program_scores[order(as.integer(cluster), program_rank)]

top_overlap <- merge(top50_annotation[, .(cluster, gene, effect_specificity_score)],
                     program_map, by = "gene", allow.cartesian = TRUE)
top_program <- top_overlap[, .(
  top50_overlap_number = uniqueN(gene),
  top50_overlap_genes = paste(gene[order(-effect_specificity_score)], collapse = ","),
  top50_overlap_score = sum(effect_specificity_score)
), by = .(cluster, program)]
program_scores <- merge(program_scores, top_program, by = c("cluster", "program"), all.x = TRUE)
program_scores[is.na(top50_overlap_number), `:=`(top50_overlap_number = 0L,
                                                 top50_overlap_genes = "NONE",
                                                 top50_overlap_score = 0)]
program_scores <- program_scores[order(as.integer(cluster), program_rank)]
write_tsv(program_scores, audit_path("04_tables", "FOUR_CLUSTER_LINEAGE_PROGRAM_SCORES.tsv"))

author_map <- rbindlist(lapply(names(author_panels), function(panel) {
  data.table(author_panel = panel, gene = intersect(author_panels[[panel]], present))
}))
author_stats <- merge(gene_stats, author_map, by = "gene", allow.cartesian = TRUE)
setcolorder(author_stats, c("cluster", "author_panel", "gene", setdiff(names(author_stats), c("cluster", "author_panel", "gene"))))
write_tsv(author_stats, audit_path("04_tables", "FOUR_CLUSTER_AUTHOR_MARKER_STATS.tsv"))

stage_message("Rebuilding reference PCA only and mapping four target clusters to GSE217727 P7")
ref <- readRDS(reference_path)
DefaultAssay(ref) <- "RNA"
if (ncol(ref) != 7711L || nrow(ref) != 19285L) stop("Reference invariant failed")
confident <- as.logical(ref$reference_signature_confident)
label_counts <- table(as.character(ref$reference_subtype_reconstructed[confident]))
eligible_labels <- names(label_counts[label_counts >= 20L])
ref_cells <- colnames(ref)[confident & as.character(ref$reference_subtype_reconstructed) %in% eligible_labels]
ref_train <- subset(ref, cells = ref_cells)
ref_train$reference_subtype <- factor(as.character(ref_train$reference_subtype_reconstructed), levels = eligible_labels)
ref_train$reference_broad <- factor(as.character(ref_train$reference_broad_reconstructed))
ref_features <- intersect(VariableFeatures(ref), rownames(obj))
if (length(ref_features) < 2000L) stop("Insufficient shared reference variable features")
VariableFeatures(ref_train) <- ref_features
ref_train <- ScaleData(ref_train, features = ref_features, verbose = FALSE)
ref_train <- RunPCA(ref_train, features = ref_features, npcs = 30, seed.use = 42, verbose = FALSE)

query_cells <- colnames(obj)[cluster_vec %in% target_clusters]
query <- subset(obj, cells = query_cells)
if (ncol(query) != 4726L || !"pca" %in% Reductions(query)) stop("Query subset invariant failed")
anchors <- FindTransferAnchors(
  reference = ref_train, query = query,
  normalization.method = "LogNormalize", reference.assay = "RNA", query.assay = "RNA",
  reduction = "pcaproject", features = ref_features, dims = 1:30,
  k.anchor = 5, k.filter = 200, k.score = 30, verbose = TRUE
)
pred_subtype <- TransferData(anchorset = anchors, refdata = ref_train$reference_subtype,
                             weight.reduction = query[["pca"]], dims = 1:30,
                             k.weight = 50, verbose = TRUE)
pred_broad <- TransferData(anchorset = anchors, refdata = ref_train$reference_broad,
                           weight.reduction = query[["pca"]], dims = 1:30,
                           k.weight = 50, verbose = TRUE)

pred_subtype <- as.data.table(as.data.frame(pred_subtype), keep.rownames = "barcode")
pred_broad <- as.data.table(as.data.frame(pred_broad), keep.rownames = "barcode")
setnames(pred_subtype, setdiff(names(pred_subtype), "barcode"),
         paste0("subtype_", setdiff(names(pred_subtype), "barcode")))
setnames(pred_broad, setdiff(names(pred_broad), "barcode"),
         paste0("broad_", setdiff(names(pred_broad), "barcode")))
reference_cell <- merge(pred_subtype, pred_broad, by = "barcode", all = TRUE, sort = FALSE)
reference_cell[, cluster := cluster_vec[match(barcode, colnames(obj))]]

score_margin <- function(dt, prefix) {
  score_cols <- grep(paste0("^", prefix, "_prediction\\.score\\."), names(dt), value = TRUE)
  score_cols <- setdiff(score_cols, paste0(prefix, "_prediction.score.max"))
  mat <- as.matrix(dt[, ..score_cols])
  sorted <- t(apply(mat, 1, sort, decreasing = TRUE))
  sorted[, 1] - sorted[, 2]
}
reference_cell[, subtype_prediction_raw := as.character(subtype_predicted.id)]
reference_cell[, subtype_score := subtype_prediction.score.max]
reference_cell[, subtype_margin := score_margin(reference_cell, "subtype")]
reference_cell[, subtype_rejected := subtype_prediction_raw %in% c("Unknown", "UNRESOLVED") |
                 subtype_score < 0.50 | subtype_margin < 0.10]
reference_cell[, subtype_prediction_rejected := fifelse(subtype_rejected, "UNRESOLVED", subtype_prediction_raw)]
reference_cell[, broad_prediction_raw := as.character(broad_predicted.id)]
reference_cell[, broad_score := broad_prediction.score.max]
reference_cell[, broad_margin := score_margin(reference_cell, "broad")]
reference_cell[, broad_rejected := broad_prediction_raw == "UNRESOLVED" |
                 broad_score < 0.50 | broad_margin < 0.10]
reference_cell[, broad_prediction_rejected := fifelse(broad_rejected, "UNRESOLVED", broad_prediction_raw)]
reference_cell[, query_order__ := match(barcode, query_cells)]
setorder(reference_cell, query_order__)
reference_cell[, query_order__ := NULL]
fwrite(reference_cell, audit_path("04_tables", "FOUR_CLUSTER_REFERENCE_MAPPING_CELL_LEVEL.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

reference_summary <- reference_cell[, {
  subtype_raw_tab <- sort(table(subtype_prediction_raw), decreasing = TRUE)
  subtype_acc <- .SD[subtype_rejected == FALSE]
  subtype_acc_tab <- sort(table(subtype_acc$subtype_prediction_rejected), decreasing = TRUE)
  subtype_top <- if (length(subtype_acc_tab)) names(subtype_acc_tab)[1] else "UNRESOLVED"
  subtype_rows <- subtype_acc[subtype_prediction_rejected == subtype_top]
  broad_raw_tab <- sort(table(broad_prediction_raw), decreasing = TRUE)
  broad_acc <- .SD[broad_rejected == FALSE]
  broad_acc_tab <- sort(table(broad_acc$broad_prediction_rejected), decreasing = TRUE)
  broad_top <- if (length(broad_acc_tab)) names(broad_acc_tab)[1] else "UNRESOLVED"
  broad_rows <- broad_acc[broad_prediction_rejected == broad_top]
  list(
    cell_number = .N,
    subtype_top_raw = names(subtype_raw_tab)[1],
    subtype_top_raw_fraction = unname(subtype_raw_tab[1]) / .N,
    subtype_top_accepted = subtype_top,
    subtype_top_accepted_fraction_total = if (length(subtype_acc_tab)) unname(subtype_acc_tab[1]) / .N else 0,
    subtype_top_accepted_fraction_among_accepted = if (nrow(subtype_acc)) unname(subtype_acc_tab[1]) / nrow(subtype_acc) else 0,
    subtype_rejected_fraction = mean(subtype_rejected),
    subtype_median_score_top = if (nrow(subtype_rows)) median(subtype_rows$subtype_score) else NA_real_,
    subtype_median_margin_top = if (nrow(subtype_rows)) median(subtype_rows$subtype_margin) else NA_real_,
    broad_top_raw = names(broad_raw_tab)[1],
    broad_top_raw_fraction = unname(broad_raw_tab[1]) / .N,
    broad_top_accepted = broad_top,
    broad_top_accepted_fraction_total = if (length(broad_acc_tab)) unname(broad_acc_tab[1]) / .N else 0,
    broad_top_accepted_fraction_among_accepted = if (nrow(broad_acc)) unname(broad_acc_tab[1]) / nrow(broad_acc) else 0,
    broad_rejected_fraction = mean(broad_rejected),
    broad_median_score_top = if (nrow(broad_rows)) median(broad_rows$broad_score) else NA_real_,
    broad_median_margin_top = if (nrow(broad_rows)) median(broad_rows$broad_margin) else NA_real_
  )
}, by = cluster]
reference_summary <- reference_summary[order(as.integer(cluster))]
write_tsv(reference_summary, audit_path("04_tables", "FOUR_CLUSTER_REFERENCE_MAPPING_SUMMARY.tsv"))

reference_fractions <- reference_cell[, .N, by = .(cluster, subtype_prediction_rejected)]
reference_fractions[, fraction := N / sum(N), by = cluster]
reference_fractions <- reference_fractions[order(as.integer(cluster), -fraction)]
write_tsv(reference_fractions, audit_path("04_tables", "FOUR_CLUSTER_REFERENCE_SUBTYPE_FRACTIONS.tsv"))

stage_message("Technical and library composition audit")
meta <- as.data.table(obj[[]], keep.rownames = "barcode")
meta[, cluster := cluster_vec]
technical_rows <- list()
library_rows <- list()
for (cl in target_clusters) {
  inside <- meta[cluster == cl]
  outside <- meta[cluster != cl]
  lib_in <- inside[, .N, by = library_id][, fraction := N / sum(N)]
  lib_out <- outside[, .N, by = library_id][, fraction := N / sum(N)]
  lib_in[, `:=`(cluster = cl, group = "CLUSTER")]
  lib_out[, `:=`(cluster = cl, group = "REST")]
  library_rows[[length(library_rows) + 1L]] <- rbind(lib_in, lib_out)
  in_frac <- setNames(lib_in$fraction, lib_in$library_id)
  out_frac <- setNames(lib_out$fraction, lib_out$library_id)
  libs <- sort(unique(c(names(in_frac), names(out_frac))))
  in_vec <- in_frac[libs]; out_vec <- out_frac[libs]
  in_vec[is.na(in_vec)] <- 0; out_vec[is.na(out_vec)] <- 0
  robust_z <- function(in_x, out_x) {
    denom <- mad(out_x, constant = 1.4826, na.rm = TRUE)
    (median(in_x, na.rm = TRUE) - median(out_x, na.rm = TRUE)) / max(denom, 1e-12)
  }
  technical_rows[[length(technical_rows) + 1L]] <- data.table(
    cluster = cl, cell_number = nrow(inside),
    LIB1_fraction = in_frac[["LIB1"]], LIB2_fraction = in_frac[["LIB2"]], LIB3_fraction = in_frac[["LIB3"]],
    dominant_library = names(which.max(in_frac)), max_library_fraction = max(in_frac),
    library_total_variation_vs_rest = 0.5 * sum(abs(in_vec - out_vec)),
    median_nCount_RNA = median(inside$nCount_RNA), median_nFeature_RNA = median(inside$nFeature_RNA),
    median_percent_mt = median(inside$percent.mt),
    robust_z_median_nCount_vs_rest = robust_z(inside$nCount_RNA, outside$nCount_RNA),
    robust_z_median_nFeature_vs_rest = robust_z(inside$nFeature_RNA, outside$nFeature_RNA),
    robust_z_median_percent_mt_vs_rest = robust_z(inside$percent.mt, outside$percent.mt)
  )
}
technical <- rbindlist(technical_rows)
technical[, technical_warning := fifelse(
  max_library_fraction >= 0.90 | library_total_variation_vs_rest >= 0.50 |
    abs(robust_z_median_nCount_vs_rest) >= 3 |
    abs(robust_z_median_nFeature_vs_rest) >= 3 |
    abs(robust_z_median_percent_mt_vs_rest) >= 3,
  "POTENTIAL_LIBRARY_OR_QC_DRIVER", "NONE"
)]
library_composition <- rbindlist(library_rows)
write_tsv(technical, audit_path("04_tables", "FOUR_CLUSTER_TECHNICAL_AUDIT.tsv"))
write_tsv(library_composition, audit_path("04_tables", "FOUR_CLUSTER_LIBRARY_COMPOSITION.tsv"))

# Cluster-level evidence bundle for manual final adjudication.
marker_strings <- top50_annotation[, .(
  top20_effect_markers = paste(head(gene, 20L), collapse = ","),
  top10_specificity_markers = paste(head(gene[order(-pct_diff, -avg_log2FC)], 10L), collapse = ",")
), by = cluster]
program_strings <- program_scores[order(program_rank), .(
  top_program_candidates = paste0(head(program, 5L), "[", head(supported_genes, 5L), "]", collapse = ";"),
  top_program_score_text = paste0(head(program, 5L), "=", sprintf("%.4f", head(program_effect_score, 5L)), collapse = ";")
), by = cluster]
sgn_author <- dcast(author_stats[author_panel == "SGN_author"], cluster ~ gene, value.var = "pct_in")
setnames(sgn_author, intersect(c("Car1", "Ces2g"), names(sgn_author)),
         paste0(intersect(c("Car1", "Ces2g"), names(sgn_author)), "_pct_in"))
sgn_canonical <- program_scores[program == "SGN_canonical_neuronal",
  .(cluster, canonical_neuronal_supported_gene_number = supported_gene_number,
    canonical_neuronal_supported_genes = supported_genes,
    canonical_neuronal_program_score = program_effect_score)]
evidence_bundle <- Reduce(function(x, y) merge(x, y, by = "cluster", all = TRUE),
  list(target_counts[, .(cluster, cell_number = N)], marker_strings, program_strings,
       sgn_author, sgn_canonical, reference_summary, technical))
evidence_bundle <- evidence_bundle[order(as.integer(cluster))]
write_tsv(evidence_bundle, audit_path("04_tables", "FOUR_CLUSTER_EVIDENCE_BUNDLE.tsv"))

# Figure: marker-program dot plot (pct in and pct difference; full numeric table is retained).
display_programs <- c("Erythroid", "Hematopoietic_progenitor", "SGN_canonical_neuronal", "Schwann", "Satellite_glia", "Chondrocyte",
                      "Fibrocyte_mesenchymal", "Supporting_epithelial", "Hair_sensory",
                      "Endothelial", "Pericyte_mural", "SMC", "Immune_myeloid")
display_genes <- unique(unlist(lapply(display_programs, function(program_name) {
  genes <- programs[[program_name]]
  scored <- program_gene_stats[program == program_name & gene %in% genes,
                               .(max_effect = max(pmax(pct_diff, 0) * pmax(mean_log_normalized_diff, 0))), by = gene]
  scored[order(-max_effect), head(gene, 6L)]
}), use.names = FALSE))
dot_dt <- program_gene_stats[program %in% display_programs & gene %in% display_genes]
dot_dt[, gene := factor(gene, levels = rev(display_genes))]
dot_dt[, cluster := factor(cluster, levels = target_clusters)]
p_program_dot <- ggplot(dot_dt, aes(cluster, gene, size = 100 * pct_in, colour = pct_diff)) +
  geom_point() +
  scale_size_area(max_size = 6, limits = c(0, 100), breaks = c(10, 25, 50, 75)) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0) +
  labs(x = "Cluster", y = NULL, size = "% expressing", colour = "pct.1-pct.2",
       title = "Expanded lineage-marker audit",
       subtitle = "Dot size shows within-cluster detection; colour shows specificity versus all other cells") +
  theme_branch()
save_audit_dual(p_program_dot, "FOUR_CLUSTER_LINEAGE_MARKER_DOTPLOT", 11, 20)

score_plot <- copy(program_scores[program %in% display_programs])
score_plot[, score_z := safe_zscore(program_effect_score), by = program]
score_plot[, cluster := factor(cluster, levels = target_clusters)]
score_plot[, program := factor(program, levels = rev(display_programs))]
p_score <- ggplot(score_plot, aes(cluster, program, fill = score_z)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = supported_gene_number), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(x = "Cluster", y = NULL, fill = "program-wise z",
       title = "Lineage program evidence",
       subtitle = "Numbers are supported genes; colour is relative program effect across the four audited clusters") +
  theme_branch()
save_audit_dual(p_score, "FOUR_CLUSTER_LINEAGE_SCORE_HEATMAP", 9, 8)

sgn_genes <- intersect(unique(c(author_panels$SGN_author, programs$SGN_canonical_neuronal,
                                "Ligp1", "Mgst3", "Hist1h1a", "Arhgdib")), present)
sgn_dt <- gene_stats[gene %in% sgn_genes]
sgn_dt[, gene := factor(gene, levels = rev(sgn_genes))]
sgn_dt[, cluster := factor(cluster, levels = target_clusters)]
p_sgn <- ggplot(sgn_dt, aes(cluster, gene, size = 100 * pct_in, colour = pct_diff)) +
  geom_point() + scale_size_area(max_size = 7, limits = c(0, 100)) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0) +
  labs(x = "Cluster", y = NULL, size = "% expressing", colour = "pct.1-pct.2",
       title = "Author SGN markers versus canonical neuronal program",
       subtitle = "Car1/Ces2g positivity is not sufficient without the multi-gene neuronal program") +
  theme_branch()
save_audit_dual(p_sgn, "FOUR_CLUSTER_SGN_AUTHOR_VS_CANONICAL_DOTPLOT", 9, 9)

ref_plot <- copy(reference_fractions)
ref_plot[, cluster := factor(cluster, levels = target_clusters)]
subtype_levels <- unique(ref_plot[order(-fraction), subtype_prediction_rejected])
ref_plot[, subtype_prediction_rejected := factor(subtype_prediction_rejected, levels = rev(subtype_levels))]
p_ref <- ggplot(ref_plot, aes(cluster, subtype_prediction_rejected, fill = fraction)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(fraction >= 0.03, sprintf("%.1f%%", 100 * fraction), "")), size = 2.8) +
  scale_fill_viridis_c(option = "C", limits = c(0, 1), labels = scales::percent_format()) +
  labs(x = "Target cluster", y = "GSE217727 prediction after rejection", fill = "Cell fraction",
       title = "Mouse P7 cochlea reference mapping",
       subtitle = "Closed-world predictions with score <0.50 or margin <0.10 are UNRESOLVED") +
  theme_branch()
save_audit_dual(p_ref, "FOUR_CLUSTER_REFERENCE_MAPPING_HEATMAP", 10, 9)

library_composition[, cluster := factor(cluster, levels = target_clusters)]
library_composition[, group := factor(group, levels = c("CLUSTER", "REST"))]
technical[, cluster := factor(cluster, levels = target_clusters)]
p_lib <- ggplot(library_composition, aes(group, fraction, fill = library_id)) +
  geom_col() + facet_wrap(~cluster, nrow = 1) +
  scale_fill_manual(values = c(LIB1 = "#0072B2", LIB2 = "#D55E00", LIB3 = "#009E73")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Library fraction", fill = NULL,
       title = "Library composition: cluster versus rest") + theme_branch()
tech_long <- melt(technical,
  id.vars = "cluster",
  measure.vars = c("robust_z_median_nCount_vs_rest", "robust_z_median_nFeature_vs_rest",
                   "robust_z_median_percent_mt_vs_rest"),
  variable.name = "metric", value.name = "robust_z")
p_tech <- ggplot(tech_long, aes(cluster, robust_z, colour = metric, group = metric)) +
  geom_hline(yintercept = c(-3, 3), linetype = 2, colour = "grey50") +
  geom_point(size = 2.8) + geom_line() +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Cluster", y = "Robust median shift vs rest (MAD units)", colour = NULL,
       title = "QC-covariate audit") + theme_branch()
save_audit_dual(p_lib / p_tech, "FOUR_CLUSTER_LIBRARY_QC_AUDIT", 12, 10)

top_heat_genes <- unique(top50_annotation[, head(gene, 10L), by = cluster]$V1)
mean_rows <- lapply(target_clusters, function(cl) {
  cells <- which(cluster_vec == cl)
  data.table(cluster = cl, gene = top_heat_genes,
             mean_log_expression = as.numeric(Matrix::rowMeans(data_mat[top_heat_genes, cells, drop = FALSE])))
})
top_heat <- rbindlist(mean_rows)
top_heat[, gene_z := safe_zscore(mean_log_expression), by = gene]
top_heat[, cluster := factor(cluster, levels = target_clusters)]
top_heat[, gene := factor(gene, levels = rev(top_heat_genes))]
p_top_heat <- ggplot(top_heat, aes(cluster, gene, fill = gene_z)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(x = "Cluster", y = NULL, fill = "gene-wise z",
       title = "Top effect/specificity markers") + theme_branch()
save_audit_dual(p_top_heat, "FOUR_CLUSTER_TOP_MARKER_HEATMAP", 9, 14)

parameters <- data.table(
  parameter = c("target_object", "target_clusters", "fixed_resolution", "FindMarkers_test",
                "FindMarkers_min_pct", "FindMarkers_logfc_threshold", "marker_ranking",
                "reference", "reference_mapping", "broad_reject_score", "broad_reject_margin",
                "subtype_reject_score", "subtype_reject_margin", "target_PCA_rerun",
                "target_clustering_rerun", "target_object_modified", "condition_DE"),
  value = c(target_path, paste(target_clusters, collapse = ","), "0.6", "wilcox",
            "0.05", "0", "effect_specificity_score=avg_log2FC*max(pct.1-pct.2,0)",
            "GSE217727_P7_C57BL6J_whole_cochlea", "Seurat_FindTransferAnchors_TransferData_pcaproject_30PC",
            "0.50", "0.10", "0.50", "0.10", "NO", "NO", "NO", "NOT_PERFORMED")
)
write_tsv(parameters, audit_path("00_manifest", "FOUR_CLUSTER_AUDIT_PARAMETERS.tsv"))
file.copy(reference_provenance_path, audit_path("00_manifest", "REFERENCE_GSE217727_PROVENANCE.tsv"), overwrite = TRUE)
capture.output(sessionInfo(), file = audit_path("07_logs", "sessionInfo.txt"))
stage_message("Evidence stage complete")
