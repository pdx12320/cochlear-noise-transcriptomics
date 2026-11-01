source(file.path(Sys.getenv("GSE252285_SCRIPT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "03_scripts")), "common.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(patchwork)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
})

stage_message("Stage 02: Seurat object, QC, and per-library doublet audit")
counts <- readRDS(p("02_processed", "01_counts_matrix.rds"))
features <- fread(p("00_manifest", "FEATURE_NAME_MAP.tsv.gz"), showProgress = FALSE)
map <- fread(p("00_manifest", "library_map_template.tsv"), showProgress = FALSE)
if (!all(c("suffix", "library_id", "condition") %in% names(map))) {
  stop("library_map_template.tsv must contain suffix, library_id, and condition")
}

barcodes <- colnames(counts)
suffix <- barcode_suffix(barcodes)
meta <- data.frame(
  barcode = barcodes,
  barcode_suffix = suffix,
  library_id = map$library_id[match(suffix, map$suffix)],
  condition = map$condition[match(suffix, map$suffix)],
  row.names = barcodes,
  check.names = FALSE
)
if (anyNA(meta$library_id) || anyNA(meta$condition)) {
  stop("One or more barcode suffixes are absent from library_map_template.tsv")
}

obj <- CreateSeuratObject(counts = counts, project = "GSE252285",
                          assay = "RNA", min.cells = 0, min.features = 0,
                          meta.data = meta)
feature_md <- data.frame(
  feature_id = features$feature_id,
  original_symbol = features$gene_symbol,
  feature_type = features$feature_type,
  row.names = features$feature_name,
  check.names = FALSE
)
obj[["RNA"]] <- SeuratObject::AddMetaData(obj[["RNA"]], metadata = feature_md)
saveRDS(obj, p("06_objects", "01_raw_seurat.rds"), compress = FALSE)
rm(counts)
gc(verbose = FALSE)

mt_features <- rownames(obj)[grepl("^mt-", rownames(obj), ignore.case = TRUE)]
obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt_features)
qc <- as.data.table(obj[[]])

thresholds <- qc[, {
  nf <- robust_upper(nFeature_RNA)
  nc <- robust_upper(nCount_RNA)
  mt <- robust_upper(percent.mt)
  list(
    initial_cells = .N,
    min_nFeature_RNA = 200,
    max_nFeature_RNA = floor(nf$limit),
    max_nCount_RNA = floor(nc$limit),
    max_percent_mt = max(10, min(30, mt$limit)),
    nFeature_median = nf$median,
    nFeature_MAD = nf$mad,
    nFeature_q99.5 = nf$q995,
    nFeature_median_plus_6MAD = nf$mad_limit,
    nCount_median = nc$median,
    nCount_MAD = nc$mad,
    nCount_q99.5 = nc$q995,
    nCount_median_plus_6MAD = nc$mad_limit,
    percent_mt_median = mt$median,
    percent_mt_MAD = mt$mad,
    percent_mt_q99.5 = mt$q995,
    percent_mt_median_plus_6MAD = mt$mad_limit,
    rule = "nFeature>=200; high nFeature/nCount=min(q99.5,median+6*MAD); mt=min(q99.5,median+6*MAD), bounded 10-30%"
  )
}, by = library_id]
write_tsv(thresholds, p("04_tables", "QC_THRESHOLDS.tsv"))

qc <- merge(qc, thresholds[, .(library_id, min_nFeature_RNA, max_nFeature_RNA,
                                max_nCount_RNA, max_percent_mt)],
            by = "library_id", all.x = TRUE, sort = FALSE)
qc[, `:=`(
  pass_min_features = nFeature_RNA >= min_nFeature_RNA,
  pass_max_features = nFeature_RNA <= max_nFeature_RNA,
  pass_max_counts = nCount_RNA <= max_nCount_RNA,
  pass_mito = percent.mt <= max_percent_mt
)]
qc[, qc_pass := pass_min_features & pass_max_features & pass_max_counts & pass_mito]

scenario_defs <- data.table(
  scenario = c("PERMISSIVE", "MODERATE", "STRINGENT", "MAIN_ROBUST"),
  min_features = c(200, 300, 500, NA_real_),
  max_mito = c(20, 15, 10, NA_real_)
)
sensitivity <- rbindlist(lapply(seq_len(nrow(scenario_defs)), function(i) {
  s <- scenario_defs[i]
  qc[, {
    th <- thresholds[library_id == .BY$library_id]
    keep <- if (s$scenario == "MAIN_ROBUST") {
      qc_pass
    } else {
      nFeature_RNA >= s$min_features & nFeature_RNA <= th$max_nFeature_RNA &
        nCount_RNA <= th$max_nCount_RNA & percent.mt <= s$max_mito
    }
    .(initial_cells = .N, retained_cells = sum(keep), retained_fraction = mean(keep))
  }, by = library_id][, scenario := s$scenario]
}))
setcolorder(sensitivity, c("scenario", "library_id", "initial_cells", "retained_cells", "retained_fraction"))
write_tsv(sensitivity, p("09_qc", "QC_SENSITIVITY.tsv"))

plot_dt <- copy(qc)
p_nf <- VlnPlot(obj, features = "nFeature_RNA", group.by = "library_id", pt.size = 0) + theme_sc()
p_nc <- VlnPlot(obj, features = "nCount_RNA", group.by = "library_id", pt.size = 0) + theme_sc()
p_mt <- VlnPlot(obj, features = "percent.mt", group.by = "library_id", pt.size = 0) + theme_sc()
save_dual(p_nf, "QC_nFeature_VIOLIN", 8, 5)
save_dual(p_nc, "QC_nCount_VIOLIN", 8, 5)
save_dual(p_mt, "QC_percent_mt_VIOLIN", 8, 5)

p_sc1 <- ggplot(plot_dt, aes(nCount_RNA, nFeature_RNA, colour = library_id)) +
  geom_point(alpha = 0.20, size = 0.25) + scale_x_log10() + scale_y_log10() +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") + theme_sc()
p_sc2 <- ggplot(plot_dt, aes(nCount_RNA, percent.mt, colour = library_id)) +
  geom_point(alpha = 0.20, size = 0.25) + scale_x_log10() +
  labs(x = "nCount_RNA (log10)", y = "percent.mt") + theme_sc()
save_dual(p_sc1, "QC_nFeature_vs_nCount_SCATTER", 8, 6)
save_dual(p_sc2, "QC_nCount_vs_percent_mt_SCATTER", 8, 6)

p_den <- ggplot(plot_dt, aes(nFeature_RNA, colour = library_id, fill = library_id)) +
  geom_density(alpha = 0.12) + scale_x_log10() + theme_sc()
p_hist <- ggplot(plot_dt, aes(percent.mt, fill = library_id)) +
  geom_histogram(bins = 60, alpha = 0.65, position = "identity") +
  facet_wrap(~library_id, scales = "free_y") + theme_sc()
save_dual(p_den / p_hist, "QC_DENSITY_HISTOGRAM", 10, 10)

lib_initial <- plot_dt[, .N, by = library_id]
p_lib <- ggplot(lib_initial, aes(library_id, N, fill = library_id)) +
  geom_col(show.legend = FALSE) + geom_text(aes(label = N), vjust = -0.25) +
  labs(x = NULL, y = "Initial cells") + theme_sc()
save_dual(p_lib, "QC_LIBRARY_CELL_COUNTS", 7, 5)

obj$qc_pass <- qc$qc_pass[match(colnames(obj), qc$barcode)]
qc_obj <- subset(obj, subset = qc_pass)
rm(obj)
gc(verbose = FALSE)

stage_message("Running scDblFinder independently for each library")
counts_qc <- SeuratObject::LayerData(qc_obj[["RNA"]], layer = "counts")
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_qc))
doublet_list <- vector("list", length(unique(qc_obj$library_id)))
libs <- sort(unique(as.character(qc_obj$library_id)))
for (i in seq_along(libs)) {
  lib <- libs[i]
  idx <- which(qc_obj$library_id == lib)
  stage_message("scDblFinder: ", lib, " (", length(idx), " cells)")
  set.seed(seed)
  sce_lib <- scDblFinder::scDblFinder(
    sce[, idx, drop = FALSE], samples = rep(lib, length(idx)), clusters = FALSE,
    BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
  )
  doublet_list[[i]] <- data.table(
    barcode = colnames(sce_lib),
    library_id = lib,
    scDblFinder_score = as.numeric(SummarizedExperiment::colData(sce_lib)$scDblFinder.score),
    scDblFinder_class = as.character(SummarizedExperiment::colData(sce_lib)$scDblFinder.class)
  )
  rm(sce_lib)
  gc(verbose = FALSE)
}
doublets <- rbindlist(doublet_list)
write_tsv(doublets, p("04_tables", "DOUBLET_RESULTS.tsv.gz"))

qc_obj$scDblFinder_score <- doublets$scDblFinder_score[match(colnames(qc_obj), doublets$barcode)]
qc_obj$scDblFinder_class <- doublets$scDblFinder_class[match(colnames(qc_obj), doublets$barcode)]

doublet_rate <- doublets[, .(
  qc_retained = .N,
  doublets = sum(scDblFinder_class == "doublet"),
  doublet_rate = mean(scDblFinder_class == "doublet")
), by = library_id]
p_doublet <- ggplot(doublet_rate, aes(library_id, doublet_rate, fill = library_id)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * doublet_rate, doublets)), vjust = -0.2) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "scDblFinder doublet rate") + theme_sc()
save_dual(p_doublet, "DOUBLET_RATE_BY_LIBRARY", 7, 5)

final_obj <- subset(qc_obj, subset = scDblFinder_class == "singlet")
initial_n <- nrow(qc)
qc_n <- ncol(qc_obj)
doublet_n <- sum(doublets$scDblFinder_class == "doublet")
final_n <- ncol(final_obj)
qc_flow <- data.table(
  step = c("INITIAL_CELLS", "QC_RETAINED", "DOUBLETS_REMOVED", "FINAL_CELLS"),
  n_cells = c(initial_n, qc_n, doublet_n, final_n),
  removed_from_previous = c(0L, initial_n - qc_n, doublet_n, 0L),
  definition = c("official processed matrix barcodes", "robust per-library QC pass",
                 "scDblFinder class=doublet, run independently per library",
                 "QC-passed scDblFinder singlets")
)
write_tsv(qc_flow, p("04_tables", "QC_FLOW.tsv"))

qc[, scDblFinder_score := doublets$scDblFinder_score[match(barcode, doublets$barcode)]]
qc[, scDblFinder_class := doublets$scDblFinder_class[match(barcode, doublets$barcode)]]
write_tsv(qc, p("04_tables", "QC_CELL_METADATA.tsv.gz"))

final_lib <- as.data.table(table(final_obj$library_id))
setnames(final_lib, c("library_id", "cell_number"))
write_tsv(final_lib, p("04_tables", "CELLS_PER_LIBRARY_FINAL.tsv"))

warnings <- rbindlist(lapply(seq_len(nrow(thresholds)), function(i) {
  th <- thresholds[i]
  retained <- sensitivity[scenario == "MAIN_ROBUST" & library_id == th$library_id, retained_fraction]
  dbl_rate <- doublet_rate[library_id == th$library_id, doublet_rate]
  msgs <- character()
  if (retained < 0.80) msgs <- c(msgs, "QC_RETAINED_BELOW_80_PERCENT")
  if (abs(th$max_percent_mt - 10) < 1e-8) msgs <- c(msgs, "MITO_THRESHOLD_LOWER_BOUND_10_PERCENT")
  if (abs(th$max_percent_mt - 30) < 1e-8) msgs <- c(msgs, "MITO_THRESHOLD_UPPER_BOUND_30_PERCENT")
  if (length(dbl_rate) && is.finite(dbl_rate) && dbl_rate > 0.10) {
    msgs <- c(msgs, sprintf("HIGH_SCDOUBLET_RATE_GT_10_PERCENT(%.2f%%)", 100 * dbl_rate))
  }
  if (!length(msgs)) msgs <- "NONE"
  data.table(library_id = th$library_id, warning = msgs)
}))
write_tsv(warnings, p("09_qc", "QC_WARNINGS.tsv"))

saveRDS(final_obj, p("06_objects", "02_qc_seurat.rds"), compress = FALSE)
rm(qc_obj, final_obj, sce, counts_qc)
gc(verbose = FALSE)
stage_message("Stage 02 complete: initial=", initial_n, ", QC=", qc_n,
              ", doublets=", doublet_n, ", final=", final_n)
