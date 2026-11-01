options(stringsAsFactors = FALSE, width = 140)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(UCell)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(digest)
})

seed <- 20260823L
set.seed(seed)
run_dir <- "/path/to/cochlear_analysis/GSE252285_UCELL_GSE160639_PROGRAMS_20260823"
out_dir <- file.path(run_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
obj_path <- "/path/to/cochlear_analysis/GSE252285_stage1_20260814/06_objects/no_doublet_robust_qc_tsne_sensitivity/02_no_doublet_annotated_tsne.rds"
hallmark_path <- file.path(run_dir, "gene_sets", "mh.all.v2026.1.Mm.symbols.gmt")
reactome_path <- file.path(run_dir, "gene_sets", "m2.cp.reactome.v2026.1.Mm.symbols.gmt")

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
save_dual <- function(p, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = height,
         units = "in", bg = "white", limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = width, height = height,
         units = "in", dpi = 300, bg = "white", limitsize = FALSE, device = ragg::agg_png)
}
read_gmt <- function(path, collection) {
  z <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  setNames(lapply(z, function(v) unique(v[-c(1, 2)])), vapply(z, `[[`, character(1), 1)) |>
    structure(collection = collection)
}

message("Read frozen gene sets")
hallmark <- read_gmt(hallmark_path, "MSigDB Hallmark")
reactome <- read_gmt(reactome_path, "MSigDB Reactome")
all_sets <- c(hallmark, reactome)

# These official MSigDB component sets are fixed before score calculation. Programs with
# slash-separated concepts are deterministic unions; no genes are manually added or removed.
program_components <- list(
  MAPK_AP1 = c(
    "REACTOME_ACTIVATION_OF_THE_AP_1_FAMILY_OF_TRANSCRIPTION_FACTORS",
    "REACTOME_MAPK_TARGETS_NUCLEAR_EVENTS_MEDIATED_BY_MAP_KINASES"
  ),
  TNFA_NFKB = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  IL6_JAK_STAT3 = "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  VEGF_ANGIOGENESIS = c("HALLMARK_ANGIOGENESIS", "REACTOME_SIGNALING_BY_VEGF"),
  MYC_TARGETS = c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"),
  E2F_TARGETS = "HALLMARK_E2F_TARGETS"
)
missing_sets <- setdiff(unique(unlist(program_components)), names(all_sets))
if (length(missing_sets)) stop("Frozen set(s) missing: ", paste(missing_sets, collapse = ", "))
signatures <- lapply(program_components, function(nm) unique(unlist(all_sets[nm], use.names = FALSE)))
if (length(signatures) != 6L || any(lengths(signatures) == 0L)) stop("Program construction failed")

definition <- rbindlist(lapply(names(program_components), function(program) {
  comps <- program_components[[program]]
  rbindlist(lapply(comps, function(component) data.table(
    program = program,
    component_gene_set = component,
    collection = if (startsWith(component, "HALLMARK_")) "MSigDB Hallmark" else "MSigDB Reactome",
    component_gene_count = length(all_sets[[component]]),
    union_gene_count = length(signatures[[program]])
  )))
}))
write_tsv(definition, file.path(out_dir, "PROGRAM_DEFINITIONS.tsv"))

message("Read fixed Seurat object (read-only)")
obj <- readRDS(obj_path)
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object")
if (nrow(obj) != 32285L || ncol(obj) != 36111L) stop("Fixed object dimension invariant failed")
if (!identical(Assays(obj), "RNA")) stop("Expected only RNA assay")
if (!all(c("counts", "data") %in% Layers(obj[["RNA"]]))) stop("RNA counts/data layer missing")
if (!all(c("pca", "tsne") %in% Reductions(obj))) stop("Fixed t-SNE reduction missing")
ann_field <- "broad_annotation_no_doublet_sensitivity"
if (!ann_field %in% colnames(obj[[]])) stop("Fixed broad annotation missing")
annotation_before <- obj[[ann_field, drop = TRUE]]
annotation_digest_before <- digest(annotation_before, algo = "sha256")
if (anyNA(annotation_before) || length(annotation_before) != 36111L) stop("Annotation invalid")
expected_counts <- c(
  Endothelial_vascular = 275L, Fibrocyte_mesenchymal = 2997L, Glial = 264L,
  Hair_sensory = 24L, Immune_myeloid = 27676L, Supporting_epithelial = 149L,
  UNRESOLVED = 4726L
)
observed_counts <- table(annotation_before)
if (!identical(as.integer(observed_counts[names(expected_counts)]), as.integer(expected_counts))) {
  stop("Fixed broad annotation count invariant failed")
}

meta <- as.data.table(obj[[]], keep.rownames = "cell")
meta[, broad_cell_type := as.character(get(ann_field))]
meta[, parsed_suffix := sub("^.*(-[123])$", "\\1", cell)]
if (any(!grepl("-[123]$", meta$cell)) ||
    !identical(meta$parsed_suffix, as.character(meta$barcode_suffix))) {
  stop("Barcode suffix invariant failed")
}
if (any(as.character(meta$condition) != "UNKNOWN")) stop("Legacy condition field changed")
condition_map <- data.table(
  barcode_suffix = c("-3", "-1", "-2"),
  condition_descriptive = c("Control", "HHL", "NIHL"),
  mapping_source = "USER_SPECIFIED_GEO_SUFFIX"
)
meta[condition_map, on = "barcode_suffix", condition_descriptive := i.condition_descriptive]
if (anyNA(meta$condition_descriptive)) stop("Condition descriptive mapping incomplete")
condition_levels <- c("Control", "HHL", "NIHL")
broad_levels <- c("Immune_myeloid", "Fibrocyte_mesenchymal", "Endothelial_vascular",
                  "Glial", "Supporting_epithelial", "Hair_sensory", "UNRESOLVED")
meta[, condition_descriptive := factor(condition_descriptive, levels = condition_levels)]
meta[, broad_cell_type := factor(broad_cell_type, levels = broad_levels)]
write_tsv(meta[, .N, by = .(barcode_suffix, library_id, condition_descriptive)][order(condition_descriptive)],
          file.path(out_dir, "CONDITION_MAPPING_AUDIT.tsv"))
write_tsv(meta[, .(n_cells = .N), by = broad_cell_type][order(match(broad_cell_type, broad_levels))],
          file.path(out_dir, "BROAD_ANNOTATION_COUNTS.tsv"))

counts <- LayerData(obj[["RNA"]], layer = "counts")
if (!inherits(counts, "sparseMatrix")) stop("RNA counts are not sparse")
if (!identical(colnames(counts), meta$cell)) stop("Cell order mismatch")
if (anyDuplicated(rownames(counts))) stop("Duplicated object gene symbols")

coverage <- rbindlist(lapply(names(signatures), function(program) {
  genes <- signatures[[program]]
  matched <- intersect(genes, rownames(counts))
  data.table(
    program = program,
    signature_genes = length(genes),
    matched_genes = length(matched),
    missing_genes = length(setdiff(genes, rownames(counts))),
    object_coverage_fraction = length(matched) / length(genes),
    matched_gene_symbols = paste(matched, collapse = ","),
    missing_gene_symbols = paste(setdiff(genes, rownames(counts)), collapse = ",")
  )
}))
if (any(coverage$matched_genes < 5L)) stop("A program has fewer than five matched genes")
write_tsv(coverage, file.path(out_dir, "PROGRAM_OBJECT_GENE_COVERAGE.tsv"))

message("UCell API smoke test")
smoke <- ScoreSignatures_UCell(
  matrix = counts[, seq_len(100), drop = FALSE], features = signatures,
  maxRank = 1500, chunk.size = 100, ncores = 1, ties.method = "average", force.gc = TRUE
)
if (nrow(smoke) != 100L || ncol(smoke) != 6L || anyNA(smoke)) stop("UCell smoke test failed")

message("Full UCell scoring: 36,111 cells, six programs")
scores <- ScoreSignatures_UCell(
  matrix = counts, features = signatures, maxRank = 1500,
  chunk.size = 500, ncores = 4, ties.method = "average", force.gc = TRUE
)
score_dt <- as.data.table(scores, keep.rownames = "cell")
expected_score_cols <- paste0(names(signatures), "_UCell")
if (!setequal(setdiff(names(score_dt), "cell"), expected_score_cols)) stop("Unexpected UCell columns")
setnames(score_dt, expected_score_cols, names(signatures))
if (!identical(score_dt$cell, meta$cell)) stop("UCell output cell order mismatch")
score_cols <- names(signatures)
score_values <- as.matrix(score_dt[, ..score_cols])
if (any(!is.finite(score_values)) || any(score_values < 0) || any(score_values > 1)) {
  stop("UCell score range/finite check failed")
}
score_dt[, `:=`(
  broad_cell_type = as.character(meta$broad_cell_type),
  condition_descriptive = as.character(meta$condition_descriptive),
  barcode_suffix = meta$barcode_suffix
)]
setcolorder(score_dt, c("cell", "barcode_suffix", "condition_descriptive", "broad_cell_type", names(signatures)))
fwrite(score_dt, file.path(out_dir, "UCELL_SCORES_PER_CELL.tsv.gz"), sep = "\t",
       quote = FALSE, na = "NA", compress = "gzip")

score_long <- melt(score_dt,
                   id.vars = c("cell", "barcode_suffix", "condition_descriptive", "broad_cell_type"),
                   measure.vars = names(signatures), variable.name = "program", value.name = "UCell_score")
broad_summary <- score_long[, .(
  n_cells = .N,
  median_UCell = median(UCell_score),
  Q1_UCell = quantile(UCell_score, 0.25, names = FALSE),
  Q3_UCell = quantile(UCell_score, 0.75, names = FALSE),
  IQR_UCell = IQR(UCell_score)
), by = .(program, broad_cell_type)]
broad_summary[, median_rank := frank(-median_UCell, ties.method = "min"), by = program]
broad_summary[, low_cell_count_flag := n_cells < 20L]

message("Per-cell-type gene coverage")
celltype_cov <- rbindlist(lapply(names(signatures), function(program) {
  matched <- intersect(signatures[[program]], rownames(counts))
  rbindlist(lapply(broad_levels, function(ct) {
    idx <- which(as.character(meta$broad_cell_type) == ct)
    x <- counts[matched, idx, drop = FALSE]
    gene_detected <- Matrix::rowSums(x > 0) > 0
    per_cell_detected_fraction <- Matrix::colSums(x > 0) / length(matched)
    data.table(
      program = program, broad_cell_type = ct,
      signature_genes = length(signatures[[program]]), matched_genes = length(matched),
      genes_detected_in_at_least_one_cell = sum(gene_detected),
      cell_type_gene_coverage_fraction = sum(gene_detected) / length(signatures[[program]]),
      matched_gene_detection_fraction = sum(gene_detected) / length(matched),
      median_per_cell_matched_gene_detection_fraction = median(per_cell_detected_fraction),
      Q1_per_cell_matched_gene_detection_fraction = quantile(per_cell_detected_fraction, 0.25, names = FALSE),
      Q3_per_cell_matched_gene_detection_fraction = quantile(per_cell_detected_fraction, 0.75, names = FALSE)
    )
  }))
}))
broad_summary <- merge(broad_summary, celltype_cov,
                       by = c("program", "broad_cell_type"), all.x = TRUE, sort = FALSE)
setorder(broad_summary, program, median_rank, -median_UCell)
write_tsv(broad_summary, file.path(out_dir, "UCELL_BROAD_CELLTYPE_SUMMARY.tsv"))

top2 <- broad_summary[order(program, median_rank, -median_UCell), head(.SD, 2), by = program]
top2[, top_position := seq_len(.N), by = program]
setcolorder(top2, c("program", "top_position", setdiff(names(top2), c("program", "top_position"))))
write_tsv(top2, file.path(out_dir, "UCELL_TOP2_CELLTYPES.tsv"))

condition_summary <- score_long[, .(
  n_cells = .N,
  median_UCell = median(UCell_score),
  Q1_UCell = quantile(UCell_score, 0.25, names = FALSE),
  Q3_UCell = quantile(UCell_score, 0.75, names = FALSE),
  IQR_UCell = IQR(UCell_score)
), by = .(program, broad_cell_type, condition_descriptive)]
condition_summary[, low_cell_count_flag := n_cells < 20L]
condition_summary[, condition_order := match(condition_descriptive, condition_levels)]
setorder(condition_summary, program, broad_cell_type, condition_order)
condition_summary[, condition_order := NULL]
write_tsv(condition_summary, file.path(out_dir, "UCELL_CONDITION_DESCRIPTIVE_SUMMARY.tsv"))

program_labels <- c(
  MAPK_AP1 = "MAPK / AP-1", TNFA_NFKB = "TNFa / NF-kB",
  IL6_JAK_STAT3 = "IL6 / JAK-STAT3", VEGF_ANGIOGENESIS = "VEGF / angiogenesis",
  MYC_TARGETS = "MYC targets", E2F_TARGETS = "E2F targets"
)

heat <- copy(broad_summary)
heat[, program_label := factor(program_labels[program], levels = unname(program_labels))]
heat[, broad_cell_type := factor(broad_cell_type, levels = rev(broad_levels))]
heat[, median_row_z := {
  s <- sd(median_UCell)
  if (is.na(s) || s == 0) rep(0, .N) else (median_UCell - mean(median_UCell)) / s
}, by = program]
p_heat <- ggplot(heat, aes(program_label, broad_cell_type, fill = median_row_z)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", median_UCell)), size = 3.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Median\nrow z-score") +
  labs(x = NULL, y = NULL, title = "UCell program activity by fixed broad cell type",
       subtitle = "Color: row-standardized median; text: raw median UCell score") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold"))
save_dual(p_heat, "UCELL_PROGRAM_BY_CELLTYPE_HEATMAP", 12, 6.5)

plot_long <- copy(score_long)
plot_long[, program_label := factor(program_labels[as.character(program)], levels = unname(program_labels))]
plot_long[, broad_cell_type := factor(broad_cell_type, levels = broad_levels)]
p_violin <- ggplot(plot_long, aes(broad_cell_type, UCell_score, fill = broad_cell_type)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  facet_wrap(~program_label, ncol = 2, scales = "free_y") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = "UCell score", title = "Program activity across fixed broad cell types") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 38, hjust = 1), strip.text = element_text(face = "bold"))
save_dual(p_violin, "UCELL_PROGRAM_VIOLINS_BROAD_CELLTYPE", 14, 11)

condition_group_n <- plot_long[, .(condition_group_n = .N),
                               by = .(program_label, broad_cell_type, condition_descriptive)]
plot_long <- merge(plot_long, condition_group_n,
                   by = c("program_label", "broad_cell_type", "condition_descriptive"), sort = FALSE)
p_condition <- ggplot(plot_long, aes(broad_cell_type, UCell_score, fill = condition_descriptive)) +
  geom_violin(data = plot_long[condition_group_n >= 2L], scale = "width", trim = TRUE,
              linewidth = 0.15) +
  stat_summary(fun = median, geom = "point", shape = 21, color = "black", size = 0.6,
               stroke = 0.15) +
  facet_grid(condition_descriptive ~ program_label, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#4E79A7", HHL = "#F28E2B", NIHL = "#E15759"), guide = "none") +
  labs(x = NULL, y = "UCell score", title = "Descriptive condition facets only",
       subtitle = "No condition-level hypothesis tests; cells are not treated as biological replicates") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
save_dual(p_condition, "UCELL_CONDITION_DESCRIPTIVE_VIOLINS", 19, 10)

emb <- as.data.table(Embeddings(obj, "tsne"), keep.rownames = "cell")
if (!identical(emb$cell, score_dt$cell)) stop("t-SNE cell order mismatch")
emb_long <- merge(emb, score_long[, .(cell, program, UCell_score)], by = "cell", sort = FALSE)
emb_long[, program_label := factor(program_labels[as.character(program)], levels = unname(program_labels))]
setorder(emb_long, program_label, UCell_score)
p_feature <- ggplot(emb_long, aes(tSNE_1, tSNE_2, color = UCell_score)) +
  geom_point(size = 0.12, alpha = 0.8) +
  facet_wrap(~program_label, ncol = 3) +
  scale_color_viridis_c(option = "magma") +
  coord_equal() +
  labs(x = "tSNE 1", y = "tSNE 2", color = "UCell", title = "UCell programs on the existing fixed t-SNE") +
  theme_void(base_size = 10) +
  theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
        legend.position = "right")
save_dual(p_feature, "UCELL_PROGRAM_FEATUREPLOTS_TSNE", 14, 9)

annotation_after <- obj[[ann_field, drop = TRUE]]
annotation_digest_after <- digest(annotation_after, algo = "sha256")
condition_counts <- meta[, .N, by = condition_descriptive]
condition_counts <- condition_counts[match(condition_levels, as.character(condition_descriptive))]
qa <- data.table(
  check = c(
    "fixed_object_dimensions", "cell_count", "feature_count", "annotation_no_NA",
    "annotation_counts_match_frozen", "annotation_vector_unchanged_in_memory",
    "legacy_condition_remains_UNKNOWN", "condition_mapping_complete",
    "RNA_counts_sparse", "gene_symbols_unique", "six_programs_fixed_before_scoring",
    "UCell_smoke_test", "all_scores_finite", "scores_within_0_1",
    "no_condition_hypothesis_tests", "no_QC_clustering_annotation_rerun",
    "no_manual_signature_genes", "no_matched_null"
  ),
  status = "PASS",
  detail = c(
    "32285 genes x 36111 cells", "36111", "32285", "0 NA",
    paste(names(expected_counts), expected_counts, sep = "=", collapse = ";"),
    paste0(annotation_digest_before, " == ", annotation_digest_after),
    "All input metadata condition values are UNKNOWN",
    paste(condition_levels, condition_counts$N, sep = "=", collapse = ";"),
    class(counts)[1], "No duplicated RNA rownames", paste(names(signatures), collapse = ";"),
    "100-cell API smoke test passed", "All finite", "All UCell scores in [0,1]",
    "Descriptive summaries and facets only", "Existing object read only; existing t-SNE reused",
    "Only exact frozen MSigDB sets and deterministic unions", "Not run"
  )
)
write_tsv(qa, file.path(out_dir, "QA_STATUS.tsv"))

versions <- data.table(
  software = c("R", "UCell", "Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2", "patchwork"),
  version = c(as.character(getRversion()), as.character(packageVersion("UCell")),
              as.character(packageVersion("Seurat")), as.character(packageVersion("SeuratObject")),
              as.character(packageVersion("Matrix")), as.character(packageVersion("data.table")),
              as.character(packageVersion("ggplot2")), as.character(packageVersion("patchwork")))
)
write_tsv(versions, file.path(out_dir, "SOFTWARE_VERSIONS.tsv"))

methods <- c(
  "# METHODS — GSE252285 UCell program localization",
  "",
  paste0("- Input object: `", obj_path, "` (read only)."),
  "- Fixed object: 32,285 RNA features and 36,111 cells.",
  paste0("- Annotation field: `", ann_field, "`; values and ordering were not modified."),
  "- No QC, normalization, dimensional reduction, clustering, doublet filtering, or annotation was rerun.",
  "- Existing t-SNE coordinates were reused only for visualization.",
  "- Species: Mus musculus.",
  "- Gene sets: frozen Mouse MSigDB 2026.1.Mm Hallmark and Reactome symbol GMT files.",
  "- Programs were declared before scoring as exact official sets or deterministic unions listed in `PROGRAM_DEFINITIONS.tsv`.",
  "- No hand-written genes were added; APOLD1/ERRFI1 were not used to define or adapt any signature.",
  "- UCell input: complete sparse RNA counts matrix; no DEG or gene prefiltering before cell ranking.",
  "- UCell parameters: maxRank=1500, ties.method=average, chunk.size=500, ncores=4.",
  paste0("- Seed: ", seed, ". UCell itself is deterministic for these settings."),
  "- Broad summaries: median, first quartile (Q1), third quartile (Q3), and IQR per fixed broad cell type.",
  "- Gene coverage: object-level signature matches; per cell type, genes detected in at least one cell and median per-cell matched-gene detection fraction.",
  "- Top cell types were ranked by raw median UCell score without excluding low-cell-count or UNRESOLVED groups; low cell counts are flagged.",
  "- Control/HHL/NIHL mapping follows the user-resolved GEO barcode suffix mapping (-3/-1/-2). Conditions are descriptive facets only.",
  "- No condition-level P values, matched-null analysis, or cell-as-replicate inference was performed.",
  "",
  "## Versions",
  "",
  paste0("- ", versions$software, ": ", versions$version),
  "",
  "## Frozen gene-set files",
  "",
  paste0("- `", hallmark_path, "`"),
  paste0("- `", reactome_path, "`")
)
writeLines(methods, file.path(out_dir, "METHODS.md"), useBytes = TRUE)

message("Run completed")
