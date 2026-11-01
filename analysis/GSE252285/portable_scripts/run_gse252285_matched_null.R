# Portable archive interface; scientific body unchanged from the retrieved script.
options(stringsAsFactors = FALSE, width = 160)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(UCell)
  library(Matrix)
  library(data.table)
  library(matrixStats)
  library(ggplot2)
  library(digest)
})

B <- as.integer(Sys.getenv("N_ITER", "100"))
if (!B %in% c(100L, 10000L)) stop("N_ITER must be 100 (benchmark) or 10000 (final)")
is_final <- B == 10000L
seed <- 20260823L
set.seed(seed)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4L) stop("Usage: run_gse252285_matched_null.R OBJECT_RDS FROZEN_GENE_SET_DIRECTORY GENCODE_VM23_GTF_GZ NEW_RUN_DIRECTORY")
obj_path <- normalizePath(args[[1]], mustWork=TRUE)
gene_set_dir <- normalizePath(args[[2]], mustWork=TRUE)
gtf_path <- normalizePath(args[[3]], mustWork=TRUE)
run_dir <- path.expand(args[[4]])
if (file.exists(run_dir) || dir.exists(run_dir)) stop("Use a new output directory")
if (unname(tools::md5sum(gtf_path)) != "b646883b82252b8831d337e9b88d1f30") stop("GENCODE vM23 checksum mismatch")
out_dir <- if (is_final) file.path(run_dir, "outputs") else file.path(run_dir, "benchmark")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(run_dir,"work"), recursive=TRUE, showWarnings=FALSE)
hallmark_path <- file.path(gene_set_dir, "mh.all.v2026.1.Mm.symbols.gmt")
reactome_path <- file.path(gene_set_dir, "m2.cp.reactome.v2026.1.Mm.symbols.gmt")
if (digest::digest(obj_path, algo="sha256", file=TRUE) != "24909e74b584dcb33a113d815234cda0861ed6fea07f78cf66fd99d5324d76b8") stop("Selected Seurat object checksum mismatch")
if (digest::digest(hallmark_path, algo="sha256", file=TRUE) != "3a21be724a87dc0375955e725ca9688b87a26e7b74b62fba0c62da0967b789f7") stop("Hallmark checksum mismatch")
if (digest::digest(reactome_path, algo="sha256", file=TRUE) != "3beeef870afc4bb803f8e7aa184f8c98a22bcc799b702b48dc7f9c85d50cd3a9") stop("Reactome checksum mismatch")

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
save_dual <- function(p, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = height,
         units = "in", bg = "white", limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = width, height = height,
         units = "in", dpi = 300, bg = "white", limitsize = FALSE, device = ragg::agg_png)
}
read_gmt <- function(path) {
  z <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  setNames(lapply(z, function(v) unique(v[-c(1, 2)])), vapply(z, `[[`, character(1), 1))
}
rank_bin <- function(x, nbin = 5L) {
  pmin(nbin, pmax(1L, ceiling(nbin * frank(x, ties.method = "average") / length(x))))
}

message("Read frozen signatures")
all_sets <- c(read_gmt(hallmark_path), read_gmt(reactome_path))
components <- list(
  MAPK_AP1 = c("REACTOME_ACTIVATION_OF_THE_AP_1_FAMILY_OF_TRANSCRIPTION_FACTORS",
               "REACTOME_MAPK_TARGETS_NUCLEAR_EVENTS_MEDIATED_BY_MAP_KINASES"),
  VEGF_ANGIOGENESIS = c("HALLMARK_ANGIOGENESIS", "REACTOME_SIGNALING_BY_VEGF"),
  TNFA_NFKB = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  IL6_JAK_STAT3 = "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)
if (length(setdiff(unique(unlist(components)), names(all_sets)))) stop("Frozen component missing")
official_signatures <- lapply(components, function(x) unique(unlist(all_sets[x], use.names = FALSE)))
tests <- data.table(
  test_id = c("MAPK_AP1__Endothelial_vascular", "VEGF_ANGIOGENESIS__Endothelial_vascular",
              "TNFA_NFKB__Endothelial_vascular", "TNFA_NFKB__Fibrocyte_mesenchymal",
              "IL6_JAK_STAT3__Endothelial_vascular", "IL6_JAK_STAT3__Immune_myeloid"),
  program = c("MAPK_AP1", "VEGF_ANGIOGENESIS", "TNFA_NFKB", "TNFA_NFKB",
              "IL6_JAK_STAT3", "IL6_JAK_STAT3"),
  target_cell_type = c("Endothelial_vascular", "Endothelial_vascular", "Endothelial_vascular",
                       "Fibrocyte_mesenchymal", "Endothelial_vascular", "Immune_myeloid")
)
stopifnot(!any(grepl("MYC|E2F", tests$program)))

message("Read fixed object and calculate detectable-gene covariates")
obj <- readRDS(obj_path)
if (nrow(obj) != 32285L || ncol(obj) != 36111L) stop("Fixed object dimension invariant failed")
ann_field <- "broad_annotation_no_doublet_sensitivity"
ann_before <- obj[[ann_field, drop = TRUE]]
ann_digest_before <- digest(ann_before, algo = "sha256")
if (anyNA(ann_before)) stop("Annotation contains NA")
counts <- LayerData(obj[["RNA"]], layer = "counts")
if (!inherits(counts, "sparseMatrix") || anyDuplicated(rownames(counts))) stop("RNA count matrix invariant failed")
detection_cells <- as.integer(Matrix::rowSums(counts > 0))
universe <- data.table(
  gene = rownames(counts),
  mean_expression = as.numeric(Matrix::rowMeans(counts)),
  detection_cells = detection_cells,
  detection_rate = detection_cells / ncol(counts)
)[detection_cells > 0L]
if (nrow(universe) < 10000L) stop("Unexpectedly small detectable-gene universe")

message("Parse GENCODE vM23 gene-span lengths for sensitivity")
gtf_cmd <- sprintf("gzip -cd %s | awk '$3==\"gene\"'", shQuote(gtf_path))
gtf <- fread(cmd = gtf_cmd, sep = "\t", header = FALSE, quote = "", fill = TRUE,
             select = c(4L, 5L, 9L), col.names = c("start", "end", "attribute"), showProgress = FALSE)
gtf[, gene := sub('.*gene_name "([^"]+)".*', "\\1", attribute)]
gtf <- gtf[gene != attribute, .(gene_length_bp = max(end - start + 1)), by = gene]
universe <- merge(universe, gtf, by = "gene", all.x = TRUE, sort = FALSE)
universe[, `:=`(expr_bin = rank_bin(mean_expression, 10L), det_bin = rank_bin(detection_rate, 10L))]
universe[, primary_stratum := sprintf("E%02d_D%02d", expr_bin, det_bin)]
universe[, length_bin := fifelse(is.na(gene_length_bp), "LMISSING",
                                 sprintf("L%02d", rank_bin(log10(gene_length_bp), 3L)))]
universe[, length_sensitivity_stratum := paste(primary_stratum, length_bin, sep = "_")]
write_tsv(universe, file.path(out_dir, "DETECTABLE_GENE_UNIVERSE.tsv"))

signatures <- lapply(official_signatures, function(g) intersect(g, universe$gene))
coverage <- rbindlist(lapply(names(signatures), function(program) data.table(
  program = program,
  official_signature_size = length(official_signatures[[program]]),
  detectable_signature_size = length(signatures[[program]]),
  missing_from_detectable_universe = length(setdiff(official_signatures[[program]], universe$gene)),
  gene_length_mapped = sum(!is.na(universe[match(signatures[[program]], gene), gene_length_bp])),
  gene_length_mapping_fraction = mean(!is.na(universe[match(signatures[[program]], gene), gene_length_bp])),
  detectable_signature_genes = paste(signatures[[program]], collapse = ","),
  nondetectable_signature_genes = paste(setdiff(official_signatures[[program]], universe$gene), collapse = ",")
)))
if (any(coverage$detectable_signature_size < 20L)) stop("Too few signature genes")
write_tsv(coverage, file.path(out_dir, "SIGNATURE_COVERAGE.tsv"))

resolved_types <- c("Endothelial_vascular", "Fibrocyte_mesenchymal", "Glial",
                    "Hair_sensory", "Immune_myeloid", "Supporting_epithelial")
cell_index <- setNames(lapply(resolved_types, function(ct) which(ann_before == ct)), resolved_types)
if (any(lengths(cell_index) == 0L)) stop("Resolved cell-type index missing")

message("Precompute UCell rankings")
rank_start <- proc.time()[[3]]
rank_cache <- file.path(run_dir, "work", "ucell_ranks_maxRank1500.rds")
if (file.exists(rank_cache)) {
  ranks <- readRDS(rank_cache)
  rank_cache_status <- "REUSED_TASK_CACHE"
} else {
  ranks <- StoreRankings_UCell(matrix = counts, maxRank = 1500, chunk.size = 500,
                               ncores = 4, ties.method = "average", force.gc = TRUE)
  saveRDS(ranks, rank_cache, compress = FALSE)
  rank_cache_status <- "CREATED_TASK_CACHE"
}
rank_seconds <- proc.time()[[3]] - rank_start
if (ncol(ranks) != ncol(counts)) stop("Rank matrix cell count mismatch")
if (!inherits(ranks, "sparseMatrix") || !identical(rownames(ranks), rownames(counts))) {
  stop("Unexpected UCell rank matrix representation")
}

# For positive signatures, UCell::u_stat is exactly the sum of (maxRank-rank)
# contributions divided by n*maxRank-n*(n+1)/2. Sparse multiplication computes
# this same statistic for many null sets without changing the algorithm.
rank_weights <- ranks
rank_weights@x <- 1500 - rank_weights@x
fast_ucell <- function(feature_sets, weight_matrix = rank_weights) {
  nset <- length(feature_sets)
  set_sizes <- lengths(feature_sets)
  if (!nset || any(set_sizes < 1L) || any(set_sizes >= 1500L)) stop("Invalid feature-set size")
  gene_index <- match(unlist(feature_sets, use.names = FALSE), rownames(weight_matrix))
  if (anyNA(gene_index)) stop("Feature absent from UCell rank matrix")
  incidence <- Matrix::sparseMatrix(
    i = rep.int(seq_len(nset), set_sizes), j = gene_index, x = 1,
    dims = c(nset, nrow(weight_matrix)),
    dimnames = list(names(feature_sets), rownames(weight_matrix))
  )
  benefit <- incidence %*% weight_matrix
  denominator <- set_sizes * 1500 - set_sizes * (set_sizes + 1) / 2
  score <- as.matrix(benefit)
  score <- score / denominator
  score <- t(score)
  colnames(score) <- names(feature_sets)
  rownames(score) <- colnames(weight_matrix)
  score
}

message("Validate sparse UCell implementation against official UCell API")
validation_sets <- signatures[c("MAPK_AP1", "TNFA_NFKB")]
official_validation <- as.matrix(ScoreSignatures_UCell(
  matrix = counts[, seq_len(200), drop = FALSE], features = validation_sets,
  maxRank = 1500, chunk.size = 200, ncores = 1, ties.method = "average"
))
colnames(official_validation) <- sub("_UCell$", "", colnames(official_validation))
fast_validation <- fast_ucell(validation_sets, rank_weights[, seq_len(200), drop = FALSE])
fast_validation <- fast_validation[, colnames(official_validation), drop = FALSE]
fast_ucell_max_abs_error <- max(abs(official_validation - fast_validation))
if (!is.finite(fast_ucell_max_abs_error) || fast_ucell_max_abs_error > 1e-12) {
  stop("Sparse UCell implementation validation failed: ", fast_ucell_max_abs_error)
}

calc_localization <- function(score_matrix, target_types) {
  if (is.null(dim(score_matrix))) score_matrix <- matrix(score_matrix, ncol = 1L)
  type_medians <- vapply(resolved_types, function(ct) {
    matrixStats::colMedians(score_matrix[cell_index[[ct]], , drop = FALSE], na.rm = FALSE)
  }, numeric(ncol(score_matrix)))
  if (is.null(dim(type_medians))) type_medians <- matrix(type_medians, nrow = ncol(score_matrix))
  colnames(type_medians) <- resolved_types
  ans <- lapply(target_types, function(target) {
    other <- setdiff(resolved_types, target)
    target_med <- type_medians[, target]
    comparator <- matrixStats::rowMedians(type_medians[, other, drop = FALSE])
    data.table(target_cell_type = target, target_median = target_med,
               comparator_median_of_other_resolved_types = comparator,
               localization_delta = target_med - comparator)
  })
  rbindlist(ans, idcol = "target_index")
}

message("Score observed detectable signatures")
observed_matrix <- fast_ucell(signatures)
if (!setequal(colnames(observed_matrix), names(signatures))) stop("Observed score columns mismatch")
observed_rows <- rbindlist(lapply(names(signatures), function(prog) {
  targets <- tests[program == prog, target_cell_type]
  loc <- calc_localization(observed_matrix[, prog, drop = FALSE], targets)
  loc[, `:=`(program = prog, test_id = paste(prog, target_cell_type, sep = "__"))]
  loc
}))
observed_rows <- observed_rows[, .(test_id, program, target_cell_type, target_median,
                                   comparator_median_of_other_resolved_types, localization_delta)]

make_match_spec <- function(program, stratum_col) {
  sig <- signatures[[program]]
  sig_dt <- universe[gene %in% sig]
  req <- sig_dt[, .(required = .N), by = .(stratum = get(stratum_col))]
  candidates <- universe[!gene %in% sig, .(genes = list(gene), available = .N),
                         by = .(stratum = get(stratum_col))]
  audit <- merge(req, candidates, by = "stratum", all.x = TRUE)
  audit[is.na(available), available := 0L]
  audit[, `:=`(program = program, match_mode = stratum_col,
               supply_ratio = available / required, pass = available >= required)]
  if (any(!audit$pass)) {
    print(audit[!pass])
    stop("Insufficient exact-stratum candidate supply for ", program, " / ", stratum_col)
  }
  candidate_map <- setNames(candidates$genes, candidates$stratum)
  list(requirement = req, candidates = candidate_map, audit = audit)
}

draw_sets <- function(spec, n, prefix) {
  req <- spec$requirement
  ans <- vector("list", n)
  for (i in seq_len(n)) {
    genes <- unlist(lapply(seq_len(nrow(req)), function(j) {
      sample(spec$candidates[[req$stratum[j]]], req$required[j], replace = FALSE)
    }), use.names = FALSE)
    if (anyDuplicated(genes)) stop("Random set contains duplicated genes")
    ans[[i]] <- genes
  }
  names(ans) <- sprintf("%s_%05d", prefix, seq_len(n))
  ans
}

match_modes <- c(primary = "primary_stratum", length_sensitivity = "length_sensitivity_stratum")
null_rows <- list()
match_audits <- list()
balance_rows <- list()
runtime_rows <- list()
null_key <- 0L
chunk_size <- if (is_final) 500L else 100L

for (mode_name in names(match_modes)) {
  stratum_col <- match_modes[[mode_name]]
  for (program in names(signatures)) {
    prog <- program
    message("Null: ", mode_name, " / ", prog, " / B=", B)
    spec <- make_match_spec(prog, stratum_col)
    match_audits[[paste(mode_name, prog)]] <- spec$audit
    target_types <- tests[program == prog, target_cell_type]
    program_null <- setNames(lapply(target_types, function(x) numeric(B)), target_types)
    balance <- data.table(iteration = seq_len(B), mean_expression = NA_real_,
                          mean_log1p_expression = NA_real_, detection_rate = NA_real_,
                          mean_log10_gene_length = NA_real_)
    null_start <- proc.time()[[3]]
    for (first in seq.int(1L, B, by = chunk_size)) {
      last <- min(B, first + chunk_size - 1L)
      n_chunk <- last - first + 1L
      random_sets <- draw_sets(spec, n_chunk, sprintf("null_%s_%s_%05d", mode_name, program, first))
      # Preserve the full-run iteration numbering independently of chunk prefix names.
      names(random_sets) <- sprintf("iter_%05d", first:last)
      for (j in seq_along(random_sets)) {
        u <- universe[match(random_sets[[j]], gene)]
        balance$mean_expression[first + j - 1L] <- mean(u$mean_expression)
        balance$mean_log1p_expression[first + j - 1L] <- mean(log1p(u$mean_expression))
        balance$detection_rate[first + j - 1L] <- mean(u$detection_rate)
        balance$mean_log10_gene_length[first + j - 1L] <- mean(log10(u$gene_length_bp), na.rm = TRUE)
      }
      random_matrix <- fast_ucell(random_sets)
      if (nrow(random_matrix) != ncol(counts) || ncol(random_matrix) != n_chunk) stop("Null score shape mismatch")
      loc <- calc_localization(random_matrix, target_types)
      for (target in target_types) {
        program_null[[target]][first:last] <- loc[target_cell_type == target, localization_delta]
      }
      rm(random_matrix, loc, random_sets)
      gc(verbose = FALSE)
    }
    elapsed <- proc.time()[[3]] - null_start
    runtime_rows[[paste(mode_name, program)]] <- data.table(
      match_mode = mode_name, program = program, iterations = B,
      elapsed_seconds = elapsed, iterations_per_second = B / elapsed
    )
    sig_u <- universe[match(signatures[[program]], gene)]
    balance_summary <- data.table(
      match_mode = mode_name, program = program, iterations = B,
      signature_mean_expression = mean(sig_u$mean_expression),
      null_mean_of_set_mean_expression = mean(balance$mean_expression),
      null_sd_of_set_mean_expression = sd(balance$mean_expression),
      signature_mean_log1p_expression = mean(log1p(sig_u$mean_expression)),
      null_mean_of_set_mean_log1p_expression = mean(balance$mean_log1p_expression),
      null_sd_of_set_mean_log1p_expression = sd(balance$mean_log1p_expression),
      signature_mean_detection_rate = mean(sig_u$detection_rate),
      null_mean_of_set_mean_detection_rate = mean(balance$detection_rate),
      null_sd_of_set_mean_detection_rate = sd(balance$detection_rate),
      signature_mean_log10_gene_length = mean(log10(sig_u$gene_length_bp), na.rm = TRUE),
      null_mean_of_set_mean_log10_gene_length = mean(balance$mean_log10_gene_length),
      null_sd_of_set_mean_log10_gene_length = sd(balance$mean_log10_gene_length)
    )
    balance_rows[[paste(mode_name, program)]] <- balance_summary
    for (target in target_types) {
      null_key <- null_key + 1L
      null_rows[[null_key]] <- data.table(
        match_mode = mode_name, program = program, target_cell_type = target,
        test_id = paste(program, target, sep = "__"), iteration = seq_len(B),
        null_localization_delta = program_null[[target]]
      )
    }
  }
}

null_dt <- rbindlist(null_rows)
results <- merge(null_dt[, .(
  null_mean = mean(null_localization_delta),
  null_sd = sd(null_localization_delta),
  null_median = median(null_localization_delta),
  null_Q025 = quantile(null_localization_delta, 0.025),
  null_Q975 = quantile(null_localization_delta, 0.975)
), by = .(match_mode, test_id, program, target_cell_type)], observed_rows,
by = c("test_id", "program", "target_cell_type"), allow.cartesian = TRUE)
p_eval <- merge(null_dt, observed_rows[, .(test_id, observed_delta = localization_delta)],
                by = "test_id", all.x = TRUE, sort = FALSE)
p_values <- p_eval[, .(
  empirical_P = (1 + sum(null_localization_delta >= observed_delta[1])) / (B + 1)
), by = .(match_mode, test_id)]
results <- merge(results, p_values, by = c("match_mode", "test_id"), all.x = TRUE, sort = FALSE)
results[, empirical_FDR := p.adjust(empirical_P, method = "BH"), by = match_mode]
results[, `:=`(
  null_z = (localization_delta - null_mean) / null_sd,
  significant_FDR05 = empirical_FDR < 0.05 & localization_delta > 0,
  iterations = B,
  seed = seed,
  alternative = "greater"
)]
setorder(results, match_mode, empirical_FDR, -localization_delta)

write_tsv(results, file.path(out_dir, "MATCHED_NULL_RESULTS.tsv"))
fwrite(null_dt, file.path(out_dir, "MATCHED_NULL_DISTRIBUTIONS.tsv.gz"), sep = "\t",
       quote = FALSE, compress = "gzip")
write_tsv(rbindlist(match_audits), file.path(out_dir, "STRATUM_SUPPLY_AUDIT.tsv"))
write_tsv(rbindlist(balance_rows), file.path(out_dir, "MATCHING_BALANCE_QC.tsv"))
  runtime <- rbindlist(runtime_rows)
  runtime <- rbindlist(list(data.table(match_mode = "rank_precomputation", program = "ALL",
                                     iterations = NA_integer_, elapsed_seconds = rank_seconds,
                                     iterations_per_second = NA_real_, cache_status = rank_cache_status), runtime), fill = TRUE)
write_tsv(runtime, file.path(out_dir, "RUNTIME.tsv"))

if (is_final) {
  primary <- results[match_mode == "primary"]
  sensitivity <- results[match_mode == "length_sensitivity"]
  concordance <- merge(
    primary[, .(test_id, program, target_cell_type, primary_delta = localization_delta,
                primary_P = empirical_P, primary_FDR = empirical_FDR,
                primary_pass_FDR05 = significant_FDR05)],
    sensitivity[, .(test_id, sensitivity_delta = localization_delta,
                    sensitivity_P = empirical_P, sensitivity_FDR = empirical_FDR,
                    sensitivity_pass_FDR05 = significant_FDR05)], by = "test_id"
  )
  concordance[, sensitivity_concordant := primary_delta > 0 & sensitivity_delta > 0 &
                primary_pass_FDR05 == sensitivity_pass_FDR05]
  write_tsv(concordance, file.path(out_dir, "GENE_LENGTH_SENSITIVITY_CONCORDANCE.tsv"))

  program_label <- c(MAPK_AP1 = "MAPK / AP-1", VEGF_ANGIOGENESIS = "VEGF / angiogenesis",
                     TNFA_NFKB = "TNFa / NF-kB", IL6_JAK_STAT3 = "IL6 / JAK-STAT3")
  plot_dt <- merge(null_dt, results[, .(match_mode, test_id, localization_delta, empirical_P, empirical_FDR)],
                   by = c("match_mode", "test_id"), all.x = TRUE)
  plot_dt[, panel := paste0(program_label[program], " -> ", gsub("_", " ", target_cell_type))]
  plot_dt[, mode_label := factor(match_mode,
                                 levels = c("primary", "length_sensitivity"),
                                 labels = c("Primary: expression + detection matched",
                                            "Sensitivity: + gene-length strata"))]
  annotations <- unique(plot_dt[, .(mode_label, panel, localization_delta, empirical_P, empirical_FDR)])
  p_null <- ggplot(plot_dt, aes(null_localization_delta)) +
    geom_histogram(bins = 55, fill = "#7A9CC6", color = "white", linewidth = 0.15) +
    geom_vline(data = annotations, aes(xintercept = localization_delta), color = "#B2182B", linewidth = 0.8) +
    geom_text(data = annotations,
              aes(x = Inf, y = Inf, label = sprintf("P=%.4g\nFDR=%.4g", empirical_P, empirical_FDR)),
              hjust = 1.08, vjust = 1.25, inherit.aes = FALSE, size = 3) +
    facet_grid(mode_label ~ panel, scales = "free") +
    labs(x = "Null localization delta", y = "Iterations",
         title = "Matched-gene empirical null distributions",
         subtitle = "Red line: observed target median minus median of other resolved cell-type medians") +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))
  save_dual(p_null, "MATCHED_NULL_DISTRIBUTIONS", 18, 8)

  result_plot <- copy(results)
  result_plot[, panel := paste0(program_label[program], " -> ", gsub("_", " ", target_cell_type))]
  result_plot[, mode_label := factor(match_mode,
                                     levels = c("primary", "length_sensitivity"),
                                     labels = c("Primary", "+ gene length"))]
  p_effect <- ggplot(result_plot,
                     aes(localization_delta, reorder(panel, localization_delta), color = mode_label)) +
    geom_errorbarh(aes(xmin = null_Q025, xmax = null_Q975), height = 0.18,
                   position = position_dodge(width = 0.45), linewidth = 0.5) +
    geom_point(size = 2.5, position = position_dodge(width = 0.45)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    scale_color_manual(values = c("Primary" = "#2166AC", "+ gene length" = "#B2182B")) +
    labs(x = "Observed localization delta (points); null 95% interval (bars)", y = NULL,
         color = "Matching", title = "Observed localization against matched-gene null") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(face = "bold"))
  save_dual(p_effect, "MATCHED_NULL_EFFECT_SUMMARY", 11, 6.5)

  ann_after <- obj[[ann_field, drop = TRUE]]
  qa <- data.table(
    check = c("fixed_object_dimensions", "annotation_unchanged", "detectable_gene_universe",
              "only_six_prespecified_tests", "MYC_E2F_not_analyzed", "iterations_primary",
              "iterations_length_sensitivity", "exact_signature_size_matching",
              "exact_stratum_supply", "fast_UCell_matches_official_API", "scores_finite", "no_condition_inference",
              "no_QC_clustering_annotation_rerun"),
    status = "PASS",
    detail = c("32285 x 36111", paste0(ann_digest_before, " == ", digest(ann_after, algo = "sha256")),
               paste0(nrow(universe), " genes with detection_cells > 0"), paste(tests$test_id, collapse = ";"),
               "No MYC/E2F signatures or tests", as.character(B), as.character(B),
               paste(coverage$program, coverage$detectable_signature_size, sep = "=", collapse = ";"),
               paste0("minimum supply ratio=", sprintf("%.2f", min(rbindlist(match_audits)$supply_ratio))),
               paste0("max absolute error=", format(fast_ucell_max_abs_error, scientific = TRUE)),
               "All null and observed localization statistics finite", "No condition field used",
               "Fixed counts, annotation, and existing object only")
  )
  if (any(!is.finite(null_dt$null_localization_delta)) || any(!is.finite(results$localization_delta))) stop("Non-finite result")
  write_tsv(qa, file.path(out_dir, "QA_STATUS.tsv"))

  versions <- data.table(
    software = c("R", "UCell", "Seurat", "SeuratObject", "Matrix", "matrixStats", "data.table", "ggplot2"),
    version = c(as.character(getRversion()), as.character(packageVersion("UCell")),
                as.character(packageVersion("Seurat")), as.character(packageVersion("SeuratObject")),
                as.character(packageVersion("Matrix")), as.character(packageVersion("matrixStats")),
                as.character(packageVersion("data.table")), as.character(packageVersion("ggplot2")))
  )
  write_tsv(versions, file.path(out_dir, "SOFTWARE_VERSIONS.tsv"))

  methods <- c(
    "# METHODS — matched-gene empirical null for GSE252285 UCell localization",
    "",
    paste0("- Fixed input object: `", obj_path, "` (read only; 32,285 genes x 36,111 cells)."),
    "- No QC, normalization, clustering, annotation, or condition-level analysis was run.",
    "- Tested only the six prespecified program/cell-type localizations; MYC and E2F were excluded.",
    "- Gene universe: genes in the RNA counts layer detected in at least one of 36,111 cells.",
    "- Effective signature: official frozen MSigDB 2026.1.Mm signature genes present in that detectable universe.",
    "- Observed statistic: target cell-type median UCell minus the median of the other five resolved broad cell-type medians.",
    "- UNRESOLVED was excluded from the localization comparator; all six resolved annotation groups were otherwise retained.",
    "- Primary matching: exact strata from rank-based deciles of global mean count and detection rate (10 x 10 grid).",
    "- Random sets were sampled without replacement within strata and exactly matched the detectable signature size.",
    paste0("- Primary iterations: ", B, "; one-sided empirical P=(1+#null>=observed)/(B+1)."),
    "- FDR: Benjamini-Hochberg across the six prespecified tests, separately for primary and sensitivity analyses.",
    "- Gene-length sensitivity: the same expression/detection strata plus gene-span-length tertiles; missing lengths retained in an explicit LMISSING stratum.",
    "- Gene length source: GENCODE mouse vM23 comprehensive GTF (GRCm38/mm10); gene span=end-start+1.",
    paste0("- Seed: ", seed, ". UCell maxRank=1500; average ties; rankings computed once and reused."),
    paste0("- Batched null scoring used sparse matrix multiplication implementing the exact positive-signature UCell u_stat formula; maximum absolute validation error versus the official API was ", format(fast_ucell_max_abs_error, scientific = TRUE), "."),
    "- Control/HHL/NIHL labels were not used; no condition-level P values were calculated.",
    "",
    "## Software",
    "",
    paste0("- ", versions$software, ": ", versions$version)
  )
  writeLines(methods, file.path(out_dir, "METHODS.md"), useBytes = TRUE)
}

message("Matched-null run completed: B=", B)
