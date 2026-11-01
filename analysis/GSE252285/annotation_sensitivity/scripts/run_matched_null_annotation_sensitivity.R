# Optional existing fallback library is appended after the original runtime libraries.
# Used only to resolve an already installed data-only package; nothing is installed.
fallback_library <- Sys.getenv("EXISTING_DATA_LIBRARY", "")
if(nzchar(fallback_library)) .libPaths(c(.libPaths(),normalizePath(fallback_library,mustWork=TRUE)))
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
if (length(args) != 6L) stop("Usage: script.R OBJECT_RDS FROZEN_GENE_SET_DIRECTORY GENCODE_VM23_GTF_GZ NEW_RUN_DIRECTORY EXISTING_RANK_CACHE_RDS ORIGINAL_NULL_OUTPUTS")
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
cluster <- as.character(obj$seurat_clusters_no_doublet)
stopifnot(sum(cluster=="19")==319L, all(ann_before[cluster=="19"]=="UNRESOLVED"))
ann_revised <- ann_before
ann_revised[cluster=="19"] <- "Fibrocyte_mesenchymal"
original_cell_index <- setNames(lapply(resolved_types, function(ct) which(ann_before == ct)), resolved_types)
cell_index <- setNames(lapply(resolved_types, function(ct) which(ann_revised == ct)), resolved_types)
if (any(lengths(cell_index) == 0L)) stop("Resolved cell-type index missing")

message("Precompute UCell rankings")
rank_start <- proc.time()[[3]]
rank_cache <- normalizePath(args[[5]],mustWork=TRUE)
ranks <- readRDS(rank_cache)
rank_cache_status <- "REUSED_READ_ONLY_ORIGINAL_CACHE"
stopifnot(identical(colnames(ranks),colnames(counts)))
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
  original_medians <- type_medians
  original_medians[,"Fibrocyte_mesenchymal"] <- matrixStats::colMedians(
    score_matrix[original_cell_index[["Fibrocyte_mesenchymal"]],,drop=FALSE],na.rm=FALSE)
  summarize_medians <- function(medians) {
    rbindlist(lapply(target_types, function(target) {
      other <- setdiff(resolved_types, target)
      target_med <- medians[, target]
      comparator <- matrixStats::rowMedians(medians[, other, drop = FALSE])
      data.table(target_cell_type = target, target_median = target_med,
        comparator_median_of_other_resolved_types = comparator,
        localization_delta = target_med - comparator)
    }),idcol="target_index")
  }
  ans <- summarize_medians(type_medians)
  attr(ans,"original_localization") <- summarize_medians(original_medians)
  ans
}
original_observed_rows <- list()

message("Score observed detectable signatures")
observed_matrix <- fast_ucell(signatures)
if (!setequal(colnames(observed_matrix), names(signatures))) stop("Observed score columns mismatch")
observed_rows <- rbindlist(lapply(names(signatures), function(prog) {
  targets <- tests[program == prog, target_cell_type]
  loc <- calc_localization(observed_matrix[, prog, drop = FALSE], targets)
  original_loc <- copy(attr(loc,"original_localization"))
  original_loc[, `:=`(program=prog,test_id=paste(prog,target_cell_type,sep="__"))]
  original_observed_rows[[prog]] <<- original_loc
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
original_null_rows <- list()
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
    original_program_null <- setNames(lapply(target_types, function(x) numeric(B)), target_types)
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
      original_loc <- attr(loc,"original_localization")
      for (target in target_types) {
        program_null[[target]][first:last] <- loc[target_cell_type == target, localization_delta]
        original_program_null[[target]][first:last] <- original_loc[target_cell_type == target, localization_delta]
      }
      rm(random_matrix, loc, random_sets)
      gc(verbose = FALSE)
    }
    checkpoint <- rbindlist(lapply(target_types,function(target) data.table(
      match_mode=mode_name,program=program,target_cell_type=target,iteration=seq_len(B),
      revised_delta=program_null[[target]],original_delta=original_program_null[[target]])))
    checkpoint_dir <- file.path(out_dir,"checkpoints")
    dir.create(checkpoint_dir,showWarnings=FALSE)
    fwrite(checkpoint,file.path(checkpoint_dir,paste0(mode_name,"__",program,".tsv.gz")),sep="\t",quote=FALSE,compress="gzip")
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
      original_null_rows[[null_key]] <- data.table(
        match_mode=mode_name,program=program,target_cell_type=target,
        test_id=paste(program,target,sep="__"),iteration=seq_len(B),
        null_localization_delta=original_program_null[[target]])
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

# The same null draws are evaluated under the original labels for a direct numerical
# reproduction check. No expression, rank, gene-set, matching, statistic or BH rule changes.
original_null <- rbindlist(original_null_rows)
original_observed <- rbindlist(original_observed_rows)[,.(test_id,program,target_cell_type,
  target_median,comparator_median_of_other_resolved_types,localization_delta)]
write_tsv(original_observed,file.path(out_dir,"ORIGINAL_OBSERVED_REPRODUCED.tsv"))
fwrite(original_null,file.path(out_dir,"ORIGINAL_NULL_REPRODUCED.tsv.gz"),sep="\t",quote=FALSE,compress="gzip")
reference_dir <- normalizePath(args[[6]],mustWork=TRUE)
archived <- fread(file.path(reference_dir,"MATCHED_NULL_DISTRIBUTIONS.tsv.gz"))
keys <- c("match_mode","test_id","iteration")
if(!is_final) {
 original_for_check <- original_null[match_mode=="primary" & program=="MAPK_AP1"]
 archived <- archived[match_mode=="primary" & program=="MAPK_AP1" & iteration<=B]
} else original_for_check <- original_null
check <- merge(original_for_check[,c(keys,"null_localization_delta"),with=FALSE],
 archived[,c(keys,"null_localization_delta"),with=FALSE],by=keys,suffixes=c("_reproduced","_archived"))
max_null_error <- max(abs(check$null_localization_delta_reproduced-check$null_localization_delta_archived))
stopifnot(nrow(check)==nrow(original_for_check),max_null_error<1e-12)
archived_results <- fread(file.path(reference_dir,"MATCHED_NULL_RESULTS.tsv"))
observed_check <- merge(original_observed,unique(archived_results[,.(test_id,localization_delta)]),by="test_id")
max_observed_error <- max(abs(observed_check$localization_delta.x-observed_check$localization_delta.y))
stopifnot(max_observed_error<1e-12,identical(ann_before,obj[[ann_field,drop=TRUE]]))
comparison <- merge(results,archived_results,by=c("match_mode","test_id","program","target_cell_type"),suffixes=c("_c19_revised","_original"))
write_tsv(comparison,file.path(out_dir,"C19_VS_ORIGINAL_MATCHED_NULL.tsv"))
qa <- data.table(check=c("source_object_read_only","original_annotation_unchanged","only_cluster19_reassigned",
 "same_null_draws_reproduce_original","original_observed_delta_reproduced","six_tests_BH_per_matching_mode",
 "same_gene_sets_matching_UCell_and_statistic","all_requested_iterations_completed"),status="PASS",
 detail=c("No saveRDS of original object",as.character(identical(ann_before,obj[[ann_field,drop=TRUE]])),
 "319 UNRESOLVED cells to Fibrocyte_mesenchymal; other cells unchanged",
 paste0(nrow(check)," archived null values; max_abs_error=",format(max_null_error,scientific=TRUE)),
 format(max_observed_error,scientific=TRUE),"6 tests separately in primary and length_sensitivity",
 "seed=20260823; maxRank=1500; average ties; fixed strata and signature sizes",as.character(B)))
write_tsv(qa,file.path(out_dir,"VALIDATION.tsv"))
capture.output(sessionInfo(),file=file.path(out_dir,"sessionInfo.txt"))
message("Completed paired original versus cluster19 annotation sensitivity; B=",B)
