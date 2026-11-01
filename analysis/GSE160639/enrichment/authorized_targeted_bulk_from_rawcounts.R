#!/usr/bin/env Rscript
# AUTHORIZED 2026-11-01: targeted enrichment chained to the new raw-count DESeq2 rebuild.
# Scope: the existing nine targeted pathways x three unchanged bulk contrasts only.
# The 12 samples, 19673-gene universe, gene-set memberships and enrichment parameters are fixed.
# Upstream result tables and reference caches remain read-only.
# Input hashes are fixed to the completed 2026-11-01 public-count DESeq2 reconstruction.

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5L || !args[[1L]] %in% c("--preflight-only", "--execute-approved")) {
    stop(paste("NOT EXECUTED. Requires prior user approval and:",
      "Rscript --vanilla authorized_targeted_bulk_from_rawcounts.R --preflight-only|--execute-approved",
      "NEW_DE_RESULTS_DIRECTORY LOCKED_SAMPLE_MANIFEST EXISTING_MSIGDB_CACHE NEW_OUTPUT_DIRECTORY"))
  }
  input_root <- normalizePath(args[[2L]], mustWork = TRUE)
  manifest_path <- normalizePath(args[[3L]], mustWork = TRUE)
  cache_root <- normalizePath(args[[4L]], mustWork = TRUE)
  output_dir <- path.expand(args[[5L]])
  if (file.exists(output_dir) || dir.exists(output_dir)) {
    stop("Output already exists; refuse to overwrite any previous result")
  }
  if (as.character(getRversion()) != "4.5.3") stop("Requires existing R 4.5.3")
  if (!requireNamespace("fgsea", quietly = TRUE) ||
      as.character(utils::packageVersion("fgsea")) != "1.36.2") {
    stop("Requires existing fgsea 1.36.2; this script never installs packages")
  }
  if (!nzchar(Sys.which("sha256sum"))) stop("Existing sha256sum required")
  options(stringsAsFactors = FALSE, width = 180)
  Sys.setlocale("LC_COLLATE", "C")
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  sha256 <- function(path) {
    ans <- system2(Sys.which("sha256sum"), shQuote(path), stdout = TRUE)
    if (!is.null(attr(ans, "status"))) stop("Could not hash ", path)
    sub(" .*", "", ans[[1L]])
  }
  write_tsv <- function(x, name) {
    utils::write.table(x, file.path(output_dir, name), sep = "\t",
      quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  }
  contrasts <- data.frame(
    contrast = c("94dB_vs_70dB", "105dB_vs_70dB", "105dB_vs_94dB"),
    numerator = c("94dB", "105dB", "105dB"),
    denominator = c("70dB", "70dB", "94dB"),
    sha256 = c(
      "8e888afbdf32c31b6e87f8a054624584d247b955f369fa596c46da0b74f26dde",
      "02ad6ad60a08cdc457fac67983e9daaf1064601264c17c8f41ad8c03d9bfedb7",
      "5cbc2ec718c846130fcfa06c49dd1337034648125f78e693788ddd29cea9ae4b")
  )
  contrasts$path <- file.path(input_root,
    paste0("NEW_DESeq2_", contrasts$contrast, "_all_genes.csv.gz"))
  go_path <- file.path(cache_root, "msigdb.2026.1.Mm.M5.rds")
  hall_path <- file.path(cache_root, "msigdb.2026.1.Mm.MH.rds")
  input_hashes <- rbind(
    data.frame(input = contrasts$contrast, path = contrasts$path,
               expected_sha256 = contrasts$sha256),
    data.frame(input = c("locked_sample_manifest", "native_mouse_GO", "native_mouse_Hallmark"),
      path = c(manifest_path, go_path, hall_path),
      expected_sha256 = c(
        "9a31de34c1e36e06d70068e72c22498306f7ea662c7a94e541dd3e016ab54164",
        "ae55ef78392cbf58fb728f0fc67618a03bdfab87f76bbe78edd6c5f5c01bea13",
        "1309ce698e76c7e698bca729b69f53ccf354171b880f79f775a69d6abd8f6c82")))
  if (!all(file.exists(input_hashes$path))) stop("A frozen input/cache file is missing")
  input_hashes$observed_sha256 <- unname(vapply(input_hashes$path, sha256, character(1)))
  if (!all(input_hashes$observed_sha256 == input_hashes$expected_sha256)) {
    stop("A frozen input hash differs; no analysis or output creation")
  }
  sample_manifest <- utils::read.csv(manifest_path, check.names = FALSE)
  if (nrow(sample_manifest) != 12L || anyDuplicated(sample_manifest$gsm) ||
      !identical(as.integer(table(factor(sample_manifest$dose_db, levels = c(70, 94, 105)))),
                 c(4L, 4L, 4L))) stop("Expected 12 locked libraries, four per dose")

  # Read ONLY existing native-mouse cache files. No msigdbr download/API call.
  go <- readRDS(go_path)
  hall <- readRDS(hall_path)
  need <- c("gs_name", "db_gene_symbol", "gs_collection", "gs_subcollection",
            "db_version", "db_target_species")
  if (!all(need %in% names(go)) || !all(need %in% names(hall))) stop("Cache schema differs")
  if (!all(go$db_version == "2026.1.Mm") || !all(hall$db_version == "2026.1.Mm") ||
      !all(go$db_target_species == "MM") || !all(hall$db_target_species == "MM")) {
    stop("Requires frozen native-mouse MSigDB 2026.1.Mm")
  }
  wanted <- c(
    "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY", "GOBP_PROTEIN_FOLDING",
    "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS", "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
    "HALLMARK_APOPTOSIS", "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
  members <- rbind(
    go[go$gs_collection == "M5" & go$gs_subcollection == "GO:BP" & go$gs_name %in% wanted, need],
    hall[hall$gs_collection == "MH" & hall$gs_name %in% wanted, need])
  members <- unique(members[!is.na(members$db_gene_symbol) & nzchar(members$db_gene_symbol), ])
  if (!setequal(unique(members$gs_name), wanted)) stop("One of the exact nine pathways is absent")
  pathways <- setNames(lapply(wanted, function(w) sort(unique(
    as.character(members$db_gene_symbol[members$gs_name == w])))), wanted)
  ranks <- vector("list", nrow(contrasts))
  duplicates <- vector("list", nrow(contrasts))
  qc <- vector("list", nrow(contrasts))
  overlaps <- vector("list", nrow(contrasts))
  for (i in seq_len(nrow(contrasts))) {
    con <- gzfile(contrasts$path[[i]], open = "rt")
    raw <- tryCatch(utils::read.csv(con, check.names = FALSE), finally = close(con))
    if (!all(c("ensembl_id", "gene_symbol", "stat", "numerator", "denominator") %in% names(raw))) {
      stop("Frozen DE table schema differs")
    }
    if (nrow(raw) != 19673L || anyNA(raw$numerator) || anyNA(raw$denominator) ||
        !all(raw$numerator == contrasts$numerator[[i]]) ||
        !all(raw$denominator == contrasts$denominator[[i]])) stop("Input/direction mismatch")
    valid <- is.finite(raw$stat) & !is.na(raw$gene_symbol)
    z <- raw[valid, c("gene_symbol", "ensembl_id", "stat")]
    z$gene_symbol <- trimws(as.character(z$gene_symbol))
    z <- z[nzchar(z$gene_symbol), ]
    before_dedup <- nrow(z)
    # Fixed rule carried forward from the authorized targeted-enrichment specification.
    z <- z[order(z$gene_symbol, -abs(z$stat), -z$stat, z$ensembl_id, method = "radix"), ]
    dupe <- duplicated(z$gene_symbol) | duplicated(z$gene_symbol, fromLast = TRUE)
    duplicates[[i]] <- z[dupe, ]
    z <- z[!duplicated(z$gene_symbol), ]
    z <- z[order(-z$stat, z$gene_symbol, z$ensembl_id, method = "radix"), ]
    # No jitter, expression/DEG/P/FDR/LFC filtering, or alternative ranking is allowed.
    if (anyDuplicated(z$gene_symbol) || !all(diff(z$stat) < 0)) {
      stop("Non-unique symbols or tied statistics: stop for review, without changing the rule")
    }
    ranks[[i]] <- data.frame(rank = seq_len(nrow(z)), z, row.names = NULL)
    sizes <- vapply(pathways, function(g) length(intersect(g, z$gene_symbol)), integer(1))
    if (any(sizes < 10L | sizes > 500L)) stop("A target falls outside fixed 10–500 overlap limits")
    overlaps[[i]] <- data.frame(contrast = contrasts$contrast[[i]], pathway = wanted,
      source_size = lengths(pathways), overlap_size = sizes, row.names = NULL)
    qc[[i]] <- data.frame(contrast = contrasts$contrast[[i]], input_rows = nrow(raw),
      finite_stat_and_valid_symbol_rows = before_dedup, ranked_symbols = nrow(z),
      duplicate_rows_removed = before_dedup - nrow(z), positive = sum(z$stat > 0),
      negative = sum(z$stat < 0), zeros = sum(z$stat == 0), tied_statistics = 0)
  }
  # A dry run validates inputs/ranks/membership without writing or running enrichment.
  if (args[[1L]] == "--preflight-only") {
    print(input_hashes[c("input", "observed_sha256")], row.names = FALSE)
    print(do.call(rbind, qc), row.names = FALSE)
    print(do.call(rbind, overlaps), row.names = FALSE)
    message("PREFLIGHT_PASS: hashes, versions, 12 samples, all ranks and 27 target overlaps checked; no enrichment run")
    return(invisible(NULL))
  }
  # All frozen-input and pathway preflight checks precede output creation.
  if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) stop("Cannot create new output")
  completed <- FALSE
  on.exit(writeLines(if (completed) "COMPLETED_TARGETED_ENRICHMENT_FROM_NEW_DESEQ2" else "INCOMPLETE_DO_NOT_USE",
                     file.path(output_dir, "STATUS.txt")), add = TRUE)
  writeLines(c("Current targeted enrichment from the explicitly reconstructed public-count DESeq2 tables.",
    "Nine fixed gene sets, three contrasts and the unchanged global 27-test BH family.",
    "Upstream tables, sample metadata, labels and reference files are read-only."),
    file.path(output_dir, "README.txt"))
  write_tsv(input_hashes, "INPUT_SHA256.tsv")
  write_tsv(members, "TARGETED_PATHWAY_MEMBERSHIP.tsv")
  write_tsv(do.call(rbind, overlaps), "PATHWAY_OVERLAP_AUDIT.tsv")
  write_tsv(do.call(rbind, qc), "RANKING_QC.tsv")
  parameters <- data.frame(parameter = c("status", "R", "fgsea", "MSigDB", "ranking",
      "duplicate_rule", "seed_each_contrast", "sampleSize", "nPermSimple", "minSize", "maxSize",
      "eps", "scoreType", "gseaParam", "nproc", "RNG", "FDR_family", "FDR_criterion"),
    value = c("TARGETED_ENRICHMENT_FROM_NEW_DESEQ2", "4.5.3", "1.36.2", "2026.1.Mm",
      "signed DESeq2 Wald stat descending; all finite symbol-mapped genes",
      "maximum abs(stat), maximum signed stat, lexical Ensembl ID", "20260822", "101", "1000",
      "10", "500", "0", "std", "1", "1", "Mersenne-Twister/Inversion/Rejection",
      "BH over exactly nine pathways times three contrasts=27 tests", "global_FDR < 0.05"))
  write_tsv(parameters, "PARAMETERS.tsv")
  results <- vector("list", nrow(contrasts))
  objects <- vector("list", nrow(contrasts))
  warnings_seen <- character()
  start <- proc.time()[["elapsed"]]
  for (i in seq_len(nrow(contrasts))) {
    label <- contrasts$contrast[[i]]
    write_tsv(ranks[[i]], paste0("RANKING_", label, ".tsv"))
    write_tsv(duplicates[[i]], paste0("DUPLICATE_SYMBOL_AUDIT_", label, ".tsv"))
    message("Running authorized contrast: ", label)
    gene_stats <- setNames(ranks[[i]]$stat, ranks[[i]]$gene_symbol)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(20260822L)
    fit <- withCallingHandlers(fgsea::fgseaMultilevel(pathways = pathways, stats = gene_stats,
      sampleSize = 101L, nPermSimple = 1000L, minSize = 10L, maxSize = 500L,
      eps = 0, scoreType = "std", gseaParam = 1, nproc = 1L), warning = function(w) {
        warnings_seen <<- c(warnings_seen, paste(label, conditionMessage(w)))
      })
    fit <- as.data.frame(fit)
    if (nrow(fit) != 9L || !setequal(fit$pathway, wanted) ||
        any(!is.finite(fit$pval)) || any(!is.finite(fit$NES)) || any(!is.finite(fit$ES)) ||
        any(!is.finite(fit$log2err))) {
      stop("Incomplete/nonfinite enrichment output; no alternative method or threshold attempted")
    }
    fit <- fit[match(wanted, fit$pathway), ]
    if (any(fit$pval < 0 | fit$pval > 1) ||
        !all(fit$size == overlaps[[i]]$overlap_size) ||
        !all(vapply(fit$leadingEdge, function(g) all(g %in% names(gene_stats)), logical(1))) ||
        !all(mapply(function(g, pid) all(g %in% pathways[[pid]]), fit$leadingEdge, fit$pathway))) {
      stop("Output invariant failed")
    }
    objects[[i]] <- fit
    results[[i]] <- data.frame(contrast = label, numerator = contrasts$numerator[[i]],
      denominator = contrasts$denominator[[i]], pathway = fit$pathway, NES = fit$NES,
      ES = fit$ES, P = fit$pval, within_contrast_FDR = p.adjust(fit$pval, method = "BH"),
      size = fit$size, log2err = fit$log2err,
      leadingEdge = vapply(fit$leadingEdge, paste, collapse = ",", character(1)),
      analysis_status = "TARGETED_ENRICHMENT_FROM_NEW_DESEQ2")
  }
  combined <- do.call(rbind, results)
  stopifnot(nrow(combined) == 27L, !anyDuplicated(combined[c("contrast", "pathway")]))
  combined$global_27_test_FDR <- p.adjust(combined$P, method = "BH")
  combined$significant_global_FDR05 <- combined$global_27_test_FDR < 0.05
  write_tsv(combined, "TARGETED_27_RESULTS.tsv")
  names(objects) <- contrasts$contrast
  saveRDS(objects, file.path(output_dir, "FGSEA_OBJECTS.rds"))
  writeLines(warnings_seen, file.path(output_dir, "WARNINGS.txt"))
  capture.output(sessionInfo(), file = file.path(output_dir, "SESSIONINFO.txt"))
  write_tsv(data.frame(elapsed_seconds = proc.time()[["elapsed"]] - start,
                      statistical_tests = 27L, requested_processes = 1L), "RUNTIME.tsv")
  # Recheck input hashes to verify all upstream data/cache files remain unchanged.
  after <- unname(vapply(input_hashes$path, sha256, character(1)))
  if (!identical(after, input_hashes$observed_sha256)) stop("Input changed during run")
  completed <- TRUE
  message("Completed targeted enrichment from NEW DESeq2 tables in NEW directory: ", output_dir)
  message("Original results and manuscript were not altered; new statistics require review.")
}
main()
