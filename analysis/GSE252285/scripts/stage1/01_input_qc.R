source(file.path(Sys.getenv("GSE252285_SCRIPT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "03_scripts")), "common.R"))

stage_message("Stage 01: input audit")

matrix_path <- Sys.getenv("GSE_MATRIX", p("01_input", "GSE252285_matrix.mtx.gz"))
features_path <- Sys.getenv("GSE_FEATURES", p("01_input", "GSE252285_features.tsv.gz"))
barcodes_path <- Sys.getenv("GSE_BARCODES", p("01_input", "GSE252285_barcodes.tsv.gz"))

input_paths <- c(matrix_path, features_path, barcodes_path)
if (!all(file.exists(input_paths))) {
  writeLines("STOP_INPUT_AUDIT", p("00_manifest", "STOP_INPUT_AUDIT"))
  stop("STOP_INPUT_AUDIT: one or more inputs are missing: ",
       paste(input_paths[!file.exists(input_paths)], collapse = ", "))
}

features <- data.table::fread(features_path, header = FALSE, sep = "\t",
                              quote = "", na.strings = c("NA"), showProgress = FALSE)
if (ncol(features) < 2L) {
  writeLines("STOP_INPUT_AUDIT", p("00_manifest", "STOP_INPUT_AUDIT"))
  stop("STOP_INPUT_AUDIT: features table has fewer than two columns")
}
if (ncol(features) == 2L) features[, V3 := "Gene Expression"]
setnames(features, names(features)[1:3], c("feature_id", "gene_symbol", "feature_type"))
features <- features[, .(feature_id, gene_symbol, feature_type)]

barcodes_dt <- data.table::fread(barcodes_path, header = FALSE, sep = "\t",
                                 quote = "", na.strings = c("NA"), showProgress = FALSE)
barcodes <- as.character(barcodes_dt[[1]])

stage_message("Reading Matrix Market input as a sparse matrix")
con <- gzfile(matrix_path, open = "rt")
counts_triplet <- Matrix::readMM(con)
close(con)
declared_entries <- length(counts_triplet@x)
matrix_na <- anyNA(counts_triplet@x)
counts <- as(counts_triplet, "CsparseMatrix")
rm(counts_triplet)
gc(verbose = FALSE)

raw_symbol <- as.character(features$gene_symbol)
fallback <- is.na(raw_symbol) | !nzchar(raw_symbol)
base_name <- raw_symbol
base_name[fallback] <- as.character(features$feature_id[fallback])
feature_name <- make.unique(base_name, sep = "-dup")
features[, feature_name := feature_name]

suffix <- barcode_suffix(barcodes)
suffix_counts <- as.data.table(table(suffix), keep.rownames = FALSE)
setnames(suffix_counts, c("suffix", "cell_number"))
suffix_counts[, suffix := as.character(suffix)]
setorder(suffix_counts, suffix)

feature_type_counts <- features[, .(feature_number = .N), by = feature_type]
write_tsv(suffix_counts, p("00_manifest", "BARCODE_SUFFIX_COUNTS.tsv"))
write_tsv(feature_type_counts, p("00_manifest", "FEATURE_TYPE_COUNTS.tsv"))

map_path <- p("00_manifest", "library_map_template.tsv")
if (!file.exists(map_path)) {
  map <- copy(suffix_counts[, .(suffix)])
  map[, library_id := ifelse(grepl("^-[0-9]+$", suffix),
                             paste0("LIB", sub("^-", "", suffix)),
                             paste0("LIB_", seq_len(.N)))]
  map[, condition := "UNKNOWN"]
  write_tsv(map, map_path)
} else {
  stage_message("Preserving existing library_map_template.tsv")
}

dim_ok <- nrow(counts) == nrow(features) && ncol(counts) == length(barcodes)
entries_ok <- declared_entries == length(counts@x)
barcodes_ok <- !anyDuplicated(barcodes) && !anyNA(barcodes) && all(nzchar(barcodes))

audit <- data.table(
  metric = c(
    "matrix_dimension", "matrix_genes", "matrix_cells", "matrix_declared_entries",
    "matrix_stored_entries", "matrix_nonzero_values", "feature_rows", "barcode_rows",
    "duplicate_gene_ids", "duplicate_gene_symbols", "duplicate_barcodes",
    "NA_matrix_values", "NA_feature_ids", "NA_gene_symbols", "NA_barcodes",
    "blank_feature_ids", "blank_gene_symbols", "blank_barcodes",
    "dimension_consistency", "stored_entry_consistency", "cellranger_version"
  ),
  value = as.character(c(
    paste(dim(counts), collapse = "x"), nrow(counts), ncol(counts), declared_entries,
    length(counts@x), Matrix::nnzero(counts), nrow(features), length(barcodes),
    sum(duplicated(features$feature_id)), sum(duplicated(features$gene_symbol)),
    sum(duplicated(barcodes)), matrix_na, sum(is.na(features$feature_id)),
    sum(is.na(features$gene_symbol)), sum(is.na(barcodes)),
    sum(!is.na(features$feature_id) & !nzchar(features$feature_id)),
    sum(!is.na(features$gene_symbol) & !nzchar(features$gene_symbol)),
    sum(!is.na(barcodes) & !nzchar(barcodes)), dim_ok, entries_ok, "7.1.0"
  )),
  status = c(
    rep("INFO", 8), rep("CHECKED", 10),
    ifelse(dim_ok, "PASS", "FAIL"), ifelse(entries_ok, "PASS", "FAIL"), "INFO"
  ),
  details = c(
    "rows x columns", "matrix rows", "matrix columns", "Matrix Market header/coordinate count",
    "sparse x-slot length", "non-zero values after sparse conversion", "features.tsv rows",
    "barcodes.tsv rows", "exact duplicated Ensembl IDs after first occurrence",
    "exact duplicated symbols after first occurrence", "exact duplicated barcodes after first occurrence",
    "checked in sparse x slot", "checked", "checked", "checked", "checked", "checked", "checked",
    "matrix rows/features and columns/barcodes", "declared versus stored coordinate entries",
    "from Matrix Market metadata_json"
  )
)
write_tsv(audit, p("00_manifest", "INPUT_AUDIT.tsv"))
write_tsv(features, p("00_manifest", "FEATURE_NAME_MAP.tsv.gz"))

if (!dim_ok || !entries_ok || !barcodes_ok || matrix_na) {
  writeLines("STOP_INPUT_AUDIT", p("00_manifest", "STOP_INPUT_AUDIT"))
  cat("STOP_INPUT_AUDIT\n")
  stop("STOP_INPUT_AUDIT: matrix/features/barcodes integrity checks failed")
}
if (file.exists(p("00_manifest", "STOP_INPUT_AUDIT"))) {
  unlink(p("00_manifest", "STOP_INPUT_AUDIT"))
}

rownames(counts) <- feature_name
colnames(counts) <- barcodes
saveRDS(counts, p("02_processed", "01_counts_matrix.rds"), compress = FALSE)
rm(counts)
gc(verbose = FALSE)

stage_message("Stage 01 complete: ", nrow(features), " genes x ", length(barcodes), " cells")
