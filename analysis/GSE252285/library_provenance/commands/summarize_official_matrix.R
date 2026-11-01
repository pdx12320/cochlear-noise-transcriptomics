#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Matrix))

project <- Sys.getenv(
  "PROJECT",
  "/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815"
)
matrix_path <- file.path(project, "official", "GSE252285_matrix.mtx.gz")
barcode_path <- file.path(project, "official", "GSE252285_barcodes.tsv.gz")
feature_path <- file.path(project, "official", "GSE252285_features.tsv.gz")
output_path <- file.path(project, "comparison", "official_matrix_summary.tsv")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
barcodes <- readLines(gzfile(barcode_path))
features <- read.delim(
  gzfile(feature_path), header = FALSE, sep = "\t", quote = "",
  stringsAsFactors = FALSE
)
matrix <- readMM(gzfile(matrix_path))

stopifnot(nrow(matrix) == 32285L)
stopifnot(ncol(matrix) == 39092L)
stopifnot(length(barcodes) == ncol(matrix))
stopifnot(nrow(features) == nrow(matrix))
stopifnot(length(matrix@x) == 89879619L)

genes_per_cell <- tabulate(matrix@j + 1L, nbins = ncol(matrix))
umis_per_cell <- as.numeric(colSums(matrix))
suffix <- sub(".*(-[123])$", "\\1", barcodes)
stopifnot(identical(as.integer(table(factor(suffix, levels = c("-1", "-2", "-3")))), c(17720L, 10028L, 11344L)))

summarize_group <- function(label, indices) {
  data.frame(
    suffix = label,
    cell_count = length(indices),
    median_genes_per_cell = median(genes_per_cell[indices]),
    median_umi_per_cell = median(umis_per_cell[indices]),
    total_umis = sum(umis_per_cell[indices]),
    stringsAsFactors = FALSE
  )
}

summary <- rbind(
  summarize_group("ALL", seq_along(barcodes)),
  summarize_group("-1", which(suffix == "-1")),
  summarize_group("-2", which(suffix == "-2")),
  summarize_group("-3", which(suffix == "-3"))
)
write.table(summary, output_path, sep = "\t", row.names = FALSE, quote = FALSE)

versions_path <- file.path(project, "comparison", "official_matrix_software_versions.tsv")
write.table(
  data.frame(
    software = c("R", "Matrix"),
    version = c(R.version.string, as.character(packageVersion("Matrix")))
  ),
  versions_path, sep = "\t", row.names = FALSE, quote = FALSE
)

cat("OFFICIAL_MATRIX_SUMMARY_COMPLETE", output_path, "\n")
