#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop("Usage: run_scdblfinder_per_sample.R <root>")
root <- normalizePath(args[[1]], mustWork=TRUE)
manifest <- read.delim(file.path(root, "00_manifest", "S1_SAMPLE_MANIFEST.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
qc <- read.delim(file.path(root, "02_qc", "S1_QC_CELLS.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
qc$keep <- tolower(as.character(qc$keep)) == "true"
out_dir <- file.path(root, "03_doublet", "per_sample")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
base_seed <- 20260825L

read_lines_gz <- function(path) {
  con <- gzfile(path, open="rt")
  on.exit(close(con))
  readLines(con)
}

for (i in seq_len(nrow(manifest))) {
  m <- manifest[i, ]
  message("scDblFinder ", m$GSM, " (", i, "/", nrow(manifest), ")")
  barcodes <- read_lines_gz(m$barcode_file)
  global_barcodes <- paste0(m$GSM, "_", barcodes)
  keep_ids <- qc$barcode[qc$GSM == m$GSM & qc$keep]
  keep_index <- match(keep_ids, global_barcodes)
  if (anyNA(keep_index)) stop("QC barcode mapping failed for ", m$GSM)
  con <- gzfile(m$matrix_file, open="rt")
  counts <- Matrix::readMM(con)
  close(con)
  if (!all(dim(counts) == c(m$n_genes_raw, m$n_cells_raw))) stop("Matrix dimensions changed for ", m$GSM)
  counts <- as(counts[, keep_index, drop=FALSE], "dgCMatrix")
  if (length(counts@x) && any(counts@x != floor(counts@x))) stop("Non-integer count detected for ", m$GSM)
  colnames(counts) <- keep_ids
  rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
  sce <- SingleCellExperiment(list(counts=counts))
  sample_seed <- base_seed + i - 1L
  set.seed(sample_seed)
  sce <- scDblFinder(sce, BPPARAM=SerialParam(), verbose=TRUE)
  result <- data.frame(
    barcode=colnames(sce),
    GSM=m$GSM,
    sample_id=m$sample_id,
    doublet_score=colData(sce)$scDblFinder.score,
    doublet_class=as.character(colData(sce)$scDblFinder.class),
    keep=as.character(colData(sce)$scDblFinder.class) == "singlet",
    seed=sample_seed,
    stringsAsFactors=FALSE
  )
  write.table(result, file.path(out_dir, paste0(m$GSM, "_scDblFinder.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  rm(counts, sce); invisible(gc())
}

session <- capture.output(sessionInfo())
writeLines(session, file.path(root, "14_methods", "R_SCDblFinder_SESSIONINFO.txt"))
