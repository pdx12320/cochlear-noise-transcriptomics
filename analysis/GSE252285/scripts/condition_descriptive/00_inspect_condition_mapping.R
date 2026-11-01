suppressPackageStartupMessages({
  library(SeuratObject)
  library(data.table)
})

project_dir <- Sys.getenv("GSE252285_PROJECT_ROOT", unset = ".")
object_path <- file.path(
  project_dir, "06_objects", "no_doublet_robust_qc_tsne_sensitivity",
  "02_no_doublet_annotated_tsne.rds"
)
obj <- readRDS(object_path)
meta <- as.data.table(obj[[]], keep.rownames = "cell_name")
meta[, parsed_suffix := sub("^.*(-[123])$", "\\1", cell_name)]
meta[!grepl("-[123]$", cell_name), parsed_suffix := NA_character_]

cat("OBJECT_DIMENSION\t", nrow(obj), "\t", ncol(obj), "\n", sep = "")
cat("REDUCTIONS\t", paste(Reductions(obj), collapse = ","), "\n", sep = "")
cat("METADATA_COLUMNS\t", paste(colnames(meta), collapse = ","), "\n", sep = "")
cat("CELL_NAME_EXAMPLES\n")
print(meta[, .(cell_name, parsed_suffix)][1:12])

candidate_fields <- intersect(
  c("orig.ident", "library_id", "barcode_suffix", "condition",
    "RNA_snn_res.0.6", "broad_annotation_no_doublet_sensitivity"),
  colnames(meta)
)
for (field in candidate_fields) {
  cat("UNIQUE_FIELD\t", field, "\t", sep = "")
  cat(paste(sort(unique(as.character(meta[[field]]))), collapse = ","), "\n", sep = "")
}

cross_fields <- intersect(c("parsed_suffix", "barcode_suffix", "library_id", "orig.ident", "condition"), colnames(meta))
cat("MAPPING_CROSSTAB\n")
print(meta[, .N, by = cross_fields][order(parsed_suffix, library_id)])
cat("CELL_COUNT\t", nrow(meta), "\n", sep = "")
