options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

obj_path <- "/path/to/cochlear_analysis/GSE252285_stage1_20260814/06_objects/no_doublet_robust_qc_tsne_sensitivity/02_no_doublet_annotated_tsne.rds"
obj <- readRDS(obj_path)
cat("class\t", paste(class(obj), collapse = ","), "\n", sep = "")
cat("dimensions\t", nrow(obj), "\t", ncol(obj), "\n", sep = "")
cat("assays\t", paste(Assays(obj), collapse = ","), "\n", sep = "")
cat("default_assay\t", DefaultAssay(obj), "\n", sep = "")
for (a in Assays(obj)) {
  cat("layers_", a, "\t", paste(Layers(obj[[a]]), collapse = ","), "\n", sep = "")
}
cat("reductions\t", paste(Reductions(obj), collapse = ","), "\n", sep = "")
cat("meta_columns\t", paste(colnames(obj[[]]), collapse = ","), "\n", sep = "")
ann <- obj$broad_annotation_no_doublet_sensitivity
stopifnot(length(ann) == ncol(obj), !anyNA(ann))
cat("annotation_class\t", paste(class(ann), collapse = ","), "\n", sep = "")
cat("annotation_levels\t", paste(levels(as.factor(ann)), collapse = "|"), "\n", sep = "")
print(table(ann, useNA = "ifany"))
cat("cell_head\t", paste(head(colnames(obj)), collapse = "|"), "\n", sep = "")
cat("cell_tail\t", paste(tail(colnames(obj)), collapse = "|"), "\n", sep = "")
for (red in Reductions(obj)) {
  em <- Embeddings(obj, red)
  cat("embedding_", red, "\t", nrow(em), "\t", ncol(em), "\t", paste(colnames(em), collapse = ","), "\n", sep = "")
}
cat("feature_head\t", paste(head(rownames(obj)), collapse = "|"), "\n", sep = "")
cat("feature_unique\t", !anyDuplicated(rownames(obj)), "\n", sep = "")
cat("session\n")
print(sessionInfo())
