suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})
project <- Sys.getenv("GSE252285_PROJECT_ROOT", unset = ".")
ref_path <- file.path(project, "02_processed", "reference_final_audit", "GSE217727",
                      "GSE217727_P7_reconstructed_reference.rds")
target_path <- file.path(project, "06_objects", "no_doublet_robust_qc_tsne_sensitivity",
                         "02_no_doublet_annotated_tsne.rds")
ref <- readRDS(ref_path)
cat("REF_DIM\t", nrow(ref), "\t", ncol(ref), "\n", sep = "")
cat("REF_LAYERS\t", paste(Layers(ref[["RNA"]]), collapse = ","), "\n", sep = "")
cat("REF_REDUCTIONS\t", paste(Reductions(ref), collapse = ","), "\n", sep = "")
cat("REF_VAR_FEATURES\t", length(VariableFeatures(ref)), "\n", sep = "")
cat("REF_META\t", paste(colnames(ref[[]]), collapse = ","), "\n", sep = "")
for (field in intersect(c("reference_subtype", "reference_broad",
                          "reference_subtype_reconstructed", "reference_broad_reconstructed",
                          "reference_signature_confident"), colnames(ref[[]]))) {
  tab <- sort(table(ref[[field, drop = TRUE]], useNA = "ifany"), decreasing = TRUE)
  cat("FIELD\t", field, "\t", paste(names(tab), as.integer(tab), sep = "=", collapse = ";"), "\n", sep = "")
}
target <- readRDS(target_path)
cluster <- as.character(target[["RNA_snn_res.0.6", drop = TRUE]])
cat("TARGET_DIM\t", nrow(target), "\t", ncol(target), "\n", sep = "")
cat("TARGET_REDUCTIONS\t", paste(Reductions(target), collapse = ","), "\n", sep = "")
cat("TARGET_CLUSTER_COUNTS\n")
print(data.table(cluster = cluster)[cluster %in% c("6", "9", "13", "19"), .N, by = cluster][order(as.integer(cluster))])
cat("TARGET_LAYERS\t", paste(Layers(target[["RNA"]]), collapse = ","), "\n", sep = "")
cat("SEURAT_VERSION\t", as.character(packageVersion("Seurat")), "\n", sep = "")
