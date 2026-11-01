source(file.path(Sys.getenv("GSE252285_SCRIPT_DIR", file.path(Sys.getenv("PROJECT_DIR", "."), "03_scripts")), "common.R"))
suppressPackageStartupMessages(library(Seurat))

stage_message("Stage 06: descriptive composition and release summary")
obj <- readRDS(p("06_objects", "04_annotated_seurat.rds"))
meta <- as.data.table(obj[[]])

composition <- meta[, .(cell_count = .N), by = .(library_id, final_broad_label)]
composition[, fraction_within_library := cell_count / sum(cell_count), by = library_id]
composition[, analysis_scope := "DESCRIPTIVE_ONLY_NO_CONDITION_P_VALUES"]
setorder(composition, library_id, final_broad_label)
write_tsv(composition, p("04_tables", "CELLTYPE_COMPOSITION_DESCRIPTIVE.tsv"))

p_comp <- ggplot(composition,
                 aes(library_id, fraction_within_library, fill = final_broad_label)) +
  geom_col() + scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Library", y = "Fraction within library", fill = "Broad cell type",
       subtitle = "Descriptive only; no condition-level inference") + theme_sc()
save_dual(p_comp, "CELLTYPE_COMPOSITION_STACKED_BAR", 10, 7)

audit <- fread(p("00_manifest", "INPUT_AUDIT.tsv"))
flow <- fread(p("04_tables", "QC_FLOW.tsv"))
res <- fread(p("04_tables", "CLUSTER_RESOLUTION_AUDIT.tsv"))
lib_counts <- fread(p("04_tables", "CELLS_PER_LIBRARY_FINAL.tsv"))
annotation <- fread(p("04_tables", "BROAD_ANNOTATION.tsv"))
top20_audit <- fread(p("04_tables", "TOP20_FEATURE_AUDIT.tsv"))
qc_warnings <- fread(p("09_qc", "QC_WARNINGS.tsv"))
confound <- fread(p("09_qc", "LIBRARY_CONFOUNDING_AUDIT.tsv"))
map <- fread(p("00_manifest", "library_map_template.tsv"))

input_dim <- audit[metric == "matrix_dimension", value]
initial_n <- flow[step == "INITIAL_CELLS", n_cells]
qc_n <- flow[step == "QC_RETAINED", n_cells]
doublet_n <- flow[step == "DOUBLETS_REMOVED", n_cells]
final_n <- flow[step == "FINAL_CELLS", n_cells]
clusters_main <- res[as.character(resolution) == "0.6", number_of_clusters]
broad_counts <- meta[, .N, by = final_broad_label][order(final_broad_label)]
unresolved_clusters <- annotation[final_broad_label == "UNRESOLVED", cluster]
found_genes <- top20_audit[status != "NOT_FOUND", gene]
not_found_genes <- top20_audit[status == "NOT_FOUND", gene]
warn_text <- qc_warnings[warning != "NONE",
                         paste0(library_id, ":", warning, collapse = "; ")]
if (!length(warn_text) || !nzchar(warn_text)) warn_text <- "NONE"
confound_status <- confound[metric == "status", value]
condition_status <- if (any(map$condition == "UNKNOWN"))
  "CONDITION_MAPPING_UNRESOLVED" else "MAPPED_BY_USER"

summary_lines <- c(
  "# GSE252285 First-stage Single-cell Module Summary",
  "",
  paste0("- INPUT_DIMENSION: ", input_dim),
  paste0("- INITIAL_CELLS: ", initial_n),
  paste0("- QC_RETAINED: ", qc_n),
  paste0("- DOUBLETS_REMOVED: ", doublet_n),
  paste0("- FINAL_CELLS: ", final_n),
  paste0("- NUMBER_OF_CLUSTERS_RES0.6: ", clusters_main),
  paste0("- CELLS_PER_LIBRARY: ", paste0(lib_counts$library_id, "=", lib_counts$cell_number,
                                           collapse = "; ")),
  paste0("- BROAD_CELLTYPE_COUNTS: ", paste0(broad_counts$final_broad_label, "=", broad_counts$N,
                                               collapse = "; ")),
  paste0("- UNRESOLVED_CLUSTERS: ", ifelse(length(unresolved_clusters),
                                            paste(unresolved_clusters, collapse = ","), "NONE")),
  paste0("- TOP20_FOUND: ", length(found_genes), "/20 [", paste(found_genes, collapse = ","), "]"),
  paste0("- TOP20_NOT_FOUND: ", length(not_found_genes), "/20 [",
         ifelse(length(not_found_genes), paste(not_found_genes, collapse = ","), "NONE"), "]"),
  paste0("- MAIN_QC_WARNINGS: ", warn_text),
  paste0("- LIBRARY_CONFOUNDING_WARNING: ", confound_status),
  "- ANNOTATION_STATUS: PROVISIONAL_BROAD_ANNOTATION",
  paste0("- CONDITION_MAPPING_STATUS: ", condition_status),
  "",
  "## Reconstruction audit against public author baseline",
  "",
  paste0("- Initial cells: observed ", initial_n, " vs public baseline ~39092 (difference ", initial_n - 39092, ")"),
  paste0("- QC-passed comparison: pipeline post-QC ", qc_n,
         " vs public baseline ~36474 (difference ", qc_n - 36474,
         "); no parameters were tuned to reproduce this number"),
  paste0("- Clusters: resolution 0.6 observed ", clusters_main,
         " vs public baseline ~31 (difference ", clusters_main - 31,
         "); resolution was pre-specified and not tuned"),
  "",
  "No cell-level Control/HHL/NIHL differential expression or differential-abundance inference was performed."
)
writeLines(summary_lines, p("10_release", "MODULE_SUMMARY.md"), useBytes = TRUE)
stage_message("Stage 06 complete")
