suppressPackageStartupMessages(library(data.table))

project_dir <- normalizePath(
  Sys.getenv("PROJECT_DIR", Sys.getenv("GSE252285_PROJECT_ROOT", unset = ".")),
  mustWork = TRUE
)
branch <- "no_doublet_robust_qc_tsne_sensitivity"
audit <- "cluster_6_9_13_19_annotation_audit"
audit_path <- function(section, ...) file.path(project_dir, section, branch, audit, ...)
for (section in c("00_manifest", "04_tables", "05_figures", "07_logs", "08_methods", "09_qc", "10_release")) {
  dir.create(audit_path(section), recursive = TRUE, showWarnings = FALSE)
}
write_tsv <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")

target_path <- file.path(project_dir, "06_objects", branch, "02_no_doublet_annotated_tsne.rds")
table_path <- function(stem) audit_path("04_tables", stem)
required <- c(
  target_path,
  table_path("FOUR_CLUSTER_FINDMARKERS_ALL.tsv.gz"),
  table_path("FOUR_CLUSTER_TOP50_ANNOTATION_MARKERS.tsv"),
  table_path("FOUR_CLUSTER_PROGRAM_GENE_STATS.tsv"),
  table_path("FOUR_CLUSTER_LINEAGE_PROGRAM_SCORES.tsv"),
  table_path("FOUR_CLUSTER_AUTHOR_MARKER_STATS.tsv"),
  table_path("FOUR_CLUSTER_REFERENCE_MAPPING_CELL_LEVEL.tsv.gz"),
  table_path("FOUR_CLUSTER_REFERENCE_MAPPING_SUMMARY.tsv"),
  table_path("FOUR_CLUSTER_TECHNICAL_AUDIT.tsv")
)
if (any(!file.exists(required))) stop("Missing finalization input: ", paste(required[!file.exists(required)], collapse = ", "))

markers <- fread(table_path("FOUR_CLUSTER_FINDMARKERS_ALL.tsv.gz"), showProgress = FALSE)
top50 <- fread(table_path("FOUR_CLUSTER_TOP50_ANNOTATION_MARKERS.tsv"))
gene_stats <- fread(table_path("FOUR_CLUSTER_PROGRAM_GENE_STATS.tsv"))
program_scores <- fread(table_path("FOUR_CLUSTER_LINEAGE_PROGRAM_SCORES.tsv"))
author_stats <- fread(table_path("FOUR_CLUSTER_AUTHOR_MARKER_STATS.tsv"))
reference_cell <- fread(table_path("FOUR_CLUSTER_REFERENCE_MAPPING_CELL_LEVEL.tsv.gz"), showProgress = FALSE)
reference_summary <- fread(table_path("FOUR_CLUSTER_REFERENCE_MAPPING_SUMMARY.tsv"))
technical <- fread(table_path("FOUR_CLUSTER_TECHNICAL_AUDIT.tsv"))

target_clusters <- c("6", "9", "13", "19")
for (x in list(markers, top50, gene_stats, program_scores, author_stats, reference_cell,
               reference_summary, technical)) x[, cluster := as.character(cluster)]

if (!identical(sort(unique(markers$cluster)), sort(target_clusters))) stop("Marker cluster invariant failed")
if (!all(top50[, .N, by = cluster]$N == 50L)) stop("Top50 row-count invariant failed")
if (nrow(reference_cell) != 4726L || sum(technical$cell_number) != 4726L) stop("Target cell conservation failed")
if (!identical(technical[order(as.integer(cluster)), cell_number], c(2359L, 1210L, 838L, 319L))) {
  stop("Target cluster cell counts changed")
}

# Consolidate within-cluster detection percentages from all independently generated evidence tables.
pct_sources <- rbindlist(list(
  markers[, .(cluster, gene, pct_in = pct.1, source_priority = 3L)],
  unique(gene_stats[, .(cluster, gene, pct_in, source_priority = 1L)]),
  unique(author_stats[, .(cluster, gene, pct_in, source_priority = 2L)])
), use.names = TRUE)
setorder(pct_sources, cluster, gene, source_priority)
pct_lookup <- pct_sources[, .SD[1L], by = .(cluster, gene)]
pct_text <- function(cluster_id, genes) {
  z <- merge(data.table(gene = genes), pct_lookup[cluster == cluster_id, .(gene, pct_in)],
             by = "gene", all.x = TRUE, sort = FALSE)
  z[, paste0(gene, "=", fifelse(is.na(pct_in), "NA", sprintf("%.1f%%", 100 * pct_in)))] |>
    paste(collapse = ";")
}

decision <- data.table(
  cluster = target_clusters,
  candidate_label = c("UNRESOLVED", "UNRESOLVED", "Immune_hematopoietic", "Fibrocyte_mesenchymal"),
  candidate_lineage_hypothesis = c(
    "Erythroid-like",
    "Erythroid-like",
    "Hematopoietic progenitor-like; not myeloid-resolved",
    "Type I OMC-like/perivascular fibrocyte; subtype not forced"
  ),
  confidence = c("UNRESOLVED", "UNRESOLVED", "MODERATE", "HIGH"),
  reference_consistent = c("NO_REFERENCE_CLASS", "NO_DUE_TO_CLOSED_WORLD_CONFLICT", "YES_BROAD", "YES_BROAD_AND_SUBTYPE"),
  supporting_markers = c(
    pct_text("6", c("Klf1", "Gypa", "Gata1", "Slc4a1", "Alas2", "Car1", "Ces2g")),
    pct_text("9", c("Alas2", "Gypa", "Slc4a1", "Hbb-bt", "Hba-a1", "Hba-a2")),
    pct_text("13", c("Cd34", "Kit", "Flt3", "Myb", "Hoxa9", "Meis1", "Msi2", "Dntt", "Bcl11a", "Mef2c", "Ptprc", "Ptprcap")),
    pct_text("19", c("Igfbp5", "Coch", "Otor", "Pou3f4", "Col1a2", "Pdgfra", "Tbx18", "Rgs2"))
  ),
  conflicting_markers = c(
    pct_text("6", c("Tubb3", "Snap25", "Syt1", "Rbfox3", "Nefh", "Prph", "Ptprc", "Ctss")),
    pct_text("9", c("Tubb3", "Snap25", "Syt1", "Rbfox3", "Nefh", "Prph", "Ptprc", "Ctss")),
    pct_text("13", c("Egfl7", "Pecam1", "Adgrl4", "Pdgfrb", "Lyz2", "C1qa", "Cd68")),
    pct_text("19", c("Emcn", "Rbp1", "Abcc9", "Pdgfrb", "Kcnj8", "Cald1", "Sox9", "Epha4", "Bhlhe22", "Syt1"))
  ),
  reason = c(
    "Nine-of-nine erythroid program with large pct.1-pct.2 and effect size; Car1/Ces2g are positive but the canonical neuronal program has zero supported genes. GSE217727 lacks an erythroid class and rejects 91.6% of cells, so the strict program-plus-reference rule requires UNRESOLVED.",
    "Seven-of-nine erythroid program with near-universal Alas2/Slc4a1/globin detection and no canonical neuronal support. The P7 reference has no erythroid class and therefore maps the closed-world query to Immune despite immune/myeloid genes being less detected than in the rest; retain UNRESOLVED.",
    "Fifteen-of-seventeen hematopoietic-progenitor program plus 100% accepted P7 broad Immune prediction. Strong Egfl7/Pecam1 and Pdgfrb create endothelial/hemogenic ambiguity, and the mature myeloid program is not enriched; name only the broad hematopoietic lineage with MODERATE confidence.",
    "Eight-gene fibrocyte/mesenchymal program plus 99.7% TypeI_OMCs reference prediction (median score 0.966, margin 0.932). Endothelial/pericyte/chondrocyte/glial-associated genes indicate a specialized perivascular mesenchymal profile, not sufficient to force a finer subtype. Partial neuronal-associated genes lack the pan-neuronal core."
  )
)

ref_cols <- c("cluster", "subtype_top_raw", "subtype_top_raw_fraction", "subtype_top_accepted",
              "subtype_top_accepted_fraction_total", "subtype_rejected_fraction",
              "subtype_median_score_top", "subtype_median_margin_top",
              "broad_top_accepted", "broad_top_accepted_fraction_total")
tech_cols <- c("cluster", "cell_number", "LIB1_fraction", "LIB2_fraction", "LIB3_fraction",
               "library_total_variation_vs_rest", "median_nCount_RNA", "median_nFeature_RNA",
               "median_percent_mt", "robust_z_median_nCount_vs_rest",
               "robust_z_median_nFeature_vs_rest", "robust_z_median_percent_mt_vs_rest",
               "technical_warning")
decision <- merge(decision, reference_summary[, ..ref_cols], by = "cluster", all.x = TRUE)
decision <- merge(decision, technical[, ..tech_cols], by = "cluster", all.x = TRUE)
decision <- decision[order(as.integer(cluster))]

decision[, SGN_recovered := "NO"]
decision[, author_marker_conflict := fifelse(cluster == "6", "YES_Car1_Ces2g_WITHOUT_CANONICAL_NEURONAL_PROGRAM", "NO")]
write_tsv(decision, table_path("FOUR_CLUSTER_FINAL_ANNOTATION_AUDIT.tsv"))

flags <- data.table(
  flag = c("SGN_RECOVERED", "AUTHOR_LABEL_CONFLICT", "UNRESOLVED_REMAINING",
           "TARGET_CLUSTER_COUNT", "TARGET_CELL_COUNT", "UPSTREAM_OBJECT_MODIFIED",
           "TARGET_PCA_RERUN", "TARGET_CLUSTERING_RERUN", "CONDITION_MAPPING_STATUS"),
  value = c("NO", "YES", "6,9", "4", "4726", "NO", "NO", "NO", "CONDITION_MAPPING_UNRESOLVED")
)
write_tsv(flags, table_path("FOUR_CLUSTER_AUDIT_FLAGS.tsv"))

summary_lines <- c(
  "# GSE252285 clusters 6/9/13/19 annotation audit",
  "",
  "STATUS=SUCCESS",
  "",
  "This focused audit is read-only with respect to the 36,111-cell Seurat object. Resolution 0.6, all 28 clusters, PCA and t-SNE were held fixed. No QC, PCA, UMAP, t-SNE or clustering was rerun, and no annotation metadata was written back.",
  "",
  "## Adjudication",
  "",
  "| cluster | cells | candidate label | confidence | key evidence | reference evidence |",
  "|---:|---:|---|---|---|---|",
  "| 6 | 2359 | UNRESOLVED (erythroid-like hypothesis) | UNRESOLVED | 9/9 erythroid genes; Car1/Ces2g high but canonical neuronal support 0/17 | 91.6% rejected; P7 reference lacks erythroid |",
  "| 9 | 1210 | UNRESOLVED (erythroid-like hypothesis) | UNRESOLVED | 7/9 erythroid genes; no canonical neuronal support | closed-world Immune call conflicts with erythroid program; reference lacks erythroid |",
  "| 13 | 838 | Immune_hematopoietic (progenitor-like; not myeloid-resolved) | MODERATE | 15/17 hematopoietic-progenitor genes; endothelial/hemogenic conflict | 100% accepted broad Immune; median score 0.843, margin 0.687 |",
  "| 19 | 319 | Fibrocyte_mesenchymal (Type I OMC-like/perivascular) | HIGH broad | 8/13 fibrocyte genes; mixed vascular/pericyte/chondrocyte programs | 99.7% TypeI_OMCs; median score 0.966, margin 0.932 |",
  "",
  "SGN_RECOVERED=NO",
  "AUTHOR_LABEL_CONFLICT=YES",
  "UNRESOLVED_REMAINING=6,9",
  "",
  "## SGN decision",
  "",
  "Cluster 6 expresses the paper's SGN-associated Car1 (72.3%) and Ces2g (74.9%), but its dominant program is erythroid (Klf1/Gypa/Gata1/Slc4a1/Alas2) and the canonical neuronal program has 0/17 supported genes. Cluster 19 has four partial neuronal-associated genes, but lacks the pan-neuronal Tubb3/Snap25/Rbfox3/Nefh/Prph core and maps to Type I OMCs. SGN is therefore not recovered.",
  "",
  "## Technical warnings",
  "",
  "- No audited cluster crossed the prespecified library/QC-driver flag (library fraction >=0.90, total-variation >=0.50, or absolute robust median shift >=3 MAD).",
  "- Library composition is nevertheless imbalanced: LIB1 is only 6.2% of cluster 6 and 2.1% of cluster 9; cluster 19 is 73.4% LIB3. CONDITION_MAPPING_UNRESOLVED prevents biological interpretation of these shifts.",
  "- GSE217727 is a P7 whole-cochlea closed-world reference and lacks an erythroid class. Its Immune prediction for cluster 9 is explicitly treated as a reference-space failure, not lineage proof.",
  "- Cluster 13 has a strong progenitor program and broad reference agreement, but its endothelial/hemogenic signature prevents a mature myeloid label.",
  "- FindMarkers p-values are descriptive because clustering and testing use the same data. Labels use effect size, pct.1, pct.2, pct.1-pct.2, multi-gene coherence and reference evidence.",
  "- Exact author cluster/cell metadata were not found in the project; author marker comparisons are not claimed as cell-level author label transfer.",
  "",
  "## Reference",
  "",
  "GSE217727, P7 C57BL/6J whole mouse cochlea, reconstructed 15-subtype reference; Seurat anchor transfer with score <0.50 or top-two margin <0.10 rejected. Reference PCA was rebuilt; target PCA was not rerun.",
  "",
  "## Key files",
  "",
  "- `04_tables/.../FOUR_CLUSTER_FINAL_ANNOTATION_AUDIT.tsv`",
  "- `04_tables/.../FOUR_CLUSTER_TOP50_ANNOTATION_MARKERS.tsv`",
  "- `04_tables/.../FOUR_CLUSTER_LINEAGE_PROGRAM_SCORES.tsv`",
  "- `04_tables/.../FOUR_CLUSTER_REFERENCE_MAPPING_SUMMARY.tsv`",
  "- `04_tables/.../FOUR_CLUSTER_TECHNICAL_AUDIT.tsv`",
  "- `05_figures/.../FOUR_CLUSTER_LINEAGE_MARKER_DOTPLOT.pdf/png`",
  "- `05_figures/.../FOUR_CLUSTER_SGN_AUTHOR_VS_CANONICAL_DOTPLOT.pdf/png`",
  "- `05_figures/.../FOUR_CLUSTER_REFERENCE_MAPPING_HEATMAP.pdf/png`"
)
writeLines(summary_lines, audit_path("10_release", "FOUR_CLUSTER_ANNOTATION_AUDIT_SUMMARY.md"))

methods_lines <- c(
  "# Focused annotation-audit methods",
  "",
  "The existing no-doublet robust-QC t-SNE sensitivity object was read only. The fixed RNA_snn_res.0.6 identities were used for clusters 6, 9, 13 and 19. No upstream reduction or graph calculation was run.",
  "",
  "For each cluster, Seurat FindMarkers compared the cluster with all remaining cells using the Wilcoxon test, only.pos=TRUE, min.pct=0.05 and logfc.threshold=0. Ranking prioritized avg_log2FC multiplied by max(pct.1-pct.2,0). P-values were retained but were not used alone for annotation. The fixed Top20 candidate genes were excluded from annotation-marker ranking.",
  "",
  "Canonical multi-gene programs covered SGN/neuronal, Schwann, satellite glia, chondrocyte, fibrocyte/mesenchymal, supporting epithelium, hair cells, endothelium, pericyte/mural, SMC, immune/myeloid, erythroid, hematopoietic progenitor and strial intermediate/melanocyte lineages. The hematopoietic-progenitor program was added after the unbiased Top markers nominated that lineage; it was evaluated as a multi-gene module rather than by a single marker.",
  "",
  "Independent reference mapping used the existing reconstructed GSE217727 P7 whole-cochlea reference. Reference PCA was rebuilt from 2,914 shared variable features and 30 PCs. Target cells retained their existing PCA. Seurat FindTransferAnchors/TransferData used pcaproject; predictions with score below 0.50 or top-two margin below 0.10 were rejected. Closed-world limitations, especially the absence of erythroid cells, were explicitly considered.",
  "",
  "A candidate label was accepted only when a coherent multi-gene program and reference evidence agreed. Conflicting or missing reference evidence resulted in UNRESOLVED. Library fractions and robust median shifts in nCount_RNA, nFeature_RNA and percent.mt were audited without altering thresholds or identities."
)
writeLines(methods_lines, audit_path("08_methods", "METHODS_FOUR_CLUSTER_ANNOTATION_AUDIT.md"))

# Validate immutable inputs and figure/file structure.
sha_before <- audit_path("00_manifest", "INPUT_SHA256_BEFORE.tsv")
sha_after_candidates <- c(
  audit_path("00_manifest", "INPUT_SHA256_AFTER_FINAL_QA_RERUN.tsv"),
  audit_path("00_manifest", "INPUT_SHA256_AFTER_FIGURE_QA_RERUN.tsv"),
  audit_path("00_manifest", "INPUT_SHA256_AFTER.tsv")
)
sha_after <- sha_after_candidates[file.exists(sha_after_candidates)][1]
if (is.na(sha_after) || !identical(readLines(sha_before), readLines(sha_after))) stop("Immutable input checksum mismatch")

figure_stems <- c("FOUR_CLUSTER_LINEAGE_MARKER_DOTPLOT", "FOUR_CLUSTER_LINEAGE_SCORE_HEATMAP",
                  "FOUR_CLUSTER_SGN_AUTHOR_VS_CANONICAL_DOTPLOT", "FOUR_CLUSTER_REFERENCE_MAPPING_HEATMAP",
                  "FOUR_CLUSTER_LIBRARY_QC_AUDIT", "FOUR_CLUSTER_TOP_MARKER_HEATMAP")
figure_files <- unlist(lapply(figure_stems, function(stem) audit_path("05_figures", paste0(stem, c(".pdf", ".png")))))
if (any(!file.exists(figure_files)) || any(file.info(figure_files)$size <= 1000)) stop("Figure output missing or empty")
pdf_files <- grep("\\.pdf$", figure_files, value = TRUE)
if (any(vapply(pdf_files, function(p) rawToChar(readBin(p, "raw", 4L)) != "%PDF", logical(1)))) stop("Invalid PDF signature")
png_files <- grep("\\.png$", figure_files, value = TRUE)
png_sig <- as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
if (any(vapply(png_files, function(p) !identical(readBin(p, "raw", 8L), png_sig), logical(1)))) stop("Invalid PNG signature")
if (length(list.files(audit_path("06_objects"), recursive = TRUE, pattern = "\\.rds$", full.names = TRUE))) stop("Unexpected audit RDS output")

validation <- data.table(
  check = c("target_clusters_exact", "target_cells_conserved", "top50_exact", "cell_reference_rows",
            "decision_rows", "SGN_recovered_no", "source_sha256_unchanged", "pdf_png_pairs",
            "figures_nonempty_valid_signature", "no_new_RDS", "condition_DE_not_run"),
  status = "PASS",
  detail = c("6,9,13,19", "4726=2359+1210+838+319", "50 per cluster", "4726",
             "4", "NO", basename(sha_after), paste(length(figure_stems), "pairs"),
             "PDF and PNG signatures valid", "no RDS under audit", "NOT_PERFORMED")
)
write_tsv(validation, audit_path("09_qc", "FOUR_CLUSTER_AUDIT_VALIDATION.tsv"))

# Reproducible manifest after all final outputs except the manifest itself exist.
roots <- c("00_manifest", "04_tables", "05_figures", "07_logs", "08_methods", "09_qc", "10_release")
all_files <- unlist(lapply(roots, function(section) list.files(audit_path(section), recursive = TRUE, full.names = TRUE)))
manifest_path <- audit_path("00_manifest", "FOUR_CLUSTER_AUDIT_MANIFEST_SHA256.tsv")
all_files <- sort(setdiff(all_files[file.info(all_files)$isdir == FALSE], manifest_path))
sha_cmd <- system2("sha256sum", all_files, stdout = TRUE)
manifest <- tstrsplit(sha_cmd, "  ", fixed = TRUE)
manifest <- data.table(sha256 = manifest[[1]], path = manifest[[2]])
manifest[, relative_path := sub(paste0("^", project_dir, "/"), "", path)]
manifest[, path := NULL]
write_tsv(manifest, manifest_path)

cat("SUCCESS\n")
