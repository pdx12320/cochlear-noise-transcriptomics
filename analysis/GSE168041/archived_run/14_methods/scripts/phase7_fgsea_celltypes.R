#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop("Usage: phase7_fgsea_celltypes.R <root>")
root <- normalizePath(args[[1]], mustWork=TRUE)
seed_value <- 20260825L
min_size <- 10L
max_size <- 500L
gsea_dir <- file.path(root, "10_gsea")
dir.create(gsea_dir, recursive=TRUE, showWarnings=FALSE)

safe_name <- function(x) gsub("(^_+|_+$)", "", gsub("[^A-Za-z0-9]+", "_", x))

message("[1/5] Freezing Mouse MSigDB 2026.1.Mm Hallmark + Reactome membership")
collections <- as.data.table(msigdbr_collections(db_species="MM"))
hall <- as.data.table(msigdbr(db_species="MM", species="Mus musculus", collection="MH"))
react <- as.data.table(msigdbr(db_species="MM", species="Mus musculus", collection="M2", subcollection="CP:REACTOME"))
db_versions <- unique(c(hall$db_version, react$db_version))
if (length(db_versions) != 1L || db_versions[[1]] != "2026.1.Mm") stop("MSigDB version mismatch")
term2gene <- unique(rbindlist(list(
  unique(hall[, .(pathway_id=paste0("Hallmark::", gs_name), gene_symbol)]),
  unique(react[, .(pathway_id=paste0("Reactome::", gs_name), gene_symbol)])
)))
term2gene <- term2gene[!is.na(gene_symbol) & gene_symbol != ""]
set_sizes_source <- term2gene[, .(source_size=uniqueN(gene_symbol)), by=pathway_id]
if (uniqueN(term2gene[grepl("^Hallmark::", pathway_id), pathway_id]) != 50L) stop("Hallmark count mismatch")
if (uniqueN(term2gene[grepl("^Reactome::", pathway_id), pathway_id]) != 1333L) stop("Reactome count mismatch")
setorder(term2gene, pathway_id, gene_symbol)
membership_path <- file.path(gsea_dir, "MSIGDB_2026.1.Mm_Hallmark_Reactome_TERM2GENE.tsv.gz")
fwrite(term2gene, membership_path, sep="\t")
membership_sha <- system2("sha256sum", membership_path, stdout=TRUE)
writeLines(membership_sha, file.path(gsea_dir, "MSIGDB_TERM2GENE.sha256"))
fwrite(collections, file.path(gsea_dir, "MSIGDB_COLLECTIONS_SNAPSHOT.tsv"), sep="\t")
pathways <- split(term2gene$gene_symbol, term2gene$pathway_id)

focused_map <- data.table(
  program=c("TNFa/NFkB", "IL6/JAK-STAT3", "VEGF/angiogenesis", "VEGF/angiogenesis", "MAPK/AP-1", "MAPK/AP-1"),
  pathway_id=c(
    "Hallmark::HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "Hallmark::HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "Hallmark::HALLMARK_ANGIOGENESIS",
    "Reactome::REACTOME_SIGNALING_BY_VEGF",
    "Reactome::REACTOME_MAPK_FAMILY_SIGNALING_CASCADES",
    "Reactome::REACTOME_ACTIVATION_OF_THE_AP_1_FAMILY_OF_TRANSCRIPTION_FACTORS"
  ),
  hypothesis_role=c("PRIMARY", "PRIMARY", "PRIMARY", "PRIMARY_CORROBORATIVE", "SECONDARY", "SECONDARY")
)
if (!all(focused_map$pathway_id %in% names(pathways))) stop("Focused pathway missing from frozen membership")
fwrite(focused_map, file.path(gsea_dir, "FOCUSED_PATHWAY_MAPPING_FROZEN.tsv"), sep="\t")

summary <- fread(file.path(root, "09_differential", "S1_DS_SUMMARY.tsv"))
if (nrow(summary) != 13L || any(summary$status != "PASS")) stop("DE summary eligibility mismatch")

all_results <- list()
qc_rows <- list()
focus_results <- list()

message("[2/5] Running cell-type-specific preranked fgseaMultilevel")
for (i in seq_len(nrow(summary))) {
  s <- summary[i]
  tissue <- s$tissue
  cell_type <- s$cell_type
  message(sprintf("GSEA %d/%d: %s / %s", i, nrow(summary), tissue, cell_type))
  de <- fread(s$result_file)
  n_input <- nrow(de)
  valid <- de[!is.na(gene_symbol) & trimws(gene_symbol) != "" & is.finite(stat)]
  valid[, gene_symbol := trimws(gene_symbol)]
  valid[, abs_stat := abs(stat)]
  valid[, duplicate_group_size := .N, by=gene_symbol]
  setorderv(valid, c("gene_symbol", "abs_stat", "stat", "gene_id"), c(1L, -1L, -1L, 1L))
  valid[, duplicate_rank := seq_len(.N), by=gene_symbol]
  duplicate_audit <- valid[duplicate_group_size > 1L,
                           .(gene_id, gene_symbol, stat, abs_stat, duplicate_group_size,
                             duplicate_rank, retained=duplicate_rank == 1L)]
  rank_dt <- valid[duplicate_rank == 1L, .(gene_id, gene_symbol, stat)]
  setorderv(rank_dt, c("stat", "gene_symbol", "gene_id"), c(-1L, 1L, 1L))
  gene_list <- rank_dt$stat
  names(gene_list) <- rank_dt$gene_symbol
  if (anyDuplicated(names(gene_list)) || any(!is.finite(gene_list))) stop("Invalid rank vector")
  overlap_rank_genes <- uniqueN(term2gene[gene_symbol %in% names(gene_list), gene_symbol])
  overlap_fraction <- overlap_rank_genes / length(gene_list)
  overlap_counts <- term2gene[gene_symbol %in% names(gene_list), .(overlap_size=uniqueN(gene_symbol)), by=pathway_id]
  set_sizes <- merge(set_sizes_source, overlap_counts, by="pathway_id", all.x=TRUE)
  set_sizes[is.na(overlap_size), overlap_size := 0L]
  eligible_sets <- set_sizes[overlap_size >= min_size & overlap_size <= max_size]
  if (overlap_rank_genes < 1000L || overlap_fraction < 0.10 || nrow(eligible_sets) < 500L) {
    stop(sprintf("STOP: abnormal gene-ID mapping for %s/%s", tissue, cell_type))
  }
  out_dir <- file.path(gsea_dir, tissue, safe_name(cell_type))
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  fwrite(rank_dt, file.path(out_dir, "RANKING_WALD_STAT.tsv"), sep="\t")
  fwrite(duplicate_audit, file.path(out_dir, "DUPLICATE_SYMBOL_AUDIT.tsv"), sep="\t")
  fwrite(set_sizes, file.path(out_dir, "GENE_SET_OVERLAP_AUDIT.tsv"), sep="\t")
  set.seed(seed_value)
  gse <- fgseaMultilevel(
    pathways=pathways, stats=gene_list, minSize=min_size, maxSize=max_size,
    eps=0, scoreType="std", nproc=1
  )
  out <- as.data.table(gse)
  out <- merge(out, set_sizes, by.x="pathway", by.y="pathway_id", all.x=TRUE)
  setnames(out, c("pathway", "pval", "padj"), c("pathway_id", "P", "FDR"))
  out[, source := sub("::.*$", "", pathway_id)]
  out[, pathway := sub("^[^:]+::", "", pathway_id)]
  out[, leadingEdge := vapply(leadingEdge, paste, collapse=",", FUN.VALUE=character(1))]
  out[, gene_coverage := size / source_size]
  out[, `:=`(tissue=tissue, cell_type=cell_type, direction=ifelse(NES > 0, "Noise", "Naive"))]
  setcolorder(out, c("tissue", "cell_type", "source", "pathway_id", "pathway", "NES", "P", "FDR",
                     "size", "source_size", "gene_coverage", "leadingEdge", "ES", "log2err", "direction"))
  out[, abs_NES_sort := abs(NES)]
  setorderv(out, c("FDR", "abs_NES_sort", "pathway_id"), c(1L, -1L, 1L))
  out[, abs_NES_sort := NULL]
  fwrite(out, file.path(out_dir, "S1_GSEA_ALL.tsv"), sep="\t")
  fwrite(out[FDR <= 0.05], file.path(out_dir, "S1_GSEA_FDR05.tsv"), sep="\t")
  saveRDS(gse, file.path(out_dir, "S1_GSEA_fgsea_object.rds"))
  all_results[[i]] <- out
  focus <- merge(focused_map, out[, .(pathway_id, tissue, cell_type, NES, P, FDR, size, source_size, gene_coverage, leadingEdge, direction)],
                 by="pathway_id", all.x=TRUE)
  focus_results[[i]] <- focus
  qc_rows[[i]] <- data.table(
    tissue=tissue, cell_type=cell_type, de_rows=n_input, valid_symbol_finite_stat=nrow(valid),
    duplicate_symbols=uniqueN(valid[duplicate_group_size > 1L, gene_symbol]),
    unique_rank_genes=length(gene_list), positive_stats=sum(gene_list > 0),
    negative_stats=sum(gene_list < 0), zero_stats=sum(gene_list == 0),
    tied_stat_rows=sum(duplicated(gene_list) | duplicated(gene_list, fromLast=TRUE)),
    rank_genes_matching_any_set=overlap_rank_genes, rank_gene_set_union_overlap_fraction=overlap_fraction,
    eligible_gene_sets=nrow(eligible_sets), tested_gene_sets=nrow(out), fdr05_gene_sets=sum(out$FDR <= 0.05),
    status="PASS"
  )
}

all_out <- rbindlist(all_results, fill=TRUE)
focus_out <- rbindlist(focus_results, fill=TRUE)
qc_out <- rbindlist(qc_rows, fill=TRUE)
fwrite(all_out, file.path(gsea_dir, "S1_GSEA_ALL.tsv"), sep="\t")
fwrite(all_out[FDR <= 0.05], file.path(gsea_dir, "S1_GSEA_FDR05.tsv"), sep="\t")
fwrite(focus_out, file.path(gsea_dir, "S1_GSEA_FOCUSED_PROGRAMS.tsv"), sep="\t")
fwrite(qc_out, file.path(gsea_dir, "S1_GSEA_QA.tsv"), sep="\t")

message("[3/5] Focused lateral-wall program heatmap")
focus_plot <- focus_out[tissue == "Lateral_wall" & cell_type %in% c("Endothelial", "Fibrocyte / mesenchymal", "Immune / myeloid")]
focus_plot[, cell_type := factor(cell_type, levels=c("Endothelial", "Fibrocyte / mesenchymal", "Immune / myeloid"))]
focus_plot[, pathway_label := factor(
  paste0(program, " — ", sub("^[^:]+::", "", pathway_id)),
  levels=rev(unique(paste0(focused_map$program, " — ", sub("^[^:]+::", "", focused_map$pathway_id))))
)]
p <- ggplot(focus_plot, aes(x=cell_type, y=pathway_label, fill=NES)) +
  geom_tile(color="white", linewidth=0.4) +
  geom_text(aes(label=ifelse(is.na(NES), "NT", sprintf("%.2f\nq=%.2g", NES, FDR))), size=2.7) +
  scale_fill_gradient2(low="#3B6FB6", mid="white", high="#C33D3D", midpoint=0, na.value="grey85") +
  labs(title="GSE168041 lateral-wall sample-aware pathway validation", x=NULL, y=NULL, fill="NES\nNoise vs Naive") +
  theme_bw(base_size=9) + theme(axis.text.x=element_text(angle=25, hjust=1), panel.grid=element_blank())
ggsave(file.path(root, "12_figures", "FOCUSED_PROGRAM_NES_HEATMAP.png"), p, width=9.5, height=6.5, dpi=260)
ggsave(file.path(root, "12_figures", "FOCUSED_PROGRAM_NES_HEATMAP.pdf"), p, width=9.5, height=6.5)
fwrite(focus_plot, file.path(root, "13_source_data", "FIG_FOCUSED_PROGRAM_NES_HEATMAP.tsv"), sep="\t")

message("[4/5] Top signed NES plots for prespecified major cell types")
for (ct in c("Endothelial", "Fibrocyte / mesenchymal", "Immune / myeloid")) {
  d <- copy(all_out[tissue == "Lateral_wall" & cell_type == ct])
  pos <- head(d[order(-NES)], 10L)
  neg <- head(d[order(NES)], 10L)
  top <- unique(rbindlist(list(pos, neg)))
  top[, label := factor(pathway, levels=rev(top[order(NES)]$pathway))]
  gp <- ggplot(top, aes(x=NES, y=label, color=FDR <= 0.05, size=-log10(pmax(FDR, 1e-300)))) +
    geom_point() + geom_vline(xintercept=0, linewidth=0.3) +
    scale_color_manual(values=c(`TRUE`="#D95F02", `FALSE`="grey55")) +
    labs(title=paste("Lateral wall", ct, "top signed NES"), x="NES (positive = Noise)", y=NULL,
         color="FDR <= 0.05", size="-log10 FDR") + theme_bw(base_size=8)
  stem <- paste0("GSEA_TOP_NES_Lateral_wall_", safe_name(ct))
  ggsave(file.path(root, "12_figures", paste0(stem, ".png")), gp, width=8.5, height=6.5, dpi=260)
  ggsave(file.path(root, "12_figures", paste0(stem, ".pdf")), gp, width=8.5, height=6.5)
  fwrite(top, file.path(root, "13_source_data", paste0("FIG_", stem, ".tsv")), sep="\t")
}

message("[5/5] Recording versions and parameters")
params <- data.table(
  parameter=c("R", "fgsea", "msigdbr", "MSigDB", "collections", "rank_metric", "duplicate_rule",
              "minSize", "maxSize", "eps", "scoreType", "nproc", "seed", "direction"),
  value=c(R.version.string, as.character(packageVersion("fgsea")), as.character(packageVersion("msigdbr")),
          db_versions[[1]], "MH + M2/CP:REACTOME", "PyDESeq2 Wald stat descending",
          "max abs(stat), then max signed stat, then lexical gene_id", min_size, max_size, 0, "std", 1, seed_value,
          "positive NES = Noise; negative NES = Naive")
)
fwrite(params, file.path(root, "14_methods", "GSEA_PARAMETERS.tsv"), sep="\t")
cat("Phase 7 complete\n")
