#!/usr/bin/env Rscript
# AUTHORIZED 2026-11-01: isolated runtime installation and explicit new reconstruction approved by user.
# Explicit new reconstruction from public counts, not historical 03b/03d recovery.
# Frozen source tables, samples, gene universe and annotations are read-only.

main <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args)!=4L || !args[[1]] %in% c("--preflight-only","--execute-approved")) {
    stop("Usage after approval: Rscript --vanilla proposed_bulk_deseq2_reestimate.R MODE PUBLIC_COUNTS_GZ LOCKED_METADATA_DIRECTORY NEW_OUTPUT_DIRECTORY")
  }
  counts_path <- normalizePath(args[[2]],mustWork=TRUE)
  source_dir <- normalizePath(args[[3]],mustWork=TRUE)
  output_dir <- path.expand(args[[4]])
  if(file.exists(output_dir)||dir.exists(output_dir)) stop("New output directory must not already exist")
  if(as.character(getRversion())!="4.6.1") stop("Requires proposed R 4.6.1; no silent version substitution")
  if(!requireNamespace("DESeq2",quietly=TRUE) || as.character(utils::packageVersion("DESeq2"))!="1.52.0") stop("Requires DESeq2 1.52.0; no dependency installation by this script")
  if(!nzchar(Sys.which("sha256sum"))) stop("Requires existing sha256sum")
  Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
  Sys.setlocale("LC_COLLATE","C")
  options(stringsAsFactors=FALSE,width=160)
  hash <- function(p) {
    z <- system2(Sys.which("sha256sum"),shQuote(p),stdout=TRUE)
    if(!is.null(attr(z,"status"))||length(z)!=1L) stop("Hash failure: ",p)
    sub(" .*","",z[[1]])
  }
  contrasts <- c("94dB_vs_70dB","105dB_vs_70dB","105dB_vs_94dB","overall_dose_LRT")
  manifest_path <- file.path(source_dir,"03b_GSE160639_sample_manifest_locked.csv")
  annotation_path <- file.path(source_dir,"LOCKED_GENE_UNIVERSE_ANNOTATION.tsv")
  inputs <- data.frame(input=c("public_raw_counts","locked_sample_manifest","locked_gene_universe_annotation"),
    path=c(counts_path,manifest_path,annotation_path),
    expected_sha256=c("5f4e4cfa3fa9cf1169fa13acedcf0000961814beea57a3ad94b49dc94eb4df44",
      "9a31de34c1e36e06d70068e72c22498306f7ea662c7a94e541dd3e016ab54164",
      "c26d7bf3aecd2541504bec004b9a55d6f8cbdea87379cfd083f8c7a4a79619f9"))
  inputs$observed_sha256 <- unname(vapply(inputs$path,hash,character(1)))
  if(!all(inputs$expected_sha256==inputs$observed_sha256)) stop("Frozen input hash mismatch")
  annotation <- utils::read.delim(annotation_path,check.names=FALSE)
  if(!identical(names(annotation),c("ensembl_id","gene_symbol","entrez_id"))) stop("Locked annotation schema differs")
  if(nrow(annotation)!=19673L||anyNA(annotation$ensembl_id)||anyDuplicated(annotation$ensembl_id)) stop("Frozen gene universe is invalid")
  con <- gzfile(counts_path,"rt")
  raw <- tryCatch(utils::read.delim(con,check.names=FALSE),finally=close(con))
  sample_order <- c(paste0("A",1:4,"_70dB"),paste0("B",1:4,"_94dB"),"C1_105dB","C2_105dB","C4_105dB","C5_105dB")
  if(nrow(raw)!=35496L || !identical(names(raw),c("gene_id",sample_order)) || anyDuplicated(raw$gene_id)) stop("Public matrix dimensions/IDs/order differ")
  values <- as.matrix(raw[sample_order])
  if(!is.numeric(values)||anyNA(values)||any(!is.finite(values))||any(values<0)||any(values!=floor(values))||any(values>.Machine$integer.max)) stop("Invalid raw integer counts")
  if(!all(annotation$ensembl_id %in% raw$gene_id)) stop("Frozen fit genes absent from public matrix")
  # Preserve the exact historical fitted gene membership; infer no historical prefilter rule.
  count_matrix <- values[match(annotation$ensembl_id,raw$gene_id),,drop=FALSE]
  storage.mode(count_matrix)<-"integer"
  rownames(count_matrix)<-annotation$ensembl_id
  meta <- utils::read.csv(manifest_path,check.names=FALSE)
  if(nrow(meta)!=12L||anyDuplicated(meta$gsm)||!setequal(meta$matrix_column,sample_order)) stop("Locked metadata mismatch")
  meta <- meta[match(sample_order,meta$matrix_column),,drop=FALSE]
  expected_dose <- rep(c(70,94,105),each=4)
  if(!identical(as.numeric(meta$dose_db),expected_dose)) stop("Dose labels do not match locked columns")
  coldata <- data.frame(dose_factor=factor(paste0(meta$dose_db,"dB"),levels=c("70dB","94dB","105dB")),row.names=sample_order)
  if(!identical(colnames(count_matrix),rownames(coldata))) stop("Count/metadata order mismatch")
  if(args[[1]]=="--preflight-only") {
    print(inputs[c("input","observed_sha256")],row.names=FALSE)
    message("PREFLIGHT_PASS: 35496 public genes; 19673 fixed fit genes; twelve samples, four per dose; no normalization/fitting/output creation")
    return(invisible(NULL))
  }
  if(!dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)) stop("Cannot create new output")
  complete <- FALSE
  on.exit(writeLines(if(complete) "COMPLETED_NEW_DESEQ2_RECONSTRUCTION" else "INCOMPLETE_DO_NOT_USE",file.path(output_dir,"STATUS.txt")),add=TRUE)
  write_tsv <- function(x,f) utils::write.table(x,file.path(output_dir,f),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
  write_csv_gz <- function(x,f) {
    con <- gzfile(file.path(output_dir,f),"wt")
    tryCatch(utils::write.csv(x,con,row.names=FALSE,na="NA"),finally=close(con))
  }
  write_tsv(inputs,"INPUT_SHA256.tsv")
  write_tsv(meta,"LOCKED_SAMPLES.tsv")
  write_tsv(annotation,"FROZEN_GENE_UNIVERSE_AND_ANNOTATION.tsv")
  writeLines(c("NEW reconstruction from public raw gene counts with explicit new median-ratio normalization.",
    "Not recovery or exact reproduction of the unavailable historical 03b/03d code, object, size factors or prefilter rule.",
    "Exact historical 19673-gene membership and twelve samples retained; no new gene-level prefilter or sample exclusion.",
    "Gene symbols/Entrez IDs come from the checksum-locked membership/annotation input, independently checked against all four frozen tables; old statistical outputs are not runtime inputs.",
    "GSE168041 models and GSE252285 fixed programs/localization/matched-null results remain unchanged.",
    "These new DE outputs require independent comparison before manuscript adoption and targeted-enrichment chaining."),file.path(output_dir,"README.txt"))
  pars <- data.frame(parameter=c("status","R","DESeq2","gene_universe","design","reference","normalization","fitType","betaPrior","minReplicatesForReplace","useT","minmu","alpha","independentFiltering","cooksCutoff","pAdjustMethod","lfcThreshold","lfcShrink","parallel","LRT_reduced","seed"),
    value=c("NEW_RECONSTRUCTION_NOT_HISTORICAL_REPRODUCTION","4.6.1","1.52.0","exact19673frozenIDs_no_new_prefilter","~dose_factor","70dB","estimateSizeFactors(type=ratio,locfunc=median) on exact19673genes","parametric","FALSE","7","FALSE","0.5","0.05","TRUE","TRUE (DESeq2 default Cook cutoff)","BH","0","not_applied","FALSE","~1","20261101"))
  write_tsv(pars,"PARAMETERS.tsv")
  warnings_seen <- character()
  warning_record <- function(w) warnings_seen <<- c(warnings_seen,conditionMessage(w))
  start <- proc.time()[["elapsed"]]
  set.seed(20261101L)
  dds0 <- DESeq2::DESeqDataSetFromMatrix(countData=count_matrix,colData=coldata,design=~dose_factor)
  # This normalization is newly specified and must not be described as the old size factors.
  dds0 <- withCallingHandlers(DESeq2::estimateSizeFactors(dds0,type="ratio",locfunc=stats::median),warning=warning_record)
  sf <- DESeq2::sizeFactors(dds0)
  if(length(sf)!=12L||any(!is.finite(sf))||any(sf<=0)) stop("Invalid newly estimated size factors")
  write_tsv(data.frame(sample_id=sample_order,dose_db=meta$dose_db,raw_library_size=colSums(count_matrix),size_factor=sf),"NEW_SIZE_FACTORS.tsv")
  saveRDS(dds0,file.path(output_dir,"NEW_UNFITTED_NORMALIZED_DDS.rds"))
  wald <- withCallingHandlers(DESeq2::DESeq(dds0,test="Wald",fitType="parametric",sfType="ratio",betaPrior=FALSE,minReplicatesForReplace=7,parallel=FALSE,useT=FALSE,minmu=0.5,quiet=FALSE),warning=warning_record)
  if(!identical(attr(DESeq2::dispersionFunction(wald),"fitType"),"parametric")) stop("Unexpected dispersion fallback; require review, no silent method switch")
  lrt <- withCallingHandlers(DESeq2::DESeq(dds0,test="LRT",reduced=~1,fitType="parametric",sfType="ratio",betaPrior=FALSE,minReplicatesForReplace=7,parallel=FALSE,minmu=0.5,quiet=FALSE),warning=warning_record)
  if(!identical(attr(DESeq2::dispersionFunction(lrt),"fitType"),"parametric")) stop("Unexpected LRT dispersion fallback")
  comparison_plan <- data.frame(name=contrasts,numerator=c("94dB","105dB","105dB",NA),denominator=c("70dB","70dB","94dB",NA))
  summaries <- list()
  for(i in seq_len(nrow(comparison_plan))) {
    row<-comparison_plan[i,]
    fit<-if(i<=3) wald else lrt
    call <- list(object=fit,alpha=.05,independentFiltering=TRUE,cooksCutoff=TRUE,pAdjustMethod="BH",parallel=FALSE)
    if(i<=3) call$contrast<-c("dose_factor",row$numerator,row$denominator)
    ans <- withCallingHandlers(do.call(DESeq2::results,call),warning=warning_record)
    output <- as.data.frame(ans)
    if(!identical(rownames(output),annotation$ensembl_id)) stop("Output gene order differs")
    output<-cbind(annotation,output)
    output$comparison<-row$name;output$numerator<-row$numerator;output$denominator<-row$denominator
    output$analysis_status<-"NEW_DESEQ2_RECONSTRUCTION"
    write_csv_gz(output,paste0("NEW_DESeq2_",row$name,"_all_genes.csv.gz"))
    sig<-!is.na(output$padj)&output$padj<.05
    summaries[[i]]<-data.frame(comparison=row$name,genes=nrow(output),significant=sum(sig),positive=if(i<=3)sum(sig&output$log2FoldChange>0) else NA,negative=if(i<=3)sum(sig&output$log2FoldChange<0) else NA,pvalue_NA=sum(is.na(output$pvalue)),padj_NA=sum(is.na(output$padj)),filterThreshold=as.numeric(S4Vectors::metadata(ans)$filterThreshold))
  }
  write_tsv(do.call(rbind,summaries),"NEW_DESEQ2_SUMMARY.tsv")
  write_tsv(data.frame(ensembl_id=rownames(wald),as.data.frame(S4Vectors::mcols(wald))),"WALD_GENE_DIAGNOSTICS.tsv")
  write_tsv(data.frame(ensembl_id=rownames(wald),SummarizedExperiment::assays(wald)[["cooks"]],check.names=FALSE),"WALD_COOKS_DISTANCE.tsv")
  normalized<-DESeq2::counts(wald,normalized=TRUE)
  write_tsv(data.frame(ensembl_id=rownames(wald),normalized,check.names=FALSE),"NEW_NORMALIZED_COUNTS.tsv")
  saveRDS(wald,file.path(output_dir,"NEW_FITTED_WALD_DDS.rds"))
  saveRDS(lrt,file.path(output_dir,"NEW_FITTED_LRT_DDS.rds"))
  capture.output(utils::sessionInfo(),file=file.path(output_dir,"SESSIONINFO.txt"))
  writeLines(unique(warnings_seen),file.path(output_dir,"WARNINGS.txt"))
  write_tsv(data.frame(elapsed_seconds=proc.time()[["elapsed"]]-start,requested_processes=1L),"RUNTIME.tsv")
  after<-unname(vapply(inputs$path,hash,character(1)))
  if(!identical(after,inputs$observed_sha256)) stop("An upstream input changed during run")
  complete<-TRUE
  message("Completed new DESeq2 reconstruction; originals untouched. Independent review and new targeted-GSEA chaining remain required.")
}
main()
