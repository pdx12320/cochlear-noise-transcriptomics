options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(Matrix);library(data.table);library(digest);library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: script.R ORIGINAL_OBJECT_RDS EXISTING_UCELL_SCORES_TSV_GZ NEW_OUTPUT_DIRECTORY")
obj_path <- normalizePath(args[1],mustWork=TRUE)
score_path <- normalizePath(args[2],mustWork=TRUE)
out <- args[3]
if(file.exists(out)) stop("Output directory already exists")
dir.create(out,recursive=TRUE);dir.create(file.path(out,"figure_inputs"))
write_tsv <- function(x,n) fwrite(x,file.path(out,n),sep="\t",quote=FALSE,na="NA")
source_sha <- digest(obj_path,algo="sha256",file=TRUE)
stopifnot(source_sha=="24909e74b584dcb33a113d815234cda0861ed6fea07f78cf66fd99d5324d76b8")
obj <- readRDS(obj_path)
stopifnot(ncol(obj)==36111L,nrow(obj)==32285L)
original <- as.character(obj$broad_annotation_no_doublet_sensitivity)
original_digest <- digest(original,algo="sha256")
clusters <- as.character(obj$seurat_clusters_no_doublet)
stopifnot(identical(clusters,as.character(obj[["RNA_snn_res.0.6",drop=TRUE]])))
c19 <- clusters=="19"
stopifnot(sum(c19)==319L,all(original[c19]=="UNRESOLVED"))
revised <- original;revised[c19] <- "Fibrocyte_mesenchymal"
stopifnot(sum(revised!=original)==319L,all(revised[!c19]==original[!c19]))
suffix <- sub(".*(-[123])$","\\1",colnames(obj))
condition <- c("-3"="Control","-1"="HHL","-2"="NIHL")[suffix]
stopifnot(!anyNA(condition))
meta <- data.table(cell=colnames(obj),cluster=clusters,original_broad=original,c19_revised_broad=revised,
                   barcode_suffix=suffix,condition=unname(condition))
write_tsv(meta,"ANNOTATION_SENSITIVITY_METADATA.tsv")
counts <- LayerData(obj[["RNA"]],layer="counts")
norm <- LayerData(obj[["RNA"]],layer="data")
genes <- c("Dusp1","Fos","Nfkbia","Arrdc2","Apold1","Sik1","Per1","Ccn1","Cebpb","Fosb","Rasd1","Gadd45g","Depp1","Serpine1","Errfi1","Nr4a1","Sgk1","Rgs1","Egr1","Jun")
# Panel order is checked against the executed condition table by the downstream review.
additional <- c("Coch","Otor","Pou3f4","Col1a2","Pdgfra","Pecam1","Kdr","Vwf")
stopifnot(all(c(genes,additional)%in%rownames(counts)))
summarize_expression <- function(labels,variant) {
  rbindlist(lapply(c("ALL","Control","HHL","NIHL"),function(cond) {
    rbindlist(lapply(sort(unique(labels)),function(ct) {
      ii <- which(labels==ct & (cond=="ALL" | meta$condition==cond))
      x <- counts[genes,ii,drop=FALSE];y <- norm[genes,ii,drop=FALSE]
      data.table(variant=variant,condition=cond,broad_cell_type=ct,gene=genes,n_cells=length(ii),
        n_expressing=as.integer(Matrix::rowSums(x>0)),pct_expressing=100*Matrix::rowMeans(x>0),
        mean_raw_count=as.numeric(Matrix::rowMeans(x)),mean_log_normalized=as.numeric(Matrix::rowMeans(y)))
    }))
  }))
}
expr <- rbindlist(list(summarize_expression(original,"original"),summarize_expression(revised,"c19_fibrocyte")))
expr[,rank_mean:=frank(-mean_log_normalized,ties.method="min"),by=.(variant,condition,gene)]
expr[,rank_pct:=frank(-pct_expressing,ties.method="min"),by=.(variant,condition,gene)]
write_tsv(expr,"TOP20_EXPRESSION_BOTH_ANNOTATIONS.tsv")
write_tsv(expr[rank_mean==1|rank_pct==1],"TOP20_TOP_CELLTYPE_MEAN_OR_PCT.tsv")
cell_counts <- rbindlist(list(meta[,.(n_cells=.N),by=.(condition,broad_cell_type=original_broad)][,variant:="original"],
  meta[,.(n_cells=.N),by=.(condition,broad_cell_type=c19_revised_broad)][,variant:="c19_fibrocyte"]))
write_tsv(cell_counts,"CELL_COUNTS_BY_CONDITION.tsv")
cg <- unique(c(genes,additional))
cluster_marker <- rbindlist(lapply(c("ALL","Control","HHL","NIHL"),function(cond){
 ii<-which(c19 & (cond=="ALL"|meta$condition==cond))
 data.table(cluster="19",condition=cond,gene=cg,n_cells=length(ii),n_expressing=as.integer(Matrix::rowSums(counts[cg,ii,drop=FALSE]>0)),
  pct_expressing=100*Matrix::rowMeans(counts[cg,ii,drop=FALSE]>0),mean_log_normalized=as.numeric(Matrix::rowMeans(norm[cg,ii,drop=FALSE])))
}))
write_tsv(cluster_marker,"CLUSTER19_MARKER_EXPRESSION.tsv")
scores <- fread(score_path)
stopifnot(nrow(scores)==nrow(meta),uniqueN(scores$cell)==nrow(meta))
scores <- scores[match(meta$cell,cell)]
stopifnot(identical(scores$cell,meta$cell),identical(as.character(scores$broad_cell_type),original),
          identical(as.character(scores$condition_descriptive),meta$condition))
programs <- c("MAPK_AP1","TNFA_NFKB","IL6_JAK_STAT3","VEGF_ANGIOGENESIS","MYC_TARGETS","E2F_TARGETS")
stopifnot(all(programs%in%names(scores)))
score_hash <- digest(as.matrix(scores[,..programs]),algo="sha256")
ucell_summary <- rbindlist(lapply(c("original","c19_fibrocyte"),function(v){
 labels<-if(v=="original")original else revised
 rbindlist(lapply(c("ALL","Control","HHL","NIHL"),function(cond){
   rbindlist(lapply(sort(unique(labels)),function(ct){
    ii<-which(labels==ct & (cond=="ALL"|meta$condition==cond))
    rbindlist(lapply(programs,function(pr){z<-scores[[pr]][ii];data.table(variant=v,condition=cond,broad_cell_type=ct,
      program=pr,n_cells=length(ii),median_UCell=median(z),mean_UCell=mean(z),Q1_UCell=quantile(z,.25,names=FALSE),Q3_UCell=quantile(z,.75,names=FALSE))}))
   }))
 }))
}))
ucell_summary[,median_rank:=frank(-median_UCell,ties.method="min"),by=.(variant,condition,program)]
write_tsv(ucell_summary,"UCELL_SUMMARY_BOTH_ANNOTATIONS.tsv")
embed <- Embeddings(obj,"tsne")
stopifnot(identical(rownames(embed),meta$cell))
figure <- copy(meta);figure[,`:=`(tSNE_1=embed[,1],tSNE_2=embed[,2])]
for(g in cg) set(figure,j=paste0("expr__",g),value=as.numeric(norm[g,]))
for(pr in programs) set(figure,j=paste0("UCell__",pr),value=scores[[pr]])
fwrite(figure,file.path(out,"figure_inputs","CELL_TSNE_ANNOTATION_EXPRESSION_UCELL.tsv.gz"),sep="\t",quote=FALSE,compress="gzip")
stopifnot(identical(original,as.character(obj$broad_annotation_no_doublet_sensitivity)),
 digest(original,algo="sha256")==original_digest,digest(as.matrix(scores[,..programs]),algo="sha256")==score_hash)
qa <- list(source_object_sha256=source_sha,source_object_written=FALSE,original_annotation_unchanged=TRUE,
 cells=36111L,changed_cells=319L,changed_cluster="19",change="UNRESOLVED -> Fibrocyte_mesenchymal",
 revised_fibrocyte_cells=sum(revised=="Fibrocyte_mesenchymal"),revised_unresolved_cells=sum(revised=="UNRESOLVED"),
 score_cache_sha256=digest(score_path,algo="sha256",file=TRUE),cached_scores_recomputed=FALSE,
 scores_unchanged=TRUE,normalized_expression_layer="RNA data",detection_layer="RNA counts > 0",
 figure_input_columns=names(figure),R=as.character(getRversion()))
write_json(qa,file.path(out,"VALIDATION.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
message("Completed metadata sensitivity and figure inputs")
