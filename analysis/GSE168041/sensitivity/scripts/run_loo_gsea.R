#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(fgsea)})
a<-commandArgs(TRUE);stopifnot(length(a)==4L,as.character(getRversion())=='4.5.3',as.character(packageVersion('fgsea'))=='1.36.2')
root<-normalizePath(a[1],mustWork=TRUE);tissue<-a[2];cell<-a[3];omit<-a[4]
src<-file.path(root,'inputs',tissue,cell);fit<-file.path(root,'fits',tissue,cell,omit);out<-file.path(fit,'gsea')
stopifnot(!dir.exists(out),readLines(file.path(fit,'STATUS.txt'))=='COMPLETED_LOO_DE')
stopifnot(identical(formals(fgseaMultilevel)$sampleSize,101),identical(formals(fgseaMultilevel)$nPermSimple,1000),identical(formals(fgseaMultilevel)$gseaParam,1))
sha<-function(p)sub(' .*','',system2('sha256sum',shQuote(p),stdout=TRUE)[1])
mf<-fread(file.path(root,'FROZEN_INPUT_SHA256.tsv'))
stopifnot(all(vapply(file.path(root,'inputs',mf$file),sha,character(1))==mf$sha256))
read_gz_tsv<-function(path) { con<-gzfile(path,'rt');on.exit(close(con));as.data.table(read.delim(con,check.names=FALSE)) }
term<-read_gz_tsv(file.path(root,'inputs','MSIGDB_2026.1.Mm_Hallmark_Reactome_TERM2GENE.tsv.gz'))
family<-fread(file.path(src,'FROZEN_TESTING_FAMILY.tsv'));setnames(family,'size','baseline_overlap_size')
de<-read_gz_tsv(file.path(fit,'LOO_DE_ALL.tsv.gz'));valid<-de[!is.na(gene_symbol)&trimws(gene_symbol)!=''&is.finite(stat)]
valid[,gene_symbol:=trimws(gene_symbol)];valid[,abs_stat:=abs(stat)];valid[,duplicate_group_size:=.N,by=gene_symbol]
setorderv(valid,c('gene_symbol','abs_stat','stat','gene_id'),c(1L,-1L,-1L,1L));valid[,duplicate_rank:=seq_len(.N),by=gene_symbol]
ranks<-valid[duplicate_rank==1L,.(gene_id,gene_symbol,stat)];setorderv(ranks,c('stat','gene_symbol','gene_id'),c(-1L,1L,1L))
v<-ranks$stat;names(v)<-ranks$gene_symbol;stopifnot(!anyDuplicated(names(v)),all(is.finite(v)))
pathways<-split(term[order(pathway_id,gene_symbol)]$gene_symbol,term[order(pathway_id,gene_symbol)]$pathway_id)
pathways<-pathways[family$pathway_id];stopifnot(!any(vapply(pathways,is.null,logical(1))))
family[,overlap_size:=vapply(pathways,function(x)length(intersect(x,names(v))),integer(1))]
family[,eligible_after_omission:=overlap_size>=10L&overlap_size<=500L]
dir.create(out);writeLines('INCOMPLETE_DO_NOT_USE',file.path(out,'STATUS.txt'))
fwrite(ranks,file.path(out,'RANKING_WALD_STAT.tsv.gz'),sep='\t')
fwrite(valid[duplicate_group_size>1L,.(gene_id,gene_symbol,stat,duplicate_group_size,duplicate_rank)],file.path(out,'DUPLICATE_SYMBOL_AUDIT.tsv'),sep='\t')
fwrite(family,file.path(out,'FIXED_FAMILY_OVERLAP_AUDIT.tsv'),sep='\t')
set.seed(20260825L);start<-proc.time()[['elapsed']];warnings_seen<-character()
gse<-withCallingHandlers(fgseaMultilevel(pathways=pathways,stats=v,minSize=10,maxSize=500,eps=0,scoreType='std',nproc=1,gseaParam=1,sampleSize=101,nPermSimple=1000),warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')})
gse<-as.data.table(gse);setnames(gse,c('pathway','pval','padj'),c('pathway_id','P','FDR_returned'))
gse[,leadingEdge:=vapply(leadingEdge,paste,collapse=',',FUN.VALUE=character(1))]
res<-merge(family,gse,by='pathway_id',all.x=TRUE,sort=FALSE)
res[,FDR:=NA_real_];res[!is.na(P),FDR:=p.adjust(P,method='BH',n=nrow(family))]
res[,`:=`(tissue=tissue,safe_cell_type=cell,omitted_sample_id=omit,original_family_size=nrow(family),direction=ifelse(is.na(NES),'UNAVAILABLE',ifelse(NES>0,'Noise','Naive')),test_status=ifelse(is.na(P),'UNAVAILABLE','TESTED'))]
stopifnot(nrow(res)==nrow(family),all(is.na(res$P[!res$eligible_after_omission])))
fwrite(res,file.path(out,'LOO_GSEA_ALL.tsv.gz'),sep='\t')
fwrite(res[pathway_id %in% c('Hallmark::HALLMARK_TNFA_SIGNALING_VIA_NFKB','Reactome::REACTOME_MAPK_FAMILY_SIGNALING_CASCADES')],file.path(out,'FOCUSED_RESULTS.tsv'),sep='\t')
writeLines(unique(warnings_seen),file.path(out,'WARNINGS.txt'))
capture.output(sessionInfo(),file=file.path(out,'SESSIONINFO.txt'))
fwrite(data.table(tissue=tissue,safe_cell_type=cell,omitted_sample_id=omit,original_family_size=nrow(family),tested=sum(!is.na(res$P)),unavailable=sum(is.na(res$P)),ranked_symbols=length(v),tied_stat_rows=sum(duplicated(v)|duplicated(v,fromLast=TRUE)),elapsed_seconds=proc.time()[['elapsed']]-start,seed=20260825L,nproc=1),file.path(out,'GSEA_RECORD.tsv'),sep='\t')
writeLines('COMPLETED_LOO_GSEA',file.path(out,'STATUS.txt'));print(res[pathway_id=='Hallmark::HALLMARK_TNFA_SIGNALING_VIA_NFKB',.(tissue,safe_cell_type,omitted_sample_id,NES,P,FDR,original_family_size)])
