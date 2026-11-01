#!/usr/bin/env Rscript
# Create publication copies only; do not fit models or alter original RDS files.
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(TRUE)
stopifnot(length(args)==2L, !file.exists(args[2]), as.character(packageVersion("DESeq2"))=="1.52.0")
dir.create(args[2],recursive=TRUE)
hash <- function(p) sub(" .*","",system2("sha256sum",shQuote(p),stdout=TRUE)[1])
report <- list()
for (name in c("NEW_UNFITTED_NORMALIZED_DDS.rds","NEW_FITTED_WALD_DDS.rds","NEW_FITTED_LRT_DDS.rds")) {
  source <- file.path(args[1],name)
  x <- readRDS(source); y <- x
  f <- design(y); environment(f) <- baseenv(); slot(y,"design") <- f
  d <- dispersionFunction(y)
  fitted <- grepl("FITTED_(WALD|LRT)",name)
  if (fitted) {
    stopifnot(identical(attr(d,"fitType"),"parametric"),identical(paste(deparse(body(d)),collapse=""),"coefs[1] + coefs[2]/q"))
    e <- new.env(parent=baseenv()); e$coefs <- get("coefs",environment(d),inherits=FALSE)
    environment(d) <- e
  } else {
    stopifnot(identical(body(d),NULL))
    environment(d) <- baseenv()
  }
  slot(y,"dispersionFunction") <- d
  unchanged <- setdiff(slotNames(x),c("design","dispersionFunction"))
  stopifnot(all(vapply(unchanged,function(s)identical(slot(x,s),slot(y,s)),logical(1))))
  extra <- setdiff(names(attributes(x)),slotNames(x))
  stopifnot(identical(attributes(x)[extra],attributes(y)[extra]))
  stopifnot(identical(as.character(design(x)),as.character(design(y))),identical(attributes(dispersionFunction(x)),attributes(d)))
  stopifnot(identical(counts(x),counts(y)),identical(counts(x,normalized=TRUE),counts(y,normalized=TRUE)))
  disp_error <- 0; result_error <- 0
  if (fitted) {
    means <- mcols(x)$baseMean
    stopifnot(identical(dispersionFunction(x)(means),dispersionFunction(y)(means)))
    plans <- if (grepl("WALD",name)) list(c("dose_factor","94dB","70dB"),c("dose_factor","105dB","70dB"),c("dose_factor","105dB","94dB")) else list(NULL)
    for (contrast in plans) {
      a <- list(object=x,alpha=.05,independentFiltering=TRUE,cooksCutoff=TRUE,pAdjustMethod="BH",parallel=FALSE)
      if(!is.null(contrast)) a$contrast <- contrast
      b <- a; b$object <- y
      rx <- do.call(results,a); ry <- do.call(results,b)
      stopifnot(identical(as.data.frame(rx),as.data.frame(ry)),identical(metadata(rx),metadata(ry)))
    }
  }
  validObject(y)
  out <- file.path(args[2],name); saveRDS(y,out)
  z <- readRDS(out)
  stopifnot(all(vapply(unchanged,function(s)identical(slot(x,s),slot(z,s)),logical(1))))
  # XDR serialized bytes are searched for the private machine's root and account.
  blob <- serialize(z,NULL,ascii=TRUE)
  text <- rawToChar(blob)
  stopifnot(!grepl("/data/|/home/|/Users/",text))
  report[[name]] <- data.frame(file=name,source_sha256=hash(source),public_sha256=hash(out),source_bytes=file.info(source)$size,public_bytes=file.info(out)$size,unchanged_slots=paste(unchanged,collapse=","),all_assays_row_col_metadata_identical=TRUE,normalized_counts_identical=TRUE,design_expression_identical=TRUE,dispersion_attributes_identical=TRUE,dispersion_outputs_identical=if(fitted)TRUE else NA,regenerated_results_and_metadata_identical=if(fitted)TRUE else NA,serialized_private_path_scan_pass=TRUE)
}
write.table(do.call(rbind,report),file.path(args[2],"PUBLIC_OBJECT_VALIDATION.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c("Publication copies of the current R4.6.1/DESeq2 1.52.0 analysis objects.","Only design formula environments and dispersion-function environments were replaced.","Design expressions, function bodies/attributes, every other slot, assays, row/column metadata and extra object attributes are identical.","The design environment now points to baseenv; fitted parametric dispersion closures retain only their original coefficients with baseenv parent; the unfitted empty closure uses baseenv.","Dispersion predictions at every fitted gene baseMean and all three Wald/one LRT regenerated results and result metadata were identical.","This removes captured machine-local execution arguments and redundant objects. Original execution RDS files remain unchanged.","PUBLIC_OBJECT_VALIDATION.tsv records source/public SHA-256 and checks. sanitize_public_dds.R documents the exact transformation."),file.path(args[2],"README.txt"))
message("PUBLIC_OBJECT_VALIDATION_PASS")
