#!/usr/bin/env Rscript
# Portable optional restore helper. Not used to create the original runtime.
# Inputs: source-tarball directory, SHA manifest, NEW R library directory.
# Base R only; no repository access or package-version substitution.
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=3L) stop("Usage: Rscript --vanilla restore_locked_r_packages.R PACKAGE_SOURCES PACKAGE_SOURCE_SHA256.tsv NEW_LIBRARY")
stopifnot(as.character(getRversion())=="4.6.1")
.libPaths(.Library,include.site=FALSE)
sources<-normalizePath(a[[1]],mustWork=TRUE)
manifest<-utils::read.delim(a[[2]],check.names=FALSE,stringsAsFactors=FALSE)
lib<-path.expand(a[[3]])
if(dir.exists(lib)||file.exists(lib)) stop("Refusing to replace an existing library")
stopifnot(identical(names(manifest),c("file","sha256")),!anyDuplicated(manifest$file))
sha<-function(p) {
 z<-system2(Sys.which("sha256sum"),shQuote(p),stdout=TRUE)
 if(!is.null(attr(z,"status"))) stop("Hash command failed")
 sub(" .*","",z[[1]])
}
if(!nzchar(Sys.which("sha256sum"))) stop("Requires sha256sum")
files<-file.path(sources,manifest$file)
if(!all(file.exists(files))||!all(unname(vapply(files,sha,character(1)))==manifest$sha256)) stop("Source archive hash mismatch")
read_description<-function(f) {
 entries<-utils::untar(f,list=TRUE)
 target<-entries[grepl("^[^/]+/DESCRIPTION$",entries)]
 if(length(target)!=1L) stop("Expected one package DESCRIPTION")
 tmp<-tempfile("locked_r_metadata_");dir.create(tmp)
 on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
 utils::untar(f,files=target,exdir=tmp)
 as.list(read.dcf(file.path(tmp,target))[1,])
}
desc<-lapply(files,read_description)
names(desc)<-vapply(desc,function(x)x$Package,character(1))
if(anyDuplicated(names(desc))) stop("More than one source archive for a package")
deps<-lapply(desc,function(x) {
 fields<-unlist(x[intersect(c("Depends","Imports","LinkingTo"),names(x))],use.names=FALSE)
 parts<-unlist(strsplit(paste(fields,collapse=","),",",fixed=TRUE))
 unique(setdiff(trimws(sub("\\s*\\(.*\\)$","",parts)),c("","R")))
})
pending<-names(desc);order<-character()
while(length(pending)) {
 ready<-pending[vapply(pending,function(x)!any(deps[[x]] %in% pending),logical(1))]
 if(!length(ready)) stop("Dependency cycle or malformed source metadata")
 order<-c(order,sort(ready));pending<-setdiff(pending,ready)
}
outside<-setdiff(unique(unlist(deps)),names(desc))
if(length(outside)) {
 missing<-outside[!vapply(outside,requireNamespace,logical(1),quietly=TRUE)]
 if(length(missing)) stop("Required base/runtime dependencies absent: ",paste(missing,collapse=","))
}
dir.create(lib,recursive=TRUE)
.libPaths(c(lib,.Library),include.site=FALSE)
Sys.setenv(MAKEFLAGS="-j4",OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
for(pkg in order) {
 f<-files[[match(pkg,names(desc))]]
 utils::install.packages(f,lib=lib,repos=NULL,type="source",dependencies=FALSE,Ncpus=1)
 if(!requireNamespace(pkg,quietly=TRUE)||as.character(utils::packageVersion(pkg,lib.loc=lib))!=desc[[pkg]]$Version) stop("Exact package restore failed: ",pkg)
}
capture.output(sessionInfo(),file=file.path(lib,"RESTORE_SESSIONINFO.txt"))
utils::write.table(data.frame(package=order,version=vapply(desc[order],function(x)x$Version,character(1))),file.path(lib,"RESTORED_PACKAGE_ORDER.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
message("Exact source package restore completed; analysis has not been run")
