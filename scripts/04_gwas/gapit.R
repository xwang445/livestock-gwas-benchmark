#!/usr/bin/env Rscript
# Runs both historically logged models; Kansas/NYC are emitted by GAPIT itself.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5L) stop("Usage: gapit.R GAPIT_SOURCE PHENO_TSV PLINK_RAW GM_TSV NEW_OUT_DIR")
library(data.table)
source_file <- normalizePath(args[1],mustWork=TRUE)
expected_md5 <- "d20fd577f8634f508c1bfb8313c8dbad"
if(unname(tools::md5sum(source_file))!=expected_md5) stop("GAPIT source differs from audited retained library; review version provenance")
# The retained source auto-installs missing dependencies. Require them beforehand.
needed <- c("gplots","genetics","ape","compiler","grid","bigmemory","EMMREML","scatterplot3d","lme4","multtest","snpStats")
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Install dependencies explicitly first: ",paste(missing,collapse=", "))
y <- fread(args[2],data.table=FALSE);raw <- fread(args[3],data.table=FALSE)
gm <- fread(args[4],data.table=FALSE)
stopifnot(all(c("FID","IID",paste0("pheno",1:10)) %in% names(y)),!anyDuplicated(y$IID),nrow(y)==20000)
y <- y[,c("IID",paste0("pheno",1:10))];names(y)[1] <- "Taxa"
gd <- data.frame(Taxa=raw$IID,raw[,-(1:6)],check.names=FALSE)
idx <- match(names(gd)[-1],gm$SNP)
stopifnot(!anyNA(idx),!anyDuplicated(gm$SNP),!anyDuplicated(gd$Taxa),all(y$Taxa %in% gd$Taxa))
gm <- gm[idx,];stopifnot(identical(names(gd)[-1],gm$SNP),nrow(gm)==135358)
if(dir.exists(args[5])) stop("Output directory exists")
dir.create(args[5],recursive=TRUE)
source(source_file)
# GAPIT emits fixed filenames in its current directory; isolate all writes here.
previous <- getwd();setwd(normalizePath(args[5]))
tryCatch(GAPIT(Y=y,GD=gd,GM=gm,model=c("FarmCPU","BLINK"),Multiple_analysis=TRUE),finally=setwd(previous))
