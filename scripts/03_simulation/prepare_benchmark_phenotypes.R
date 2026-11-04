#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: prepare_benchmark_phenotypes.R PHENO_SIM_CSV NEW_OUT_DIR")
library(data.table)
if(dir.exists(args[2])) stop("Output exists")
p <- fread(args[1],select=c("iid",as.character(1:10)))
stopifnot(nrow(p)==31301, !anyDuplicated(p$iid))
p <- p[1:20000]
dir.create(args[2],recursive=TRUE)
dir.create(file.path(args[2],"fastGWA_data"))
for(i in 1:10) fwrite(setNames(data.table(FID=0,IID=p$iid,P=p[[as.character(i)]]),c("FID","IID",paste0("pheno",i))),
  file.path(args[2],"fastGWA_data",paste0("pheno.",i,".txt")),sep="\t",col.names=TRUE)
y <- data.table(FID=0,IID=p$iid)
for(i in 1:10) y[[paste0("pheno",i)]] <- p[[as.character(i)]]
fwrite(y,file.path(args[2],"pheno_combined.txt"),sep="\t")
# Historical edit_phenofile.sh writes FID=0 and keeps the IID unchanged.
fwrite(data.table(IID=p$iid),file.path(args[2],"IID.txt"))
fwrite(y[,-1],file.path(args[2],"pheno_slemm.csv"))
