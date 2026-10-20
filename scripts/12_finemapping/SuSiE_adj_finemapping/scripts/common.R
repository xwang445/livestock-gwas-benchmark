suppressPackageStartupMessages({library(data.table); library(susieR); library(ggplot2)})
setDTthreads(4)
project <- normalizePath("..")
geno_source <- Sys.getenv("GENO_SOURCE")
regenie_geno <- Sys.getenv("REGENIE_GENO_SOURCE")
snp_info_source <- Sys.getenv("SNP_INFO_SOURCE")
read <- function(p, ...) fread(p, showProgress=FALSE, ...)
write <- function(d,p) fwrite(d,p,sep="\t",na="NA")
read_bim <- function(p) read(paste0(p,".bim"),header=FALSE,
  col.names=c("chr","SNP","cm","position","A1","A2"))
read_fam <- function(p) read(paste0(p,".fam"),header=FALSE,colClasses="character",
  col.names=c("FID","IID","PAT","MAT","SEX","PHENO"))
params <- list(h2=0.5,L=5L,coverage=0.95,min_abs_corr=0.5,max_iter=1000L,
               tol=0.001,seed=20260929L,permutations=10000L,half_window=1000000L)
stopifnot(as.character(packageVersion("susieR")) == "0.14.2")

stopifnot(nzchar(geno_source), nzchar(regenie_geno), nzchar(snp_info_source))
