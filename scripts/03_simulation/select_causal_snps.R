#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: select_causal_snps.R QC_CHIP_SNPLIST CHROM_SIZES NEW_OUT_DIR")
library(data.table)
if(dir.exists(args[3])) stop("Output directory already exists")
s <- fread(args[1], header=FALSE)
sizes <- fread(args[2], header=FALSE, col.names=c("chr","length"))
sizes[, chr := as.integer(sub("chr", "", chr))]
sizes <- sizes[chr %in% 1:18]
stopifnot(nrow(sizes)==18, !anyDuplicated(sizes$chr), !anyDuplicated(s$V1))
s[, c("chr","pos") := tstrsplit(V1,"_",keep=1:2)]
s[, `:=`(chr=as.integer(chr),pos=as.integer(pos))]
s <- s[chr %in% 1:18]
eligible <- s[pos %% 6000000 < 2000000]
stopifnot(nrow(s)==31062, nrow(eligible)==10461)
# R >=3.6 uses Rejection sampling; historical validation used R4.1.2.
RNGkind("Mersenne-Twister","Inversion","Rejection")
set.seed(1)
causal <- sample(eligible$V1,10000)
regions <- rbindlist(lapply(1:18,function(ch) {
  len <- sizes[chr==ch,length];start <- seq(0,len,by=6000000)
  data.table(chr=ch,start=start,end=pmin(start+2000000,len))
}))[start<end]
dir.create(args[3],recursive=TRUE)
fwrite(data.table(SNP=causal),file.path(args[3],"snp.info.chr1_18.10k.csv"))
fwrite(regions,file.path(args[3],"eligible_regions.csv"))
fwrite(sizes,file.path(args[3],"chr1_18_lengths.csv"))
