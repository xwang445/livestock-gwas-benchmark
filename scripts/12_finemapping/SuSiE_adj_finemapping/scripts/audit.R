source("scripts/common.R")
sink("qc/input_audit.txt",split=TRUE)
cat("Input audit:",format(Sys.time()),"\n")
for (folder in c("BOLTLMM","REGENIE")) {
  cat("\nDirectory:",folder,"\n"); print(list.files(file.path(project,folder)))
}
bim <- read_bim(geno_source); fam <- read_fam(geno_source)
stopifnot(nrow(bim)==135358,nrow(fam)==31301,!anyDuplicated(bim$SNP),!anyDuplicated(fam$IID))
ids <- read(file.path(project,"simu_pheno/IID.txt"),colClasses="character")
stopifnot(identical(names(ids),"IID"),nrow(ids)==20000,!anyDuplicated(ids$IID),all(ids$IID %in% fam$IID))
keep <- fam[match(ids$IID,IID),.(FID,IID)]
fwrite(keep,"genotype/keep.txt",sep="\t",col.names=FALSE)
cat("IID.txt: IID only, header present; exact unique samples:",nrow(keep),"\n")
ph <- read(file.path(project,"simu_pheno/pheno_combined.txt"))
stopifnot(setequal(as.character(ph$IID),ids$IID),nrow(ph)==20000,
          all(vapply(ph[,paste0("pheno",1:10),with=FALSE],function(x) all(is.finite(x)),logical(1))))
rb <- read_bim(regenie_geno); rf <- read_fam(regenie_geno)
stopifnot(identical(fam,rf),nrow(rb)==nrow(bim),!anyDuplicated(rb$SNP),
          identical(bim[,.(chr,position,A1,A2)],rb[,.(chr,position,A1,A2)]))
hash <- function(p) strsplit(system2("sha256sum",shQuote(p),stdout=TRUE)," ")[[1]][1]
bh <- hash(paste0(geno_source,".bed")); rh <- hash(paste0(regenie_geno,".bed"))
stopifnot(identical(bh,rh))
cat("BOLT and REGENIE genotypes identical BED SHA256:",bh,"\n")
info <- read(snp_info_source)
stopifnot("SNP" %in% names(info),!anyDuplicated(info$SNP))
info_ok <- setequal(info$SNP,bim$SNP) && nrow(info)==nrow(bim)
cat("External SNP info matches full panel:",info_ok,"; rows:",nrow(info),"\n")
# Always derive the actual all-marker GRM list directly from the verified BIM.
fwrite(bim[,.(SNP)],"grm/snp_info.csv")
bim[, regenie_ID := rb$SNP]
write(bim,"qc/variants.tsv")
write(fam[IID %in% ids$IID],"qc/expected_samples.tsv")
saveRDS(params,"qc/parameters.rds")
stats <- list(); all_loci <- list()
for (i in 1:10) for (method in c("BOLT-LMM","REGENIE-LOCO")) {
  phenotype <- paste0("pheno",i)
  path <- if(method=="BOLT-LMM") file.path(project,"BOLTLMM",paste0("bolt_",phenotype,".stats")) else
    file.path(project,"REGENIE",paste0("regenie_",phenotype,"_",phenotype,".regenie"))
  x <- read(path)
  if(i==1) {cat("\n",path,"\n");print(head(x,3))}
  if(method=="BOLT-LMM") {
    stopifnot(all(c("SNP","CHR","BP","ALLELE1","ALLELE0","BETA","SE","P_BOLT_LMM_INF") %in% names(x)),!anyDuplicated(x$SNP))
    idx <- match(bim$SNP,x$SNP); stopifnot(!anyNA(idx)); x<-x[idx]
    d<-data.table(SNP=bim$SNP,chr=x$CHR,position=x$BP,effect_allele=x$ALLELE1,
                  other_allele=x$ALLELE0,beta=x$BETA,SE=x$SE,P=x$P_BOLT_LMM_INF,
                  logP=-log10(x$P_BOLT_LMM_INF),frequency=x$A1FREQ)
  } else {
    stopifnot(all(c("ID","CHROM","GENPOS","ALLELE1","ALLELE0","BETA","SE","LOG10P","CHISQ","N","TEST") %in% names(x)))
    x<-x[TEST=="ADD"];stopifnot(!anyDuplicated(x$ID),all(x$N==20000))
    idx<-match(bim$regenie_ID,x$ID);stopifnot(!anyNA(idx));x<-x[idx]
    d<-data.table(SNP=bim$SNP,chr=x$CHROM,position=x$GENPOS,effect_allele=x$ALLELE1,
                  other_allele=x$ALLELE0,beta=x$BETA,SE=x$SE,P=10^(-x$LOG10P),
                  logP=x$LOG10P,frequency=x$A1FREQ)
  }
  stopifnot(nrow(x)==nrow(bim),all(d$chr==bim$chr & d$position==bim$position))
  same<-d$effect_allele==bim$A1 & d$other_allele==bim$A2
  flip<-d$effect_allele==bim$A2 & d$other_allele==bim$A1
  d[,reason := fifelse(!same & !flip,"allele_mismatch",
    fifelse(!is.finite(beta)|!is.finite(SE)|SE<=0,"invalid_beta_SE",
    fifelse(is.na(logP)|logP<0,"invalid_P","OK")))]
  d[, `:=`(z=beta/SE,allele_flip=flip,
            ambiguous_pair=paste0(effect_allele,other_allele) %in% c("AT","TA","CG","GC"))]
  # Same BED and explicit allele labels establish strand, including palindromes.
  expected_logp <- -pchisq(d$z^2,df=1,lower.tail=FALSE,log.p=TRUE)/log(10)
  valid<-d$reason=="OK" & is.finite(expected_logp) & is.finite(d$logP)
  delta<-max(abs(expected_logp[valid]-d$logP[valid]))
  stopifnot(delta<0.1,all(sign(d$z[valid])==sign(d$beta[valid])))
  if(method=="REGENIE-LOCO") stopifnot(max(abs(d$z[valid]^2-x$CHISQ[valid])/pmax(1,x$CHISQ[valid]))<0.001)
  d[, `:=`(phenotype=phenotype,GWAS_method=method)]
  write(d,sprintf("qc/%s_%s.sumstats.tsv",method,phenotype))
  stats[[length(stats)+1]]<-data.table(phenotype,method,n=nrow(d),invalid=sum(d$reason!="OK"),
    flips=sum(flip),ambiguous_resolved=sum(d$ambiguous_pair),max_logP_difference=delta,
    significant=sum(d$reason=="OK" & d$logP>=-log10(5e-8)))
  # Greedy GWAS leads, >=1 Mb separation; +/-1 Mb flanks. Partition overlap
  # at adjacent lead midpoints, yielding disjoint loci within each analysis.
  loci<-list()
  for(ch in sort(unique(bim$chr))) {
    sig<-d[chr==ch & reason=="OK" & logP>=-log10(5e-8)]
    setorder(sig,-logP,position,SNP);leads<-sig[0]
    while(nrow(sig)) {
      lead<-sig[1];leads<-rbind(leads,lead)
      sig<-sig[abs(position-lead$position)>params$half_window]
    }
    if(!nrow(leads)) next
    setorder(leads,position)
    start<-pmax(1,leads$position-params$half_window)
    end<-pmin(max(bim[chr==ch,position]),leads$position+params$half_window)
    if(nrow(leads)>1) for(j in 1:(nrow(leads)-1)) if(end[j]>=start[j+1]) {
      split<-floor((leads$position[j]+leads$position[j+1])/2)
      end[j]<-split;start[j+1]<-split+1
    }
    loci[[length(loci)+1]]<-data.table(phenotype,method,chr=ch,start,end,
      lead_SNP=leads$SNP,lead_position=leads$position,lead_logP=leads$logP)
  }
  loc<-rbindlist(loci)
  if(nrow(loc)) {
    loc[,locus_id:=sprintf("chr%s_%s_%s",chr,start,end)]
    loc[,job_id:=paste(method,phenotype,locus_id,sep="__")]
    loc[,requested_snps:=vapply(seq_len(.N),function(k) bim[chr==loc$chr[k] & position>=loc$start[k] & position<=loc$end[k],.N],integer(1))]
    # Ensure every significant marker lies in exactly one selected interval.
    sig<-d[reason=="OK" & logP>=-log10(5e-8)]
    covered<-integer(nrow(sig))
    for(k in seq_len(nrow(loc))) covered<-covered+as.integer(sig$chr==loc$chr[k] & sig$position>=loc$start[k] & sig$position<=loc$end[k])
    stopifnot(all(covered==1))
    all_loci[[length(all_loci)+1]]<-loc
  }
}
write(rbindlist(stats),"qc/gwas_audit.tsv")
loci<-rbindlist(all_loci);stopifnot(nrow(loci)>0,!anyDuplicated(loci$job_id))
write(loci,"loci/manifest.tsv")
# Pilot selection uses only phenotype-1 BOLT GWAS significance.
pilot<-loci[method=="BOLT-LMM" & phenotype=="pheno1"][order(-lead_logP,chr,start)][1]
write(pilot,"loci/pilot.tsv")
print(rbindlist(stats));print(loci[,.(loci=.N),by=.(method,phenotype)])
cat("Pilot selected by strongest BOLT phenotype-1 peak:\n");print(pilot)
cat("Causal list and effect-window labels were not read.\n")
print(sessionInfo());sink()
writeLines("PASS", "status/audit.PASS")
