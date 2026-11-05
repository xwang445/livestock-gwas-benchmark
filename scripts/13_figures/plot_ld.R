#!/usr/bin/env Rscript
# Same constrained mean spline as historical ld_decay_common.R; no median fit.
args <- commandArgs(TRUE)
if(length(args)!=2) stop("Usage: plot_ld.R DERIVED_DATA FIGURE_DIR")
library(data.table);library(ggplot2);library(patchwork)
d <- fread(file.path(args[1],"ld_summary/ld_decay_10kb.tsv"))
setorder(d,bin_id)
stopifnot(nrow(d)==1000L,identical(d$bin_id,seq_len(1000L)),
          all(abs(d$bin_start-(d$bin_id-1)/100)<1e-12),
          all(abs(d$bin_end-d$bin_id/100)<1e-12),
          all(abs(d$mid_mb-(d$bin_start+d$bin_end)/2)<1e-12),
          all(d$n_pairs>0),all(is.finite(d$mean_r2)))
basis <- splines::ns(log10(d$mid_mb),df=20,intercept=TRUE)
constraint <- basis[which.min(d$mid_mb),,drop=FALSE]
system <- rbind(cbind(crossprod(basis),t(constraint)),cbind(constraint,0))
coef <- solve(system,c(crossprod(basis,d$mean_r2),d$mean_r2[1]))[seq_len(ncol(basis))]
x <- sort(unique(c(d$mid_mb,seq(min(d$mid_mb),max(d$mid_mb),length.out=5001),10^seq(log10(min(d$mid_mb)),log10(max(d$mid_mb)),length.out=2001))))
x <- unique(pmin(max(d$mid_mb),pmax(min(d$mid_mb),x)))
curve <- data.table(mid_mb=x,mean_r2=as.numeric(predict(basis,log10(x)) %*% coef))
stopifnot(abs(curve$mean_r2[1]-d$mean_r2[1])<1e-10)
curve_file <- file.path(args[1],"ld_summary/ld_decay_fitted_mean.tsv")
# Figure commands never overwrite the source dataset. Reuse its canonical
# coordinates when the fit matches within numerical solver roundoff.
if(file.exists(curve_file)) {
  saved <- fread(curve_file)
  # R versions may collapse a different number of almost-identical grid points.
  # Check the fitted function at the canonical coordinates, not row alignment.
  expected_saved <- as.numeric(predict(basis,log10(saved$mid_mb)) %*% coef)
  stopifnot(identical(names(saved),names(curve)),
            max(abs(saved$mean_r2-expected_saved))<1e-10)
  curve <- saved
} else {
  dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
  fwrite(curve,file.path(args[2],"ld_decay_fitted_mean.tsv"),sep="\t")
}
plots <- lapply(c(2,10),function(limit) ggplot(curve,aes(mid_mb,mean_r2))+geom_line(size=.7)+
  coord_cartesian(xlim=c(0,limit),ylim=c(0,.6),expand=FALSE)+theme_bw(base_size=12)+
  labs(x="Distance (Mb)",y=expression("Mean "*R^2),title=paste("LD decay to",limit,"Mb")))
stopifnot(identical(ggplot_build(plots[[1]])$data[[1]][,c("x","y")],
                    ggplot_build(plots[[2]])$data[[1]][,c("x","y")]))
p <- plots[[1]]+plots[[2]]+plot_annotation(tag_levels="a")
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
for(ext in c("pdf","png","tiff")) {
  options <- list(filename=file.path(args[2],paste0("figure1_ld.",ext)),plot=p,width=10,height=4,dpi=600)
  if(ext=="tiff") options$compression <- "lzw"
  do.call(ggsave,options)
}
message(sprintf("Shared 10-kb bins; first midpoint %.3f Mb, mean R2 %.12f; plotted coordinates identical",d$mid_mb[1],d$mean_r2[1]))
