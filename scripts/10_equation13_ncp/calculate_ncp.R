#!/usr/bin/env Rscript
# Compatibility entry point for the recovered Figure 7 downstream analysis.
# The older exploratory script is preserved in Git history; it did not implement
# Equation 12 and incorrectly called observed chi-square the theoretical NCP.
args <- commandArgs(TRUE)
if(length(args)!=3L) stop("Usage: calculate_ncp.R RESTRICTED_ROOT CACHED_THEORY_TSV NEW_OUT_DIR. Upstream Equation 12 is not implemented.")
script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script <- normalizePath(sub("^--file=","",script_arg))
core <- file.path(dirname(dirname(script)),"benchmark_core.py")
python <- Sys.getenv("PYTHON",unset="python3")
status <- system2(python,c(shQuote(core),"equation13","--root",shQuote(args[1]),"--theory",shQuote(args[2]),"--out",shQuote(args[3])))
quit(status=status)
