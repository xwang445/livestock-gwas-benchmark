# Nine recovered GWAS configurations

Every shell runner requires a new OUTPUT_DIR and authorized GENO_DIR prefix; it refuses an existing output directory. Defaults preserve the recovered scientific command choices. These runners were syntax-checked and selected arguments were tested with mocks; expensive fitting was NOT repeated.

| Runner | Required phenotype input | Configuration |
|---|---|---|
| run_boltlmm.sh | PHENO_FILE combined TSV | BOLT --lmm --LDscoresUseChip; ten traits |
| run_regenie.sh | PHENO_FILE combined TSV | Step1/2 --qt --bsize 1000; default LOCO |
| run_regenie_wo_loco.sh | PHENO_FILE combined TSV | Step1 --print-prs; Step2 --use-prs |
| run_fastGWA.sh | PHENO_DIR containing pheno.1.txt..10 | Full GRM, sparse cutoff 0.05, --fastGWA-mlm; final extension .fastgwa.fastGWA |
| run_slemm_gwas.sh | PHENO_FILE SLEMM CSV; SNP_INFO_FILE | Full null fit, gamma calibration, chr1–18 association; genotype prefix must have both bed and pgen companions |
| gapit.R | Five positional arguments: retained source, combined TSV, PLINK raw, matching GM, new output directory | Both FarmCPU and BLINK; Kansas is pre-local refinement and NYC post-refinement |

THREADS defaults to 16; executable paths are configurable through BOLT_BIN, REGENIE_BIN, GCTA_BIN, SLEMM_BIN, SLEMM_GAMMA_BIN and SLEMM_GWA_BIN. No covariate option was recovered in these runs; none is added. The fastGWA absolute output prefix fixes the earlier relative-path nesting bug. Per-chromosome SLEMM outputs are retained for audit.

GAPIT source is required externally and checked against the retained library hash. The library identifies itself as 4.0 while the manuscript says 3.5; logs establish execution of both models but do not establish which library bytes executed. The pin enables a reproducible retained-source run, **not a claim that the historical software version is resolved**. Required R packages must already be present; the wrapper does not silently install dependencies.
