# GWAS

The benchmark contains nine configurations across ten simulated phenotypes: 90 saved GWAS panels. To analyze these results, use `python scripts/workflow.py summarize --profile mac --out ../benchmark-rerun`. The source filenames are in [methods.yaml](../../config/methods.yaml); formats are in [input_file_formats.md](../../docs/input_file_formats.md).

## Fitting commands

These runners are for a new, separately allocated run. Each shell runner requires a genotype prefix in `GENO_DIR` and a new `OUTPUT_DIR` and refuses to overwrite existing output. No new fits were needed for the published downstream tables.

| Runner | Phenotype input | Model settings |
|---|---|---|
| `run_boltlmm.sh` | `PHENO_FILE`: combined TSV | BOLT-LMM 2.5; `--lmm --LDscoresUseChip`; infinitesimal output |
| `run_regenie.sh` | `PHENO_FILE`: combined TSV | REGENIE 4.1.2; quantitative two-step model, block size 1,000, LOCO |
| `run_regenie_wo_loco.sh` | `PHENO_FILE`: combined TSV | Step 1 `--print-prs`; step 2 `--use-prs` |
| `run_fastGWA.sh` | `PHENO_DIR`: `pheno.1.txt` through `pheno.10.txt` | GCTA 1.94.1; sparse GRM cutoff 0.05; `--fastGWA-mlm` |
| `run_slemm_gwas.sh` | `PHENO_FILE`: SLEMM CSV; `SNP_INFO_FILE` | SLEMM 0.90.1; full null fit, gamma calibration, chromosome-wise testing |
| `gapit.R` | Retained GAPIT source, combined TSV, PLINK RAW, matching marker map, new output directory | FarmCPU and BLINK; Kansas before local refinement, NYC after refinement |

SLEMM needs BED and PGEN companions for its genotype prefix. `THREADS` defaults to 16. Executables are selected through `BOLT_BIN`, `REGENIE_BIN`, `GCTA_BIN`, `SLEMM_BIN`, `SLEMM_GAMMA_BIN`, and `SLEMM_GWA_BIN`. Phenotypes are matched by individual ID; genotype-file row order is not assumed to equal simulation order. No additional covariates were used beyond the intercept.

The retained GAPIT library is externally supplied and hash-checked. Its version string is 4.0; the manuscript reports 3.5. This version distinction is recorded in [methods.yaml](../../config/methods.yaml) rather than silently changing the historical claim.

## Result definitions

Significant-SNP/locus counts use P < 5 × 10⁻⁸ and chromosome-specific 500 kb bins. The REGENIE_LOCO maximum is 20,461 significant SNPs; 21,484 is the BOLT-LMM maximum. The full-range maximum −log10(P) is 107.517 for REGENIE_LOCO phenotype 8. Both were independently checked from raw result columns in [manuscript_number_checks.tsv](../../audit/manuscript_number_checks.tsv).
