# GWAS Benchmarking Scripts

This directory contains scripts for running various Genome-Wide Association Study (GWAS) tools on livestock data.

## Scripts Overview

### 1. **run_boltlmm.sh**
Runs BOLT-LMM (Bayesian and Likelihood Optimization Techniques for Linear Mixed Models) for GWAS analysis.

**Key Configuration Variables:**
- `GENO_DIR`: Path to genotype data files
- `PHENO_FILE`: Path to phenotype file
- `OUTPUT_DIR`: Output directory for results
- `THREADS`: Number of threads to use

### 2. **run_fastGWA.sh**
Runs fastGWA (GCTA-fastGWA) for fast genome-wide association analysis using sparse GRM.

**Key Configuration Variables:**
- `GENO_DIR`: Path to genotype data files (PLINK binary format)
- `PHENO_DIR`: Path to phenotype data files
- `OUTPUT_DIR`: Output directory for results
- `PREFIX`: Prefix for output files
- `THREADS`: Number of threads to use

**Steps:**
1. Creates genetic relationship matrix (GRM)
2. Creates sparse GRM
3. Runs fastGWA-MLM for association analysis

### 3. **run_regenie.sh**
Runs REGENIE (Regression-based Efficient Genome-wide Association Tests with Efficient Interaction) with LOCO (Leave-One-Chromosome-Out) validation.

**Key Configuration Variables:**
- `GENO_DIR`: Path to genotype data files (PLINK binary format)
- `PHENO_FILE`: Path to phenotype file
- `OUTPUT_DIR`: Output directory for results
- `THREADS`: Number of threads to use
- `BSIZE`: Batch size for model fitting

**Steps:**
1. Step 1: Fits null LMM model
2. Step 2: Tests association with LOCO cross-validation

### 4. **run_regenie_wo_loco.sh**
Runs REGENIE without LOCO (Leave-One-Chromosome-Out) validation, using polygenic risk scores (PRS) instead.

**Key Configuration Variables:**
- `GENO_DIR`: Path to genotype data files (PLINK binary format)
- `PHENO_FILE`: Path to phenotype file
- `OUTPUT_DIR`: Output directory for results
- `THREADS`: Number of threads to use
- `BSIZE`: Batch size for model fitting

### 5. **run_slemm_gwas.sh**
Runs SLEMM (Spatial and Linear/nonlinear Effects Mixed Model) for GWAS analysis with multiple steps.

**Key Configuration Variables:**
- `GENO_DIR`: Path to genotype data files (PLINK binary format)
- `PHENO_FILE`: Path to phenotype file (CSV format)
- `SNP_INFO_FILE`: Path to SNP information file
- `OUTPUT_DIR`: Output directory for results
- `THREADS`: Number of threads to use

**Steps:**
1. Fits null LMM model
2. Computes GRAMMAR-Gamma
3. Runs per-chromosome GWAS
4. Merges results across chromosomes

### 6. **gapit.R**
Runs GAPIT (Genome Association and Prediction Integrated Tool) with BLINK (Bayesian-information and Linkage-disequilibrium Iteratively Nested Keyway) method.

**Key Configuration Variables:**
- `R_LIB_PATH`: Path to R library installation
- `WORK_DIR`: Working directory for output
- `GAIT_FUNCTIONS_FILE`: Path to GAPIT functions source file
- `PHENO_FILE`: Path to phenotype file
- `GENO_FILE`: Path to genotype file (PLINK recode A format)
- `GM_FILE`: Path to genetic map file

## Usage Instructions

### Prerequisites
- BOLT-LMM, GCTA, REGENIE, SLEMM, or R/GAPIT installed
- Genotype data in appropriate format (PLINK binary or other formats as specified)
- Phenotype data in required format

### Basic Setup
1. Open the desired script
2. Modify the path variables at the top of the script:
   ```bash
   GENO_DIR=/your/path/to/genotype
   PHENO_FILE=/your/path/to/phenotype
   OUTPUT_DIR=./results
   ```
3. Make the script executable:
   ```bash
   chmod +x run_boltlmm.sh
   ```
4. Run the script:
   ```bash
   ./run_boltlmm.sh
   ```

### Configuration Examples
Replace placeholder paths with your actual data locations:
- `/path/to/genotype/` → Your PLINK binary genotype files (without .bed/.bim/.fam extension)
- `/path/to/phenotype/` → Your phenotype data file
- `./results` → Your desired output directory

## Output Files

Each tool produces different output formats:
- **BOLT-LMM**: `.stats` files with association statistics
- **fastGWA**: `.fastgwa` files with association results
- **REGENIE**: `.regenie` files for each phenotype and step
- **SLEMM**: `.gwas.txt` files for each chromosome
- **GAPIT**: GAPIT output directory with association results

## Notes

- All scripts are configured to run analysis for phenotypes 1-10
- Modify the loop `for i in {1..10}` in each script to change the number of phenotypes
- Update `THREADS` variable based on your system's CPU availability
- Ensure sufficient disk space for intermediate files, especially for sparse GRM calculations

## References

- [BOLT-LMM](https://alkesgroup.broadinstitute.org/BOLT-LMM/)
- [GCTA](http://cnsgenomics.com/software/gcta/)
- [REGENIE](https://rgcgithub.github.io/regenie/)
- [SLEMM](https://github.com/niekwit/slemm)
- [GAPIT](http://zzlab.net/GAPIT/)
