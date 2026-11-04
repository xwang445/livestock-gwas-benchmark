# Restricted input contracts

All file names below are relative to an authorized benchmark root. No individual-level input is distributed.

| Input | Required content |
|---|---|
| Final PLINK prefix | .bed/.bim/.fam; 135,358 SNPs; historical full-resource sample count 31,301; original sample labels/order must match simulation |
| QC chip list | One SNP ID per line, original order; 31,062 QC-passing chip SNPs, 10,461 eligible in effect windows |
| Chromosome sizes | Two columns without header: chromosome and length; chromosomes 1–18 |
| `eligible_regions.csv` | chr,start,end; historical 2 Mb effect intervals including terminal truncations |
| `snp.info.chr1_18.10k.csv` | SNP column; 10,000 unique chr_position-based identifiers, shared across all phenotypes; restricted |
| `pheno.sim.csv` | MPH iid, columns 1..1000, followed by gv/egv; first 20,000 rows and first ten phenotype columns selected |
| `simu_pheno/pheno_combined.txt` | FID,IID,pheno1..pheno10; historical FID=0; tab-separated |
| `simu_pheno/fastGWA_data/pheno.N.txt` | Header FID,IID,phenoN; historical FID=0 |
| SLEMM phenotype | IID,pheno1..pheno10 CSV, derived by removing FID; SNP annotation required separately |
| GAPIT genotype | PLINK recode-A table with first six metadata columns; accompanying GM SNP names exactly match the remaining column names and order (including any allele suffix) |
| Equation 13 cache | SNP/chr/pos, causal_ldscore,c_loco,ncp_framework,expected_chi2,is_central_null; 135,358 unique coordinates |

`benchmark_core.py input_path` documents the exact 90 GWAS filenames. GAPIT files may be in FarmCPU/BLINK subdirectories or GAPIT; conflicting copies cause failure. The harmonizer requires SNP ID, chromosome, position and P/LOG10P. Mixed-model outputs also provide BETA/SE; REGENIE and SLEMM provide native CHISQ. The common schema is phenotype,method,snp,chr,bp,effect_allele,other_allele,beta,se,p,chisq, with extra columns recording statistic definitions and source orientation. Missing GAPIT effects/SE/alleles remain NA; SLEMM REF/ALT are retained separately pending effect-orientation confirmation. Every per-variant export remains restricted.

No imputation or empirical field-trait collection is implemented here. The imputed pgens are upstream restricted inputs. Array name and imputation software/version remain review items; neither absent missingness thresholds nor unrecorded covariates are invented.
