# Genotype preparation

`prepare_panel.sh` recovers INFO ≥ 0.8, MAF ≥ 0.01, HWE P ≥ 1e-9, 20 kb spacing, and the union with QC-passing chip SNPs. Supply SOURCE_PGEN_DIR, CHIP_IDS, OUTPUT_DIR (new), optional PLINK2_BIN and THREADS. It retains intermediates and checks 135,358 SNPs / 31,301 animals. It has not been executed against proprietary source pgens in this reconstruction.

The original chip-list bytes and historical genotype identity/order are unresolved; satisfying counts is insufficient proof of identical inputs. GWAS uses the first 20,000 individuals through phenotype selection. No unrecorded missingness threshold or imputation algorithm is added. See docs/input_file_formats.md.
