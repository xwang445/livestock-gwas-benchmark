# Genotype preparation

This upstream stage constructs the 135,358-SNP, 31,301-animal panel. It is separate from the saved-result workflow in [docs/workflow.md](../../docs/workflow.md).

`prepare_panel.sh` requires `SOURCE_PGEN_DIR`, `CHIP_IDS`, and a new `OUTPUT_DIR`. `PLINK2_BIN` and `THREADS` select the executable and resources. The script preserves INFO ≥ 0.8, MAF ≥ 0.01, HWE P ≥ 10⁻⁹, 20 kb thinning, and the union with QC-passing chip SNPs. It retains intermediate files and checks panel dimensions.

The simulation used the retained chip dataset, not the imputed GWAS panel. The first 20,000 saved simulation IDs define the analysis sample; genotype rows are aligned by IID. See [sample provenance](../../docs/scientific_provenance.md). Panel construction was not rerun in this update. Array/imputation metadata and the exact historical GAPIT version remain author metadata items.
