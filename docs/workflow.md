# Workflow

1. Prepare or thin the genotype panel.
   - Script: scripts/01_genotype_preparation/README.md
   - Inputs: genotype resources, SNP list, phenotype sample IDs
   - Outputs: thinned genotype panel and sample subset
   - Requirements: proprietary data required

2. Calculate LD decay.
   - Script: scripts/02_ld_analysis/README.md
   - Inputs: genotype panel
   - Outputs: LD-decay summaries
   - Requirements: proprietary data required

3. Define effect and null windows.
   - Script: scripts/03_simulation/README.md
   - Inputs: effect-window definitions
   - Outputs: effect and null window tables
   - Requirements: no proprietary data if provided as a public example

4. Select causal variants.
   - Script: scripts/03_simulation/README.md
   - Inputs: analysis panel and window definitions
   - Outputs: causal SNP list
   - Requirements: proprietary data required

5. Simulate phenotypes.
   - Script: scripts/03_simulation/simulate_polygenic_phenotypes.sh
   - Inputs: causal variant list, genotype panel
   - Outputs: phenotype replicates
   - Requirements: proprietary data required

6. Prepare the 20,000-animal subset.
   - Script: scripts/03_simulation/README.md
   - Inputs: phenotype file and genotype sample list
   - Outputs: subset genotype/phenotype files
   - Requirements: proprietary data required

7. Run all GWAS methods.
   - Scripts: scripts/04_gwas/*
   - Inputs: prepared genotype and phenotype data
   - Outputs: GWAS summary statistics
   - Requirements: proprietary data required

8. Harmonize GWAS outputs.
   - Script: scripts/06_significant_loci/summarize_significant_snps_and_loci.R
   - Inputs: method-specific GWAS outputs
   - Outputs: harmonized summary tables
   - Requirements: proprietary data required

9. Summarize significant SNPs and 500 kb loci.
   - Script: scripts/06_significant_loci/summarize_significant_snps_and_loci.R
   - Inputs: harmonized GWAS results
   - Outputs: significant SNP and locus counts
   - Requirements: proprietary data required

10. Calculate method overlap.
    - Script: scripts/09_locus_overlap/README.md
    - Inputs: significant loci by method
    - Outputs: Jaccard summaries
    - Requirements: proprietary data required

11. Calculate spillover-locus rates.
    - Script: scripts/07_spillover/figure6A_spillover_locus_500kb.R
    - Inputs: GWAS results and window definitions
    - Outputs: spillover summaries
    - Requirements: proprietary data required

12. Calculate central-null lambda GC.
    - Script: scripts/08_lambda_gc/lambda_gc_range.R
    - Inputs: GWAS results and null windows
    - Outputs: lambda GC summaries
    - Requirements: proprietary data required

13. Calculate 50 kb observed-to-expected ratios.
    - Script: scripts/07_spillover/fpr_50kb_window.R
    - Inputs: GWAS results and null windows
    - Outputs: FPR summaries
    - Requirements: proprietary data required

14. Calculate cross-replicate locus stability.
    - Script: scripts/09_locus_overlap/README.md
    - Inputs: significant loci across 10 replicates
    - Outputs: overlap and stability summaries
    - Requirements: proprietary data required

15. Calculate realized maximum PVE and power.
    - Scripts: scripts/05_effect_size_power/calculate_max_snp_pve.sh and scripts/05_effect_size_power/power_analysis.R
    - Inputs: causal variant lists and GWAS outputs
    - Outputs: realized PVE and power summaries
    - Requirements: proprietary data required

16. Calculate Equation 12 LOCO coefficients.
    - Script: scripts/10_equation13_ncp/calculate_ncp.R
    - Inputs: GWAS results and causal LD summaries
    - Outputs: LOCO coefficient summaries
    - Requirements: proprietary data required

17. Calculate causal LD scores.
    - Script: scripts/10_equation13_ncp/calculate_ncp.R
    - Inputs: genotype data and causal SNP list
    - Outputs: LD-score summaries
    - Requirements: proprietary data required

18. Calculate Equation 13 NCP.
    - Script: scripts/10_equation13_ncp/calculate_ncp.R
    - Inputs: causal LD scores and effect sizes
    - Outputs: non-centrality summaries
    - Requirements: proprietary data required

19. Compare predicted and observed chi-square.
    - Script: scripts/10_equation13_ncp/calculate_ncp.R
    - Inputs: theoretical NCP and GWAS results
    - Outputs: chi-square comparison tables
    - Requirements: proprietary data required

20. Generate manuscript figures.
    - Scripts: scripts/11_figures/README.md
    - Inputs: derived summary tables
    - Outputs: figure files
    - Requirements: proprietary data required for some panels
