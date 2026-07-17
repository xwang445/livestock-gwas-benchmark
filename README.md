# livestock-gwas-polygenic-benchmark

Repository for the manuscript: Widely used GWAS methods can produce misleading localization under diffuse polygenic architecture in livestock.

GitHub repository: https://github.com/xwang445/livestock-gwas-polygenic-benchmark

## Project title

Widely used GWAS methods can produce misleading localization under diffuse polygenic architecture in livestock

## Scientific question

This repository contains analysis code for benchmarking GWAS localization under diffuse polygenic architecture in Duroc pigs. The focus is on how commonly used GWAS methods can overstate the precision of causal localization when the signal is distributed across many loci.

## Simulation design summary

The analysis workflow was designed around a large-scale Duroc pig resource with:

- 31,301 genotyped animals as the genomic resource;
- 20,000 phenotyped animals used for GWAS;
- a 135,358-SNP analysis panel;
- 10,000 causal SNPs;
- heritability $h^2 = 0.5$;
- 10 phenotype replicates;
- alternating 2 Mb effect windows and 4 Mb null windows;
- central-null regions;
- 500 kb locus windows; and
- a genome-wide significance threshold of $P \le 5 \times 10^{-8}$.

## Methods evaluated

- BOLT-LMM
- REGENIE with LOCO
- REGENIE without LOCO
- fastGWA
- FarmCPU (Kansas)
- FarmCPU (NYC)
- BLINK (Kansas)
- BLINK (NYC)
- SLEMM

## Repository structure

- scripts/: analysis scripts grouped by workflow stage
- derived_data/: non-identifying derived summary templates
- docs/: workflow and manuscript mapping documentation
- environment/: software and package notes
- config/: example configuration file for local paths
- example/: minimal usage examples

## Required software

- PLINK / PLINK2
- REGENIE v4.1.2
- BOLT-LMM v2.5
- GCTA v1.94.1
- GAPIT v3.5
- SLEMM v0.90.1
- R with data.table, ggplot2, scales, dplyr, and stringr
- Python with numpy, pandas, scipy, and matplotlib

## Configuring paths

Copy the example configuration in config/example_paths.yaml and edit the placeholders to match your local execution environment. The scripts are written to accept command-line arguments or environment variables in place of hard-coded paths.

## Reproduction order

1. Prepare or thin the genotype panel.
2. Calculate LD decay.
3. Define effect and null windows.
4. Simulate phenotypes.
5. Prepare the 20,000-animal subset.
6. Run all GWAS methods.
7. Harmonize GWAS outputs and SNP identifiers.
8. Summarize significant SNPs and 500 kb loci.
9. Calculate spillover-locus rates and genomic inflation.
10. Calculate realized PVE and power.
11. Calculate Equation 12 / Equation 13 non-centrality summaries.
12. Generate manuscript figures.

## Figure-to-script map

- Figure 6A: scripts/07_spillover/figure6A_spillover_locus_500kb.R
- Figure 6B: scripts/08_lambda_gc/lambda_gc_range.R
- Figure 6C / 6D: scripts/07_spillover/fpr_50kb_window.R
- Equation 13 / Figure 8: scripts/10_equation13_ncp/calculate_ncp.R
- Significant SNP summaries: scripts/06_significant_loci/summarize_significant_snps_and_loci.R
- Power and realized PVE: scripts/05_effect_size_power/power_analysis.R and scripts/05_effect_size_power/calculate_max_snp_pve.sh

## Data availability

The individual-level Duroc pig genotype and phenotype data used in this study are proprietary to Smithfield Premium Genetics and are not publicly available. Access to these data may be available from the corresponding author and Smithfield Premium Genetics upon reasonable request and subject to a data-sharing agreement and any applicable institutional approvals.

This repository contains the public-safe analysis code, documentation, and non-identifying derived summary templates needed to inspect the workflow and adapt the analysis for similar projects. Complete reruns require access to the underlying Duroc data and authorized computational resources.

## Citation

Please cite the manuscript when using this repository. A citation file is included in CITATION.cff.

## Contact

Please contact the corresponding author of the manuscript for data access requests and collaboration inquiries. For repository-related questions, please open an issue or contact the repository maintainer through the GitHub project page.

## Notes on data sharing

This repository does not include individual-level genotype or phenotype data. Scripts can be inspected and reused, but complete reruns require authorized access to the Duroc data.
