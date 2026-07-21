# livestock-gwas-benchmark

Repository for the manuscript: Widely used GWAS methods can produce misleading localization under diffuse polygenic architecture in livestock.

GitHub repository: https://github.com/xwang445/livestock-gwas-benchmark

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

```text
livestock-gwas-benchmark/
├── config/                 # Example local-path configuration
├── derived_data/           # Non-identifying derived summary templates
├── docs/                   # Workflow, input-format, and manuscript mapping documentation
├── environment/            # Software versions and package requirements
├── example/                # Minimal usage examples
└── scripts/                # Analysis scripts grouped by workflow stage
    ├── 01_genotype_preparation/
    ├── 02_ld_analysis/
    ├── 03_simulation/       # MPH-based polygenic phenotype simulation
    ├── 04_gwas/             # Runners for BOLT-LMM, fastGWA, REGENIE, SLEMM, and GAPIT/BLINK
    ├── 05_effect_size_power/
    ├── 06_significant_loci/
    ├── 07_spillover/
    ├── 08_lambda_gc/
    ├── 09_locus_overlap/
    ├── 10_equation13_ncp/
    └── 11_figures/
```

## Running the GWAS scripts

The method-specific runners are in [`scripts/04_gwas/`](scripts/04_gwas/). Each script contains placeholder paths that must be changed for the local genotype, phenotype, software, and output locations before execution.

| Script | Method | Main workflow |
| --- | --- | --- |
| [`run_boltlmm.sh`](scripts/04_gwas/run_boltlmm.sh) | BOLT-LMM | Runs the linear mixed-model association test for phenotypes 1–10 |
| [`run_fastGWA.sh`](scripts/04_gwas/run_fastGWA.sh) | fastGWA | Builds a GRM and sparse GRM, then runs fastGWA-MLM |
| [`run_regenie.sh`](scripts/04_gwas/run_regenie.sh) | REGENIE with LOCO | Fits the step-1 model and performs step-2 association testing |
| [`run_regenie_wo_loco.sh`](scripts/04_gwas/run_regenie_wo_loco.sh) | REGENIE without LOCO | Fits the model and performs association testing using PRS predictions |
| [`run_slemm_gwas.sh`](scripts/04_gwas/run_slemm_gwas.sh) | SLEMM | Fits the null LMM, computes GRAMMAR-Gamma, runs chromosomes 1–18, and merges results |
| [`gapit.R`](scripts/04_gwas/gapit.R) | GAPIT/BLINK | Loads phenotype, genotype, and map files; validates marker alignment; and runs BLINK |

See the [GWAS scripts README](scripts/04_gwas/README.md) for required variables, input expectations, output types, and example commands. The supplied runners default to 10 phenotype replicates; adjust their phenotype loops if a different number of replicates is required.

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
4. Simulate phenotypes with [`simulate_polygenic_phenotypes.sh`](scripts/03_simulation/simulate_polygenic_phenotypes.sh).
5. Prepare the 20,000-animal subset.
6. Run the method-specific workflows in [`scripts/04_gwas/`](scripts/04_gwas/).
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
