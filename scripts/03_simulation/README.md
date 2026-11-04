# Historical simulation

```bash
Rscript scripts/03_simulation/select_causal_snps.R /path/to/QC_chip.snplist /path/to/chrom_sizes /path/to/new/architecture
bash scripts/03_simulation/simulate_polygenic_phenotypes.sh /path/to/historical/genotype /path/to/new/architecture/snp.info.chr1_18.10k.csv /path/to/new/simulation
Rscript scripts/03_simulation/prepare_benchmark_phenotypes.R /path/to/new/simulation/pheno.sim.csv /path/to/new/phenotypes
```

The commands require verified original sample/SNP identity and order; do not point them at a known replaced historical path. No original simulation was rerun here.

Causal selection uses R Mersenne-Twister/Inversion/Rejection, seed 1, chr1–18 and position modulo 6 Mb < 2 Mb. MPH v0.55.1 uses an equal-weight standardized-genotype GRM, listed with variance weight 0.5. Its seed argument 0 seeds mt19937 at 54321. The historical command generates 1,000 replicates; Gaussian genetic draws occupy the first 1,000 random-vector columns and residuals the next 1,000. Reducing the generation count changes the residual stream even for the first ten retained traits.

Versioned MPH source (`13ffe633e0170d35e01e6164fc2bb4e01f35e2f2`, src/option.cpp) computes genetic values from Cholesky(0.5 G + 0.00005 I), adds residual variance 0.5, and does not rescale realized phenotypes. It does not output explicitly drawn per-SNP effects. These are implementation facts, not a newly executed simulation validation.

The phenotype converter preserves the first 20,000 IID rows and columns 1–10, writes FID=0 as in the historical converter, and creates the combined, fastGWA and SLEMM formats. Private IDs/phenotypes must remain excluded from Git. Software details: https://github.com/jiang18/mph/tree/13ffe633e0170d35e01e6164fc2bb4e01f35e2f2
