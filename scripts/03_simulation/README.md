# MPH simulation

The historical simulation used MPH v0.55.1, a fixed set of 10,000 causal chip SNPs, and target heritability 0.5. MPH drew genetic values from N(0, 0.5G + 0.00005I) and independent residuals from N(0, 0.5I), then added them. It did not explicitly draw or save SNP effects or rescale the resulting phenotypes.

The command generated 1,000 replicates with seed argument 0, which initializes mt19937 with 54321. The benchmark uses the first ten phenotypes and first 20,000 saved simulation rows. Keep the generation count at 1,000 when reproducing the random stream.

| Script | Inputs | Output |
|---|---|---|
| `select_causal_snps.R` | Ordered QC chip list, chromosome sizes, new output directory | Fixed causal list and effect windows |
| `simulate_polygenic_phenotypes.sh` | Verified retained chip prefix, causal list, new output directory | MPH GRM and 1,000 simulated phenotypes |
| `prepare_benchmark_phenotypes.R` | Saved MPH simulation CSV, new output directory | First-ten/first-20,000 combined, fastGWA and SLEMM files |

Causal selection uses R seed 1 and samples 10,000 of 10,461 eligible QC chip SNPs within the repeating 2 Mb effect windows. The same positions are used across phenotypes.

The retained `FineMapping/spg` genotype order matches the saved simulation and checked GRM entries. The files now at the old simulation source path do not. Always align analysis phenotypes to genotype rows by IID. The full [provenance record](../../docs/scientific_provenance.md) and [read-only checker](../provenance/verify_sample_order.py) document this chain. No new simulation was run during this update.
