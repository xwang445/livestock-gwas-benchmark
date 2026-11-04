# Recovered execution workflow

The selected manuscript is the scientific specification. The shared `scripts/benchmark_core.py` implements stages 05–11 in one validated pass to avoid nine separate parsers and repeated loading of all 90 files. Older numbered directories remain for compatibility; do not assume their numeric prefixes describe the current order.

| Stage | Input → implementation → output | Next stage |
|---|---|---|
| 01 Genotype preparation | Restricted chromosome pgens + original chip list → `01_genotype_preparation/prepare_panel.sh` → 31,301-animal / 135,358-SNP panel | LD, simulation, GWAS |
| 02 LD | Full resource PLINK panel → `02_ld_analysis/run_plink_ld.sh`, `verify_ld_means.py` → restricted LD / public 10 kb means | Figure 1 |
| 03 Simulation | Ordered QC chip SNPs + chromosome lengths → `select_causal_snps.R`; historical ordered genotype → `simulate_polygenic_phenotypes.sh`; `prepare_benchmark_phenotypes.R` → restricted ten-trait subset | GWAS |
| 04 GWAS | Authorized genotypes + ten traits → five shell runners and `04_gwas/gapit.R` → 90 method-specific outputs | Harmonization; fine-mapping |
| 05 Harmonization | Saved outputs → `benchmark_core.py harmonize` → restricted common schema | Downstream statistics |
| 06 Power/FPR | Shared truth, windows, harmonized P → `benchmark_core.py summarize` → `power_summary/` | Tables/figures |
| 07 Significant loci | Strict P and chromosome-position bins → same command → `method_summary/` | Figure 3, overlap |
| 08 Spillover | Inclusive P and central-null window IDs → same command → `spillover_summary/` | Figure 6A |
| 09 Lambda | Historical statistic and central-null mask → same command → `lambda_gc_summary/` | Figure 6B |
| 10 Overlap | Strict significant-locus sets → same command → `overlap_summary/` | Figure 5 |
| 11 Equation 13 | Cached theory + four methods × ten traits → `benchmark_core.py equation13` → `equation13_summary/` | Figure 7; upstream Equation 12 unresolved |
| 12 Fine-mapping | Restricted genotypes/GWAS/LOCO inputs → installed `12_finemapping` pipelines → restricted fits; `verify_saved_results.py` → aggregate JSON | Figure 8 |
| 13 Figures | Verified aggregates → `13_figures/plot_aggregates.py`, `plot_ld.R`; restricted GWAS → `plot_manhattan.py` | Reviewable PDF/PNG |

The first four stages and new fine-mapping fits are prepared but were not re-executed. Current verification starts from saved GWAS, existing pairwise LD and converged fitted PIPs. Old exploratory NCP code was replaced by an explicit downstream compatibility wrapper; its prior content remains in Git history and historical source directories.
