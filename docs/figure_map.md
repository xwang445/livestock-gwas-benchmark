# Figure reproduction map

All 19 figures were generated as PDF and PNG and visually inspected. Artwork is a redraw: colors, panel arrangements and typography need not be pixel-identical. No new GWAS/SuSiE model fitting is implied.

| Figure | Analysis input | Public derived table | Generator | Final base path under figures/reconstructed |
|---|---|---|---|---|
| 1 | Saved pairwise LD, 31,301 animals | ld_summary/ld_decay_10kb.tsv, ld_decay_fitted_mean.tsv | plot_ld.R | figure1_ld |
| 2 | Recovered 2/4 Mb architecture | config/analysis_parameters.yaml defines constants; geometry is shown schematically | plot_aggregates.py | figure2_windows |
| 3 | 90 harmonized GWAS outputs | method_summary/significant_snps_by_pheno.tsv | plot_aggregates.py | figure3_counts |
| 4 | 90 restricted GWAS outputs; strict significance | manhattan_manifest.tsv contains counts/maxima only | plot_manhattan.py | manhattan/figure4_manhattan |
| 5 | Strict locus sets | overlap_summary/locus_overlap_by_pheno.tsv | plot_aggregates.py | figure5_jaccard |
| 6 | Central-null and spatial summaries | spillover_summary, lambda_gc_summary, power_summary/spatial_enrichment_summary.tsv | plot_aggregates.py | figure6_spillover_inflation |
| 7 | Cached prediction and observed chi-square | equation13_summary/equation13_summary.tsv | plot_aggregates.py | figure7_equation13 |
| 8 | Independently recounted saved PIPs | finemapping_summary/high_pip_summary.tsv | plot_aggregates.py | figure8_finemapping |
| S1–S10 | All SNPs in each phenotype's nine GWAS files | manhattan_manifest.tsv records row count and maximum | plot_manhattan.py | manhattan/figureS01_manhattan through figureS10_manhattan |
| Extra central-null | Re-binned cached theory and observed values | equation13_summary/equation13_summary.tsv | plot_aggregates.py | supplement_equation13_central |

Figure 4 uses the retained source's 30+1 display limit; the caption's approximately 32 is approximate. All ten full-range supplements include the complete observed statistic range. Maximum −log10(P) is 107.517, not approximately 102. P=0 would be plotted at a labeled computational floor (1e-300), but the real panels' zero-P counts are recorded separately. This display operation does not alter statistical results.

Figure 6 means/SDs use historical finite-only lambda definitions, not the stable P-based sensitivity estimate. Figure 7's interval unit is the phenotype replicate, not individual SNPs. Figure 8 presents the primary causal-precision/MAF comparison and additionally visualizes the reported null-window proportions. Its pipeline comparison changes LD, effective sample size and residual-variance estimation together.
