# Executed aggregate results

This directory contains **27 populated TSV files** from the executed reconstruction of the saved livestock GWAS benchmark, plus JSON provenance records. The former placeholder templates are preserved separately in [example/schemas](../example/schemas/README.md).

## Main result tables

| Analysis | Aggregate result | Per-phenotype or supporting result |
|---|---|---|
| Significant SNPs and 500 kb loci | [Method summary: 9 methods](method_summary/method_summary.tsv) | [SNP counts: 90 rows](method_summary/significant_snps_by_pheno.tsv); [locus counts: 90 rows](method_summary/significant_loci_by_pheno.tsv) |
| Power | [306 method/threshold summaries](power_summary/power_summary.tsv) | [3,060 phenotype/threshold rows](power_summary/power_by_pheno.tsv) |
| False-positive rates | [612 method/region/threshold summaries](power_summary/fpr_summary.tsv) | [6,120 phenotype/region/threshold rows](power_summary/fpr_by_pheno.tsv) |
| Spatial enrichment | [360 method/distance summaries](power_summary/spatial_enrichment_summary.tsv) | [3,600 phenotype/distance rows](power_summary/spatial_enrichment_by_pheno.tsv) |
| Central-null spillover | [9 method summaries](spillover_summary/spillover_summary.tsv) | [90 phenotype/method rows](spillover_summary/spillover_by_pheno.tsv) |
| Genomic inflation | [9 method summaries](lambda_gc_summary/lambda_gc_summary.tsv) | [90 historical lambda values and stable P-based diagnostics](lambda_gc_summary/lambda_gc_by_pheno.tsv) |
| Locus overlap | [81 method-pair summaries](overlap_summary/locus_overlap_summary.tsv) | [810 phenotype-specific pairs](overlap_summary/locus_overlap_by_pheno.tsv); [900 cross-replicate comparisons](overlap_summary/cross_replicate_stability.tsv) |
| Equation 13 downstream | [160 method/region/bin summaries](equation13_summary/equation13_summary.tsv) | [1,600 replicate-bin means](equation13_summary/equation13_by_pheno.tsv) |
| Fine-mapping | [4 pooled high-PIP summaries](finemapping_summary/high_pip_summary.tsv) | [40 phenotype/pipeline rows](finemapping_summary/high_pip_by_pheno.tsv); [null/effect counts](finemapping_summary/null_vs_effect_summary.tsv); [converged loci](finemapping_summary/converged_loci.tsv) |
| LD decay | [1,000 independently recounted 10 kb bins](ld_summary/ld_decay_10kb.tsv) | [7,998 fitted curve coordinates](ld_summary/ld_decay_fitted_mean.tsv) |
| Heritability | [20 estimates recovered from fitting logs](method_summary/heritability_from_logs.tsv) | Means/ranges are checked in the validation ledger |

For example, the pooled high-PIP table contains 1,974 calls with 1,153 in null windows for BOLT-LMM/SuSiE-adj, and 2,101 calls with 1,274 in null windows for REGENIE_LOCO/SuSiE-adj. These are independently recounted saved results, not template values.

## Provenance and limits

[Input fingerprints](input_manifest.tsv) identify all 90 source GWAS panels. [Run metadata](run_manifest.json), [cached-theory checks](equation13_summary/cached_theory_validation.tsv), and [fine-mapping verification metadata](finemapping_summary/verification_manifest.json) describe the executed calculations. The [numerical validation ledger](../audit/manuscript_validation.tsv) retains 95 agreements, two manuscript discrepancies and explicit unresolved items.

The study uses real livestock genotypes and **simulated phenotypes**. These tables recompute downstream results from saved GWAS, cached theoretical predictions, saved converged fits and pairwise LD. They are not fresh model refits. Upstream Equation 12 and the manuscript maximum-PVE provenance remain unresolved. Legitimate NA values in overlap comparisons indicate an empty locus set, as defined in the manuscript; they are not missing template data.

No animal IDs, individual phenotypes, genotypes, raw GWAS, GRMs, raw LD matrices or per-SNP PIPs are distributed here. Synthetic code-validation data are separate in [example/synthetic](../example/synthetic/README.md). See the [validation report](../repository_validation_report.md) and [data availability](../docs/data_availability.md).
