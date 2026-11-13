# Benchmark results

The root summary folders contain the current results: all ten GWAS configurations and ten phenotypes, updated SLEMM, completed GCTA MLMA, and inclusive P ≤ 5 × 10⁻⁸. The previous root tables are preserved byte for byte in [historical_2026-11-04](historical_2026-11-04/README.md). The [November 13 summary update](../docs/revision_2026-11-13.md) records this migration.

## Current result tables

| Analysis | Summary | Phenotype-level or supporting result |
|---|---|---|
| Significant SNPs and 500 kb loci | [Ten configurations](method_summary/method_summary.tsv) | [100 SNP-count rows](method_summary/significant_snps_by_pheno.tsv); [100 locus-count rows](method_summary/significant_loci_by_pheno.tsv) |
| Genomic inflation | [Ten configuration summaries](lambda_gc_summary/lambda_gc_summary.tsv) | [100 current lambda rows, including GCTA MLMA](lambda_gc_summary/lambda_gc_by_pheno.tsv) |
| Central-null spillover | [Ten summaries](spillover_summary/spillover_summary.tsv) | [100 phenotype/configuration rows](spillover_summary/spillover_by_pheno.tsv) |
| Power | [Method/threshold summaries](power_summary/power_summary.tsv) | [Current causal-SNP rates and exact thresholds](power_summary/power_by_pheno.tsv) |
| Null-region SNP rates | [Region/threshold summaries](power_summary/fpr_summary.tsv) | [Outside-effect and central-null rates](power_summary/fpr_by_pheno.tsv) |
| Spatial enrichment | [400 configuration/distance summaries](power_summary/spatial_enrichment_summary.tsv) | [4,000 phenotype/configuration/distance rows](power_summary/spatial_enrichment_by_pheno.tsv) |
| Locus overlap | [100 method-pair summaries](overlap_summary/locus_overlap_summary.tsv) | [1,000 phenotype-specific pairs](overlap_summary/locus_overlap_by_pheno.tsv); [1,000 cross-replicate comparisons](overlap_summary/cross_replicate_stability.tsv) |
| Heritability | [30 estimates from fastGWA, current SLEMM and GCTA MLMA](method_summary/heritability_from_logs.tsv) | GCTA ratios retain rounded-log precision; no unrounded values or standard errors are inferred |
| Fine-mapping | [Four current GWAS/pipeline combinations](finemapping_summary/high_pip_summary.tsv) | [40 phenotype-specific rows](finemapping_summary/high_pip_by_pheno.tsv); [null/effect counts](finemapping_summary/null_vs_effect_summary.tsv); [converged loci](finemapping_summary/converged_loci.tsv) |
| Equation 13 / current Figure 7 | [Retained BOLT-LMM and REGENIE_LOCO bin summaries](equation13_summary/equation13_summary.tsv) | [Phenotype/bin observations](equation13_summary/equation13_by_pheno.tsv); [scope and provenance](equation13_summary/cached_theory_validation.tsv) |
| LD | [1,000 verified 10 kb bins](ld_summary/ld_decay_10kb.tsv) | [Shared fitted spline used by current artwork](ld_summary/ld_decay_fitted_mean.tsv) |

`n_significant_snps` and `n_significant_loci` are current inclusive counts. Explicit `strict_*` fields retain comparisons at the older strict cutoff. The root files retain the established method key `REGENIE_woLOCO`; dated association sources use `REGENIE_WO_LOCO`. GCTA's stable key is `GCTA_MLMA`.

`lambda_gc` retains the central-null statistic estimator used in Figure 6, with the current inputs for each method. Its `statistic_definition` column identifies the estimator; `lambda_gc_p_based` provides the stable P-based alternative. The word "historical" in an estimator name describes that calculation, rather than an old SLEMM or incomplete GCTA fit. Method summaries use sample SD across ten phenotypes.

Power/null-rate tables retain the exact audited threshold grid and its hexadecimal representation. Spatial enrichment retains the established closed-endpoint 50 kb bins. Empty locus-set comparisons remain undefined where specified; overlap tables distinguish the number of paired phenotypes from the number with defined Jaccard values.

Fine-mapping was performed for BOLT-LMM and REGENIE_LOCO inputs only. Both current pipelines estimate residual variance; RSS high-PIP counts are 2,021 and 2,072. No GCTA fine-mapping results are claimed. Current Equation 13 exports cover the two inputs shown in Figure 7 using retained genome-wide bins. The original four-method/central-null theory audit is archived; the upstream theoretical generator remains unverified.

## Dated canonical analyses and submission sources

[revision_2026-11-05/association](revision_2026-11-05/association) supplies the detailed current counts, lambda, overlap, spatial and sensitivity sources. [revision_2026-11-05/susie](revision_2026-11-05/susie) retains both final TRUE and paired historical FALSE fine-mapping analyses. [revision_2026-11-11](revision_2026-11-11/README.md) contains recovered-effect lead/PVE summaries and separately labeled GWAS-coefficient PVE. [submission_2026-11-12](submission_2026-11-12/README.md) contains figure aggregates and the current pooled four-column Table S1. Their dated definitions remain intact.

## Rebuild and verify

```bash
python scripts/workflow.py refresh-results
python scripts/workflow.py check
```

The first command exports the root views from committed aggregate sources and refreshes result/artifact fingerprints. The check compares root exports with those sources, requires complete ten-configuration GWAS coverage, verifies the historical snapshot checksums and runs the numerical tests. [Current input fingerprints](input_manifest.tsv) and [run metadata](run_manifest.json) describe all 100 GWAS panels.

The portable [support derivation](../scripts/revision_2026-11-05/association/derive_public_summary_support.py) produces aggregate QC, cross-replicate overlap and relative source fingerprints from authorized private records. Detailed loci, SNP IDs and individual-level data remain excluded.

## Historical audit

The [historical snapshot manifest](historical_2026-11-04/snapshot_manifest.json) records the original bytes and checksums of 29 files, including 27 TSVs and two metadata files. Its nine-configuration strict association results, original SLEMM, RSS FALSE results and cached four-method theory comparison support the historical manuscript ledger. Use explicit `--historical` commands to reproduce that audit. Historical files never supply current primary summary rows.

Individual-level inputs, raw GWAS, GRMs, raw LD, causal variant identifiers and per-SNP PIPs remain restricted. [Synthetic examples](../example/synthetic/README.md) are separate from study results.
