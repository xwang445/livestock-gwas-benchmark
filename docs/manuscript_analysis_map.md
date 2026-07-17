# Manuscript analysis map

| manuscript section or figure | analysis purpose | source script | expected input | expected output | public data included? | proprietary input required? |
| --- | --- | --- | --- | --- | --- | --- |
| Figure 6A | Spillover-locus rate in central-null 500 kb windows | scripts/07_spillover/figure6A_spillover_locus_500kb.R | GWAS results, window definitions | spillover-rate table and figure | No | Yes |
| Figure 6B | Central-null lambda GC | scripts/08_lambda_gc/lambda_gc_range.R | GWAS results | lambda GC summary | No | Yes |
| Figure 6C / 6D | 50 kb null-window false-positive rates | scripts/07_spillover/fpr_50kb_window.R | GWAS results, null windows | FPR summary and figure | No | Yes |
| Equation 13 / Figure 8 | Theoretical NCP and comparison with observed chi-square | scripts/10_equation13_ncp/calculate_ncp.R | effect sizes, LD scores, GWAS statistics | NCP diagnostic summary | No | Yes |
| Significant SNP summaries | Count significant SNPs and loci | scripts/06_significant_loci/summarize_significant_snps_and_loci.R | harmonized GWAS outputs | counts by method and phenotype | No | Yes |
| Power and realized PVE | Max-SNP PVE and power calculations | scripts/05_effect_size_power/calculate_max_snp_pve.sh and scripts/05_effect_size_power/power_analysis.R | causal SNP list, phenotype data, GWAS outputs | realized PVE and power summary | No | Yes |
