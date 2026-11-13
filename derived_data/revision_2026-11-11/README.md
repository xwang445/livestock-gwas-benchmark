# November 11 recovered-effect and PVE aggregates

These tables were regenerated with the portable
[November 11 calculators](../../scripts/revision_2026-11-11/). See the
[scientific definitions and source limitations](../../docs/recovered_effects_2026-11-11.md).

| Directory / file | Contents |
|---|---|
| `recovered_effects/recovered_lead_pve_summary.tsv` | Four current methods' lead counts, causal shares and recovered-PVE ratios/percentiles |
| `recovered_effects/recovered_effect_summary.tsv` | Supplied recovered-PVE summaries per phenotype and across 100,000 causal pairs |
| `recovered_effects/pve_normalization_check.tsv` | Equality of supplied alpha²/PVE and saved 31,301-animal phenotype variance |
| `finemapping/pve_groups_all_runs.tsv` | PVE-group rates and replicate CIs for adjusted TRUE, RSS TRUE and historical RSS FALSE |
| `finemapping/pve_groups_by_phenotype.tsv` | Group counts/calls/rates for every phenotype and run |
| `finemapping/independently_recounted_pooled_metrics.tsv` | Independently recounted evaluated and high-PIP class counts |
| `finemapping/claim_evidence.tsv` | Numerical final-fit rate ranges and supplied recovered-PVE limits |
| `gcta_estimated_pve/estimated_pve_summary.tsv` | GCTA association-coefficient PVE summaries by group/phenotype and pooled |
| `gcta_estimated_pve/phenotype_variances.tsv` | Analysis phenotype variances and aggregate GWAS/lead counts |

Manifests record source roles and hashes without local paths. All outputs omit
SNP IDs, genotype/phenotype records and detailed PIPs. Current Figure 8 display
data and final-fit PVE-group summaries are indexed separately in
[`submission_2026-11-12`](../submission_2026-11-12/).
