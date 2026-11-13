# Public submission figure sources, 12 November 2026

These are non-identifying aggregate source tables for the latest publication artwork. No animal IDs, SNP IDs, per-variant GWAS/PIP rows, genotypes, or individual phenotypes are included.

- `figure01_*` are retained LD distance bins and the fitted shared curve. Each SNP pair has equal weight within a distance bin; the fitted spline uses equally weighted bin means.
- `figure02_*` records the schematic window design.
- `Figure_03_*`, `Figure_05_*`, and `Figure_06_*` contain configuration/phenotype counts, 1,000 overlap cells, spillover and calibration summaries, and distance-bin ratios. Figure 6 retains the documented historical closed-endpoint spatial definition and historical central-null lambda statistic.
- `figure07_*` contains retained theoretical bin means and observed phenotype-level means; upstream variant predictions and assignments were not regenerated.
- `figure08_pr_display_vertices.tsv` is a reduced display representation. The exact AP values in `figure08_pooled_metrics.tsv` were calculated from the full SNP–phenotype ranking, not these vertices. The four current pipeline/GWAS combinations all use ERV_TRUE.
- `figure08_pve_group_summary.tsv` contains phenotype-mean PIP ≥0.9 rates and t9 intervals based on the supervisor's supplied recovered effects. PVE is a fraction of phenotype variance and the design expectation is 0.00005 (0.005%). `power`, `lo`, and `hi` are fractions; `pooled` is a different estimator. The group interval calculation retains the supplied rounded t critical value 2.262157. See `derived_data/revision_2026-11-11/` for independently rederived per-phenotype results and provenance.
- `manhattan_counts.tsv` verifies all 100 restricted-source panels; it contains only counts and display-range summaries.
- `tables/Table_S1.tsv` is the displayed four-column table. Supporting tables retain full pooled counts, all phenotype-level regional counts, genomic denominators, and central-null spillover intervals. Pooled percentages use detected locus–phenotype pairs as their denominator; zero-detection phenotypes contribute zero counts. The table generator recomputes these values from phenotype aggregates. Supporting spillover intervals use sample SD and full-precision Student t(0.975,9).

The supporting `display_method` column contains publication labels; `configuration` is the stable method identifier. Original analysis and figure sources are preserved in their dated directories. Source hashes and the restricted-input derivation workflow are documented with the corresponding revision records.
