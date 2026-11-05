# SuSiE-RSS residual-variance comparison

These sources reproduce the November 5 rerun. The fitter retains the original loci, GWAS statistics, in-sample LD, sample size 20,000, allele matching, L=10, coverage=0.95, minimum absolute correlation=0.5 and maximum 1,000 iterations. The changed fitting option is `estimate_residual_variance=TRUE`. The old ERV_FALSE outputs remain separate under the original `SuSiE_finemapping` directory. The unchanged SuSiE-adj/LOCO fits provide the other arm of the original fine-mapping comparison.

Use the existing R environment with susieR 0.14.2, data.table and parallel. Set `R_LIBS_USER` to the authorized existing library if needed. Python comparison/evaluation uses numpy, pandas, scipy, scikit-learn and matplotlib. All output directories are private: detailed outputs include per-variant PIPs and causal identifiers and must not be committed.

```bash
export GWAS_BENCHMARK_ROOT=/path/to/verified/benchmark
export SUSIE_REVISION_ROOT=/path/to/new/private/susie_erv_true
Rscript scripts/revision_2026-11-05/susie/02_susie_rss_erv_true.R 20
python scripts/revision_2026-11-05/susie/compare_erv.py \
  --benchmark "$GWAS_BENCHMARK_ROOT" --new "$SUSIE_REVISION_ROOT"
python scripts/revision_2026-11-05/susie/05_pip_maf_matched.py
python scripts/revision_2026-11-05/susie/06_localization_precision.py
python scripts/revision_2026-11-05/susie/07_high_pip_window_percentages.py
python scripts/revision_2026-11-05/susie/08_precision_threshold.py
```

The fitter can resume completed loci in this dedicated directory. A settings record checks the fitting options, LD root, locus manifest, marker table and GWAS fingerprints before resumption; it refuses an unverified existing fit directory. Preserve the separately recorded LD checksums because this guard does not rehash every dense LD matrix. It checks every LD matrix, sample count and allele mapping. The comparison validates SNP/locus membership, effect alleles and z-scores against the old run before calculating quantitative differences. It reports both average precision and trapezoidal PR area and the exact random-ranking baseline for each evaluation universe. These quantities must not be conflated.

Use the new ERV_TRUE RSS tables consistently when preparing the revised primary comparison. Label old ERV_FALSE results as the paired sensitivity baseline. Do not merge per-SNP calls between the runs. Only reviewed aggregate comparison tables and figures belong in the public package.

The manuscript's Figure 8 pools SNP–phenotype pairs into one ranking per GWAS input and pipeline. This is distinct from the sensitivity table's mean phenotype-specific AP. Recompute the exact pooled metrics from verified private aggregate PIPs, then generate the figure:

```bash
python scripts/revision_2026-11-05/susie/pooled_finemapping_metrics.py \
  --benchmark "$GWAS_BENCHMARK_ROOT" --rss-run "$SUSIE_REVISION_ROOT" \
  --output "$SUSIE_REVISION_ROOT/results"
python scripts/revision_2026-11-05/susie/plot_figure8_pooled.py \
  --results "$SUSIE_REVISION_ROOT/results" --output NEW_FIGURE_DIRECTORY
```

For a redraw using only distributed aggregate data, use `python finemapping/pr_curve/plot_pooled.py --output NEW_FIGURE_DIRECTORY`. The pooled metric table states each fraction's numerator and denominator, while average precision and trapezoidal PR area remain separate ranking integrals. Source hashes identify the original fixed RSS, new estimated RSS and unchanged LOCO-adjusted inputs. Original realized-PVE strata remain unresolved; reconstructed coefficient deciles are sensitivity-only and do not supply generating SNP effects.
