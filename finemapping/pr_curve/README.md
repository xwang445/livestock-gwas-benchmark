# Fine-mapping precision and recall

This is the publishable entry point for the November 5 fine-mapping revision. The original requested workstation `finemapping/pr_curve/` directory could not be accessed in the execution environment; these sources were reconstructed from the verified retained RSS evaluation and the newly executed paired rerun. They are not represented as a recovered original directory.

The full rerun and exact metric calculator are in `scripts/revision_2026-11-05/susie/`. Public aggregate sensitivity tables retain both ERV_FALSE and ERV_TRUE, the phenotype-specific evaluation universes, average precision, trapezoidal PR area and random-ranking baselines. Ordinary in-sample LD, N=20,000, GWAS inputs and loci are identical across this pair. The final revised RSS arm uses ERV_TRUE; the old FALSE run is the paired baseline.

Final manuscript Figure 8 uses a single pooled ranking of SNP–phenotype pairs for each GWAS input and fine-mapping pipeline. Each SNP is counted once per phenotype, retaining its maximum PIP across overlapping loci. The two panels show BOLT-LMM and REGENIE-LOCO, using final ERV_TRUE RSS and unchanged ERV_TRUE SuSiE-adj. Exact pooled average precision and explicit numerator/denominator counts are in `derived_data/revision_2026-11-05/susie/pooled_finemapping_metrics.tsv`. The paired old RSS ERV_FALSE values are identified separately in that table.

Public display vertices contain only pipeline/method labels, precision, recall and PIP cutoffs; they contain no sample identifiers or per-SNP records. The full pooled curve is used for exact average precision; the reduced display table is not an estimator of reported AUPRC. A zero-call endpoint is omitted from display because empirical precision is undefined there. Redraw the final figure without restricted inputs:

```bash
python finemapping/pr_curve/plot_pooled.py --output NEW_FIGURE_DIRECTORY
```

The distributed final figure, caption and source manifest are in `figures/revision_2026-11-05/figure8/`. The original realized-PVE panel is preserved privately; its generating-effect inputs remain unverified and it is not presented as an updated result.

The residual-variance sensitivity analysis uses a different summary: the mean of ten phenotype-specific AP values and phenotype-specific PR curves. `export_and_plot.py` exports every 200th vertex and each endpoint from those sensitivity curves. Its plotted curves are phenotype means interpolated to a common recall grid; they do not reproduce pooled manuscript Figure 8. Use it only for the explicitly labeled sensitivity display:

```bash
python finemapping/pr_curve/export_and_plot.py \
  --full-curves "$SUSIE_REVISION_ROOT/results/precision_recall_curves.tsv.gz" \
  --display-curves NEW_DISPLAY_TABLE.tsv \
  --summary "$SUSIE_REVISION_ROOT/results/sensitivity_comparison_summary.tsv" \
  --out NEW_FIGURE_DIRECTORY
```

For a public sensitivity redraw, omit `--full-curves` and use the distributed sensitivity display table and exact summary paths in the revision record. Detailed fits, truth identifiers, full LD and per-SNP PIPs remain restricted.
