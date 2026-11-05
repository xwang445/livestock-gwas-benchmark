# Fine-mapping precision and recall

This is the publishable entry point for the November 5 RSS comparison. The original requested workstation `finemapping/pr_curve/` directory could not be accessed in the execution environment; these sources were reconstructed from the verified retained RSS evaluation and the newly executed paired rerun. They are not represented as a recovered original directory.

The full rerun and exact metric calculator are in `scripts/revision_2026-11-05/susie/`. Public aggregate sensitivity tables retain both ERV_FALSE and ERV_TRUE, the phenotype-specific evaluation universes, average precision, trapezoidal PR area and random-ranking baselines. Ordinary in-sample LD, N=20,000, GWAS inputs and loci are identical across this pair. The final revised RSS arm uses ERV_TRUE; the old FALSE run is the paired baseline.

`export_and_plot.py` exports every 200th PR vertex and each endpoint without SNP IDs or thresholds. Its curves are for display only; exact areas come from the full original vertices used by `compare_erv.py`. Plotted curves are phenotype means interpolated to a common recall grid. Do not calculate reported AUPRC from this reduced display table.

```bash
python finemapping/pr_curve/export_and_plot.py \
  --full-curves "$SUSIE_REVISION_ROOT/results/precision_recall_curves.tsv.gz" \
  --display-curves NEW_DISPLAY_TABLE.tsv \
  --summary "$SUSIE_REVISION_ROOT/results/sensitivity_comparison_summary.tsv" \
  --out NEW_FIGURE_DIRECTORY
```

For a public redraw, omit `--full-curves` and use the distributed display table and exact summary paths in the revision record. Detailed fits, truth identifiers and per-SNP PIPs remain restricted.
