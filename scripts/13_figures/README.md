# Figures

Use `python scripts/workflow.py figures --data derived_data --out ../benchmark-figures` for public aggregate figures and the LD curve. Use `python scripts/workflow.py manhattan --profile mac --out ../benchmark-manhattan` for Figure 4 and S1–S10 from restricted GWAS inputs.

`plot_aggregates.py` consumes aggregate tables. `plot_ld.R` uses data.table, ggplot2 and patchwork. `plot_manhattan.py` plots all tested SNPs in the supplements; the main Figure 4 redraw has a display cap of 31. The uncapped maximum is 107.517 in REGENIE_LOCO phenotype 8.

Each figure is written as PDF and PNG. Use a new output directory. The [figure map](../../docs/figure_map.md) records each source table and final figure path. Equation 13 figures use the cached theoretical predictions.
