# Current figures

`python scripts/workflow.py figures --data derived_data --out NEW_DIRECTORY` draws Figures 3/5/6 and the power supplement from the revised canonical association tables, plus the verified LD curve. `plot_aggregates.py` defaults to this revision source. The explicit `--historical` option selects the old strict figure set.

`workflow.py manhattan --profile hpcc --gcta-directory VERIFIED_DIRECTORY --out NEW_DIRECTORY` delegates to the current native-scale inclusive plotting code, including updated SLEMM. It requires complete coverage by default. Use `--allow-incomplete` for labeled review panels. The main plot retains the disclosed display cap; supplements show full ranges. Exact-threshold points are included. Historical strict plots require `--historical`.

Fine-mapping figures use the separate ERV_TRUE canonical tables and [fine-mapping sources](../revision_2026-11-05/susie/README.md). Revised Figure 7 remains dependent on the verified theoretical input. The [figure map](../../docs/figure_map.md) separates current results from historical snapshots.
