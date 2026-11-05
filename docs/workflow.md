# Supported current workflow

`python scripts/workflow.py` defaults to the November 5 revision. All active association decisions use **P ≤ 5 × 10⁻⁸**. Native REGENIE LOG10P is compared directly with −log10(5 × 10⁻⁸), without exponentiating first or adding a tolerance. The revised SLEMM qf90/window3000 files and ten-configuration registry supply current summaries and plots.

## Public checks and figures

```bash
python scripts/workflow.py check
python scripts/workflow.py figures --data derived_data --out NEW_FIGURE_DIRECTORY
```

The check validates current configuration/table consistency, exact-threshold regressions, syntax, aggregate schemas and the separately labeled historical manuscript ledger. The default figure command reads `derived_data/revision_2026-11-05/association`, draws Figures 3/5/6 and the power supplement, and regenerates Figure 1 from its verified LD table. It does not silently redraw the unavailable revised Figure 7 or substitute old RSS results. Fine-mapping figures have their own [canonical workflow](../scripts/revision_2026-11-05/susie/README.md).

## Restricted-input reanalysis

Install an ignored machine profile as described in [configuration](../config/README.md). Use new output directories outside the source tree. Supply verified external inputs explicitly; no workstation paths are inferred.

```bash
python scripts/workflow.py summarize --profile hpcc \
  --simulation "$ORIGINAL_SIMULATION_DIRECTORY" \
  --genotype-prefix "$VERIFIED_GENOTYPE_PREFIX" \
  --gcta-directory "$GCTA_OUTPUT_DIRECTORY" --out NEW_SUMMARY_DIRECTORY
python scripts/workflow.py manhattan --profile hpcc \
  --gcta-directory "$GCTA_OUTPUT_DIRECTORY" --out NEW_MANHATTAN_DIRECTORY
python scripts/workflow.py ld --profile hpcc --out NEW_LD_CHECK_DIRECTORY
python scripts/workflow.py fine-mapping --profile hpcc \
  --rss-directory "$COMPLETED_ERV_TRUE_RSS_RUN" --out NEW_FINEMAPPING_CHECK_DIRECTORY
```

Current summaries are written under `NEW_SUMMARY_DIRECTORY/tables`. Missing GCTA outputs are recorded explicitly, and ten-phenotype means are withheld until all ten are present. Manhattan plots require complete configurations by default; `--allow-incomplete` permits clearly labeled diagnostic panels. Fine-mapping verification requires the completed TRUE run and checks its fit-mode field; it preserves the original SuSiE locus manifest and does not refit models.

Add `--individual` to the current `manhattan` command to write one PDF and PNG for each completed configuration × phenotype, in addition to the composite figures. Once all ten GCTA phenotypes validate, this produces 100 individual PDF/PNG pairs plus the eleven existing composite PDF/PNG pairs. Each individual uses the same loaded SNP coordinates, native inclusive decisions, chromosome colors and threshold line; the vertical axis shows that configuration’s full range. `manhattan_manifest.tsv` records relative `individual_pdf` and `individual_png` paths. Pending rows have no individual files.

The existing snapshot has 100 planned manifest rows: 91 completed GWAS panels and nine pending GCTA panels. Those are panel counts, rather than counts of composite image files. The ten supplementary composite figures contain ten method rows each; pooled Figure 4 is a separate composite. Use the complete command without `--allow-incomplete` for the final set of 100 individual analyses.

`equation13 --profile mac --out NEW_DIRECTORY` reads the verified cached theory and updated SLEMM input. The theoretical source remains unavailable in this execution environment; the command cannot supply missing inputs or reconstruct Equation 12.

## Explicit historical audit

The old strict Figure 3/Figure 5 results and original SLEMM/RSS fits remain immutable snapshots. To reproduce them deliberately:

```bash
python scripts/workflow.py summarize --historical --profile hpcc --out NEW_HISTORICAL_AUDIT
python scripts/workflow.py figures --historical --data derived_data --out NEW_HISTORICAL_FIGURES
python scripts/workflow.py manhattan --historical --profile hpcc --out NEW_HISTORICAL_MANHATTAN
```

Only this historical summary supports `--resume`, using its input/code/design hash-checked cache. Its strict cutoff is an audit definition and never supplies current primary columns. The legacy core refuses a summary invocation without explicit historical opt-in. Existing historical outputs are not overwritten or relabeled as current results.

Upstream fitting commands remain separate from these saved-result stages. See [statistical definitions](statistical_definitions.md), [revision status](revision_2026-11-05.md), and [figure map](figure_map.md) for the scientific scope and outstanding dependencies.
