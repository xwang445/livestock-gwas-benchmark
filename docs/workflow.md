# Supported current workflow

`python scripts/workflow.py` defaults to the November 5 revision. All active association decisions use **P ≤ 5 × 10⁻⁸**. Native REGENIE LOG10P is compared directly with −log10(5 × 10⁻⁸), without exponentiating first or adding a tolerance. The revised SLEMM qf90/window3000 files and ten-configuration registry supply current summaries and plots.

## Public checks and figures

The root summary folders now contain all ten configurations. Rebuild these views from committed source aggregates with `python scripts/workflow.py refresh-results`. The check verifies their contents and complete method/phenotype coverage before validating the configuration. See [the root-summary update](revision_2026-11-13.md).

```bash
python scripts/workflow.py check
python scripts/workflow.py figures --data derived_data --out NEW_FIGURE_DIRECTORY
```

The check validates current configuration/table consistency, exact-threshold regressions, syntax, aggregate schemas and the separately labeled historical manuscript ledger. The default figure command reads `derived_data/revision_2026-11-05/association`, draws Figures 3/5/6 and the power supplement, and regenerates Figure 1 from its verified LD table. This command retains the November 5 diagnostic layouts. The November 12 submission stage below redraws Figure 7 from retained bin summaries and Figure 8 from final RSS plus recovered-PVE aggregates. Fine-mapping figures have their own [canonical workflow](../scripts/revision_2026-11-05/susie/README.md).

## Submission artwork from public aggregates

```bash
python scripts/workflow.py submission --data derived_data --out NEW_SUBMISSION_DIRECTORY
```

This stage redraws Figures 1–3 and 5–8, S11–S20, and the pooled four-column Table S1 from [committed aggregate sources](../derived_data/submission_2026-11-12). The distributed [submission artwork](../figures/submission_2026-11-12) also contains Figure 4/S1–S10. Regenerating those Manhattan figures requires the explicit restricted-input command in the [submission script instructions](../scripts/submission_2026-11-12/README.md). Their raw variant-level caches are excluded from Git. The older `figures` stage retains its November 5 diagnostic layouts.

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

For publication at the root summary paths, the aggregate support exporter supplies QC counts, source fingerprints and cross-replicate locus-set overlaps:

```bash
python scripts/revision_2026-11-05/association/derive_public_summary_support.py \
  --qc "$PRIVATE_ASSOCIATION_QC_TSV" --detected-loci "$PRIVATE_DETECTED_LOCI_TSV" \
  --input-manifest "$PRIVATE_ASSOCIATION_INPUT_MANIFEST_TSV" \
  --benchmark-root "$GWAS_BENCHMARK_ROOT" --gcta-root "$GCTA_OUTPUT_DIRECTORY" \
  --output "$NEW_AGGREGATE_SUPPORT_DIRECTORY"
```

Reviewed support files are stored with the dated canonical association tables. `refresh-results` reads those public aggregates and never opens restricted inputs.

Add `--individual` to the current `manhattan` command to write one PDF and PNG for each completed configuration × phenotype, in addition to the composite figures. With all ten GCTA phenotypes validated, this produces 100 individual PDF/PNG pairs plus the eleven existing composite PDF/PNG pairs. Each individual uses the same loaded SNP coordinates, native inclusive decisions, chromosome colors and threshold line; the vertical axis shows that configuration’s full range. `manhattan_manifest.tsv` records relative `individual_pdf` and `individual_png` paths. Pending rows have no individual files.

The existing snapshot has 100 planned manifest rows: 100 completed GWAS panels and zero pending GCTA panels. Those are panel counts, rather than counts of composite image files. The ten supplementary composite figures contain ten method rows each; pooled Figure 4 is a separate composite. Use the complete command without `--allow-incomplete` for the final set of 100 individual analyses.

`equation13 --profile mac --out NEW_DIRECTORY` reads the verified cached theory and updated SLEMM input. The theoretical source remains unavailable in this execution environment; the command cannot supply missing inputs or reconstruct Equation 12.

## Explicit historical audit

The old strict Figure 3/Figure 5 results and original SLEMM/RSS fits remain immutable snapshots in `derived_data/historical_2026-11-04/`. Historical figure commands given the root `derived_data` path resolve that archive automatically, including its LD curve. To reproduce them deliberately:

```bash
python scripts/workflow.py summarize --historical --profile hpcc --out NEW_HISTORICAL_AUDIT
python scripts/workflow.py figures --historical --data derived_data --out NEW_HISTORICAL_FIGURES
python scripts/workflow.py manhattan --historical --profile hpcc --out NEW_HISTORICAL_MANHATTAN
```

Only this historical summary supports `--resume`, using its input/code/design hash-checked cache. Its strict cutoff is an audit definition and never supplies current primary columns. The legacy core refuses a summary invocation without explicit historical opt-in. Existing historical outputs are not overwritten or relabeled as current results.

Upstream fitting commands remain separate from these saved-result stages. See [statistical definitions](statistical_definitions.md), [revision status](revision_2026-11-05.md), and [figure map](figure_map.md) for the scientific scope and outstanding dependencies.
