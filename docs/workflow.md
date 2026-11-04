# Supported workflow

Start with `python scripts/workflow.py --help`. This is the supported entry point. It uses the checked configuration and calls the analysis modules below. Historical alternatives are preserved in `archive/legacy_workflow` and are not part of this workflow.

## Check the published results

```bash
python scripts/workflow.py check
```

This checks configuration, source-table fingerprints, schemas, script syntax, synthetic tests, and the corrected manuscript claims. It needs Python dependencies from `environment/requirements-validation.txt` and R for syntax checks. It does not require restricted inputs or fit a model.

## Recalculate from saved inputs

Install the private machine profile described in [configuration](../config/README.md). Run each command on the machine holding the files. Choose new output directories outside the historical data tree.

| Command | Profile | Inputs | Outputs |
|---|---|---|---|
| `summarize` | Mac or HPCC | 90 saved GWAS panels, causal SNP list, effect windows | SNP/locus counts, power/FPR, spillover, lambda, overlap and stability |
| `equation13` | Mac | Saved GWAS and cached theoretical predictions | Replicate-bin means, intervals and cache checks |
| `ld` | HPCC | Saved pairwise LD | 1,000 distance-bin counts and means |
| `fine-mapping` | HPCC | Saved converged PIPs, locus diagnostics and MAF annotations | Per-phenotype and pooled fine-mapping summaries |
| `figures` | Either | Public aggregate tables | Aggregate figures and LD curve |
| `manhattan` | Mac or HPCC | Saved GWAS panels | Figure 4 and S1–S10 |

```bash
python scripts/workflow.py summarize --profile mac --out ../benchmark-rerun
python scripts/workflow.py equation13 --profile mac --out ../benchmark-equation13
python scripts/workflow.py figures --data derived_data --out ../benchmark-figures
python scripts/workflow.py manhattan --profile mac --out ../benchmark-manhattan
```

On HPCC, with the local profile installed:

```bash
python scripts/workflow.py ld --profile hpcc --out ../benchmark-ld-check
python scripts/workflow.py fine-mapping --profile hpcc --out ../benchmark-finemapping-check
```

Only `summarize` accepts `--resume`. Its cache is reused only when input, code and design hashes match. The `.restricted_cache` stays local. The workflow refuses outputs inside the historical source tree and does not replace committed result tables.

`benchmark_core.py` implements harmonization and summary definitions; the wrapper does not create a second statistical implementation. Equation 13 uses the existing prediction file. Equation 12 and maximum-SNP PVE provenance are deferred to Professor Jicai Jiang at the author's request.

## Upstream study stages

The retained fitting commands are documented in [genotype preparation](../scripts/01_genotype_preparation/README.md), [simulation](../scripts/03_simulation/README.md), [GWAS](../scripts/04_gwas/README.md), and [fine-mapping](../scripts/12_finemapping/README.md). These are for a separately allocated new run. They are not invoked by the saved-result workflow and were not refitted during this update.

Use [statistical definitions](statistical_definitions.md) for thresholds and boundaries, [input formats](input_file_formats.md) for schemas, and [figure map](figure_map.md) for table-to-figure links.
