# Downstream analysis

Use `python scripts/workflow.py summarize --profile mac --out ../benchmark-rerun` from the repository root. The same command supports `--profile hpcc` on HPCC. Add `--resume` to continue a verified summary run.

`benchmark_core.py` reads the nine observed output schemas and calculates SNP/locus counts, power/FPR, spatial enrichment, spillover, lambda, Jaccard overlap and replicate stability in one pass. Outputs are aggregate tables under the new output directory's `derived_data/`. Per-variant harmonized data and the resumable cache remain restricted.

Use `python scripts/workflow.py equation13 --profile mac --out ../benchmark-equation13` for the existing cached-theory comparison. This calculates 20 prediction-only bins and intervals across the ten replicate means. It does not regenerate Equation 12.

Definitions are in [statistical_definitions.md](../../docs/statistical_definitions.md); result links are in [derived_data](../../derived_data/README.md). Legacy R analyses and inverse-regression PVE utilities are archived and are not alternate production commands.
