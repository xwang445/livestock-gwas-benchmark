# Reproducibility

The supported commands are in [workflow.md](workflow.md). Use `scripts/workflow.py` for public checks, saved-result reanalysis and figure generation. Configuration is in `config/`; actual Mac/HPCC roots are in the Git-ignored `config/local_paths.yaml`.

## Recorded execution

The public tables were regenerated from 90 saved GWAS panels, cached Equation 13 predictions, saved converged fine-mapping fits and pairwise LD. Python 3.10.4 with the pinned scientific dependencies produced the aggregate results; R 4.3.2 produced the LD curve. No GWAS or SuSiE model was refitted during reconstruction or consolidation.

The [run registry](../audit/run_registry.tsv), [input hashes](../derived_data/input_manifest.tsv), [artifact hashes](../audit/artifact_manifest.tsv), and [current validation ledger](../audit/manuscript_validation.tsv) record the execution. The [original manuscript ledger](../audit/historical_manuscript_validation.tsv) preserves the two discrepancies before text correction.

## Reanalysis and resumption

Use a new output directory outside the source tree. The summary stage accepts `--resume`; unchanged jobs use verified cache entries. Cache reuse checks input SHA256, window/truth SHA256, implementation SHA256 and SNP count. Retain restricted caches locally and publish only reviewed aggregate outputs.

After a new analysis has been reviewed and its tables installed, regenerate the public snapshot and hashes:

```bash
python scripts/build_result_config.py
python scripts/build_artifact_manifest.py
python scripts/workflow.py check
```

## Provenance

The simulation description and sample chain are documented in [scientific provenance](scientific_provenance.md). The retained chip dataset matches the saved simulation order and checked GRM entries. GWAS and fine-mapping rows are aligned by IID. Original source data and the October manuscript remain unchanged; corrections are in a separate local manuscript copy.

Equation 12 and maximum-SNP PVE provenance are deferred to Professor Jicai Jiang. The cached Equation 13 comparison remains available, with that scope explicit. GAPIT's retained source identifies version 4.0 while the manuscript reports 3.5; this version question remains in the author review record.
