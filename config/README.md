# Verified configuration

`analysis_parameters.yaml` defines one current rule, P ≤ 5 × 10⁻⁸, with direct native-log comparison. Its separate `historical_audit` section applies only to explicit legacy reproduction. `methods.yaml` contains the ten current configurations, including updated SLEMM and partial GCTA MLMA; `historical_output_patterns` retains the original SLEMM paths. `revision_methods.tsv` defines the method-family grouping and plotting colors.

`paths.yaml` points primary result and figure paths to the current revision. The previous 27-table paths remain under `historical_result_paths`. `verified_results.yaml` is generated from the actual canonical tables: primary blocks describe the current revision, and `historical_*` blocks preserve the old strict snapshot. Planned configurations and validated completed panels are reported separately.

The ignored `config/local_paths.yaml` contains `repository_root` and `profiles.mac` / `profiles.hpcc`. Each profile requires `benchmark_root` and can supply these additional verified revision inputs:

```yaml
simulation_directory: VERIFIED_ORIGINAL_SIMULATION_DIRECTORY
genotype_prefix: VERIFIED_FULL_PANEL_PREFIX
gcta_directory: VALIDATED_GCTA_OUTPUT_DIRECTORY
rss_directory: COMPLETED_ERV_TRUE_RSS_RUN_DIRECTORY
```

These symbolic labels illustrate keys only; install actual authorized paths privately. Command-line arguments override these optional keys. The simulation directory contains the saved chromosome lengths and `data_20000` phenotype files. The genotype prefix is supplied explicitly rather than substituting an unverified panel. Current fine-mapping verification refuses the old FALSE run unless `--historical` is supplied.

```bash
python scripts/check_config.py
python scripts/check_config.py --profile hpcc --check-inputs
python scripts/build_result_config.py --check
```

The checker does not launch fits or SSH. Restricted stages check their relevant inputs rather than failing on unrelated unavailable resources. Missing GCTA outputs are optional for independent reanalysis and are reported as missing, never counted as completed. Publication-quality complete Manhattan comparisons still require every planned panel; diagnostic plots require explicit `--allow-incomplete`.

Regenerate the result snapshot and artifact hashes after validated canonical tables change. Never edit a result to match an expected conclusion. The original provenance questions and historical ledgers remain separate scientific evidence.
