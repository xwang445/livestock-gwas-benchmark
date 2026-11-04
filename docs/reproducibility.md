# Reproducibility and exact commands

Use a new output workspace; treat the historical benchmark as read-only. Paths below are placeholders. The new Python summary command requires the original directory/file naming described in input_file_formats.md. It writes only aggregate tables and a private resumable cache. Do not commit the cache, harmonized variant exports or raw inputs.

```bash
export RESTRICTED_ROOT=/path/to/authorized/W2Mb_N20k
export REBUILD=/path/to/new/reconstruction
mkdir -p "$REBUILD"
python scripts/benchmark_core.py summarize --root "$RESTRICTED_ROOT" --out "$REBUILD/derived_data" --cache "$REBUILD/private_cache"
python scripts/benchmark_core.py equation13 --root "$RESTRICTED_ROOT" --theory "$RESTRICTED_ROOT/equation13_analysis/framework_ncp_variant.tsv" --out "$REBUILD/derived_data/equation13_summary"
python scripts/12_finemapping/verify_saved_results.py --root "$RESTRICTED_ROOT" > "$REBUILD/finemapping_verification.json"
python scripts/export_verification.py --kind finemapping --json "$REBUILD/finemapping_verification.json" --out "$REBUILD/derived_data/finemapping_summary"
python scripts/02_ld_analysis/verify_ld_means.py --input "$RESTRICTED_ROOT/geno_data/ld_analysis/ld_10Mb.ld.gz" > "$REBUILD/ld_verification.json"
python scripts/export_verification.py --kind ld --json "$REBUILD/ld_verification.json" --out "$REBUILD/derived_data/ld_summary"
Rscript scripts/13_figures/plot_ld.R "$REBUILD/derived_data" "$REBUILD/figures"
python scripts/13_figures/plot_aggregates.py --data "$REBUILD/derived_data" --out "$REBUILD/figures"
python scripts/13_figures/plot_manhattan.py --root "$RESTRICTED_ROOT" --out "$REBUILD/manhattan"
```

LD/fine-mapping verification should run on the host holding those files. During this audit scripts were streamed over authorized SSH and only aggregate JSON was returned. No raw LD/genotype was copied to the Mac. Source inputs were not modified.

The committed `derived_data` tables were generated with Python 3.10.4, NumPy 2.2.6, pandas 2.3.3, SciPy 1.14.1, matplotlib 3.10.9 and PyYAML 6.0.3; LD plotting used R 4.3.2. Python dependencies are pinned for validation. Production tool versions are separate in config/methods.yaml and the fine-mapping README. The system SciPy binary failed to import; the working environment used for this audit is recorded privately in PROJECT_STATE.md. Dependencies in other projects were not modified.

Cache reuse requires the source-file SHA256, window/truth SHA256, implementation SHA256 and expected SNP count to match. Original results are immutable. Changing code invalidates caches intentionally. Saved aggregates can be checked without restricted inputs using:

```bash
python scripts/validate_manuscript.py
python tests/check_validation_status.py
```

`audit/reference_values` contains historical non-identifying comparison tables with source hashes; these are references, not newly regenerated results. `audit/manuscript_validation.tsv` records expected/observed values and tolerances. Exact counts use zero tolerance, historical floating summaries use 1e-10 (spatial ratios 1e-7), Equation 13 bins/intervals 1e-9, and manuscript rounded percentages use half the displayed unit. Alpha grids are joined by log10 threshold rounded to ten decimal places solely to reconcile text serialization.

## Limits of the evidence

The same raw outputs can reproduce a historical coding convention without proving that convention is statistically ideal. In particular finite-only GAPIT lambda discards very small P tails after rounding, central-null endpoints differ between analyses, and Figure 3 differs from spillover at P exactly 5e-8. Both primary conventions and diagnostics are explicit.

Equation 13 reproduction verifies cached NCP algebra, joins, all 1,600 replicate-bin means and their aggregate intervals. It does not recover the Equation 12 covariance computation. No exact upstream rerun command can honestly be provided until its generating formula implementation and inputs are identified. The current `calculate_ncp.R` wrapper only performs the recovered downstream calculation and explicitly says so.

The finite-SNP PVE claim is unresolved. Recovered LSMR reconstruction of effects from genetic values is not a recovery of original sampled causal effects. Do not substitute it to force agreement.

Manhattan plots are visual redraws; colors/layout differ from historical artwork. Their counts are checked against regenerated strict-threshold counts, and full-range supplements prevent the main display cap from hiding the existence of extreme results. No editable manuscript or author approval was created by this task.

## Verified configurations

`config/paths.yaml` resolves the exact result files and restricted-source filenames; `config/methods.yaml` lists the distinct Mac and HPCC GWAS layouts. `config/verified_results.yaml` is generated directly from executed tables, with source hashes and the validation status. Run `python scripts/build_result_config.py --check` to detect stale results configuration.

The actual machine roots are installed locally in Git-ignored `config/local_paths.yaml`. `python scripts/check_config.py --profile mac --check-inputs` resolves and checks the 90 local panels and shared inputs. `--profile hpcc` resolves HPCC paths without SSH or computation. The Equation 13 cache is verified on the Mac; phenotype/LD/fine-mapping inputs are verified on HPCC. See config/README.md and audit/config_path_verification.tsv. Files absent on a machine are explicitly excluded from that profile.
