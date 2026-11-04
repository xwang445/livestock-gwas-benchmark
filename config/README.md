# Verified project configuration

| File | Contents |
|---|---|
| [analysis_parameters.yaml](analysis_parameters.yaml) | Recovered scientific design, thresholds and seeds. Unknown historical settings remain explicitly null with explanations. |
| [methods.yaml](methods.yaml) | All nine configurations, recovered versions/settings, and verified source filenames for all ten phenotypes on Mac and HPCC. |
| [paths.yaml](paths.yaml) | Real repository-relative paths to all 27 populated TSV files and evidence, plus exact restricted-input filenames and their verified machine availability. |
| [verified_results.yaml](verified_results.yaml) | Actual method counts, spillover, lambda, h2 summaries, fine-mapping totals, LD and Equation 13 checks generated from executed tables, with source hashes. |
| `local_paths.yaml` (Git-ignored) | Installed locally with the actual Mac/HPCC source roots, SSH host and repository location. No credentials or dataset contents are stored. |

The old `/path/to/...` example has moved to [example/schemas/example_paths.yaml](../example/schemas/example_paths.yaml). Published code resolves real filenames relative to each machine profile, so it does not depend on a particular user's absolute paths.

## Check and use

```bash
python scripts/check_config.py
python scripts/build_result_config.py --check
python scripts/check_config.py --profile mac --check-inputs
python scripts/check_config.py --profile hpcc
```

The Mac check confirms 90 GWAS files, causal/window metadata and the cached theory file. The HPCC profile resolves its exact paths without attempting SSH; run `--check-inputs` on the host holding that root for a fresh filesystem check. This reconstruction verified HPCC paths through read-only SSH on 2026-11-04. [Verification records](../audit/config_path_verification.tsv) distinguish present and absent inputs.

The cached Equation 13 table is on the Mac. Phenotypes, raw pairwise LD, fine-mapping saved fits and the 20,000-animal / 135,358-SNP fine-mapping genotype subset are on HPCC. The two profiles are deliberately different. No original simulation genotype prefix is activated because its historical identity/order remains unresolved.

To regenerate the results configuration after a genuinely new executed analysis:

```bash
python scripts/build_result_config.py
python scripts/build_artifact_manifest.py
```

Do not edit result values to match a preferred conclusion. The generated configuration retains the two manuscript discrepancies and upstream Equation 12 limitation. It is a results snapshot, not a set of tuning targets. The checker reads metadata only and does not run a model, transfer data or launch remote jobs.
