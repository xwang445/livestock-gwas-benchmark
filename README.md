# Livestock GWAS benchmark

We compare how well widely used GWAS methods locate causal signals under diffuse polygenic architecture in livestock. The benchmark uses Duroc pig genotypes and simulated phenotypes, followed by SuSiE-adj and SuSiE-RSS fine-mapping.

This repository contains the analysis scripts, result tables, figures, and configuration files. The [November 12 update](docs/revision_2026-11-12.md) adds the submission figures, pooled Table S1, and author-supplied recovered-effect PVE analyses. All ten SLEMM and GCTA MLMA phenotypes and all 100 configuration–phenotype panels are complete. Both final SuSiE pipelines estimate residual variance. The original November 4 results and [November 5–9 analysis revision](docs/revision_2026-11-05.md) remain available in their established paths.

[Publication](#publication) · [Result tables](derived_data/README.md) · [Current figures](figures/submission_2026-11-12) · [Configuration](config/README.md) · [Analysis commands](docs/reproducibility.md)

## Publication

Wang, X., Wang, J., Tiezzi, F., Huang, Y., Huang, W., Maltecca, C., & Jiang, J. (2026). **[Widely used GWAS methods can be poorly suited to SNP-level localization under diffuse polygenic architecture in livestock](https://www.biorxiv.org/content/10.64898/2026.07.20.739621v1)**. *bioRxiv* preprint. DOI: [10.64898/2026.07.20.739621](https://doi.org/10.64898/2026.07.20.739621).

Please cite this paper if you use the code or results from this repository.

## Study design

| Item | Setting |
|---|---|
| Population | Duroc pigs |
| Animals in the GWAS benchmark | 20,000 |
| SNPs | 135,358 across 18 autosomes |
| Causal SNPs | 10,000 fixed positions |
| Analyzed phenotypes | 10 simulated replicates |
| Target heritability | 0.5 |
| Effect windows | 2 Mb |
| Null gaps | 4 Mb, with the central 2 Mb used for central-null analyses |
| Locus size | 500 kb |

The original comparison includes nine GWAS configurations:

- BOLT-LMM
- REGENIE with and without leave-one-chromosome-out (LOCO) analysis
- fastGWA
- SLEMM
- FarmCPU: Kansas and NYC outputs
- BLINK: Kansas and NYC outputs

Kansas and NYC refer to the GAPIT outputs before and after local refinement. The nine configurations across ten phenotypes give 90 GWAS result files. Full settings are in [analysis_parameters.yaml](config/analysis_parameters.yaml) and [methods.yaml](config/methods.yaml).

The revision registry adds GCTA MLMA, giving seven method families and ten configurations. All 100 configuration–phenotype summaries are validated in [the registry and canonical tables](docs/revision_2026-11-05.md).

## Workflow

Use `python scripts/workflow.py --help`. Default commands use the current canonical revision, updated SLEMM and **P ≤ 5 × 10⁻⁸**, with direct comparison of native LOG10P. The old strict workflow requires explicit `--historical`; it never supplies current primary results.

The [threshold audit](docs/threshold_consistency.md) checks equality-boundary calls, canonical counts and the unchanged original SuSiE locus manifest.

```text
Saved GWAS + causal SNPs + windows → counts, power/FPR, spillover, lambda, overlap
Saved GWAS + cached theory        → Equation 13 comparisons
Saved LD and fine-mapping fits     → LD and fine-mapping summaries
Aggregate tables                  → main and supplementary figures
```

| Step | Command | Where |
|---|---|---|
| Check the public package | `python scripts/workflow.py check` | Any configured Python/R environment |
| Redraw submission plots and Table S1 from public aggregates | `submission --data derived_data --out NEW_DIRECTORY` | Any configured Python environment; Manhattan sources remain restricted |
| Summarize the saved GWAS panels | `summarize --profile hpcc --simulation VERIFIED_SIMULATION --genotype-prefix VERIFIED_PREFIX --gcta-directory VERIFIED_GCTA --out NEW_DIRECTORY` | Machine holding the verified inputs |
| Compare cached theory and observed statistics | `equation13 --profile mac --out ../benchmark-equation13` | Mac |
| Verify saved LD / fine-mapping results | `ld`; or `fine-mapping --rss-directory COMPLETED_ERV_TRUE_RUN`, with profile/output arguments | HPCC |
| Draw revised association figures and LD | `figures --data derived_data --out NEW_DIRECTORY` | Mac or HPCC |
| Draw Manhattan plots | `manhattan --profile hpcc --gcta-directory VERIFIED_GCTA --out NEW_DIRECTORY` | Complete configurations required; diagnostic `--allow-incomplete` is explicit |

For rows after the first, prefix the command with `python scripts/workflow.py`. The [workflow guide](docs/workflow.md) gives complete commands and inputs. Upstream genotype preparation, simulation and fitting scripts are retained for a separately allocated new run. Superseded scripts are preserved in [archive](archive/README.md) and are not part of the active workflow.

MPH generated genetic values from the fixed causal-SNP relationship matrix and added independent residuals. The benchmark uses the first ten of 1,000 generated phenotypes and the first 20,000 simulation rows. GWAS and fine-mapping files are aligned by individual ID. The [provenance record](docs/scientific_provenance.md) documents the source code and sample-chain checks.

## Current and historical results

Current association and residual-variance comparisons are in [derived_data/revision_2026-11-05](derived_data/revision_2026-11-05). The [November 11 recovered-effect results](derived_data/revision_2026-11-11) and [November 12 submission sources](derived_data/submission_2026-11-12) add lead-PVE summaries, the updated Figure 8, and pooled Table S1. Table S1 reports SLEMM effect/null percentages of 93.75/6.25 and GCTA MLMA 93.33/6.67; central-null percentages are zero for both. These pooled ratios differ from the phenotype-mean percentages retained in the earlier analysis tables.

`config/verified_results.yaml` derives its primary blocks from committed tables and fingerprints their sources. Its `historical_*` blocks preserve the older strict snapshot. The [workflow guide](docs/workflow.md) describes the current and audit commands.

The historical root set in [derived_data](derived_data/README.md) contains 27 populated TSV files; current revision tables are separately indexed above. The historical set covers:

- Significant SNPs and loci, power, false-positive rates, and spatial enrichment
- Central-null spillover and genomic inflation
- Locus overlap between methods and stability across replicates
- Equation 13 comparisons, fine-mapping summaries, LD decay, and heritability estimates from model logs

The historical [figures/reconstructed](figures/reconstructed) set contains Figures 1–8, Manhattan supplements S1–S10, and the central-null Equation 13 supplement in PDF and PNG formats. The [figure map](docs/figure_map.md) links each figure to its inputs and script.

The preserved strict-audit REGENIE_LOCO maximum is 20,461 significant SNPs. The largest −log10(P) is 107.517 for REGENIE_LOCO phenotype 8. [Independent historical checks](audit/manuscript_number_checks.tsv) support both values; current inclusive counts come from the revision tables.

The machine-readable [results summary](config/verified_results.yaml) records the historical calculated values and source-file hashes. [Statistical definitions](docs/statistical_definitions.md) document the thresholds, window boundaries, and denominators used in each analysis. The [November 5 revision record](docs/revision_2026-11-05.md) distinguishes new analyses from these preserved historical results.

## Set up and check the repository

Use Python 3.10 or later and R. The checks below use public tables and synthetic test data.

```bash
git clone https://github.com/xwang445/livestock-gwas-benchmark.git
cd livestock-gwas-benchmark
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/requirements-validation.txt

python scripts/workflow.py check
```

The LD figure script uses the R packages `data.table`, `ggplot2`, and `patchwork`. Software versions and dependencies are listed in [environment](environment).

## Paths and configuration

| File | Purpose |
|---|---|
| [analysis_parameters.yaml](config/analysis_parameters.yaml) | Study design and analysis settings |
| [methods.yaml](config/methods.yaml) | Method settings and GWAS filenames for Mac and HPCC |
| [paths.yaml](config/paths.yaml) | Result locations and source filenames |
| [verified_results.yaml](config/verified_results.yaml) | Calculated results and source hashes |
| `config/local_paths.yaml` | Local Mac/HPCC roots; excluded from Git |

On a machine with the local path profile installed, check the Mac inputs or display the HPCC paths with:

```bash
python scripts/check_config.py --profile mac --check-inputs
python scripts/check_config.py --profile hpcc
```

See [configuration instructions](config/README.md), [input formats](docs/input_file_formats.md), and [HPCC execution](docs/hpcc_execution.md) for working with the source files.

## Repository structure

```text
config/                 Study settings, paths, and calculated results
derived_data/           Aggregate result tables
scripts/                Analysis and plotting code
figures/reconstructed/  Main and supplementary figures
docs/                   Methods and execution guides
environment/            Dependencies and software versions
tests/                  Code and result checks
example/                Synthetic test data and schema examples
audit/                  Source records and numerical validation
archive/                Superseded workflow, retained for reference
```

## Data and citation

Individual-level genotypes, phenotypes, animal IDs, raw GWAS results, GRMs, pairwise LD files, and per-SNP fine-mapping results are restricted and are not distributed here. The public repository provides code and aggregate results. See [data availability](docs/data_availability.md).

Citation information is in [CITATION.cff](CITATION.cff). Repository code is distributed under the [MIT license](LICENSE).

The [workflow documentation](docs/workflow.md) and [manuscript-to-code map](docs/manuscript_to_code_map.md) provide the detailed analysis record.
