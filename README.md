# Livestock GWAS benchmark

This repository supports a benchmark of SNP localization under diffuse polygenic architecture using restricted Duroc pig genotypes and simulated phenotypes. The design analyzes 20,000 animals, 135,358 SNPs, ten phenotypes, 10,000 fixed causal positions and target heritability 0.5. It compares BOLT-LMM, REGENIE with/without LOCO, fastGWA, SLEMM, and FarmCPU/BLINK Kansas and NYC outputs, followed by SuSiE-adj and SuSiE-RSS.

The November 2026 reconstruction follows the author-supplied October 22 manuscript. **Downstream results have been independently regenerated from saved GWAS and fine-mapping outputs. Full model fitting and the upstream Equation 12 calculation have not been reproduced.** Two confirmed manuscript numerical discrepancies and unresolved methodological provenance are retained in the [validation report](repository_validation_report.md).

## Start here

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r environment/requirements-validation.txt
python -m unittest discover -s tests -v
python tests/check_static.py
python scripts/validate_manuscript.py
python tests/check_validation_status.py
```

Static checks also require R on PATH; no R packages are required merely to parse scripts. Figure 1 rendering requires data.table, ggplot2 and patchwork. The unit tests need no proprietary inputs and fit no external GWAS model. CI is configured in [.github/workflows/validation.yml](.github/workflows/validation.yml); its local equivalent was executed, and hosted run status is available on the repository Actions page.

## Contents

| Area | Content |
|---|---|
| `scripts/` | Public analysis code, harmonization, recovered fitting commands and figure generators |
| `config/` | Scientific constants and example machine paths; unresolved parameters remain null |
| `derived_data/` | Non-identifying counts, rates, binned summaries and input fingerprints; real results alongside explicitly named schema templates |
| `figures/reconstructed/` | Eight main-figure redraws, S1–S10 Manhattan supplements, and a central-null Equation 13 supplement, each PDF and PNG |
| `audit/` | Provenance, numerical claim ledger, discrepancies and reference aggregates |
| `example/synthetic/` | Generator for artificial code-validation fixtures; not benchmark evidence |
| `tests/` | Boundary/statistic/schema tests and mock-command tests |

Restricted genotypes, animal phenotypes/IDs, raw GWAS, GRMs, raw LD, fitted SNP PIPs and the manuscript DOCX are not included. Public-data eligibility still requires the authors' data-agreement review before publication; no data-rights approval is implied by this local package.

## Scientific distinctions

Historical Figure 3 uses strict P < 5e-8; spillover uses P ≤ 5e-8. Loci are chromosome-specific floor(bp/500000) bins, not LD-clumped peaks. Central-null spillover uses 1,480 windows. Historical lambda calculations use beta/SE statistics for mixed models and finite-only inverse-CDF statistics for GAPIT; a stable P-based diagnostic is separate. See [exact definitions](docs/statistical_definitions.md).

The simulation executed MPH Gaussian genetic-value draws from a causal-SNP GRM. It generated 1,000 replicates, then retained the first ten phenotypes and first 20,000 animals in historical order. Explicit per-SNP effect draws were not saved. The corrected wrapper retains the original 0.5 GRM weight, seed argument 0 and generation count. It must not be used against the replaced historical genotype path without verifying identity and order.

## Reproduction and limitations

Follow [workflow](docs/workflow.md), [reproducibility](docs/reproducibility.md), [input schemas](docs/input_file_formats.md), [HPCC execution](docs/hpcc_execution.md), [figure map](docs/figure_map.md) and [manuscript-to-code map](docs/manuscript_to_code_map.md).

Unresolved: the original Equation 12 covariance generator, the source of the manuscript maximum-PVE range, the executed GAPIT version (manuscript 3.5 versus retained source identifying 4.0), original genotype/sample ordering and array/imputation metadata. These limitations prevent a claim of complete end-to-end reproducibility. No expensive GWAS/SuSiE fitting was repeated, and the reviewed package is now published on branch `reconstruct/20261104` following explicit user authorization.

Citation metadata from the existing repository is preserved in [CITATION.cff](CITATION.cff). Authors should reconcile the cited preprint title/version with the current manuscript. The existing MIT license applies as recorded; this audit does not establish rights to third-party software or restricted data. See [license review](docs/license_review.md).
