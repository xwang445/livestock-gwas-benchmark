# Project state

Status: reconstruction/review package published to GitHub main following the user request for real results on the default branch on 2026-11-04. End-to-end scientific reproduction remains PARTIAL. No manuscript deposit, submission, or external message. Original data, original manuscript, original local main, and historical HPCC outputs were preserved.

## Completed work and outputs

- Audited 590 local and 502,804 HPCC files; inspected the author-confirmed October manuscript and selected historical code.
- Recomputed all 90 GWAS downstream panels and 40 Equation 13 jobs. Final source-write-guard rebuild produced 20 byte-identical TSVs.
- Recounted 40 fine-mapping method/pipeline/phenotype groups on HPCC; streamed and independently checked 1,000 LD bins. No expensive model fits were rerun.
- Numerical ledger: 95 PASS, 2 FAIL, 4 AMBIGUOUS, 3 NOT_RUN, 1 REQUIRES_RECOMPUTATION. Two FAIL rows describe manuscript count/statistic discrepancies, not test failures.
- Generated 19 figures (eight main, S1–S10, central-null Equation 13 supplement), each PDF and PNG; inspected their layout and verified PDF/PNG files.
- All 12 unit/synthetic/mock tests and shell/Python/R/YAML/schema checks passed, including from this installed clone. GitHub-hosted CI is configured; its current status is available on the repository Actions page. Local validation results are recorded above.
- Detailed outputs: derived_data/, figures/reconstructed/, audit/manuscript_validation.tsv, audit/run_registry.tsv, docs/stage_status.md, repository_validation_report.md, docs/author_review.md.

## Local commits

Public baseline 14ff11e; original local main ec6cb32 remains preserved.

- 32ffb12: initial inventory and manuscript provenance.
- 00e08b0: recovered runners, harmonization and downstream implementation.
- 1deaa16: verified aggregate results and 19 figure redraws.
- b8524cf: synthetic example, tests and CI.
- d49a978: documentation and scientific review records.
- Publication status: the user authorized the push after reviewing the completion report. The five reviewed commits through d49a978 were pushed successfully; this status record is a documentation-only follow-up. The remote review branch is reconstruct/20261104.

Active jobs: NONE. Remote verification jobs and local rendering/rebuilds exited successfully.

## Failures and corrections

System Python SciPy import failed; an existing working scientific environment was used read-only. Git 2.22 lacks switch; checkout was used. Renderer fallback used pypdf instead of unavailable fitz. Initial comparison labels were canonicalized explicitly. One mock test initially compared /var and /private/var literally; resolving paths fixed the test. TSV line endings/trailing empty fields were normalized without changing numeric results. Historical lambda and simulation packaging errors were corrected using source evidence.

## Scientific blockers and human decisions

Original Equation 12 generator, maximum-PVE claim source, exact executed GAPIT library, original genotype/sample ordering and chip/imputation metadata remain unresolved. Authors must review the two numeric discrepancies, simulation description, interval/threshold/lambda distinctions, rights to publish aggregates/figures, citation and ethics/access metadata. These are listed in docs/author_review.md. No author approval is inferred.

## Recovery

The durable review clone is this repository. Full path inventories are here under audit/ but Git-ignored. Private source snapshots, original manuscript snapshot/render, execution logs and caches are retained at /private/tmp/livestock-reconstruction-20261104; do not publish them. Original research directories remain the source of truth if temporary storage expires.

Exact next review commands from this repository root:

```bash
git log --oneline 14ff11e..HEAD
git diff --stat 14ff11e..HEAD
python tests/check_validation_status.py
```

Continue from this state. Do not rerun completed experiments simply because a session resumes. Use docs/reproducibility.md only when a new verification is needed; use new output directories. The reviewed package push was explicitly authorized and completed. Any future publication outside that scope requires the applicable authorization.

## Default-branch results update

The earlier push published only reconstruct/20261104, leaving main at the template-only baseline. The user requested real derived_data on GitHub. All 27 populated TSV files were verified against their fingerprints; the five schema templates were relocated intact to example/schemas. The results index now links every analysis. No numerical result was altered or invented.

Publication: fast-forward remote main and the review branch to this update, without force. The default-branch update includes all 27 populated TSV files and preserves the former templates under example/schemas. The original local main remains preserved at ec6cb32. After publication, verify with:

```bash
git ls-remote --heads origin main reconstruct/20261104
git status --short --branch
```

No analysis jobs are active. Existing scientific limitations and two manuscript discrepancies remain unchanged.

Verification for this update: all 12 tests and static/schema checks passed; the numerical ledger is unchanged at 95 PASS and two documented manuscript discrepancies. No numeric analysis table changed.

## Configuration update

Replaced the active placeholder path example with config/paths.yaml and moved the original example to example/schemas. Added config/verified_results.yaml generated from executed result tables, versioned Mac/HPCC method filename patterns, and path validation. The local machine roots are installed in Git-ignored config/local_paths.yaml. Public configuration contains real repository-relative and source-relative paths, with no personal machine paths or credentials.

Read-only filesystem checks confirmed all 90 GWAS outputs on both machines. The cached Equation 13 predictions are present only in the checked Mac location; phenotype/LD/fine-mapping inputs are present on HPCC. The fine-mapping genotype files contain 20,000 animals and 135,358 SNPs. The original simulation genotype identity remains unresolved and has no active default.

Exact next configuration checks:

```bash
python scripts/check_config.py
python scripts/check_config.py --profile mac --check-inputs
python scripts/build_result_config.py --check
```

Configuration verification passed: all 15 tests, script/schema checks, result-source equality and the 90-file Mac input check. HPCC source presence was verified separately through read-only SSH. See audit/configuration_checks.txt. No analysis results were changed and no model jobs were launched. Publication targets are main and reconstruct/20261104 under the user's continuing authorization; verify both remote refs after the push.

Publication verified on 2026-11-04: configuration commit 020832eb241a52116f5f300e645d3728af1a5717 reached both remote main and reconstruct/20261104. The private profile was installed locally and remains ignored by Git. Installed verification passed for all 90 Mac GWAS files, all 94 artifact hashes and exact result-snapshot equality. The tracked working tree was clean; no jobs remain active. The commands above are the next checks for any later configuration change.

## README maintenance — 2026-11-04

Rewrote README.md in direct research prose with study design, an eight-step workflow, result locations, setup commands, configuration, repository structure, and data/citation details. Removed audit commentary from the main page while retaining links to the detailed records. The source of the published tables and use of simulated phenotypes remain explicit. All 34 local links and all Bash command blocks passed checks. No code, data, configuration, or scientific results changed; no jobs are active. Publication targets remain main and reconstruct/20261104. Next check: `git ls-remote --heads origin main reconstruct/20261104` followed by `git status --short --branch`.
