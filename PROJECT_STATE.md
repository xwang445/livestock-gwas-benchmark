# Project state

Updated 2026-11-04. Workflow consolidation, the authorized provenance checks, manuscript correction, installation and GitHub publication are complete. Both remote main and reconstruct/20261104 were verified at content commit 997d219c5fb10e876e235d140834c30758599f8e; this state-only follow-up records that verification. No jobs are active. The installed repository is clean, its private profile remains ignored, and all 90 Mac GWAS input paths passed the installed configuration check. The previous state and publication history are preserved in Git at commit 0182c6e.

## Completed

- One supported entry point: `python scripts/workflow.py`. The root README, result index, workflow guide, validation report, configuration notes and stage READMEs now describe the same commands.
- Eleven superseded files moved intact to `archive/legacy_workflow`; their original hashes and replacements are in `audit/legacy_relocation.tsv`.
- MPH v0.55.1 source checked: genetic-value covariance 0.5G + 0.00005I, residual variance 0.5, 1,000 generated replicates, default seed mapping and no phenotype rescaling.
- Retained simulation/chip/GRM sample order verified. All ten phenotype subsets exactly match the original simulation values. GWAS and fine-mapping inputs contain the correct ID sets in different row orders; alignment is by IID.
- All 31,301 GRM diagonals and 4,096 deterministic off-diagonals checked against the retained causal genotypes. This was a partial matrix check, not a full original-file comparison.
- Independent raw-column scan of all 90 GWAS panels confirmed REGENIE_LOCO's significant-SNP maximum of 20,461 and the overall untruncated maximum −log10(P) of 107.517 in REGENIE_LOCO phenotype 8. The source of the historical approximately 102 wording was not established.
- A separate corrected local DOCX contains the two numerical corrections and consistent MPH/sample-selection wording. The original manuscript remains unchanged. Installed local delivery location: `../manuscript_review_20261104/`, containing the corrected DOCX and full edit log; manuscript files are excluded from Git.

## Evidence and validation

The current numerical ledger contains 108 PASS, 3 NOT_RUN, 2 AMBIGUOUS and 1 REQUIRES_RECOMPUTATION records. The original 95-PASS/two-FAIL manuscript ledger is preserved separately. Nineteen unit, synthetic and command-routing tests pass. Static/schema checks, result-configuration equality and corrected-claim checks pass through `python scripts/workflow.py check`.

The original 27 numerical result tables and all reconstructed figures are byte-identical to commit 0182c6e. This update adds provenance evidence and corrects descriptions; it does not replace prior numerical analyses. See `audit/run_registry.tsv`, `audit/consolidation_checks.txt`, `audit/manuscript_corrections.json`, and `docs/scientific_provenance.md`.

## Deferred and remaining author review

Equation 12/theory and maximum-SNP PVE provenance were deferred at the user's request because those inputs came from Professor Jicai Jiang. Their manuscript statements remain unchanged. The broad source search was interrupted after this instruction; it produced no final evidence and was not used to support a conclusion.

The exact historical GAPIT version, chip/imputation production metadata, data-access terms and final manuscript approval remain author items. PVE wording still needs reconciliation with the implemented MPH sampling procedure when the deferred source is reviewed. New GWAS, SuSiE fitting and the complete simulation were not rerun. No manuscript submission, preprint deposit or external messages were made.

## Recovery

Active jobs: NONE after final checks. Failed/aborted jobs: broad provenance search interrupted on user deferral (exit 255); no inputs or results modified. The isolated working checkout was `/private/tmp/livestock-consolidation-20261104`; private source inspection and manuscript rendering used `/private/tmp/livestock-provenance-20261104`. Original research data and the installed repository remain durable sources if temporary files expire.

The installed branch is `reconstruct/20261104`; original local `main` remains at ec6cb32. Publication targets are remote `main` and `reconstruct/20261104`, under the user's existing push authorization. Do not force-push or stage the ignored private path profile.

Next verification commands from the repository root:

```bash
python scripts/workflow.py check
git ls-remote --heads origin main reconstruct/20261104
git status --short --branch
```

Continue from verified state. Do not restart completed analyses on resumption. Deferred provenance requires the author's source material and renewed instruction, not another unconstrained search.
