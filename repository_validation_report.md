# Repository validation report

The repository now has one supported execution guide, [docs/workflow.md](docs/workflow.md), and one entry point, `scripts/workflow.py`. Superseded downstream scripts are preserved outside the active code tree in `archive/legacy_workflow`. File hashes in [legacy_relocation.tsv](audit/legacy_relocation.tsv) show that archival moves preserved their contents.

## Result validation

The current ledger contains 108 PASS, 3 NOT_RUN, 2 AMBIGUOUS and 1 REQUIRES_RECOMPUTATION records, with no failed corrected numerical claims. These statuses describe the specific comparisons, not an end-to-end refit.

| Analysis | Executed evidence |
|---|---|
| GWAS downstream summaries | 90 panels, each containing 135,358 SNPs; counts, power/FPR, spillover, lambda and overlap checked against retained references |
| Equation 13 | 1,600 replicate-bin rows and 160 aggregate rows checked against cached theoretical inputs |
| Fine-mapping | 40 method/pipeline/phenotype groups independently recounted from saved converged fits |
| LD | 1,000 distance-bin means and counts independently recounted |
| Figures | 8 main figures, S1–S10, and central-null Equation 13 supplement |
| Manuscript corrections | Independent raw-column scan reproduced all 90 counts and maxima |
| Sample provenance | Exact phenotype chain and retained chip/GRM ID order; all GRM diagonals and 4,096 off-diagonal entries checked |

The [current ledger](audit/manuscript_validation.tsv) records values, tolerances and sources. The [original October-manuscript ledger](audit/historical_manuscript_validation.tsv) preserves 95 PASS and two FAIL records before correction; those original discrepancies were not erased.

## Corrected manuscript values

REGENIE_LOCO's significant-SNP maximum is **20,461**, not 21,484. The latter is the BOLT-LMM maximum. The complete range for REGENIE_LOCO is 12,738–20,461.

The Figure 4 uncapped maximum is **107.517**, from REGENIE_LOCO phenotype 8. It is present directly in the raw LOG10P column. None of the 90 panel maxima falls between 100 and 103. The plot's display cap does not account for the old approximately 102 statement. The corrected caption states the exact value and panel; the original approximation's source was not established.

The original DOCX was preserved. A separate local copy contains these two corrections and the verified MPH/sample-order wording. The manuscript itself is excluded from the public repository.

## Scientific provenance

MPH v0.55.1 generated Gaussian genetic values from 0.5G + 0.00005I and residuals with variance 0.5, with no explicit saved SNP-effect draws or post-simulation rescaling. It generated 1,000 replicates; the benchmark uses the first ten and first 20,000 saved simulation rows. Versioned source hashes, phenotype comparisons and GRM checks are in [scientific_provenance.md](docs/scientific_provenance.md).

The retained chip file shares the saved simulation's 31,301-ID order. GWAS and fine-mapping inputs have the correct ID sets in different row orders and are aligned by IID. The old simulation source path now contains different labels/order and must not be reused as if unchanged.

Equation 12 and maximum-SNP PVE provenance were deferred at the author's explicit request because those inputs came from Professor Jicai Jiang. No replacement calculation was substituted. The retained GAPIT source-version discrepancy remains in the author review record. New GWAS, SuSiE fitting and simulation were not run during this update.

## Checks and recovery

Run `python scripts/workflow.py check`. This exercises static syntax, schemas, configuration, synthetic/unit tests, source-table equality and corrected-claim status. Restricted files are not needed for public checks.

[audit/run_registry.tsv](audit/run_registry.tsv) records executed jobs. [PROJECT_STATE.md](PROJECT_STATE.md) records completion, output locations and the next command. Restricted inputs and original results were not modified. Only code, documentation and non-identifying aggregate evidence are published.
