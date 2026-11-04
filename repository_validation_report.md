# Repository reconstruction and validation report

Review package prepared on 2026-11-04 from the author-confirmed October manuscript. Original Mac/HPCC files and the manuscript were preserved. **The downstream reconstruction is complete; end-to-end scientific reproducibility remains partial. Following explicit user authorization on 2026-11-04, the reviewed commits were pushed to GitHub branch `reconstruct/20261104`.**

Publication note: the remaining report records the pre-push audit and review findings. Statements that remote CI had not run or that publication authorization was pending describe that audit snapshot; scientific limitations remain unresolved.

The numerical ledger contains 95 PASS checks, two FAIL, four AMBIGUOUS, three NOT_RUN and one REQUIRES_RECOMPUTATION. PASS means numeric agreement with a stated manuscript value or a historical manuscript-supporting aggregate, at the recorded tolerance. It does not mean a new GWAS fit succeeded. Every PASS gives expected, observed, difference, tolerance and source in [audit/manuscript_validation.tsv](audit/manuscript_validation.tsv).

| Stage | Status | Source | Validated against manuscript | GitHub ready |
|---|---|---|---|---|
| Genotype QC | PARTIAL | HPCC QC script and prior chip/INFO audits | Counts supported by historical evidence; no new QC fit | Code ready; source identity unresolved |
| LD decay | COMPLETE | HPCC pairwise LD and harmonized binning code | 1,000 counts/means independently matched | Yes, aggregate review candidate |
| Architecture / simulation | AMBIGUOUS | Historical R/MPH commands, versioned MPH source | Design recovered; new simulation not run | Code ready with implementation caveats |
| Nine GWAS configurations | PARTIAL | HPCC runners/logs and 90 local results | All 90 panels validated downstream; no new fits | Code ready; GAPIT version unresolved |
| Harmonization | COMPLETE | Nine observed output schemas | 90 × 135,358 rows and common coordinates checked | Yes; per-variant exports remain restricted |
| Significant SNPs/loci | COMPLETE | Local historical Figure 3 logic | All 90 SNP and 90 locus cells matched | Yes |
| Power / FPR | COMPLETE | Local pwer_analysis.R and 90 outputs | 2,970 historical rows for each curve matched | Yes |
| Spillover | COMPLETE | Local Figure 6A code | 90 rows; nine rounded rates; denominator 1,480 | Yes |
| Lambda GC | COMPLETE | Local historical statistic definitions | All 90 historical values matched | Yes; numerical tail limitation disclosed |
| Overlap / stability | COMPLETE | Local Jaccard code, manuscript caption | 30 displayed spot checks; full tables recomputed | Yes; stability is supporting output |
| Equation 12 / upstream theory | REQUIRES_RECOMPUTATION | Cached prediction table; generator not found | NOT_RUN upstream | No end-to-end claim |
| Equation 13 downstream | COMPLETE | Cached table and local production R | 1,600 replicate bins, 160 means/intervals matched | Yes with upstream limitation |
| Fine-mapping fitting | PARTIAL | HPCC primary LOCO/RSS code and saved fits | Fit counts recounted; fitting not repeated | Code ready with restricted inputs |
| Fine-mapping evaluation | COMPLETE | Saved converged PIPs and MAF annotations | All four counts, precision and MAF expectations match | Yes; MAF not recalculated from genotypes |
| Figures / supplements | COMPLETE | Verified aggregates plus restricted GWAS for Manhattan | 8 main + S1–S10 + central-null supplement rendered | Yes after data-agreement review |
| Synthetic / tests / CI | COMPLETE | New artificial fixtures and tests | Local tests/static checks passed; remote CI not run | Yes |
| Maximum SNP PVE | AMBIGUOUS | Recovered LSMR outputs differ from manuscript | Manuscript range not reproduced | No validated claim |
| Metadata / publication | PARTIAL | Author review required | Data rights, array/imputation and citation unresolved | No push authorized |

## 1. Local files discovered

The original local W2Mb_N20k inventory contains 590 files (about 3.29 GB), including all 90 GWAS output panels, power/FPR/spillover/lambda tables, the 135,358-row cached theory table and historical figures. Additional manuscript-folder inventory contains 87 files, including production plotting code. The actual local Git clone was found clean at ec6cb32, behind public 14ff11e. Full inventories remain local and Git-ignored.

## 2. HPCC files discovered

Read-only inventory found 502,804 files (about 173 GB), 190 small code/document candidates and no traversal errors. Sources include historical GWAS, simulation/QC evidence, LD, fine-mapping branches, fitted results and prior provenance audits. Only selected small code/log/aggregate evidence was returned; no large genotype/GRM/LD dataset was copied.

## 3. Authoritative scripts identified

Independent numeric comparisons support the selected local Figure 3, spillover, power/FPR, lambda and Equation 13 production definitions and HPCC harmonized LD binning. Kansas/NYC logic is supported by retained library branches and output schemas, but exact executed source version is unresolved. All 208 code candidates are classified in audit/script_provenance.tsv; current classification counts are {'AUTHORITATIVE': 7, 'UNKNOWN': 115, 'EXPLORATORY': 24, 'SUPERSEDED': 13, 'LIKELY_AUTHORITATIVE': 47, 'FAILED': 2}.

## 4. Superseded scripts identified

Thirteen candidates were classified SUPERSEDED, including pre-harmonization LD code, older fine-mapping branches and older overlap conventions. Exploratory NCP/PVE work is distinguished from production evidence. No historical file was deleted or changed. The misleading public calculate_ncp.R entry now calls only the recovered downstream implementation; its old content remains in Git history.

## 5. GitHub files added

New public-review code covers genotype preparation, causal selection, phenotype formatting, streaming LD verification, common GWAS harmonization/statistics, independent fine-mapping recount, aggregate exports, scientific validation, figure rendering, synthetic fixtures, automated tests and CI. New audit/reference/derived tables and 19 PDF/PNG figure pairs provide evidence rather than placeholders. Exact added paths and hashes appear in audit/change_manifest.tsv and audit/artifact_manifest.tsv.

## 6. GitHub files modified

Corrected simulation generation count/GRM weight/seed; safe and configurable five GWAS shell runners; GAPIT wrapper invoking both models with an explicit retained-source hash; Equation 13 compatibility entry; README, stage documentation, data/rights limits and current manuscript map. Original local main remains preserved; changes are on reconstruct/20261104.

## 7. Files deliberately excluded

Genotypes, individual phenotypes/IDs, raw GWAS, causal SNP lists, GRMs, raw LD, detailed SNP PIPs, full manuscript DOCX/render, private caches, full path inventories and third-party source snapshots are excluded. Synthetic fixtures are artificial and generated on demand. Large files were not staged automatically.

## 8. Analyses recovered from local Mac

All nine output schemas and 90 GWAS files, strict SNP/locus counts, historical central-null lambda, power/FPR, 50 kb spatial enrichment, spillover, pairwise/cross-replicate overlap, cached theory and Equation 13 plotting logic, and main/supplementary Manhattan definitions.

## 9. Analyses recovered from HPCC

QC and simulation command provenance; replaced-genotype-path evidence; full-resource LD and its binning/spline code; variance-component logs; primary LOCO SuSiE-adj/RSS pipelines, converged PIPs and MAF pools; earlier audit/reference tables used for comparison. Independent verification read saved LD and fitted results, without new fitting.

## 10. Analyses reconstructed from manuscript

Manuscript-defined stage map, empty-set Jaccard convention, rounded spillover/fine-map claims, expected chi-square = 1+NCP, confidence intervals across replicate bin means, and the illustrative effect/null diagram. The idealized approximately 1.2 causal detections calculation is an analytic normal-mixture check and is explicitly separate from the unresolved maximum-PVE estimate. Missing Equation 12 algebra was not invented.

## 11. Analyses requiring recomputation

Upstream Equation 12 requires recovery of the exact original generator or an author-reviewed implementation from the verified source equation and original genotype ordering. The exact command cannot be specified honestly until those are resolved. Maximum-PVE provenance also needs recovery before an appropriate rerun can be chosen. Full genotype QC, phenotype simulation, GWAS and SuSiE fitting are NOT_RUN in this reconstruction; commands are prepared but resource/scheduler and input identity require confirmation.

## 12. Manuscript results successfully reproduced

All 90 significant-SNP and 90 significant-locus cells; all 90 historical lambda values; nine rounded spillover means; 2,970 historical power rows and 2,970 central-null FPR rows; 3,600 spatial-enrichment rows; 30 displayed Jaccard spot checks; 1,600 Equation 13 replicate-bin means and 160 summary rows/intervals; 1,000 LD counts/means; fastGWA/SLEMM h2 summary rounding; 6,732/6,569 converged fine-map locus counts in each pipeline; four high-PIP/null/causal counts and rounded MAF expectations/enrichments. Figures 1–8 and S1–S10 were redrawn, plus a central-null Equation 13 supplement.

## 13. Results that do not match

- Results text assigns REGENIE_LOCO maximum 21,484 significant SNPs; raw counts and Figure 3 show 20,461. The 21,484 value belongs to BOLT-LMM.
- Figure 4 says maximum −log10(P) approximately 102; the raw maximum is 107.517 in REGENIE_LOCO phenotype 8.
- Recovered PVE output gives .209–.433% on 31,301 animals, versus .063–.093% in the manuscript. This is AMBIGUOUS because the original estimator/input combination has not been identified.
- Retained GAPIT source identifies 4.0 versus manuscript 3.5, without evidence proving the exact historical executed version.

The complete ledger is audit/discrepancies.tsv. No definition was changed merely to force agreement.

## 14. Privacy and restriction concerns

All tracked derived tables contain aggregates, with no animal IDs, SNP identifiers, per-variant effects/PIPs or raw matrix entries. Detailed Manhattan plots and aggregate tables still require the authors' data-agreement review before public release. Existing license/citation metadata were preserved; no new rights or approvals are inferred.

## 15. Remaining TODO items

See docs/author_review.md: two text corrections, Equation 12 generator, maximum-PVE evidence, GAPIT version, original genotype/sample order and array/imputation metadata, exact simulation description, threshold/boundary/lambda disclosure, access/license/citation/ethics decisions, and execution allocation for any future expensive rerun. These unresolved items are not marked complete.

## 16. Git diff summary

See audit/change_manifest.tsv for exact added/modified file inventory and audit/git_review.txt for the local diff/commit summary. Original source directories were unchanged and unrelated existing work was preserved. Private inventories are installed locally but ignored. No remote update was performed.

## 17. Local commits

Changes are split into logical audit, implementation, derived-results/figures, tests/CI and documentation commits. The initial audit commit is 32ffb12. Implementation/results/test commit IDs are recorded in audit/git_review.txt and PROJECT_STATE.md; the final documentation commit is identified by branch HEAD. Review these commits before any later authorized push.

## 18. Recommended next steps

Review the scientific discrepancy ledger and metadata decisions, confirm publication rights for the candidate aggregates/figures, recover missing Equation 12/PVE/version evidence, and inspect the local Git diff. Only after separate explicit approval should the review branch be pushed or manuscript claims revised. The current repository is usable for public downstream validation and synthetic tests, with its limits visible.
