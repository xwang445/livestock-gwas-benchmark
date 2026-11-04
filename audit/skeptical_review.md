# Skeptical code and scientific review

The first section records the initial reconstruction review. Its historical discrepancy statements are superseded by the consolidation review below; the original ledger remains available.

Scope: a second pass over the completed implementation, source definitions, observed discrepancies, argument safety and rendered figures. Performed by the same reconstruction agent; this is not an independent human peer review.

| Finding | Resolution / remaining limit |
|---|---|
| Previous public simulation omitted GRM variance weight and generated only ten replicates | Restored weight 0.5 and 1,000 draws; mock-command test asserts weight/seed/count. Full simulation remains unexecuted. |
| Initial phenotype adapter assumed FID=IID | Historical converter writes FID=0; adapter corrected before any use. Original input files were untouched. |
| fastGWA relative output risk | Absolute output prefix; mock executable verifies all 12 command destinations and refusal to overwrite. |
| Historical threshold/interval differences | Tests cover strict versus inclusive, 500 kb boundaries, central-null endpoints and shared 50 kb boundaries. No unification applied. |
| Stable P-based lambda does not reproduce historical estimator | Both estimators exposed; primary finite-only beta/SE or qchisq(1-P) definition matched all 90 historical values. |
| Missing GAPIT effect orientation and stale NYC effect column | NA beta/SE/alleles; no invented cross-method effect comparison. SLEMM orientation remains explicit NA. |
| Source directories could accidentally be selected as outputs | Added a guard prohibiting summary/cache/Equation13 output beneath the historical input root; tested. Final numerical rebuild verifies the guard did not change results. |
| Misleading exploratory NCP code | Replaced public compatibility entry point with recovered downstream command; upstream generator remains missing. Prior exploratory code is preserved in Git history and originals. |
| Cached-theory agreement could be oversold | Every document and Figure 7 disclose the missing upstream Equation 12 validation. |
| Fine-map comparison could be described as LD-only | Documentation names LD, effective sample size and residual-variance estimation changes; saved MAF annotations were not recomputed from genotypes. |
| Main Manhattan cap could hide full ranges | Generated all ten untruncated supplements; detected 107.517 vs approximately 102 caption discrepancy. |
| Tests initially compared /var with /private/var literally | Corrected test to compare resolved paths; command behavior was already correct. |
| CI could appear to accept all scientific claims | Known two manuscript FAIL rows and unresolved statuses remain explicit; CI checks they are unchanged, not that all manuscript claims pass. |
| Provenance may outlive private temporary workspace | Full metadata inventories installed locally and ignored; tracked code/source hashes, aggregate references, validation ledger and run registry survive independently. Original sources remain the source of truth. |

Remaining risks: data-access conditions and software provenance require author decisions; no arbitrary expensive refits can establish historical identity after input replacement. Scientific generalization remains limited to this population/simulation architecture. FPR here measures localization outside simulated effect regions, not proof of zero LD-mediated marginal association.


## Consolidation review — 2026-11-04

Performed as a separate pass by the same agent; no independent human approval is claimed.

| Review point | Outcome |
|---|---|
| Multiple apparent workflows | One supported CLI delegates to the existing statistical implementation. Eleven old files moved outside active scripts with byte-for-byte hash checks. |
| Output safety | CLI rejects repository and historical-source output paths; existing output requires explicit summary resume. Tests cover these guards and profile dispatch. |
| Manuscript numbers | Separate raw-column scan agrees with all 90 panel counts/maxima. The corrected manuscript and current ledger use 20,461 and 107.517. Historical discrepancies remain recorded. |
| Overstated sample identity | Order and set equality are reported separately. GRM verification covers every diagonal and a deterministic off-diagonal sample; full matrix identity is not claimed. |
| Simulation versus effect-model interpretation | MPH source establishes genetic-value draws with a ridge. Abstract, Methods and Discussion wording were reconciled. Deferred maximum-PVE statements remain for author review. |
| Evidence versus rerun | Current validation does not claim new GWAS, fine-mapping or complete phenotype simulation. Cached-theory analysis does not establish Equation 12 provenance. |
| Preservation | Existing numerical result tables, figures and original DOCX are unchanged. Restricted files and the corrected manuscript are excluded from the Git update. |
| Manuscript layout | Corrected DOCX rendered to 32 pages. Every page inspected; unchanged pages compared against the prior render. No new clipping or overlap found. |
| User deferral | Equation 12 and maximum-SNP PVE search stopped. Neither question marked resolved. |

The retained GAPIT source identifies version 4.0 while the manuscript states 3.5. That historical identity and the deferred PVE/theory source remain author review items.
