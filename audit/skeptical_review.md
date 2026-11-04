# Separate skeptical code and scientific review

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
