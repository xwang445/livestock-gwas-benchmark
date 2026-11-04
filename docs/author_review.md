# Decisions requiring human author review

1. Correct the Results REGENIE_LOCO significant-SNP upper bound to 20,461 if the selected Figure 3 remains authoritative; 21,484 is BOLT-LMM.
2. Correct or explain the Figure 4 approximately 102 maximum: independently read raw REGENIE_LOCO phenotype 8 LOG10P is 107.517.
3. Recover the original Equation 12 generator and verify its covariance projection, variance scaling, genotype standardization, allele/sample order and causal LD score calculations. Cached downstream agreement is insufficient.
4. Identify the exact source for maximum SNP PVE .063–.093%, corresponding NCP and marginal-power range. The retained LSMR estimate uses a different sample/estimator and gives .209–.433%; it must not replace the original silently.
5. Reconcile GAPIT 3.5 manuscript wording with retained source self-identifying 4.0. Historical logs prove both models ran but not the exact library bytes.
6. Describe the executed MPH genetic-value Cholesky simulation accurately, including 1,000 generated replicates, first-ten/first-20,000 selection, GRM weight, diagonal ridge and lack of original explicit SNP-effect draws.
7. Confirm historical genotype identities/order and the correct source prefix. A prior audit found the historical path's labels/order had changed. Confirm array model, QC chip list, imputation production software/version and INFO interpretation.
8. Clarify phase-specific strict/inclusive significance and interval-boundary conventions. Decide whether to retain historical finite-only GAPIT lambda or add the stable diagnostic as sensitivity; a change would be a documented analysis revision, not a reconstruction correction.
9. Confirm rights to publish aggregate tables and detailed Manhattan figures under the governing data agreement, and confirm author/ethics/access/citation metadata. No such approval is recorded here.
10. Identify the actual scheduler/allocation before any expensive new fitting. Full pipeline refitting remains NOT_RUN. Review resource concurrency before submitting the prepared template.

No preferred model was used as an acceptance target. Counts were not altered to force agreement, and unresolved claims remain marked as such. No manuscript submission, preprint deposit, email, or GitHub push was performed.
