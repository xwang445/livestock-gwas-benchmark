# Initial audit workflow (historical planning record)

1. Recover original genotype QC and causal-selection commands from the upstream paths referenced by HPCC logs. Pin input hashes and sample order. Do not rerun an altered genotype source.
2. Preserve GRM-based MPH generation with 1,000 draws and the historical first-ten / first-20,000 subset; verify the software RNG before claiming exact regeneration.
3. Retain all nine historical GWAS outputs, using the recovered GAPIT source to define Kansas and NYC. Repair wrappers independently of expensive GWAS execution.
4. Harmonize with a versioned schema, preserve original P/statistic values and explicitly label missing alleles. Compute counts, spillover, lambda, power/FPR and overlaps into public aggregates. Stage-specific interval and threshold conventions remain explicit.
5. Recover Equation13 downstream binning from the original variant prediction file. Verify its algebra; keep Equation12 provenance unresolved until the generator is found or verified from the source paper.
6. Independently recount saved fine-mapping calls on the server; export only aggregate counts, not calls, fitted objects, LD or GRMs. Keep primary LOCO fits separate from exploratory whole-genome/sensitivity runs.
7. Compare independently regenerated aggregates to the manuscript and saved outputs, make figures from aggregates, and add synthetic/static tests and CI. Disagreements remain in the discrepancy ledger rather than being tuned away.

The supplied server has no scheduler executable on PATH. No expensive jobs are warranted yet: GWAS and fine-mapping outputs already exist. Full fit regeneration is distinct from aggregate validation and must not be implied by a successful count comparison.

The implementation and verification are now recorded in docs/workflow.md; this file preserves the initial audit decisions, not current execution status.
