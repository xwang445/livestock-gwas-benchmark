# Derived aggregate results

These are real regenerated downstream summaries, not simulated demonstration numbers. They derive from the saved historical GWAS, cached theoretical predictions, saved converged fine-mapping fits and pairwise LD. The synthetic generator is separate in example/synthetic.

`input_manifest.tsv` fingerprints all 90 GWAS files and their coordinate panels. `run_manifest.json` records the summary implementation hash. `equation13_summary/cached_theory_validation.tsv` records the cached-theory hash and algebra errors, explicitly marking upstream Equation 12 NOT_RUN. Fine-mapping verification_manifest.json fingerprints the saved PIPs/diagnostics/MAF inputs and records limits. `method_summary/heritability_from_logs.tsv` contains estimates recovered from fitting logs, not newly refitted values.

Subdirectories contain per-phenotype and aggregate method counts, power/FPR, spatial enrichment, spillover, lambda, overlap/stability, Equation 13 bins, fine-mapping counts and LD means/fitted curves. Historical comparison aggregates are separately stored in audit/reference_values with source fingerprints. Templates ending `_template.tsv` remain schema examples only.

All files are classified as non-identifying aggregate review candidates. Authors must confirm publication rights before any push; raw inputs and per-variant outputs remain excluded. See docs/data_availability.md and repository_validation_report.md.
