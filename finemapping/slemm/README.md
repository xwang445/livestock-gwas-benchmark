# SLEMM lead-SNP summaries

The requested workstation `finemapping/slemm/` source directory is not accessible in this execution environment. This directory documents the reproducible replacement for lead-SNP causal status and effect summaries; it does not claim the missing original pipeline was recovered.

The canonical association derivation under `scripts/revision_2026-11-05/` uses the verified SLEMM qf90/window3000 outputs for all ten phenotypes. Lead SNPs are the minimum-P SNP within each detected chromosome-specific 500 kb bin, with deterministic tie ordering. Causal membership uses the fixed verified 10,000-SNP truth. Summary tables include causal lead counts and absolute-effect distributions, with source allele/orientation limitations recorded separately. These are lead-SNP association summaries, rather than a new SLEMM posterior fine-mapping model.

The private detailed lead table contains variant identifiers and is excluded from the public package. The aggregate canonical association tables and old/new SLEMM reconciliation table are the publishable outputs. The revision record gives their exact paths and commands. Effect-size-stratified fine-mapping detection based on a verified historical SNP-effect definition remains dependent on the missing PI-supplied definition; reconstructed inverse-regression effects must not be labeled original simulated effects.
