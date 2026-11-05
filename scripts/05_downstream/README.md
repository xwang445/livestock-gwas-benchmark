# Current downstream analysis

Use the [supported workflow](../../docs/workflow.md). Default `workflow.py summarize` calls `revision_2026-11-05/association/derive_association_tables.py` and its power-grid derivation with the verified simulation, genotype and GCTA locations. It writes one canonical set of tables under the new private output directory.

All current SNP and locus decisions use P ≤ 5 × 10⁻⁸; REGENIE native LOG10P is compared directly. Updated SLEMM and the ten-configuration manifest supply counts, window distributions, spillover, lambda, overlap and plots. Incomplete GCTA coverage is explicit. Per-variant detailed outputs remain private.

`benchmark_core.py summarize --historical` is a gated historical audit, including the old strict cutoff and original SLEMM source. It is not an alternate current estimator. Existing strict outputs and archived R scripts remain available for scientific comparison. Current Equation 13 uses the updated SLEMM file but still requires the unavailable verified prediction cache.
