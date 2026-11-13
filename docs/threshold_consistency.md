# Consistent association threshold

The current revision uses **P ≤ 5 × 10⁻⁸** for association significance throughout SNP counts, 500 kb loci, overlap/Jaccard sets, regional/spillover summaries, Manhattan selection and annotations, and fine-mapping locus selection. For native REGENIE output the equivalent decision is **LOG10P ≥ −log10(5 × 10⁻⁸)**, evaluated before exponentiation. No epsilon or approximate-equality tolerance is added.

The shared Python implementation is `scripts/revision_2026-11-05/association/significance.py`. Current configuration and the default workflow select the inclusive revision. Strict P < 5 × 10⁻⁸ survives only in explicitly labeled historical reproduction and audit columns; the old result files remain preserved. This prevents old strict Figure 3/5 output from silently becoming the current result.

The initial November 5 independent audit checked 101 source outputs, 91 current Manhattan panels, 829 pairwise overlap records, 10302 sensitivity-rate rows, and 36 Table S1 statistics, with no primary boundary disagreements. BOLT-LMM has 212 SNP records exactly on the cutoff, contributing 19 additional phenotype–locus detections relative to the historical strict rule. Its inclusive maximum is 21500 significant SNPs in phenotype 3, compared with 21484 using the strict cutoff. The revised canonical primary tables already used the inclusive rule, so their primary numbers remain unchanged by this enforcement.

The cutoff-sensitivity tables retain all 10302 numerator/denominator pairs. Exact cutoff coordinates, 17-digit decimal serialization and hexadecimal float representations now prevent a reread from changing the represented cutoff. Classification uses native source statistics; the plotted P coordinate is a display value.

All 13,301 fine-mapping loci reconstruct exactly under the inclusive selection rule. Existing same-locus SuSiE reruns remain valid; no refitting is needed for this correction. PIP ≥ 0.9 is a separate posterior-call criterion. LambdaGC uses the defined full central-null statistic pool rather than selecting significant SNPs.

Figures were regenerated from canonical tables. The pooled Manhattan plots now leave space below the threshold line so boundary points are visible. All ten GCTA phenotypes subsequently completed and validated; the current manifest contains 100 completed panels and zero pending panels. The [November 12 submission exports](revision_2026-11-12.md) retain the inclusive rule.

Evidence: `derived_data/revision_2026-11-05/association/inclusive_threshold_audit.json`, `inclusive_threshold_by_input.tsv`, and `power_precision_fix_validation.json`.
