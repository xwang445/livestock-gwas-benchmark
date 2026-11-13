# Statistical definitions

## Current revision

All current association counts, Figure 3/Table S1 loci, Figure 5 overlaps, spillover, power and Manhattan selections use **P ≤ 5 × 10⁻⁸**. For native −log10(P) columns, compare directly with −log10(5 × 10⁻⁸); do not exponentiate before the decision or widen the cutoff with a numerical tolerance. `significance.py` is the shared policy. Equality-boundary regressions cover both source scales. Current default commands use updated SLEMM and the ten-configuration registry.

The original SuSiE locus manifest already used the inclusive threshold and remains unchanged; no new fine-mapping loci or fits are generated solely to unify this rule. Central-null fractions retain their distinct denominators, and current grid/boundary definitions are recorded in the revision captions and inventories.

The [threshold audit](threshold_consistency.md) records the raw-source checks, including exact-threshold calls and preservation of the fine-mapping manifest.

The root lambda tables retain the Figure 6 statistic estimator for current method fits, with an explicit `statistic_definition` and a stable P-based alternative. All ten configurations are included. Current root SNP/locus counts use the inclusive rule; strict comparison columns are labeled separately.

## Explicit historical audit definitions

The definitions below document immutable pre-revision results in `derived_data/historical_2026-11-04/`. Their strict cutoff applies only when `--historical` is explicitly selected and does not supply current primary results.

- Figure 3 significant SNPs: strict P<5e-8. Loci are chromosome-specific floor(bp/500000) bins; no distance-based clumping or lead-SNP selection. Inclusive sensitivity counts are also written.
- Figure 5 uses those strict loci. Jaccard is intersection/union, with NA whenever either set is empty, as specified by the October caption. Cross-replicate stability is a separate table. Empty comparisons must not be silently converted to zero.
- Spillover: inclusive P<=5e-8. Only full null gaps between consecutive effect intervals are used. The historical integer central interval is [null_start+1000001,null_end-1000000], with four closed 500 kb windows. Denominator 1,480 includes empty windows. This differs at boundaries from the simulation's half-open rule.
- Historical lambdaGC: median finite chi-square / qchisq(.5,1). For mixed models the source script uses (beta/SE)^2; GAPIT uses qchisq(1-P,1) and drops nonfinite statistics. That inverse-CDF calculation loses extreme tails to floating-point rounding. The primary historical result and a stable P-based inverse-survival diagnostic are both retained and labeled. Neither is applied as a correction to original GWAS P values.
- Power is the fraction of the 10,000 causal positions with P<=alpha among nonmissing tested causal SNPs. FPR is the corresponding fraction in central-null or outside-effect regions, reported separately. It is a localization metric under this architecture, not an assertion that LD-induced marginal associations have no signal.
- Figure 6C–D: full 4 Mb null gaps split into 80 bins, mirrored distance indices 1..40. The original closed-interval join can count an exact boundary SNP in both adjacent 50 kb bins. Preserve that convention. Compute significant/(tested*alpha) per nonempty window, then average across windows and phenotypes. Pooled ratio is a separate diagnostic.
- Equation 13 observed statistic: BOLT (beta/SE)^2; REGENIE and SLEMM native CHISQ. Theory-only 20 equal-count bins retain input order to resolve ties; the central-null subset is independently binned with chromosome-position key tie order. Intervals use t(9,.975)*SD/sqrt(10) across replicate bin means. This is downstream reproduction from cached predictions; it does not validate the missing upstream Equation 12 calculation.
- Fine-mapping uses a different locus definition: +/-1Mb around GWAS-selected leads, overlaps split at lead midpoints. Count unique high-PIP SNPs per method/pipeline/phenotype, PIP>=.9, from converged fits. Null is the entire 4 Mb null region. Percentages are ratios of pooled counts, not averages of phenotype percentages. SuSiE-adj and RSS differ in LD, sample-size input and residual-variance estimation.

Harmonization does not infer GAPIT beta/SE or effect alleles: local refinement may replace P without replacing effect. SLEMM source REF/ALT are retained, but harmonized effect-allele orientation remains NA pending independent confirmation. P values are never recalibrated. Zero P remains zero; missingness, duplicates, invalid coordinates and invalid ranges are checked.
