# Recovered-effect PVE and lead summaries, November 11

The authors supplied a table of recovered causal effects and per-phenotype PVE,
together with scripts for fine-mapping PR curves, recovered-PVE groups, spatial
high-PIP summaries and SLEMM leads. The supplied table matches all 10,000 fixed
causal SNP IDs and covers ten phenotypes. Its recovery generator was not supplied.
These analyses therefore describe **author-supplied recovered-effect PVE**;
they do not certify the coefficients as the original independent generating draws.

The [portable calculators](../scripts/revision_2026-11-11/) were executed against
the authorized inputs. Their [public outputs](../derived_data/revision_2026-11-11/)
contain aggregates and source fingerprints. Raw recovered effects, causal IDs,
per-SNP PIPs, detailed leads, private paths and the supervisor's correspondence
remain outside this repository. The current Figure 8 uses the
[November 12 final-fit summaries](../derived_data/submission_2026-11-12/).

## Current lead results

| Method | Loci | Causal leads | Causal share | Median recovered PVE / expected PVE | Median causal percentile |
|---|---:|---:|---:|---:|---:|
| SLEMM | 16 | 11 | 68.75% | 2.36464 | 89.576 |
| GCTA MLMA | 15 | 10 | 66.67% | 2.38336 | 89.702 |
| BOLT-LMM | 16,002 | 1,519 | 9.49% | 0.66341 | 61.843 |
| REGENIE_LOCO | 15,450 | 1,432 | 9.27% | 0.61948 | 60.1935 |

These use current inclusive association decisions, updated SLEMM, and completed
GCTA outputs. Canonical leads are the minimum-P SNPs within chromosome-specific
fixed 500 kb bins; ties use position then SNP ID. They are association leads,
not a new SLEMM posterior fine-mapping analysis. The shared sheet's older LOCO
comparison values are superseded by these current counts.

The expected PVE is the design reference `h²/m = 0.5/10,000 = 0.00005`, or
0.005% of phenotypic variance. Each supplied PVE column equals the corresponding
recovered standardized coefficient squared divided by the saved phenotype
variance among **31,301 simulation animals**. The normalization check verifies
this equality; it does not rederive the missing recovery algorithm or rescale
the supplied values. The pooled median recovered-PVE/expected ratio is
0.3948388894, and the pooled mean ratio is 0.8967950609. This observed mean is
approximately 0.90 rather than exactly one. The supplied data also flag 947 causal SNPs in perfect-LD
groups. Their flags are retained in aggregate lead counts; attribution among
indistinguishable variants depends on the unavailable recovery procedure.

Percentiles rank **all 100,000 supplied causal SNP–phenotype ratios together**.
For tied values, use their average rank, divided by 100,000 and multiplied by 100.
The reported median percentile is the median of the causal leads' individual
percentiles, rather than a new ranking of the leads or the percentile of their
median PVE. Reported quartiles use pandas' linear interpolation.

## Fine-mapping estimator and supervisor baseline

Each pipeline retains one SNP per phenotype and GWAS input, taking the maximum
PIP across overlapping loci. A call is `PIP ≥ 0.9`. Causal groups are left-closed
intervals of supplied PVE divided by `h²/m`: `[0,0.5)`, `[0.5,1)`, `[1,2)`,
`[2,4)` and `[4,∞)`. The evaluated pool comprises selected fine-mapped loci;
group counts are evaluated SNP–phenotype pairs, rather than all causal SNPs in
the genome. The separate noncausal group has no recovered PVE.

The supervisor's `powbin_rep.py` estimator is reproduced before updating RSS:
compute a call fraction within each phenotype/group, take the equal-weight mean
of ten fractions, and use `mean ± 2.262157 × sample_SD / sqrt(10)` for its 95% CI.
The rounded Student's t critical value is preserved exactly. These are replicate
intervals; the pooled fraction `calls/n` is a separate statistic. No clipping of
interval endpoints is applied. Every group must contain all ten phenotypes.

The original shared Figure 8 and group table use the older RSS residual-variance
FALSE fits: 2,009 BOLT-LMM and 2,064 REGENIE_LOCO high-PIP calls. The final RSS
TRUE rerun has **2,021 and 2,072** calls. Adjusted fits remain TRUE. The public
all-run and phenotype tables label these modes explicitly. Maximum discrepancies
when reproducing the original 24-row group table are below `1e-12`; canonical
pooled call counts were independently recounted for all three runs.

The final pooled AP values are approximately 9.5624%/9.5301% for adjusted/RSS
BOLT-LMM and 9.0968%/9.0356% for adjusted/RSS REGENIE_LOCO, compared with causal
prevalence of 8.6200% and 8.4835%. AP is stepwise average precision from the full
pooled ranking, grouping equal PIPs at one threshold. It differs from trapezoidal
PR area and from the mean of ten phenotype-specific AP values. Reduced plotting
vertices are display data and do not estimate the reported exact AP.

The shared `c66.py` uses the total high-PIP call pool as the denominator for null,
central-null, and distance-qualified null shares; its `pool_*` statistics use the
total evaluated pool. These are distinct from fractions conditioned on being in
null windows. The portable recount uses the verified benchmark design encoded by
`chromosome_position_alleles` IDs: position modulo 6 Mb is null at `[2,6)` Mb and
central-null at `[3,5)` Mb. This boundary convention matches the canonical tables.

The portable `summarize_null_distance.py` also preserves `c66.py`'s null-call and
evaluated-pool shares at distance **≥200, ≥500 and ≥1,000 kb**, using the supplied
framework's `dist_to_boundary` column. The original variant-level framework is
not available locally, so those distance-qualified results have not been
published or substituted with reconstructed distances. The script requires
`--framework` explicitly. Its tests check inclusive PIP/distance boundaries and
the full-call/full-pool denominators. Original framework classification and
terminal chromosome boundary treatment remain dependencies for this analysis.

## GWAS-estimated GCTA PVE is a different quantity

The GCTA calculator computes `b_MLMA² × Var(G, ddof=1) / Var(y, ddof=1)` on the
matched **20,000 analysis animals**. Genotype variance is read from their actual
BED dosages. FAM/phenotype intersections, reference BIM, allele frequencies,
completed logs and canonical lead statistics are checked. Missing genotype calls
cause the benchmark run to stop. The BED helper's general calculation handles
mean-imputation and excludes byte-padding slots; its synthetic test checks both.

Across 15 leads, median GWAS-estimated PVE is 0.581439% (range 0.360087–0.968487%).
Among ten causal leads it is 0.525328%. These selected association coefficients
can be inflated by selection and LD. They do not replace recovered-effect PVE,
generating-effect ratios or causal-effect percentiles. Components from correlated
SNPs must not be summed.

Figure 7 still uses retained theoretical caches; supplying recovered PVE does not
recover or verify the missing upstream Equation 12/13 generators. Manuscript,
response-letter and citation-manager editing are tracked separately from this
analysis repository.
