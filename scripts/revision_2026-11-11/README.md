# Recovered-effect and association-PVE calculators

All source paths are explicit arguments. Inputs are restricted; outputs contain
only aggregate tables and source fingerprints with role labels, byte counts and
SHA-256 hashes. No detailed SNP or individual records are written. Use a new
output directory for each analysis.

```bash
python scripts/revision_2026-11-11/calculate_recovered_lead_pve.py \
  --effects "$RECOVERED_EFFECTS_TSV" --loci "$CANONICAL_DETECTED_LOCI_TSV" \
  --causal "$FIXED_CAUSAL_CSV" --reference-variances "$SIMULATION_VARIANCE_TSV" \
  --output "$NEW_LEAD_SUMMARY_DIRECTORY"

python scripts/revision_2026-11-11/check_finemapping_paragraph.py \
  --effects "$RECOVERED_EFFECTS_TSV" --causal "$FIXED_CAUSAL_CSV" \
  --adjusted "$ADJUSTED_SNP_LEVEL_TSV" --rss "$RSS_TRUE_SNP_LEVEL_TSV" \
  --rss-historical "$RSS_FALSE_SNP_LEVEL_TSV" \
  --pooled-metrics "$POOLED_FINEMAPPING_METRICS_TSV" \
  --supervisor-groups "$ORIGINAL_POWBIN_REP_TSV" \
  --output "$NEW_PVE_GROUP_DIRECTORY"

python scripts/revision_2026-11-11/calculate_gcta_mlma_pve.py \
  --genotype-prefix "$MATCHED_20K_BED_PREFIX" --reference-bim "$ORIGINAL_GCTA_BIM" \
  --causal "$FIXED_CAUSAL_CSV" --loci "$CANONICAL_DETECTED_LOCI_TSV" \
  --phenotype-dir "$FASTGWA_FORMAT_PHENOTYPE_DIRECTORY" \
  --gcta-dir "$COMPLETED_GCTA_DIRECTORY" --output "$NEW_GCTA_PVE_DIRECTORY"

python scripts/revision_2026-11-11/summarize_null_distance.py \
  --framework "$ORIGINAL_FRAMEWORK_NCP_VARIANT_TSV" \
  --adjusted "$ADJUSTED_SNP_LEVEL_TSV" --rss "$RSS_TRUE_SNP_LEVEL_TSV" \
  --rss-historical "$RSS_FALSE_SNP_LEVEL_TSV" --output "$NEW_DISTANCE_SUMMARY_DIRECTORY"
```

`--supervisor-groups` is optional; when present, it checks exact reproduction of
the original adjusted TRUE/RSS FALSE grouping before recalculating final RSS TRUE.
The fine-mapping inputs contain `phenotype,GWAS_method,SNP,PIP` columns. Recovered
effects contain `snp`, `alpha1..alpha10`, `pve1..pve10` and `in_perfectLD_group`;
causal truth is CSV with a unique `SNP` column. Canonical loci use the existing
association derivation schema. Reference variance uses `PHENOTYPE,VAR_Y`.

The standalone null-distance calculator ports the supplied `c66.py` estimator.
It requires the original framework columns `snp,is_effect_window,is_null_window,
is_central_null,dist_to_boundary`. It preserves that framework's boundaries and
uses inclusive distance cutoffs. That source framework is currently unavailable;
no distance-qualified table is distributed and no reconstructed distance input
is substituted. Overall and central-null counts are already independently
verified in the other fine-mapping calculator.

GCTA phenotype files are named `pheno.1.txt` through `pheno.10.txt` with FID, IID
and the corresponding named phenotype column. Outputs/logs use
`spg.even.100k.phenoN.mlma/.log`. The calculator enforces the benchmark's 20,000
animals, 135,358 SNPs, ten phenotypes, 10,000 causal SNPs and 15 GCTA leads.

The [scientific record](../../docs/recovered_effects_2026-11-11.md) defines the
estimators, source limitations and distinction between recovered-effect and
GWAS-estimated PVE. No fitting or phenotype simulation is launched by these scripts.
