# LD verification for the November 5 revision

The pairwise-LD input and PLINK log were read independently of the figure caption. PLINK v1.9.0-b.8 used 31,301 animals and 135,358 markers across 18 autosomes, with `--r2 gz --ld-window 999999 --ld-window-kb 10000 --ld-window-r2 0`. These LD summaries describe the full genotype panel; the GWAS subset contains 20,000 animals.

The raw file contains 78,082,232 eligible pairs. Independent streaming recalculation reproduced all 1,000 bin counts exactly; the maximum difference in bin means was 5.56 × 10⁻¹⁶. The [machine-readable verification](../audit/ld_verification_2026-11-05.json) records the source checksum. Only non-identifying aggregate LD summaries are distributed.

Both panels use one table of 10 kb bins across 0–10 Mb. Bins are left-closed/right-open, with the final 10 Mb endpoint included. Exactly 2 Mb belongs to [2.00, 2.01) Mb rather than being reassigned to the preceding bin in the shorter view. The shorter panel changes only the coordinate limits. Every eligible SNP pair receives equal weight within a bin. The mean curve is a natural cubic regression spline in log10(distance), with 20 basis columns, equal weight across bins, and a constraint at the observed first-bin mean. The plotted coordinate arrays are tested for identity across panels.

The first observed bin, [0, 10 kb), contains 31,209 pairs and mean R² 0.592525238699, plotted at 0.005 Mb. There is no observation at distance zero. Mean R² in [0.99, 1.00), [1.99, 2.00), and [9.99, 10.00] Mb is approximately 0.2330, 0.1609, and 0.0396. These are observed bin means, rather than fitted values.

The original pre-harmonization scripts differed: 10 kb bins in the 2 Mb panel versus 100 kb bins in the 10 Mb panel, with independent LOESS fits. Their first-bin means were 0.592525 and 0.485459; their first smoothed values were approximately 0.557520 and 0.464157. A prior harmonization corrected this discrepancy. The November 5 regeneration verifies the harmonized source and uses lowercase panel labels and separate vector PDF, 600 dpi PNG and TIFF files.

Reproduce from the repository root:

```bash
python scripts/02_ld_analysis/verify_ld_means.py --input "$LD_PAIR_FILE" > ld_recount.json
Rscript scripts/13_figures/plot_ld.R derived_data NEW_FIGURE_DIRECTORY
```

`LD_PAIR_FILE` is the authorized restricted `ld_10Mb.ld.gz` input. Compare the recomputed bin counts/means with `derived_data/ld_summary/ld_decay_10kb.tsv`. The raw file and animal IDs remain excluded from Git.
