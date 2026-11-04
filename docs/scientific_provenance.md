# Simulation and sample provenance

## MPH implementation

The historical command used MPH v0.55.1 and one relationship matrix built from 10,000 fixed causal chip SNPs. With standardized genotype matrix Z and m = 10,000, G = ZZᵀ/m. The GRM weight was 0.5. MPH generated genetic values by Cholesky factorization of 0.5G + 0.00005I and added independent residuals with variance 0.5. It did not explicitly draw or save SNP effects and did not rescale the simulated phenotypes.

The default seed argument 0 initializes mt19937 with 54321. The command generated 1,000 genetic random vectors followed by 1,000 residual vectors. The benchmark selected phenotypes 1–10 and the first 20,000 simulation rows. Changing the generation count changes the residual random stream even for the first ten phenotypes.

These statements follow the [versioned MPH source](https://github.com/jiang18/mph/tree/13ffe633e0170d35e01e6164fc2bb4e01f35e2f2), especially `src/option.cpp` lines 548–572. Source hashes and settings are in [mph_source_verification.json](../audit/mph_source_verification.json). The additive model with independent normal SNP effects is a distributional interpretation of the covariance, apart from the numerical ridge; it is not the executed sampling procedure.

## Historical sample chain

A read-only HPCC check on 2026-11-04 established:

| Comparison | Result |
|---|---|
| Saved simulation IDs versus saved MPH GRM IDs | Same 31,301 IDs in the same order |
| Saved simulation versus retained `FineMapping/spg` chip dataset | Same IDs and order |
| Saved simulation versus GWAS panel | Same 31,301-ID set; different row order |
| First 20,000 simulation rows versus combined phenotypes | Same IDs and order |
| First 20,000 simulation rows versus fine-mapping genotypes | Same 20,000-ID set; different row order |
| All ten original subsets, fastGWA files and combined phenotype columns | Exact numeric agreement with simulation columns 1–10 |
| Current files at the old simulation source path | Different labels/order; not the retained simulation source |

Genotypes must be joined to phenotypes by IID. The fine-mapping implementation does this explicitly in `L2_loco_ld.py`, including reordering REGENIE LOCO prediction columns to the genotype FAM order. The old statement that every file shared a single row order was incorrect.

To check the retained chip file against the saved MPH GRM, the audit decoded its 10,000 causal SNPs, standardized them using all 31,301 animals, and recomputed all 31,301 diagonal entries and 4,096 deterministic off-diagonal entries. Maximum absolute errors were 1.081 × 10⁻⁷ and 4.757 × 10⁻⁸, respectively, below the prespecified 10⁻⁴ tolerance. This supports the retained genotype/order provenance; it is not a full matrix or original-file bytewise comparison. No missing calls occurred in the checked causal genotypes.

[Verification code](../scripts/provenance/verify_sample_order.py) and [aggregate evidence](../audit/sample_order_verification.json) retain the calculation and hashes without exporting IDs or genotypes. Realized Var(g)/Var(y) in the analyzed animals ranges from 0.4732 to 0.5285; the design target remains 0.5.

## Deferred inputs

The author identified Professor Jicai Jiang as the source of the Equation 12/theory and maximum-SNP PVE inputs and requested that their provenance investigation be skipped for now. These questions are deferred, not resolved by new calculations. Existing tables and manuscript PVE wording remain unchanged. The inverse-regression PVE scripts are archived because they estimate a different quantity and cannot establish the original realized-effect claim.

## Manuscript corrections

Independent raw-column checks reproduce the REGENIE_LOCO maximum of 20,461 significant SNPs. The former upper bound, 21,484, is the BOLT-LMM maximum.

The largest −log10(P) across all 90 panels is 107.517, directly recorded in REGENIE_LOCO phenotype 8. None of the 90 panel maxima lies between 100 and 103. The plotting code applies a display limit but does not transform 107.517 into approximately 102. The corrected caption reports the full value and identifies its method/phenotype; the origin of the former approximation is not established.

[Raw-column checks](../audit/manuscript_number_checks.tsv) are independent of the shared harmonizer. Both the current corrected-claim ledger and the original discrepancy ledger are retained. The corrected manuscript is a separate local file; the original October DOCX was not overwritten.
