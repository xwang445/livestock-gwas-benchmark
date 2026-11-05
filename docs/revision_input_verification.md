# Revision input verification, 5 November 2026

The revision retains the established 135,358-marker, 31,301-animal panel and the
20,000-animal analysis subset. All 111 checks in the private verification passed.
Aggregate results and file hashes are supplied in
[`derived_data/revision_2026-11-05`](../derived_data/revision_2026-11-05/).
Individual identifiers, phenotypes and genotype calls are excluded from this repository.

The ten combined, SLEMM, fastGWA/GCTA and GAPIT phenotype inputs agree exactly with
original simulation columns 1–10 for the first 20,000 simulation rows. Phenotypes are
selected by name: the combined table places `pheno10` before `pheno1`. Simulation
order and panel FAM order differ; matching by IID, rather than row position, is required.
The PGEN PSAM and BED FAM have identical order; PVAR/BIM SNP IDs, positions and
alleles agree (BED A1 is ALT and A2 is REF).

A complete PGEN-to-BED roundtrip using the logged historical PLINK2 version also
reproduced BED, BIM and FAM byte for byte, across all 31,301 animals and 135,358
markers. This verifies hardcall equivalence; it does not equate original imputation
dosages with hardcalls. The temporary genotype outputs remain private.

The actual analysis subset contains 2,707,160,000 genotype calls and zero missing
calls. All ten phenotype columns are finite. The 10,000 fixed causal SNPs are present
in both retained chip and analysis panels, with matching alleles and positions. All
200,000,000 causal genotype calls agree after sample and allele alignment. The retained
historical chip panel is required: the present file at the old logged source path has
different sample labels and cannot be substituted merely because its pathname matches
the old command. See [scientific provenance](scientific_provenance.md) for the recovered
sample chain and MPH source verification.

This validation does not imply access permission to distribute the individual-level
inputs. It verifies the inputs used for the revision and preserves their established
values; it does not rerun simulation or redefine the fixed truth. Exact source commands,
paths and private audit records are retained in the PI handoff outside this public tree.
