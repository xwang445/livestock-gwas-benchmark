# Data availability and classification

| Class | Included or excluded | Treatment |
|---|---|---|
| A PUBLIC CODE | Analysis scripts, configurations, tests, documentation | Included for local author review; existing license retained |
| B NON-IDENTIFYING DERIVED RESULTS | Counts, rates, lambda summaries, binned predictions, LD bin means, pooled fine-mapping metrics, manuscript figure redraws, comparison aggregates | Included locally; authors must confirm the governing data agreement permits publication before any push |
| C RESTRICTED INPUT/OUTPUT | Animal genotypes, individual phenotypes and IDs, raw GWAS, causal variant identifiers, GRMs, raw LD, detailed PIPs, manuscript DOCX | Excluded from Git; source files remain in authorized locations |
| D TEMPORARY/HISTORICAL | Caches, inventories containing private paths, old outputs, renderer scratch files, third-party source snapshots | Retained privately where relevant; not committed |

The public aggregate tables have no animal IDs, SNP identifiers, individual phenotype values, genotype dosages, GRM entries, raw LD pairs or per-variant association/PIP rows. Source fingerprints and relative filenames identify provenance without distributing the underlying data. Full local/HPCC inventories are Git-ignored.

Synthetic data are generated independently and labeled `SYNTHETIC EXAMPLE FOR CODE VALIDATION ONLY`. Synthetic IDs cannot be confused with the restricted sample labels.

Access to proprietary data requires the authors and data owners; this repository does not grant permission. Array model, imputation production details, institutional approval language and availability terms require human confirmation. No public deposit or external transfer was performed during reconstruction.
