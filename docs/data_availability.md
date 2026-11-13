# Data availability and classification

| Class | Included or excluded | Treatment |
|---|---|---|
| A PUBLIC CODE | Analysis scripts, configurations, tests, documentation | Included for local author review; existing license retained |
| B NON-IDENTIFYING DERIVED RESULTS | Counts, rates, lambda summaries, binned predictions, LD bin means, pooled fine-mapping metrics, recovered-effect/PVE-group summaries, figure redraws, editable aggregate tables and captions, comparison aggregates | Included after review for identifying fields; underlying access terms remain an author-provided item |
| C RESTRICTED INPUT/OUTPUT | Animal genotypes, individual phenotypes and IDs, raw GWAS, causal variant identifiers, GRMs, raw LD, detailed PIPs, manuscript DOCX | Excluded from Git; source files remain in authorized locations |
| D TEMPORARY/HISTORICAL | Caches, inventories containing private paths, old outputs, renderer scratch files, third-party source snapshots | Retained privately where relevant; not committed |

The public aggregate tables have no animal IDs, SNP identifiers, individual phenotype values, genotype dosages, GRM entries, raw LD pairs or per-variant association/PIP rows. Source fingerprints and relative filenames identify provenance without distributing the underlying data. Full local/HPCC inventories are Git-ignored.

Synthetic data are generated independently and labeled `SYNTHETIC EXAMPLE FOR CODE VALIDATION ONLY`. Synthetic IDs cannot be confused with the restricted sample labels.

Access to proprietary data requires the authors and data owners; this repository does not grant permission. Array model, imputation production details, institutional approval language and availability terms require author confirmation. The revision distributes no individual-level genotypes or identifying sample information.

Fine-mapping is explicitly covered: the [RSS fitter and exact sensitivity calculator](../scripts/revision_2026-11-05/susie/README.md), [precision–recall display source](../finemapping/pr_curve/README.md), and [SLEMM lead-summary entry point](../finemapping/slemm/README.md) are included. Actual aggregate tables and separate historical/revised fit settings are indexed in [the revision record](revision_2026-11-05.md). Dense LD, truth SNP identifiers, per-SNP PIPs, credible-set membership and detailed lead-SNP calls remain restricted.
