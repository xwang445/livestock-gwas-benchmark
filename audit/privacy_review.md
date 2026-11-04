# Local publication-file review

Class A: code, configuration, synthetic generator, tests and prose. Class B: aggregate tables, reference counts, binned summaries and 19 figure redraws. Class C excluded: genotypes, real animal IDs/phenotypes, raw GWAS, SNP-specific truth/PIPs, raw GRMs/LD, manuscript source. Class D excluded: private caches, snapshots and full local/HPCC inventories.

Automated aggregate-column checks reject SNP/animal ID, beta/SE, genotype and PIP columns. All public derived tables passed. Figure text contains chromosome labels and summary annotations, not animal/SNP labels. Input provenance uses relative filenames and hashes. All publication files are below 5 MB each. No executable script contains a hard-coded personal workstation/HPCC source path.

This checks content, not legal permission. The authors must confirm the data agreement allows these aggregate derivatives and detailed figures before a separately authorized push. Existing MIT/citation metadata are preserved; the GAPIT library and production executables are not vendored. Full inventories are installed locally but Git-ignored.
