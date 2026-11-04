# Locus overlap

Recovered production logic is in `scripts/benchmark_core.py summarize`. It writes phenotype-specific Jaccard/intersection/union and cross-replicate stability tables. Strict P < 5e-8 and floor(bp/500000) define loci. Jaccard is NA whenever either set is empty, matching the October Figure 5 caption. Thirty displayed values were independently checked against historical figure references. All pairwise tables are regenerated; empty comparisons are tested.
