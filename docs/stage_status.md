# Stage status

“Ready” means technically prepared for review, not permission to publish restricted derivatives.

| Stage | Status | Source | Validated against manuscript | GitHub ready |
|---|---|---|---|---|
| Genotype QC | PARTIAL | HPCC QC script and prior chip/INFO audits | Counts supported by historical evidence; no new QC fit | Code ready; source identity unresolved |
| LD decay | COMPLETE | HPCC pairwise LD and harmonized binning code | 1,000 counts/means independently matched | Yes, aggregate review candidate |
| Architecture / simulation | AMBIGUOUS | Historical R/MPH commands, versioned MPH source | Design recovered; new simulation not run | Code ready with implementation caveats |
| Nine GWAS configurations | PARTIAL | HPCC runners/logs and 90 local results | All 90 panels validated downstream; no new fits | Code ready; GAPIT version unresolved |
| Harmonization | COMPLETE | Nine observed output schemas | 90 × 135,358 rows and common coordinates checked | Yes; per-variant exports remain restricted |
| Significant SNPs/loci | COMPLETE | Local historical Figure 3 logic | All 90 SNP and 90 locus cells matched | Yes |
| Power / FPR | COMPLETE | Local pwer_analysis.R and 90 outputs | 2,970 historical rows for each curve matched | Yes |
| Spillover | COMPLETE | Local Figure 6A code | 90 rows; nine rounded rates; denominator 1,480 | Yes |
| Lambda GC | COMPLETE | Local historical statistic definitions | All 90 historical values matched | Yes; numerical tail limitation disclosed |
| Overlap / stability | COMPLETE | Local Jaccard code, manuscript caption | 30 displayed spot checks; full tables recomputed | Yes; stability is supporting output |
| Equation 12 / upstream theory | REQUIRES_RECOMPUTATION | Cached prediction table; generator not found | NOT_RUN upstream | No end-to-end claim |
| Equation 13 downstream | COMPLETE | Cached table and local production R | 1,600 replicate bins, 160 means/intervals matched | Yes with upstream limitation |
| Fine-mapping fitting | PARTIAL | HPCC primary LOCO/RSS code and saved fits | Fit counts recounted; fitting not repeated | Code ready with restricted inputs |
| Fine-mapping evaluation | COMPLETE | Saved converged PIPs and MAF annotations | All four counts, precision and MAF expectations match | Yes; MAF not recalculated from genotypes |
| Figures / supplements | COMPLETE | Verified aggregates plus restricted GWAS for Manhattan | 8 main + S1–S10 + central-null supplement rendered | Yes after data-agreement review |
| Synthetic / tests / CI | COMPLETE | New artificial fixtures and tests | Local tests/static checks passed; remote CI not run | Yes |
| Maximum SNP PVE | AMBIGUOUS | Recovered LSMR outputs differ from manuscript | Manuscript range not reproduced | No validated claim |
| Metadata / publication | PARTIAL | Author review required | Data rights, array/imputation and citation unresolved | No push authorized |
