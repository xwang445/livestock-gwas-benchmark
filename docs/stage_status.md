# Stage status

| Stage | Status | Evidence / entry point |
|---|---|---|
| Public checks | Implemented and locally tested | `workflow.py check` |
| Saved GWAS summaries | Complete | 90 panels; `workflow.py summarize` |
| Cached Equation 13 comparisons | Complete downstream | `workflow.py equation13` |
| Saved LD summaries | Complete | 1,000 bins; `workflow.py ld` |
| Saved fine-mapping evaluation | Complete | 40 groups; `workflow.py fine-mapping` |
| Figures | Complete | 19 PDF/PNG pairs; `workflow.py figures` / `manhattan` |
| MPH implementation description | Verified and corrected | Versioned source and corrected manuscript copy |
| Historical sample chain | Verified for retained files | Exact ID/phenotype checks; all diagonal and sampled off-diagonal GRM entries |
| Manuscript numerical corrections | Complete in separate local copy | REGENIE 20,461; maximum −log10(P) 107.517 |
| Equation 12 / maximum SNP PVE provenance | Deferred by author | Professor Jicai Jiang's inputs; no substitution |
| New genotype preparation, simulation, GWAS or SuSiE fitting | Not run in this update | Upstream runners retained |
| GAPIT executed version | Author review | Retained source 4.0 versus manuscript 3.5 |

See [workflow.md](workflow.md) for supported commands and [scientific_provenance.md](scientific_provenance.md) for evidence and scope. Public repository updates remain authorized; manuscript submission and preprint deposition are outside this task.
