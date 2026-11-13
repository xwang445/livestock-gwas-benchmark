# Analysis code

Run `python scripts/workflow.py --help` from the repository root. [docs/workflow.md](../docs/workflow.md) is the authoritative execution guide.

| Directory or module | Role |
|---|---|
| `workflow.py` | Supported entry point for checks, saved-result analysis and figures |
| `publish_current_summaries.py` | Rebuild/check current root views from committed aggregate sources |
| `benchmark_core.py` | Harmonization and downstream statistical implementation |
| `01_genotype_preparation/` | Upstream genotype panel construction |
| `02_ld_analysis/` | LD generation and saved-LD verification |
| `03_simulation/` | Causal selection, MPH simulation and phenotype formatting |
| `04_gwas/` | Original GWAS runners; the completed GCTA MLMA runner is in `revision_2026-11-05/run_gcta_mlma.py` |
| `05_downstream/` | Guide to the shared downstream calculations |
| `12_finemapping/` | Saved-fit evaluation and retained fitting implementation |
| `13_figures/` | Figure generators |
| `revision_2026-11-11/` | Recovered-effect lead/PVE and GWAS-coefficient PVE summaries |
| `submission_2026-11-12/` | Portable submission figures and pooled Table S1 |
| `provenance/` | Read-only sample-chain and manuscript-number checks |

Older stage scripts are preserved in [archive](../archive/README.md), outside the supported workflow. The directory numbers preserve established paths; the workflow guide defines execution order.
