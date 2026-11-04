# Analysis code

Run `python scripts/workflow.py --help` from the repository root. [docs/workflow.md](../docs/workflow.md) is the authoritative execution guide.

| Directory or module | Role |
|---|---|
| `workflow.py` | Supported entry point for checks, saved-result analysis and figures |
| `benchmark_core.py` | Harmonization and downstream statistical implementation |
| `01_genotype_preparation/` | Upstream genotype panel construction |
| `02_ld_analysis/` | LD generation and saved-LD verification |
| `03_simulation/` | Causal selection, MPH simulation and phenotype formatting |
| `04_gwas/` | Nine GWAS configurations |
| `05_downstream/` | Guide to the shared downstream calculations |
| `12_finemapping/` | Saved-fit evaluation and retained fitting implementation |
| `13_figures/` | Figure generators |
| `provenance/` | Read-only sample-chain and manuscript-number checks |

Older stage scripts are preserved in [archive](../archive/README.md), outside the supported workflow. The directory numbers preserve established paths; the workflow guide defines execution order.
