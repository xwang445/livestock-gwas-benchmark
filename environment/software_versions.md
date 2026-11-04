# Software evidence

| Tool | Version evidence / scope |
|---|---|
| BOLT-LMM | 2.5, historical method materials/logs |
| REGENIE | 4.1.2, historical method materials/logs |
| GCTA | 1.94.1, historical fastGWA logs |
| SLEMM | 0.90.1, historical method materials/logs |
| GAPIT | Manuscript says 3.5; retained source says 2025.10.03 GAPIT 4.0; executed source version unresolved |
| MPH | 0.55.1, retained simulation log and inspected source commit 13ffe633e0170d35e01e6164fc2bb4e01f35e2f2 |
| PLINK LD | 1.9.0-b.8 (22 Oct 2024), ld_10Mb.log |
| PLINK2 QC | Exact executed version not established in this reconstruction |
| susieR | 0.14.2, both primary fitters assert packageVersion |
| Local R | 4.3.2; figure packages data.table, ggplot2, patchwork |
| Local Python | 3.10.4; NumPy 2.2.6, pandas 2.3.3, SciPy 1.14.1, matplotlib 3.10.9, PyYAML 6.0.3 |

Validation dependencies are pinned in requirements-validation.txt. This is not a complete lockfile of the original HPCC software environment. CI targets Python 3.11 and R 4.3.2; it was prepared but not run remotely.
