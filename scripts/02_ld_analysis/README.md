# LD decay

Run `python scripts/workflow.py ld --profile hpcc --out ../benchmark-ld-check` from the repository root on HPCC. It streams the saved pairwise LD file and exports only 1,000 10 kb bin counts and mean r² values. Bins are left-closed/right-open, with the final 10 Mb endpoint included. The independent recount matched the saved summaries within 10⁻¹⁰.

The source LD calculation used all 31,301 animals. `run_plink_ld.sh` preserves the logged PLINK options for a separately allocated new run; raw LD remains restricted.

The public mean-bin table feeds `scripts/13_figures/plot_ld.R`, which applies the historical endpoint-constrained spline in log10 distance. Use `workflow.py figures` to generate the LD curve with the other aggregate figures.
