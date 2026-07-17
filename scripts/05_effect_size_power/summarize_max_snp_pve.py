#!/usr/bin/env python3
"""Summarize max-SNP PVE outputs from the benchmark workflow.

Purpose
    Summarize maximum single-SNP PVE and realized PVE from the causal-SNP analysis.
Manuscript section or figure
    Realized per-SNP PVE and power
Required software
    Python3 with pandas and matplotlib
Required inputs
    Files emitted by calculate_max_snp_pve.sh
Outputs
    Summary report and optional figure
Example command
    python3 scripts/05_effect_size_power/summarize_max_snp_pve.py max_snp_pve_results
Note
    Proprietary genotype and phenotype inputs are not distributed with this repository.
"""

#!/usr/bin/env python3
"""
summarize_max_snp_pve.py

Summarizes the output of calculate_max_snp_pve.sh:
  - pve_and_max_snp_pve_by_pheno.tsv
  - overall_max_snp_pve.tsv
  - all_phenotypes_snp_effects_pve.tsv.gz

Writes a text report and a bar chart of max-SNP PVE per phenotype.

Usage:
    python3 summarize_max_snp_pve.py [results_dir] [--top N]

Defaults to ./max_snp_pve_results if results_dir is omitted.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_inputs(results_dir: Path):
    by_pheno_file = results_dir / "pve_and_max_snp_pve_by_pheno.tsv"
    overall_file = results_dir / "overall_max_snp_pve.tsv"
    all_effects_file = results_dir / "all_phenotypes_snp_effects_pve.tsv.gz"

    for f in (by_pheno_file, overall_file, all_effects_file):
        if not f.exists():
            raise FileNotFoundError(f"Expected output file not found: {f}")

    by_pheno = pd.read_csv(by_pheno_file, sep="\t")
    overall = pd.read_csv(overall_file, sep="\t")
    all_effects = pd.read_csv(all_effects_file, sep="\t", compression="gzip")

    return by_pheno, overall, all_effects


def build_report(by_pheno: pd.DataFrame, overall: pd.DataFrame,
                  all_effects: pd.DataFrame, top_n: int) -> str:
    lines = []
    add = lines.append

    n_pheno = len(by_pheno)
    n_snps = int(by_pheno["N_SNPS"].iloc[0])
    n_ind = int(by_pheno["N_INDIVIDUALS"].iloc[0])

    add("=" * 72)
    add("MAX SNP PVE SUMMARY")
    add("=" * 72)
    add(f"Phenotypes analyzed:   {n_pheno}")
    add(f"Individuals:           {n_ind:,}")
    add(f"SNPs (causal set):     {n_snps:,}")
    add("")

    add("-" * 72)
    add("Total PVE (phenotypic variance explained by all causal SNPs)")
    add("-" * 72)
    total_pve = by_pheno["OBSERVED_PVE_VAR_G_OVER_VAR_Y"]
    add(f"  mean:   {total_pve.mean():.6f}")
    add(f"  sd:     {total_pve.std(ddof=1):.6f}")
    add(f"  min:    {total_pve.min():.6f}  (phenotype {int(by_pheno.loc[total_pve.idxmin(), 'PHENOTYPE'])})")
    add(f"  max:    {total_pve.max():.6f}  (phenotype {int(by_pheno.loc[total_pve.idxmax(), 'PHENOTYPE'])})")
    add("")

    add("-" * 72)
    add("Maximum single-SNP PVE, per phenotype")
    add("-" * 72)
    view = by_pheno[
        ["PHENOTYPE", "MAX_SNP", "MAX_SNP_BETA_STANDARDIZED", "MAX_SNP_PVE_PERCENT"]
    ].sort_values("MAX_SNP_PVE_PERCENT", ascending=False)
    add(view.to_string(index=False))
    add("")

    add("-" * 72)
    add("Overall maximum across all phenotypes")
    add("-" * 72)
    row = overall.iloc[0]
    add(f"  Phenotype:        {int(row['MAXIMUM_ACROSS_PHENOTYPES_PHENOTYPE'])}")
    add(f"  SNP:              {row['MAXIMUM_ACROSS_PHENOTYPES_SNP']}")
    add(f"  Beta (standardized): {row['MAXIMUM_ACROSS_PHENOTYPES_BETA_STANDARDIZED']:.6f}")
    add(f"  PVE:              {row['MAXIMUM_ACROSS_PHENOTYPES_PVE']:.6f} "
        f"({row['MAXIMUM_ACROSS_PHENOTYPES_PVE_PERCENT']:.4f}%)")
    add(f"  Mean of per-phenotype max PVE: {row['MEAN_OF_PHENOTYPE_MAX_PVE']:.6f}")
    add(f"  SD of per-phenotype max PVE:   {row['SD_OF_PHENOTYPE_MAX_PVE']:.6f}")
    add("")

    add("-" * 72)
    add(f"Top {top_n} SNPs overall, ranked by single-SNP PVE (any phenotype)")
    add("-" * 72)
    top = all_effects.sort_values("PVE_COMPONENT", ascending=False).head(top_n)
    top_view = top[
        ["PHENOTYPE", "SNP", "CHR", "BP", "ALLELE_FREQUENCY",
         "BETA_STANDARDIZED", "PVE_COMPONENT"]
    ].copy()
    top_view["PVE_COMPONENT"] = top_view["PVE_COMPONENT"].round(6)
    top_view["ALLELE_FREQUENCY"] = top_view["ALLELE_FREQUENCY"].round(4)
    top_view["BETA_STANDARDIZED"] = top_view["BETA_STANDARDIZED"].round(4)
    add(top_view.to_string(index=False))
    add("")

    add("-" * 72)
    add("SNPs that are the top hit for more than one phenotype")
    add("-" * 72)
    repeats = by_pheno["MAX_SNP"].value_counts()
    repeats = repeats[repeats > 1]
    if repeats.empty:
        add("  (none)")
    else:
        for snp, count in repeats.items():
            phenos = by_pheno.loc[by_pheno["MAX_SNP"] == snp, "PHENOTYPE"].tolist()
            add(f"  {snp}: top hit in {count} phenotypes {phenos}")
    add("")

    add("=" * 72)

    return "\n".join(lines)


def make_plot(by_pheno: pd.DataFrame, out_path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ordered = by_pheno.sort_values("PHENOTYPE")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(
        ordered["PHENOTYPE"].astype(str),
        ordered["MAX_SNP_PVE_PERCENT"],
        color="steelblue",
    )
    for x, (snp, pve) in enumerate(
        zip(ordered["MAX_SNP"], ordered["MAX_SNP_PVE_PERCENT"])
    ):
        ax.text(x, pve, snp, ha="center", va="bottom", fontsize=7, rotation=90)

    ax.set_xlabel("Phenotype")
    ax.set_ylabel("Maximum single-SNP PVE (%)")
    ax.set_title("Maximum SNP PVE by phenotype")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "results_dir",
        nargs="?",
        default="max_snp_pve_results",
        help="Directory containing calculate_max_snp_pve.sh output (default: %(default)s)",
    )
    parser.add_argument(
        "--top", type=int, default=20,
        help="Number of top SNPs to list overall (default: %(default)s)",
    )
    parser.add_argument(
        "--no-plot", action="store_true",
        help="Skip generating the bar chart",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    by_pheno, overall, all_effects = load_inputs(results_dir)

    report = build_report(by_pheno, overall, all_effects, args.top)
    print(report)

    report_path = results_dir / "summary_report.txt"
    report_path.write_text(report + "\n")
    print(f"\nSaved report to: {report_path}")

    if not args.no_plot:
        plot_path = results_dir / "max_snp_pve_by_phenotype.png"
        make_plot(by_pheno, plot_path)
        print(f"Saved plot to:   {plot_path}")


if __name__ == "__main__":
    sys.exit(main())
