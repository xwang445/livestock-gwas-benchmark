#!/usr/bin/env bash
set -euo pipefail

# Purpose
#   Estimate maximum single-SNP PVE from a causal-SNP list and phenotype table.
# Manuscript section or figure
#   Realized per-SNP PVE and power
# Required software
#   PLINK2, Python3 with numpy/pandas/scipy
# Required inputs
#   genotype prefix, phenotype file, causal SNP list
# Outputs
#   max_snp_pve_results/*.tsv, summary report, optional plot
# Example command
#   bash scripts/05_effect_size_power/calculate_max_snp_pve.sh
# Note
#   Proprietary genotype and phenotype inputs are not distributed with this repository.

#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# calculate_max_snp_pve.sh
#
# Purpose
#   1. Extract the 10,000 simulated causal SNPs from the PLINK
#      genotype dataset.
#   2. Keep individuals present in pheno.sim.csv.
#   3. Export additive genotype dosages.
#   4. Standardize each SNP as:
#          Z = (X - 2p) / sqrt(2pq)
#   5. For phenotype 1-10, pair phenotype r with breeding
#      value gvr (for example, phenotype 1 with gv1).
#   6. Estimate true effects for each phenotype:
#          beta_r = (Z'Z)^(-1) Z'g_r
#      using a numerically stable least-squares solver.
#   7. Calculate total PVE for each phenotype:
#          PVE_total,r = Var(Z beta_r) / Var(y_r)
#   8. Calculate each SNP's component PVE:
#          PVE_jr = Var(Z_j beta_jr) / Var(y_r)
#   9. Report the maximum single-SNP PVE for each phenotype
#      and the maximum across all 10 phenotypes.
#
# Run from any directory:
#   bash calculate_max_snp_pve.sh
#
# Optional environment variables:
#   PLINK2=/path/to/plink2
#   N_REPLICATES=10
#   THREADS=8
# ============================================================

# ----------------------------
# User paths
# ----------------------------

WORKDIR="${PROJECT_ROOT}"

GENO_PREFIX="${GENO_PREFIX}"

PHENO_FILE="${PHENO_FILE}"
CAUSAL_SNP_FILE="${CAUSAL_SNP_FILE}"

OUTDIR="${WORKDIR}/max_snp_pve_results"
EXTRACT_PREFIX="${OUTDIR}/causal10k_simulation_samples"
DOSAGE_FILE="${EXTRACT_PREFIX}.raw"

N_REPLICATES="${N_REPLICATES:-10}"
THREADS="${THREADS:-8}"
PLINK2="${PLINK2:-plink2}"

mkdir -p "${OUTDIR}"

LOG_FILE="${OUTDIR}/calculate_max_snp_pve_20k.log"
exec > >(tee "${LOG_FILE}") 2>&1

echo "============================================================"
echo "Starting maximum SNP PVE calculation using all 31,301 simulation individuals and the exact MPH genotype dataset"
echo "Date: $(date)"
echo "Working directory: ${WORKDIR}"
echo "Output directory: ${OUTDIR}"
echo "Phenotypes analyzed: 1-${N_REPLICATES}"
echo "============================================================"

# ----------------------------
# Validate input files
# ----------------------------

for file in \
    "${GENO_PREFIX}.bed" \
    "${GENO_PREFIX}.bim" \
    "${GENO_PREFIX}.fam" \
    "${PHENO_FILE}" \
    "${CAUSAL_SNP_FILE}"
do
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Required input file is missing or empty:"
        echo "  ${file}"
        exit 1
    fi
done

if ! command -v "${PLINK2}" >/dev/null 2>&1; then
    echo "ERROR: Cannot find PLINK 2 executable: ${PLINK2}"
    echo "Load the appropriate module or set PLINK2=/path/to/plink2"
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is not available."
    exit 1
fi

python3 - <<'PY'
required = ["numpy", "pandas", "scipy"]
missing = []

for package in required:
    try:
        __import__(package)
    except ImportError:
        missing.append(package)

if missing:
    raise SystemExit(
        "ERROR: Missing Python packages: " + ", ".join(missing)
    )

print("Required Python packages are available.")
PY

# ----------------------------
# Prepare causal-SNP list
# ----------------------------

CAUSAL_LIST="${OUTDIR}/causal_10k_snps.txt"

awk -F',' '
NR == 1 {
    if ($1 != "SNP") {
        print "ERROR: Expected SNP header in causal SNP file." > "/dev/stderr"
        exit 1
    }
    next
}
NF >= 1 && $1 != "" {
    gsub(/\r/, "", $1)
    print $1
}
' "${CAUSAL_SNP_FILE}" > "${CAUSAL_LIST}"

N_CAUSAL_REQUESTED=$(wc -l < "${CAUSAL_LIST}")

echo "Causal SNPs requested: ${N_CAUSAL_REQUESTED}"

if [[ "${N_CAUSAL_REQUESTED}" -ne 10000 ]]; then
    echo "WARNING: Expected 10,000 causal SNPs but found ${N_CAUSAL_REQUESTED}."
fi

# ----------------------------
# Prepare individual keep file
# ----------------------------

IID_LIST="${OUTDIR}/simulation_iids.txt"
KEEP_FILE="${OUTDIR}/simulation_individuals.keep"

awk -F',' '
NR >= 2 && NR <= 20001 && $1 != "" {
    gsub(/\r/, "", $1)
    print $1
}
' "${PHENO_FILE}" > "${IID_LIST}"

awk '
NR == FNR {
    keep[$1] = 1
    next
}
($2 in keep) {
    print $1, $2
}
' "${IID_LIST}" "${GENO_PREFIX}.fam" > "${KEEP_FILE}"

N_PHENO=$(wc -l < "${IID_LIST}")

if [[ "${N_PHENO}" -ne 20000 ]]; then
    echo "ERROR: Expected 20,000 phenotype IDs but found ${N_PHENO}."
    exit 1
fi
N_KEEP=$(wc -l < "${KEEP_FILE}")

echo "Individuals in phenotype file: ${N_PHENO}"
echo "Individuals matched in PLINK FAM: ${N_KEEP}"

if [[ "${N_KEEP}" -ne 20000 ]]; then
    echo "ERROR: Expected 20,000 matched individuals but found ${N_KEEP}."
    exit 1
fi

# ----------------------------
# Extract genotypes
# ----------------------------

echo "Extracting causal SNPs and matched individuals with PLINK 2..."

"${PLINK2}" \
    --bfile "${GENO_PREFIX}" \
    --extract "${CAUSAL_LIST}" \
    --keep "${KEEP_FILE}" \
    --make-bed \
    --threads "${THREADS}" \
    --out "${EXTRACT_PREFIX}"

N_SNPS_EXTRACTED=$(wc -l < "${EXTRACT_PREFIX}.bim")
N_SAMPLES_EXTRACTED=$(wc -l < "${EXTRACT_PREFIX}.fam")

echo "SNPs extracted: ${N_SNPS_EXTRACTED}"
echo "Individuals extracted: ${N_SAMPLES_EXTRACTED}"

if [[ "${N_SNPS_EXTRACTED}" -eq 0 ]]; then
    echo "ERROR: No causal SNPs were extracted."
    exit 1
fi

if [[ "${N_SNPS_EXTRACTED}" -ne "${N_CAUSAL_REQUESTED}" ]]; then
    echo "WARNING: Requested ${N_CAUSAL_REQUESTED} SNPs but extracted ${N_SNPS_EXTRACTED}."
    echo "Check ${EXTRACT_PREFIX}.log for unmatched variant IDs."
fi

# ----------------------------
# Export additive dosages
# ----------------------------

echo "Exporting additive genotype dosages..."

"${PLINK2}" \
    --bfile "${EXTRACT_PREFIX}" \
    --export A \
    --threads "${THREADS}" \
    --out "${EXTRACT_PREFIX}"

if [[ ! -s "${DOSAGE_FILE}" ]]; then
    echo "ERROR: PLINK dosage file was not created:"
    echo "  ${DOSAGE_FILE}"
    exit 1
fi

# ----------------------------
# Estimate effects and PVE
# ----------------------------

echo "Estimating SNP effects and PVE..."

export PHENO_FILE
export DOSAGE_FILE
export BIM_FILE="${EXTRACT_PREFIX}.bim"
export OUTDIR
export N_REPLICATES

python3 - <<'PY'
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse.linalg import lsmr


pheno_file = Path(os.environ["PHENO_FILE"])
dosage_file = Path(os.environ["DOSAGE_FILE"])
bim_file = Path(os.environ["BIM_FILE"])
outdir = Path(os.environ["OUTDIR"])
n_replicates = int(os.environ["N_REPLICATES"])

replicates = list(range(1, n_replicates + 1))

# Numerical solver settings
lsmr_atol = 1e-8
lsmr_btol = 1e-8
lsmr_maxiter = 5000

print("Reading phenotype and genetic-value columns...")

phenotype_columns = [str(r) for r in replicates]
genetic_value_columns = [f"gv{r}" for r in replicates]
required_pheno_columns = ["iid"] + phenotype_columns + genetic_value_columns

pheno = pd.read_csv(
    pheno_file,
    usecols=required_pheno_columns,
    dtype={"iid": str},
)

if pheno["iid"].duplicated().any():
    examples = pheno.loc[
        pheno["iid"].duplicated(),
        "iid",
    ].head().tolist()
    raise ValueError(
        "Duplicated IID values in pheno.sim.csv. "
        f"Examples: {examples}"
    )

print(f"Phenotype records read: {len(pheno):,}")

print("Reading PLINK additive dosage file...")

geno = pd.read_csv(
    dosage_file,
    sep=r"\s+",
    dtype={"FID": str, "IID": str},
)

metadata_columns = {
    "FID",
    "IID",
    "PAT",
    "MAT",
    "SEX",
    "PHENOTYPE",
}

missing_metadata = metadata_columns.difference(geno.columns)

if missing_metadata:
    raise ValueError(
        "Missing expected PLINK metadata columns: "
        + ", ".join(sorted(missing_metadata))
    )

snp_columns = [
    column for column in geno.columns
    if column not in metadata_columns
]

print(f"Genotype individuals read: {len(geno):,}")
print(f"Genotype SNPs read: {len(snp_columns):,}")

# Align phenotype rows to PLINK genotype order.
pheno_indexed = pheno.set_index("iid")

matched_mask = geno["IID"].isin(pheno_indexed.index)

if not matched_mask.all():
    n_unmatched = int((~matched_mask).sum())
    print(
        f"WARNING: Removing {n_unmatched:,} genotype individuals "
        "without phenotype records."
    )
    geno = geno.loc[matched_mask].reset_index(drop=True)

pheno_aligned = pheno_indexed.loc[geno["IID"]].reset_index()

if not np.array_equal(
    geno["IID"].to_numpy(),
    pheno_aligned["iid"].to_numpy(),
):
    raise RuntimeError("Individual alignment failed.")

print(f"Matched individuals used: {len(geno):,}")

# Construct genotype dosage matrix.
X = geno[snp_columns].to_numpy(dtype=np.float32, copy=True)
del geno

# Mean-impute missing dosages.
dosage_means = np.nanmean(X, axis=0, dtype=np.float64)
has_data = np.isfinite(dosage_means)

if not has_data.all():
    n_removed = int((~has_data).sum())
    print(f"Removing {n_removed:,} SNPs with no observed genotypes.")
    X = X[:, has_data]
    dosage_means = dosage_means[has_data]
    snp_columns = [
        snp for snp, keep in zip(snp_columns, has_data)
        if keep
    ]

missing_rows, missing_cols = np.where(np.isnan(X))

if len(missing_rows) > 0:
    print(f"Mean-imputing {len(missing_rows):,} missing dosages.")
    X[missing_rows, missing_cols] = dosage_means[missing_cols]

# Standardize as (X - 2p) / sqrt(2pq).
p = X.mean(axis=0, dtype=np.float64) / 2.0
q = 1.0 - p
two_pq = 2.0 * p * q

polymorphic = (
    np.isfinite(two_pq)
    & (two_pq > 1e-12)
    & (p > 0.0)
    & (p < 1.0)
)

if not polymorphic.all():
    n_removed = int((~polymorphic).sum())
    print(f"Removing {n_removed:,} monomorphic or near-monomorphic SNPs.")
    X = X[:, polymorphic]
    p = p[polymorphic]
    q = q[polymorphic]
    two_pq = two_pq[polymorphic]
    snp_columns = [
        snp for snp, keep in zip(snp_columns, polymorphic)
        if keep
    ]

Z = X
Z -= (2.0 * p).astype(np.float32)
Z /= np.sqrt(two_pq).astype(np.float32)
del X

n_individuals, n_snps = Z.shape

print(
    f"Final Z matrix: {n_individuals:,} individuals x "
    f"{n_snps:,} SNPs"
)

# Calculate observed sample variance for each standardized SNP.
snp_variance = Z.var(axis=0, ddof=1, dtype=np.float64)

# Read BIM annotation.
bim = pd.read_csv(
    bim_file,
    sep=r"\s+",
    header=None,
    names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
    dtype={
        "CHR": str,
        "SNP": str,
        "CM": float,
        "BP": int,
        "A1": str,
        "A2": str,
    },
)

# PLINK --export A appends the counted allele to each SNP ID.
# Example:
#   BIM SNP:       1_273618278_C
#   dosage column: 1_273618278_C_C
marker_ids = [
    re.sub(r"_[ACGT]$", "", column)
    for column in snp_columns
]

summary_records = []
all_effect_tables = []

for replicate in replicates:
    phenotype_column = str(replicate)
    genetic_value_column = f"gv{replicate}"

    print()
    print("=" * 72)
    print(f"Analyzing phenotype {replicate} with gv{replicate}")
    print("=" * 72)

    y = pheno_aligned[phenotype_column].to_numpy(dtype=np.float64)
    g = pheno_aligned[genetic_value_column].to_numpy(dtype=np.float64)

    if not np.isfinite(y).all():
        raise ValueError(
            f"Missing or non-finite phenotype values in replicate {replicate}."
        )

    if not np.isfinite(g).all():
        raise ValueError(
            f"Missing or non-finite genetic values in replicate {replicate}."
        )

    y_centered = y - y.mean()
    g_centered = g - g.mean()

    var_y = float(np.var(y_centered, ddof=1))
    var_g = float(np.var(g_centered, ddof=1))

    if var_y <= 0:
        raise ValueError(
            f"Phenotypic variance is non-positive in replicate {replicate}."
        )

    # This least-squares solution is algebraically equivalent to
    # beta = (Z'Z)^(-1) Z'g when Z has full column rank, while
    # avoiding explicit construction and inversion of Z'Z.
    solution = lsmr(
        Z,
        g_centered,
        atol=lsmr_atol,
        btol=lsmr_btol,
        maxiter=lsmr_maxiter,
        show=False,
    )

    beta = np.asarray(solution[0], dtype=np.float64)
    stopping_code = int(solution[1])
    iterations = int(solution[2])
    residual_norm = float(solution[3])
    condition_estimate = float(solution[6])

    g_hat = np.asarray(Z @ beta, dtype=np.float64)
    var_g_hat = float(np.var(g_hat, ddof=1))

    reconstruction_correlation = float(
        np.corrcoef(g_centered, g_hat)[0, 1]
    )
    reconstruction_r2 = reconstruction_correlation ** 2
    reconstruction_mse = float(
        np.mean((g_centered - g_hat) ** 2)
    )

    observed_pve = var_g / var_y
    reconstructed_pve = var_g_hat / var_y

    # Per-SNP PVE based on each SNP's own component:
    #   Var(Z_j beta_j) / Var(y)
    #
    # These values are useful for identifying the largest SNP,
    # but they are not additive when SNPs are correlated by LD.
    component_variance = (beta ** 2) * snp_variance
    pve_component = component_variance / var_y

    max_index = int(np.nanargmax(pve_component))
    max_snp = marker_ids[max_index]
    max_snp_pve = float(pve_component[max_index])
    max_beta = float(beta[max_index])

    beta_per_allele = beta / np.sqrt(two_pq)

    effect_table = pd.DataFrame({
        "PHENOTYPE": replicate,
        "SNP_COLUMN": snp_columns,
        "SNP": marker_ids,
        "ALLELE_FREQUENCY": p,
        "TWO_PQ": two_pq,
        "Z_VARIANCE": snp_variance,
        "BETA_STANDARDIZED": beta,
        "BETA_PER_ALLELE": beta_per_allele,
        "COMPONENT_VARIANCE": component_variance,
        "PVE_COMPONENT": pve_component,
    })

    effect_table = effect_table.merge(
        bim,
        on="SNP",
        how="left",
        sort=False,
    )

    effect_table["ABS_BETA_STANDARDIZED"] = (
        effect_table["BETA_STANDARDIZED"].abs()
    )

    effect_table = effect_table.sort_values(
        "PVE_COMPONENT",
        ascending=False,
    ).reset_index(drop=True)

    effect_output = outdir / f"snp_pve_first20k_pheno{replicate}.tsv"
    effect_table.to_csv(
        effect_output,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    all_effect_tables.append(
        effect_table[
            [
                "PHENOTYPE",
                "SNP",
                "CHR",
                "BP",
                "A1",
                "A2",
                "ALLELE_FREQUENCY",
                "BETA_STANDARDIZED",
                "BETA_PER_ALLELE",
                "PVE_COMPONENT",
            ]
        ].copy()
    )

    summary_records.append({
        "PHENOTYPE": replicate,
        "BREEDING_VALUE_COLUMN": genetic_value_column,
        "N_INDIVIDUALS": n_individuals,
        "N_SNPS": n_snps,
        "VAR_Y": var_y,
        "VAR_G": var_g,
        "OBSERVED_PVE_VAR_G_OVER_VAR_Y": observed_pve,
        "VAR_G_RECONSTRUCTED": var_g_hat,
        "RECONSTRUCTED_TOTAL_PVE": reconstructed_pve,
        "G_RECONSTRUCTION_CORRELATION":
            reconstruction_correlation,
        "G_RECONSTRUCTION_R2": reconstruction_r2,
        "G_RECONSTRUCTION_MSE": reconstruction_mse,
        "MAX_SNP": max_snp,
        "MAX_SNP_BETA_STANDARDIZED": max_beta,
        "MAX_SNP_PVE": max_snp_pve,
        "MAX_SNP_PVE_PERCENT": 100.0 * max_snp_pve,
        "LSMR_STOPPING_CODE": stopping_code,
        "LSMR_ITERATIONS": iterations,
        "LSMR_RESIDUAL_NORM": residual_norm,
        "LSMR_CONDITION_ESTIMATE": condition_estimate,
    })

    print(f"Var(y): {var_y:.10g}")
    print(f"Var(g): {var_g:.10g}")
    print(f"Observed PVE Var(g)/Var(y): {observed_pve:.10g}")
    print(f"Reconstructed total PVE: {reconstructed_pve:.10g}")
    print(f"Reconstruction R2: {reconstruction_r2:.10g}")
    print(f"Maximum-PVE SNP: {max_snp}")
    print(f"Maximum SNP PVE: {max_snp_pve:.10g}")
    print(f"Maximum SNP PVE (%): {100.0 * max_snp_pve:.8g}")
    print(f"Saved replicate table: {effect_output}")

summary = pd.DataFrame(summary_records)
summary_output = outdir / "pve_summary_first20k.tsv"
summary.to_csv(summary_output, sep="\t", index=False)

all_effects = pd.concat(all_effect_tables, ignore_index=True)
all_effects_output = outdir / "all_snp_pve_first20k.tsv.gz"
all_effects.to_csv(
    all_effects_output,
    sep="\t",
    index=False,
    compression="gzip",
)

overall_row = summary.loc[summary["MAX_SNP_PVE"].idxmax()].copy()

overall_summary = pd.DataFrame([{
    "MAXIMUM_ACROSS_PHENOTYPES_PHENOTYPE":
        int(overall_row["PHENOTYPE"]),
    "MAXIMUM_ACROSS_PHENOTYPES_SNP":
        overall_row["MAX_SNP"],
    "MAXIMUM_ACROSS_PHENOTYPES_BETA_STANDARDIZED":
        overall_row["MAX_SNP_BETA_STANDARDIZED"],
    "MAXIMUM_ACROSS_PHENOTYPES_PVE":
        overall_row["MAX_SNP_PVE"],
    "MAXIMUM_ACROSS_PHENOTYPES_PVE_PERCENT":
        overall_row["MAX_SNP_PVE_PERCENT"],
    "MEAN_OF_PHENOTYPE_MAX_PVE":
        summary["MAX_SNP_PVE"].mean(),
    "SD_OF_PHENOTYPE_MAX_PVE":
        summary["MAX_SNP_PVE"].std(ddof=1),
    "MIN_OF_PHENOTYPE_MAX_PVE":
        summary["MAX_SNP_PVE"].min(),
    "MAX_OF_PHENOTYPE_MAX_PVE":
        summary["MAX_SNP_PVE"].max(),
}])

overall_output = outdir / "overall_max_snp_pve_first20k.tsv"
overall_summary.to_csv(overall_output, sep="\t", index=False)

print()
print("=" * 72)
print("PVE and maximum SNP PVE using the exact first 20,000 GWAS individuals")
print("=" * 72)
print(
    summary[
        [
            "PHENOTYPE",
            "BREEDING_VALUE_COLUMN",
            "MAX_SNP",
            "MAX_SNP_BETA_STANDARDIZED",
            "MAX_SNP_PVE",
            "MAX_SNP_PVE_PERCENT",
            "OBSERVED_PVE_VAR_G_OVER_VAR_Y",
            "RECONSTRUCTED_TOTAL_PVE",
            "G_RECONSTRUCTION_R2",
        ]
    ].to_string(index=False)
)

print()
print("=" * 72)
print("Overall maximum")
print("=" * 72)
print(overall_summary.to_string(index=False))
print()
print(f"Summary by phenotype: {summary_output}")
print(f"Overall maximum: {overall_output}")
print(f"All SNP results: {all_effects_output}")
PY

echo
echo "============================================================"
echo "Analysis completed successfully"
echo "Date: $(date)"
echo "Main output files:"
echo "  ${OUTDIR}/pve_summary_first20k.tsv"
echo "  ${OUTDIR}/overall_max_snp_pve_first20k.tsv"
echo "  ${OUTDIR}/all_snp_pve_first20k.tsv.gz"
echo "  ${LOG_FILE}"
echo "============================================================"
