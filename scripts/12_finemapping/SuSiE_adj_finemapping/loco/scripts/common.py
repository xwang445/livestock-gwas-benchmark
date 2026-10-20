"""Shared paths and helpers for the LOCO SuSiE-adj analysis (W2Mb_N20k)."""
import os
from pathlib import Path
import re
import numpy as np
import pandas as pd

W = str(Path(__file__).resolve().parents[2])
L = f"{W}/loco"
BENCH = str(Path(W).parent)
GENO = f"{W}/genotype/spg.20k"             # 20,000 GWAS animals x 135,358 SNPs (= BOLT/REGENIE SNP set)
MPH_GRM = f"{W}/grm/spg.20k"               # existing non-LOCO MPH GRM (QC comparison only)
MANIFEST = f"{W}/loci/manifest.tsv"        # existing locus definition (unchanged)
REGIONS = f"{BENCH}/eligible_regions.csv"
CAUSAL = f"{BENCH}/snp.info.chr1_18.10k.csv"
PHENOS = [f"pheno{i}" for i in range(1, 11)]
METHOD_DIR = {"BOLT-LMM": "BOLT_LMM", "REGENIE-LOCO": "REGENIE"}


def read_bim():
    return pd.read_csv(f"{GENO}.bim", sep="\t", header=None, names=["chr", "SNP", "cm", "pos", "A1", "A2"])


def read_fam():
    return pd.read_csv(f"{GENO}.fam", sep=r"\s+", header=None, names=["FID", "IID", "P", "M", "S", "Y"], dtype=str)


_LUT = np.array([2.0, np.nan, 1.0, 0.0])   # PLINK .bed 2-bit codes -> count of .bim A1 allele


def read_bed_cols(idx, n):
    """Genotypes (n x len(idx), float64, A1 counts) for SNP row indices `idx` of the SNP-major .bed."""
    nb = (n + 3) // 4
    bed = np.memmap(f"{GENO}.bed", dtype=np.uint8, mode="r", offset=3)
    raw = np.asarray(bed.reshape(-1, nb)[np.asarray(idx)])                      # (m, nb)
    codes = np.stack([(raw >> s) & 3 for s in (0, 2, 4, 6)], axis=2).reshape(len(idx), -1)[:, :n]
    return _LUT[codes].T.copy()                                                 # (n, m)


def standardize(X):
    """ld_adjuster v1.1.1 convention: center, divide by SD with n-1 denominator (in place)."""
    if np.isnan(X).any():
        raise ValueError("missing genotypes not expected in spg.20k")
    X -= X.mean(0)
    sd = np.sqrt((X ** 2).sum(0) / (X.shape[0] - 1))
    ok = sd > 1e-6
    X[:, ok] /= sd[ok]
    return X, ok


def bolt_params():
    """Per-phenotype BOLT-LMM variance parameters parsed from BOLTLMM/boltlmm.log."""
    txt = open(f"{BENCH}/BOLTLMM/boltlmm.log").read()
    rows = []
    for blk in txt.split("--phenoCol ")[1:]:
        ph = blk.split()[0]
        h2 = float(re.search(r"Estimated \(pseudo-\)heritability: h2g = ([\d.]+)", blk).group(1))
        m = re.search(r"Variance params: sigma\^2_K = ([\d.eE+-]+), logDelta = ([\d.eE+-]+)", blk)
        cal = re.search(r"Calibration: ([\d.]+) \(([\d.]+)\)", blk)
        lg = re.search(r"lambdaGC: ([\d.]+)", blk)
        chunks = re.search(r"Assigning SNPs to (\d+) chunks for leave-out", blk)
        ld = float(m.group(2))
        rows.append(dict(phenotype=ph, h2g_bolt=h2, sigma2_K=float(m.group(1)), logDelta=ld,
                         eta=float(np.exp(-ld)), h2_from_delta=1 / (1 + np.exp(ld)),
                         calibration=float(cal.group(1)), calibration_se=float(cal.group(2)),
                         lambdaGC=float(lg.group(1)), n_loco_chunks=int(chunks.group(1))))
    d = pd.DataFrame(rows)
    assert sorted(d.phenotype) == sorted(PHENOS), d.phenotype.tolist()
    return d
