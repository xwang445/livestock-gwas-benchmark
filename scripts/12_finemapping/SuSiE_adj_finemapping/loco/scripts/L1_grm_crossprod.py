#!/usr/bin/env python3
"""L1: per-chromosome standardized genotype cross-products S_k = Z_k Z_k' (20,000 x 20,000).

SNP set = all 135,358 panel SNPs: the BOLT-LMM model SNPs ("135358 SNPs to include in model (i.e., GRM)")
and the REGENIE Step-1 markers ("145 blocks for 135358 variants") are this same set.
Animals = the 20,000 phenotyped animals used by both GWAS (Nused = 20000), in spg.20k.fam order.
Standardization: centre and divide by the empirical SD (n-1), as BOLT ("projecting and normalizing
SNPs"), REGENIE ("residualizing and scaling genotypes") and ld_adjuster do.

Outputs: grm/S_chr{k}.npy (float64), grm/S_all.npy, grm/snp_counts.tsv, qc/L1_grm_qc.txt
(comparison of S_all/M with the existing MPH GRM), qc/bolt_variance_params.tsv.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

t0 = time.time()
bim, fam = read_bim(), read_fam()
n = len(fam)
assert n == 20000 and len(bim) == 135358
bolt_params().to_csv(f"{L}/qc/bolt_variance_params.tsv", sep="\t", index=False)

S_all = np.zeros((n, n))
counts = []
for k in range(1, 19):
    out = f"{L}/grm/S_chr{k}.npy"
    idx = np.flatnonzero(bim.chr.to_numpy() == k)
    if os.path.exists(out):
        S = np.load(out)
    else:
        S = np.zeros((n, n))
        for s in range(0, len(idx), 4000):                   # blocks keep memory modest
            Z, ok = standardize(read_bed_cols(idx[s:s + 4000], n))
            assert ok.all(), "monomorphic SNP"
            S += Z @ Z.T
        np.save(out + ".tmp.npy", S); os.replace(out + ".tmp.npy", out)
    S_all += S
    counts.append(dict(chr=k, M_chr=len(idx), M_loco=len(bim) - len(idx)))
    print(f"chr{k}: {len(idx)} SNPs  {time.time() - t0:.0f}s", flush=True)
np.save(f"{L}/grm/S_all.npy", S_all)
pd.DataFrame(counts).to_csv(f"{L}/grm/snp_counts.tsv", sep="\t", index=False)

# QC against the existing MPH GRM (same 20k animals, all SNPs); MPH lower triangle, column-major
G = S_all / len(bim)
with open(f"{MPH_GRM}.grm.bin", "rb") as f:
    npp = np.frombuffer(f.read(4), np.int32)[0]; nm = np.frombuffer(f.read(4), np.float32)[0]
    tri = np.frombuffer(f.read(), np.float32)
iid = [l.strip() for l in open(f"{MPH_GRM}.grm.iid")]
assert npp == n and iid == fam.IID.tolist()
Mph = np.zeros((n, n), np.float32)
Mph[np.triu_indices(n)] = tri          # column-major lower triangle == row-major upper triangle
Mph = Mph.T / nm
il = np.tril_indices(n, -1)
off_ours, off_mph = G[il], Mph[il].astype(np.float64)
with open(f"{L}/qc/L1_grm_qc.txt", "w") as f:
    print(f"N={n} M={len(bim)}  MPH sum-of-weights={nm}", file=f)
    print(f"diag mean ours={np.diag(G).mean():.6f} mph={np.diag(Mph).mean():.6f}", file=f)
    print(f"diag corr={np.corrcoef(np.diag(G), np.diag(Mph))[0, 1]:.8f}  max|diff|={np.abs(np.diag(G) - np.diag(Mph)).max():.4g}", file=f)
    print(f"offdiag corr={np.corrcoef(off_ours, off_mph)[0, 1]:.8f}  max|diff|={np.abs(off_ours - off_mph).max():.4g}"
          f"  mean ours={off_ours.mean():.3g} mph={off_mph.mean():.3g}", file=f)
    print(f"offdiag SD ours={off_ours.std():.5f}", file=f)
print(open(f"{L}/qc/L1_grm_qc.txt").read())
print(f"done {time.time() - t0:.0f}s")
