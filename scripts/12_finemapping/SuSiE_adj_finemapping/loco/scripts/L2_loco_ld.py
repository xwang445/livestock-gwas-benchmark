#!/usr/bin/env python3
"""L2 (one chromosome k per call): LOCO GRM, its eigendecomposition, REGENIE-matched eta, and the
relatedness-adjusted LD for every manifest locus on chromosome k (both GWAS methods).

LOCO GRM:  G_-k = (S_all - S_k) / M_-k    (standardized SNPs of the 17 other chromosomes, equal weights)
           G_-k = U diag(lam) U'          (float64 eigendecomposition, once per chromosome)
Adjustment (ld_adjuster v1.1.1 algebra, no covariates):
           V = eta * G_-k + I,  X standardized (n-1 SD),  R* = X' V^-1 X = (U'X)' diag(1/(1+eta*lam)) (U'X)
           n_eff = mean(diag R*),  R_adj = R* scaled to unit diagonal.
eta:
  BOLT-LMM     eta_p = 1/delta_p = exp(-logDelta_p), BOLT's own REML estimate for phenotype p, held fixed
               across the 18 leave-out chunks (BOLT estimates variance parameters once per phenotype).
  REGENIE-LOCO eta_{p,k} fitted so that the BLUP-type LOCO prediction eta G_-k V^-1 y best matches REGENIE's
               own Step-1 LOCO prediction for chromosome k (regenie_step1_pheno{p}_1.loco, row k), where y is
               the phenotype residualized on the intercept and scaled, as in REGENIE. Rationale: if
               yhat_-k = eta G V^-1 y, then y - yhat_-k = V^-1 y, so REGENIE's Step-2 OLS numerator x'(y - yhat_-k)
               equals x'V^-1 y and its z-scores have correlation X'V^-1 X, i.e. the ld_adjuster form.

Usage: python3 L2_loco_ld.py CHR
Outputs: grm/loco_chr{k}.grm.{bin,iid} (MPH format), grm/eig_chr{k}.{lam,U}.npy,
         qc/regenie_eta_chr{k}.tsv, ld/<METHOD>/<pheno>/<locus_id>/{R.bin,snps.tsv,summary.tsv}, status/L2_chr{k}.DONE
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import scipy.linalg as sl
from scipy.optimize import minimize_scalar

k = int(sys.argv[1]); t0 = time.time()
log = lambda *a: print(f"[chr{k} {time.time() - t0:6.0f}s]", *a, flush=True)
if os.path.exists(f"{L}/status/L2_chr{k}.DONE"):
    log("already DONE"); sys.exit(0)
bim, fam = read_bim(), read_fam(); n = len(fam)
cnt = pd.read_csv(f"{L}/grm/snp_counts.tsv", sep="\t").set_index("chr")
M_loco = int(cnt.loc[k, "M_loco"]); assert M_loco == len(bim) - (bim.chr == k).sum()

# ---------------------------------------------------------------- LOCO GRM + eigendecomposition
eig = f"{L}/grm/eig_chr{k}"
if os.path.exists(f"{eig}.U.npy"):                               # restart: reuse the saved decomposition
    lam, U = np.load(f"{eig}.lam.npy"), np.load(f"{eig}.U.npy"); log("loaded saved eigendecomposition")
else:
  G = np.load(f"{L}/grm/S_all.npy"); G -= np.load(f"{L}/grm/S_chr{k}.npy"); G /= M_loco
  log(f"G_-k built: M_-k={M_loco} diag mean={np.diag(G).mean():.5f}")
  with open(f"{L}/grm/loco_chr{k}.grm.bin", "wb") as f:           # MPH format, for ld_adjuster validation/reuse
      f.write(np.int32(n).tobytes()); f.write(np.float32(M_loco).tobytes())
      for s in range(0, n, 2000):                                   # column-major lower triangle, scaled by nm
          f.write(np.concatenate([G[j, j:] * M_loco for j in range(s, min(n, s + 2000))]).astype(np.float32).tobytes())
  with open(f"{L}/grm/loco_chr{k}.grm.iid", "w") as f:
      f.write("\n".join(fam.IID) + "\n")
  lam, U = sl.eigh(G, driver="evd", overwrite_a=True, check_finite=False); del G
  log(f"eigh done: lam min={lam.min():.3g} max={lam.max():.4g} n(lam<1e-8)={(lam < 1e-8).sum()}")
  lam = np.clip(lam, 0, None)
  np.save(f"{eig}.lam.npy", lam); np.save(f"{eig}.U.npy.tmp.npy", U); os.replace(f"{eig}.U.npy.tmp.npy", f"{eig}.U.npy")

# ---------------------------------------------------------------- REGENIE-matched eta (per phenotype)
bp = bolt_params().set_index("phenotype")
ph = pd.read_csv(f"{BENCH}/simu_pheno/pheno_combined.txt", sep="\t", dtype={"IID": str}).set_index("IID")
shrink = lambda eta: eta * lam / (1 + eta * lam)
rows = []
for p in PHENOS:
    y = ph.loc[fam.IID, p].to_numpy(float); y = (y - y.mean()) / y.std(ddof=1)
    lo = pd.read_csv(f"{BENCH}/REGENIE/regenie_step1_{p}_1.loco", sep=" ", index_col=0)
    lo = lo.loc[:, ~lo.columns.str.startswith("Unnamed")]          # trailing separator in the .loco file
    # rows 1-23 are written; 19-23 (no SNPs in this panel) are the whole-genome prediction
    assert list(lo.index[:18]) == list(range(1, 19)) and lo.shape[1] == n and set(lo.columns) == set("0_" + fam.IID)
    lo = lo[["0_" + i for i in fam.IID]]
    yhat = lo.loc[k].to_numpy(float)
    a, b = U.T @ y, U.T @ yhat
    sse = lambda le: float(((b - shrink(np.exp(le)) * a) ** 2).sum())
    grid = np.linspace(-8, 6, 141); g0 = grid[np.argmin([sse(g) for g in grid])]
    opt = minimize_scalar(sse, bounds=(g0 - 0.1, g0 + 0.1), method="bounded", options={"xatol": 1e-6})
    eta = float(np.exp(opt.x)); fit = U @ (shrink(eta) * a); fb = U @ (shrink(bp.loc[p, "eta"]) * a)
    rows.append(dict(phenotype=p, chr=k, eta_regenie=eta, h2_equiv=eta / (1 + eta), at_grid_edge=bool(g0 in (grid[0], grid[-1])),
                     var_yhat=yhat.var(), R2_fit=1 - opt.fun / (yhat ** 2).sum(), corr_fit=np.corrcoef(yhat, fit)[0, 1],
                     slope_yhat_on_fit=float(fit @ yhat / (fit @ fit)),
                     eta_bolt=bp.loc[p, "eta"], R2_at_bolt_eta=1 - float(((yhat - fb) ** 2).sum()) / (yhat ** 2).sum(),
                     corr_at_bolt_eta=np.corrcoef(yhat, fb)[0, 1]))
eta_tab = pd.DataFrame(rows); eta_tab.to_csv(f"{L}/qc/regenie_eta_chr{k}.tsv", sep="\t", index=False)
log("REGENIE eta:", " ".join(f"{r.phenotype}={r.eta_regenie:.3f}(R2={r.R2_fit:.3f})" for r in eta_tab.itertuples()))
eta_reg = eta_tab.set_index("phenotype").eta_regenie

# ---------------------------------------------------------------- adjusted LD per locus
man = pd.read_csv(MANIFEST, sep="\t"); man = man[man.chr == k]
bk = bim[bim.chr == k].reset_index()                              # 'index' = row in .bed
Xc, ok = standardize(read_bed_cols(bk["index"].to_numpy(), n)); assert ok.all()
pos = bk.pos.to_numpy()
UtX = U.T @ Xc; del Xc, U                                         # (n x M_k): all SNPs of chromosome k at once
for i, loc in enumerate(man.itertuples()):
    od = f"{L}/ld/{METHOD_DIR[loc.method]}/{loc.phenotype}/{loc.locus_id}"
    if os.path.exists(f"{od}/summary.tsv"):
        continue
    os.makedirs(od, exist_ok=True)
    sel = np.flatnonzero((pos >= loc.start) & (pos <= loc.end)); assert len(sel) == loc.requested_snps
    eta = float(bp.loc[loc.phenotype, "eta"]) if loc.method == "BOLT-LMM" else float(eta_reg[loc.phenotype])
    Wm = UtX[:, sel]
    Rs = Wm.T @ (Wm / (1 + eta * lam)[:, None])
    Rs = (Rs + Rs.T) / 2
    n_eff = float(np.diag(Rs).mean()); dsq = np.sqrt(np.diag(Rs)); R = Rs / dsq[:, None] / dsq[None, :]
    np.fill_diagonal(R, 1.0)
    R.astype(np.float64).tofile(f"{od}/R.bin")
    bk.iloc[sel][["SNP", "pos", "A1", "A2"]].assign(ld_index=np.arange(1, len(sel) + 1)).rename(
        columns={"A1": "counted", "A2": "other"})[["ld_index", "SNP", "pos", "counted", "other"]].to_csv(f"{od}/snps.tsv", sep="\t", index=False)
    pd.DataFrame([dict(job_id=loc.job_id, method=loc.method, phenotype=loc.phenotype, chr=k, n_snps=len(sel), n=n, M_loco=M_loco,
                       eta=eta, h2_equiv=eta / (1 + eta), eta_source="BOLT_1/delta" if loc.method == "BOLT-LMM" else "REGENIE_loco_fit",
                       n_eff=n_eff, min_diag_Rstar=float(np.diag(Rs).min()), max_diag_Rstar=float(np.diag(Rs).max()))]
                 ).to_csv(f"{od}/summary.tsv.tmp", sep="\t", index=False)
    os.replace(f"{od}/summary.tsv.tmp", f"{od}/summary.tsv")
    if i % 200 == 0: log(f"locus {i + 1}/{len(man)} {loc.job_id} m={len(sel)} n_eff={n_eff:.1f}")
open(f"{L}/status/L2_chr{k}.DONE", "w").write(f"{len(man)} loci\n")
log(f"DONE {len(man)} loci")
