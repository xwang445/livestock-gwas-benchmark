#!/usr/bin/env python3
"""L5: post hoc evaluation of the LOCO SuSiE-adj fits (run after L4).

Region truth:  eligible_regions.csv = the 2-Mb polygenic-effect windows, 0-based half-open [start, end).
               Null window = gap between consecutive effect windows on a chromosome (exactly 4 Mb), or from the
               last effect window to the chromosome end (chr1_18_lengths.csv of the simulation).
               Central null = [null_start + 1 Mb, null_start + 3 Mb), i.e. >= 1 Mb from both neighbouring effect windows.
               500-kb bins: cycle position (pos - effect_start_of_cycle) // 500 kb: 0-3 effect, 4-11 null
               (4-5 and 10-11 null flank, 6-9 central null).
Causal truth:  snp.info.chr1_18.10k.csv (10,000 SNP IDs; the same causal set for every phenotype; per-SNP effects
               were drawn implicitly by `mph --simulate` from a GRM of these SNPs and are not stored).
Loci are classified by the position of their lead SNP (never by GWAS significance or PIP).
Primary evaluation set = loci whose primary fit (estimate_residual_variance = TRUE) succeeded and converged.
"""
import glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from scipy.stats import mannwhitneyu, fisher_exact

N_PERM, SEED = 10000, 20260930
SIMDIR = os.environ["SIMULATION_DIR"]
RES, EV = f"{L}/results", f"{L}/evaluation"
os.makedirs(RES, exist_ok=True); os.makedirs(EV, exist_ok=True)
RT_ORDER = ["effect", "null_flank", "null_central"]

# ------------------------------------------------------------------ region annotation
reg = pd.read_csv(REGIONS); reg[["start", "end"]] = reg[["start", "end"]].astype(int)
reg = reg.sort_values(["chr", "start"]).reset_index(drop=True)
chr_len = pd.read_csv(f"{SIMDIR}/chr1_18_lengths.csv").set_index("chr").length


def annotate(chrom, pos):
    chrom, pos = np.asarray(chrom), np.asarray(pos)
    rtype = np.empty(len(pos), dtype=object); window = np.empty(len(pos), dtype=object); cyc = np.full(len(pos), -1)
    for c in np.unique(chrom):
        r = reg[reg.chr == c]; st, en = r.start.to_numpy(), r.end.to_numpy()
        i = np.flatnonzero(chrom == c); p = pos[i]
        j = np.searchsorted(st, p, side="right") - 1
        assert (j >= 0).all(), "SNP before first effect window"
        eff = p < en[j]
        off = p - en[j]
        rtype[i] = np.where(eff, "effect", np.where((off >= 1_000_000) & (off < 3_000_000), "null_central", "null_flank"))
        window[i] = [f"chr{c}:{st[a]}" for a in j]            # cycle id = its effect window
        cyc[i] = (p - st[j]) // 500_000
    assert (cyc >= 0).all() and (cyc <= 11).all()
    return rtype, window, cyc


# ------------------------------------------------------------------ aggregate fits
def cat(pattern):
    fs = sorted(glob.glob(f"{L}/susie_results/*/pheno*/*/{pattern}"))
    fs = [f for f in fs if os.path.getsize(f) > 1]
    return pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True) if fs else pd.DataFrame()

man = pd.read_csv(MANIFEST, sep="\t")
fits = cat("fit_summary.tsv"); snp = cat("snp_pip.tsv"); cs = cat("credible_sets.tsv"); sens = cat("sens_ERV_FALSE_fit_summary.tsv")
fits = man.merge(fits.drop(columns=["method", "phenotype", "chr"]), left_on="job_id", right_on="locus_id", how="left")
fits["fit_status"] = fits.fit_status.fillna("not_run")
fits["converged"] = fits.converged.fillna(False).astype(bool)
fits["region_type"], fits["cycle_window"], fits["lead_cycle_bin"] = annotate(fits.chr, fits.lead_position)
fits["region_class"] = np.where(fits.region_type == "effect", "effect", "null")
fits["primary"] = (fits.fit_status == "ok") & fits.converged
snp["region_type"], snp["cycle_window"], snp["cycle_bin"] = annotate(snp.chr, snp.position)
snp.to_csv(f"{RES}/snp_level_all.tsv.gz", sep="\t", index=False, compression="gzip")
cs.to_csv(f"{RES}/credible_sets_all.tsv", sep="\t", index=False)
cs = cs.merge(fits[["job_id", "method", "phenotype", "region_type", "primary"]], left_on="locus_id", right_on="job_id")

# QC: region rule == modular simulation rule (pos %% 6e6 < 2e6)
bim = read_bim(); rt_b, _, _ = annotate(bim.chr, bim.pos)
assert ((rt_b == "effect") == (bim.pos % 6_000_000 < 2_000_000)).all()
causal = pd.read_csv(CAUSAL); assert list(causal.columns) == ["SNP"] and causal.SNP.is_unique and len(causal) == 10000
assert causal.SNP.isin(bim.SNP).all()
assert (rt_b[bim.SNP.isin(causal.SNP).to_numpy()] == "effect").all(), "causal SNP outside effect windows"
cset = set(causal.SNP)
snp["causal"] = snp.SNP.isin(cset)

# ------------------------------------------------------------------ per-locus diagnostics
g = snp.groupby("locus_id")
pl = g.agg(n_SNPs_analyzed=("PIP", "size"), max_PIP=("PIP", "max"), total_PIP=("PIP", "sum"),
           n_PIP_gt_0_1=("PIP", lambda x: int((x > 0.1).sum())), n_PIP_gt_0_5=("PIP", lambda x: int((x > 0.5).sum())),
           n_causal_in_locus=("causal", "sum"), n_effect_window_SNPs=("region_type", lambda x: int((x == "effect").sum())))
pl["top1_share"] = pl.max_PIP / pl.total_PIP
pl["eff_n_SNPs"] = g.PIP.apply(lambda x: 1 / ((x / x.sum()) ** 2).sum() if x.sum() > 0 else np.nan)
pl["PIP_on_effect_window_SNPs"] = snp.assign(w=snp.PIP * (snp.region_type == "effect")).groupby("locus_id").w.sum()
csz = cs.groupby("locus_id").cs_size.agg(median_CS_size="median", total_CS_SNPs="sum")
D = fits.merge(pl, left_on="job_id", right_index=True, how="left").merge(csz, left_on="job_id", right_index=True, how="left")
D["failure_reason"] = np.where(D.fit_status == "ok", np.where(D.converged, "", "not_converged"),
                               D.failure_category.fillna("") + ": " + D.message.fillna(D.fit_status))
diag_cols = ["phenotype", "method", "chr", "job_id", "start", "end", "lead_SNP", "lead_position", "lead_logP", "region_type",
             "requested_snps", "n_SNPs_analyzed", "fit_status", "converged", "failure_category", "failure_reason", "n_cs", "max_PIP",
             "total_PIP", "top1_share", "eff_n_SNPs", "n_PIP_gt_0_1", "n_PIP_gt_0_5", "median_CS_size", "sigma2", "niter", "n_eff",
             "eta", "eta_source", "s_rss", "max_abs_z", "n_causal_in_locus", "n_effect_window_SNPs", "PIP_on_effect_window_SNPs", "warnings"]
D = D.rename(columns={"n_cs": "n_CS"}); diag_cols[diag_cols.index("n_cs")] = "n_CS"
D[diag_cols].rename(columns={"job_id": "locus"}).to_csv(f"{RES}/per_locus_diagnostics.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ A. model fitting
def qd(x, pre):
    x = pd.Series(x).dropna()
    return {f"{pre}_median": x.median(), f"{pre}_q25": x.quantile(.25), f"{pre}_q75": x.quantile(.75),
            f"{pre}_min": x.min(), f"{pre}_max": x.max()} if len(x) else {}

def fit_summary(d, cs_d):
    ok = d.fit_status == "ok"; err = d.fit_status == "error"; P = d[d.primary]
    out = dict(n_loci_attempted=len(d), n_fit_ok=int(ok.sum()), n_converged=int(d.primary.sum()),
               n_not_converged=int((ok & ~d.converged).sum()), n_failed=int(err.sum()),
               n_other_status=int((~ok & ~err).sum()), failure_rate=err.mean(), nonconvergence_rate=(ok & ~d.converged).mean(),
               n_fail_residual_variance=int((err & (d.failure_category == "residual_variance_estimation")).sum()),
               n_fail_other=int((err & (d.failure_category != "residual_variance_estimation")).sum()),
               n_CS=int(P.n_CS.sum()), loci_with_CS=int((P.n_CS > 0).sum()), frac_loci_with_CS=(P.n_CS > 0).mean(),
               mean_CS_per_locus=P.n_CS.mean(), frac_max_PIP_gt_0_5=(P.max_PIP > .5).mean(), frac_max_PIP_gt_0_9=(P.max_PIP > .9).mean())
    for k_, v in [("CS_size", cs_d.cs_size), ("max_PIP", P.max_PIP), ("total_PIP", P.total_PIP), ("sigma2", P.sigma2),
                  ("n_eff", d.n_eff), ("s_rss", P.s_rss), ("niter", P.niter), ("top1_share", P.top1_share), ("eff_n_SNPs", P.eff_n_SNPs),
                  ("lead_logP", d.lead_logP), ("n_SNPs", P.n_SNPs_analyzed)]:
        out.update(qd(v, k_))
    out["CS_size_mean"] = cs_d.cs_size.mean(); out["frac_CS_size1"] = (cs_d.cs_size == 1).mean() if len(cs_d) else np.nan
    return out

rows = []
for (meth, ph), d in [((m, "all"), D[D.method == m]) for m in sorted(D.method.unique())] + \
        [((m, p), D[(D.method == m) & (D.phenotype == p)]) for m in sorted(D.method.unique()) for p in PHENOS]:
    rows.append(dict(method=meth, phenotype=ph, **fit_summary(d, cs[cs.locus_id.isin(d.job_id[d.primary])])))
pd.DataFrame(rows).to_csv(f"{RES}/A_fit_diagnostics.tsv", sep="\t", index=False)
bad = D[~D.primary][["method", "phenotype", "job_id", "chr", "lead_SNP", "lead_logP", "region_type", "fit_status", "converged",
                     "failure_category", "message", "niter", "sigma2", "n_eff", "s_rss", "max_abs_z", "warnings"]]
bad.to_csv(f"{RES}/A_failed_or_nonconverged_loci.tsv", sep="\t", index=False)
D.assign(failure_category=D.failure_category.fillna(D.fit_status)).groupby(["method", "fit_status", "failure_category"]).size() \
    .rename("n").reset_index().to_csv(f"{RES}/A_failure_categories.tsv", sep="\t", index=False)
if len(sens):
    sens.to_csv(f"{RES}/A_sensitivity_ERV_FALSE_for_failed_loci.tsv", sep="\t", index=False)
wn = D[D.warnings.notna() & (D.warnings != "")]
wn.groupby(["method", "warnings"]).size().rename("n").reset_index().to_csv(f"{RES}/A_warnings.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ B. effect vs null regions
rows, tests = [], []
for m in sorted(D.method.unique()):
    for rt, sel in [("effect", D.region_type == "effect"), ("null_all", D.region_type != "effect"),
                    ("null_flank", D.region_type == "null_flank"), ("null_central", D.region_type == "null_central")]:
        d = D[(D.method == m) & sel]
        r = fit_summary(d, cs[cs.locus_id.isin(d.job_id[d.primary])])
        P = d[d.primary]
        r.update(frac_loci_with_causal_SNP=(P.n_causal_in_locus > 0).mean(), mean_PIP_on_effect_window_SNPs=P.PIP_on_effect_window_SNPs.mean(),
                 frac_total_PIP_on_effect_window_SNPs=P.PIP_on_effect_window_SNPs.sum() / P.total_PIP.sum() if len(P) else np.nan,
                 frac_SNPs_in_effect_windows=P.n_effect_window_SNPs.sum() / P.n_SNPs_analyzed.sum() if len(P) else np.nan)
        rows.append(dict(method=m, region=rt, **r))
    P = D[(D.method == m) & D.primary]
    for a, b in [("effect", "null_all"), ("effect", "null_central")]:
        A_ = P[P.region_type == "effect"]; B_ = P[P.region_type != "effect"] if b == "null_all" else P[P.region_type == b]
        for v in ["max_PIP", "n_CS", "top1_share", "eff_n_SNPs", "total_PIP", "median_CS_size", "lead_logP", "sigma2", "s_rss"]:
            x, y = A_[v].dropna(), B_[v].dropna()
            if len(x) and len(y):
                tests.append(dict(method=m, contrast=f"{a}_vs_{b}", variable=v, n_a=len(x), n_b=len(y), median_a=x.median(), median_b=y.median(),
                                  mean_a=x.mean(), mean_b=y.mean(), MannWhitney_P=mannwhitneyu(x, y, alternative="two-sided").pvalue))
        Dm = D[D.method == m]
        Da = Dm[Dm.region_type == "effect"]; Db = Dm[Dm.region_type != "effect"] if b == "null_all" else Dm[Dm.region_type == b]
        for v, fa, fb in [("has_CS", (A_.n_CS > 0), (B_.n_CS > 0)), ("failed", Da.fit_status == "error", Db.fit_status == "error"),
                          ("not_converged_or_failed", ~Da.primary, ~Db.primary)]:
            t = [[int(fa.sum()), int((~fa).sum())], [int(fb.sum()), int((~fb).sum())]]
            tests.append(dict(method=m, contrast=f"{a}_vs_{b}", variable=f"prop_{v}", n_a=len(fa), n_b=len(fb), mean_a=fa.mean(), mean_b=fb.mean(),
                              Fisher_P=fisher_exact(t)[1] if min(len(fa), len(fb)) else np.nan))
pd.DataFrame(rows).to_csv(f"{RES}/B_effect_vs_null_by_method.tsv", sep="\t", index=False)
pd.DataFrame(tests).to_csv(f"{RES}/B_effect_vs_null_tests.tsv", sep="\t", index=False)
# stratified by GWAS strength (lead -log10P quartiles pooled over both methods)
P = D[D.primary].copy()
P["lead_logP_bin"] = pd.qcut(P.lead_logP, 4, duplicates="drop").astype(str)
P.groupby(["method", "lead_logP_bin", "region_type"]).agg(n=("job_id", "size"), median_max_PIP=("max_PIP", "median"),
      mean_n_CS=("n_CS", "mean"), frac_with_CS=("n_CS", lambda x: (x > 0).mean()), median_top1_share=("top1_share", "median"),
      median_eff_n_SNPs=("eff_n_SNPs", "median")).reset_index().to_csv(f"{RES}/B_effect_vs_null_by_leadP_quartile.tsv", sep="\t", index=False)
B_ph = []
for (m, p, rt), d in D.groupby(["method", "phenotype", "region_type"]):
    Pp = d[d.primary]
    B_ph.append(dict(method=m, phenotype=p, region=rt, n_loci=len(d), n_primary=len(Pp), failure_rate=(d.fit_status == "error").mean(),
                     n_CS=int(Pp.n_CS.sum()), frac_with_CS=(Pp.n_CS > 0).mean(), median_max_PIP=Pp.max_PIP.median()))
pd.DataFrame(B_ph).to_csv(f"{RES}/B_effect_vs_null_by_phenotype.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ C. causal localization (effect-window lead loci only) + permutations
E_ids = set(D.job_id[D.primary & (D.region_type == "effect")])
S = snp[snp.locus_id.isin(E_ids)].copy()
rng = np.random.default_rng(SEED)
loc_rows, csnp_rows = [], []
pooled = {}                                  # method -> {stat: null draws summed over loci}; also per (method, phenotype)
PSTATS = ["sum_PIP", "n_in_CS", "any_CS_hit", "n_CS_with", "top_is_causal", "n_PIP_gt_0_1", "n_PIP_gt_0_5"]
for lid, d in S.sort_values("position").groupby("locus_id", sort=True):
    meth, ph = d.GWAS_method.iloc[0], d.phenotype.iloc[0]
    pip = d.PIP.to_numpy(); c = d.causal.to_numpy(); m_ = len(pip); K = int(c.sum()); tot = pip.sum()
    labels = d.credible_set.fillna("").astype(str).to_numpy()
    names = sorted({x for lab in labels for x in lab.split(",") if x})
    memb = np.array([[nm in lab.split(",") for nm in names] for lab in labels], dtype=bool).reshape(m_, len(names))
    in_cs = memb.any(1)
    rank = pd.Series(-pip).rank(method="average").to_numpy()          # 1 = highest PIP
    top = pip == pip.max()
    r = dict(locus_id=lid, method=meth, phenotype=ph, chr=d.chr.iloc[0], n_SNPs=m_, total_PIP=tot, n_CS=len(names), n_CS_SNPs=int(in_cs.sum()),
             K_causal_present=K, A_causal_present=K > 0,
             causal_SNPs=";".join(d.SNP[c]), causal_positions=";".join(map(str, d.position[c])),
             causal_PIPs=";".join(f"{x:.4g}" for x in pip[c]))
    obs = dict(sum_PIP=pip[c].sum(), n_in_CS=int((c & in_cs).sum()), any_CS_hit=int(memb[c].any()) if K else 0,
               n_CS_with=int(memb[c].any(0).sum()) if K else 0, top_is_causal=(c & top).sum() / top.sum(),
               n_PIP_gt_0_1=int((c & (pip > .1)).sum()), n_PIP_gt_0_5=int((c & (pip > .5)).sum()))
    r.update(sum_PIP_causal=obs["sum_PIP"], frac_PIP_causal=obs["sum_PIP"] / tot if tot > 0 else np.nan,
             max_PIP_causal=pip[c].max() if K else np.nan, best_causal_rank=rank[c].min() if K else np.nan,
             top_SNP_is_causal=obs["top_is_causal"], expected_top_is_causal=K / m_,
             B_n_causal_in_CS=obs["n_in_CS"], B_any_CS_contains_causal=bool(obs["any_CS_hit"]), n_CS_containing_causal=obs["n_CS_with"],
             C_n_causal_PIP_gt_0_1=obs["n_PIP_gt_0_1"], C_n_causal_PIP_gt_0_5=obs["n_PIP_gt_0_5"],
             CS_sizes_with_causal=";".join(str(int(memb[:, j].sum())) for j in range(len(names)) if K and memb[c, j].any()),
             CS_sizes_without_causal=";".join(str(int(memb[:, j].sum())) for j in range(len(names)) if not (K and memb[c, j].any())))
    for j, i in enumerate(np.flatnonzero(c)):
        csnp_rows.append(dict(locus_id=lid, method=meth, phenotype=ph, SNP=d.SNP.iloc[i], chr=d.chr.iloc[i], position=d.position.iloc[i],
                              z=d.z.iloc[i], log10P=d.log10P.iloc[i], PIP=pip[i], PIP_rank=rank[i], in_CS=bool(in_cs[i]),
                              credible_set=labels[i] or "", n_SNPs_locus=m_, K_locus=K))
    if K:
        sel = np.argpartition(rng.random((N_PERM, m_)), K - 1, axis=1)[:, :K] if K < m_ else np.tile(np.arange(m_), (N_PERM, 1))
        ps = pip[sel]
        nd = dict(sum_PIP=ps.sum(1), n_in_CS=in_cs[sel].sum(1),
                  any_CS_hit=memb[sel].any(axis=(1, 2)).astype(int) if len(names) else np.zeros(N_PERM),
                  n_CS_with=memb[sel].any(1).sum(1) if len(names) else np.zeros(N_PERM),
                  top_is_causal=top[sel].sum(1) / top.sum(), n_PIP_gt_0_1=(ps > .1).sum(1), n_PIP_gt_0_5=(ps > .5).sum(1))
        for key in [(meth, "all"), (meth, ph)]:
            acc = pooled.setdefault(key, {s: np.zeros(N_PERM) for s in PSTATS})
            for s_ in PSTATS: acc[s_] += nd[s_]
        f = nd["sum_PIP"] / tot
        r.update(perm_mean_frac=f.mean(), perm_median_frac=np.median(f), perm_q025_frac=np.quantile(f, .025), perm_q975_frac=np.quantile(f, .975),
                 perm_mean_sum_PIP=nd["sum_PIP"].mean(), perm_P_sum_PIP=(1 + (nd["sum_PIP"] >= obs["sum_PIP"] - 1e-12).sum()) / (N_PERM + 1),
                 perm_mean_n_in_CS=nd["n_in_CS"].mean(), perm_P_n_in_CS=(1 + (nd["n_in_CS"] >= obs["n_in_CS"]).sum()) / (N_PERM + 1),
                 perm_mean_any_CS_hit=nd["any_CS_hit"].mean())
    loc_rows.append(r)
LC = pd.DataFrame(loc_rows); LC.to_csv(f"{EV}/C_effect_loci_causal_localization.tsv", sep="\t", index=False)
CS_ = pd.DataFrame(csnp_rows); CS_.to_csv(f"{EV}/C_effect_loci_causal_SNPs.tsv", sep="\t", index=False)
cs_E = cs[cs.locus_id.isin(E_ids)].copy()
cs_E["has_causal"] = cs_E.SNPs.str.split(",").apply(lambda v: any(x in cset for x in v))
cs_E.to_csv(f"{EV}/C_effect_loci_credible_sets.tsv", sep="\t", index=False)

perm_rows = []
for (meth, ph), nd in sorted(pooled.items()):
    lc = LC[(LC.method == meth) & ((LC.phenotype == ph) if ph != "all" else True) & LC.A_causal_present]
    tot = lc.total_PIP.sum()
    obs = dict(sum_PIP=lc.sum_PIP_causal.sum(), n_in_CS=lc.B_n_causal_in_CS.sum(), any_CS_hit=lc.B_any_CS_contains_causal.sum(),
               n_CS_with=lc.n_CS_containing_causal.sum(), top_is_causal=lc.top_SNP_is_causal.sum(),
               n_PIP_gt_0_1=lc.C_n_causal_PIP_gt_0_1.sum(), n_PIP_gt_0_5=lc.C_n_causal_PIP_gt_0_5.sum())
    for s_ in PSTATS:
        v = nd[s_]; o = float(obs[s_]); den = tot if s_ == "sum_PIP" else 1.0
        perm_rows.append(dict(method=meth, phenotype=ph, statistic=s_ if s_ != "sum_PIP" else "causal_PIP_fraction", n_loci=len(lc),
                              observed=o / den, perm_mean=v.mean() / den, perm_median=np.median(v) / den,
                              perm_q025=np.quantile(v, .025) / den, perm_q975=np.quantile(v, .975) / den,
                              ratio_obs_to_perm=o / v.mean() if v.mean() > 0 else np.nan,
                              perm_P_greater=(1 + (v >= o - 1e-9).sum()) / (N_PERM + 1), perm_P_less=(1 + (v <= o + 1e-9).sum()) / (N_PERM + 1)))
PERM = pd.DataFrame(perm_rows); PERM.to_csv(f"{EV}/C_permutation_pooled.tsv", sep="\t", index=False)
np.savez_compressed(f"{EV}/C_permutation_null_draws.npz", **{f"{m}__{p}__{s}": v for (m, p), d in pooled.items() for s, v in d.items()})

acc_rows = []
for m in sorted(LC.method.unique()):
    lc = LC[LC.method == m]; lk = lc[lc.A_causal_present]; cm = CS_[CS_.method == m]; ce = cs_E[cs_E.method == m]
    pp = PERM[(PERM.method == m) & (PERM.phenotype == "all")].set_index("statistic")
    acc_rows.append(dict(
        method=m, n_effect_loci_analyzed=len(lc), frac_with_causal_present=lc.A_causal_present.mean(),
        K_mean=lc.K_causal_present.mean(), K_median=lc.K_causal_present.median(), K_min=lc.K_causal_present.min(), K_max=lc.K_causal_present.max(),
        n_causal_SNP_instances=int(lc.K_causal_present.sum()), frac_SNPs_causal=lc.K_causal_present.sum() / lc.n_SNPs.sum(),
        mean_causal_SNP_PIP=cm.PIP.mean(), median_causal_SNP_PIP=cm.PIP.median(),
        mean_noncausal_SNP_PIP=(lc.total_PIP.sum() - lc.sum_PIP_causal.sum()) / (lc.n_SNPs.sum() - lc.K_causal_present.sum()),
        mean_sum_PIP_causal_per_locus=lk.sum_PIP_causal.mean(), median_sum_PIP_causal_per_locus=lk.sum_PIP_causal.median(),
        pooled_causal_PIP_fraction=lk.sum_PIP_causal.sum() / lk.total_PIP.sum(), mean_per_locus_causal_PIP_fraction=lk.frac_PIP_causal.mean(),
        perm_pooled_causal_PIP_fraction=pp.loc["causal_PIP_fraction", "perm_mean"],
        perm_pooled_q025=pp.loc["causal_PIP_fraction", "perm_q025"], perm_pooled_q975=pp.loc["causal_PIP_fraction", "perm_q975"],
        perm_pooled_P=pp.loc["causal_PIP_fraction", "perm_P_greater"],
        median_max_causal_PIP=lk.max_PIP_causal.median(), frac_max_causal_PIP_gt_0_1=(lk.max_PIP_causal > .1).mean(),
        frac_max_causal_PIP_gt_0_5=(lk.max_PIP_causal > .5).mean(),
        median_best_causal_rank=lk.best_causal_rank.median(), median_best_causal_rank_rel=(lk.best_causal_rank / lk.n_SNPs).median(),
        frac_top_SNP_causal=lk.top_SNP_is_causal.mean(), expected_frac_top_SNP_causal=lk.expected_top_is_causal.mean(),
        n_CS=len(ce), frac_CS_containing_causal=ce.has_causal.mean() if len(ce) else np.nan,
        perm_frac_CS_containing_causal=pp.loc["n_CS_with", "perm_mean"] / len(ce) if len(ce) else np.nan,
        frac_loci_with_CS_containing_causal=lk.B_any_CS_contains_causal.mean(),
        perm_frac_loci_with_CS_containing_causal=pp.loc["any_CS_hit", "perm_mean"] / len(lk),
        causal_SNP_CS_coverage=lk.B_n_causal_in_CS.sum() / lk.K_causal_present.sum(),
        perm_causal_SNP_CS_coverage=pp.loc["n_in_CS", "perm_mean"] / lk.K_causal_present.sum(),
        frac_all_SNPs_in_CS=lc.n_CS_SNPs.sum() / lc.n_SNPs.sum(),
        median_CS_size_with_causal=ce[ce.has_causal].cs_size.median(), median_CS_size_without_causal=ce[~ce.has_causal].cs_size.median(),
        n_causal_PIP_gt_0_1=int(lk.C_n_causal_PIP_gt_0_1.sum()), n_causal_PIP_gt_0_5=int(lk.C_n_causal_PIP_gt_0_5.sum()),
        frac_causal_PIP_gt_0_1=lk.C_n_causal_PIP_gt_0_1.sum() / lk.K_causal_present.sum(),
        frac_causal_PIP_gt_0_5=lk.C_n_causal_PIP_gt_0_5.sum() / lk.K_causal_present.sum(),
        frac_loci_perm_P_lt_0_05=(lk.perm_P_sum_PIP < .05).mean()))
pd.DataFrame(acc_rows).T.to_csv(f"{EV}/C_causal_localization_summary.tsv", sep="\t", header=False)

# ------------------------------------------------------------------ 500-kb bins (downstream spatial summary)
Sp = snp[snp.locus_id.isin(set(D.job_id[D.primary]))].copy()
Sp["bin_start"] = Sp.position // 500_000 * 500_000
cs_top = cs[cs.primary].copy()
top_snp = Sp.set_index(["locus_id", "SNP"]).PIP
def cs_top_pos(row):
    v = row.SNPs.split(","); p = [top_snp.get((row.locus_id, x), -1) for x in v]; return v[int(np.argmax(p))]
cs_top["top_SNP"] = cs_top.apply(cs_top_pos, axis=1)
cs_top["position"] = cs_top.top_SNP.str.split("_").str[1].astype(int); cs_top["chr"] = cs_top.top_SNP.str.split("_").str[0].astype(int)
cs_top["bin_start"] = cs_top.position // 500_000 * 500_000
leads = D[D.primary][["method", "phenotype", "chr", "lead_position"]].assign(bin_start=lambda x: x.lead_position // 500_000 * 500_000)
# all genomic bins (from the panel)
b = bim.assign(bin_start=bim.pos // 500_000 * 500_000)
b["region_type"], _, b["cycle_bin"] = annotate(b.chr, b.pos)
gb = b.groupby(["chr", "bin_start"]).agg(n_panel_SNPs=("SNP", "size"), cycle_bin=("cycle_bin", "first"),
                                         region_type=("region_type", "first"), n_causal=("SNP", lambda x: int(x.isin(cset).sum()))).reset_index()
bins = []
for (m, p), d in Sp.groupby(["GWAS_method", "phenotype"]):
    a = d.groupby(["chr", "bin_start"]).agg(n_SNPs_analyzed=("PIP", "size"), sum_PIP=("PIP", "sum"),
                                            n_PIP_gt_0_1=("PIP", lambda x: int((x > .1).sum())), n_PIP_gt_0_5=("PIP", lambda x: int((x > .5).sum())),
                                            sum_PIP_causal=("PIP", lambda x: float(x[d.loc[x.index, "causal"]].sum())))
    c_ = cs_top[(cs_top.method == m) & (cs_top.phenotype == p)].groupby(["chr", "bin_start"]).size().rename("n_CS")
    l_ = leads[(leads.method == m) & (leads.phenotype == p)].groupby(["chr", "bin_start"]).size().rename("n_leads")
    x = gb.merge(a, on=["chr", "bin_start"], how="left").merge(c_, on=["chr", "bin_start"], how="left").merge(l_, on=["chr", "bin_start"], how="left")
    x = x.fillna({"n_SNPs_analyzed": 0, "sum_PIP": 0, "n_PIP_gt_0_1": 0, "n_PIP_gt_0_5": 0, "sum_PIP_causal": 0, "n_CS": 0, "n_leads": 0})
    bins.append(x.assign(method=m, phenotype=p))
BINS = pd.concat(bins, ignore_index=True)
BINS.to_csv(f"{RES}/D_bins500kb_genome.tsv.gz", sep="\t", index=False, compression="gzip")
BINS["bin_class"] = np.where(BINS.cycle_bin < 4, "effect", np.where(BINS.cycle_bin.isin([6, 7, 8, 9]), "null_central", "null_flank"))
agg = dict(n_bins=("sum_PIP", "size"), frac_bins_covered=("n_SNPs_analyzed", lambda x: (x > 0).mean()), mean_sum_PIP=("sum_PIP", "mean"),
           mean_n_CS=("n_CS", "mean"), mean_n_leads=("n_leads", "mean"), mean_n_PIP_gt_0_5=("n_PIP_gt_0_5", "mean"),
           total_PIP=("sum_PIP", "sum"), total_CS=("n_CS", "sum"), total_leads=("n_leads", "sum"), total_panel_SNPs=("n_panel_SNPs", "sum"),
           total_analyzed_SNPs=("n_SNPs_analyzed", "sum"), total_causal=("n_causal", "sum"))
cyc = BINS.groupby(["method", "cycle_bin", "bin_class"]).agg(**agg).reset_index()
cyc["PIP_per_analyzed_SNP"] = cyc.total_PIP / cyc.total_analyzed_SNPs
cyc.to_csv(f"{RES}/D_bins500kb_by_cycle_position.tsv", sep="\t", index=False)
cls = BINS.groupby(["method", "bin_class"]).agg(**agg).reset_index()
for m in cls.method.unique():
    t = cls.method == m
    for c_ in ["total_PIP", "total_CS", "total_leads", "total_panel_SNPs", "total_analyzed_SNPs"]:
        cls.loc[t, "share_" + c_.replace("total_", "")] = cls.loc[t, c_] / cls.loc[t, c_].sum()
cls.to_csv(f"{RES}/D_bins500kb_by_class.tsv", sep="\t", index=False)
BINS.groupby(["method", "phenotype", "bin_class"]).agg(**agg).reset_index().to_csv(f"{RES}/D_bins500kb_by_class_phenotype.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ LOCO adjustment parameters
eta = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{L}/qc/regenie_eta_chr*.tsv"))])
eta.to_csv(f"{RES}/loco_eta_REGENIE_by_pheno_chr.tsv", sep="\t", index=False)
pd.read_csv(f"{L}/qc/bolt_variance_params.tsv", sep="\t").to_csv(f"{RES}/loco_eta_BOLT_by_pheno.tsv", sep="\t", index=False)

print(pd.read_csv(f"{RES}/A_fit_diagnostics.tsv", sep="\t").query("phenotype == 'all'").T.to_string())
print(pd.read_csv(f"{RES}/B_effect_vs_null_by_method.tsv", sep="\t")[["method", "region", "n_loci_attempted", "n_converged", "failure_rate",
      "n_CS", "frac_loci_with_CS", "max_PIP_median", "CS_size_median", "top1_share_median", "total_PIP_median", "lead_logP_median"]].to_string())
print(pd.read_csv(f"{EV}/C_causal_localization_summary.tsv", sep="\t", header=None).to_string())
print(PERM[PERM.phenotype == "all"].to_string())
print(cls.to_string())
