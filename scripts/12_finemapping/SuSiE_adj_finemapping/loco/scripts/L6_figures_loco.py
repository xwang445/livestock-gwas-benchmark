#!/usr/bin/env python3
"""L6: figures for the LOCO SuSiE-adj analysis (reads L5 outputs). PDF + PNG (300 dpi)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES, EV, FIG = f"{L}/results", f"{L}/evaluation", f"{L}/figures"
os.makedirs(FIG, exist_ok=True)
COL = {"BOLT-LMM": "#2a78d6", "REGENIE-LOCO": "#eb6834"}
GRAY, INK, INK2 = "#9a9993", "#0b0b0b", "#52514e"
METHODS = ["BOLT-LMM", "REGENIE-LOCO"]
RT = ["effect", "null_flank", "null_central"]
RT_LAB = {"effect": "Effect\n(2 Mb)", "null_flank": "Null\nflank", "null_central": "Central\nnull"}
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "axes.edgecolor": INK2, "axes.linewidth": 0.6,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "figure.dpi": 150, "savefig.dpi": 300, "font.family": "DejaVu Sans"})

def save(fig, name):
    fig.savefig(f"{FIG}/{name}.pdf", bbox_inches="tight"); fig.savefig(f"{FIG}/{name}.png", bbox_inches="tight"); plt.close(fig)

def tag(ax, s): ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")

def grouped_box(ax, D, var, ylabel, log=False):
    for j, m in enumerate(METHODS):
        data = [D[(D.method == m) & (D.region_type == r)][var].dropna().to_numpy() for r in RT]
        pos = np.arange(len(RT)) + (j - 0.5) * 0.36
        bp = ax.boxplot(data, positions=pos, widths=0.3, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1))
        for b in bp["boxes"]: b.set(facecolor=COL[m], alpha=0.55, edgecolor=COL[m], lw=0.8)
        for w in bp["whiskers"] + bp["caps"]: w.set(color=COL[m], lw=0.8)
    ax.set_xticks(range(len(RT))); ax.set_xticklabels([RT_LAB[r] for r in RT]); ax.set_ylabel(ylabel)
    if log: ax.set_yscale("log")

D = pd.read_csv(f"{RES}/per_locus_diagnostics.tsv", sep="\t").rename(columns={"locus": "job_id"})
P = D[(D.fit_status == "ok") & (D.converged == True)]
A = pd.read_csv(f"{RES}/A_fit_diagnostics.tsv", sep="\t")

# ---------------------------------------------------------------- Fig 1: model fitting
fig, ax = plt.subplots(1, 4, figsize=(11, 2.7), gridspec_kw=dict(wspace=0.45))
cats = [("converged", lambda d: (d.fit_status == "ok") & (d.converged == True)),
        ("not\nconverged", lambda d: (d.fit_status == "ok") & (d.converged != True)),
        ("σ² est.\nfailure", lambda d: (d.fit_status == "error") & (d.failure_category == "residual_variance_estimation")),
        ("other\nfailure", lambda d: ~d.fit_status.isin(["ok"]) & (d.failure_category != "residual_variance_estimation"))]
for j, m in enumerate(METHODS):
    d = D[D.method == m]
    for i, (lab, f) in enumerate(cats):
        v = int(f(d).sum())
        ax[0].bar(i + (j - 0.5) * 0.38, v, 0.34, color=COL[m], label=m if i == 0 else None)
        ax[0].text(i + (j - 0.5) * 0.38, v * 1.15 if v else 0.3, f"{v:,}", ha="center", va="bottom", fontsize=5.5, color=INK2, rotation=90 if v > 999 else 0)
ax[0].set_xticks(range(4)); ax[0].set_xticklabels([c[0] for c in cats], fontsize=6.5); ax[0].set_ylabel("Loci"); ax[0].set_yscale("symlog", linthresh=10); ax[0].set_ylim(0, 1e5)
ax[0].legend(loc="upper right", fontsize=7); ax[0].set_title("Fit status (ERV = TRUE)"); tag(ax[0], "a")
for k, (var, lab, lg) in enumerate([("sigma2", "Estimated residual variance σ²", False), ("n_eff", "LOCO-adjusted n_eff", False),
                                    ("s_rss", "z/LD inconsistency s (estimate_s_rss)", True)]):
    a = ax[k + 1]
    src = P if var != "n_eff" else D
    for m in METHODS:
        x = src[src.method == m][var].dropna()
        if lg: x = np.log10(x.clip(lower=1e-9))
        a.hist(x, bins=60, color=COL[m], alpha=0.6, label=m)
    a.set_xlabel(("log10 " if lg else "") + lab); a.set_ylabel("Loci"); tag(a, "bcd"[k])
ax[1].set_title("Residual variance"); ax[2].set_title("Effective sample size"); ax[3].set_title("z/LD consistency")
save(fig, "fig1_model_fitting")

# ---------------------------------------------------------------- Fig 2: effect vs null
fig, ax = plt.subplots(2, 3, figsize=(10, 5.6), gridspec_kw=dict(wspace=0.38, hspace=0.55))
ax = ax.ravel()
grouped_box(ax[0], P, "max_PIP", "Maximum PIP per locus"); ax[0].set_title("Maximum PIP"); tag(ax[0], "a")
for j, m in enumerate(METHODS):
    fr = [(P[(P.method == m) & (P.region_type == r)].n_CS > 0).mean() for r in RT]
    mn = [P[(P.method == m) & (P.region_type == r)].n_CS.mean() for r in RT]
    ax[1].bar(np.arange(3) + (j - 0.5) * 0.38, mn, 0.34, color=COL[m], label=m)
    for i, v in enumerate(mn): ax[1].text(i + (j - 0.5) * 0.38, v, f"{v:.2f}", ha="center", va="bottom", fontsize=6, color=INK2)
    ax[2].bar(np.arange(3) + (j - 0.5) * 0.38, fr, 0.34, color=COL[m])
    for i, v in enumerate(fr): ax[2].text(i + (j - 0.5) * 0.38, v, f"{v:.2f}", ha="center", va="bottom", fontsize=6, color=INK2)
for a, t, yl in [(ax[1], "Credible sets per locus", "Mean # 95% CS"), (ax[2], "Loci with ≥1 credible set", "Fraction of loci")]:
    a.set_xticks(range(3)); a.set_xticklabels([RT_LAB[r] for r in RT]); a.set_title(t); a.set_ylabel(yl)
ax[1].legend(fontsize=7, loc="upper left", bbox_to_anchor=(0, 1.0), ncol=2); ax[1].set_ylim(0, ax[1].get_ylim()[1] * 1.2); tag(ax[1], "b"); tag(ax[2], "c")
cs = pd.read_csv(f"{RES}/credible_sets_all.tsv", sep="\t").merge(P[["job_id", "method", "region_type"]], left_on="locus_id", right_on="job_id")
grouped_box(ax[3], cs, "cs_size", "SNPs per 95% CS", log=True); ax[3].set_title("Credible-set size"); tag(ax[3], "d")
grouped_box(ax[4], P, "top1_share", "Max PIP / total locus PIP"); ax[4].set_title("PIP concentration"); tag(ax[4], "e")
q = pd.read_csv(f"{RES}/B_effect_vs_null_by_leadP_quartile.tsv", sep="\t")
qs = sorted(q.lead_logP_bin.unique(), key=lambda s: float(s.strip("([]").split(",")[0]))
for m in METHODS:
    for r, ls in [("effect", "-"), ("null_central", "--")]:
        z = q[(q.method == m) & (q.region_type == r)].set_index("lead_logP_bin").reindex(qs)
        ax[5].plot(range(len(qs)), z.median_max_PIP, ls, color=COL[m], lw=2, marker="o", ms=4,
                   label=f"{m.split('-')[0]}, {'effect' if r == 'effect' else 'central null'}")
ax[5].set_xticks(range(len(qs))); ax[5].set_xticklabels([f"Q{i + 1}" for i in range(len(qs))])
ax[5].set_xlabel("Lead −log10 P quartile"); ax[5].set_ylabel("Median max PIP"); ax[5].set_title("Controlling for GWAS strength")
ax[5].legend(fontsize=6); tag(ax[5], "f")
save(fig, "fig2_effect_vs_null")

# ---------------------------------------------------------------- Fig 3: causal localization + permutation
LC = pd.read_csv(f"{EV}/C_effect_loci_causal_localization.tsv", sep="\t"); LC = LC[LC.A_causal_present]
PERM = pd.read_csv(f"{EV}/C_permutation_pooled.tsv", sep="\t"); PERM = PERM[PERM.phenotype == "all"].set_index(["method", "statistic"])
nd = np.load(f"{EV}/C_permutation_null_draws.npz")
fig, ax = plt.subplots(1, 4, figsize=(12, 2.9), gridspec_kw=dict(wspace=0.42))
for m in METHODS:
    d = LC[LC.method == m]
    ax[0].scatter(d.perm_mean_frac, d.frac_PIP_causal, s=4, alpha=0.35, color=COL[m], lw=0, label=m)
lim = [0, max(1e-3, LC.frac_PIP_causal.max(), LC.perm_mean_frac.max()) * 1.02]
ax[0].plot(lim, lim, color=INK2, lw=0.8, ls="--"); ax[0].set_xlim(lim); ax[0].set_ylim(lim)
ax[0].set_xlabel("Matched-permutation mean\n(K random SNPs, same locus)"); ax[0].set_ylabel("Observed fraction of locus PIP\non true causal SNPs")
ax[0].set_title("Per locus"); ax[0].legend(fontsize=6, markerscale=3); tag(ax[0], "a")
for m in METHODS:
    tot = LC[LC.method == m].total_PIP.sum(); v = nd[f"{m}__all__sum_PIP"] / tot; o = PERM.loc[(m, "causal_PIP_fraction"), "observed"]
    ax[1].hist(v, bins=60, color=COL[m], alpha=0.35); ax[1].axvline(o, color=COL[m], lw=2)
    pv = PERM.loc[(m, 'causal_PIP_fraction'), 'perm_P_greater']; pv = "P<1e-4" if pv <= 1.5 / 10001 else f"P={pv:.2g}"
    ax[1].text(0.93, 0.55 - 0.1 * METHODS.index(m), f"{m.split('-')[0]}: obs {o:.3f}\n{pv}" if False else f"{m.split('-')[0]}: {o:.3f}, {pv}",
               transform=ax[1].transAxes, ha="right", va="top", color=INK, fontsize=6)
    ax[1].plot([], [], color=COL[m], lw=2)
ax[1].set_xlabel("Pooled causal-PIP fraction"); ax[1].set_ylabel("Permutation draws"); ax[1].set_title("Pooled vs matched null (10,000)"); tag(ax[1], "b")
stats = [("any_CS_hit", "Locus CS\nhas causal", "n_loci"), ("n_in_CS", "Causal\nin CS", "K"),
         ("n_PIP_gt_0_1", "Causal\nPIP>0.1", "K"), ("top_is_causal", "Top SNP\ncausal", "n_loci")]
for j, m in enumerate(METHODS):
    d = LC[LC.method == m]; den = {"n_loci": len(d), "K": d.K_causal_present.sum()}
    for i, (s, lab, dn) in enumerate(stats):
        o = PERM.loc[(m, s), "observed"] / den[dn]; e = PERM.loc[(m, s), "perm_mean"] / den[dn]
        lo_, hi_ = PERM.loc[(m, s), "perm_q025"] / den[dn], PERM.loc[(m, s), "perm_q975"] / den[dn]
        x = i + (j - 0.5) * 0.38
        ax[2].bar(x, o, 0.34, color=COL[m], label=m if i == 0 else None)
        ax[2].errorbar(x, e, yerr=[[e - lo_], [hi_ - e]], fmt="_", color=INK, ms=10, mew=1.5, lw=0.8, label="perm. mean (95%)" if i == 0 and j == 0 else None)
ax[2].set_xticks(range(len(stats))); ax[2].set_xticklabels([s[1] for s in stats], fontsize=6.5); ax[2].set_ylabel("Fraction")
ax[2].set_title("Observed vs matched null"); ax[2].legend(fontsize=6); tag(ax[2], "c")
for m in METHODS:
    ax[3].hist(LC[LC.method == m].perm_P_sum_PIP, bins=np.linspace(0, 1, 21), color=COL[m], alpha=0.5, label=m)
ax[3].legend(fontsize=6); ax[3].set_xlabel("Per-locus permutation P (causal PIP)"); ax[3].set_ylabel("Loci"); ax[3].set_title("Per-locus enrichment P"); tag(ax[3], "d")
save(fig, "fig3_causal_localization_permutation")

# ---------------------------------------------------------------- Fig 4: 500-kb cycle profile
cyc = pd.read_csv(f"{RES}/D_bins500kb_by_cycle_position.tsv", sep="\t")
fig, ax = plt.subplots(1, 3, figsize=(11, 2.7), gridspec_kw=dict(wspace=0.35))
for a, v, yl in [(ax[0], "mean_sum_PIP", "Mean Σ PIP per 500-kb bin"), (ax[1], "mean_n_CS", "Mean # CS (top SNP) per bin"),
                 (ax[2], "mean_n_leads", "Mean # GWAS lead SNPs per bin")]:
    a.axvspan(-0.5, 3.5, color="#cde2fb", alpha=0.5, lw=0); a.axvspan(5.5, 9.5, color="#f0efec", lw=0)
    for m in METHODS:
        z = cyc[cyc.method == m].sort_values("cycle_bin")
        a.plot(z.cycle_bin, z[v], color=COL[m], lw=2, marker="o", ms=4, label=m)
    a.set_xticks(range(12)); a.set_xticklabels([f"{i * 0.5:.1f}" for i in range(12)], fontsize=6)
    a.set_xlabel("Position within 6-Mb cycle (Mb)"); a.set_ylabel(yl); a.set_xlim(-0.5, 11.5)
ax[0].text(1.5, ax[0].get_ylim()[1], "effect", ha="center", va="top", fontsize=7, color=INK2)
ax[0].text(7.5, ax[0].get_ylim()[1], "central null", ha="center", va="top", fontsize=7, color=INK2)
ax[2].legend(fontsize=6, loc="upper right"); [tag(a, s) for a, s in zip(ax, "abc")]
save(fig, "fig4_bins500kb_cycle_profile")
print("figures written to", FIG)
