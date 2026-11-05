#!/usr/bin/env python3
"""Plot empirical causal precision and retained calls across saved PIP cutoffs."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/susie_precision_threshold_mpl')
from pathlib import Path
import os
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

W = Path(os.environ["SUSIE_REVISION_ROOT"]).resolve()
B = Path(os.environ["GWAS_BENCHMARK_ROOT"]).resolve()
OUT = W / 'publication/precision_threshold'
OUT.mkdir(parents=True, exist_ok=True)
ROOTS = {'SuSiE-RSS': W / 'results',
         'SuSiE-adj': B / 'SuSiE_adj_finemapping/loco/results'}
TRUTH_PATH = B / 'snp.info.chr1_18.10k.csv'
truth = set(pd.read_csv(TRUTH_PATH).SNP)
assert len(truth) == 10000
# Extra resolution in the high-confidence tail; avoid an exact-PIP=1 endpoint
# whose inclusion depends on numerical saturation/serialization.
thresholds = np.unique(np.concatenate([np.arange(90) / 100,
                                      np.arange(900, 1000) / 1000,
                                      [.9995, .9999]]))
sources = [TRUTH_PATH]
rows, qc = [], []
pool_sets = {}
for analysis, root in ROOTS.items():
    sp, dp = root / 'snp_level_all.tsv.gz', root / 'per_locus_diagnostics.tsv'
    sources.extend([sp, dp])
    diag = pd.read_csv(dp, sep='\t')
    idcol = 'job_id' if 'job_id' in diag else ('locus_id' if 'locus_id' in diag else 'locus')
    good = diag.loc[diag.fit_status.eq('ok') & diag.converged.eq(True), idcol]
    assert len(good) == len(diag) == 13301
    s = pd.read_csv(sp, sep='\t', usecols=['phenotype', 'GWAS_method', 'SNP', 'PIP', 'locus_id'])
    assert s.locus_id.isin(good).all()
    assert s.PIP.notna().all() and s.PIP.between(0, 1).all()
    for (method, pheno), q in s.groupby(['GWAS_method', 'phenotype']):
        original = len(q)
        q = q.sort_values('PIP', ascending=False).drop_duplicates('SNP')
        pool_sets[(analysis, method, pheno)] = set(q.SNP)
        pip = q.PIP.to_numpy()
        causal = q.SNP.isin(truth).to_numpy()
        qc.append(dict(fine_mapping=analysis, GWAS_method=method, phenotype=pheno,
                       n_pool=len(q), n_causal_pool=int(causal.sum()),
                       duplicates_removed=original-len(q)))
        causal_sum = np.concatenate([[0], np.cumsum(causal)])
        pip_sum = np.concatenate([[0], np.cumsum(pip)])
        # Descending PIP; count all ties at the cutoff.
        counts = np.searchsorted(-pip, -thresholds, side='right')
        for t, n in zip(thresholds, counts):
            rows.append(dict(fine_mapping=analysis, GWAS_method=method, phenotype=pheno,
                             threshold=t, n_calls=int(n), n_causal=int(causal_sum[n]),
                             sum_PIP=float(pip_sum[n])))
    print('Read and checked', analysis, flush=True)
for (analysis, method, pheno), ids in pool_sets.items():
    assert ids == pool_sets[('SuSiE-adj', method, pheno)]
per = pd.DataFrame(rows)
keys = ['GWAS_method', 'fine_mapping', 'threshold']
pooled = per.groupby(keys, as_index=False)[['n_calls', 'n_causal', 'sum_PIP']].sum()
for frame in [per, pooled]:
    den = frame.n_calls.replace(0, np.nan)
    frame['precision_pct'] = 100 * frame.n_causal / den
    frame['mean_PIP_pct'] = 100 * frame.sum_PIP / den
    frame['noncausal_pct'] = 100 - frame.precision_pct
for _, q in pooled.groupby(['GWAS_method', 'fine_mapping']):
    q = q.sort_values('threshold')
    assert (np.diff(q.n_calls) <= 0).all()
    assert (np.diff(q.n_causal) <= 0).all()
    assert (q.n_causal <= q.n_calls).all()
    assert (q.mean_PIP_pct >= 100*q.threshold-1e-7).all()
# Independent reconciliation against the existing publication table at 0.90.
reference = pd.read_csv(W / 'publication/localization_precision/Table_3_numeric.csv')
check = pooled[pooled.threshold.eq(.9)].merge(
    reference, on=['GWAS_method', 'fine_mapping'], suffixes=('_new', '_old'), validate='one_to_one')
assert len(check) == 4
for col in ['n_calls', 'n_causal', 'precision_pct']:
    assert np.allclose(check[col+'_new'], check[col+'_old'])
assert np.allclose(check.mean_PIP_pct, check.mean_assigned_PIP_pct)
per.to_csv(OUT / 'precision_by_threshold_phenotype.csv', index=False)
pooled.to_csv(OUT / 'precision_by_threshold_pooled.csv', index=False)
pd.DataFrame(qc).to_csv(OUT / 'pool_qc.csv', index=False)
selected = pooled[pooled.threshold.isin([.5, .9, .95, .99, .999, .9999])]
selected.to_csv(OUT / 'selected_thresholds.csv', index=False)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.labelsize': 10, 'axes.titlesize': 11,
                     'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'svg.fonttype': 'none'})
colors = {'BOLT-LMM': '#007A91', 'REGENIE-LOCO': '#BC572A'}
styles = {'SuSiE-RSS': '-', 'SuSiE-adj': (0, (4, 2))}
order = [(m, a) for m in colors for a in styles]
fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.8), gridspec_kw={'width_ratios': [1.15, 1, 1]})
ax, zoom, count = axes
ax.plot([0, 1], [0, 100], color='#8B9298', lw=1.3, ls=(0, (2, 3)), zorder=1)
ax.text(.43, 48, 'PIP threshold reference', color='#70777D', rotation=39,
        fontsize=8.5, rotation_mode='anchor')
legend = []
for method, analysis in order:
    q = pooled[(pooled.GWAS_method == method) & (pooled.fine_mapping == analysis)]
    q = q.sort_values('threshold')
    kw = dict(color=colors[method], ls=styles[analysis], lw=1.8)
    ax.plot(q.threshold, q.precision_pct, **kw)
    high = q[q.threshold >= .9]
    zoom.plot(high.threshold, high.precision_pct, **kw)
    count.plot(high.threshold, high.n_calls, **kw)
    marked = q[q.threshold.isin([.9, .95, .99])]
    marker = 'o' if analysis == 'SuSiE-RSS' else 's'
    for target, ycol in [(zoom, 'precision_pct'), (count, 'n_calls')]:
        target.plot(marked.threshold, marked[ycol], ls='none', marker=marker,
                    ms=4, markerfacecolor='white', markeredgecolor=colors[method], zorder=4)
    legend.append(Line2D([], [], **kw, label=method+' / '+analysis))
ax.set(xlim=(0, 1), ylim=(0, 103), xlabel='PIP threshold, t',
       ylabel='Causal precision (%)', title='a   Full threshold range')
ax.set_xticks(np.arange(0, 1.01, .2))
ax.set_yticks(np.arange(0, 101, 20))
ax.text(.035, .96, f'At PIP ≥ 0.90\nMean assigned PIP: {reference.mean_assigned_PIP_pct.min():.1f}–{reference.mean_assigned_PIP_pct.max():.1f}%\nObserved precision: {reference.precision_pct.min():.1f}–{reference.precision_pct.max():.1f}%',
        transform=ax.transAxes, va='top', fontsize=9, linespacing=1.5,
        bbox=dict(boxstyle='round,pad=.45', facecolor='#F1F4F5', edgecolor='none'))
zoom.set(xlim=(.898, 1.001), ylim=(0, 20), xlabel='PIP threshold, t',
         ylabel='Causal precision (%)', title='b   High-PIP detail')
zoom.set_yticks(np.arange(0, 21, 5))
count.set(xlim=(.898, 1.001), ylim=(0, 2300), xlabel='PIP threshold, t',
          ylabel='Retained SNP–phenotype calls', title='c   High-PIP call counts')
count.set_yticks(np.arange(0, 2001, 500))
for target in [zoom, count]:
    target.set_xticks([.90, .95, .99])
    target.set_xticklabels(['0.90', '0.95', '0.99'])
for target in axes:
    target.grid(axis='y', color='#E5E8EA', lw=.6)
    target.set_axisbelow(True)
    target.tick_params(length=3)
fig.suptitle('Causal precision across PIP thresholds', x=.065, ha='left', y=.99,
             fontsize=16, weight='bold')
fig.text(.065, .92, 'LOCO GWAS inputs • diffuse polygenic simulation • pooled across 10 phenotypes',
         fontsize=10, color='#4D5963')
fig.legend(handles=legend, loc='lower center', bbox_to_anchor=(.52, .04), ncol=2,
           frameon=False, fontsize=9, handlelength=3.5, columnspacing=2.5)
fig.text(.065, .015, 'Precision = causal calls / all calls with PIP ≥ t. Curves are descriptive; no independence-based confidence intervals.',
         fontsize=8, color='#59636B')
fig.subplots_adjust(left=.065, right=.985, bottom=.26, top=.80, wspace=.34)
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(OUT / ('precision_vs_threshold.'+ext), dpi=300, facecolor='white')
fig.savefig(OUT / 'precision_vs_threshold.tiff', dpi=600, facecolor='white',
            pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig)

caption = 'Causal precision across PIP thresholds. (a) Pooled precision over the full threshold range. The gray diagonal is the PIP threshold expressed as a percentage, a reference rather than observed precision. (b) Precision for thresholds from 0.90 to 0.9999. (c) Retained calls over the same high-PIP range. Markers show thresholds 0.90, 0.95 and 0.99. Colors identify BOLT-LMM and REGENIE-LOCO inputs; solid and dashed lines identify standard SuSiE-RSS and LOCO SuSiE-adj. Precision is causal calls divided by all calls with PIP at least the threshold. Counts pool ten phenotypes, with each SNP counted once per phenotype, GWAS method and pipeline, retaining its maximum PIP across overlapping loci. All 13,301 fits converged in both pipelines. Truth comprises the 10,000 simulated causal SNPs and enters only the evaluation. Both pipelines estimate residual variance. Standard RSS uses ordinary in-sample LD and N=20,000; SuSiE-adj uses LOCO-adjusted LD and locus-specific effective sample size. Other settings are L=10, 95% credible sets, minimum absolute correlation 0.5 and maximum 1,000 iterations. Thus the comparison changes LD and sample-size input together. Thresholds use 0.01 increments below 0.90, 0.001 increments from 0.90 to 0.999, plus 0.9995 and 0.9999. Exact PIP=1 is excluded because numerical saturation affects that endpoint. Curves are descriptive and do not assume independent SNP calls.\n'
(OUT / 'figure_caption.txt').write_text(caption)
def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()
(OUT / 'provenance.json').write_text(json.dumps(dict(
    script=str(Path(__file__).resolve()), thresholds=thresholds.tolist(),
    sources_sha256={str(p): digest(p) for p in sources},
    aggregation='ratio of sums across phenotypes; one call per SNP per method/pipeline/phenotype',
    validation='all fits converged; same SNP pools; threshold monotonicity of counts; PIP 0.90 matches existing Table 3',
    uncertainty='descriptive; no independence-based confidence intervals'), indent=2)+'\n')
print(selected.to_string(index=False))
print('Saved figure and source tables to', OUT)
