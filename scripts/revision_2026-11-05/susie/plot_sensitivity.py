#!/usr/bin/env python3
"""Plot only canonical paired sensitivity tables; no fitting or manual statistics."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/gwas_revision_susie_mpl')
from pathlib import Path
import os
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

W=Path(os.environ["SUSIE_REVISION_ROOT"]).resolve()
OUT=W/'figures';OUT.mkdir(exist_ok=True)
src=W/'results/sensitivity_comparison_by_phenotype.tsv'
diag=W/'results/fit_diagnostics_comparison.tsv'
d=pd.read_csv(src,sep='\t');d=d[d.universe=='analyzed_snps']
f=pd.read_csv(diag,sep='\t')
methods=['BOLT-LMM','REGENIE-LOCO'];modes=['ERV_FALSE','ERV_TRUE']
colors={'BOLT-LMM':'#0072B2','REGENIE-LOCO':'#D55E00'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,
    'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,2,figsize=(7.2,6.1))
for ax,letter,title in zip(axs.ravel(),'abcd',
        ['Precision–recall ranking','High-PIP causal precision','Credible-set coverage','Residual variance estimates']):
    ax.text(-.16,1.07,letter,transform=ax.transAxes,fontweight='bold',fontsize=12)
    ax.set_title(title,pad=11,fontsize=10)
    ax.grid(axis='y',lw=.5,color='#E7E7E7');ax.set_axisbelow(True)
    ax.set_xticks([0,1]);ax.set_xticklabels(['BOLT-LMM','REGENIE\n(LOCO)'])
    ax.set_xlim(-.5,1.5)
for i,m in enumerate(methods):
    for j,mode in enumerate(modes):
        q=d[(d.method==m)&(d.fit_mode==mode)].copy()
        q['number']=q.phenotype.str.replace('pheno','',regex=False).astype(int)
        q=q.sort_values('number');x=i+(j-.5)*.25
        col=colors[m];marker='o' if j==0 else '^'
        for ax,key,scale in [(axs[0,0],'average_precision',1),(axs[1,0],'cs_coverage',100)]:
            vals=q[key].to_numpy()*scale
            ax.scatter(x+np.linspace(-.035,.035,10),vals,s=13,marker=marker,
                facecolors='white' if j==0 else col,edgecolors=col,linewidths=.6,zorder=3)
            ax.plot([x-.06,x+.06],[vals.mean()]*2,color=col,lw=1.8)
        if j==0:
            b=q.random_precision_baseline.mean()
            axs[0,0].plot([i-.30,i+.30],[b,b],ls=':',lw=1.2,color='#777777')
        totals=q[['n_high_pip','n_high_pip_causal']].sum()
        val=100*totals.n_high_pip_causal/totals.n_high_pip
        axs[0,1].bar(x,val,width=.20,color='white' if j==0 else col,edgecolor=col,lw=1)
        axs[0,1].text(x,val+.4,f'{val:.2f}',ha='center',fontsize=7)
        z=f[(f.method==m)&(f.fit_mode==mode)]
        lo,med,hi=z.sigma2_min.min(),z.sigma2_median.median(),z.sigma2_max.max()
        axs[1,1].errorbar(x,med,yerr=[[med-lo],[hi-med]],fmt=marker,
            color=col,markerfacecolor='white' if j==0 else col,markersize=4,lw=1,capsize=3)
axs[0,0].set_ylabel('Average precision')
axs[0,1].set_ylabel('Causal SNPs among\nhigh-PIP calls (%)');axs[0,1].set_ylim(0,16)
axs[1,0].set_ylabel('Credible sets containing\na causal SNP (%)')
axs[1,1].set_ylabel('Residual variance');axs[1,1].set_ylim(.97,1.004)
fig.legend(handles=[Line2D([],[],marker='o',color='#444444',markerfacecolor='white',ls='',label='Fixed residual variance'),
    Line2D([],[],marker='^',color='#444444',ls='',label='Estimated residual variance'),
    Line2D([],[],color='#777777',ls=':',label='Random-ranking precision baseline')],
    loc='lower center',bbox_to_anchor=(.5,.015),ncol=1,frameon=False,fontsize=8)
fig.subplots_adjust(left=.105,right=.98,bottom=.23,top=.93,wspace=.47,hspace=.58)
for ext in ['png','pdf','svg','tiff']:
    kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
    fig.savefig(OUT/f'SuSiE_residual_variance_sensitivity.{ext}',dpi=600,**kw)
plt.close(fig)
(OUT/'SuSiE_residual_variance_sensitivity_caption.txt').write_text(
    'Residual-variance sensitivity of standard SuSiE-RSS. '
    '(a) Average precision among unique analyzed SNPs; dotted segments indicate the mean causal prevalence '
    'in each GWAS-specific evaluation pool. (b) Pooled precision of SNP calls with PIP ≥0.9 across ten phenotypes. '
    '(c) Percentage of nominal 95% credible sets containing at least one causal SNP. Points in a and c '
    'represent ten phenotypes, and horizontal segments represent means. (d) Median of the ten phenotype-level '
    'median residual variances, with the full locus-level minimum–maximum range. Fixed variance is one. '
    'Both runs use identical loci, GWAS statistics, allele alignment, ordinary in-sample LD, N=20,000, L=10, '
    '95% credible sets, minimum absolute correlation 0.5, maximum 1,000 iterations, and seed 20261014. '
    'Only estimate_residual_variance changes. All 6,732 BOLT-LMM and 6,569 REGENIE-LOCO fits converged '
    'under both settings without warnings. Each SNP is counted once per phenotype and GWAS method, using '
    'its maximum PIP across overlapping loci. These are descriptive comparisons; correlated calls are '
    'not treated as independent replicates.\n')
(OUT/'source_manifest.json').write_text(json.dumps({str(p.relative_to(W)):
    hashlib.sha256(p.read_bytes()).hexdigest() for p in [src,diag]},indent=2)+'\n')
