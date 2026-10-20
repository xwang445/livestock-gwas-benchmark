#!/usr/bin/env python3
"""Aggregate saved fits and evaluate CS coverage/power/size after fitting; no model tuning."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/standard_susie_mpl')
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B=Path(__file__).resolve().parents[2]
W=B/'SuSiE_finemapping'
OLD=B/'SuSiE_adj_finemapping'
RES=W/'results'; FIG=W/'figures'
METHODS=['BOLT-LMM','REGENIE-LOCO']
MDIR={'BOLT-LMM':'BOLT_LMM','REGENIE-LOCO':'REGENIE'}
man=pd.read_csv(OLD/'loci/manifest.tsv',sep='\t')
truth=set(pd.read_csv(B/'snp.info.chr1_18.10k.csv').SNP)
assert len(truth)==10000

def read_fit(row):
    root=W/'susie_results'/MDIR[row.method]/row.phenotype/row.locus_id
    assert (root/'DONE').exists(), str(root)
    f=pd.read_csv(root/'fit_summary.tsv',sep='\t')
    assert len(f)==1 and f.locus_id.iloc[0]==row.job_id
    s=c=None
    if f.fit_status.iloc[0]=='ok':
        s=pd.read_csv(root/'snp_pip.tsv',sep='\t')
        if (root/'credible_sets.tsv').stat().st_size>1:
            c=pd.read_csv(root/'credible_sets.tsv',sep='\t')
        assert len(s)==f.matched.iloc[0]
        assert f.n_cs.iloc[0]==(0 if c is None else len(c))
    return f,s,c

with ThreadPoolExecutor(max_workers=12) as pool:
    inputs=list(pool.map(read_fit,man.itertuples()))
F=pd.concat([x[0] for x in inputs],ignore_index=True)
S=pd.concat([x[1] for x in inputs if x[1] is not None],ignore_index=True)
C=pd.concat([x[2] for x in inputs if x[2] is not None],ignore_index=True)
del inputs
assert len(F)==len(man)==13301 and F.locus_id.is_unique
assert F.n_eff.eq(20000).all() and F.fit_mode.eq('ERV_FALSE').all()
assert F.dropped.eq(0).all(), 'Unexpected SNP filtering; investigate before comparing'
assert S.PIP.between(0,1).all() and not S.duplicated(['locus_id','SNP']).any()
assert not C.duplicated(['locus_id','cs']).any()
D=man.merge(F.drop(columns=['method','phenotype','chr']),left_on='job_id',right_on='locus_id',validate='one_to_one')
D.to_csv(RES/'per_locus_diagnostics.tsv',sep='\t',index=False)
S.to_csv(RES/'snp_level_all.tsv.gz',sep='\t',index=False)
C.to_csv(RES/'credible_sets_all.tsv',sep='\t',index=False)
F.groupby(['method','fit_status','converged'],dropna=False).size().rename('n_loci').reset_index().to_csv(RES/'fit_status_summary.tsv',sep='\t',index=False)

oldD=pd.read_csv(OLD/'loco/results/per_locus_diagnostics.tsv',sep='\t').rename(columns={'locus':'job_id'})
oldS=pd.read_csv(OLD/'loco/results/snp_level_all.tsv.gz',sep='\t')
oldC=pd.read_csv(OLD/'loco/results/credible_sets_all.tsv',sep='\t')
# Ensure the requested comparison really keeps the same harmonized input statistics.
check=S[['locus_id','SNP','effect_allele','z']].merge(oldS[['locus_id','SNP','effect_allele','z']],on=['locus_id','SNP'],suffixes=('_new','_old'),how='outer',indicator=True,validate='one_to_one')
assert check._merge.eq('both').all() and check.effect_allele_new.eq(check.effect_allele_old).all()
max_z_diff=float(np.abs(check.z_new-check.z_old).max())
assert max_z_diff < 1e-10
(W/'qc/comparison_input_check.json').write_text(json.dumps({'n_SNP_locus_rows':len(check),'max_abs_z_difference':max_z_diff,'alleles_identical':True,'all_manifest_loci_attempted':True},indent=2)+'\n')
del check

def metrics(d,s,c,analysis,restrict=None):
    d=d[(d.fit_status=='ok') & (d.converged==True)]
    if restrict is not None: d=d[d.job_id.isin(restrict)]
    s=s[s.locus_id.isin(d.job_id)]
    c=c[c.locus_id.isin(d.job_id)].copy()
    c['members']=c.SNPs.str.split(',').map(set)
    assert (c.members.map(len)==c.cs_size).all()
    c['has_causal']=c.members.map(lambda x:bool(x & truth))
    if analysis=='Standard SuSiE-RSS' and restrict is None:
        c.drop(columns='members').to_csv(RES/'credible_sets_truth.tsv.gz',sep='\t',index=False)
    rows=[]
    for (method,phenotype),q in d.groupby(['method','phenotype']):
        ss=s[s.locus_id.isin(q.job_id)]
        cc=c[c.locus_id.isin(q.job_id)]
        causal_analyzed=set(ss.SNP)&truth
        recovered=set().union(*cc.members.tolist())&truth
        assert recovered<=causal_analyzed
        high=ss[ss.PIP>=.9].drop_duplicates('SNP')
        rows.append(dict(analysis=analysis,method=method,phenotype=phenotype,n_loci=len(q),n_CS=len(cc),n_CS_with_causal=int(cc.has_causal.sum()),
            n_causal_analyzed=len(causal_analyzed),n_causal_recovered=len(recovered),coverage=100*cc.has_causal.mean(),
            power=100*len(recovered)/len(causal_analyzed) if causal_analyzed else np.nan,median_CS_size=cc.cs_size.median(),
            genome_wide_recovery=100*len(recovered)/len(truth),n_PIP_ge_0_9=len(high),n_PIP_ge_0_9_causal=int(high.SNP.isin(truth).sum())))
    return pd.DataFrame(rows)

M=metrics(D,S,C,'Standard SuSiE-RSS')
O=metrics(oldD,oldS,oldC,'LOCO-adjusted SuSiE')
assert len(M)==len(O)==20
M.to_csv(RES/'metrics_by_phenotype.csv',index=False)
A=pd.concat([M,O],ignore_index=True)
A.to_csv(RES/'comparison_by_phenotype.csv',index=False)
keys=['coverage','power','median_CS_size']
A.groupby(['analysis','method'])[keys].agg(['mean','min','max']).to_csv(RES/'comparison_summary.csv')
common=set(D.loc[(D.fit_status=='ok') & (D.converged==True),'job_id'])&set(oldD.loc[(oldD.fit_status=='ok') & (oldD.converged==True),'job_id'])
paired=pd.concat([metrics(D,S,C,'Standard SuSiE-RSS',common),metrics(oldD,oldS,oldC,'LOCO-adjusted SuSiE',common)])
paired.to_csv(RES/'comparison_common_converged_loci.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':13,'axes.labelsize':11,
    'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
colors=['#2878A6','#D47936']

def plot(data,comparison=False):
    fig,axes=plt.subplots(1,3,figsize=(11.5,4.5))
    for ax,key,title,ylabel,letter in zip(axes,keys,['Credible-set coverage','Causal-variant recovery','Credible-set size'],
        ['Sets containing a causal SNP (%)','Causal SNPs recovered (%)','Median SNPs per set'],'ABC'):
        for i,method in enumerate(METHODS):
            analyses=['Standard SuSiE-RSS','LOCO-adjusted SuSiE'] if comparison else ['Standard SuSiE-RSS']
            for j,analysis in enumerate(analyses):
                vals=data[(data.method==method)&(data.analysis==analysis)].sort_values('phenotype')[key].to_numpy()
                x=i+(-.15 if j==0 else .15) if comparison else i
                width=.08 if comparison else .10
                ax.scatter(x+np.linspace(-width,width,len(vals)),vals,s=30,
                    facecolors=colors[i] if j==0 else 'white',edgecolors=colors[i],linewidths=.9,alpha=.8,zorder=3)
                ax.plot([x-width-.03,x+width+.03],[vals.mean()]*2,color=colors[i],lw=2.6,ls='-' if j==0 else '--',zorder=4)
                if not comparison: ax.annotate(f'{vals.mean():.1f}',(x,vals.max()),xytext=(0,9),textcoords='offset points',ha='center',fontweight='bold',color=colors[i])
        ax.set_xticks([0,1]);ax.set_xticklabels(METHODS,fontsize=10);ax.set_xlim(-.55,1.55)
        ax.set_title(title,pad=18);ax.set_ylabel(ylabel)
        ax.text(-.18,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=15)
        ax.yaxis.grid(True,color='#E9E9E9',lw=.7);ax.set_axisbelow(True)
        ax.set_ylim(0,105 if key=='coverage' else data[key].max()*1.35)
        if key=='coverage':
            ax.axhline(95,color='#777777',ls='--',lw=1);ax.text(1.45,96.8,'95% target',ha='right',fontsize=9,color='#666666')
    title='Standard versus LOCO-adjusted SuSiE' if comparison else 'Standard SuSiE-RSS | Unadjusted LD'
    fig.suptitle(title,fontsize=16,y=.98)
    if comparison:
        handles=[Line2D([0],[0],marker='o',color='#555555',label='Standard SuSiE-RSS',ls='-'),
                 Line2D([0],[0],marker='o',markerfacecolor='white',color='#555555',label='LOCO-adjusted SuSiE',ls='--')]
        fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.055),ncol=2,frameon=False,fontsize=10)
        footer='Same converged loci; residual variance fixed (standard) versus estimated (LOCO-adjusted)'
    else: footer='Each dot: one phenotype (10 total)   |   Horizontal lines: means   |   95% credible sets'
    fig.text(.5,.025,footer,ha='center',fontsize=9,color='#555555')
    fig.subplots_adjust(left=.08,right=.98,bottom=.25 if comparison else .21,top=.78,wspace=.50)
    name='susie_vs_loco_comparison' if comparison else 'susie_rss_summary'
    for ext in ['png','pdf','svg']:fig.savefig(FIG/f'{name}.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
plot(M)
plot(paired,True)
(FIG/'figure_caption.txt').write_text('Standard SuSiE-RSS with ordinary Pearson LD from the 20,000 GWAS animals, n=20,000, L=10, coverage=0.95, min_abs_corr=0.5, estimate_residual_variance=FALSE, and max_iter=1000. Same BOLT-LMM/REGENIE-LOCO summary statistics and GWAS-defined loci as the previous analysis. All successful converged loci are included, including null-region loci. Each dot is one phenotype; lines are means across 10 phenotypes. A: empirical coverage, the percentage of credible sets containing at least one causal SNP. B: conditional power, unique causal SNPs in at least one credible set divided by unique causal SNPs analyzed in converged loci. This denominator excludes untested loci; genome-wide recovery using all 10,000 causal SNPs is provided in the CSV. C: median set size per phenotype. Truth is SNP membership in snp.info.chr1_18.10k.csv and is read only after fitting.\nComparison figure uses loci converged in both pipelines. The previous LOCO analysis used adjusted LD, locus-specific effective sample size, and estimated residual variance; the new standard analysis uses ordinary LD, N=20,000 and fixed residual variance. Thus this comparison is not an isolated test of LD adjustment alone. The retained mixed-model GWAS statistics may be less compatible with ordinary LD; the s_rss diagnostic is saved for each locus.\nMetric references: Wang et al. (2020), https://doi.org/10.1111/rssb.12388; Zou et al. (2022), https://doi.org/10.1371/journal.pgen.1010299.\n')
print(F.groupby(['method','fit_status','converged']).size().to_string())
print(A.groupby(['analysis','method'])[keys].mean().to_string())
(W/'status/EVALUATION.COMPLETE').write_text('Verified inputs, evaluated fits, and generated figures.\n')
