#!/usr/bin/env python3
"""Redraw manuscript figures from verified aggregate tables, not restricted SNP data."""
import argparse,sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from benchmark_core import METHODS
LABELS={m:m.replace('_',' ') for m in METHODS}
COLORS=dict(zip(METHODS,['#0072B2','#D55E00','#E69F00','#009E73','#CC79A7','#7B3294','#56B4E9','#999933','#222222']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})

def draw(data,out):
    data,out=Path(data),Path(out);out.mkdir(parents=True,exist_ok=True)
    def read(name):return pd.read_csv(data/name,sep='\t')
    def save(fig,name):
        fig.savefig(out/(name+'.pdf'),bbox_inches='tight');fig.savefig(out/(name+'.png'),dpi=180,bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,2.8))
    for x in [0,6,12]:ax.add_patch(Rectangle((x,0),2,1,facecolor='#0072B2',alpha=.65))
    for x in [2,8]:
        ax.add_patch(Rectangle((x,0),4,1,facecolor='#eeeeee',edgecolor='white'))
        ax.add_patch(Rectangle((x+1,.1),2,.8,facecolor='#E69F00',alpha=.5))
        for k in range(1,4):ax.plot([x+1+k*.5]*2,[.1,.9],color='white')
    ax.set(xlim=(0,14),ylim=(-.1,1.3),xlabel='Chromosome position (Mb)',yticks=[],title='Figure 2. Repeating 2 Mb effect / 4 Mb null design')
    ax.text(1,.5,'Effect',ha='center',va='center',color='white');ax.text(4,.5,'Central null\n4 × 500 kb',ha='center',va='center')
    ax.text(4,1.12,'Full null region',ha='center');save(fig,'figure2_windows')
    c=read('method_summary/significant_snps_by_pheno.tsv')
    a=c.pivot(index='method',columns='phenotype',values='n_significant_loci').reindex(METHODS)
    snp=c.pivot(index='method',columns='phenotype',values='n_significant_snps').reindex(METHODS)
    fig,ax=plt.subplots(figsize=(13,7));im=ax.imshow(a,cmap='Blues',aspect='auto')
    for i in range(9):
        for j in range(10):ax.text(j,i,f'{int(a.iloc[i,j]):,}\n({int(snp.iloc[i,j]):,})',ha='center',va='center',fontsize=8,color='white' if a.iloc[i,j]>1000 else '#111111')
    ax.set(xticks=range(10),xticklabels=range(1,11),yticks=range(9),yticklabels=[LABELS[m] for m in METHODS],xlabel='Simulated phenotype',title='Figure 3. Significant 500 kb loci (significant SNPs in parentheses)')
    fig.colorbar(im,ax=ax,label='Loci');save(fig,'figure3_counts')
    ov=read('overlap_summary/locus_overlap_by_pheno.tsv')
    fig,axs=plt.subplots(5,2,figsize=(13,23),layout='constrained')
    short=['BOLT','REG-L','REG-W','fast','Farm-N','Farm-K','BL-N','BL-K','SLEMM']
    for ph,ax in zip(range(1,11),axs.flat):
        a=ov[ov.phenotype==ph].pivot(index='method_a',columns='method_b',values='jaccard').reindex(index=METHODS,columns=METHODS)
        im=ax.imshow(a,vmin=0,vmax=1,cmap='YlGnBu')
        for i in range(9):
            for j in range(9):
                v=a.iloc[i,j];ax.text(j,i,'NA' if np.isnan(v) else f'{v:.2f}',ha='center',va='center',fontsize=7,color='white' if v>.6 else 'black')
        ax.set(xticks=range(9),xticklabels=short,yticks=range(9),yticklabels=short,title=f'Phenotype {ph}')
        ax.tick_params(axis='x',rotation=55,labelsize=8);ax.tick_params(axis='y',labelsize=8)
    fig.suptitle('Figure 5. Jaccard overlap of strict-threshold 500 kb loci\nREG-L/W: REGENIE with/without LOCO; N/K: NYC/Kansas',fontsize=13)
    fig.colorbar(im,ax=list(axs.flat),shrink=.3,label='Jaccard index');save(fig,'figure5_jaccard')
    s=read('spillover_summary/spillover_summary.tsv').set_index('method').reindex(METHODS)
    lam=read('lambda_gc_summary/lambda_gc_summary.tsv').set_index('method').reindex(METHODS)
    sp=read('power_summary/spatial_enrichment_summary.tsv')
    fig,axs=plt.subplots(2,2,figsize=(14,11),layout='constrained')
    for ax,frame,col,title,factor in [(axs[0,0],s,'spillover_rate','A  Central-null spillover loci',100),(axs[0,1],lam,'lambda_gc','B  Historical central-null λGC',1)]:
        ax.bar(range(9),factor*frame['mean_'+col],yerr=factor*frame['sd_'+col],color=[COLORS[m] for m in METHODS],capsize=3)
        ax.set(xticks=range(9),xticklabels=[LABELS[m] for m in METHODS],title=title,ylabel='Percent of 1,480 windows' if factor==100 else 'λGC')
        ax.tick_params(axis='x',rotation=60,labelsize=8)
        if factor==100:
            for i,value in enumerate(frame['mean_'+col]):ax.text(i,100*value+100*frame['sd_'+col].iloc[i]+.3,f'{100*value:.1f}%',ha='center',fontsize=8)
    for ax,methods,title in [(axs[1,0],METHODS,'C  All configurations'),(axs[1,1],METHODS[2:],'D  Lower-enrichment configurations')]:
        for m in methods:
            q=sp[sp.method==m];ax.plot((q.distance_index-.5)*.05,q.mean_observed_expected/1e6,label=LABELS[m],color=COLORS[m])
        ax.set(xlabel='Distance from nearest effect boundary (Mb)',ylabel='Observed / expected significant SNPs (×10⁶)',title=title,xlim=(0,2))
        ax.legend(fontsize=7,frameon=False,ncol=2)
    save(fig,'figure6_spillover_inflation')
    e=read('equation13_summary/equation13_summary.tsv');methods=['BOLT-LMM','REGENIE_LOCO','REGENIE_woLOCO','SLEMM']
    for region,name in [('Genome-wide','figure7_equation13'),('Central-null','supplement_equation13_central')]:
        fig,axs=plt.subplots(2,2,figsize=(10,8),layout='constrained')
        for m,ax in zip(methods,axs.flat):
            q=e[(e.region==region)&(e.method==m)].sort_values('bin');ax.errorbar(q.predicted_chisq,q.mean_observed_chisq,yerr=[q.mean_observed_chisq-q.ci_lower,q.ci_upper-q.mean_observed_chisq],fmt='o',ms=4,color=COLORS[m],capsize=2)
            lim=max(q.predicted_chisq.max(),q.ci_upper.max())*1.04;ax.plot([0,lim],[0,lim],'--',color='#888888',lw=1)
            ax.set(xlabel='Predicted χ² = 1 + cached NCP',ylabel='Mean observed χ²',title=LABELS[m],xlim=(0,lim),ylim=(0,lim))
        fig.suptitle(f'Equation 13: {region}; 20 bins; 95% intervals across 10 phenotypes\nCached upstream predictions; Equation 12 generator unresolved',fontsize=11);save(fig,name)
    f=read('finemapping_summary/high_pip_summary.tsv')
    fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained');x=np.arange(4)
    labels=[r.method.replace('_LOCO','')+'\n'+r.pipeline for r in f.itertuples()]
    axs[0].bar(x-.18,f.causal_percent,width=.36,label='Observed causal precision',color='#0072B2');axs[0].bar(x+.18,f.maf_expected_percent,width=.36,label='MAF-matched expectation',color='#bbbbbb')
    for i,r in enumerate(f.itertuples()):axs[0].text(i-.18,r.causal_percent+.3,f'{r.causal_percent:.1f}%',ha='center',fontsize=8)
    axs[0].set(ylabel='Percent of high-PIP calls',ylim=(0,17),title='A  Causal precision at PIP ≥ 0.9');axs[0].legend(fontsize=8,frameon=False)
    axs[1].bar(x,f.null_percent,color='#E69F00')
    for i,r in enumerate(f.itertuples()):axs[1].text(i,r.null_percent+1,f'{r.n_null:,}/{r.n_high_pip:,}\n{r.null_percent:.2f}%',ha='center',fontsize=9)
    axs[1].set(ylabel='Percent of high-PIP calls',ylim=(0,78),title='B  Calls in entire null windows')
    for ax in axs:ax.set_xticks(x,labels)
    fig.suptitle('Figure 8. Recount of saved converged fine-mapping fits',fontsize=12);save(fig,'figure8_finemapping')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--out',required=True);a=p.parse_args();draw(a.data,a.out)
