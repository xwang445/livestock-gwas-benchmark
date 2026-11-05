#!/usr/bin/env python3
"""Draw revised Manhattan panels from restricted files; publish plots/counts only.

Missing configurations are labeled pending. A pooled comparison includes only
configurations with all ten phenotypes; --require-complete enforces a final run.
"""
import argparse
import os
from pathlib import Path
import numpy as np
import pandas as pd
os.environ.setdefault('MPLCONFIGDIR','/tmp/gwas_revision_mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from association.significance import ALPHA, significant
from association.derive_association_tables import config as association_source

METHODS = ('BOLT-LMM','REGENIE_LOCO','REGENIE_WO_LOCO','fastGWA',
           'FarmCPU_NYC','FarmCPU_Kansas','BLINK_NYC','BLINK_Kansas','SLEMM','GCTA_MLMA')
COLORS = dict(zip(METHODS,('#0072B2','#D55E00','#E69F00','#009E73','#CC79A7',
                         '#7B3294','#56B4E9','#999933','#222222','#A6761D')))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})

def source(root, gcta, method, ph):
    return association_source(root,method,ph,gcta)[0]

def load(path,method):
    if method.startswith(('FarmCPU','BLINK')):
        d=pd.read_csv(path);cols=['Chr','Pos','SNP','P.value']
    else:
        d=pd.read_csv(path,sep=r'\s+')
        cols=(['CHR','BP','SNP','P_BOLT_LMM_INF'] if method=='BOLT-LMM' else
              ['CHROM','GENPOS','ID','LOG10P'] if method.startswith('REGENIE') else
              ['#CHROM','POS','ID','P'] if method=='SLEMM' else
              ['Chr','bp','SNP','p'] if method=='GCTA_MLMA' else ['CHR','POS','SNP','P'])
    v=d[cols].copy();v.columns=['chr','bp','snp','p']
    native=v.p.to_numpy(copy=True)
    is_log10=method.startswith('REGENIE')
    v['significant']=significant(native,is_log10=is_log10)
    v['minus_log10p']=native if is_log10 else -np.log10(np.maximum(native,1e-300))
    if is_log10:v['p']=10.**(-native)
    v['chr']=pd.to_numeric(v.chr.astype(str).str.replace(r'^chr','',regex=True))
    assert len(v)==135358 and v.chr.between(1,18).all() and v.bp.gt(0).all()
    assert not v.snp.duplicated().any() and not v[['chr','bp']].duplicated().any()
    assert v.p.dropna().between(0,1).all()
    return v

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--gcta-dir',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--require-complete',action='store_true');a=p.parse_args()
    paths={(m,ph):source(a.root,a.gcta_dir,m,ph) for m in METHODS for ph in range(1,11)}
    missing=[(m,ph) for (m,ph),path in paths.items() if not path.is_file()]
    if missing and a.require_complete:raise ValueError(f'Incomplete analyses: {missing}')
    if a.out.exists():raise ValueError('Use a new output directory to preserve figures')
    a.out.mkdir(parents=True)
    complete=[m for m in METHODS if all(paths[m,ph].is_file() for ph in range(1,11))]
    ref=load(paths['BOLT-LMM',1],'BOLT-LMM');lengths=ref.groupby('chr').bp.max().sort_index()
    offsets=lengths.cumsum().shift(fill_value=0);centers=offsets+lengths/2
    mainfig,mainaxes=plt.subplots(2,1,figsize=(16,9),sharex=True,constrained_layout=True)
    records=[]
    for ph in range(1,11):
        fig,axes=plt.subplots(len(METHODS),1,figsize=(16,23),sharex=True,constrained_layout=True)
        ymax=0
        for idx,(m,ax) in enumerate(zip(METHODS,axes)):
            path=paths[m,ph]
            ax.set_title(f'{chr(97+idx)}  {m.replace("_"," ")}',loc='left',fontsize=10)
            if not path.is_file():
                ax.text(.5,.5,'Analysis pending: no result reported',transform=ax.transAxes,ha='center')
                records.append(dict(configuration=m,phenotype=ph,status='PENDING',n_plotted=np.nan,n_significant=np.nan,max_minus_log10p=np.nan))
                continue
            d=load(path,m);x=(d.bp+d.chr.map(offsets))/1e6
            y=d.minus_log10p;ymax=max(ymax,float(y.max()));sig=d.significant
            for parity,color in [(0,'#0072B2'),(1,'#595959')]:
                use=d.chr.mod(2).eq(parity);ax.scatter(x[use],y[use],s=.65,color=color,alpha=.65,linewidths=0,rasterized=True)
            ax.axhline(-np.log10(ALPHA),ls='--',lw=.6,color='#D55E00');ax.set_ylabel('−log₁₀(P)',fontsize=9)
            if m in complete:
                group=int(m not in ('BOLT-LMM','REGENIE_LOCO'))
                mainaxes[group].scatter(x[sig],y[sig],s=.8,color=COLORS[m],alpha=.45,linewidths=0,rasterized=True,label=m.replace('_',' ') if ph==1 else None)
            records.append(dict(configuration=m,phenotype=ph,status='COMPLETED',n_plotted=len(d),n_significant=int(sig.sum()),max_minus_log10p=float(y.max())))
        for ax in axes:ax.set_ylim(0,np.ceil(ymax*1.03));ax.set_xlim(0,(offsets.iloc[-1]+lengths.iloc[-1])/1e6)
        axes[-1].set_xticks(centers/1e6);axes[-1].set_xticklabels(centers.index);axes[-1].set_xlabel('Chromosome')
        fig.suptitle(f'Phenotype {ph}: all tested SNPs\nShared vertical scale; dashed line P = 5 × 10⁻⁸'+ ('; provisional: incomplete GCTA analysis' if missing else ''),fontsize=12)
        for ext in ('pdf','png'):fig.savefig(a.out/f'figureS{ph:02}_manhattan.{ext}',dpi=300,bbox_inches='tight')
        plt.close(fig);print(f'Phenotype {ph} rendered',flush=True)
    for idx,ax in enumerate(mainaxes):
        ax.set_ylim(-np.log10(ALPHA)-.15,31)
        ax.axhline(-np.log10(ALPHA),ls='--',lw=.6,color='#595959')
        ax.set_ylabel('−log₁₀(P)');ax.legend(loc='upper right',markerscale=5,fontsize=8,ncol=3,frameon=False)
        ax.set_title(chr(97+idx),loc='left')
    mainaxes[-1].set_xticks(centers/1e6);mainaxes[-1].set_xticklabels(centers.index);mainaxes[-1].set_xlabel('Chromosome')
    mainfig.suptitle(f'Significant SNPs pooled across ten phenotypes ({len(complete)} complete configurations)\nP ≤ 5 × 10⁻⁸; display cap 31; full ranges in S1–S10'+ ('; GCTA pending, excluded from pooled comparison' if missing else ''),fontsize=12)
    for ext in ('pdf','png'):mainfig.savefig(a.out/f'figure4_manhattan.{ext}',dpi=300,bbox_inches='tight')
    plt.close(mainfig)
    pd.DataFrame(records).to_csv(a.out/'manhattan_manifest.tsv',sep='\t',index=False,na_rep='NA')

if __name__=='__main__':main()
