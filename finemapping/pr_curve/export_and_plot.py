#!/usr/bin/env python3
"""Export display-only PR vertices and plot paired RSS fits without variant IDs.

Exact AUPRC/average precision always comes from the full comparison table.
Sampling vertices for display never supplies the reported numerical areas.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

KEYS=['fit_mode','method','phenotype','universe']
COLORS={'BOLT-LMM':'#0072B2','REGENIE-LOCO':'#D55E00'}

def export(source,target):
    rows=[];seen={};last={}
    for chunk in pd.read_csv(source,sep='\t',chunksize=200000):
        for key,g in chunk.groupby(KEYS,sort=False):
            start=seen.get(key,0);indices=np.arange(start,start+len(g))
            rows.append(g.iloc[np.flatnonzero(indices%200==0)][KEYS+['precision','recall']])
            seen[key]=start+len(g);last[key]=g.iloc[[-1]][KEYS+['precision','recall']]
    sampled=pd.concat(rows+list(last.values()),ignore_index=True).drop_duplicates()
    target.parent.mkdir(parents=True,exist_ok=True)
    sampled.to_csv(target,sep='\t',index=False)

def main():
    p=argparse.ArgumentParser();p.add_argument('--full-curves',type=Path)
    p.add_argument('--display-curves',type=Path,required=True)
    p.add_argument('--summary',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.full_curves:
        if a.display_curves.exists():raise ValueError('Display curve export exists; preserve it and choose a new path')
        export(a.full_curves,a.display_curves)
    d=pd.read_csv(a.display_curves,sep='\t');s=pd.read_csv(a.summary,sep='\t')
    assert len(d.groupby(KEYS))==80 and d.precision.between(0,1).all() and d.recall.between(0,1).all()
    a.out.mkdir(parents=True,exist_ok=True)
    grid=np.linspace(0,1,201)
    for universe in d.universe.unique():
        fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
        for idx,(method,ax) in enumerate(zip(COLORS,axes)):
            q=d[(d.method==method)&(d.universe==universe)]
            for mode,style in [('ERV_FALSE','--'),('ERV_TRUE','-')]:
                curves=[]
                for _,g in q[q.fit_mode==mode].groupby('phenotype'):
                    g=g.sort_values('recall');curves.append(np.interp(grid,g.recall,g.precision))
                mean=np.mean(curves,axis=0)
                ax.plot(grid,mean,style,color=COLORS[method],label=f'Residual variance estimated: {mode=="ERV_TRUE"}')
            base=float(s[(s.method==method)&(s.universe==universe)].random_precision_baseline_mean.mean())
            ax.axhline(base,color='#777777',linestyle=':',label=f'Random ranking ({base:.4f})')
            ax.set(title=f'{chr(97+idx)}  {method}',xlabel='Recall',ylabel='Causal precision',xlim=(0,1),ylim=(0,1))
            ax.legend(fontsize=8,frameon=False)
        fig.suptitle('SuSiE-RSS: mean of ten phenotype-specific display curves\n'+universe.replace('_',' '),fontsize=11)
        for ext in ('pdf','png'):fig.savefig(a.out/f'pr_sensitivity_{universe}.{ext}',dpi=300,bbox_inches='tight')
        plt.close(fig)

if __name__=='__main__':main()
