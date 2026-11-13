#!/usr/bin/env python3
"""Current inclusive Manhattan entry point; --historical selects the old strict audit.

The active implementation uses native P<=5e-8 (or native LOG10P >= -log10(alpha))
and the updated SLEMM fit. The function below is an explicitly historical strict
redraw. No SNP or animal IDs appear in generated plots.
"""
import argparse,sys,subprocess
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from benchmark_core import METHODS,ALPHA,harmonize,input_path,output_table
from plot_aggregates import COLORS,LABELS

def draw_historical(root,out):
    root,out=Path(root),Path(out);out.mkdir(parents=True,exist_ok=True)
    if (out/'manhattan_manifest.tsv').exists():raise ValueError('Completed output exists; use a new output directory')
    ref=harmonize(input_path(root,'BOLT-LMM',1,historical=True),'BOLT-LMM',1,135358)
    lengths=ref.groupby('chr').bp.max().sort_index();offset=lengths.cumsum().shift(fill_value=0)
    centers=offset+lengths/2
    main,mainaxes=plt.subplots(2,1,figsize=(16,8),sharex=True,layout='constrained');records=[]
    for ph in range(1,11):
        fig,axes=plt.subplots(9,1,figsize=(16,20),sharex=True,layout='constrained');maximum=0.
        for method,ax in zip(METHODS,axes):
            d=harmonize(input_path(root,method,ph,historical=True),method,ph,135358)
            x=d.bp+d.chr.map(offset);y=-np.log10(np.maximum(d.p,1e-300));maximum=max(maximum,float(y.max()))
            for parity,color in [(0,'#0072B2'),(1,'#595959')]:
                use=d.chr.mod(2).eq(parity);ax.scatter(x[use]/1e6,y[use],s=.65,color=color,alpha=.65,linewidths=0,rasterized=True)
            ax.axhline(-np.log10(ALPHA),ls='--',lw=.6,color='#D55E00');ax.set_title(LABELS[method],loc='left',fontsize=10);ax.set_ylabel('−log₁₀(P)',fontsize=9)
            sig=d.p.lt(ALPHA);group=0 if method in ['BOLT-LMM','REGENIE_LOCO'] else 1
            mainaxes[group].scatter(x[sig]/1e6,y[sig],s=.8,color=COLORS[method],alpha=.45,linewidths=0,rasterized=True,label=LABELS[method] if ph==1 else None)
            records.append(dict(phenotype=ph,method=method,n_plotted=len(d),n_significant=int(sig.sum()),max_minus_log10p=float(y.max()),source='restricted historical GWAS panel; see derived_data/historical_2026-11-04/input_manifest.tsv'))
        for ax in axes:ax.set_ylim(0,np.ceil(maximum*1.03));ax.set_xlim(0,(offset.iloc[-1]+lengths.iloc[-1])/1e6)
        axes[-1].set_xticks(centers/1e6,centers.index);axes[-1].set_xlabel('Chromosome')
        fig.suptitle(f'Supplementary Figure S{ph}. All tested SNPs, phenotype {ph}\nShared vertical scale; dashed line P = 5 × 10⁻⁸',fontsize=13)
        for ext in ['png','pdf']:fig.savefig(out/f'figureS{ph:02}_manhattan.{ext}',dpi=150,bbox_inches='tight')
        plt.close(fig);print(f'Manhattan phenotype {ph} complete',flush=True)
    for ax in mainaxes:
        ax.set_ylim(-np.log10(ALPHA),31);ax.set_ylabel('−log₁₀(P)');ax.legend(loc='upper right',markerscale=5,fontsize=8,ncol=4,frameon=False)
    mainaxes[-1].set_xticks(centers/1e6,centers.index);mainaxes[-1].set_xlabel('Chromosome')
    main.suptitle('Figure 4. Significant SNPs pooled across 10 phenotypes\nDisplay truncated at 31; full ranges are shown in S1–S10',fontsize=13)
    for ext in ['png','pdf']:main.savefig(out/f'figure4_manhattan.{ext}',dpi=180,bbox_inches='tight')
    plt.close(main);output_table(pd.DataFrame(records),out/'manhattan_manifest.tsv')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',required=True)
    p.add_argument('--historical',action='store_true');p.add_argument('--gcta-dir',type=Path);p.add_argument('--require-complete',action='store_true');p.add_argument('--individual',action='store_true');a=p.parse_args()
    if a.historical and a.individual:p.error('--individual applies to current plots only')
    if a.historical:draw_historical(a.root,a.out)
    else:
        source=Path(__file__).resolve().parents[1]/'revision_2026-11-05/plot_manhattan.py'
        args=[sys.executable,str(source),'--root',str(a.root),'--gcta-dir',str(a.gcta_dir or a.root/'GCTA_mlma'),'--out',a.out]
        if a.require_complete:args.append('--require-complete')
        if a.individual:args.append('--individual')
        subprocess.run(args,check=True)
