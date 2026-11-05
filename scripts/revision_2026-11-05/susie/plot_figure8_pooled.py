#!/usr/bin/env python3
"""Final Figure 8: pooled PR curves, with data-derived labels and safe vertices."""
import argparse,hashlib,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/gwas_figure8_pooled_mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser()
p.add_argument('--results',type=Path,help='Private aggregate results directory with exact pooled curves')
p.add_argument('--display-vertices',type=Path,help='Public display vertices for an aggregate-only redraw')
p.add_argument('--metrics',type=Path,help='Exact pooled metric table; defaults to --results/pooled_finemapping_metrics.tsv')
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
if a.results is None and (a.display_vertices is None or a.metrics is None):
    p.error('Provide --results, or both --display-vertices and --metrics')
OUT=a.output.resolve();OUT.mkdir(parents=True,exist_ok=True)
msrc=a.metrics if a.metrics is not None else a.results/'pooled_finemapping_metrics.tsv'
M=pd.read_csv(msrc,sep='\t');M=M[M.fit_mode.eq('ERV_TRUE')]
if a.display_vertices is None:
    src=a.results/'pooled_precision_recall_curves.tsv.gz'
    parts=[]
    for d in pd.read_csv(src,sep='\t',chunksize=250000):
        # A zero-call endpoint has undefined empirical precision; omit its conventional value.
        parts.append(d[d.fit_mode.eq('ERV_TRUE')&d.threshold.notna()])
    D=pd.concat(parts,ignore_index=True)
else:
    src=a.display_vertices
    D=pd.read_csv(src,sep='\t')
    assert D.fit_mode.eq('ERV_TRUE').all() and D.threshold.notna().all()
display=[]
for (pipeline,method),q in D.groupby(['pipeline','method']):
    q=q.sort_values('threshold').reset_index(drop=True)
    assert np.isfinite(q[['precision','recall','threshold']]).all().all()
    assert q.precision.between(0,1).all() and q.recall.between(0,1).all()
    assert np.all(np.diff(q.recall)<=0)
    if a.display_vertices is None and len(q)>3000:
        target=np.linspace(q.recall.min(),q.recall.max(),1500)
        idx=np.minimum(np.searchsorted(-q.recall.to_numpy(),-target,side='left'),len(q)-1)
        ix=np.unique(np.r_[np.linspace(0,len(q)-1,1500).astype(int),idx,
            0,len(q)-1,q.precision.idxmin(),q.precision.idxmax()])
        q=q.iloc[ix]
    display.append(q)
V=pd.concat(display,ignore_index=True)
assert len(V)<15000
vertices=OUT/'Figure_8_display_vertices.tsv';V.to_csv(vertices,sep='\t',index=False,float_format='%.15g')
M.to_csv(OUT/'Figure_8_pooled_metrics.tsv',sep='\t',index=False,float_format='%.15g')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,
    'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,2,figsize=(7.2,3.35),sharey=True)
color={'SuSiE-adj':'#267B9B','SuSiE-RSS':'#C66B35'}
top=max(.16,float(V.precision.max())*1.10)
for ax,method,label,title in zip(axes,['BOLT-LMM','REGENIE-LOCO'],'ab',['BOLT-LMM','REGENIE (LOCO)']):
    for pipeline in ['SuSiE-adj','SuSiE-RSS']:
        q=V[V.method.eq(method)&V.pipeline.eq(pipeline)].sort_values('threshold')
        row=M[M.method.eq(method)&M.pipeline.eq(pipeline)].iloc[0]
        ax.plot(q.recall,100*q.precision,lw=1.5,color=color[pipeline],
            label='%s (AP %.2f%%)'%(pipeline,100*row.average_precision))
    baseline=float(M[M.method.eq(method)].random_precision_baseline.iloc[0])
    ax.axhline(100*baseline,color='#666666',ls=':',lw=1,
        label='Random baseline (%.2f%%)'%(100*baseline))
    ax.set_title(title,fontsize=10,pad=10)
    ax.text(-.15,1.07,label,transform=ax.transAxes,fontweight='bold',fontsize=12)
    ax.set(xlim=(0,1),ylim=(0,100*top),xlabel='Recall')
    ax.set_xticks(np.linspace(0,1,6));ax.grid(axis='y',color='#E7E7E7',lw=.5);ax.set_axisbelow(True)
    ax.legend(loc='upper right',fontsize=7,frameon=False,handlelength=2)
axes[0].set_ylabel('Precision (%)')
fig.text(.5,.02,'SNP–phenotype pairs pooled across ten phenotypes; both pipelines estimate residual variance',
    ha='center',fontsize=7,color='#4D5963')
fig.subplots_adjust(left=.09,right=.985,bottom=.20,top=.86,wspace=.25)
for ext in ['png','pdf','svg','tiff']:
    kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
    fig.savefig(OUT/('Figure_8.'+ext),dpi=600,**kw)
plt.close(fig)
caption=('Figure 8. Pooled precision–recall performance of fine-mapping LOCO associations. '
    '(a) BOLT-LMM and (b) REGENIE-LOCO inputs analyzed using SuSiE-adj and standard SuSiE-RSS. '
    'Each curve ranks SNP–phenotype pairs pooled across ten simulated phenotypes by posterior inclusion '
    'probability (PIP). Each SNP contributes once per phenotype, GWAS method and pipeline, retaining its '
    'maximum PIP across overlapping loci. The dotted line is the causal fraction of that analyzed pool '
    'and represents the precision baseline for random ranking. Average precision (AP) is the stepwise '
    'area under the full pooled curve; it is computed before selecting display vertices and differs '
    'from the mean of ten phenotype-specific AP values. Both pipelines estimate residual variance. '
    'Standard RSS uses ordinary in-sample LD and N=20,000; SuSiE-adj retains LOCO-adjusted LD and '
    'locus-specific effective sample sizes. Both use L=10, 95% credible sets, minimum absolute '
    'correlation 0.5 and maximum 1,000 iterations. All 6,732 BOLT-LMM and 6,569 REGENIE-LOCO locus '
    'fits converged in both pipelines. Loci were selected at P ≤5×10⁻⁸. Display vertices are sampled '
    'from the exact pooled curves, retaining endpoints and precision extrema; the zero-call endpoint '
    'is omitted because empirical precision is undefined. The comparison changes LD and sample-size '
    'inputs together. These descriptive curves do not treat correlated SNP calls as independent.\n')
(OUT/'Figure_8_caption.txt').write_text(caption)
assert len(caption.split())<=300
manifest=dict(status='complete',primary_RSS='estimate_residual_variance=TRUE',
    n_display_rows=len(V),display_data_contains_no_individual_identifiers=True,
    figure_caption_words=len(caption.split()),y_axis_max_percent=100*top,
    source_sha256={str(f.name):hashlib.sha256(f.read_bytes()).hexdigest() for f in [src,msrc,vertices]},
    note='Original embedded Figure8 panels preserved separately; unverified realized-PVE panel is not presented as an updated result.')
(OUT/'Figure_8_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
