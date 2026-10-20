#!/usr/bin/env python3
"""Summarize empirical precision of high-PIP calls under simulated diffuse polygenicity."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/susie_precision_mpl')
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
W=Path(__file__).resolve().parents[1]
SRC=W/'publication/maf_matched'; OUT=W/'publication/localization_precision'
OUT.mkdir(exist_ok=True)
d=pd.read_csv(SRC/'counts_by_phenotype.csv')
calls=pd.read_csv(SRC/'high_pip_calls.tsv.gz',sep='\t')
assert calls.PIP.ge(.9).all()
keys=['GWAS_method','fine_mapping']
classes=['Causal','Non-causal in effect windows','Null windows']
order=[(m,a) for m in ['BOLT-LMM','REGENIE-LOCO'] for a in ['SuSiE-adj','SuSiE-RSS']]
rows=[]; per=[]
for method,analysis in order:
    q=d[(d.GWAS_method==method)&(d.fine_mapping==analysis)]
    c=calls[(calls.GWAS_method==method)&(calls.fine_mapping==analysis)]
    assert not c.duplicated(['phenotype','SNP']).any()
    totals=q.groupby('class')['observed'].sum()
    n=int(totals.sum()); causal=int(totals['Causal'])
    assert n==len(c) and causal==c['class'].eq('Causal').sum()
    cq=q[q['class']=='Causal'].copy()
    for r in cq.itertuples():
        cc=c[c.phenotype==r.phenotype]
        assert len(cc)==r.n_calls
        per.append(dict(GWAS_method=method,fine_mapping=analysis,phenotype=r.phenotype,
                        n_calls=r.n_calls,n_causal=r.observed,precision_pct=100*r.observed/r.n_calls,
                        maf_expected_precision_pct=100*r.expected/r.n_calls,mean_PIP_pct=100*cc.PIP.mean()))
    precision=100*cq.observed/cq.n_calls
    expected=cq.expected.sum()
    rows.append(dict(GWAS_method=method,fine_mapping=analysis,n_calls=n,n_causal=causal,
        n_noncausal_effect=int(totals[classes[1]]),n_null=int(totals[classes[2]]),
        precision_pct=100*causal/n,noncausal_pct=100*(n-causal)/n,
        noncausal_effect_pct=100*totals[classes[1]]/n,null_pct=100*totals[classes[2]]/n,
        mean_assigned_PIP_pct=100*c.PIP.mean(),posterior_expected_causal=c.PIP.sum(),
        maf_expected_causal=expected,maf_expected_precision_pct=100*expected/n,
        causal_enrichment=causal/expected,phenotype_precision_mean=precision.mean(),
        phenotype_precision_sd=precision.std(),phenotype_precision_min=precision.min(),
        phenotype_precision_max=precision.max()))
s=pd.DataFrame(rows); p=pd.DataFrame(per)
s.to_csv(OUT/'Table_3_numeric.csv',index=False)
p.to_csv(OUT/'precision_by_phenotype.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
causal_color='#247D91'
noncausal_color='#B8BEC5'
pip_color='#465364'
labels=[m+'\n'+a for m,a in order]
# GSE half-page width; embedded TrueType sans-serif font.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
                     'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.6,
                     'xtick.major.width':.6,'ytick.major.width':.6})
fig,ax=plt.subplots(figsize=(85/25.4,115/25.4))
x=np.array([0.,1.,2.3,3.3])
width=.48
ax.bar(x,s.precision_pct,width=width,color=causal_color,edgecolor='white',linewidth=.4)
ax.bar(x,s.noncausal_pct,bottom=s.precision_pct,width=width,color=noncausal_color,edgecolor='white',linewidth=.4)
assert np.allclose(s.precision_pct+s.noncausal_pct,100)
for xx,r in zip(x,s.itertuples()):
    ax.text(xx,r.precision_pct+r.noncausal_pct/2,f'{r.noncausal_pct:.1f}',ha='center',va='center',color='#202A34',fontsize=7)
    ax.text(xx,r.precision_pct/2,f'{r.precision_pct:.1f}',ha='center',va='center',color='white',fontsize=7)
    expected=r.maf_expected_precision_pct
    ax.plot([xx-width/2+.025,xx+width/2-.025],[expected,expected],color='white',lw=.7,ls=(0,(2,1.4)),zorder=4)
    ax.annotate(f'Exp.\n{expected:.1f}%',xy=(xx,expected),xytext=(xx,23),ha='center',va='center',fontsize=7,color='#202A34',arrowprops=dict(arrowstyle='-',color='#202A34',lw=.6,shrinkA=3,shrinkB=1),zorder=6)
    ax.text(xx,103,f'{r.causal_enrichment:.2f}×',ha='center',fontsize=7,color='#202A34')
ax.set_xlim(-.55,3.85)
ax.set_xticks(x);ax.set_xticklabels(['SuSiE-adj','SuSiE-RSS']*2,fontsize=7,rotation=45,ha='right',rotation_mode='anchor')
for xx,label in [(.5,'BOLT-LMM'),(2.8,'REGENIE-LOCO')]:
    ax.text(xx,-.19,label,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=7.5)
ax.set_ylim(0,109);ax.set_yticks(np.arange(0,101,20))
ax.set_ylabel('High-PIP SNP calls (%)',labelpad=3)
ax.tick_params(axis='both',length=2,pad=2)
ax.yaxis.grid(True,color='#E7E7E7',lw=.4);ax.set_axisbelow(True)
fig.legend(handles=[Patch(facecolor=causal_color,label='Causal SNPs (precision)'),
                    Patch(facecolor=noncausal_color,label='Non-causal SNPs'),
                    Line2D([],[],ls=(0,(2,1.4)),lw=.7,color='#58626C',label='MAF-matched expectation')],
           loc='lower center',bbox_to_anchor=(.54,.015),ncol=1,frameon=False,fontsize=7.5,
           handlelength=1.7,labelspacing=.35)
fig.subplots_adjust(left=.17,right=.98,top=.97,bottom=.31)
for ext in ['png','pdf','svg']:fig.savefig(OUT/('Figure_3.'+ext),dpi=600,facecolor='white')
fig.savefig(OUT/'Figure_3.tiff',dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
# Four-row publication table, with counts and empirical precision shown together.
t=[]
for r in s.itertuples():
    t.append([r.GWAS_method,r.fine_mapping,f'{r.n_calls:,}',f'{r.n_causal:,} ({r.precision_pct:.1f}%)',
              f'{r.n_noncausal_effect:,} ({r.noncausal_effect_pct:.1f}%)',f'{r.n_null:,} ({r.null_pct:.1f}%)',
              f'{r.mean_assigned_PIP_pct:.1f}%',f'{r.maf_expected_precision_pct:.1f}%',f'{r.causal_enrichment:.2f}'])
cols=['GWAS method','Fine-mapping','Total calls','Causal, n (%)','Non-causal effect, n (%)','Null windows, n (%)','Mean PIP','MAF-expected causal (%)','Causal O/E']
pd.DataFrame(t,columns=cols).to_csv(OUT/'Table_3.csv',index=False)
fig,ax=plt.subplots(figsize=(15.5,3.9));ax.axis('off')
head=['GWAS\nmethod','Fine-mapping','Total\ncalls','Causal\nn (%)','Non-causal in effect\nwindows, n (%)','Null windows\nn (%)','Mean\nassigned PIP','MAF-expected\ncausal fraction','Causal\nO/E']
tab=ax.table(cellText=t,colLabels=head,loc='center',cellLoc='center',colWidths=[.105,.095,.07,.12,.15,.13,.11,.13,.09])
tab.auto_set_font_size(False);tab.set_fontsize(10)
for (r,c),cell in tab.get_celld().items():
    cell.set_height(.19 if r==0 else .16);cell.set_edgecolor('#D6DCE0');cell.set_linewidth(.6)
    if r==0:cell.set_facecolor('#233E4C');cell.set_text_props(color='white',weight='bold')
    else:cell.set_facecolor('#EFF4F6' if r<=2 else 'white')
ax.set_title('High-PIP calls from LOCO inputs under diffuse polygenicity',fontsize=15,pad=14)
fig.text(.5,.08,'Counts pooled across 10 phenotypes; SNPs counted once per phenotype. O/E = observed / MAF-matched expected causal count.',ha='center',fontsize=10)
fig.subplots_adjust(left=.025,right=.975,top=.8,bottom=.19)
for ext in ['png','pdf']:fig.savefig(OUT/('Table_3.'+ext),dpi=300,facecolor='white')
plt.close(fig)
base=('Analysis is restricted to BOLT-LMM and REGENIE-LOCO inputs under the diffuse polygenic architecture simulated here. '
      'All PIP ≥ 0.9 SNPs from the 13,301 converged locus fits in each pipeline were included; no lead-window stratification '
      'or top-SNP restriction was applied. Loci were selected using genome-wide significant lead SNPs (P ≤ 5 × 10⁻⁸). '
      'SNPs were counted once within each GWAS method, fine-mapping approach, and phenotype; counts across the 10 phenotypes '
      'are SNP–phenotype calls, not globally unique variants. Empirical precision is the number of simulated causal calls '
      'divided by all high-PIP calls. Classes use each SNP’s own causal status and position: causal; non-causal in effect '
      'windows; or null windows. The latter two classes together are non-causal. The mean assigned PIP is the average '
      'saved posterior probability among selected calls, shown as a model-implied causal fraction for comparison with '
      'the observed causal fraction; it is not an independent empirical estimate. MAF-matched expectations use 10 pool-specific '
      'MAF quantile bins and E_c = sum_k n_k p_k(c), with MAF from the 20,000 GWAS animals and pools restricted to SNPs '
      'analyzed for the same fine-mapping approach, GWAS method, and phenotype. Expected counts are summed across '
      'phenotypes before dividing by total calls. ')
(OUT/'Figure_3_caption.txt').write_text('Figure 3. Observed precision and MAF-matched expectation for high-PIP SNP calls from LOCO inputs under the simulated diffuse polygenic architecture. Each GWAS–fine-mapping combination has one vertical stacked bar showing observed causal SNPs (teal) and non-causal SNPs (gray). Observed precision is the simulated causal count divided by all PIP ≥ 0.9 calls. White dashed lines within the bars mark the MAF-matched expected causal percentage at its numerical y-axis position; callouts above the causal segments prefixed Exp. give the expected percentages and point to the dashed lines. Numbers in the teal segments give observed causal percentages, and numbers in gray segments give observed non-causal percentages. All in-bar numbers are percentages. Dashed lines are reference markers, not additional bar segments. Non-causal SNPs include all selected SNPs outside the simulated causal set, irrespective of genomic window. The expected causal count is E_c = sum_k n_k p_k(c), computed separately for each fine-mapping approach, GWAS method, and phenotype. Pools contain unique SNPs analyzed in the corresponding fine-mapped loci, partitioned into 10 MAF quantile bins using MAF from the 20,000 GWAS animals. In bin k, n_k is the number of high-PIP calls and p_k(c) is the causal fraction of pool SNPs. This is the exact analytical expectation of MAF-matched draws, without Monte Carlo sampling. Observed and expected counts are summed across the 10 phenotypes before conversion to percentages; a SNP can contribute once in each phenotype. Enrichment above each bar is the ratio of pooled observed to expected causal counts. All calls in converged loci are included, without lead-window stratification or top-SNP restriction. Percentages are descriptive pooled estimates; no uncertainty intervals or significance tests are shown.\n')
(OUT/'Table_3_caption.txt').write_text('Table 3. '+base+'Counts and percentages are pooled across phenotypes. Causal O/E is the ratio of summed observed causal counts to summed MAF-matched expected causal counts. Rounding uses unrounded underlying values. No hypothesis tests or independence-based confidence intervals are reported.\n')
(OUT/'results_text.txt').write_text('Under the diffuse polygenic architecture examined here, fine-mapping of genome-wide significant loci using BOLT-LMM and REGENIE-LOCO inputs yielded limited SNP-level precision. Across SuSiE-adj and SuSiE-RSS, only 12.2–13.3% of calls with PIP ≥ 0.9 were simulated causal variants, whereas 86.7–87.8% were non-causal (Figure 3; Table 3). The mean assigned PIP among these calls was 98.2–98.3%, substantially exceeding the observed causal fraction. Causal calls were enriched 1.40–1.51-fold above the MAF-matched expectation, indicating some causal discrimination but low absolute precision. Thus, in this simulation setting, strong LOCO association signals did not translate into reliable SNP-level causal identification, and high fine-mapping posterior probabilities frequently selected non-causal variants.\n')
(OUT/'provenance.json').write_text(json.dumps(dict(script=str(Path(__file__).resolve()),
    inputs={str(SRC/f):hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['counts_by_phenotype.csv','high_pip_calls.tsv.gz']},
    aggregation='sum counts across phenotypes; ratios of sums',n_phenotypes=10,pip_threshold=.9,
    scope='LOCO GWAS inputs under the simulated diffuse polygenic architecture'),indent=2)+'\n')
print(s.to_string(index=False))
