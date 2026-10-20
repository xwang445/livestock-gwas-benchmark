#!/usr/bin/env python3
"""Three-class PIP >= 0.9 counts and exact MAF-matched expectations."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/susie_maf_mpl')
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

W = Path(__file__).resolve().parents[1]
B = W.parent
OUT = W/'publication/maf_matched'
OUT.mkdir(parents=True, exist_ok=True)
GENO = B/'SuSiE_adj_finemapping/genotype/spg.20k'
ROOTS = {'SuSiE-adj': B/'SuSiE_adj_finemapping/loco/results', 'SuSiE-RSS': W/'results'}
CLASSES = ['Causal', 'Non-causal in effect windows', 'Null windows']
METHODS = ['BOLT-LMM', 'REGENIE-LOCO']

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''): h.update(block)
    return h.hexdigest()

# Read packed BED directly; sum four A1 dosages per byte without expanding N x M.
bim = pd.read_csv(str(GENO)+'.bim', sep=r'\s+', header=None,
                  names=['chr','SNP','cm','position','A1','A2'])
n = sum(1 for _ in open(str(GENO)+'.fam'))
assert n == 20000 and n % 4 == 0 and bim.SNP.is_unique
with open(str(GENO)+'.bed','rb') as f: assert f.read(3) == b'\x6c\x1b\x01'
raw = np.memmap(str(GENO)+'.bed', mode='r', dtype=np.uint8, offset=3).reshape(len(bim), n//4)
codes = (np.arange(256)[:,None] >> np.array([0,2,4,6])) & 3
sums = np.array([2,0,1,0])[codes].sum(axis=1).astype(np.uint8)
missing = (codes == 1).sum(axis=1).astype(np.uint8)
a1count = np.zeros(len(bim), dtype=np.int64)
for start in range(0,len(bim),2048):
    chunk = raw[start:start+2048]
    assert not missing[chunk].any(), 'Unexpected missing genotypes'
    a1count[start:start+len(chunk)] = sums[chunk].sum(axis=1)
bim['A1_frequency'] = a1count/(2*n)
bim['MAF'] = np.minimum(bim.A1_frequency,1-bim.A1_frequency)
assert bim.MAF.between(0,0.5).all()
# Independent check against GWAS frequencies (printed at limited precision).
gw = pd.read_csv(B/'BOLTLMM/bolt_pheno1.stats',sep='\t').set_index('SNP')
f = gw.loc[bim.SNP,'A1FREQ'].to_numpy()
a = gw.loc[bim.SNP,'ALLELE1'].to_numpy()
assert np.all((a==bim.A1)|(a==bim.A2))
f = np.where(a==bim.A1,f,1-f)
assert np.max(np.abs(f-bim.A1_frequency)) < 0.00006
truth_path = B/'snp.info.chr1_18.10k.csv'
regions_path = B/'eligible_regions.csv'
truth = set(pd.read_csv(truth_path).SNP)
assert len(truth)==10000 and truth <= set(bim.SNP)
regions = pd.read_csv(regions_path).sort_values(['chr','start'])
effect = np.zeros(len(bim),dtype=bool)
for chrom, q in bim.groupby('chr'):
    r = regions[regions.chr==chrom]
    starts,ends = r.start.to_numpy(),r.end.to_numpy()
    assert (starts[1:] >= ends[:-1]).all()
    j = np.searchsorted(starts,q.position.to_numpy(),side='right')-1
    assert (j>=0).all()
    effect[q.index] = q.position.to_numpy() < ends[j]
# Same [start,end) position convention as the existing simulation evaluation.
assert np.array_equal(effect,bim.position.to_numpy()%6000000<2000000)
causal = bim.SNP.isin(truth).to_numpy()
assert effect[causal].all()
bim['class'] = np.where(causal,CLASSES[0],np.where(effect,CLASSES[1],CLASSES[2]))
bim[['SNP','chr','position','A1','A2','A1_frequency','MAF','class']].to_csv(OUT/'snp_annotations.tsv.gz',sep='\t',index=False)
print('Computed genotype MAF and three-class labels for',len(bim),'SNPs',flush=True)

rows=[]; bins=[]; calls=[]; qc=[]; pool_sets={}; source_paths=[truth_path,regions_path]
for analysis,root in ROOTS.items():
    sp=root/'snp_level_all.tsv.gz'; dp=root/'per_locus_diagnostics.tsv'
    source_paths.extend([sp,dp])
    d=pd.read_csv(dp,sep='\t')
    idcol='job_id' if 'job_id' in d else 'locus'
    good=d.loc[d.fit_status.eq('ok') & d.converged.eq(True),idcol]
    s=pd.read_csv(sp,sep='\t',usecols=['phenotype','GWAS_method','SNP','chr','position','PIP','locus_id'])
    assert s.locus_id.isin(good).all(), 'Nonconverged or unsuccessful input fit'
    s=s.merge(bim[['SNP','MAF','class']],on='SNP',how='left',validate='many_to_one')
    assert s.MAF.notna().all() and s.PIP.between(0,1).all()
    for (method,pheno),q in s.groupby(['GWAS_method','phenotype']):
        original=len(q)
        # If a SNP occurs in overlapping loci, count once and call if any PIP >= .9.
        q=q.sort_values('PIP',ascending=False).drop_duplicates('SNP').copy()
        pool_sets[(analysis,method,pheno)]=set(q.SNP)
        q['bin'],edges=pd.qcut(q.MAF,q=10,labels=False,retbins=True,duplicates='drop')
        assert q.bin.notna().all()
        q['call']=q.PIP.ge(.9)
        high=q[q.call].copy(); high.insert(0,'fine_mapping',analysis); calls.append(high)
        qc.append(dict(fine_mapping=analysis,GWAS_method=method,phenotype=pheno,
                       n_loci=q.locus_id.nunique(),pool_size=len(q),duplicate_rows_removed=original-len(q),
                       n_calls=int(q.call.sum()),n_bins=q.bin.nunique()))
        expected=dict.fromkeys(CLASSES,0.)
        for k,g in q.groupby('bin'):
            nk=int(g.call.sum())
            for cls in CLASSES:
                size=int(g['class'].eq(cls).sum()); e=nk*size/len(g)
                expected[cls]+=e
                bins.append(dict(fine_mapping=analysis,GWAS_method=method,phenotype=pheno,bin=int(k)+1,
                                 lower_edge=edges[k],upper_edge=edges[k+1],pool_n=len(g),calls_n=nk,
                                 **{'class':cls},pool_class_n=size,pool_class_fraction=size/len(g),expected=e))
        assert np.isclose(sum(expected.values()),q.call.sum())
        for cls in CLASSES:
            obs=int((q.call & q['class'].eq(cls)).sum())
            rows.append(dict(fine_mapping=analysis,GWAS_method=method,phenotype=pheno,
                             **{'class':cls},observed=obs,expected=expected[cls],
                             pool_size=len(q),n_calls=int(q.call.sum())))
    print('Counted',analysis,flush=True)
for method in METHODS:
    for p in range(1,11):
        assert pool_sets[('SuSiE-adj',method,'pheno%d'%p)]==pool_sets[('SuSiE-RSS',method,'pheno%d'%p)]
res=pd.DataFrame(rows); assert len(res)==120
res.to_csv(OUT/'counts_by_phenotype.csv',index=False)
pd.DataFrame(bins).to_csv(OUT/'maf_bin_expectations.csv',index=False)
pd.DataFrame(qc).to_csv(OUT/'pool_qc.csv',index=False)
pd.concat(calls,ignore_index=True).to_csv(OUT/'high_pip_calls.tsv.gz',sep='\t',index=False)
summary=res.groupby(['GWAS_method','fine_mapping','class'],sort=False).agg(
    observed_total=('observed','sum'),expected_total=('expected','sum'),
    observed_mean=('observed','mean'),observed_sd=('observed','std'),
    expected_mean=('expected','mean'),expected_sd=('expected','std'),n_phenotypes=('phenotype','nunique')).reset_index()
summary['observed_expected_ratio']=summary.observed_total/summary.expected_total
assert summary.n_phenotypes.eq(10).all()
summary.to_csv(OUT/'Table_2_numeric.csv',index=False)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
colors={'SuSiE-adj':'#267B9B','SuSiE-RSS':'#C66B35'}
fig,axes=plt.subplots(2,3,figsize=(11,6.6))
for i,method in enumerate(METHODS):
    for j,cls in enumerate(CLASSES):
        ax=axes[i,j]
        for x,analysis in enumerate(ROOTS):
            z=res[(res.GWAS_method==method)&(res.fine_mapping==analysis)&(res['class']==cls)]
            color=colors[analysis]
            for dx,col,hatch in [(-.18,'observed',None),(.18,'expected','///')]:
                ax.bar(x+dx,z[col].mean(),width=.32,color=color if col=='observed' else 'white',
                       edgecolor=color,hatch=hatch,linewidth=1.3,zorder=3)
                ax.errorbar(x+dx,z[col].mean(),yerr=z[col].std(),fmt='none',ecolor='#333333',
                            capsize=3,lw=1,zorder=4)
        ax.set_xticks([0,1]); ax.set_xticklabels(list(ROOTS))
        ax.set_ylim(bottom=0); ax.set_xlim(-.6,1.6)
        ax.yaxis.grid(True,color='#E5E7EB',lw=.7); ax.set_axisbelow(True)
        ax.set_title(cls if j!=1 else 'Non-causal in\neffect windows',fontsize=11)
        ax.text(-.12,1.06,'ABCDEF'[i*3+j],transform=ax.transAxes,fontweight='bold',fontsize=13)
        if j==0: ax.set_ylabel(method+'\nSNP count per phenotype',labelpad=10)
fig.suptitle('High-confidence SNPs (PIP ≥ 0.9): observed and MAF-matched expectation',fontsize=14,y=.98)
fig.legend(handles=[Patch(facecolor='#777777',label='Observed'),
                    Patch(facecolor='white',edgecolor='#777777',hatch='///',label='MAF-matched expectation')],
           loc='lower center',bbox_to_anchor=(.5,.035),ncol=2,frameon=False)
fig.text(.5,.015,'Bars: mean across 10 phenotypes; error bars: sample SD. Panel y-axis scales differ.',ha='center',fontsize=9)
fig.subplots_adjust(left=.10,right=.98,bottom=.16,top=.86,hspace=.55,wspace=.32)
for ext in ['png','pdf','svg']: fig.savefig(OUT/('Figure_2.'+ext),dpi=400,facecolor='white')
plt.close(fig)
# Compact table: totals across phenotypes; each SNP counted once within a phenotype.
table_rows=[]
for method in METHODS:
    for analysis in ROOTS:
        for cls in CLASSES:
            r=summary[(summary.GWAS_method==method)&(summary.fine_mapping==analysis)&(summary['class']==cls)].iloc[0]
            table_rows.append([method,analysis,cls,f'{r.observed_total:,}',f'{r.expected_total:,.1f}',f'{r.observed_expected_ratio:.2f}'])
cols=['GWAS method','Fine-mapping','SNP class','Observed','Expected','Obs./Exp.']
pd.DataFrame(table_rows,columns=cols).to_csv(OUT/'Table_2.csv',index=False)
fig,ax=plt.subplots(figsize=(12,5.1));ax.axis('off')
t=ax.table(cellText=table_rows,colLabels=cols,cellLoc='left',colLoc='left',loc='center',
           colWidths=[.15,.13,.30,.12,.17,.13])
t.auto_set_font_size(False);t.set_fontsize(10);t.scale(1,1.65)
for (r,c),cell in t.get_celld().items():
    cell.set_edgecolor('#D5DADF');cell.set_linewidth(.5)
    if r==0: cell.set_facecolor('#223C4B');cell.set_text_props(color='white',weight='bold')
    else: cell.set_facecolor('#F0F4F6' if ((r-1)//3)%2==0 else 'white')
ax.set_title('PIP ≥ 0.9 SNPs by class: totals across 10 phenotypes',fontsize=14,pad=18)
fig.text(.5,.035,'Expected counts use phenotype-specific MAF-matched pools; Obs./Exp. is the ratio of summed counts.',ha='center',fontsize=10)
fig.subplots_adjust(left=.02,right=.98,bottom=.10,top=.88)
for ext in ['png','pdf']:fig.savefig(OUT/('Table_2.'+ext),dpi=300,facecolor='white')
plt.close(fig)
caption=('SNPs with PIP ≥ 0.9 were classified as causal, non-causal in effect windows, or in null windows, '
         'using each SNP’s own position, without lead-window stratification or restriction to top-PIP SNPs. '
         'Causal SNPs are the 10,000 IDs in snp.info.chr1_18.10k.csv. Effect windows are the [start,end) '
         'intervals in eligible_regions.csv, following the existing simulation evaluation convention; all other '
         'positions are null windows. All causal SNPs lie in effect windows. All 13,301 fits converged in both '
         'pipelines (6,732 BOLT-LMM; 6,569 REGENIE-LOCO). Each SNP is counted once per fine-mapping approach, '
         'GWAS method, and phenotype; repeated appearances across phenotypes are counted separately. '
         'The background pool consists of all unique SNPs analyzed in that same group, including high-PIP calls. '
         'MAF was computed directly from the 20,000 GWAS animals. Each pool was divided into 10 MAF quantile '
         'bins, retaining tied MAF values together and dropping duplicate bin edges if necessary. '
         'For class c, the exact expectation is E_c = sum_k n_k p_k(c), where n_k is the number of high-PIP '
         'calls in bin k and p_k(c) is the class fraction among pool SNPs in that bin. This is the analytical '
         'mean of MAF-matched draws; no Monte Carlo approximation or significance test was used. ')
(OUT/'Figure_2_caption.txt').write_text('Figure 2. '+caption+'Bars show mean counts across 10 phenotypes, and error bars show sample SD across phenotypes for both observed and expected counts; they are not confidence intervals or random-draw uncertainty. Rows show GWAS methods; columns show SNP classes. Panel y-axis scales differ.\n')
(OUT/'Table_2_caption.txt').write_text('Table 2. '+caption+'Observed and expected counts are summed across 10 phenotypes. Expected counts are displayed to one decimal place, and observed/expected ratios use unrounded totals.\n')
(OUT/'provenance.json').write_text(json.dumps(dict(script=str(Path(__file__).resolve()),
    threshold=.9,n_gwas_animals=n,requested_maf_bins=10,lead_window_stratification=False,
    expectation='sum_k calls_in_bin_k * pool_class_count_k / pool_count_k',
    pool='unique SNPs in all successful converged loci for each fine-mapping/GWAS method/phenotype',
    source_sha256={str(p):sha(p) for p in source_paths},
    genotype_sha256={str(GENO)+s:sha(Path(str(GENO)+s)) for s in ['.bed','.bim','.fam']},
    matching_pools_between_approaches=True),indent=2)+'\n')
print(summary.to_string(index=False))
print('Saved outputs to',OUT,flush=True)
