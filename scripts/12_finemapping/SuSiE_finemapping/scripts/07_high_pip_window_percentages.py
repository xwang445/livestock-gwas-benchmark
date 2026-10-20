#!/usr/bin/env python3
"""Summarize saved PIP >= 0.9 calls by each SNP's own genomic window."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
W=Path(__file__).resolve().parents[1]
OUT=W/'results/high_pip_window_percentages'
OUT.mkdir(parents=True,exist_ok=True)
source=W/'publication/maf_matched/high_pip_calls.tsv.gz'
d=pd.read_csv(source,sep='\t')
truth_path=W.parent/'snp.info.chr1_18.10k.csv'
region_path=W.parent/'eligible_regions.csv'
truth=set(pd.read_csv(truth_path).SNP)
regions=pd.read_csv(region_path)
assert len(truth)==10000 and d.PIP.ge(.9).all()
assert not d.duplicated(['GWAS_method','fine_mapping','phenotype','SNP']).any()
d['effect_window']=False
for ch,indices in d.groupby('chr').groups.items():
    r=regions[regions.chr==ch].sort_values('start')
    p=d.loc[indices,'position'].to_numpy()
    ix=np.searchsorted(r.start.to_numpy(),p,side='right')-1
    valid=ix>=0
    effect=np.zeros(len(p),dtype=bool)
    effect[valid]=p[valid]<r.end.to_numpy()[ix[valid]]
    d.loc[indices,'effect_window']=effect
# Independent check against the simulation's alternating 2-Mb effect / 4-Mb null architecture.
assert np.array_equal(d.effect_window,d.position%6000000<2000000)
d['true_causal']=d.SNP.isin(truth)
assert d.loc[d.true_causal,'effect_window'].all()
assert np.array_equal(d.true_causal,d['class'].eq('Causal'))
assert np.array_equal(~d.effect_window,d['class'].eq('Null windows'))

def summarize(keys):
    rows=[]
    for k,g in d.groupby(keys):
        if not isinstance(k,tuple):k=(k,)
        n=len(g);e=int(g.effect_window.sum());t=int(g.true_causal.sum())
        rows.append(dict(zip(keys,k),n_high_pip=n,n_effect=e,n_null=n-e,
            effect_pct=100*e/n,null_pct=100*(n-e)/n,
            n_true_causal=t,n_true_causal_effect=t,n_true_causal_null=0,
            true_causal_effect_pct=100. if t else np.nan,true_causal_null_pct=0. if t else np.nan,
            n_noncausal_effect=e-t,true_causal_pct_of_all_calls=100*t/n))
    return pd.DataFrame(rows)
pooled=summarize(['GWAS_method','fine_mapping'])
per=summarize(['GWAS_method','fine_mapping','phenotype'])
assert len(pooled)==4 and len(per)==40
assert np.allclose(pooled.effect_pct+pooled.null_pct,100)
reference=pd.read_csv(W/'publication/localization_precision/Table_3_numeric.csv')
qc=pooled.merge(reference,on=['GWAS_method','fine_mapping'],validate='one_to_one',suffixes=('','_previous'))
assert (qc.n_high_pip==qc.n_calls).all() and (qc.n_true_causal==qc.n_causal).all()
assert (qc.n_null==qc.n_null_previous).all()
assert (qc.n_noncausal_effect==qc.n_noncausal_effect_previous).all()
pooled.to_csv(OUT/'pooled_by_method.csv',index=False)
per.to_csv(OUT/'by_phenotype.csv',index=False)
d.to_csv(OUT/'annotated_high_pip_calls.tsv.gz',sep='\t',index=False)
(OUT/'analysis_notes.txt').write_text(
    'High PIP means PIP >= 0.9, consistent with the existing evaluation.\n'
    'SuSiE-adj is the LOCO-adjusted pipeline used in the existing SuSiE-RSS comparison.\n'
    'Calls are from successful converged fits, using the existing saved call table.\n'
    'Pooled percentages are ratios of summed counts across ten phenotypes, not means of percentages.\n'
    'A SNP is counted once within GWAS method, fine-mapping pipeline and phenotype.\n'
    'A SNP selected in multiple phenotypes contributes once in each phenotype.\n'
    'Region membership uses each SNP position and the original [start,end) effect intervals;\n'
    'null means the complementary 4-Mb null regions, not only their central 2 Mb.\n'
    'All high-PIP calls: denominator is all selected SNP-phenotype calls within a method/pipeline.\n'
    'True-causal subset: denominator is only selected SNPs in the simulated causal list.\n'
    'Every selected true-causal SNP lies in an effect window; null-region high-PIP calls are noncausal.\n'
    'Counts independently match the existing Table_3_numeric.csv for all four comparisons.\n'
    'No GWAS or fine-mapping was rerun.\n')
(OUT/'provenance.json').write_text(json.dumps({
    'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,truth_path,region_path]},
    'threshold':.9,'n_calls':len(d),'validated_against_existing_table':True},indent=2)+'\n')
print(pooled[['GWAS_method','fine_mapping','n_high_pip','n_effect','effect_pct','n_null','null_pct','n_true_causal']].to_string(index=False))
