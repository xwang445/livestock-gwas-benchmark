#!/usr/bin/env python3
"""Reproduce the manuscript's pooled SNP–phenotype precision–recall estimator.

Each SNP contributes once per phenotype, retaining max PIP across overlapping
loci. A single ranking is then made over the ten phenotypes for each GWAS input
and pipeline. This differs from the mean of ten phenotype-specific AP values.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score,precision_recall_curve,auc

p=argparse.ArgumentParser()
p.add_argument('--benchmark',type=Path,required=True)
p.add_argument('--rss-run',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
B=a.benchmark.resolve();W=a.rss_run.resolve();OUT=a.output.resolve()
OUT.mkdir(parents=True,exist_ok=True)
truth=set(pd.read_csv(B/'snp.info.chr1_18.10k.csv').SNP)
assert len(truth)==10000
roots=[('SuSiE-RSS','ERV_FALSE',B/'SuSiE_finemapping/results'),
       ('SuSiE-RSS','ERV_TRUE',W/'results'),
       ('SuSiE-adj','ERV_TRUE',B/'SuSiE_adj_finemapping/loco/results')]
rows=[];curves=[];sources={};universes={}
for pipeline,mode,root in roots:
    dp=root/'per_locus_diagnostics.tsv';sp=root/'snp_level_all.tsv.gz'
    d=pd.read_csv(dp,sep='\t')
    ident='job_id' if 'job_id' in d else 'locus' if 'locus' in d else 'locus_id'
    good=d.loc[d.fit_status.eq('ok')&d.converged.eq(True),ident]
    assert len(good)==len(d)==13301
    s=pd.read_csv(sp,sep='\t',usecols=['phenotype','GWAS_method','SNP','chr','position','PIP','locus_id'])
    assert s.locus_id.isin(good).all() and s.PIP.between(0,1).all()
    q=s.sort_values('PIP',ascending=False).drop_duplicates(['GWAS_method','phenotype','SNP']).copy()
    q['causal']=q.SNP.isin(truth)
    q['effect']=q.position.mod(6000000).lt(2000000)
    q['central']=q.position.mod(6000000).ge(3000000)&q.position.mod(6000000).lt(5000000)
    assert q.loc[q.causal,'effect'].all()
    for method,g in q.groupby('GWAS_method'):
        assert g.phenotype.nunique()==10
        key=(pipeline,mode,method)
        universes[key]=set(zip(g.phenotype,g.SNP))
        y=g.causal.to_numpy();score=g.PIP.to_numpy()
        precision,recall,threshold=precision_recall_curve(y,score)
        high=g[g.PIP.ge(.9)]
        rows.append(dict(pipeline=pipeline,fit_mode=mode,method=method,
            aggregation='pooled SNP-phenotype pairs across ten phenotypes',n_phenotypes=10,
            n_evaluated=len(g),n_causal_evaluated=int(y.sum()),
            random_precision_baseline=float(y.mean()),
            average_precision=average_precision_score(y,score),trapezoid_auprc=auc(recall,precision),
            n_high_pip=len(high),n_high_pip_causal=int(high.causal.sum()),
            high_pip_precision=float(high.causal.mean()),
            n_high_pip_effect=int(high.effect.sum()),n_high_pip_null=int((~high.effect).sum()),
            high_pip_null_fraction=float((~high.effect).mean()),
            n_high_pip_central_null=int(high.central.sum()),high_pip_central_null_fraction=float(high.central.mean()),
            n_central_null_evaluated=int(g.central.sum()),central_null_pool_fraction=float(g.central.mean()),
            mean_assigned_high_pip=float(high.PIP.mean())))
        curves.append(pd.DataFrame(dict(pipeline=pipeline,fit_mode=mode,method=method,
            precision=precision,recall=recall,threshold=np.r_[threshold,np.nan])))
    for file in [dp,sp]:sources[pipeline+'_'+mode+'_'+file.name]=hashlib.sha256(file.read_bytes()).hexdigest()
for method in ['BOLT-LMM','REGENIE-LOCO']:
    assert universes[('SuSiE-RSS','ERV_FALSE',method)]==universes[('SuSiE-RSS','ERV_TRUE',method)]
    assert universes[('SuSiE-RSS','ERV_TRUE',method)]==universes[('SuSiE-adj','ERV_TRUE',method)]
M=pd.DataFrame(rows)
M.to_csv(OUT/'pooled_finemapping_metrics.tsv',sep='\t',index=False)
pd.concat(curves,ignore_index=True).to_csv(OUT/'pooled_precision_recall_curves.tsv.gz',sep='\t',index=False)
(OUT/'pooled_finemapping_manifest.json').write_text(json.dumps(dict(
    status='complete',primary_RSS='ERV_TRUE',historical_RSS_sensitivity='ERV_FALSE',
    adjusted_pipeline='unchanged LOCO SuSiE-adj ERV_TRUE',n_locus_fits=13301,
    metric='average precision over one pooled ranking of unique SNP-phenotype pairs',
    source_sha256=sources,
    causal_truth_sha256=hashlib.sha256((B/'snp.info.chr1_18.10k.csv').read_bytes()).hexdigest(),
    identical_analyzed_pools=True,
    definitions={
        'random_precision_baseline':'n_causal_evaluated/n_evaluated',
        'high_pip_precision':'n_high_pip_causal/n_high_pip',
        'high_pip_null_fraction':'n_high_pip_null/n_high_pip',
        'high_pip_central_null_fraction':'n_high_pip_central_null/n_high_pip',
        'central_null_pool_fraction':'n_central_null_evaluated/n_evaluated',
        'average_precision':'stepwise PR area; not a count ratio and not mean phenotype AP',
        'central_null':'position modulo6Mb in [3Mb,5Mb), nested in null [2Mb,6Mb)'}),indent=2)+'\n')
print(M.to_string(index=False))
