#!/usr/bin/env python3
"""Read-only independent recount of saved converged fits; stdout is aggregates.

No fit, GWAS, genotype, LD or original output is modified. Does not validate the
SuSiE optimizer itself. Input paths are supplied at runtime and omitted from JSON.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()

def verify(root):
    root=Path(root);sources={};rows=[];fits=[]
    truth_path=root/'snp.info.chr1_18.10k.csv';truth=set(pd.read_csv(truth_path).SNP)
    if len(truth)!=10000:raise ValueError('Truth set size')
    annpath=root/'SuSiE_finemapping/publication/maf_matched/snp_annotations.tsv.gz'
    ann=pd.read_csv(annpath,sep='\t',usecols=['SNP','chr','position','MAF'])
    if ann.SNP.duplicated().any() or not ann.MAF.between(0,.5).all():raise ValueError('Invalid annotations')
    sources['truth']=sha(truth_path);sources['saved_MAF_annotations']=sha(annpath)
    for pipeline,folder in [('SuSiE-adj','SuSiE_adj_finemapping/loco/results'),('SuSiE-RSS','SuSiE_finemapping/results')]:
        dpath=root/folder/'per_locus_diagnostics.tsv';spath=root/folder/'snp_level_all.tsv.gz'
        diagnostics=pd.read_csv(dpath,sep='\t');idcol='job_id' if 'job_id' in diagnostics else 'locus'
        good=diagnostics.loc[diagnostics.fit_status.eq('ok')&diagnostics.converged.eq(True),idcol]
        s=pd.read_csv(spath,sep='\t',usecols=['phenotype','GWAS_method','SNP','chr','position','PIP','locus_id'])
        if not s.locus_id.isin(good).all() or not s.PIP.between(0,1).all():raise ValueError('Invalid fit provenance/PIP')
        s=s.merge(ann,on=['SNP','chr','position'],how='left',validate='many_to_one')
        if s.MAF.isna().any():raise ValueError('Unmatched MAF/coordinate')
        for (method,ph),g in s.groupby(['GWAS_method','phenotype']):
            q=g.sort_values('PIP',ascending=False).drop_duplicates('SNP').copy()
            q['causal']=q.SNP.isin(truth);q['effect']=q.position.mod(6000000).lt(2000000)
            if not q.loc[q.causal,'effect'].all():raise ValueError('Causal SNP in null')
            q['bin']=pd.qcut(q.MAF,10,labels=False,duplicates='drop');q['call']=q.PIP.ge(.9)
            high=q[q.call];expected=0.
            for _,b in q.groupby('bin'):expected+=int(b.call.sum())*int(b.causal.sum())/len(b)
            rows.append(dict(method=method.replace('REGENIE-LOCO','REGENIE_LOCO'),pipeline=pipeline,
                phenotype=int(str(ph).replace('pheno','')),n_high_pip=len(high),n_causal=int(high.causal.sum()),
                n_null=int((~high.effect).sum()),n_effect=int(high.effect.sum()),expected_causal=expected,
                n_pool=len(q),n_maf_bins=int(q.bin.nunique()),n_converged_loci=int(g.locus_id.nunique())))
        for method,g in s.groupby('GWAS_method'):fits.append(dict(method=method,pipeline=pipeline,n_converged_loci=int(g.locus_id.nunique())))
        sources[pipeline+'_diagnostics']=sha(dpath);sources[pipeline+'_PIPs']=sha(spath)
    per=pd.DataFrame(rows)
    if len(per)!=40 or not per.groupby(['method','pipeline']).phenotype.nunique().eq(10).all():raise ValueError('Incomplete replicates')
    pooled=per.groupby(['method','pipeline'])[['n_high_pip','n_causal','n_null','n_effect','expected_causal']].sum().reset_index()
    pooled['null_percent']=100*pooled.n_null/pooled.n_high_pip
    pooled['causal_percent']=100*pooled.n_causal/pooled.n_high_pip
    pooled['noncausal_percent']=100-pooled.causal_percent
    pooled['maf_expected_percent']=100*pooled.expected_causal/pooled.n_high_pip
    pooled['enrichment']=pooled.n_causal/pooled.expected_causal
    return dict(status='SAVED_FIT_AGGREGATES_RECOMPUTED',per_phenotype=per.to_dict('records'),pooled=pooled.to_dict('records'),fits=fits,input_sha256=sources,
        limitations=['No GWAS or SuSiE refitting','MAF annotations recovered; genotype MAF calculation not rerun by this independent recount'])
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);args=ap.parse_args()
    print(json.dumps(verify(args.root),indent=2,allow_nan=False))
