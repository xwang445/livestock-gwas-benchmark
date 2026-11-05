#!/usr/bin/env python3
"""Extract final full-GRM REML variances from completed MLMA logs only.

GCTA MLMA prints rounded iteration variances, not a separate unrounded .hsq.
The derived h2 is explicitly a ratio of those printed estimates; no SE is inferred.
"""
import argparse, hashlib, json, re
from pathlib import Path
import numpy as np
import pandas as pd

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for s in iter(lambda:f.read(1048576),b''):h.update(s)
    return h.hexdigest()

def derive(directory,tables,qc):
    directory=Path(directory);tables=Path(tables);qc=Path(qc)
    tables.mkdir(parents=True,exist_ok=True);qc.mkdir(parents=True,exist_ok=True)
    records=[]
    for ph in range(1,11):
        log=directory/('spg.even.100k.pheno%d.log'%ph)
        mlma=directory/('spg.even.100k.pheno%d.mlma'%ph)
        r=dict(configuration='GCTA_MLMA',phenotype=ph,status='pending_completion',
               Vg_printed=np.nan,Ve_printed=np.nan,h2_from_printed_variances=np.nan,
               final_iteration=np.nan,n=20000,log_sha256='',mlma_sha256='',
               precision_note='Derived from final REML iteration printed to 5 decimal places; no unrounded estimate or SE is available')
        s=log.read_text() if log.exists() else ''
        if mlma.exists() and 'Analysis finished' in s:
            assert '20000 individuals are in common' in s,(ph,'sample correspondence not confirmed')
            assert 'Log-likelihood ratio converged.' in s,(ph,'REML not converged')
            section=s.split('Iter.\tlogL\tV(G)\tV(e)',1)
            assert len(section)==2,(ph,'variance iteration header missing')
            tail=section[1].split('Log-likelihood ratio converged.',1)[0]
            iterations=[]
            for line in tail.splitlines():
                m=re.match(r'^\s*(\d+)\s+([-+\deE.]+)\s+([-+\deE.]+)\s+([-+\deE.]+)\s*$',line)
                if m:iterations.append((int(m[1]),float(m[2]),float(m[3]),float(m[4])))
            assert iterations,(ph,'final variances not found')
            iteration,logl,vg,ve=iterations[-1]
            assert np.isfinite([logl,vg,ve]).all() and vg>=0 and ve>0
            with mlma.open() as f:
                assert next(f).split()==['Chr','SNP','bp','A1','A2','Freq','b','se','p']
                assert sum(1 for line in f if line.strip())==135358
            r.update(status='completed_log_rounded_estimate',Vg_printed=vg,Ve_printed=ve,
                     h2_from_printed_variances=vg/(vg+ve),final_iteration=iteration,
                     log_sha256=sha(log),mlma_sha256=sha(mlma))
        records.append(r)
    frame=pd.DataFrame(records);valid=frame.h2_from_printed_variances.dropna()
    complete=len(valid)==10
    summary=dict(configuration='GCTA_MLMA',n_completed=len(valid),n_required=10,
                 complete_ten_phenotypes=complete,
                 h2_mean=valid.mean() if complete else np.nan,
                 h2_sample_sd=valid.std(ddof=1) if complete else np.nan,
                 h2_min=valid.min() if complete else np.nan,
                 h2_max=valid.max() if complete else np.nan,
                 precision_note='Ratios of rounded printed Vg and Ve; aggregate withheld until all ten MLMA analyses complete')
    frame.to_csv(tables/'gcta_heritability_by_phenotype.tsv',sep='\t',index=False,na_rep='NA',float_format='%.17g')
    pd.DataFrame([summary]).to_csv(tables/'gcta_heritability_summary.tsv',sep='\t',index=False,na_rep='NA',float_format='%.17g')
    result=dict(status='PASS',completed_logs=len(valid),required=10,aggregate_available=complete,derivation_script_sha256=sha(Path(__file__)),
                source='completed GCTA MLMA final REML iteration',input_paths_published=False,
                precision='printed variances, five decimal places; no confidence interval or unrounded precision inferred')
    (qc/'gcta_heritability_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--gcta-directory',type=Path,required=True);p.add_argument('--tables',type=Path,required=True);p.add_argument('--qc',type=Path,required=True);a=p.parse_args()
    print(json.dumps(derive(a.gcta_directory,a.tables,a.qc),indent=2))
