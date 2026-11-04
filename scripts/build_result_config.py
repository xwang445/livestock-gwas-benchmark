#!/usr/bin/env python3
"""Generate the verified-result configuration from committed executed tables."""
import argparse, hashlib
from pathlib import Path
import pandas as pd
import yaml

def build(repo):
    repo=Path(repo);data=repo/'derived_data'
    def read(name):return pd.read_csv(data/name,sep='\t')
    sources={}
    for path in sorted(data.rglob('*.tsv')):
        frame=pd.read_csv(path,sep='\t')
        sources[str(path.relative_to(data))]={'path':str(path.relative_to(repo)),
            'rows':len(frame),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    counts=read('method_summary/significant_snps_by_pheno.tsv')
    methods=read('method_summary/method_summary.tsv')
    spill=read('spillover_summary/spillover_summary.tsv').set_index('method').to_dict(orient='index')
    inflation=read('lambda_gc_summary/lambda_gc_summary.tsv').set_index('method').to_dict(orient='index')
    per_method={}
    for row in methods.to_dict('records'):
        method=row.pop('method')
        per_method[method]={'significant_snps_and_loci':row,
            'spillover':spill[method],
            'lambda_gc':inflation[method]}
    h=read('method_summary/heritability_from_logs.tsv')
    h=h.groupby('method').h2.agg(['mean','min','max','count']).reset_index()
    ledger=pd.read_csv(repo/'audit/manuscript_validation.tsv',sep='\t')
    ld=read('ld_summary/ld_decay_10kb.tsv')
    result={'schema_version':1,
        'scope':'Executed downstream reanalysis of saved GWAS, cached theory, LD and converged fine-mapping results; no fresh model fits.',
        'regenerate_command':'python scripts/build_result_config.py',
        'study_checks':{'methods':int(counts.method.nunique()),'phenotypes':int(counts.phenotype.nunique()),
            'gwas_panels':len(counts),'snps_per_panel':sorted(int(x) for x in counts.n_snps.unique()),
            'equation13_replicate_bin_rows':sources['equation13_summary/equation13_by_pheno.tsv']['rows'],
            'finemapping_phenotype_pipeline_rows':sources['finemapping_summary/high_pip_by_pheno.tsv']['rows']},
        'validation':{'ledger':'audit/manuscript_validation.tsv',
            'status_counts':{k:int(v) for k,v in ledger.status.value_counts().sort_index().items()},
            'manuscript_discrepancies':ledger[ledger.status=='FAIL'][['manuscript_item','expected','observed','source','note']].to_dict('records')},
        'methods':per_method,
        'heritability_from_saved_logs':h.to_dict('records'),
        'finemapping':{'pooled':read('finemapping_summary/high_pip_summary.tsv').to_dict('records'),
            'converged_loci':read('finemapping_summary/converged_loci.tsv').to_dict('records')},
        'ld':{'n_bins':len(ld),'selected_bin_means':ld[ld.bin_id.isin([1,200,1000])].to_dict('records')},
        'equation13_cached_checks':read('equation13_summary/cached_theory_validation.tsv').to_dict('records'),
        'source_tables':sources}
    # Preserve parsed table precision and integer counts; normalize only scalar types.
    import math
    def normalize(value):
        if isinstance(value,dict):return {k:normalize(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):return [normalize(v) for v in value]
        if hasattr(value,'item'):value=value.item()
        return None if isinstance(value,float) and not math.isfinite(value) else value
    return normalize(result)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--check',action='store_true');a=p.parse_args()
    content=yaml.safe_dump(build(a.repo),sort_keys=False,allow_unicode=True)
    target=a.repo/'config/verified_results.yaml'
    if a.check:
        if not target.exists() or target.read_text()!=content:raise SystemExit('verified_results.yaml differs from executed result tables; regenerate it.')
        print('Verified results configuration matches all executed source tables.')
    else:
        target.write_text(content);print('Wrote',target)
