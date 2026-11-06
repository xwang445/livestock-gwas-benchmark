#!/usr/bin/env python3
"""Generate the verified-result configuration from committed executed tables."""
import argparse, hashlib,json
from pathlib import Path
import pandas as pd
import yaml

def build(repo):
    repo=Path(repo);data=repo/'derived_data'
    def read(name):return pd.read_csv(data/name,sep='\t')
    sources={}
    for path in sorted(list(data.rglob('*.tsv'))+list(data.rglob('*.csv'))):
        frame=pd.read_csv(path,sep=',' if path.suffix=='.csv' else '\t')
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
    result={'schema_version':2,
        'scope':'Primary blocks derive from the current inclusive revision tables refreshed November 6, updated SLEMM, ERV_TRUE RSS and actual partial GCTA coverage. Explicit historical_* blocks preserve the November 4 strict reconstruction. source_tables fingerprints both separate datasets; historical and revised results must not be mixed.',
        'regenerate_command':'python scripts/build_result_config.py',
        'historical_study_checks':{'historical_configurations':int(counts.method.nunique()),'phenotypes':int(counts.phenotype.nunique()),
            'gwas_panels':len(counts),'snps_per_panel':sorted(int(x) for x in counts.n_snps.unique()),
            'equation13_replicate_bin_rows':sources['equation13_summary/equation13_by_pheno.tsv']['rows'],
            'finemapping_phenotype_pipeline_rows':sources['finemapping_summary/high_pip_by_pheno.tsv']['rows']},
        'historical_validation':{'ledger':'audit/manuscript_validation.tsv',
            'status_counts':{k:int(v) for k,v in ledger.status.value_counts().sort_index().items()},
            'manuscript_discrepancies':ledger[ledger.status=='FAIL'][['manuscript_item','expected','observed','source','note']].to_dict('records'),
            'historical_discrepancy_ledger':'audit/historical_manuscript_validation.tsv',
            'corrected_manuscript_values':{row.manuscript_item:float(row.observed) for row in ledger.itertuples() if row.manuscript_item in ['Results REGENIE maximum SNP count','Figure4 corrected maximum minus log10 P']}},
        'historical_methods':per_method,
        'historical_heritability_from_saved_logs':h.to_dict('records'),
        'historical_finemapping':{'pooled':read('finemapping_summary/high_pip_summary.tsv').to_dict('records'),
            'converged_loci':read('finemapping_summary/converged_loci.tsv').to_dict('records')},
        'ld':{'n_bins':len(ld),'selected_bin_means':ld[ld.bin_id.isin([1,200,1000])].to_dict('records')},
        'historical_equation13_cached_checks':read('equation13_summary/cached_theory_validation.tsv').to_dict('records'),
        'source_tables':sources}
    revision=data/'revision_2026-11-05/association/association_by_phenotype.tsv'
    if revision.exists():
        registry=pd.read_csv(repo/'config/revision_methods.tsv',sep='\t')
        panels=pd.read_csv(revision,sep='\t')
        panels=panels[panels.configuration.isin(registry.configuration)]
        observed=panels.groupby('configuration').phenotype.nunique().to_dict()
        result['revision_checks']={
            'registry':'config/revision_methods.tsv',
            'planned_method_families':int(registry.method_family.nunique()),
            'planned_configurations':len(registry),
            'planned_panels':int(registry.planned_phenotypes.sum()),
            'completed_panels':len(panels),
            'completed_phenotypes_by_configuration':{m:int(observed.get(m,0)) for m in registry.configuration},
            'fully_completed_configurations':sum(observed.get(r.configuration,0)==r.planned_phenotypes for r in registry.itertuples()),
            'association_source':'derived_data/revision_2026-11-05/association/association_by_phenotype.tsv',
            'fine_mapping_source':'derived_data/revision_2026-11-05/susie/sensitivity_comparison_table.tsv',
            'status_document':'docs/revision_2026-11-05.md'}
        result['significance']={'alpha':5e-8,'rule':'P <= 5e-8','native_log10p_rule':'LOG10P >= -log10(5e-8)','historical_profile':'explicit --historical only; original snapshots preserved'}
        result['study_checks']={key:result['revision_checks'][key] for key in ['planned_method_families','planned_configurations','planned_panels','completed_panels','fully_completed_configurations']}
        current_summary=read('revision_2026-11-05/association/association_summary.tsv')
        current_summary=current_summary[current_summary.configuration.isin(registry.configuration)]
        result['methods']={row.pop('configuration'):row for row in current_summary.to_dict('records')}
        result['finemapping']={
            'rss_configuration':'estimate_residual_variance=TRUE; original loci unchanged',
            'pooled':pd.read_csv(data/'revision_2026-11-05/susie/localization_precision/Table_3_numeric.csv').to_dict('records'),
            'pooled_ranking_metrics':read('revision_2026-11-05/susie/pooled_finemapping_metrics.tsv').to_dict('records'),
            'paired_sensitivity':read('revision_2026-11-05/susie/sensitivity_comparison_table.tsv').to_dict('records')}
        result['gcta_heritability']={
            'precision':'ratio of five-decimal printed final REML variances; no unrounded estimate/SE inferred',
            'by_phenotype':read('revision_2026-11-05/association/gcta_heritability_by_phenotype.tsv').to_dict('records'),
            'summary':read('revision_2026-11-05/association/gcta_heritability_summary.tsv').to_dict('records')}
        result['validation']={'record':'audit/revision_validation.json',
            'result':json.loads((repo/'audit/revision_validation.json').read_text()),
            'scope':'Current inclusive numerical revision; the older manuscript ledger is historical_validation. The full ten-phenotype GCTA comparison and upstream Equation12/13 cache verification remain incomplete.'}
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
