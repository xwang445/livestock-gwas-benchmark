#!/usr/bin/env python3
"""Generate the verified-result configuration from committed executed tables."""
import argparse, hashlib,json
from pathlib import Path
import pandas as pd
import yaml

def build(repo):
    repo=Path(repo);data=repo/'derived_data'
    def read(name):return pd.read_csv(data/name,sep='\t')
    historical='historical_2026-11-04'
    def read_historical(name):return read(historical+'/'+name)
    sources={}
    for path in sorted(list(data.rglob('*.tsv'))+list(data.rglob('*.csv'))):
        frame=pd.read_csv(path,sep=',' if path.suffix=='.csv' else '\t')
        sources[str(path.relative_to(data))]={'path':str(path.relative_to(repo)),
            'rows':len(frame),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    counts=read_historical('method_summary/significant_snps_by_pheno.tsv')
    methods=read_historical('method_summary/method_summary.tsv')
    spill=read_historical('spillover_summary/spillover_summary.tsv').set_index('method').to_dict(orient='index')
    inflation=read_historical('lambda_gc_summary/lambda_gc_summary.tsv').set_index('method').to_dict(orient='index')
    per_method={}
    for row in methods.to_dict('records'):
        method=row.pop('method')
        per_method[method]={'significant_snps_and_loci':row,
            'spillover':spill[method],
            'lambda_gc':inflation[method]}
    h=read_historical('method_summary/heritability_from_logs.tsv')
    h=h.groupby('method').h2.agg(['mean','min','max','count']).reset_index()
    ledger=pd.read_csv(repo/'audit/manuscript_validation.tsv',sep='\t')
    ld=read('ld_summary/ld_decay_10kb.tsv')
    result={'schema_version':2,
        'scope':'Primary blocks derive from the inclusive revision tables with updated SLEMM, ERV_TRUE RSS and all ten completed GCTA phenotypes. November 11 recovered-effect results and November 12 submission aggregates are separately indexed. Explicit historical_* blocks preserve the November 4 strict reconstruction. source_tables fingerprints each dataset; historical and revised estimators must not be mixed.',
        'regenerate_command':'python scripts/build_result_config.py',
        'historical_study_checks':{'historical_configurations':int(counts.method.nunique()),'phenotypes':int(counts.phenotype.nunique()),
            'gwas_panels':len(counts),'snps_per_panel':sorted(int(x) for x in counts.n_snps.unique()),
            'data_directory':'derived_data/'+historical,
            'equation13_replicate_bin_rows':sources[historical+'/equation13_summary/equation13_by_pheno.tsv']['rows'],
            'finemapping_phenotype_pipeline_rows':sources[historical+'/finemapping_summary/high_pip_by_pheno.tsv']['rows']},
        'historical_validation':{'ledger':'audit/manuscript_validation.tsv',
            'status_counts':{k:int(v) for k,v in ledger.status.value_counts().sort_index().items()},
            'manuscript_discrepancies':ledger[ledger.status=='FAIL'][['manuscript_item','expected','observed','source','note']].to_dict('records'),
            'historical_discrepancy_ledger':'audit/historical_manuscript_validation.tsv',
            'corrected_manuscript_values':{row.manuscript_item:float(row.observed) for row in ledger.itertuples() if row.manuscript_item in ['Results REGENIE maximum SNP count','Figure4 corrected maximum minus log10 P']}},
        'historical_methods':per_method,
        'historical_heritability_from_saved_logs':h.to_dict('records'),
        'historical_finemapping':{'pooled':read_historical('finemapping_summary/high_pip_summary.tsv').to_dict('records'),
            'converged_loci':read_historical('finemapping_summary/converged_loci.tsv').to_dict('records')},
        'ld':{'n_bins':len(ld),'selected_bin_means':ld[ld.bin_id.isin([1,200,1000])].to_dict('records')},
        'historical_equation13_cached_checks':read_historical('equation13_summary/cached_theory_validation.tsv').to_dict('records'),
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
        root_counts=read('method_summary/significant_snps_by_pheno.tsv')
        root_methods=root_counts.method.replace({'REGENIE_woLOCO':'REGENIE_WO_LOCO'})
        result['current_summary_exports']={
            'data_directory':'derived_data',
            'run_manifest':'derived_data/run_manifest.json',
            'validation_record':'audit/current_summary_update_2026-11-13.json',
            'completed_panels':len(root_counts),
            'completed_configurations':int(root_methods.nunique()),
            'completed_phenotypes_by_configuration':{m:int(n) for m,n in root_counts.assign(configuration=root_methods).groupby('configuration').phenotype.nunique().items()},
            'association_rule':'Inclusive P <= 5e-8; current updated SLEMM and all ten GCTA phenotypes',
            'historical_data_directory':'derived_data/'+historical,
            'scope':'Root method, lambda, power, spillover and overlap tables are current exports. Fine-mapping and cached theory retain only the methods with executed evidence.'}
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
        submission=data/'submission_2026-11-12'
        if submission.is_dir():
            result['submission']={
                'snapshot_date':'2026-11-12',
                'status_document':'docs/revision_2026-11-12.md',
                'artwork':'figures/submission_2026-11-12',
                'figure8_pooled_metrics':read('submission_2026-11-12/figure08_pooled_metrics.tsv').to_dict('records'),
                'figure8_recovered_pve_groups':read('submission_2026-11-12/figure08_pve_group_summary.tsv').to_dict('records'),
                'table_S1_aggregation':'Ratios of pooled detected-locus counts across all ten phenotypes; central-null is a subset of null',
                'table_S1':read('submission_2026-11-12/tables/Table_S1.tsv').to_dict('records')}
        recovered=data/'revision_2026-11-11/recovered_effects/recovered_lead_pve_summary.tsv'
        if recovered.is_file():
            result['recovered_effects']={
                'definition':'Author-supplied recovered-effect PVE; original generating effects and recovery algorithm unverified',
                'status_document':'docs/recovered_effects_2026-11-11.md',
                'lead_summary':read('revision_2026-11-11/recovered_effects/recovered_lead_pve_summary.tsv').to_dict('records')}
        result['validation']={'record':'audit/revision_validation.json',
            'result':json.loads((repo/'audit/revision_validation.json').read_text()),
            'scope':'Current inclusive numerical revision; the older manuscript ledger is historical_validation. The ten-phenotype GCTA comparison is complete; upstream Equation12/13 cache verification remains incomplete.'}
        publication_validation=repo/'audit/publication_update_2026-11-13.json'
        if publication_validation.is_file():
            result['validation']['publication_update_record']='audit/publication_update_2026-11-13.json'
            result['validation']['publication_update']=json.loads(publication_validation.read_text())
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
