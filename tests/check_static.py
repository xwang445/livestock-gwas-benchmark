#!/usr/bin/env python3
"""Static syntax, configuration and public aggregate privacy/schema checks."""
from pathlib import Path
import ast, subprocess, sys
import pandas as pd
import yaml
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from benchmark_core import METHODS
for p in sorted(p for folder in ['scripts','tests','example','finemapping'] for p in (ROOT/folder).rglob('*')):
    if p.suffix=='.py':compile(p.read_text(),str(p),'exec')
    if p.suffix in ['.sh','.sbatch','.bash']:subprocess.run(['bash','-n',str(p)],check=True)
for p in ROOT.rglob('*.yaml'):yaml.safe_load(p.read_text())
for p in ROOT.rglob('*.yml'):yaml.safe_load(p.read_text())
config=yaml.safe_load((ROOT/'config/methods.yaml').read_text())
revision_methods=set(pd.read_csv(ROOT/'config/revision_methods.tsv',sep='\t').configuration)
assert set(config['methods'])==revision_methods
params=yaml.safe_load((ROOT/'config/analysis_parameters.yaml').read_text())
assert params['simulation']['causal_snps']==10000 and params['windows']['central_null_locus_denominator']==1480
assert params['significance']['alpha']==5e-8
assert all(params['significance'][key]=='le' for key in ['figure3_comparator','spillover_comparator','fine_mapping_locus_comparator','locus_overlap_comparator','manhattan_comparator'])
required_results={
    'method_summary/method_summary.tsv',
    'power_summary/power_summary.tsv',
    'power_summary/fpr_summary.tsv',
    'spillover_summary/spillover_summary.tsv',
    'lambda_gc_summary/lambda_gc_summary.tsv',
    'overlap_summary/locus_overlap_summary.tsv',
    'equation13_summary/equation13_summary.tsv',
    'finemapping_summary/high_pip_summary.tsv',
    'ld_summary/ld_decay_10kb.tsv',
}
for name in required_results:assert (ROOT/'derived_data'/name).is_file(), 'Missing executed result: '+name
def current_table(name):
    frame=pd.read_csv(ROOT/'derived_data'/name,sep='\t')
    for column in ['method','method_a','method_b']:
        if column in frame:
            frame[column]=frame[column].replace({'REGENIE_woLOCO':'REGENIE_WO_LOCO'})
    return frame

# Every applicable top-level summary must expose the completed comparison.
# Historical audits have their own directory and cannot satisfy these checks.
for folder in ['method_summary','lambda_gc_summary','spillover_summary','power_summary','overlap_summary']:
    for path in (ROOT/'derived_data'/folder).glob('*.tsv'):
        if path.name=='heritability_from_logs.tsv':continue
        frame=current_table(str(path.relative_to(ROOT/'derived_data')))
        for column in ['method','method_a','method_b']:
            if column in frame:
                assert set(frame[column])==revision_methods, 'Current root table omits a configuration: '+str(path)
counts=current_table('method_summary/significant_snps_by_pheno.tsv')
inflation=current_table('lambda_gc_summary/lambda_gc_by_pheno.tsv')
canonical=pd.read_csv(ROOT/'derived_data/revision_2026-11-05/association/association_by_phenotype.tsv',sep='\t')
canonical=canonical[canonical.configuration.isin(revision_methods)].rename(columns={'configuration':'method'})
keys=['method','phenotype']
for frame in [counts,inflation]:
    assert len(frame)==100 and not frame.duplicated(keys).any(), 'Expected all 100 current root analysis panels'
    assert frame.groupby('method').phenotype.nunique().eq(10).all()
count_check=counts.merge(canonical,on=keys,validate='one_to_one')
assert count_check.n_significant_snps.eq(count_check.significant_snps).all()
assert count_check.n_significant_loci.eq(count_check.loci).all()
lambda_check=inflation.merge(canonical,on=keys,validate='one_to_one')
assert pd.notna(lambda_check.lambda_gc).all()
assert ((lambda_check.lambda_gc-lambda_check.lambda_gc_central_null_historical).abs()<1e-12).all(), 'Current root lambda differs from canonical fit statistics'
assert ((lambda_check.lambda_gc_p_based-lambda_check.lambda_gc_central_null).abs()<1e-12).all(), 'Stable lambda diagnostic differs from canonical P-based statistics'
for name,rows in [('locus_overlap_by_pheno.tsv',1000),('locus_overlap_summary.tsv',100),('cross_replicate_stability.tsv',1000)]:
    assert len(current_table('overlap_summary/'+name))==rows, 'Incomplete current root overlap coverage: '+name
allowed_methods=set(METHODS)|revision_methods|{'SLEMM_previous','REGENIE-LOCO'}
s1=pd.read_csv(ROOT/'derived_data/revision_2026-11-05/association/supplementary_table_S1.tsv',sep='\t')
assert len(s1)==len(revision_methods) and set(s1.configuration)==revision_methods, 'S1 must contain each current configuration exactly once'
summary=pd.read_csv(ROOT/'derived_data/revision_2026-11-05/association/association_summary.tsv',sep='\t')
expected_s1=summary[summary.configuration.isin(revision_methods)][s1.columns].reset_index(drop=True)
pd.testing.assert_frame_equal(s1,expected_s1,check_exact=False,rtol=1e-12,atol=1e-12)
for p in sorted(list((ROOT/'derived_data').rglob('*.tsv'))+list((ROOT/'derived_data').rglob('*.csv'))):
    assert 'template' not in p.name, 'Keep schema templates under example/schemas: '+str(p)
    d=pd.read_csv(p,sep=',' if p.suffix=='.csv' else '\t',keep_default_na=False,na_values=['NA',''])
    assert not {'FID','IID','snp','SNP','animal_id','genotype','beta','se','PIP','pip'}.intersection(d.columns),str(p)
    for col in ['method','method_a','method_b','configuration','configuration_a','configuration_b','GWAS_method']:
        if col in d:assert set(d[col].dropna()).issubset(allowed_methods),(str(p),col)
    assert len(d)>0,str(p)
    numeric=d.select_dtypes(include='number')
    if p.name=='numerical_claims_reconciliation.tsv':
        assert 'validation_status' in d or 'status' in d, 'Reconciliation requires validation status'
    else:
        assert not numeric.empty and numeric.notna().any().any(), 'Missing executed numeric results: '+str(p)
rfiles=[str(p) for p in (ROOT/'scripts').rglob('*.R')]
subprocess.run(['Rscript','-e','for (f in commandArgs(TRUE)) { parse(file=f); }; cat("R parse checks passed\\n")',*rfiles],check=True)
print('Python compilation, shell syntax, YAML, R parsing, aggregate schemas and method labels passed.')
