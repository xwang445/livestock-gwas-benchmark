#!/usr/bin/env python3
"""Static syntax, configuration and public aggregate privacy/schema checks."""
from pathlib import Path
import ast, subprocess, sys
import pandas as pd
import yaml
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from benchmark_core import METHODS
for p in sorted(p for folder in ['scripts','tests','example'] for p in (ROOT/folder).rglob('*')):
    if p.suffix=='.py':compile(p.read_text(),str(p),'exec')
    if p.suffix in ['.sh','.sbatch','.bash']:subprocess.run(['bash','-n',str(p)],check=True)
for p in ROOT.rglob('*.yaml'):yaml.safe_load(p.read_text())
for p in ROOT.rglob('*.yml'):yaml.safe_load(p.read_text())
config=yaml.safe_load((ROOT/'config/methods.yaml').read_text())
assert set(config['methods'])==set(METHODS)
params=yaml.safe_load((ROOT/'config/analysis_parameters.yaml').read_text())
assert params['simulation']['causal_snps']==10000 and params['windows']['central_null_locus_denominator']==1480
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
for p in (ROOT/'derived_data').rglob('*.tsv'):
    assert 'template' not in p.name, 'Keep schema templates under example/schemas: '+str(p)
    d=pd.read_csv(p,sep='\t')
    assert not {'FID','IID','snp','SNP','animal_id','genotype','beta','se','PIP','pip'}.intersection(d.columns),str(p)
    for col in ['method','method_a','method_b']:
        if col in d:assert set(d[col].dropna()).issubset(METHODS),(str(p),col)
    assert len(d)>0,str(p)
    numeric=d.select_dtypes(include='number')
    assert not numeric.empty and numeric.notna().any().any(), 'Missing executed numeric results: '+str(p)
rfiles=[str(p) for p in (ROOT/'scripts').rglob('*.R')]
subprocess.run(['Rscript','-e','for (f in commandArgs(TRUE)) { parse(file=f); }; cat("R parse checks passed\\n")',*rfiles],check=True)
print('Python compilation, shell syntax, YAML, R parsing, aggregate schemas and method labels passed.')
