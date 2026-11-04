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
for p in (ROOT/'derived_data').rglob('*.tsv'):
    d=pd.read_csv(p,sep='\t')
    assert not {'FID','IID','snp','SNP','animal_id','genotype','beta','se','PIP','pip'}.intersection(d.columns),str(p)
    for col in ['method','method_a','method_b']:
        if col in d and 'template' not in p.name:assert set(d[col].dropna()).issubset(METHODS),(str(p),col)
    if 'template' not in p.name:assert len(d)>0,str(p)
rfiles=[str(p) for p in (ROOT/'scripts').rglob('*.R')]
subprocess.run(['Rscript','-e','for (f in commandArgs(TRUE)) { parse(file=f); }; cat("R parse checks passed\\n")',*rfiles],check=True)
print('Python compilation, shell syntax, YAML, R parsing, aggregate schemas and method labels passed.')
