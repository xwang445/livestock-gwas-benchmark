#!/usr/bin/env python3
"""Export non-identifying LD/fine-mapping verification JSON to aggregate tables."""
import argparse,json
from pathlib import Path
import pandas as pd
from benchmark_core import output_table
p=argparse.ArgumentParser();p.add_argument('--kind',choices=['ld','finemapping'],required=True);p.add_argument('--json',required=True);p.add_argument('--out',required=True);a=p.parse_args()
d=json.loads(Path(a.json).read_text());out=Path(a.out)
if a.kind=='ld':
    assert d['status']=='LD_BIN_MEANS_RECOMPUTED'
    output_table(pd.DataFrame(d['bins']),out/'ld_decay_10kb.tsv')
else:
    assert d['status']=='SAVED_FIT_AGGREGATES_RECOMPUTED'
    for key,name in [('per_phenotype','high_pip_by_pheno.tsv'),('pooled','high_pip_summary.tsv'),('fits','converged_loci.tsv')]:
        f=pd.DataFrame(d[key]);f.method=f.method.replace({'REGENIE-LOCO':'REGENIE_LOCO'});output_table(f,out/name)
    output_table(pd.DataFrame(d['pooled'])[['method','pipeline','n_high_pip','n_null','n_effect','null_percent']],out/'null_vs_effect_summary.tsv')
    (out/'verification_manifest.json').write_text(json.dumps({k:d[k] for k in ['status','input_sha256','limitations']},indent=2)+'\n')
