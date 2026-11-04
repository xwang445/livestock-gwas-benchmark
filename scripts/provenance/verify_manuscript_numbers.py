#!/usr/bin/env python3
"""Independently check manuscript counts/maxima from saved raw result columns."""
import argparse,csv,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

def calculate(root,repo):
 methods=yaml.safe_load((Path(repo)/'config/methods.yaml').read_text())['methods'];rows=[]
 for method,settings in methods.items():
  for phenotype in range(1,11):
   path=Path(root)/settings['output_patterns']['mac'].format(phenotype=phenotype)
   raw=pd.read_csv(path,sep=',' if method.startswith(('BLINK','FarmCPU')) else r'\s+')
   if method.startswith('REGENIE'):
    x=raw.LOG10P.to_numpy();sig=x > -np.log10(5e-8);definition='source_LOG10P'
   else:
    p=raw['P.value' if method.startswith(('BLINK','FarmCPU')) else 'P_BOLT_LMM_INF' if method=='BOLT-LMM' else 'P'].to_numpy()
    x=-np.log10(np.maximum(p,1e-300));sig=p<5e-8;definition='minus_log10_source_P'
   rows.append({'method':method,'phenotype':phenotype,'n_snps':len(raw),'n_significant_snps':int(sig.sum()),'max_minus_log10p':float(np.max(x)),'statistic_source':definition,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 return pd.DataFrame(rows)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True);p.add_argument('--out',required=True);a=p.parse_args()
 out=Path(a.out)
 if out.exists():raise SystemExit('Use a new output path.')
 out.parent.mkdir(parents=True,exist_ok=True)
 calculate(a.root,Path(__file__).resolve().parents[2]).to_csv(out,sep='\t',index=False,float_format='%.15g')
