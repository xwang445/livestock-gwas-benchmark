#!/usr/bin/env python3
"""Causal-position detection and null-region SNP association rates, same manifest."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from derive_association_tables import METHODS,config,assign,save

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--benchmark',type=Path,required=True);ap.add_argument('--gcta-directory',type=Path,required=True);ap.add_argument('--tables',type=Path,required=True);a=ap.parse_args()
    grid=pd.read_csv(a.tables/'window_inventory_500kb.tsv',sep='\t',keep_default_na=False)
    causal=pd.read_csv(a.benchmark/'snp.info.chr1_18.10k.csv').SNP.str.split('_',expand=True)
    truth=set(zip(causal[0].astype(int),causal[1].astype(int)))
    alphas=np.unique(np.r_[10.**np.arange(-10,-1.999,.25),5e-8]);rows=[]
    for m in METHODS+['SLEMM_previous']:
        for ph in range(1,11):
            path,sep,cols,islog=config(a.benchmark,m,ph,a.gcta_directory)
            if not path.is_file():continue
            d=pd.read_csv(path,sep=sep,usecols=cols[:4]).rename(columns=dict(zip(cols[:4],['chr','pos','snp','value'])))
            v=d.value.values;valid=np.isfinite(v);wid=assign(d.chr.values,d.pos.values,grid);assert (wid>=0).all()
            masks={'causal_positions':np.array([(c,p) in truth for c,p in zip(d.chr,d.pos)]),
                'central_null':grid.central_null.values[wid], 'all_null':grid.region_type.values[wid]=='null'}
            assert int(masks['causal_positions'].sum())==10000
            for alpha in alphas:
                sig=valid&(v>=-np.log10(alpha) if islog else v<=alpha)
                for region,mask in masks.items():
                    den=int((mask&valid).sum());num=int((mask&sig).sum())
                    rows.append(dict(configuration=m,phenotype=ph,alpha=alpha,region=region,n_tested=den,n_significant=num,rate=num/den if den else np.nan))
    by=pd.DataFrame(rows);save(by,a.tables,'power_and_null_snp_rate_by_phenotype.tsv')
    summary=by.groupby(['configuration','alpha','region'],sort=False).agg(mean_rate=('rate','mean'),sd_rate=('rate','std'),n_phenotypes=('phenotype','nunique'),n_tested_pooled=('n_tested','sum'),n_significant_pooled=('n_significant','sum')).reset_index()
    summary.loc[summary.n_phenotypes!=10,['mean_rate','sd_rate']]=np.nan
    summary['pooled_rate']=summary.n_significant_pooled/summary.n_tested_pooled
    summary.loc[summary.n_phenotypes!=10,'pooled_rate']=np.nan
    save(summary,a.tables,'power_and_null_snp_rate_summary.tsv')

if __name__=='__main__':main()
