#!/usr/bin/env python3
"""Refresh canonical spatial table with the verified historical closed endpoints."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from derive_association_tables import METHODS,config,assign,save,P
from significance import significant

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--benchmark',type=Path,required=True);ap.add_argument('--gcta-directory',type=Path,required=True);ap.add_argument('--tables',type=Path,required=True);a=ap.parse_args()
    ng=pd.read_csv(a.tables/'null_inventory_50kb.tsv',sep='\t');rows=[];audit=[]
    for m in METHODS:
        for ph in range(1,11):
            path,sep,cols,islog=config(a.benchmark,m,ph,a.gcta_directory)
            if not path.exists():continue
            d=pd.read_csv(path,sep=sep,usecols=cols[:4]).rename(columns=dict(zip(cols[:4],['chr','pos','snp','value'])))
            v=d.value.values;valid=np.isfinite(v);sig=significant(v,is_log10=islog)
            bid=assign(d.chr.values,d.pos.values,ng);usable=(bid>=0)&valid
            tested=np.bincount(bid[usable],minlength=len(ng));hits=np.bincount(bid[usable&sig],minlength=len(ng))
            boundary=usable&(bid>0);ii=np.flatnonzero(boundary);k=bid[ii]
            boundary[ii]=(d.pos.values[ii]==ng.start.values[k])&(d.pos.values[ii]==ng.end.values[k-1])&(d.chr.values[ii]==ng.chr.values[k-1])
            tested+=np.bincount(bid[boundary]-1,minlength=len(ng));hits+=np.bincount(bid[boundary&sig]-1,minlength=len(ng))
            ratio=np.divide(hits,tested*P,out=np.full(len(ng),np.nan),where=tested>0)
            z=ng[['distance_index']].copy();z['tested']=tested;z['hits']=hits;z['ratio']=ratio
            for ix,g in z.groupby('distance_index'):
                rows.append(dict(configuration=m,phenotype=ph,distance_index=ix,observed_expected_equal_nonempty_bins=g.ratio.mean(),nonempty_bins=int((g.tested>0).sum()),observed_expected_pooled_snps=g.hits.sum()/(g.tested.sum()*P)))
            audit.append(dict(configuration=m,phenotype=ph,exact_internal_boundary_markers=int(boundary.sum()),significant_boundary_markers=int((boundary&sig).sum()),endpoint_rule='closed, duplicating shared endpoints'))
    save(pd.DataFrame(rows),a.tables,'spatial_by_phenotype.tsv');save(pd.DataFrame(audit),a.tables,'spatial_boundary_audit.tsv')

if __name__=='__main__':main()
