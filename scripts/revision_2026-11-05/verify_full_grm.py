#!/usr/bin/env python3
"""Recompute selected off-diagonal full-GRM entries from all panel markers.

Only aggregate discrepancies are written by default. Individual relationship
entries and sample indices require an explicitly private output path.
"""
import argparse
import json
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--bfile',type=Path,required=True)
p.add_argument('--grm',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--private-details',type=Path)
a=p.parse_args()
n=sum(1 for _ in open(str(a.bfile)+'.fam'))
m=sum(1 for _ in open(str(a.bfile)+'.bim'))
nb=(n+3)//4
with open(str(a.bfile)+'.bed','rb') as f: assert f.read(3)==b'\x6c\x1b\x01'
raw=np.memmap(str(a.bfile)+'.bed',mode='r',dtype=np.uint8,offset=3).reshape(m,nb)
sel=np.unique([x for x in [0,99,999,9999,19999,25000,n-1] if x<n])
lut=np.array([2.,np.nan,1.,0.])
codes=(np.arange(256)[:,None]>>np.array([0,2,4,6]))&3
sums=np.nan_to_num(lut[codes]).sum(1)
missing=(codes==1).sum(1)
A=np.zeros((len(sel),len(sel)))
N=np.zeros_like(A)
total_missing=validmarkers=0
for st in range(0,m,1024):
    r=np.asarray(raw[st:st+1024])
    if n%4:
        tailcodes=(r[:,-1,None]>>np.arange(0,2*(n%4),2))&3
        ac=sums[r[:,:-1]].sum(1)+np.nan_to_num(lut[tailcodes]).sum(1)
        miss=missing[r[:,:-1]].sum(1)+(tailcodes==1).sum(1)
    else:
        ac=sums[r].sum(1)
        miss=missing[r].sum(1)
    total_missing+=int(miss.sum())
    freq=ac/(2*(n-miss))
    var=2*freq*(1-freq)
    ok=var>0
    x=lut[(r[:,sel//4]>>(2*(sel%4)))&3]
    observed=np.isfinite(x)&ok[:,None]
    z=np.zeros_like(x)
    np.divide(x-2*freq[:,None],np.sqrt(var[:,None]),out=z,where=observed)
    A+=z.T@z
    N+=observed.T.astype(float)@observed.astype(float)
    validmarkers+=int(ok.sum())
G=A/N
grm=np.memmap(str(a.grm)+'.grm.bin',dtype=np.float32,mode='r')
gn=np.memmap(str(a.grm)+'.grm.N.bin',dtype=np.float32,mode='r')
rows=[]
for i,ix in enumerate(sel):
    for j,jx in enumerate(sel[:i]):
        k=ix*(ix+1)//2+jx
        rows.append(dict(sample_index_i=int(ix),sample_index_j=int(jx),expected=float(G[i,j]),
            stored=float(grm[k]),abs_difference=abs(float(G[i,j])-float(grm[k])),
            n_markers_expected=int(N[i,j]),n_markers_stored=float(gn[k])))
mx=max(r['abs_difference'] for r in rows)
counts_match=all(r['n_markers_expected']==r['n_markers_stored'] for r in rows)
result=dict(status='PASS' if mx<1e-5 and counts_match else 'FAIL',n_full_samples=n,n_markers=m,
    n_valid_markers=validmarkers,full_panel_missing_genotype_calls=total_missing,
    n_offdiagonal_entries_checked=len(rows),max_abs_difference=mx,pair_marker_counts_match=counts_match,
    formula='Mean over jointly observed markers of (g_i-2p)*(g_j-2p)/(2p*(1-p)); p from all panel animals',
    limitation='Spot-check of selected off-diagonal entries, not a complete recomputation of every GRM entry')
a.output.write_text(json.dumps(result,indent=2)+'\n')
if a.private_details:a.private_details.write_text(json.dumps(dict(result,rows=rows),indent=2)+'\n')
assert result['status']=='PASS'
