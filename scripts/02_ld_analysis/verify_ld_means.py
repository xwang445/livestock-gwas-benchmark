#!/usr/bin/env python3
"""Stream restricted PLINK LD, emit 10-kb aggregate means only; no input writes."""
import argparse,json,hashlib
import numpy as np
import pandas as pd
def verify(path):
    counts=np.zeros(1000,dtype=np.int64);sums=np.zeros(1000,dtype=float)
    for d in pd.read_csv(path,sep=r'\s+',usecols=['CHR_A','BP_A','CHR_B','BP_B','R2'],chunksize=500000):
        distance=np.abs(d.BP_B.to_numpy()-d.BP_A.to_numpy());r2=d.R2.to_numpy()
        keep=(d.CHR_A.to_numpy()==d.CHR_B.to_numpy())&np.isfinite(r2)&(distance<=10000000)
        if np.any((r2[keep]<0)|(r2[keep]>1)):raise ValueError('LD outside [0,1]')
        b=np.minimum(distance[keep]//10000,999).astype(int)
        counts+=np.bincount(b,minlength=1000);sums+=np.bincount(b,weights=r2[keep],minlength=1000)
    return [dict(bin_id=i+1,bin_start=i/100,bin_end=(i+1)/100,mid_mb=(i+.5)/100,n_pairs=int(counts[i]),mean_r2=float(sums[i]/counts[i])) for i in range(1000) if counts[i]]
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);a=p.parse_args()
    print(json.dumps(dict(status='LD_BIN_MEANS_RECOMPUTED',bins=verify(a.input)),allow_nan=False))
