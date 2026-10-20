#!/usr/bin/env python3
"""Ordinary Pearson LD from the same 20k animals. No GRM or phenotype is used."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
B = Path(__file__).resolve().parents[2]
W = B/'SuSiE_finemapping'
sys.path.insert(0, str(B/'SuSiE_adj_finemapping/loco/scripts'))
from common import read_bim, read_fam, read_bed_cols, standardize
MAN = pd.read_csv(B/'SuSiE_adj_finemapping/loci/manifest.tsv', sep='\t')
BIM, FAM = read_bim(), read_fam()
N = len(FAM)
assert N == 20000
assert Path(B/'SuSiE_adj_finemapping/genotype/spg.20k.bed').open('rb').read(3) == b'\x6c\x1b\x01'
assert set(FAM.IID) == set(pd.read_csv(B/'simu_pheno/IID.txt', dtype=str).IID)

def chromosome(chrom):
    t = time.time()
    regions = MAN[MAN.chr == chrom].drop_duplicates('locus_id').sort_values(['start','end'])
    idx = np.flatnonzero(BIM.chr.to_numpy() == chrom)
    b = BIM.iloc[idx].reset_index(drop=True)
    X, ok = standardize(read_bed_cols(idx,N))
    assert ok.all(), 'Unexpected invariant SNP'
    pos = b.pos.to_numpy()
    # Compute local LD blocks once, reusing overlapping regions across phenotypes.
    blocks = {}
    checks = []
    for j, loc in enumerate(regions.itertuples()):
        out = W/'ld'/loc.locus_id
        if (out/'DONE').exists(): continue
        out.mkdir(parents=True, exist_ok=True)
        lo, hi = np.searchsorted(pos, [loc.start, loc.end], side='left')
        hi = np.searchsorted(pos, loc.end, side='right')
        assert hi-lo == loc.requested_snps
        key = lo//256
        if key not in blocks or hi > blocks[key][1]:
            start = key*256
            end = min(len(b),max(start+512,hi))
            z = X[:,start:end]
            r = z.T @ z/(N-1)
            r = (r+r.T)/2
            blocks = {key:(start,end,r)}
        start, end, block = blocks[key]
        R = block[lo-start:hi-start,lo-start:hi-start].copy()
        assert np.isfinite(R).all() and np.max(np.abs(np.diag(R)-1)) < 1e-10
        np.fill_diagonal(R,1.)
        if j == 0:
            ref = np.corrcoef(X[:,lo:hi],rowvar=False)
            err = float(np.max(np.abs(R-ref)))
            mineig = float(np.linalg.eigvalsh(R).min())
            assert err < 1e-10 and mineig > -1e-8
            checks.append(dict(chr=chrom,locus_id=loc.locus_id,n=N,n_snps=hi-lo,max_diff_numpy_corrcoef=err,min_eigenvalue=mineig))
        R.astype(np.float64).tofile(out/'R.bin')
        b.iloc[lo:hi][['SNP','pos','A1','A2']].rename(columns={'pos':'position','A1':'counted','A2':'other'}).assign(ld_index=np.arange(1,hi-lo+1)).to_csv(out/'snps.tsv',sep='\t',index=False)
        pd.DataFrame([dict(n=N,n_eff=N,n_snps=hi-lo,eta=0,eta_source='ordinary_Pearson_LD')]).to_csv(out/'summary.tsv',sep='\t',index=False)
        (out/'DONE').write_text('ordinary LD; N=20000\n')
    if checks: pd.DataFrame(checks).to_csv(W/f'qc/ld_chr{chrom}.tsv',sep='\t',index=False)
    print(f'chr{chrom}: {len(regions)} regions complete ({time.time()-t:.1f}s)',flush=True)
    return len(regions)

if __name__ == '__main__':
    chromosomes = [int(sys.argv[1])] if len(sys.argv)>1 else list(range(1,19))
    with ProcessPoolExecutor(max_workers=min(12,len(chromosomes))) as pool:
        counts=list(pool.map(chromosome,chromosomes))
    print('Completed regions:',sum(counts),flush=True)
    if len(chromosomes)==18:
        assert sum(counts)==MAN.locus_id.nunique()
        (W/'status/LD.COMPLETE').write_text(str(sum(counts))+'\n')
