#!/usr/bin/env python3
"""Read-only inclusive-threshold reconstruction of the preserved fine-map loci."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

p=argparse.ArgumentParser()
p.add_argument('--benchmark',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
B=a.benchmark.resolve()
manifest_path=B/'SuSiE_adj_finemapping/loci/manifest.tsv'
man=pd.read_csv(manifest_path,sep='\t')
bim=pd.read_csv(B/'SuSiE_adj_finemapping/genotype/spg.20k.bim',sep=r'\s+',header=None,
                names=['chr','SNP','cm','position','A1','A2'])
threshold=5e-8
log_threshold=-np.log10(threshold)
rows=[]
for ph in range(1,11):
    for method in ['BOLT-LMM','REGENIE-LOCO']:
        if method=='BOLT-LMM':
            d=pd.read_csv(B/f'BOLTLMM/bolt_pheno{ph}.stats',sep='\t',
                usecols=['SNP','CHR','BP','P_BOLT_LMM_INF','BETA','SE']).rename(
                columns={'CHR':'chr','BP':'position','P_BOLT_LMM_INF':'P'})
            d['logP']=-np.log10(d.P)
            at_threshold=d.P.eq(threshold)
            inclusive=d.P.le(threshold)
        else:
            d=pd.read_csv(B/f'REGENIE/regenie_pheno{ph}_pheno{ph}.regenie',sep=r'\s+',
                usecols=['CHROM','GENPOS','LOG10P','BETA','SE','N','TEST'])
            d=d[d.TEST=='ADD'].rename(columns={'CHROM':'chr','GENPOS':'position','LOG10P':'logP'})
            assert d.N.eq(20000).all()
            d=d.merge(bim[['chr','position','SNP']],on=['chr','position'],validate='one_to_one')
            at_threshold=d.logP.eq(log_threshold)
            inclusive=d.logP.ge(log_threshold)
        assert len(d)==135358 and d.SNP.is_unique
        assert np.isfinite(d[['BETA','SE','logP']]).all().all() and d.SE.gt(0).all()
        assert np.array_equal(inclusive,d.logP.ge(log_threshold))
        selected=[]
        for ch,q in d[inclusive].groupby('chr'):
            q=q.sort_values(['logP','position','SNP'],ascending=[False,True,True])
            remaining=q.index.to_numpy()
            positions=d.position.to_numpy()
            leads=[]
            while len(remaining):
                ix=remaining[0]
                leads.append(ix)
                remaining=remaining[np.abs(positions[remaining]-positions[ix])>1000000]
            z=d.loc[leads].sort_values('position')
            starts=np.maximum(1,z.position.to_numpy()-1000000)
            ends=np.minimum(bim.loc[bim.chr==ch,'position'].max(),z.position.to_numpy()+1000000)
            for j in range(len(z)-1):
                if ends[j]>=starts[j+1]:
                    split=(z.position.iloc[j]+z.position.iloc[j+1])//2
                    ends[j]=split
                    starts[j+1]=split+1
            for k,r in enumerate(z.itertuples()):
                selected.append(dict(chr=ch,start=int(starts[k]),end=int(ends[k]),
                    lead_SNP=r.SNP,lead_position=r.position,lead_logP=r.logP))
        reconstructed=pd.DataFrame(selected).sort_values(['chr','start']).reset_index(drop=True)
        saved=man[(man.method==method)&man.phenotype.eq(f'pheno{ph}')].sort_values(['chr','start']).reset_index(drop=True)
        assert len(reconstructed)==len(saved)
        for col in ['chr','start','end','lead_SNP','lead_position']:
            assert np.array_equal(reconstructed[col],saved[col]),(method,ph,col)
        assert np.max(np.abs(reconstructed.lead_logP-saved.lead_logP))<1e-10
        covered=np.zeros(len(d),dtype=int)
        for ch,q in saved.groupby('chr'):
            ix=np.flatnonzero(d.chr.eq(ch).to_numpy())
            k=np.searchsorted(q.start.to_numpy(),d.position.to_numpy()[ix],side='right')-1
            good=k>=0
            covered[ix[good]]=d.position.to_numpy()[ix[good]]<=q.end.to_numpy()[k[good]]
        assert covered[inclusive].sum()==inclusive.sum()
        assert covered[at_threshold].sum()==at_threshold.sum()
        rows.append(dict(method=method,phenotype=f'pheno{ph}',n_inclusive_significant_snps=int(inclusive.sum()),
            n_at_threshold=int(at_threshold.sum()),n_at_threshold_covered=int(covered[at_threshold].sum()),
            n_loci=len(saved),exact_locus_reconstruction=True,all_significant_markers_covered=True))
result=dict(status='PASS',significance='P <= 5e-8',equivalent_log_rule='-log10(P) >= -log10(5e-8)',
    n_manifest_loci=len(man),n_regions=man.locus_id.nunique(),
    manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
    total_at_threshold=sum(r['n_at_threshold'] for r in rows),
    all_threshold_boundary_markers_in_preserved_loci=True,
    original_loci_preserved=True,new_fits_required=False,
    note='Fine-mapping loci use greedy leads separated by >1 Mb, ±1 Mb flanks, clipped to the last panel marker and split at adjacent-lead midpoints. They are distinct from 500 kb association-count windows.',
    per_analysis=rows)
a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='per_analysis'},indent=2))
