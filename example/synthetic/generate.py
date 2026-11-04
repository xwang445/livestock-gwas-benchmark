#!/usr/bin/env python3
"""SYNTHETIC EXAMPLE FOR CODE VALIDATION ONLY.

32 artificial animals, 36 SNPs, two chromosomes, two phenotypes. The nine
format fixtures encode the same artificial statistics; they are NOT fits of
the nine GWAS methods and must never be used as scientific benchmark results.
"""
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import chi2
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from benchmark_core import METHODS

def generate(out):
    out=Path(out)
    if out.exists() and any(out.iterdir()):raise ValueError('Synthetic output must be empty or new')
    out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20261104)
    positions=[100,499999,500000,1999999,2000000,2000001,2999999,3000000,3000001,3499999,3500000,4999999,5000000,5000001,5999999,6000000,6000001,7000000]
    ch=np.repeat([1,2],len(positions));bp=np.tile(positions,2);ids=[f'{c}_{p}' for c,p in zip(ch,bp)]
    dosage=rng.binomial(2,rng.uniform(.1,.45,len(ids)),size=(32,len(ids)))
    animals=[f'SYNTHETIC_ANIMAL_{i:03}' for i in range(1,33)]
    pd.DataFrame(dosage,index=animals,columns=ids).rename_axis('IID').to_csv(out/'synthetic_dosages.tsv',sep='\t')
    phenotype=rng.normal(size=(32,2))+.25*dosage[:,[0]]
    pd.DataFrame(phenotype,index=animals,columns=['pheno1','pheno2']).rename_axis('IID').to_csv(out/'synthetic_phenotypes.tsv',sep='\t')
    pd.DataFrame([(c,s,s+2000000) for c in [1,2] for s in [0,6000000]],columns=['chr','start','end']).to_csv(out/'eligible_regions.csv',index=False)
    pd.DataFrame({'SNP':[ids[0],ids[18]]}).to_csv(out/'snp.info.chr1_18.10k.csv',index=False)
    for ph in [1,2]:
        p=np.tile([2e-8,5e-8,1e-6,.5,.02,.8],6);p=np.roll(p,ph-1);x=chi2.isf(p,1);beta=np.sqrt(x)*.1
        for method in METHODS:
            if method=='BOLT-LMM':
                path=out/'BOLTLMM'/f'bolt_pheno{ph}.stats'
                d=dict(SNP=ids,CHR=ch,BP=bp,ALLELE1='A',ALLELE0='G',BETA=beta,SE=.1,P_BOLT_LMM_INF=p)
            elif method.startswith('REGENIE'):
                loco=method=='REGENIE_LOCO';path=out/('REGENIE' if loco else 'REGENIE_WO_LOCO')/(f'regenie_pheno{ph}_pheno{ph}.regenie' if loco else f'regenie_pheno_wo_loco{ph}_pheno{ph}.regenie')
                d=dict(ID=ids,CHROM=ch,GENPOS=bp,ALLELE1='A',ALLELE0='G',BETA=beta,SE=.1,CHISQ=x,LOG10P=-np.log10(p))
            elif method=='fastGWA':
                path=out/'fastGWA'/f'spg.even.100k.pheno{ph}.fastgwa.fastGWA'
                d=dict(SNP=ids,CHR=ch,POS=bp,A1='A',A2='G',BETA=beta,SE=.1,P=p)
            elif method=='SLEMM':
                path=out/'SLEMM'/f'spg.even.100k.pheno{ph}.allchr.gwas.txt'
                d={'ID':ids,'#CHROM':ch,'POS':bp,'REF':'G','ALT':'A','BETA':beta,'SE':.1,'CHISQ':x,'P':p}
            else:
                m,city=method.split('_');path=out/m/f'GAPIT.Association.GWAS_Results.{m}.pheno{ph}({city}).csv'
                d=dict(SNP=ids,Chr=ch,Pos=bp,**{'P.value':p},effect=beta)
            path.parent.mkdir(parents=True,exist_ok=True)
            pd.DataFrame(d).to_csv(path,sep=',' if method.startswith(('BLINK_','FarmCPU_')) else '\t',index=False)
    (out/'README.txt').write_text(__doc__+'\n')
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True);a=p.parse_args();generate(a.out)
