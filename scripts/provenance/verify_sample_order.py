#!/usr/bin/env python3
"""Read-only sample/phenotype chain and retained MPH GRM verification.

Prints counts, hashes and error summaries only. Never exports IDs or genotypes.
The genotype check compares all GRM diagonals and 4096 fixed off-diagonal pairs;
it is not a bytewise reconstruction of the complete GRM.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

def fingerprint(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def id_hash(x):return hashlib.sha256(('\n'.join(x)+'\n').encode()).hexdigest()
def family_ids(p):return pd.read_csv(p,sep=r'\s+',header=None,dtype=str)[1].tolist()
def verify(simulation,candidate,gwas,benchmark,replaced=None,grm_check=False):
 d,b=Path(simulation),Path(benchmark)
 sources={}
 def remember(label,p):sources[label]={'sha256':fingerprint(p),'bytes':Path(p).stat().st_size}
 def compare(a,b):return {'rows_a':len(a),'rows_b':len(b),'same_order':a==b,'same_set':set(a)==set(b),'order_sha256_a':id_hash(a),'order_sha256_b':id_hash(b)}
 cols=['iid']+[str(i) for i in range(1,11)]+['gv'+str(i) for i in range(1,11)]
 sim=pd.read_csv(d/'pheno.sim.csv',usecols=cols,dtype={'iid':str})
 sim_ids=sim.iid.tolist();n=len(sim)
 iid_candidates=[d/'MPH_grm_chr1_18.grm.iid',d/'MPH_grm_chr1_18.iid']
 iid_path=next(p for p in iid_candidates if p.is_file())
 grm_ids=pd.read_csv(iid_path,sep=r'\s+',header=None,dtype=str).iloc[:,-1].tolist()
 if len(grm_ids)==n+1:grm_ids=grm_ids[1:]
 candidate_ids=family_ids(str(candidate)+'.fam');gwas_ids=family_ids(str(gwas)+'.fam')
 orders={'simulation_vs_grm':compare(sim_ids,grm_ids),'simulation_vs_candidate':compare(sim_ids,candidate_ids),'simulation_vs_gwas':compare(sim_ids,gwas_ids)}
 remember('saved_grm_ids',iid_path)
 for label,prefix in [('retained_chip',candidate),('gwas_panel',gwas)]:
  for ext in ['.fam','.bim']:remember(label+ext,str(prefix)+ext)
 if replaced:orders['simulation_vs_replaced_path']=compare(sim_ids,family_ids(str(replaced)+'.fam'))
 fm=b/'SuSiE_adj_finemapping/genotype/spg.20k.fam'
 if fm.is_file():orders['first20000_vs_finemap']=compare(sim_ids[:20000],family_ids(fm));remember('finemap_fam',fm)
 combined=pd.read_csv(b/'simu_pheno/pheno_combined.txt',sep=r'\s+',dtype={'IID':str})
 orders['first20000_vs_combined']=compare(sim_ids[:20000],combined.IID.tolist());remember('combined_phenotypes',b/'simu_pheno/pheno_combined.txt')
 chains=[]
 for i in range(1,11):
  path=b/f'simu_pheno/fastGWA_data/pheno.{i}.txt';p=pd.read_csv(path,sep=r'\s+',dtype={'IID':str});original=pd.read_csv(d/f'data_20000/pheno.{i}.csv',dtype={'ID':str})
  target=sim[str(i)].to_numpy()[:20000]
  chains.append({'phenotype':i,'ids_equal':p.IID.tolist()==sim_ids[:20000]==original.ID.tolist(),'original_max_abs_error':float(np.max(np.abs(original.iloc[:,-1].to_numpy()-target))),'fastgwa_max_abs_error':float(np.max(np.abs(p.iloc[:,-1].to_numpy()-target))),'combined_max_abs_error':float(np.max(np.abs(combined[f'pheno{i}'].to_numpy()-target)))})
  remember('fastgwa_pheno'+str(i),path)
 variances=[]
 for i in range(1,11):
  y=sim[str(i)].to_numpy()[:20000];g=sim['gv'+str(i)].to_numpy()[:20000]
  variances.append({'phenotype':i,'n':len(y),'var_y':float(np.var(y,ddof=1)),'var_g':float(np.var(g,ddof=1)),'var_e':float(np.var(y-g,ddof=1)),'var_g_over_var_y':float(np.var(g,ddof=1)/np.var(y,ddof=1))})
 result={'orders':orders,'phenotype_chain':chains,'realized_variances':variances,'sources':sources}
 if grm_check:
  assert candidate_ids==sim_ids
  causal_path=d/'snp.info.chr1_18.10k.csv';causal=set(pd.read_csv(causal_path).SNP.astype(str));bim=pd.read_csv(str(candidate)+'.bim',sep=r'\s+',header=None,dtype=str)
  indices=np.flatnonzero(bim[1].isin(causal));assert len(indices)==10000
  remember('causal_list',causal_path)
  diag=np.zeros(n);rng=np.random.default_rng(20261104);a=rng.integers(0,n,4096);c=rng.integers(0,n,4096);c=np.where(c==a,(c+1)%n,c);pair=np.zeros(len(a));missing=0
  width=(n+3)//4;lookup=np.array([2.,np.nan,1.,0.])
  with open(str(candidate)+'.bed','rb') as f:
   assert f.read(3)==b'\x6c\x1b\x01'
   for idx in indices:
    f.seek(3+int(idx)*width);raw=np.frombuffer(f.read(width),dtype=np.uint8);x=lookup[((raw[:,None]>>np.array([0,2,4,6]))&3).ravel()[:n]]
    missing+=int(np.isnan(x).sum());mean=np.nanmean(x);scale=np.sqrt(mean*(1-mean/2));assert scale>0
    z=np.nan_to_num((x-mean)/scale);diag+=z*z;pair+=z[a]*z[c]
  diag/=len(indices);pair/=len(indices)
  path=d/'MPH_grm_chr1_18.grm.bin'
  with path.open('rb') as f:
   saved_n=int(np.frombuffer(f.read(4),dtype='<i4')[0]);weight=float(np.frombuffer(f.read(4),dtype='<f4')[0])
  assert saved_n==n and weight==10000
  saved=np.memmap(path,dtype='<f4',mode='r',offset=8)
  col=np.arange(n,dtype=np.int64);diag_index=col*n-col*(col-1)//2
  low=np.minimum(a,c).astype(np.int64);high=np.maximum(a,c);pair_index=low*n-low*(low-1)//2+high-low
  errors=np.abs(diag-saved[diag_index]/weight);pair_errors=np.abs(pair-saved[pair_index]/weight)
  result['grm_check']={'n':n,'causal_snps':len(indices),'missing_calls':missing,'diagonals_checked':n,'off_diagonal_pairs_checked':len(pair),'pair_selection_seed':20261104,'max_diagonal_abs_error':float(errors.max()),'max_off_diagonal_abs_error':float(pair_errors.max()),'tolerance':0.0001,'passed':bool(errors.max()<0.0001 and pair_errors.max()<0.0001),'scope':'All row diagonals and deterministic off-diagonal sample; not full matrix equality.'}
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for name in ['simulation-root','candidate-prefix','gwas-prefix','benchmark-root']:p.add_argument('--'+name,required=True)
 p.add_argument('--replaced-prefix');p.add_argument('--grm-check',action='store_true');a=p.parse_args()
 print(json.dumps(verify(a.simulation_root,a.candidate_prefix,a.gwas_prefix,a.benchmark_root,a.replaced_prefix,a.grm_check),indent=2))
