#!/usr/bin/env python3
"""Compare regenerated aggregates with October manuscript and historical references.

PASS is numeric agreement at a stated tolerance, not evidence of a new GWAS fit.
Reference provenance is in audit/reference_values/source_manifest.tsv.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from benchmark_core import output_table

def validate(repo):
    repo=Path(repo);D=repo/'derived_data';A=repo/'audit/reference_values';rows=[]
    def read(p):
        d=pd.read_csv(p,sep='\t');d.attrs['source']=str(p.relative_to(repo));return d
    def claim(item,expected,observed,tolerance,source,note='',status=None):
        difference=abs(float(observed)-float(expected)) if isinstance(observed,(float,int,np.number)) else np.nan
        if status is None:status='PASS' if np.isfinite(difference) and difference<=tolerance else 'FAIL'
        rows.append(dict(manuscript_item=item,expected=expected,observed=observed,absolute_difference=difference,tolerance=tolerance,status=status,source=source,note=note or 'Not applicable'))
    def compare(label,actual,reference,keys,values,tolerance=1e-10):
        if actual.duplicated(keys).any() or reference.duplicated(keys).any():raise ValueError('Nonunique validation keys: '+label)
        z=reference.merge(actual,on=keys,how='left',suffixes=('_ref','_new'),validate='one_to_one',indicator=True)
        source=actual.attrs.get('source',label)+'; '+reference.attrs.get('source',label)
        claim(label+' matched rows',len(reference),int(z._merge.eq('both').sum()),0,source)
        for new,old in values.items():
            x=z[new+'_new' if new in reference else new];y=z[old+'_ref' if old in actual else old]
            delta=np.abs(x-y);bad=(x.isna()!=y.isna())
            maxdiff=float(delta.max()) if delta.notna().any() else 0.
            claim(label+' '+new,0,np.inf if bad.any() else maxdiff,tolerance,source,'Maximum absolute difference across matched rows; equal NA pairs allowed.')
    c=read(D/'method_summary/significant_snps_by_pheno.tsv')
    compare('Figure3',c,read(A/'figure3.tsv'),['method','phenotype'],{'n_significant_snps':'original_snps','n_significant_loci':'original_loci'},0)
    claim('GWAS configurations',9,c.method.nunique(),0,'method_summary/significant_snps_by_pheno.tsv')
    claim('Phenotypes',10,c.phenotype.nunique(),0,'method_summary/significant_snps_by_pheno.tsv')
    claim('All GWAS panels have 135358 SNPs',90,int(c.n_snps.eq(135358).sum()),0,'method_summary/significant_snps_by_pheno.tsv')
    man=read(repo/'figures/reconstructed/manhattan/manhattan_manifest.tsv')
    compare('Figure4 and supplements counts',man,c,['method','phenotype'],{'n_significant':'n_significant_snps','n_plotted':'n_snps'},0)
    claim('Figure4 maximum minus log10 P approximately 102',102,float(man.max_minus_log10p.max()),1,'figures/reconstructed/manhattan/manhattan_manifest.tsv','Untruncated statistic across all 90 panels.')
    s=read(D/'spillover_summary/spillover_by_pheno.tsv')
    compare('Figure6A historical spillover',s,read(A/'spillover.tsv'),['method','phenotype'],{'spillover_rate':'spillover_locus_rate','n_spillover_loci':'spillover_locus_count'})
    manuscript_spill={'BOLT-LMM':30.7,'REGENIE_LOCO':30.3,'REGENIE_woLOCO':1.1,'fastGWA':1.5,'FarmCPU_NYC':2.9,'FarmCPU_Kansas':1.9,'BLINK_NYC':1.1,'BLINK_Kansas':1.0,'SLEMM':0.0}
    for method,value in manuscript_spill.items():claim('Figure6A percent '+method,value,100*s[s.method==method].spillover_rate.mean(),.05000001,'spillover_summary/spillover_by_pheno.tsv','Manuscript rounded to one decimal percentage point.')
    lam=read(D/'lambda_gc_summary/lambda_gc_by_pheno.tsv')
    compare('Figure6B historical lambda',lam,read(A/'lambda.tsv'),['method','phenotype'],{'lambda_gc':'lambda_gc'},1e-10)
    for name in ['power','fpr']:
        a=read(D/f'power_summary/{name}_by_pheno.tsv');r=read(A/f'{name}.tsv')
        if name=='fpr':a=a[a.region=='central_null']
        # The same 33 nominal thresholds, printed at different floating precision.
        a=a.copy();r=r.copy();a['threshold']=np.round(np.log10(a.alpha),10);r['threshold']=np.round(np.log10(r.alpha),10)
        compare('Historical '+name,a,r,['method','phenotype','threshold'],{name:name},1e-10)
    compare('Figure6CD 50kb enrichment',read(D/'power_summary/spatial_enrichment_by_pheno.tsv'),read(A/'spatial.tsv'),['method','phenotype','distance_index'],{'observed_expected':'observed_expected'},1e-7)
    ov=read(D/'overlap_summary/locus_overlap_by_pheno.tsv');ref=read(A/'figure5_spot_checks.tsv')
    compare('Figure5 displayed spot checks',ov,ref,['method_a','method_b','phenotype'],{'jaccard':'original'},.00500001)
    eq=read(D/'equation13_summary/equation13_by_pheno.tsv');eqs=read(D/'equation13_summary/equation13_summary.tsv')
    for region,part in [('Genome-wide','full'),('Central-null','central')]:
        compare('Figure7 '+region+' replicate bins',eq[eq.region==region],read(A/f'eq13_{part}.tsv'),['method','phenotype','bin'],{'mean_observed_chisq':'mean_observed_chisq','n_obs':'n_obs'},1e-9)
        compare('Figure7 '+region+' means and intervals',eqs[eqs.region==region],read(A/f'eq13_{part}_summary.tsv'),['method','bin'],{'mean_observed_chisq':'mean_observed_chisq','predicted_chisq':'predicted_chisq','ci_lower':'ci_lower','ci_upper':'ci_upper'},1e-9)
    ld=read(D/'ld_summary/ld_decay_10kb.tsv')
    compare('Figure1 LD bins',ld,read(A/'ld.tsv'),['bin_id'],{'n_pairs':'n_pairs','mean_r2':'mean_r2'},1e-10)
    for index,expected in [(1,.59),(200,.16),(1000,.04)]:
        claim('Figure1 LD bin '+str(index),expected,ld.loc[ld.bin_id==index,'mean_r2'].iloc[0],.01,'ld_summary/ld_decay_10kb.tsv','Approximate manuscript value; tolerance 0.01 in R squared.')
    fm=read(D/'finemapping_summary/high_pip_summary.tsv')
    targets=[('BOLT-LMM','SuSiE-adj',1974,1153,263,58.41,13.3,8.8,1.51),('BOLT-LMM','SuSiE-RSS',2009,1170,261,58.24,13.0,8.8,1.47),('REGENIE_LOCO','SuSiE-adj',2101,1274,256,60.64,12.2,8.7,1.40),('REGENIE_LOCO','SuSiE-RSS',2064,1247,253,60.42,12.3,8.7,1.40)]
    for method,pipe,*expected in targets:
        observed=fm[(fm.method==method)&(fm.pipeline==pipe)].iloc[0]
        for column,value,tol in zip(['n_high_pip','n_null','n_causal','null_percent','causal_percent','maf_expected_percent','enrichment'],expected,[0,0,0,.00500001,.05000001,.05000001,.00500001]):
            claim('Figure8 '+method+' '+pipe+' '+column,value,float(observed[column]),tol,'finemapping_summary/high_pip_summary.tsv','Independent recount of saved converged PIPs; MAF uses saved annotations.')
    fits=read(D/'finemapping_summary/converged_loci.tsv')
    for r in fits.itertuples():claim('Figure8 converged '+r.method+' '+r.pipeline,6732 if r.method=='BOLT-LMM' else 6569,r.n_converged_loci,0,'finemapping_summary/converged_loci.tsv')
    claim('High-PIP calls approximately 87% noncausal',87,100*(1-fm.n_causal.sum()/fm.n_high_pip.sum()),1,'finemapping_summary/high_pip_summary.tsv','Pooled across four pipelines; overlaps are not independent animals or unique SNPs.')
    h=read(D/'method_summary/heritability_from_logs.tsv')
    for method,ex in [('fastGWA',[.145,.135,.153]),('SLEMM',[.529,.506,.550])]:
        v=h[h.method==method].h2
        claim(method+' h2 log estimates',10,len(v),0,'method_summary/heritability_from_logs.tsv')
        for stat,expected,observed in zip(['mean','min','max'],ex,[v.mean(),v.min(),v.max()]):claim(method+' h2 '+stat,expected,float(observed),.00050001,'method_summary/heritability_from_logs.tsv','Reaggregated saved fitting logs, not refitted.')
    expected_hits=10000*2*norm.sf(np.sqrt(chi2.isf(5e-8,1)/2))
    claim('Idealized expected causal detections',1.2,float(expected_hits),.05,'scripts/validate_manuscript.py','Analytic normal-mixture calculation only: 10000 × Pr[2 chi-square(1) > critical]. Does not verify recovered maximum PVE.')
    claim('Results REGENIE maximum SNP count',21484,int(c[c.method=='REGENIE_LOCO'].n_significant_snps.max()),0,'method_summary/significant_snps_by_pheno.tsv','Figure3 is correct at 20461; 21484 belongs to BOLT-LMM.')
    for item,status,note in [
        ('Maximum SNP PVE 0.063%-0.093%','AMBIGUOUS','Recovered LSMR estimates use 31301 animals and give 0.209%-0.433%; exact manuscript source not found.'),
        ('Equation12 covariance prediction generator','REQUIRES_RECOMPUTATION','Cached theory algebra and downstream bins verified; original generator unavailable.'),
        ('Independent explicit SNP-effect sampling','AMBIGUOUS','MPH draws independent Gaussian genetic vectors from a GRM Cholesky factor; no original explicit effect draws recovered.'),
        ('Historical 31301 to 20000 sample identity/order','AMBIGUOUS','Counts and order-selection code recovered; historical source genotype labels/order changed in place.'),
        ('GAPIT executed version 3.5','AMBIGUOUS','Retained source identifies 4.0, but historical logs do not prove which source version executed.'),
        ('New GWAS fits and SuSiE fits','NOT_RUN','Validation recounts saved results; expensive fitting was not repeated.'),
        ('10k fixed causal panel regeneration','NOT_RUN','Saved truth has 10000 distinct positions and a shared input hash; genotype-based resampling not repeated.'),
        ('Target h2=0.5 and independent residual draws','NOT_RUN','Recovered command, weight and versioned MPH source support the design; original random sequence was not re-executed.')]:
        claim(item,'See October Methods','See note',np.nan,'docs/statistical_definitions.md',note,status)
    result=pd.DataFrame(rows);output_table(result,repo/'audit/manuscript_validation.tsv')
    print(result.status.value_counts().to_string())
    print(result[result.status=='FAIL'][['manuscript_item','observed','note']].to_string(index=False))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',default=Path(__file__).resolve().parents[1]);a=p.parse_args();validate(a.repo)
