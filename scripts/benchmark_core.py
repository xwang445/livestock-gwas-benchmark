#!/usr/bin/env python3
"""Evidence-preserving downstream GWAS reconstruction.

Reads existing restricted inputs without altering them. Only aggregate tables are
written by `summarize`; optional harmonized per-variant exports remain restricted.
See docs/statistical_definitions.md for intentionally different historical rules.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2, t

METHODS = ('BOLT-LMM', 'REGENIE_LOCO', 'REGENIE_woLOCO', 'fastGWA',
           'FarmCPU_NYC', 'FarmCPU_Kansas', 'BLINK_NYC', 'BLINK_Kansas', 'SLEMM')
ALPHA = 5e-8
SCHEMA = ['phenotype','method','snp','chr','bp','effect_allele','other_allele','beta','se','p','chisq']

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()

def output_table(df, path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(path,sep='\t',index=False,na_rep='NA',float_format='%.15g')

def protect_source(root,*outputs):
    """Require outputs/caches outside the read-only historical source tree."""
    source=Path(root).resolve()
    for output in outputs:
        if Path(output).resolve().is_relative_to(source):
            raise ValueError('Output/cache must be outside the historical input tree')

def input_path(root, method, phenotype):
    root=Path(root);p=phenotype
    paths={
        'BOLT-LMM': f'BOLTLMM/bolt_pheno{p}.stats',
        'REGENIE_LOCO': f'REGENIE/regenie_pheno{p}_pheno{p}.regenie',
        'REGENIE_woLOCO': f'REGENIE_WO_LOCO/regenie_pheno_wo_loco{p}_pheno{p}.regenie',
        'fastGWA':f'fastGWA/spg.even.100k.pheno{p}.fastgwa.fastGWA',
        'SLEMM':f'SLEMM/spg.even.100k.pheno{p}.allchr.gwas.txt'}
    if method in paths:
        target=root/paths[method]
        if not target.is_file():raise FileNotFoundError(target)
        return target
    base,place=method.split('_')
    name=f'GAPIT.Association.GWAS_Results.{base}.pheno{p}({place}).csv'
    candidates=[root/base/name,root/'GAPIT'/name]
    found=[f for f in candidates if f.is_file()]
    if not found:raise FileNotFoundError(str(candidates))
    if len(found)>1 and len({sha256(f) for f in found})!=1:
        raise ValueError('Conflicting GAPIT copies: '+str(found))
    return found[0]

def validate_harmonized(d,expected_snps=None):
    if not set(SCHEMA).issubset(d):raise ValueError('Missing harmonized columns')
    if not d.method.isin(METHODS).all():raise ValueError('Unknown method label')
    if d[['phenotype','method','snp','chr','bp']].isna().any().any():raise ValueError('Missing identity/coordinate')
    if d.snp.astype(str).str.strip().eq('').any():raise ValueError('Empty SNP ID')
    if d.duplicated(['method','phenotype','snp']).any():raise ValueError('Duplicate SNP ID within run')
    if d.duplicated(['method','phenotype','chr','bp']).any():raise ValueError('Duplicate coordinates within run')
    for col in ['phenotype','chr','bp']:
        if not (np.isfinite(d[col]) & (d[col]==np.floor(d[col]))).all():raise ValueError('Non-integer '+col)
    if not d.chr.between(1,18).all() or not d.bp.gt(0).all():raise ValueError('Invalid genomic coordinate')
    if not d.phenotype.between(1,10).all():raise ValueError('Invalid phenotype label')
    p=d.p.dropna()
    if not (np.isfinite(p)&p.between(0,1)).all():raise ValueError('P outside [0,1]')
    if (d.se.dropna()<=0).any():raise ValueError('Nonpositive supplied SE')
    if (d.chisq.dropna()<0).any():raise ValueError('Negative chi-square')
    if expected_snps is not None and not d.groupby(['method','phenotype']).size().eq(expected_snps).all():
        raise ValueError('Unexpected SNP count')

def harmonize(path,method,phenotype,expected_snps=None):
    if method not in METHODS:raise ValueError(method)
    gapit=method.startswith(('BLINK_','FarmCPU_'))
    raw=pd.read_csv(path,sep=',' if gapit else r'\s+')
    fields=(['SNP','Chr','Pos','P.value'] if gapit else
            ['SNP','CHR','BP','P_BOLT_LMM_INF'] if method=='BOLT-LMM' else
            ['ID','CHROM','GENPOS','LOG10P'] if method.startswith('REGENIE') else
            ['ID','#CHROM','POS','P'] if method=='SLEMM' else ['SNP','CHR','POS','P'])
    if not set(fields).issubset(raw):raise ValueError('Missing source fields: '+str(set(fields)-set(raw)))
    d=pd.DataFrame({'snp':raw[fields[0]].astype(str),
                    'chr':pd.to_numeric(raw[fields[1]].astype(str).str.replace(r'^chr','',regex=True)),
                    'bp':pd.to_numeric(raw[fields[2]]),
                    'p':pd.to_numeric(raw[fields[3]])})
    d['phenotype']=int(phenotype);d['method']=method
    d['effect_allele']=pd.NA;d['other_allele']=pd.NA
    d['beta']=np.nan;d['se']=np.nan
    if method.startswith('REGENIE'):
        if (d.p.dropna()<0).any():raise ValueError('Negative -log10(P)')
        d['source_log10p']=d.p;d['p']=np.power(10.,-d.p)
    # GAPIT's NYC refinement replaces P values without consistently replacing
    # the original effect column. Allele-oriented beta/SE are not fabricated.
    if not gapit:
        d['beta']=pd.to_numeric(raw['BETA']);d['se']=pd.to_numeric(raw['SE'])
        if method=='BOLT-LMM' or method.startswith('REGENIE'):
            d['effect_allele']=raw['ALLELE1'];d['other_allele']=raw['ALLELE0']
        elif method=='fastGWA':
            d['effect_allele']=raw['A1'];d['other_allele']=raw['A2']
        # SLEMM REF/ALT labels are retained separately until beta orientation is
        # independently verified; no cross-method effect meta-analysis is done.
        else:d['source_ref']=raw['REF'];d['source_alt']=raw['ALT']
    if method=='BOLT-LMM' or method=='fastGWA':
        d['chisq']=(d.beta/d.se)**2;definition='beta_se_squared'
    elif method.startswith('REGENIE') or method=='SLEMM':
        d['chisq']=pd.to_numeric(raw['CHISQ']);definition='reported_CHISQ'
    else:d['chisq']=chi2.isf(d.p,1);definition='P_inverse_survival_df1'
    d['chisq_definition']=definition
    validate_harmonized(d,expected_snps)
    return d[SCHEMA+[c for c in d if c not in SCHEMA]]

def locus_keys(d,alpha=ALPHA,inclusive=False):
    sig=d.p.le(alpha) if inclusive else d.p.lt(alpha)
    return set(zip(d.loc[sig,'chr'].astype(int), (d.loc[sig,'bp']//500000).astype(int)))

def jaccard(a,b):
    # October Figure 5 explicitly treats comparisons with either empty set as NA.
    return len(a&b)/len(a|b) if a and b else np.nan

def null_regions(effects):
    rows=[]
    for ch,g in effects.groupby('chr',sort=True):
        g=g.sort_values('start')
        for end,next_start in zip(g.end.iloc[:-1],g.start.iloc[1:]):
            if next_start-end==4000000:rows.append((int(ch),int(end),int(next_start)))
    return pd.DataFrame(rows,columns=['chr','start','end'])

def central_mask_and_keys(d,nulls):
    mask=np.zeros(len(d),dtype=bool);keys=set()
    for r in nulls.itertuples():
        # Historical closed integer coordinates: (null_start+1Mb, null_end-1Mb].
        start=r.start+1000001;end=r.end-1000000
        hit=(d.chr.to_numpy()==r.chr)&(d.bp.to_numpy()>=start)&(d.bp.to_numpy()<=end)
        mask|=hit
        sig=hit&d.p.le(ALPHA).to_numpy()
        keys.update((r.chr,r.start,int(k)) for k in ((d.loc[sig,'bp']-start)//500000))
    return mask,keys

def lambda_gc(p):
    p=np.asarray(p,dtype=float);p=p[~np.isnan(p)]
    if not len(p):return np.nan
    if np.any((p<0)|(p>1)|~np.isfinite(p)):raise ValueError('Invalid lambda P input')
    return float(np.median(chi2.isf(p,1))/chi2.ppf(.5,1))

def historical_lambda(d):
    """Exact statistic choice and finite-value exclusion in pwer_analysis.R.

    GAPIT used lower-tail qchisq(1-P), so very small P values can round to 1,
    produce infinity and be excluded. This numerical limitation is preserved
    solely for historical reproduction; lambda_gc(P) is reported alongside it.
    """
    if d.method.iloc[0].startswith(('BLINK_','FarmCPU_')):
        p=d.p.to_numpy().copy();p[p<=0]=1e-300
        x=chi2.ppf(1-p,1)
        definition='historical_qchisq_1_minus_p_finite_only'
    else:
        x=(d.beta/d.se).to_numpy()**2
        definition='historical_beta_se_squared_finite_only'
    x=x[np.isfinite(x)]
    return (float(np.median(x)/chi2.ppf(.5,1)) if len(x) else np.nan,len(x),definition)

def spatial_enrichment(d,nulls):
    # Reproduce historical foverlaps closed endpoints, including both bins at a
    # shared exact boundary. Denominator counts rows even if P is missing.
    rows=[]
    for ch,g in d.groupby('chr',sort=True):
        g=g.sort_values('bp');positions=g.bp.to_numpy();p=g.p.to_numpy()
        cumulative=np.r_[0,np.cumsum(p<=ALPHA)]
        for r in nulls[nulls.chr==ch].itertuples():
            starts=r.start+50000*np.arange(80);ends=starts+50000
            lo=np.searchsorted(positions,starts,side='left');hi=np.searchsorted(positions,ends,side='right')
            total=hi-lo;sig=cumulative[hi]-cumulative[lo]
            for k in range(80):
                rows.append((min(k+1,80-k),int(total[k]),int(sig[k]),float(sig[k]/(total[k]*ALPHA)) if total[k] else np.nan))
    bins=pd.DataFrame(rows,columns=['distance_index','n_tested','n_significant','observed_expected'])
    result=bins.groupby('distance_index').agg(observed_expected=('observed_expected','mean'),
        n_nonempty_windows=('observed_expected','count'),n_tested=('n_tested','sum'),n_significant=('n_significant','sum')).reset_index()
    result['pooled_observed_expected']=result.n_significant/(result.n_tested*ALPHA)
    return result

def aggregate(d,groups,values):
    specs={}
    for name,col in values.items():
        specs.update({f'mean_{name}':(col,'mean'),f'sd_{name}':(col,'std'),f'min_{name}':(col,'min'),f'max_{name}':(col,'max')})
    specs['n_phenotypes']=('phenotype','nunique')
    return d.groupby(groups,dropna=False).agg(**specs).reset_index()

def summarize(root,out,cache,expected_snps=135358,phenotypes=range(1,11),methods=METHODS):
    protect_source(root,out,cache)
    root,out,cache=map(Path,[root,out,cache]);cache.mkdir(parents=True,exist_ok=True)
    effects=pd.read_csv(root/'eligible_regions.csv')
    if not {'chr','start','end'}.issubset(effects):raise ValueError('Invalid window schema')
    nulls=null_regions(effects)
    if not len(nulls):raise ValueError('No complete null windows')
    causal=pd.read_csv(root/'snp.info.chr1_18.10k.csv').SNP.str.split('_',expand=True)
    causal_keys=set(zip(causal[0].astype(int),causal[1].astype(int)))
    if expected_snps==135358 and len(causal_keys)!=10000:raise ValueError('Expected 10000 unique causal positions')
    inputs=[];counts=[];spills=[];lambdas=[];spatial=[];powers=[];fprs=[];loci={};reference_coords=None
    alpha_grid=np.r_[10.**np.arange(-10,-1.999,.25),ALPHA];alpha_grid=np.unique(alpha_grid)
    config_hash=hashlib.sha256((''.join(sha256(root/f) for f in ['eligible_regions.csv','snp.info.chr1_18.10k.csv'])+sha256(__file__)).encode()).hexdigest()
    for method in methods:
        for ph in phenotypes:
            path=input_path(root,method,ph);digest=sha256(path)
            cache_file=cache/f'{method}_{ph}.json'
            cached=json.loads(cache_file.read_text()) if cache_file.exists() else None
            if cached and cached.get('input_sha256')==digest and cached.get('config_sha256')==config_hash and cached.get('expected_snps')==expected_snps:
                result=cached;print(f'{method} phenotype {ph}: verified cache',flush=True)
            else:
                d=harmonize(path,method,ph,expected_snps)
                coords=sorted(zip(d.chr.astype(int),d.bp.astype(int)))
                coords_hash=hashlib.sha256(json.dumps(coords,separators=(',',':')).encode()).hexdigest()
                strict=locus_keys(d);inclusive=locus_keys(d,inclusive=True)
                central,spill=central_mask_and_keys(d,nulls)
                hist_lam,n_hist,lam_definition=historical_lambda(d.loc[central])
                hit=np.array([(int(c),int(p)) in causal_keys for c,p in zip(d.chr,d.bp)])
                effect=np.zeros(len(d),bool)
                for r in effects.itertuples():effect|=((d.chr==r.chr)&(d.bp>=r.start)&(d.bp<=r.end)).to_numpy()
                base={'method':method,'phenotype':ph}
                power=[];fpr=[]
                for a in alpha_grid:
                    valid=d.p.notna().to_numpy();sig=d.p.le(a).to_numpy()
                    den=int((hit&valid).sum())
                    power.append(dict(base,alpha=float(a),n_tested_causal=den,n_significant_causal=int((sig&hit&valid).sum()),power=float((sig&hit&valid).sum()/den) if den else np.nan))
                    for region,mask in [('central_null',central),('outside_effect',~effect)]:
                        den=int((mask&valid).sum());num=int((sig&mask&valid).sum())
                        fpr.append(dict(base,region=region,alpha=float(a),n_tested=den,n_significant=num,fpr=num/den if den else np.nan))
                sp=spatial_enrichment(d,nulls);sp['method']=method;sp['phenotype']=ph
                result=dict(input_sha256=digest,config_sha256=config_hash,expected_snps=expected_snps,coordinates_sha256=coords_hash,
                    counts=dict(base,n_snps=len(d),n_significant_snps=int(d.p.lt(ALPHA).sum()),n_significant_loci=len(strict),inclusive_snps=int(d.p.le(ALPHA).sum()),inclusive_loci=len(inclusive),at_threshold=int(d.p.eq(ALPHA).sum()),missing_p=int(d.p.isna().sum()),zero_p=int(d.p.eq(0).sum())),
                    spillover=dict(base,n_significant_loci=len(inclusive),n_spillover_loci=len(spill),n_central_windows=4*len(nulls),spillover_rate=len(spill)/(4*len(nulls))),
                    inflation=dict(base,lambda_gc=hist_lam,lambda_gc_p_based=lambda_gc(d.loc[central,'p']),statistic_definition=lam_definition,n_central_snps=int(central.sum()),n_historical_finite_statistics=n_hist,n_finite_p=int(d.loc[central,'p'].notna().sum())),
                    strict_loci=[list(k) for k in sorted(strict)],spatial=sp.to_dict('records'),power=power,fpr=fpr)
                tmp=cache_file.with_suffix('.tmp');tmp.write_text(json.dumps(result,allow_nan=True));tmp.replace(cache_file)
                print(f'{method} phenotype {ph}: {len(d)} SNPs, {len(strict)} strict loci',flush=True)
            if reference_coords is None:reference_coords=result['coordinates_sha256']
            if result['coordinates_sha256']!=reference_coords:raise ValueError('SNP coordinate panels differ across runs')
            inputs.append(dict(method=method,phenotype=ph,path=str(path.relative_to(root)),sha256=digest,rows=result['counts']['n_snps'],coordinates_sha256=result['coordinates_sha256']))
            counts.append(result['counts']);spills.append(result['spillover']);lambdas.append(result['inflation'])
            spatial.extend(result['spatial']);powers.extend(result['power']);fprs.extend(result['fpr'])
            loci[(method,ph)]={tuple(k) for k in result['strict_loci']}
    count=pd.DataFrame(counts);spill=pd.DataFrame(spills);lam=pd.DataFrame(lambdas)
    output_table(pd.DataFrame(inputs),out/'input_manifest.tsv')
    output_table(count,out/'method_summary/significant_snps_by_pheno.tsv')
    output_table(count[['method','phenotype','n_significant_loci','inclusive_loci']],out/'method_summary/significant_loci_by_pheno.tsv')
    output_table(aggregate(count,['method'],{'sig_snps':'n_significant_snps','sig_loci':'n_significant_loci'}),out/'method_summary/method_summary.tsv')
    output_table(spill,out/'spillover_summary/spillover_by_pheno.tsv')
    output_table(aggregate(spill,['method'],{'spillover_rate':'spillover_rate','spillover_loci':'n_spillover_loci'}),out/'spillover_summary/spillover_summary.tsv')
    output_table(lam,out/'lambda_gc_summary/lambda_gc_by_pheno.tsv')
    output_table(aggregate(lam,['method'],{'lambda_gc':'lambda_gc'}),out/'lambda_gc_summary/lambda_gc_summary.tsv')
    for name,data,value,groups in [('power',powers,'power',['method','alpha']),('fpr',fprs,'fpr',['method','region','alpha'])]:
        frame=pd.DataFrame(data);output_table(frame,out/f'power_summary/{name}_by_pheno.tsv')
        output_table(aggregate(frame,groups,{value:value}),out/f'power_summary/{name}_summary.tsv')
    sp=pd.DataFrame(spatial);output_table(sp,out/'power_summary/spatial_enrichment_by_pheno.tsv')
    output_table(aggregate(sp,['method','distance_index'],{'observed_expected':'observed_expected'}),out/'power_summary/spatial_enrichment_summary.tsv')
    overlap=[];stability=[]
    for ph in phenotypes:
        for a in methods:
            for b in methods:
                A,B=loci[a,ph],loci[b,ph]
                overlap.append(dict(phenotype=ph,method_a=a,method_b=b,n_loci_a=len(A),n_loci_b=len(B),n_intersection=len(A&B),n_union=len(A|B),jaccard=jaccard(A,B)))
    for m in methods:
        for a in phenotypes:
            for b in phenotypes:
                A,B=loci[m,a],loci[m,b]
                stability.append(dict(method=m,phenotype_a=a,phenotype_b=b,n_intersection=len(A&B),n_union=len(A|B),jaccard=jaccard(A,B)))
    ov=pd.DataFrame(overlap)
    output_table(ov,out/'overlap_summary/locus_overlap_by_pheno.tsv')
    output_table(aggregate(ov,['method_a','method_b'],{'jaccard':'jaccard'}),out/'overlap_summary/locus_overlap_summary.tsv')
    output_table(pd.DataFrame(stability),out/'overlap_summary/cross_replicate_stability.tsv')
    (out/'run_manifest.json').write_text(json.dumps(dict(status='AGGREGATES_REGENERATED_NOT_YET_MANUSCRIPT_VALIDATED',n_jobs=len(inputs),code_sha256=sha256(__file__),config_sha256=config_hash,expected_snps=expected_snps,central_windows=4*len(nulls)),indent=2)+'\n')

def equal_count_bins(values,n_bins=20):
    v=np.asarray(values)
    if len(v)<n_bins or np.any(~np.isfinite(v)) or np.any(v<0):raise ValueError('Invalid prediction values or too few SNPs')
    order=np.argsort(v,kind='stable');bins=np.empty(len(v),dtype=int)
    bins[order]=np.floor(np.arange(len(v))*n_bins/len(v)).astype(int)+1
    return bins

def equation13(root,theory_file,out):
    protect_source(root,out)
    root,out=Path(root),Path(out)
    theory=pd.read_csv(theory_file,sep='\t')
    required={'chr','pos','causal_ldscore','c_loco','ncp_framework','expected_chi2','is_central_null'}
    if not required.issubset(theory):raise ValueError('Missing cached-theory columns')
    if theory.duplicated(['chr','pos']).any():raise ValueError('Duplicate theoretical position')
    algebra=np.abs(theory.ncp_framework-theory.c_loco*theory.causal_ldscore)
    expected=np.abs(theory.expected_chi2-1-theory.ncp_framework)
    if algebra.max()>1e-8 or expected.max()>1e-8:raise ValueError('Cached theory does not satisfy manuscript Eq13 simplification')
    theory['bin']=equal_count_bins(theory.ncp_framework)
    central=theory.is_central_null.astype(str).str.lower().isin(['true','1'])
    tc=theory.loc[central].sort_values(['ncp_framework','snp'],kind='stable').copy()
    tc['bin']=equal_count_bins(tc.ncp_framework)
    records=[];pred=[]
    for region,frame in [('Genome-wide',theory),('Central-null',tc)]:
        pr=frame.groupby('bin').agg(n_obs=('ncp_framework','size'),mean_predicted_ncp=('ncp_framework','mean'),predicted_chisq=('expected_chi2','mean')).reset_index()
        pr['region']=region;pred.append(pr)
    for method in ['BOLT-LMM','REGENIE_LOCO','REGENIE_woLOCO','SLEMM']:
        for ph in range(1,11):
            d=harmonize(input_path(root,method,ph),method,ph,135358)
            for region,frame in [('Genome-wide',theory),('Central-null',tc)]:
                merged=frame[['chr','pos','bin']].merge(d[['chr','bp','chisq']],left_on=['chr','pos'],right_on=['chr','bp'],how='left',validate='one_to_one')
                if len(merged)!=len(frame) or not np.isfinite(merged.chisq).all():raise ValueError('Incomplete Eq13 statistic join')
                s=merged.groupby('bin').agg(n_obs=('chisq','size'),mean_observed_chisq=('chisq','mean')).reset_index()
                s['method']=method;s['phenotype']=ph;s['region']=region;records.extend(s.to_dict('records'))
            print(f'Equation13: {method} phenotype {ph}',flush=True)
    per=pd.DataFrame(records);pred=pd.concat(pred,ignore_index=True)
    summary=per.groupby(['method','region','bin']).agg(mean_observed_chisq=('mean_observed_chisq','mean'),sd_observed_chisq=('mean_observed_chisq','std'),n_phenotypes=('phenotype','nunique')).reset_index()
    summary=summary.merge(pred,on=['region','bin'],validate='many_to_one')
    error=t.ppf(.975,9)*summary.sd_observed_chisq/np.sqrt(10)
    summary['ci_lower']=summary.mean_observed_chisq-error;summary['ci_upper']=summary.mean_observed_chisq+error
    per=per.merge(pred.drop(columns='n_obs'),on=['region','bin'],validate='many_to_one')
    output_table(per,out/'equation13_by_pheno.tsv');output_table(summary,out/'equation13_summary.tsv')
    output_table(pd.DataFrame([dict(n_theory_snps=len(theory),max_ncp_algebra_error=algebra.max(),max_expectation_error=expected.max(),theory_sha256=sha256(theory_file),upstream_equation12_validation='NOT_RUN')]),out/'cached_theory_validation.tsv')

def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='command',required=True)
    s=sub.add_parser('summarize');s.add_argument('--root',required=True);s.add_argument('--out',required=True);s.add_argument('--cache',required=True)
    e=sub.add_parser('equation13');e.add_argument('--root',required=True);e.add_argument('--theory',required=True);e.add_argument('--out',required=True)
    h=sub.add_parser('harmonize');h.add_argument('--input',required=True);h.add_argument('--method',choices=METHODS,required=True);h.add_argument('--phenotype',type=int,required=True);h.add_argument('--out',required=True);h.add_argument('--expected-snps',type=int,default=135358)
    args=ap.parse_args()
    if args.command=='summarize':summarize(args.root,args.out,args.cache)
    elif args.command=='equation13':equation13(args.root,args.theory,args.out)
    else:output_table(harmonize(args.input,args.method,args.phenotype,args.expected_snps),args.out)
if __name__=='__main__':main()
