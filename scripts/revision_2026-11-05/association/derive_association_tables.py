#!/usr/bin/env python3
"""Derive one set of association tables without changing any source results.

All windows are zero-anchored [start,end), floor(position/500000), and clipped
to the saved simulation chromosome lengths. Central null is the middle 2 Mb
of complete 4 Mb null parents. No individual-level values are written.
"""
import argparse, hashlib, json, platform, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import chi2
if __package__:
    from .significance import ALPHA, RULE, significant, strict_for_historical_audit
else:
    from significance import ALPHA, RULE, significant, strict_for_historical_audit

METHODS = ['BOLT-LMM','REGENIE_LOCO','REGENIE_WO_LOCO','fastGWA',
           'FarmCPU_NYC','FarmCPU_Kansas','BLINK_NYC','BLINK_Kansas','SLEMM','GCTA_MLMA']
FAMILY = ['BOLT-LMM','REGENIE','REGENIE','fastGWA','FarmCPU','FarmCPU','BLINK','BLINK','SLEMM','GCTA']
P = ALPHA

def sha(path):
    h=hashlib.sha256()
    with open(str(path),'rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

def save(d,out,name):
    d.to_csv(out/name,sep='\t',index=False,na_rep='NA',float_format='%.15g')

def overlap_summary(overlaps):
    rows=[]
    for (a,b),g in overlaps.groupby(['configuration_a','configuration_b'],sort=False):
        complete=g.phenotype.nunique()==10
        rows.append(dict(configuration_a=a,configuration_b=b,
            n_paired_phenotypes=g.phenotype.nunique(),n_defined_phenotypes=int(g.jaccard.notna().sum()),
            complete_ten_phenotypes=complete,
            mean_jaccard=g.jaccard.mean() if complete else np.nan,
            sd_jaccard=g.jaccard.std(ddof=1) if complete else np.nan,
            pooled_intersection=int(g.intersection.sum()),pooled_union=int(g.union.sum()),
            pooled_jaccard=g.intersection.sum()/g.union.sum() if complete and g.union.sum() else np.nan))
    return pd.DataFrame(rows)

def assign(ch,pos,grid):
    ans=np.full(len(pos),-1,dtype=int)
    for c in np.unique(ch):
        ii=np.flatnonzero(ch==c);g=grid[grid.chr==c]
        if g.empty:continue
        k=np.searchsorted(g.start.values,pos[ii],side='right')-1
        valid=(k>=0)&(k<len(g));ix=ii[valid];k=k[valid]
        within=pos[ix]<g.end.values[k]
        ans[ix[within]]=g.window_id.values[k[within]]
    return ans

def grid_from_sources(bench,sim):
    effects=pd.read_csv(bench/'eligible_regions.csv').astype('int64').sort_values(['chr','start'])
    lens=pd.read_csv(sim/'chr1_18_lengths.csv').set_index('chr').length.to_dict()
    parents=[];rows=[]
    for ch,g in effects.groupby('chr'):
        cursor=0
        for r in g.itertuples():
            if cursor<r.start:parents.append((ch,cursor,r.start,'null'))
            parents.append((ch,r.start,r.end,'effect'));cursor=r.end
        if cursor<lens[ch]:parents.append((ch,cursor,lens[ch],'null'))
    for ch,st,en,typ in parents:
        full=en-st==(2000000 if typ=='effect' else 4000000)
        for k,s in enumerate(range(st,en,500000),1):
            rows.append(dict(chr=ch,start=s,end=min(s+500000,en),region_type=typ,
                parent_region_id='chr%d_%s_%d_%d'%(ch,typ,st,en),window_index=k,
                central_null=bool(typ=='null' and full and 3<=k<=6),
                parent_complete=full,window_complete=min(s+500000,en)-s==500000))
    grid=pd.DataFrame(rows).sort_values(['chr','start']).reset_index(drop=True)
    grid['window_id']=np.arange(len(grid))
    assert set(grid.chr)==set(range(1,19))
    assert (grid.start%500000==0).all()
    for ch,g in grid.groupby('chr'):
        assert g.start.iloc[0]==0 and g.end.iloc[-1]==lens[ch]
        assert np.array_equal(g.end.values[:-1],g.start.values[1:])
    cached=pd.read_csv(bench/'reviewer_peak_window_analysis/results/window_grid_500kb.tsv',sep='\t',keep_default_na=False)
    for c in grid:assert np.array_equal(grid[c],cached[c]), 'historical grid differs: '+c
    return grid,pd.DataFrame(parents,columns=['chr','start','end','region_type'])

def config(bench,m,ph,gcta):
    if m=='BOLT-LMM':return bench/('BOLTLMM/bolt_pheno%d.stats'%ph),r'\s+',['CHR','BP','SNP','P_BOLT_LMM_INF','BETA','SE'],False
    if m.startswith('REGENIE'):
        d='REGENIE' if m=='REGENIE_LOCO' else 'REGENIE_WO_LOCO'
        pref='regenie_pheno' if d=='REGENIE' else 'regenie_pheno_wo_loco'
        return bench/('%s/%s%d_pheno%d.regenie'%(d,pref,ph,ph)),r'\s+',['CHROM','GENPOS','ID','LOG10P','BETA','SE'],True
    if m=='fastGWA':return bench/('fastGWA/spg.even.100k.pheno%d.fastgwa.fastGWA'%ph),r'\s+',['CHR','POS','SNP','P','BETA','SE'],False
    if m in ['SLEMM','SLEMM_previous']:
        d='SLEMM_qf90_window3000' if m=='SLEMM' else 'SLEMM'
        return bench/('%s/spg.even.100k.pheno%d.allchr.gwas.txt'%(d,ph)),r'\s+',['#CHROM','POS','ID','P','BETA','SE'],False
    if m=='GCTA_MLMA':return gcta/('spg.even.100k.pheno%d.mlma'%ph),r'\s+',['Chr','bp','SNP','p','b','se'],False
    model,city=m.split('_')
    filename='GAPIT.Association.GWAS_Results.%s.pheno%d(%s).csv'%(model,ph,city)
    candidates=[bench/'GAPIT'/filename,bench/model/filename]
    found=[p for p in candidates if p.is_file()]
    if len(found)>1 and len(set(sha(p) for p in found))>1:
        raise ValueError('Conflicting GAPIT copies; identify authoritative input: '+str(found))
    return (found[0] if found else candidates[0]),',',['Chr','Pos','SNP','P.value','Effect' if city=='NYC' else 'effect',None],False

def sample_qc(bench,sim,geno,qc):
    fam=pd.read_csv(str(geno)+'.fam',sep=r'\s+',header=None,dtype=str)
    ph=pd.read_csv(bench/'SLEMM/pheno.slemm.csv',dtype={'IID':str})
    assert len(ph)==20000 and ph.IID.is_unique
    assert fam.iloc[:,1].is_unique and ph.IID.isin(fam.iloc[:,1]).all()
    sample_ix=pd.Series(np.arange(len(fam)),index=fam.iloc[:,1]).loc[ph.IID].values
    delta=[]
    for i in range(1,11):
        x=pd.read_csv(sim/('data_20000/pheno.%d.csv'%i),dtype={'ID':str})
        assert np.array_equal(ph.IID,x.ID)
        dif=np.abs(ph['pheno%d'%i]-x.pheno)
        assert dif.max()==0
        delta.append(dict(phenotype=i,n_samples=20000,ordered_IDs_match=True,max_phenotype_abs_difference=dif.max()))
    save(pd.DataFrame(delta),qc,'sample_phenotype_qc.tsv')
    bim=pd.read_csv(str(geno)+'.bim',sep=r'\s+',header=None,names=['chr','snp','cm','pos','a1','a2'])
    assert bim.snp.is_unique and len(bim)==135358
    info=pd.read_csv(bench.parent/'W5Mb_N10k/SLEMM/snp_info.csv')
    assert np.array_equal(info.SNP,bim.snp)
    # BED calls are independently checked by the private provenance workflow.
    existing=qc/'sample_marker_qc.json'
    if existing.is_file():
        previous=json.loads(existing.read_text())
        if previous.get('ordered_analysis_IID_sha256')==hashlib.sha256(('\n'.join(ph.IID)+'\n').encode()).hexdigest():return bim
    bpr=(len(fam)+3)//4
    with open(str(geno)+'.bed','rb') as f:assert f.read(3)==b'\x6c\x1b\x01'
    raw=np.memmap(str(geno)+'.bed',dtype='uint8',mode='r',offset=3,shape=(len(bim),bpr))
    missing=0;monomorphic=0
    for start in range(0,len(bim),1000):
        z=raw[start:start+1000,:]
        codes=np.stack([(z>>s)&3 for s in (0,2,4,6)],axis=2).reshape(len(z),-1)[:,sample_ix]
        missing+=int((codes==1).sum())
        hom0=(codes==0).sum(1);hom2=(codes==3).sum(1);het=(codes==2).sum(1)
        ac=2*hom2+het;nc=hom0+hom2+het
        monomorphic+=int(((ac==0)|(ac==2*nc)).sum())
    assert monomorphic==0
    summary=dict(panel_samples=len(fam),analysis_samples=len(ph),markers=len(bim),
        ordered_analysis_IID_sha256=hashlib.sha256(('\n'.join(ph.IID)+'\n').encode()).hexdigest(),
        first_20000_panel_order=bool(np.array_equal(ph.IID.values,fam.iloc[:20000,1].values)),
        first_20000_panel_same_ID_set=bool(set(ph.IID)==set(fam.iloc[:20000,1])),
        phenotype_IDs_matched_to_genotype_panel=True,missing_genotype_calls_in_analysis_subset=missing,
        monomorphic_markers_in_subset=monomorphic,phenotype_missing=int(ph.isna().sum().sum()),
        slemm_selected_marker_IDs_order_equal_BIM=True)
    (qc/'sample_marker_qc.json').write_text(json.dumps(summary,indent=2)+'\n')
    return bim

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--benchmark',type=Path,required=True)
    ap.add_argument('--simulation',type=Path,required=True);ap.add_argument('--genotype-prefix',type=Path,required=True)
    ap.add_argument('--gcta-directory',type=Path);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();bench=a.benchmark;out=a.output;gcta=a.gcta_directory or bench/'GCTA_mlma'
    for d in [out/'tables',out/'qc',out/'logs'] :d.mkdir(parents=True,exist_ok=True)
    tables=out/'tables';qc=out/'qc'
    grid,parents=grid_from_sources(bench,a.simulation);bim=sample_qc(bench,a.simulation,a.genotype_prefix,qc)
    save(grid,tables,'window_inventory_500kb.tsv');save(parents,tables,'parent_regions.tsv')
    causal=pd.read_csv(bench/'snp.info.chr1_18.10k.csv').SNP
    assert len(causal)==10000 and causal.is_unique
    sp=causal.str.split('_',expand=True);cw=assign(sp[0].astype(int).values,sp[1].astype(int).values,grid)
    assert (cw>=0).all() and grid.loc[cw,'region_type'].eq('effect').all()
    grid['n_causal']=np.bincount(cw,minlength=len(grid))
    marker_lookup=bim.set_index('snp');rows=[];peaks=[];spatial=[];audit=[];manifest=[];vc=[]
    full_null=parents[(parents.region_type=='null')&(parents.end-parents.start==4000000)]
    nb=[]
    for r in full_null.itertuples():
        for k in range(80):nb.append((r.chr,r.start+50000*k,r.start+50000*(k+1),min(k+1,80-k)))
    ng=pd.DataFrame(nb,columns=['chr','start','end','distance_index']);ng['window_id']=np.arange(len(ng))
    save(ng,tables,'null_inventory_50kb.tsv')
    for m in METHODS+['SLEMM_previous']:
        for ph in range(1,11):
            path,sep,cols,islog=config(bench,m,ph,gcta)
            meta=dict(configuration=m,method_family=FAMILY[METHODS.index(m)] if m in METHODS else 'SLEMM',phenotype=ph,path=str(path))
            if not path.is_file():
                meta['status']='missing';manifest.append(meta);continue
            d=pd.read_csv(path,sep=sep)
            names=dict(zip([c for c in cols if c],[n for c,n in zip(cols,['chr','pos','snp','value','beta','se']) if c]))
            d=d.rename(columns=names)
            assert d.snp.is_unique and not d.duplicated(['chr','pos']).any()
            ids_exact=bool(d.snp.isin(marker_lookup.index).all())
            d['source_snp']=d.snp
            if not ids_exact:
                coordinate_lookup=bim.set_index(['chr','pos'])
                marker=coordinate_lookup.loc[list(zip(d.chr,d.pos))]
                d['snp']=marker.snp.values
            marker=marker_lookup.loc[d.snp]
            assert np.array_equal(d.chr.values,marker.chr.values) and np.array_equal(d.pos.values,marker.pos.values)
            assert set(d.chr)==set(range(1,19))
            aa=None
            for pair in [('REF','ALT'),('ALLELE1','ALLELE0'),('A1','A2'),('ALLELE0','ALLELE1')]:
                if all(c in d for c in pair):aa=pair;break
            alleles_verified=False
            if aa:
                x,y=d[aa[0]].astype(str).values,d[aa[1]].astype(str).values
                b1,b2=marker.a1.values,marker.a2.values
                assert (((x==b1)&(y==b2))|((x==b2)&(y==b1))).all()
                alleles_verified=True
            v=d.value.values;valid=np.isfinite(v)
            assert (v[valid]>=0).all() and (islog or (v[valid]<=1).all())
            wid=assign(d.chr.values,d.pos.values,grid);assert (wid>=0).all()
            sig=significant(v,is_log10=islog)
            strict=strict_for_historical_audit(v,is_log10=islog)
            strength=v if islog else -np.log10(np.maximum(v,np.nextafter(0.,1.)))
            pval=10.**(-v) if islog else v
            d['window_id']=wid;d['strength']=strength;d['p']=pval
            if m.startswith(('FarmCPU','BLINK')):d['beta']=np.nan
            s=d[sig].copy();s['lead_is_causal']=s.snp.isin(causal)
            lead=s.sort_values(['strength','pos','snp'],ascending=[False,True,True]).drop_duplicates('window_id')
            base_columns=['chr','start','end','region_type','parent_region_id','window_index','central_null','parent_complete','window_complete','window_id','n_causal']
            lead=grid[base_columns].merge(lead[['window_id','snp','source_snp','pos','p','beta','strength','lead_is_causal']].rename(columns={'snp':'lead_snp','source_snp':'lead_source_snp','pos':'lead_pos','p':'lead_p','beta':'lead_beta'}),on='window_id')
            lead['configuration']=m;lead['phenotype']=ph
            lead['n_significant_snps']=lead.window_id.map(s.groupby('window_id').size())
            peaks.append(lead)
            n=len(lead);ne=int(lead.region_type.eq('effect').sum());nn=int(lead.region_type.eq('null').sum());nc=int(lead.central_null.sum())
            central=grid.central_null.values[wid]&valid
            histx=(d.beta.values[central]/d.se.values[central])**2 if 'se' in d else chi2.ppf(1-pval[central],1)
            if m.startswith(('FarmCPU','BLINK')):histx=chi2.ppf(1-pval[central],1)
            row=dict(configuration=m,phenotype=ph,significant_snps=int(sig.sum()),loci=n,effect_loci=ne,null_loci=nn,central_null_loci=nc,
                strict_significant_snps=int(strict.sum()),strict_loci=len(np.unique(wid[strict])),snps_at_threshold=int((sig&~strict).sum()),
                effect_locus_fraction=ne/n if n else np.nan,null_locus_fraction=nn/n if n else np.nan,central_null_locus_fraction=nc/n if n else np.nan,
                central_null_windows=int(grid.central_null.sum()),central_null_spillover_rate=nc/grid.central_null.sum(),
                lambda_gc_central_null=np.median(chi2.isf(pval[central],1))/chi2.ppf(.5,1),
                lambda_gc_central_null_historical=np.median(histx[np.isfinite(histx)])/chi2.ppf(.5,1),
                n_historical_finite_central_statistics=int(np.isfinite(histx).sum()),
                lambda_gc_all=np.median(chi2.isf(pval[valid],1))/chi2.ppf(.5,1),
                n_central_null_snps=int(central.sum()),causal_leads=int(lead.lead_is_causal.sum()),
                lead_abs_beta_mean=lead.lead_beta.abs().mean(),lead_abs_beta_median=lead.lead_beta.abs().median(),lead_abs_beta_max=lead.lead_beta.abs().max())
            assert ne+nn==n and nc<=nn
            rows.append(row)
            coverage=np.bincount(wid[valid],minlength=len(grid));grid['n_tested_snps_%s_%d'%(m,ph)]=coverage
            assert len(d)==len(bim)
            if m!='SLEMM_previous':
                bid=assign(d.chr.values,d.pos.values,ng);usable=(bid>=0)&valid
                tested=np.bincount(bid[usable],minlength=len(ng));hits=np.bincount(bid[usable&sig],minlength=len(ng))
                # Preserve the original closed-interval foverlaps convention:
                # an exact internal 50 kb boundary belongs to both neighbors.
                on_boundary=usable & (bid>0)
                bi=np.flatnonzero(on_boundary);kk=bid[bi]
                on_boundary[bi]=(d.pos.values[bi]==ng.start.values[kk])&(d.pos.values[bi]==ng.end.values[kk-1])&(d.chr.values[bi]==ng.chr.values[kk-1])
                tested+=np.bincount(bid[on_boundary]-1,minlength=len(ng));hits+=np.bincount(bid[on_boundary&sig]-1,minlength=len(ng))
                ratio=np.divide(hits,tested*P,out=np.full(len(ng),np.nan),where=tested>0)
                z=ng[['distance_index']].copy();z['tested']=tested;z['hits']=hits;z['ratio']=ratio
                for ix,g in z.groupby('distance_index'):
                    spatial.append(dict(configuration=m,phenotype=ph,distance_index=ix,
                        observed_expected_equal_nonempty_bins=g.ratio.mean(),nonempty_bins=int((g.tested>0).sum()),
                        observed_expected_pooled_snps=g.hits.sum()/(g.tested.sum()*P)))
            meta.update(status='derived_complete',sha256=sha(path),bytes=path.stat().st_size,rows=len(d),p_column=cols[3],p_scale='minus_log10_P' if islog else 'P');manifest.append(meta)
            audit.append(dict(configuration=m,phenotype=ph,rows=len(d),missing_p=int((~valid).sum()),zero_p=int((pval==0).sum()),
                chromosomes='1-18',source_marker_ID_set_matches_panel=ids_exact,canonical_marker_ID_set_matches_panel=True,marker_order_matches_panel=np.array_equal(d.snp,bim.snp),
                alleles_verified=alleles_verified,allele_columns='|'.join(aa) if aa else 'not_present',tested_windows=int((coverage>0).sum())))
            if m in ['SLEMM','SLEMM_previous']:
                vf=path.parent/('spg.even.100k.pheno%d.reml.vc.csv'%ph)
                vv=pd.read_csv(vf,header=None).set_index(0)[1].to_dict()
                vc.append(dict(configuration=m,phenotype=ph,**vv,sha256=sha(vf)))
            print('%s phenotype %d: %d SNPs; %d loci'%(m,ph,int(sig.sum()),n),flush=True)
    by=pd.DataFrame(rows);lead=pd.concat(peaks,ignore_index=True);man=pd.DataFrame(manifest)
    save(by,tables,'association_by_phenotype.tsv');save(lead,tables,'detected_loci.tsv');save(pd.DataFrame(spatial),tables,'spatial_by_phenotype.tsv')
    save(man,qc,'input_manifest.tsv');save(pd.DataFrame(audit),qc,'association_qc.tsv');save(pd.DataFrame(vc),tables,'slemm_heritability_by_phenotype.tsv')
    metrics=[c for c in by if c not in ['configuration','phenotype']]
    ss=[]
    for m,g in by.groupby('configuration',sort=False):
        r=dict(configuration=m,n_phenotypes=len(g),complete_ten_phenotypes=len(g)==10)
        for k in metrics:
            r[k+'_mean']=g[k].mean() if len(g)==10 else np.nan;r[k+'_sd']=g[k].std(ddof=1) if len(g)==10 else np.nan;r[k+'_n_defined']=int(g[k].notna().sum())
        for typ in ['effect','null','central_null']:r['pooled_'+typ+'_locus_fraction']=g[typ+'_loci'].sum()/g.loci.sum() if g.loci.sum() and len(g)==10 else np.nan
        r['pooled_loci']=g.loci.sum();r['max_loci']=g.loci.max();r['max_loci_phenotypes']=','.join(g[g.loci==g.loci.max()].phenotype.astype(str));r['zero_loci_phenotypes']=','.join(g[g.loci==0].phenotype.astype(str))
        ss.append(r)
    summary=pd.DataFrame(ss);save(summary,tables,'association_summary.tsv')
    save(summary[summary.configuration != "SLEMM_previous"][['configuration','n_phenotypes','complete_ten_phenotypes','loci_mean','loci_sd','effect_locus_fraction_mean','effect_locus_fraction_sd','effect_locus_fraction_n_defined','null_locus_fraction_mean','null_locus_fraction_sd','central_null_locus_fraction_mean','central_null_locus_fraction_sd','pooled_effect_locus_fraction','pooled_null_locus_fraction','pooled_central_null_locus_fraction']],tables,'supplementary_table_S1.tsv')
    overlaps=[]
    for ph,g in by[by.configuration!='SLEMM_previous'].groupby('phenotype'):
        ms=list(g.configuration)
        sets={m:set(lead[(lead.configuration==m)&(lead.phenotype==ph)].window_id) for m in ms}
        for m in ms:
            for n in ms:
                x,y=sets[m],sets[n]
                overlaps.append(dict(phenotype=ph,configuration_a=m,configuration_b=n,n_a=len(x),n_b=len(y),intersection=len(x&y),union=len(x|y),
                    jaccard=len(x&y)/len(x|y) if x and y else np.nan,
                    jaccard_set_theoretic=len(x&y)/len(x|y) if x|y else np.nan,
                    historical_jaccard_defined=bool(x and y)))
    overlaps=pd.DataFrame(overlaps)
    save(overlaps,tables,'pairwise_locus_overlap.tsv')
    save(overlap_summary(overlaps),tables,'pairwise_locus_overlap_summary.tsv')
    bg=[]
    for typ,mask in [('effect',grid.region_type=='effect'),('null',grid.region_type=='null'),('central_null',grid.central_null)]:
        tested=grid['n_tested_snps_BOLT-LMM_1']>0
        bg.append(dict(region=typ,n_windows=int(mask.sum()),eligible_windows=len(grid),fraction_windows=mask.mean(),
            tested_windows=int((mask&tested).sum()),total_tested_windows=int(tested.sum()),fraction_tested_windows=(mask&tested).sum()/tested.sum(),
            bp=int((grid.end[mask]-grid.start[mask]).sum()),genome_bp=int((grid.end-grid.start).sum()),fraction_bp=(grid.end[mask]-grid.start[mask]).sum()/(grid.end-grid.start).sum()))
    save(pd.DataFrame(bg),tables,'genomic_background.tsv')
    old=by[by.configuration=='SLEMM_previous'].set_index('phenotype');new=by[by.configuration=='SLEMM'].set_index('phenotype')
    comp=[]
    for ph in range(1,11):
        for k in metrics:comp.append(dict(phenotype=ph,statistic=k,previous_value=old.loc[ph,k],new_value=new.loc[ph,k],difference=new.loc[ph,k]-old.loc[ph,k]))
    save(pd.DataFrame(comp),tables,'slemm_old_new_reconciliation.tsv')
    (qc/'run_metadata.json').write_text(json.dumps(dict(command=sys.argv,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,
        completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),methods=len(set(FAMILY)),configurations=len(METHODS),
        fully_derived_configurations=int(summary[summary.configuration!='SLEMM_previous'].complete_ten_phenotypes.sum()),
        significance=RULE,comparison_scale='native source P or native minus-log10 P; no tolerance widening',intervals='[start,end), floor(position/500000), clipped at saved chromosome length',
        jaccard='historical intersection/union with NA if either set empty; set-theoretic alternative separately labeled'),indent=2)+'\n')

if __name__=='__main__':main()
