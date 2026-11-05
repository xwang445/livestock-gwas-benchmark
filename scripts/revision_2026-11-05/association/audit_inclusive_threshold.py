#!/usr/bin/env python3
"""Independent audit of reported-scale inclusive thresholds and canonical tables.

Reads source P tokens directly and reconstructs decisions independently of the
derivation's readers, then checks the shared predicate against that reconstruction.
Writes QC records only; no analysis outputs are changed.
"""
import argparse, json, math
from decimal import Decimal
from pathlib import Path
import numpy as np
import pandas as pd
from significance import significant

ALPHA=5e-8

def read(path):return pd.read_csv(path,sep='\t',keep_default_na=False,na_values=['NA'],float_precision='round_trip')

def fields(m):
    if m=='BOLT-LMM':return ['CHR','BP','P_BOLT_LMM_INF']
    if m.startswith('REGENIE'):return ['CHROM','GENPOS','LOG10P']
    if m=='fastGWA':return ['CHR','POS','P']
    if m.startswith('SLEMM'):return ['#CHROM','POS','P']
    if m=='GCTA_MLMA':return ['Chr','bp','p']
    return ['Chr','Pos','P.value']

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--analysis',type=Path,required=True);ap.add_argument('--benchmark',type=Path,required=True);ap.add_argument('--manhattan-manifest',type=Path,required=True);a=ap.parse_args()
    tables=a.analysis/'tables';qc=a.analysis/'qc';manifest=read(qc/'input_manifest.tsv');canonical=read(tables/'association_by_phenotype.tsv').set_index(['configuration','phenotype'])
    grid=read(tables/'window_inventory_500kb.tsv');lookup={(int(r.chr),int(r.start//500000)):(int(r.window_id),r.region_type,bool(r.central_null)) for r in grid.itertuples()}
    causal=pd.read_csv(a.benchmark/'snp.info.chr1_18.10k.csv').SNP.str.split('_',expand=True);truth=set(zip(causal[0].astype(int),causal[1].astype(int)))
    power=read(tables/'power_and_null_snp_rate_by_phenotype.tsv')
    cutoffs=power[['alpha','minus_log10_alpha','alpha_hex','alpha_expression']].drop_duplicates()
    serialization_disagreements=int(sum(x.alpha!=float.fromhex(x.alpha_hex) for x in cutoffs.itertuples()))
    power=power.set_index(['configuration','phenotype','alpha_hex','region'])
    tau=-math.log10(ALPHA);rows=[];sets={};failures=[];strict_details=[];power_rows=0
    for r in manifest[manifest.status=='derived_complete'].itertuples():
        cols=fields(r.configuration);gap=r.configuration.startswith(('BLINK','FarmCPU'))
        raw=pd.read_csv(r.path,sep=',' if gap else r'\s+',usecols=cols,dtype={cols[-1]:str})
        v=pd.to_numeric(raw[cols[-1]],errors='coerce').values;valid=np.isfinite(v);islog=r.p_scale=='minus_log10_P'
        p=np.power(10.,-v) if islog else v
        native=valid&(v>=tau if islog else v<=ALPHA);strict=valid&(v>tau if islog else v<ALPHA)
        if not np.array_equal(native,significant(v,is_log10=islog)):failures.append('Shared helper primary '+r.configuration+' '+str(r.phenotype))
        pp=valid&(p<=ALPHA);logp=v if islog else -np.log10(np.maximum(p,np.nextafter(0.,1.)))
        lp=valid&(logp>=tau)
        ch=raw[cols[0]].astype(int).values;pos=raw[cols[1]].astype(int).values
        keys=[(int(c),int(b//500000)) for c,b in zip(ch,pos)];assigned=[lookup[k] for k in keys]
        wid=np.array([x[0] for x in assigned]);effect=np.array([x[1]=='effect' for x in assigned]);central=np.array([x[2] for x in assigned])
        lo=set(wid[native]);ls=set(wid[strict]);c=canonical.loc[r.configuration,r.phenotype]
        check=(len(lo)==c.loci and int(native.sum())==c.significant_snps and len(ls)==c.strict_loci and int(strict.sum())==c.strict_significant_snps)
        if not check:failures.append('primary count '+r.configuration+' '+str(r.phenotype))
        cn=len(set(wid[native&central]));ef=len(set(wid[native&effect]));nu=len(set(wid[native&~effect]))
        if (cn,ef,nu)!=(c.central_null_loci,c.effect_loci,c.null_loci):failures.append('regional counts '+r.configuration+' '+str(r.phenotype))
        if not np.array_equal(native,pp):failures.append('P roundtrip primary '+r.configuration+' '+str(r.phenotype))
        if not np.array_equal(native,lp):failures.append('P log transform primary '+r.configuration+' '+str(r.phenotype))
        decimal_discordant=0
        if not islog:
            near=valid&(np.abs(v-ALPHA)<=ALPHA*1e-8)
            for ix in np.flatnonzero(near):
                exact=Decimal(str(raw.iloc[ix][cols[-1]]))<=Decimal('5e-8')
                decimal_discordant+=int(exact!=native[ix])
            if decimal_discordant:failures.append('Decimal P boundary '+r.configuration+' '+str(r.phenotype))
        scale_discordant=0;causal_mask=np.array([(int(c),int(b)) in truth for c,b in zip(ch,pos)])
        masks={'causal_positions':causal_mask,'central_null':central,'all_null':~effect}
        for cutoff in cutoffs.itertuples():
            alpha=float.fromhex(cutoff.alpha_hex)
            use=valid&(v>=cutoff.minus_log10_alpha if islog else v<=alpha)
            if not np.array_equal(use,significant(v,is_log10=islog,alpha=alpha)):failures.append('Shared helper cutoff '+r.configuration+' '+str(r.phenotype)+' '+cutoff.alpha_expression)
            converted=valid&(p<=alpha)
            scale_discordant+=int((use!=converted).sum())
            for region,mask in masks.items():
                reported=power.loc[r.configuration,r.phenotype,cutoff.alpha_hex,region]
                if int((mask&use).sum())!=reported.n_significant or int((mask&valid).sum())!=reported.n_tested:failures.append('power rate '+r.configuration+' '+str(r.phenotype)+' '+str(alpha)+' '+region)
                power_rows+=1
        row=dict(configuration=r.configuration,phenotype=r.phenotype,rows=len(raw),significant_snps=int(native.sum()),strict_significant_snps=int(strict.sum()),exact_threshold_added_snps=int((native&~strict).sum()),
            inclusive_loci=len(lo),strict_loci=len(ls),additional_loci_at_equality=len(lo-ls),inclusive_effect_loci=ef,inclusive_null_loci=nu,inclusive_central_null_loci=cn,
            additional_effect_loci_at_equality=len(set(wid[native&effect])-set(wid[strict&effect])),additional_null_loci_at_equality=len(set(wid[native&~effect])-set(wid[strict&~effect])),additional_central_loci_at_equality=len(set(wid[native&central])-set(wid[strict&central])),
            P_roundtrip_primary_disagreements=int((native!=pp).sum()),P_to_log_primary_disagreements=int((native!=lp).sum()),decimal_boundary_disagreements=decimal_discordant,cutoff_grid_roundtrip_disagreements=scale_discordant,
            nearest_reported_minus_log10P=float(logp[np.argmin(abs(logp-tau))]),minimum_absolute_log_distance_to_cutoff=float(np.nanmin(abs(logp-tau))))
        rows.append(row);sets[r.configuration,r.phenotype]=lo
        strict_details.append(dict(configuration=r.configuration,phenotype=r.phenotype,strict_effect_loci=len(set(wid[strict&effect])),strict_null_loci=len(set(wid[strict&~effect])),strict_central_loci=len(set(wid[strict&central]))))
    d=pd.DataFrame(rows);d.to_csv(qc/'inclusive_threshold_by_input.tsv',sep='\t',index=False,float_format='%.15g')
    overlaps=read(tables/'pairwise_locus_overlap.tsv')
    for r in overlaps.itertuples():
        x,y=sets[r.configuration_a,r.phenotype],sets[r.configuration_b,r.phenotype]
        if (len(x),len(y),len(x&y),len(x|y))!=(r.n_a,r.n_b,r.intersection,r.union):failures.append('Figure5 locus sets')
    man=read(a.manhattan_manifest);md=man[man.status=='COMPLETED'].merge(d,on=['configuration','phenotype'],validate='one_to_one')
    if len(md)!=len(man[man.status=='COMPLETED']) or not (md.n_significant==md.significant_snps).all():failures.append('Manhattan counts')
    s1=read(tables/'supplementary_table_S1.tsv');s1_checks=0
    for r in s1[s1.complete_ten_phenotypes].itertuples():
        g=canonical.loc[r.configuration]
        for key,calculated in [('loci_mean',g.loci.mean()),('effect_locus_fraction_mean',(g.effect_loci/g.loci.replace(0,np.nan)).mean()),('null_locus_fraction_mean',(g.null_loci/g.loci.replace(0,np.nan)).mean()),('central_null_locus_fraction_mean',(g.central_null_loci/g.loci.replace(0,np.nan)).mean())]:
            if not np.isclose(getattr(r,key),calculated,rtol=1e-12,atol=1e-12):failures.append('S1 '+r.configuration+' '+key)
            s1_checks+=1
    if serialization_disagreements:failures.append('Lossless alpha serialization')
    summary=dict(status='PASS' if not failures else 'FAIL',threshold_P=ALPHA,threshold_minus_log10P=tau,comparison_rule_P='<=',comparison_rule_minus_log10P='>=',n_GWAS_outputs=len(d),n_primary_panels_checked=len(md),n_pairwise_comparisons_checked=len(overlaps),n_power_rate_rows_checked=power_rows,n_S1_statistics_checked=s1_checks,
        exact_threshold_added_SNP_records=int(d.exact_threshold_added_snps.sum()),additional_loci_at_equality=int(d.additional_loci_at_equality.sum()),
        primary_log_or_P_roundtrip_disagreements=int(d.P_roundtrip_primary_disagreements.sum()+d.P_to_log_primary_disagreements.sum()),
        sensitivity_grid_P_roundtrip_diagnostic_disagreements=int(d.cutoff_grid_roundtrip_disagreements.sum()),
        serialized_alpha_roundtrip_disagreements=serialization_disagreements,decimal_boundary_disagreements=int(d.decimal_boundary_disagreements.sum()),failures=failures,
        limitation='Classifies reported P/LOG10P. Original unrounded REGENIE P cannot be recovered from rounded LOG10P; no artificial tolerance or epsilon is added.')
    (qc/'inclusive_threshold_audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
    if failures:raise SystemExit(1)

if __name__=='__main__':main()
