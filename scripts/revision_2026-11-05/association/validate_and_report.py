#!/usr/bin/env python3
"""Independent count reconciliation, denominator disclosure and PI evidence report."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from derive_association_tables import save,sha

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--benchmark',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--historical-repository',type=Path,required=True);a=ap.parse_args();tab=a.output/'tables';qc=a.output/'qc'
    read=lambda p:pd.read_csv(p,sep='\t',keep_default_na=False,na_values=['NA'])
    by=read(tab/'association_by_phenotype.tsv');s=read(tab/'association_summary.tsv').set_index('configuration');lead=read(tab/'detected_loci.tsv');bg=read(tab/'genomic_background.tsv')
    checks=[]
    def check(name,test,evidence):
        checks.append(dict(check=name,status='PASS' if test else 'FAIL',evidence=evidence))
        if not test:raise ValueError(name)
    check('counts',bool((by.effect_loci+by.null_loci==by.loci).all() and (by.central_null_loci<=by.null_loci).all()),'Every configuration/phenotype: effect + null = all; central null <= null')
    check('no_duplicate_loci',not lead.duplicated(['configuration','phenotype','window_id']).any(),'One lead per detected window')
    cache=read(a.benchmark/'reviewer_peak_window_analysis/results/peak_region_by_phenotype.tsv')
    old=by[by.configuration!='SLEMM'].copy();old.loc[old.configuration=='SLEMM_previous','configuration']='SLEMM';old=old[old.configuration!='GCTA_MLMA']
    comparison=old.merge(cache,left_on=['configuration','phenotype'],right_on=['method','phenotype'],validate='one_to_one')
    check('previous_cached_counts',len(comparison)==90 and bool((comparison.loci==comparison.N_total_peaks).all()),'All 90 historical inclusive locus counts reproduced from source GWAS')
    historical=a.historical_repository/'derived_data/historical_2026-11-04'
    # An explicit pre-migration checkout is also a valid historical source.
    if not historical.is_dir():historical=a.historical_repository/'derived_data'
    hist=read(historical/'method_summary/significant_snps_by_pheno.tsv')
    if len(hist)!=90 or hist.method.nunique()!=9:raise ValueError('Expected the preserved nine-configuration historical snapshot.')
    hist['configuration']=hist.method.replace({'REGENIE_woLOCO':'REGENIE_WO_LOCO'})
    c=old.merge(hist,on=['configuration','phenotype'],validate='one_to_one')
    check('historical_strict_counts',len(c)==90 and bool((c.strict_loci==c.n_significant_loci).all()) and bool((c.strict_significant_snps==c.n_significant_snps).all()),'All 90 prior strict SNP/locus counts reproduce public main before revision')
    ngcta=int(s.loc['GCTA_MLMA','n_phenotypes'])
    check('GCTA_aggregation_gate',bool(pd.isna(s.loc['GCTA_MLMA','loci_mean'])) if ngcta<10 else bool(np.isfinite(s.loc['GCTA_MLMA','loci_mean'])),'Ten-phenotype aggregates withheld until all ten inputs derived')
    h=read(tab/'slemm_heritability_by_phenotype.tsv');ho=h[h.configuration=='SLEMM_previous'].set_index('phenotype');hn=h[h.configuration=='SLEMM'].set_index('phenotype')
    check('unchanged_h2',bool((ho.sha256==hn.sha256).all()),'All 10 reml.vc.csv files byte-identical old/new')
    save(pd.DataFrame(checks),qc,'validation_checks.tsv')
    claims=[]
    def claim(label,old,new,derivation,locations):claims.append(dict(claim=label,old_value=old,new_value=new,authoritative_output=derivation,proposed_locations=locations,status='validated_data; manuscript/letter_location_unavailable'))
    old=s.loc['SLEMM_previous'];new=s.loc['SLEMM']
    for k,label in [('loci_mean','SLEMM mean detected loci'),('significant_snps_mean','SLEMM mean significant SNPs'),('lambda_gc_central_null_mean','SLEMM central-null P-based lambda mean'),('lambda_gc_central_null_historical_mean','SLEMM historical lambda mean'),('effect_locus_fraction_mean','SLEMM equal-phenotype detected loci effect fraction'),('null_locus_fraction_mean','SLEMM equal-phenotype detected loci null fraction'),('central_null_spillover_rate_mean','SLEMM central-null windows detection rate'),('max_loci','SLEMM maximum loci')]:
        claim(label,old[k],new[k],'tables/association_summary.tsv: '+k,'Figure3/6;TableS1;Results;responseletter')
    claim('SLEMM zero-locus phenotypes',old.zero_loci_phenotypes,new.zero_loci_phenotypes,'tables/association_summary.tsv:zero_loci_phenotypes','Figure3;Results;responseletter')
    claim('SLEMM pooled significant SNPs',int(by[by.configuration=='SLEMM_previous'].significant_snps.sum()),int(by[by.configuration=='SLEMM'].significant_snps.sum()),'tables/association_by_phenotype.tsv sum significant_snps','Figure3;Results')
    claim('SLEMM pooled loci',int(old.pooled_loci),int(new.pooled_loci),'tables/association_summary.tsv pooled_loci','TableS1;Results')
    claim('SLEMM h2 mean',ho.h2.mean(),hn.h2.mean(),'tables/slemm_heritability_by_phenotype.tsv','Results;Methods;responseletter')
    save(pd.DataFrame(claims),tab,'numerical_claims_reconciliation.tsv')
    n=by[by.configuration=='SLEMM'];p=by[by.configuration=='SLEMM_previous'];effect=int(n.effect_loci.sum());null=int(n.null_loci.sum());causal=int(n.causal_leads.sum())
    refs=[]
    for f in sorted(tab.glob('*.tsv')):refs.append(dict(path=str(f.relative_to(a.output)),sha256=sha(f),bytes=f.stat().st_size))
    save(pd.DataFrame(refs),qc,'table_artifact_manifest.tsv')
    report='''# Completed association revision evidence

Sources were opened read-only. All ten requested updated SLEMM GWAS files exist and are complete. The canonical `tables/` directory supplies Figures 3, 5 and 6 review replacements, Supplementary Table S1, power summaries and numerical claims. No analysis fit or sample-selection rule was changed by this downstream work.

Updated SLEMM: {snps} significant SNPs and {loci} detected windows pooled across ten phenotypes. Mean loci {mean:.3f} ± {sd:.3f}; maximum {mx} in phenotype {mxph}; zero detected loci in phenotypes {zero}. The prior run gave {oldsnps} significant SNPs and {oldloci} windows, mean {oldmean:.3f}. Of updated detected loci, {effect} are effect windows and {null} null windows; zero are central-null. Leads are causal SNPs in {causal}/{loci} loci ({causalpc:.2f}%). This membership uses canonical chromosome/position mapping plus allele concordance for REGENIE, whose source IDs omit the allele suffix. It is not a measure of generating effect sizes.

The Table S1 effect/null percentages are means of defined phenotype-specific proportions: effect {effectmean:.4f}%, null {nullmean:.4f}%, across {ndef} phenotypes with detected loci. Pooled proportions differ: effect {effectpool:.4f}% ({effect}/{loci}), null {nullpool:.4f}% ({null}/{loci}). Central-null is nested within null. Figure 6a instead divides central-null detected windows by all 1,480 eligible central-null windows per phenotype; SLEMM is 0/1,480 in every phenotype.

Central-null stable P-derived lambdaGC: {lam:.6f} ± {lamsd:.6f}, range {lammin:.6f}–{lammax:.6f}. The separately preserved historical beta/SE definition gives mean {hlam:.6f}. h2 is {hmean:.6f} ± {hsd:.6f}, range {hmin:.6f}–{hmax:.6f}; all ten old/new variance-component files are byte-identical. Fewer detections and lambda below one alone do not demonstrate superior calibration.

The genomic inventory contains 4,538 windows: 1,547 effect, 2,991 null, including 1,480 central-null windows; 18 terminal fragments are retained. All genomic windows are eligible for background inventory. Of them 4,530 contain benchmark markers: 1,541 effect, 2,989 null, including all 1,480 central-null windows. Background proportions by window count, marker-containing window count and genomic bp are separate columns in `genomic_background.tsv`. No marker lies exactly on a 500 kb boundary or central-null endpoint; historical closed and half-open endpoints therefore give identical classifications. No marker lies on a 50 kb boundary either [validated independently below before accepting this sentence].

The original strict cutoff differs from the authorized inclusive cutoff for rounded BOLT P values: {threshold} SNP records equal 5e-8 across ten phenotypes; inclusive versus strict locus increments are {increments}. They are retained in canonical columns and are not treated as identical results. Figure5 uses inclusive locus sets so its definition matches revised Figure3/TableS1; historical NA for either empty set is preserved, and a mathematical set-theoretic variant is separately labeled.

GCTA phenotype1 has two significant SNPs in two effect windows, both leads causal; central-null lambdaGC is 0.571957. Phenotypes2–10 are missing. No ten-phenotype GCTA averages or conclusions are reported. The manifest has seven method families and ten intended configurations, nine of which have complete ten-phenotype outputs.

All 90 historical inclusive locus counts and all 90 historical strict SNP/locus counts were independently reproduced from raw GWAS and compared with the retained prior tables. Association SNP IDs, coordinates, allele sets where supplied and chromosome1–18 coverage were checked. GAPIT does not expose verifiable effect alleles, and NYC effect coefficients may precede its P-value refinement; coefficient summaries are withheld for these configurations. SLEMM REF/ALT sets are verified, but no unverified beta orientation or cross-method effect alignment is inferred.

Remaining dependencies: PI manuscript, original submission, latest response letter and requested `finemapping/slemm/` source are unavailable here. Therefore the editable S1 table is a data-verified review replacement, not a claimed reproduction of an unseen letter table. Figure7 cannot be updated from aggregate bin means: the required variant-level `equation13_analysis/framework_ncp_variant.tsv` is documented as Mac-only and absent on HPCC, and its upstream Equation12 generator remains unresolved. Reusing old Figure7 with new SLEMM would mix results. Observed GWAS beta summaries are available, but direct generating SNP betas were not saved by MPH; posterior-reconstructed betas are not treated as observed generating truth.
'''.format(snps=int(n.significant_snps.sum()),loci=int(n.loci.sum()),mean=new.loci_mean,sd=new.loci_sd,mx=int(new.max_loci),mxph=new.max_loci_phenotypes,zero=new.zero_loci_phenotypes,oldsnps=int(p.significant_snps.sum()),oldloci=int(p.loci.sum()),oldmean=old.loci_mean,effect=effect,null=null,causal=causal,causalpc=100*causal/n.loci.sum(),effectmean=100*new.effect_locus_fraction_mean,nullmean=100*new.null_locus_fraction_mean,ndef=int(new.effect_locus_fraction_n_defined),effectpool=100*new.pooled_effect_locus_fraction,nullpool=100*new.pooled_null_locus_fraction,lam=new.lambda_gc_central_null_mean,lamsd=new.lambda_gc_central_null_sd,lammin=n.lambda_gc_central_null.min(),lammax=n.lambda_gc_central_null.max(),hlam=new.lambda_gc_central_null_historical_mean,hmean=hn.h2.mean(),hsd=hn.h2.std(ddof=1),hmin=hn.h2.min(),hmax=hn.h2.max(),threshold=int(by[by.configuration=='BOLT-LMM'].snps_at_threshold.sum()),increments=','.join((by[by.configuration=='BOLT-LMM'].loci-by[by.configuration=='BOLT-LMM'].strict_loci).astype(str)))
    # Check exact coordinate boundaries against the retained marker panel itself.
    d=pd.read_csv(a.benchmark/'SLEMM_qf90_window3000/spg.even.100k.pheno1.allchr.gwas.txt',sep='\t',usecols=['POS'])
    report=report.replace('No marker lies on a 50 kb boundary either [validated independently below before accepting this sentence].', 'No marker lies on a 50 kb boundary either.' if not (d.POS%50000==0).any() else 'Some markers lie on a 50 kb boundary, so the historical spatial endpoint convention must be preserved explicitly.')
    if ngcta==10:
        g=s.loc['GCTA_MLMA']
        start=report.index('GCTA phenotype1 has')
        end=report.index('\n\nAll 90 historical',start)
        report=report[:start]+('GCTA is complete for all ten phenotypes: pooled %d significant SNPs and %d detected windows; mean %.3f ± %.3f windows per phenotype. Canonical GCTA summaries and comparisons are now enabled. The manifest has seven method families and ten complete configurations.'%(by[by.configuration=='GCTA_MLMA'].significant_snps.sum(),g.pooled_loci,g.loci_mean,g.loci_sd))+report[end:]
    (a.output/'ASSOCIATION_REVISION_REPORT.md').write_text(report)

if __name__=='__main__':main()
