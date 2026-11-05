#!/usr/bin/env python3
"""Compare standard RSS fits after changing only residual-variance estimation.

Outputs use no individual-level records. AUPRC is reported as both average
precision (stepwise PR area) and trapezoidal PR area, with its exact universe.
"""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve, auc

p = argparse.ArgumentParser()
p.add_argument('--benchmark', type=Path, required=True)
p.add_argument('--new', type=Path, required=True)
a = p.parse_args()
B, W = a.benchmark.resolve(), a.new.resolve()
OLD = B/'SuSiE_finemapping'
OUT = W/'results'
OUT.mkdir(exist_ok=True)
man = pd.read_csv(B/'SuSiE_adj_finemapping/loci/manifest.tsv', sep='\t')
truth = set(pd.read_csv(B/'snp.info.chr1_18.10k.csv').SNP)
bim = pd.read_csv(B/'SuSiE_adj_finemapping/genotype/spg.20k.bim', sep=r'\s+', header=None,
                  names=['chr','SNP','cm','position','A1','A2'])
assert len(truth) == 10000 and len(bim) == 135358
MDIR = {'BOLT-LMM':'BOLT_LMM','REGENIE-LOCO':'REGENIE'}


def read_fit(row):
    root = W/'susie_results'/MDIR[row.method]/row.phenotype/row.locus_id
    assert (root/'DONE').exists(), root
    f = pd.read_csv(root/'fit_summary.tsv', sep='\t')
    assert len(f) == 1 and f.locus_id.iloc[0] == row.job_id
    s = c = None
    if f.fit_status.iloc[0] == 'ok':
        s = pd.read_csv(root/'snp_pip.tsv', sep='\t')
        assert len(s) == f.matched.iloc[0]
        if (root/'credible_sets.tsv').stat().st_size > 1:
            c = pd.read_csv(root/'credible_sets.tsv', sep='\t')
        assert f.n_cs.iloc[0] == (0 if c is None else len(c))
    return f, s, c


with ThreadPoolExecutor(max_workers=12) as pool:
    parts = list(pool.map(read_fit, man.itertuples()))
F = pd.concat([x[0] for x in parts], ignore_index=True)
S = pd.concat([x[1] for x in parts if x[1] is not None], ignore_index=True)
C = pd.concat([x[2] for x in parts if x[2] is not None], ignore_index=True)
del parts
assert len(F) == len(man) == 13301 and F.locus_id.is_unique
assert F.n_eff.eq(20000).all() and F.fit_mode.eq('ERV_TRUE').all()
assert F.dropped.eq(0).all()
assert S.PIP.between(0,1).all() and not S.duplicated(['locus_id','SNP']).any()
F['job_id'] = F.locus_id
F.to_csv(OUT/'per_locus_diagnostics.tsv', sep='\t', index=False)
S.to_csv(OUT/'snp_level_all.tsv.gz', sep='\t', index=False)
C.to_csv(OUT/'credible_sets_all.tsv', sep='\t', index=False)
oldF = pd.read_csv(OLD/'qc/fit_status_all.tsv', sep='\t')
oldS = pd.read_csv(OLD/'results/snp_level_all.tsv.gz', sep='\t')
oldC = pd.read_csv(OLD/'results/credible_sets_all.tsv', sep='\t')
assert oldF.fit_mode.eq('ERV_FALSE').all()
check = S[['locus_id','SNP','effect_allele','z']].merge(
    oldS[['locus_id','SNP','effect_allele','z']], on=['locus_id','SNP'],
    how='outer', suffixes=('_new','_old'), indicator=True, validate='one_to_one')
assert check._merge.eq('both').all() and check.effect_allele_new.eq(check.effect_allele_old).all()
max_z_diff = float(np.abs(check.z_new-check.z_old).max())
assert max_z_diff < 1e-10
(W/'qc/comparison_input_check.json').write_text(json.dumps(dict(
    n_snp_locus_rows=len(check),max_abs_z_difference=max_z_diff,
    alleles_identical=True,manifest_loci_identical=True,changed_parameter='estimate_residual_variance'),indent=2)+'\n')
del check

common = set(F.loc[(F.fit_status=='ok') & F.converged.eq(True),'locus_id']) & set(
    oldF.loc[(oldF.fit_status=='ok') & oldF.converged.eq(True),'locus_id'])
rows, strata, diag, curves = [], [], [], []
for mode, ff, ss, cc in [('ERV_FALSE',oldF,oldS,oldC),('ERV_TRUE',F,S,C)]:
    for (method,pheno), f in ff.groupby(['method','phenotype']):
        diag.append(dict(fit_mode=mode,method=method,phenotype=pheno,n_attempted=len(f),
            n_ok=int(f.fit_status.eq('ok').sum()),n_converged=int(f.converged.eq(True).sum()),
            n_warnings=int(f.warnings.fillna('').str.len().gt(0).sum()),
            sigma2_min=f.sigma2.min(),sigma2_median=f.sigma2.median(),sigma2_max=f.sigma2.max(),
            niter_max=f.niter.max(),niter_median=f.niter.median(),n_cs=f.n_cs.sum()))
        ids = set(f.locus_id) & common
        s = ss[ss.locus_id.isin(ids)].sort_values('PIP',ascending=False).drop_duplicates('SNP').copy()
        c = cc[cc.locus_id.isin(ids)].copy()
        c['members'] = c.SNPs.str.split(',').map(set)
        recovered = set().union(*c.members.tolist()) & truth if len(c) else set()
        cs_coverage = c.members.map(lambda x:bool(x & truth)).mean() if len(c) else np.nan
        high = s[s.PIP>=.9]
        true_high = set(high.SNP) & truth
        posmod = high.position % 6000000
        base = dict(fit_mode=mode,method=method,phenotype=pheno,n_matched_converged_loci=len(ids),
            n_high_pip=len(high),n_high_pip_causal=len(true_high),
            high_pip_precision=len(true_high)/len(high) if len(high) else np.nan,
            n_high_pip_effect=int((posmod<2000000).sum()),
            n_high_pip_null=int((posmod>=2000000).sum()),
            n_high_pip_central_null=int(((posmod>=3000000)&(posmod<5000000)).sum()),
            n_cs=len(c),n_cs_with_causal=int(c.members.map(lambda x:bool(x & truth)).sum()),
            cs_coverage=cs_coverage,median_cs_size=c.cs_size.median(),
            n_causal_analyzed=int(s.SNP.isin(truth).sum()),
            n_causal_cs_recovered=len(recovered),conditional_cs_recovery=len(recovered)/s.SNP.isin(truth).sum(),
            genomewide_cs_recovery=len(recovered)/len(truth))
        for universe in ['analyzed_snps','panel_with_untested_score_zero']:
            q = s if universe=='analyzed_snps' else bim[['SNP']].merge(s[['SNP','PIP']],on='SNP',how='left').fillna({'PIP':0})
            y, score = q.SNP.isin(truth).to_numpy(), q.PIP.to_numpy()
            precision,recall,threshold = precision_recall_curve(y,score)
            rows.append(dict(base,universe=universe,n_evaluated=len(q),n_causal_evaluated=int(y.sum()),
                random_precision_baseline=float(y.mean()),average_precision=average_precision_score(y,score),
                trapezoid_auprc=auc(recall,precision)))
            # Save an exact curve, including the terminal (precision=1, recall=0).
            curves.append(pd.DataFrame(dict(fit_mode=mode,method=method,phenotype=pheno,
                universe=universe,precision=precision,recall=recall,threshold=np.r_[threshold,np.nan])))
        # This is explicitly a reconstruction sensitivity, not original effect truth.
        ef = B/f'max_snp_pve_results/snp_effects_pve_{pheno}.tsv'
        if ef.exists():
            e = pd.read_csv(ef,sep='\t').sort_values(['ABS_BETA_STANDARDIZED','SNP'])
            assert len(e)==10000 and set(e.SNP)==truth
            e['effect_decile'] = np.repeat(np.arange(1,11),1000)
            for dec,g in e.groupby('effect_decile'):
                eligible=set(g.SNP) & set(s.SNP)
                strata.append(dict(fit_mode=mode,method=method,phenotype=pheno,effect_decile=int(dec),
                    effect_source='historical ridge reconstruction from breeding values; not original coefficients',
                    validation='sensitivity_only_pending_original_effect_file',
                    abs_beta_min=g.ABS_BETA_STANDARDIZED.min(),abs_beta_max=g.ABS_BETA_STANDARDIZED.max(),
                    n_causal_all=len(g),n_causal_analyzed=len(eligible),n_pip_ge_0_9=len(set(g.SNP)&true_high),
                    high_pip_rate_all=len(set(g.SNP)&true_high)/len(g),
                    n_cs_recovered=len(set(g.SNP)&recovered),cs_recovery_rate_all=len(set(g.SNP)&recovered)/len(g)))

M = pd.DataFrame(rows)
M.to_csv(OUT/'sensitivity_comparison_by_phenotype.tsv',sep='\t',index=False)
pd.DataFrame(diag).to_csv(OUT/'fit_diagnostics_comparison.tsv',sep='\t',index=False)
pd.DataFrame(strata).to_csv(OUT/'effect_strata_reconstruction_sensitivity.tsv',sep='\t',index=False)
pd.concat(curves,ignore_index=True).to_csv(OUT/'precision_recall_curves.tsv.gz',sep='\t',index=False)
numeric = [x for x in M.select_dtypes(include=np.number).columns]
summary = M.groupby(['fit_mode','method','universe'])[numeric].agg(['mean','min','max'])
summary.columns = ['_'.join(x) for x in summary.columns]
summary.reset_index().to_csv(OUT/'sensitivity_comparison_summary.tsv',sep='\t',index=False)
compact=[]
for (mode,method),q in M[M.universe=='analyzed_snps'].groupby(['fit_mode','method']):
    t=q[['n_high_pip','n_high_pip_causal','n_high_pip_null','n_high_pip_central_null','n_cs']].sum()
    compact.append(dict(fit_mode=mode,method=method,n_phenotypes=len(q),
        average_precision_mean=q.average_precision.mean(),trapezoid_auprc_mean=q.trapezoid_auprc.mean(),
        random_precision_baseline_mean=q.random_precision_baseline.mean(),
        high_pip_calls_pooled=int(t.n_high_pip),high_pip_causal_pooled=int(t.n_high_pip_causal),
        high_pip_precision_pooled=t.n_high_pip_causal/t.n_high_pip,
        high_pip_null_pooled=int(t.n_high_pip_null),high_pip_central_null_pooled=int(t.n_high_pip_central_null),
        n_cs_pooled=int(t.n_cs),cs_coverage_mean=q.cs_coverage.mean(),
        conditional_cs_recovery_mean=q.conditional_cs_recovery.mean(),
        median_cs_size_mean=q.median_cs_size.mean()))
pd.DataFrame(compact).to_csv(OUT/'sensitivity_comparison_table.tsv',sep='\t',index=False)
keys=['method','phenotype','universe']
pair=M[M.fit_mode=='ERV_FALSE'].drop(columns='fit_mode').merge(M[M.fit_mode=='ERV_TRUE'].drop(columns='fit_mode'),on=keys,suffixes=('_old','_new'),validate='one_to_one')
for k in numeric: pair[k+'_change']=pair[k+'_new']-pair[k+'_old']
pair.to_csv(OUT/'sensitivity_paired_changes.tsv',sep='\t',index=False)
(OUT/'metric_definitions.txt').write_text(
    'Paired comparison uses the same loci converged under both residual-variance settings.\n'
    'Each SNP occurs once per GWAS method/phenotype, retaining its largest PIP across overlapping loci.\n'
    'High PIP means PIP >=0.9; precision denominator is all high-PIP unique SNP calls.\n'
    'Effect: position modulo 6 Mb <2 Mb; null: remainder >=2 Mb; central null: [3,5) Mb.\n'
    'Central null is a subset of null. Genome classification uses SNP position, not lead locus.\n'
    'Average precision is the noninterpolated, stepwise area; trapezoid AUPRC is supplied separately.\n'
    'The conditional universe is all analyzed unique SNPs. The full-panel universe assigns score zero\n'
    'to SNPs not fine-mapped; this is a ranking convention, not a fitted posterior probability.\n'
    'Random precision baseline is causal SNP prevalence in the stated universe.\n'
    'CS coverage is sets containing >=1 causal SNP / all reported CSs. Conditional recovery uses\n'
    'unique analyzed causal SNPs; genome-wide recovery uses the complete 10,000 causal set.\n'
    'Explicit counts: CS coverage = n_cs_with_causal/n_cs; conditional recovery =\n'
    'n_causal_cs_recovered/n_causal_analyzed; random precision = n_causal_evaluated/n_evaluated.\n'
    'Mean summaries weight each of the ten phenotypes equally; counts are not pooled unless explicit.\n'
    'Effect-size deciles are sensitivity summaries from historical ridge-reconstructed coefficients.\n'
    'They must not be represented as verified original simulation effects.\n')
print(pd.DataFrame(diag).groupby(['fit_mode','method'])[['n_attempted','n_converged','n_warnings','n_cs']].sum().to_string())
print(M[M.universe=='analyzed_snps'].groupby(['fit_mode','method'])[['average_precision','n_high_pip','high_pip_precision','cs_coverage','conditional_cs_recovery']].mean().to_string())
(W/'qc/EVALUATION.COMPLETE').write_text('All manifest loci validated and paired; see metric_definitions.txt\n')
