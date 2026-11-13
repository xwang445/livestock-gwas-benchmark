"""Deterministic public current views for fine-mapping, cached theory and LD.

Consumed by publish_current_summaries.py. All inputs are public aggregate tables.
The historical November 4 views remain in their separate dated archive.
"""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

INPUTS = [
    'revision_2026-11-05/susie/maf_matched_counts_by_phenotype.tsv',
    'revision_2026-11-05/susie/maf_matched_pool_qc.tsv',
    'revision_2026-11-05/susie/pooled_finemapping_metrics.tsv',
    'submission_2026-11-12/figure07_cached_theory_bin_summary.tsv',
    'submission_2026-11-12/figure07_cached_observed_bins_by_phenotype.tsv',
    'submission_2026-11-12/figure01_ld_10kb_bins.tsv',
    'submission_2026-11-12/figure01_shared_spline.tsv',
]


def read(repo, relative):
    return pd.read_csv(Path(repo) / 'derived_data' / relative, sep='\t', float_precision='round_trip')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def build_auxiliary_sources(repo):
    return [Path(repo) / 'derived_data' / relative for relative in INPUTS]


def fingerprints(repo):
    return {path.relative_to(repo).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in build_auxiliary_sources(repo)}


def build_auxiliary_tables(repo):
    counts = read(repo, INPUTS[0])
    pools = read(repo, INPUTS[1])
    pooled = read(repo, INPUTS[2])
    pooled = pooled.loc[pooled.fit_mode.eq('ERV_TRUE')].copy()
    keys = ['fine_mapping', 'GWAS_method', 'phenotype']
    require(len(counts) == 120 and len(pools) == 40, 'Expected four pipeline/method pairs across ten phenotypes.')
    require(not pools.duplicated(keys).any() and not counts.duplicated(keys + ['class']).any(), 'Duplicate MAF aggregate key.')
    details = []
    for key, group in counts.groupby(keys, sort=True):
        pipeline, method, phenotype = key
        observed = group.set_index('class').observed
        expected = group.set_index('class').expected
        pool = pools.loc[(pools.fine_mapping == pipeline) & (pools.GWAS_method == method) & (pools.phenotype == phenotype)]
        require(len(pool) == 1, 'Missing MAF pool QC row.')
        pool = pool.iloc[0]
        require(set(observed.index) == {'Causal', 'Non-causal in effect windows', 'Null windows'}, 'Unexpected MAF classes.')
        require(int(observed.sum()) == int(pool.n_calls), 'MAF classes do not partition high-PIP calls.')
        require(group.pool_size.eq(pool.pool_size).all(), 'MAF pool-size mismatch.')
        details.append(dict(method=method.replace('REGENIE-LOCO', 'REGENIE_LOCO'), pipeline=pipeline,
                            phenotype=int(str(phenotype).replace('pheno', '')), fit_mode='ERV_TRUE',
                            n_high_pip=int(pool.n_calls), n_causal=int(observed['Causal']),
                            n_null=int(observed['Null windows']),
                            n_effect=int(observed['Causal'] + observed['Non-causal in effect windows']),
                            expected_causal=float(expected['Causal']), n_pool=int(pool.pool_size),
                            n_maf_bins=int(pool.n_bins), n_converged_loci=int(pool.n_loci)))
    details = pd.DataFrame(details).sort_values(['method', 'pipeline', 'phenotype']).reset_index(drop=True)
    summaries = []
    for (method, pipeline), group in details.groupby(['method', 'pipeline'], sort=True):
        require(set(group.phenotype) == set(range(1, 11)), 'Incomplete final fine-mapping phenotypes.')
        row = dict(method=method, pipeline=pipeline, fit_mode='ERV_TRUE')
        for column in ['n_high_pip', 'n_causal', 'n_null', 'n_effect', 'expected_causal']:
            row[column] = group[column].sum()
        row.update(null_percent=100 * row['n_null'] / row['n_high_pip'],
                   causal_percent=100 * row['n_causal'] / row['n_high_pip'],
                   noncausal_percent=100 * (row['n_high_pip'] - row['n_causal']) / row['n_high_pip'],
                   maf_expected_percent=100 * row['expected_causal'] / row['n_high_pip'],
                   enrichment=row['n_causal'] / row['expected_causal'])
        canonical = pooled.loc[(pooled.pipeline == pipeline) &
                               (pooled.method == method.replace('REGENIE_LOCO', 'REGENIE-LOCO'))]
        require(len(canonical) == 1, 'Missing canonical final pooled metrics.')
        for output, source in [('n_high_pip', 'n_high_pip'), ('n_causal', 'n_high_pip_causal'),
                               ('n_null', 'n_high_pip_null'), ('n_effect', 'n_high_pip_effect')]:
            require(row[output] == canonical.iloc[0][source], 'Final fine-mapping recount mismatch: ' + output)
        summaries.append(row)
    summary = pd.DataFrame(summaries)
    null_effect = summary[['method', 'pipeline', 'fit_mode', 'n_high_pip', 'n_null', 'n_effect', 'null_percent']].copy()
    converged = details.groupby(['method', 'pipeline', 'fit_mode'], as_index=False).n_converged_loci.sum()

    theory = read(repo, INPUTS[3])
    observed = read(repo, INPUTS[4])
    require(len(theory) == 40 and len(observed) == 400, 'Current Figure 7 must contain two methods and 20 bins.')
    require(set(theory.method) == {'BOLT-LMM', 'REGENIE_LOCO'} and
            set(theory.region) == {'Genome-wide'}, 'Unsupported current theoretical scope.')
    check = observed.groupby(['method', 'region', 'bin']).mean_observed_chisq.agg(['mean', 'std']).reset_index()
    check = theory.merge(check, on=['method', 'region', 'bin'], validate='one_to_one')
    require(np.allclose(check.mean_observed_chisq, check['mean'], rtol=1e-12, atol=1e-12) and
            np.allclose(check.sd_observed_chisq, check['std'], rtol=1e-12, atol=1e-12),
            'Retained Figure 7 aggregate observations disagree.')
    require(np.allclose(theory.predicted_chisq, 1 + theory.mean_predicted_ncp, rtol=1e-12, atol=1e-12),
            'Retained bin means violate expected chi-square = 1 + NCP.')
    cached = pd.DataFrame([dict(status='RETAINED_FIGURE7_BIN_CACHE', n_methods=2,
                                region='Genome-wide', n_bins_per_method=20, n_phenotypes=10,
                                n_summary_rows=40, n_observation_rows=400,
                                n_theory_snps=int(theory.loc[theory.method == 'BOLT-LMM', 'n_obs'].sum()),
                                max_expectation_error=float(np.max(np.abs(theory.predicted_chisq - 1 - theory.mean_predicted_ncp))),
                                upstream_equation12_validation='NOT_RUN', variant_level_cache_available=False)])
    bins = read(repo, INPUTS[5])
    spline = read(repo, INPUTS[6])
    require(len(bins) == 1000 and int(bins.n_pairs.sum()) == 78082232, 'Unexpected independently verified LD bins.')
    require(len(spline) >= len(bins) and np.isfinite(spline[['mid_mb', 'mean_r2']].to_numpy()).all()
            and (np.diff(spline.mid_mb) >= 0).all(), 'Invalid shared LD spline.')
    require(np.isclose(spline.mid_mb.iloc[0], bins.mid_mb.iloc[0]) and
            np.isclose(spline.mid_mb.iloc[-1], bins.mid_mb.iloc[-1]), 'Shared spline does not cover the LD bins.')
    return {
        'finemapping_summary/high_pip_by_pheno.tsv': details,
        'finemapping_summary/high_pip_summary.tsv': summary,
        'finemapping_summary/null_vs_effect_summary.tsv': null_effect,
        'finemapping_summary/converged_loci.tsv': converged,
        'equation13_summary/equation13_summary.tsv': theory,
        'equation13_summary/equation13_by_pheno.tsv': observed,
        'equation13_summary/cached_theory_validation.tsv': cached,
        'ld_summary/ld_decay_10kb.tsv': bins[['bin_id', 'bin_start', 'bin_end', 'mid_mb', 'n_pairs', 'mean_r2']],
        'ld_summary/ld_decay_fitted_mean.tsv': spline[['mid_mb', 'mean_r2']],
    }


def build_auxiliary_metadata(repo):
    return {'finemapping_summary/verification_manifest.json': dict(
        status='CURRENT_FINAL_TRUE_AGGREGATES_REFRESHED', fit_mode='ERV_TRUE',
        methods=['BOLT-LMM', 'REGENIE_LOCO'], pipelines=['SuSiE-adj', 'SuSiE-RSS'],
        n_phenotypes=10, detailed_outputs_written=False,
        input_sha256={key: value for key, value in fingerprints(repo).items() if '/susie/' in key},
        limitations=['No new fitting; final saved fits and exact MAF-matched expectations were evaluated.',
                     'No GCTA fine-mapping run is available; its GWAS summaries are separate.',
                     'Historical RSS FALSE results remain in the historical archive and labeled revision sensitivity tables.'])}
