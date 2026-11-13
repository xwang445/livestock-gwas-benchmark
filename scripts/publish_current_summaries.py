#!/usr/bin/env python3
"""Publish current root aggregate views from committed, dated scientific sources.

No restricted data are opened. The historical snapshot is immutable. Legacy
column names are retained, with explicit strict-count and precision fields.
"""
import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ALPHA = 5e-8
ALIASES = {'REGENIE_WO_LOCO': 'REGENIE_woLOCO'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_table(path):
    # Keep meaningful literal categories such as "null"; preserve audited floats.
    return pd.read_csv(path, sep='\t', keep_default_na=False, na_values=['NA', ''],
                       float_precision='round_trip')


def serialize(frame):
    buffer = io.StringIO()
    frame.to_csv(buffer, sep='\t', index=False, na_rep='NA', float_format='%.17g', lineterminator='\n')
    return buffer.getvalue().encode('utf-8')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def aggregate(frame, groups, values):
    specs = {}
    for name, field in values.items():
        for prefix, statistic in [('mean', 'mean'), ('sd', 'std'), ('min', 'min'), ('max', 'max')]:
            specs[prefix + '_' + name] = (field, statistic)
    specs['n_phenotypes'] = ('phenotype', 'nunique')
    return frame.groupby(groups, dropna=False, sort=True).agg(**specs).reset_index()


def build_tables(repo):
    """Return current views and logical source fingerprints without writing files."""
    repo = Path(repo)
    dated = repo / 'derived_data/revision_2026-11-05/association'
    historical = repo / 'derived_data/historical_2026-11-04'
    fingerprints = {}

    def load(path):
        fingerprints[str(path.relative_to(repo))] = digest(path)
        return read_table(path)

    registry = load(repo / 'config/revision_methods.tsv')
    configurations = list(registry.configuration)
    methods = [ALIASES.get(name, name) for name in configurations]
    require(len(methods) == 10 and len(set(methods)) == 10, 'Expected ten method configurations.')

    def current(frame):
        result = frame[frame.configuration.isin(configurations)].copy()
        result['method'] = result.configuration.replace(ALIASES)
        return result

    by = current(load(dated / 'association_by_phenotype.tsv'))
    require(len(by) == 100 and not by.duplicated(['method', 'phenotype']).any(), 'Expected 100 association rows.')
    require(by.groupby('method').phenotype.nunique().eq(10).all(), 'Incomplete method phenotypes.')
    qc = current(load(dated / 'association_qc.tsv'))
    require(len(qc) == 100 and not qc.duplicated(['method', 'phenotype']).any(), 'Expected 100 QC rows.')
    require(qc.rows.eq(135358).all() and qc.missing_p.eq(0).all(), 'Unexpected panel or missing P values.')
    by = by.merge(qc[['method', 'phenotype', 'rows', 'missing_p', 'zero_p']],
                  on=['method', 'phenotype'], validate='one_to_one')
    by = by.sort_values(['method', 'phenotype']).reset_index(drop=True)
    result = {}

    counts = by[['method', 'phenotype']].copy()
    for target, source in [('n_snps', 'rows'), ('n_significant_snps', 'significant_snps'),
                           ('n_significant_loci', 'loci'), ('inclusive_snps', 'significant_snps'),
                           ('inclusive_loci', 'loci'), ('at_threshold', 'snps_at_threshold'),
                           ('missing_p', 'missing_p'), ('zero_p', 'zero_p'),
                           ('strict_significant_snps', 'strict_significant_snps'), ('strict_loci', 'strict_loci')]:
        counts[target] = by[source]
    counts['significance_comparator'] = 'le'
    require((counts.n_significant_snps - counts.strict_significant_snps).eq(counts.at_threshold).all(),
            'Inclusive/strict count reconciliation failed.')
    result['method_summary/significant_snps_by_pheno.tsv'] = counts
    result['method_summary/significant_loci_by_pheno.tsv'] = counts[
        ['method', 'phenotype', 'n_significant_loci', 'inclusive_loci', 'strict_loci', 'significance_comparator']]
    result['method_summary/method_summary.tsv'] = aggregate(
        counts, ['method'], {'sig_snps': 'n_significant_snps', 'sig_loci': 'n_significant_loci'})

    spill = by[['method', 'phenotype']].copy()
    for target, source in [('n_significant_loci', 'loci'), ('n_spillover_loci', 'central_null_loci'),
                           ('n_central_windows', 'central_null_windows'), ('spillover_rate', 'central_null_spillover_rate')]:
        spill[target] = by[source]
    require(np.allclose(spill.spillover_rate, spill.n_spillover_loci / spill.n_central_windows),
            'Spillover denominator mismatch.')
    result['spillover_summary/spillover_by_pheno.tsv'] = spill
    result['spillover_summary/spillover_summary.tsv'] = aggregate(
        spill, ['method'], {'spillover_rate': 'spillover_rate', 'spillover_loci': 'n_spillover_loci'})

    lam = by[['method', 'phenotype']].copy()
    lam['lambda_gc'] = by.lambda_gc_central_null_historical
    lam['lambda_gc_p_based'] = by.lambda_gc_central_null
    lam['statistic_definition'] = np.where(
        by.method.str.startswith(('FarmCPU_', 'BLINK_')),
        'historical_qchisq_1_minus_p_finite_only', 'historical_beta_se_squared_finite_only')
    lam['n_central_snps'] = by.n_central_null_snps
    lam['n_historical_finite_statistics'] = by.n_historical_finite_central_statistics
    # Canonical central masks exclude nonfinite P; QC above verifies none missing.
    lam['n_finite_p'] = by.n_central_null_snps
    result['lambda_gc_summary/lambda_gc_by_pheno.tsv'] = lam
    result['lambda_gc_summary/lambda_gc_summary.tsv'] = aggregate(lam, ['method'], {'lambda_gc': 'lambda_gc'})

    rates = current(load(dated / 'power_and_null_snp_rate_by_phenotype.tsv'))
    require(np.array_equal(rates.alpha.to_numpy(), rates.alpha_hex.map(float.fromhex).to_numpy()),
            'Power-grid alpha decimal/hex mismatch.')
    precision = ['minus_log10_alpha', 'alpha_hex', 'alpha_expression']
    power = rates[rates.region.eq('causal_positions')][
        ['method', 'phenotype', 'alpha', 'n_tested', 'n_significant', 'rate'] + precision].rename(
            columns={'n_tested': 'n_tested_causal', 'n_significant': 'n_significant_causal', 'rate': 'power'})
    power = power.sort_values(['method', 'phenotype', 'alpha']).reset_index(drop=True)
    result['power_summary/power_by_pheno.tsv'] = power
    result['power_summary/power_summary.tsv'] = aggregate(power, ['method', 'alpha'] + precision, {'power': 'power'})
    fpr = rates[rates.region.isin(['central_null', 'all_null'])][
        ['method', 'phenotype', 'region', 'alpha', 'n_tested', 'n_significant', 'rate'] + precision].rename(columns={'rate': 'fpr'})
    fpr['source_region'] = fpr.region
    fpr['region'] = fpr.region.replace({'all_null': 'outside_effect'})
    fpr = fpr.sort_values(['method', 'phenotype', 'region', 'alpha']).reset_index(drop=True)
    result['power_summary/fpr_by_pheno.tsv'] = fpr
    result['power_summary/fpr_summary.tsv'] = aggregate(fpr, ['method', 'region', 'alpha'] + precision, {'fpr': 'fpr'})

    # Historical public denominators encode the same complete marker panel and
    # closed-endpoint 50-kb geometry. They are invariant across all 90 old jobs.
    old_spatial = load(historical / 'power_summary/spatial_enrichment_by_pheno.tsv')
    geometry = old_spatial.groupby('distance_index')[['n_tested', 'n_nonempty_windows']]
    require(geometry.nunique().eq(1).all().all(), 'Historical spatial denominators are not invariant.')
    geometry = geometry.first().reset_index()
    spatial = current(load(dated / 'spatial_by_phenotype.tsv')).merge(geometry, on='distance_index', validate='many_to_one')
    require(len(spatial) == 4000 and spatial.nonempty_bins.eq(spatial.n_nonempty_windows).all(),
            'Current/historical spatial geometry mismatch.')
    reconstructed = spatial.observed_expected_pooled_snps * spatial.n_tested * ALPHA
    require(np.max(np.abs(reconstructed - np.rint(reconstructed))) < 1e-8,
            'Spatial pooled ratios do not reconstruct integer significant counts.')
    spatial['n_significant'] = np.rint(reconstructed).astype(np.int64)
    spatial = spatial.rename(columns={'observed_expected_equal_nonempty_bins': 'observed_expected',
                                      'observed_expected_pooled_snps': 'pooled_observed_expected'})
    spatial = spatial[['distance_index', 'observed_expected', 'n_nonempty_windows', 'n_tested',
                       'n_significant', 'pooled_observed_expected', 'method', 'phenotype']].sort_values(
                           ['method', 'phenotype', 'distance_index']).reset_index(drop=True)
    result['power_summary/spatial_enrichment_by_pheno.tsv'] = spatial
    result['power_summary/spatial_enrichment_summary.tsv'] = aggregate(
        spatial, ['method', 'distance_index'], {'observed_expected': 'observed_expected'})

    overlap = load(dated / 'pairwise_locus_overlap.tsv').rename(columns={
        'configuration_a': 'method_a', 'configuration_b': 'method_b', 'n_a': 'n_loci_a',
        'n_b': 'n_loci_b', 'intersection': 'n_intersection', 'union': 'n_union'})
    overlap['method_a'] = overlap.method_a.replace(ALIASES)
    overlap['method_b'] = overlap.method_b.replace(ALIASES)
    overlap = overlap[['phenotype', 'method_a', 'method_b', 'n_loci_a', 'n_loci_b',
                       'n_intersection', 'n_union', 'jaccard', 'jaccard_set_theoretic', 'historical_jaccard_defined']]
    require(len(overlap) == 1000, 'Expected 1,000 method-pair/phenotype overlap cells.')
    result['overlap_summary/locus_overlap_by_pheno.tsv'] = overlap.sort_values(
        ['phenotype', 'method_a', 'method_b']).reset_index(drop=True)
    overlap_summary = aggregate(overlap, ['method_a', 'method_b'], {'jaccard': 'jaccard'})
    defined = overlap.groupby(['method_a', 'method_b']).jaccard.count().reset_index(name='n_defined_phenotypes')
    result['overlap_summary/locus_overlap_summary.tsv'] = overlap_summary.merge(
        defined, on=['method_a', 'method_b'], validate='one_to_one')
    crossrep = current(load(dated / 'cross_replicate_stability.tsv'))
    require(len(crossrep) == 1000 and not crossrep.duplicated(['method', 'phenotype_a', 'phenotype_b']).any(),
            'Expected 1,000 within-method phenotype-pair overlap cells.')
    cross_fields = ['method', 'phenotype_a', 'phenotype_b', 'n_intersection', 'n_union', 'jaccard']
    cross_fields += [name for name in ['jaccard_set_theoretic'] if name in crossrep]
    result['overlap_summary/cross_replicate_stability.tsv'] = crossrep[cross_fields].sort_values(
        ['method', 'phenotype_a', 'phenotype_b']).reset_index(drop=True)

    old_h2 = load(historical / 'method_summary/heritability_from_logs.tsv')
    h2_fast = old_h2[old_h2.method.eq('fastGWA')].copy()
    h2_fast['source_path'], h2_fast['source_kind'] = h2_fast.log, 'reported_fastGWA_log'
    h2_fast['precision_note'] = 'Reported log estimate retained at source precision.'
    h2_slemm = current(load(dated / 'slemm_heritability_by_phenotype.tsv'))
    h2_slemm = h2_slemm[h2_slemm.method.eq('SLEMM')][['method', 'phenotype', 'h2', 'sha256']].copy()
    h2_slemm['source_path'] = h2_slemm.phenotype.map(
        lambda ph: 'SLEMM_qf90_window3000/spg.even.100k.pheno%d.reml.vc.csv' % ph)
    h2_slemm['source_kind'] = 'reported_SLEMM_variance_component_CSV'
    h2_slemm['precision_note'] = 'Reported current qf90/window3000 CSV estimate; not the historical log source.'
    h2_gcta = current(load(dated / 'gcta_heritability_by_phenotype.tsv'))
    require(h2_gcta.status.eq('completed_log_rounded_estimate').all(), 'GCTA heritability requires completed MLMA logs.')
    h2_gcta = h2_gcta.rename(columns={'h2_from_printed_variances': 'h2', 'log_sha256': 'sha256'})
    h2_gcta['source_path'] = h2_gcta.phenotype.map(lambda ph: 'gcta/spg.even.100k.pheno%d.log' % ph)
    h2_gcta['log'], h2_gcta['source_kind'] = h2_gcta.source_path, 'ratio_of_rounded_GCTA_REML_variances'
    h2_fields = ['method', 'phenotype', 'h2', 'log', 'line', 'sha256', 'source_path', 'source_kind',
                 'precision_note', 'Vg_printed', 'Ve_printed']
    h2 = pd.concat([h2_fast, h2_slemm, h2_gcta], ignore_index=True).reindex(columns=h2_fields)
    require(len(h2) == 30 and h2.h2.notna().all(), 'Expected complete fastGWA, SLEMM and GCTA heritability.')
    result['method_summary/heritability_from_logs.tsv'] = h2.sort_values(['method', 'phenotype']).reset_index(drop=True)

    inputs = current(load(dated / 'current_input_manifest.tsv'))
    require(len(inputs) == 100 and not inputs.duplicated(['method', 'phenotype']).any(), 'Expected 100 source manifest rows.')
    require(inputs.rows.eq(135358).all() and inputs.sha256.astype(str).str.fullmatch('[0-9a-f]{64}').all(),
            'Invalid source manifest fingerprints.')
    require(inputs.path.map(lambda name: not Path(name).is_absolute() and '..' not in Path(name).parts).all(),
            'Source paths must be logical relative filenames.')
    result['input_manifest.tsv'] = inputs.drop(columns='configuration').sort_values(['method', 'phenotype']).reset_index(drop=True)

    auxiliary = repo / 'scripts/current_auxiliary_summaries.py'
    if auxiliary.exists():
        spec = importlib.util.spec_from_file_location('current_auxiliary_summaries', auxiliary)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        result.update(module.build_auxiliary_tables(repo))
        if hasattr(module, 'build_auxiliary_sources'):
            auxiliary_sources = module.build_auxiliary_sources(repo)
            if isinstance(auxiliary_sources, dict):
                fingerprints.update(auxiliary_sources)
            else:
                for path in auxiliary_sources:
                    path = Path(path)
                    fingerprints[str(path.relative_to(repo))] = digest(path)
        fingerprints['scripts/current_auxiliary_summaries.py'] = digest(auxiliary)
    fingerprints['scripts/publish_current_summaries.py'] = digest(repo / 'scripts/publish_current_summaries.py')
    return result, fingerprints


def rendered_outputs(repo):
    tables, sources = build_tables(repo)
    outputs = {name: serialize(frame) for name, frame in tables.items()}
    auxiliary = Path(repo) / 'scripts/current_auxiliary_summaries.py'
    if auxiliary.exists():
        spec = importlib.util.spec_from_file_location('current_auxiliary_metadata', auxiliary)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        if hasattr(module, 'build_auxiliary_metadata'):
            for name, metadata in module.build_auxiliary_metadata(repo).items():
                outputs[name] = (json.dumps(metadata, indent=2, sort_keys=True) + '\n').encode('utf-8')
    manifest = dict(status='CURRENT_INCLUSIVE_AGGREGATES', n_jobs=100, n_configurations=10,
                    n_phenotypes=10, expected_snps=135358, central_windows=1480,
                    significance='P <= 5e-8 on the source reporting scale',
                    current_slemm='qf90 window3000', gcta='all ten MLMA analyses complete',
                    root_aliases=ALIASES, historical_snapshot='derived_data/historical_2026-11-04',
                    lambda_definition='historical statistic estimator applied to current fits; stable P-based estimator alongside',
                    spatial_denominators='Invariant historical public closed-endpoint geometry; current pooled ratios reconstruct integer counts',
                    spatial_integer_reconstruction_tolerance=1e-8,
                    gcta_heritability_precision='Ratio of final variances printed to five decimal places; no SE inferred',
                    code_sha256=sources['scripts/publish_current_summaries.py'],
                    config_sha256=sources['config/revision_methods.tsv'],
                    validation_status='CURRENT_COUNTS_AND_DATED_AGGREGATE_SOURCES_VERIFIED',
                    source_sha256=sources,
                    output_sha256={name: hashlib.sha256(data).hexdigest() for name, data in sorted(outputs.items())})
    outputs['run_manifest.json'] = (json.dumps(manifest, indent=2, sort_keys=True) + '\n').encode('utf-8')
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=ROOT)
    parser.add_argument('--check', action='store_true', help='Verify current committed root views without writing files')
    args = parser.parse_args()
    outputs = rendered_outputs(args.repo.resolve())
    failures = []
    for name, data in sorted(outputs.items()):
        path = args.repo / 'derived_data' / name
        if args.check:
            if not path.is_file() or path.read_bytes() != data:
                failures.append(name)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
    if failures:
        parser.exit(1, 'Current summary views differ from their dated sources: ' + ', '.join(failures) + '\n')
    print(('Checked' if args.check else 'Published') + ' {} current aggregate files; all ten configurations complete.'.format(len(outputs)))


if __name__ == '__main__':
    main()
