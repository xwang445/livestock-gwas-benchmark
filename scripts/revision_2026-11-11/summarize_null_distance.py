#!/usr/bin/env python3
"""Portable c66.py estimator using an explicitly supplied restricted framework.

The original framework's dist_to_boundary column is required. Distances are never
reconstructed from genomic periodicity or inferred from aggregate tables.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import fingerprint, require, sha256


def summarize(data, pipeline, mode):
    rows = []
    for method, group in data.groupby('GWAS_method'):
        calls = group.loc[group.PIP >= 0.9]
        n_calls, n_pool = len(calls), len(group)
        row = dict(pipeline=pipeline, fit_mode=mode, method=method,
                   n_high_pip=n_calls, n_evaluated=n_pool)
        for label, column in [('null', 'is_null_window'), ('central_null', 'is_central_null')]:
            n_call = int(calls[column].eq(1).sum())
            n_evaluated = int(group[column].eq(1).sum())
            row['n_high_pip_' + label] = n_call
            row['high_pip_' + label + '_fraction'] = n_call / n_calls if n_calls else np.nan
            row['n_evaluated_' + label] = n_evaluated
            row['evaluated_' + label + '_fraction'] = n_evaluated / n_pool
        for distance in [200000, 500000, 1000000]:
            key = 'null_ge_{}kb'.format(distance // 1000)
            n_call = int((calls.is_null_window.eq(1) & calls.dist_to_boundary.ge(distance)).sum())
            n_evaluated = int((group.is_null_window.eq(1) & group.dist_to_boundary.ge(distance)).sum())
            row['n_high_pip_' + key] = n_call
            row['high_pip_' + key + '_fraction'] = n_call / n_calls if n_calls else np.nan
            row['n_evaluated_' + key] = n_evaluated
            row['evaluated_' + key + '_fraction'] = n_evaluated / n_pool
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ['framework', 'adjusted', 'rss', 'rss-historical', 'output']:
        parser.add_argument('--' + option, type=Path, required=True)
    args = parser.parse_args()
    columns = ['snp', 'is_effect_window', 'is_null_window', 'is_central_null', 'dist_to_boundary']
    framework = pd.read_csv(args.framework, sep='\t', usecols=columns)
    require(framework.snp.notna().all() and framework.snp.is_unique, 'Framework SNPs must be unique.')
    require(np.isfinite(framework.dist_to_boundary).all(), 'Framework distances must be finite.')
    for column in columns[1:4]:
        require(framework[column].isin([0, 1]).all(), 'Framework spatial flags must be binary.')
    require((framework.is_effect_window + framework.is_null_window == 1).all(),
            'Effect and null flags must partition the SNP pool.')
    require((framework.is_central_null <= framework.is_null_window).all(), 'Central-null must be a null subset.')
    sources = [('SuSiE-adj', 'ERV_TRUE', args.adjusted), ('SuSiE-RSS', 'ERV_TRUE', args.rss),
               ('SuSiE-RSS', 'ERV_FALSE', args.rss_historical)]
    tables = []
    for pipeline, mode, path in sources:
        data = pd.read_csv(path, sep='\t', usecols=['phenotype', 'GWAS_method', 'SNP', 'PIP'])
        require(data[['phenotype', 'GWAS_method', 'SNP']].notna().all().all(), 'Missing PIP key.')
        require(np.isfinite(data.PIP).all() and data.PIP.between(0, 1).all(), 'Invalid PIP.')
        data = data.groupby(['GWAS_method', 'phenotype', 'SNP'], as_index=False, sort=False).PIP.max()
        data = data.merge(framework, left_on='SNP', right_on='snp', how='left', validate='many_to_one')
        require(data.snp.notna().all(), 'Every evaluated SNP requires a supplied framework row.')
        tables.append(summarize(data, pipeline, mode))
    args.output.mkdir(parents=True, exist_ok=True)
    pd.concat(tables, ignore_index=True).to_csv(args.output / 'null_distance_summary.tsv',
                                               sep='\t', index=False, float_format='%.15g')
    manifest = dict(status='complete', estimator='Supplied c66.py with max-PIP de-duplication',
                    distance_source='Explicit supplied framework dist_to_boundary; no reconstructed distances',
                    cutoff_comparison='PIP >= 0.9; distance >= 200000, 500000, 1000000 bp',
                    denominators='All high-PIP calls for call fractions; all evaluated pairs for pool fractions',
                    detailed_outputs_written=False,
                    inputs=[fingerprint('spatial_framework', args.framework)] +
                           [fingerprint(pipeline + '_' + mode, path) for pipeline, mode, path in sources],
                    source_sha256=sha256(__file__))
    (args.output / 'null_distance_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
