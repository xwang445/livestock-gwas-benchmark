#!/usr/bin/env python3
"""Recount fine-mapping claims and recovered-PVE groups from restricted PIPs.

This ports the supplied powbin_rep.py estimator to final and historical fits.
Only method/phenotype/group aggregates and path-free fingerprints are written.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import causal_ids, effect_long, fingerprint, require, sha256

LABELS = ['<0.5×', '0.5–1×', '1–2×', '2–4×', '≥4×']
T9 = 2.262157  # Retain the supervisor's rounded t_(9, 0.975).


def pve_group_summaries(data, pipeline, mode, n_phenotypes=10):
    rows, details = [], []
    for (method, label), group in data.groupby(['GWAS_method', 'group']):
        per = group.groupby('phenotype').call.agg(['size', 'sum', 'mean'])
        require(len(per) == n_phenotypes, 'Every PVE group must cover all phenotypes.')
        per = per.reset_index().rename(columns={'size': 'n', 'sum': 'calls', 'mean': 'power'})
        per['pipeline'], per['fit_mode'], per['method'], per['group'] = pipeline, mode, method, label
        details.append(per)
        mean, sd = per.power.mean(), per.power.std(ddof=1)
        rows.append(dict(method=method, finemap=pipeline, fit_mode=mode, group=label,
                         n=len(group), calls=int(group.call.sum()), pooled=group.call.mean(),
                         power=mean, lo=mean - T9 * sd / np.sqrt(n_phenotypes),
                         hi=mean + T9 * sd / np.sqrt(n_phenotypes), min_n=int(per.n.min()),
                         max_recovered_pve=group.recovered_pve.max()))
    return pd.DataFrame(rows), pd.concat(details, ignore_index=True)


def prepare_pips(path, pve, causal):
    data = pd.read_csv(path, sep='\t', usecols=['phenotype', 'GWAS_method', 'SNP', 'PIP'])
    require(data[['phenotype', 'GWAS_method', 'SNP']].notna().all().all(), 'Missing PIP key.')
    require(np.isfinite(data.PIP).all() and data.PIP.between(0, 1).all(), 'PIP must lie in [0,1].')
    data = data.groupby(['GWAS_method', 'phenotype', 'SNP'], as_index=False, sort=False).PIP.max()
    data = data.merge(pve, on=['SNP', 'phenotype'], how='left', validate='many_to_one')
    data['causal'] = data.SNP.isin(causal)
    require(data.loc[data.causal, 'recovered_pve'].notna().all(), 'Missing recovered PVE for a causal row.')
    require(data.loc[~data.causal, 'recovered_pve'].isna().all(), 'Noncausal row has recovered PVE.')
    bins = pd.cut(data.recovered_pve / (0.5 / len(causal)), [0, 0.5, 1, 2, 4, np.inf],
                  right=False, labels=LABELS).astype(str)
    data['group'] = np.where(data.causal, bins, 'Non-causal')
    require(set(data.group) == set(LABELS + ['Non-causal']), 'Unexpected PVE group inventory.')
    data['call'] = data.PIP >= 0.9
    # Benchmark IDs encode chromosome_position_alleles. Preserve verified design boundaries.
    position = data.SNP.str.split('_', expand=True)[1].astype(np.int64)
    cycle = position % 6000000
    data['null'] = cycle >= 2000000
    data['central_null'] = (cycle >= 3000000) & (cycle < 5000000)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ['effects', 'causal', 'adjusted', 'rss', 'rss-historical', 'pooled-metrics', 'output']:
        parser.add_argument('--' + option, type=Path, required=True)
    parser.add_argument('--supervisor-groups', type=Path,
                        help='Optional shared powbin_rep.tsv; reproduce its baseline before updating it.')
    args = parser.parse_args()
    causal = causal_ids(args.causal)
    require(len(causal) == 10000, 'Expected the fixed 10,000 causal SNPs.')
    long = effect_long(pd.read_csv(args.effects, sep='\t', float_precision='round_trip'), causal)
    pve = long.rename(columns={'recovered_pve_fraction': 'recovered_pve'})
    pve.phenotype = pve.phenotype.map(lambda k: 'pheno{}'.format(k))
    metrics = pd.read_csv(args.pooled_metrics, sep='\t')
    sources = [('SuSiE-adj', 'ERV_TRUE', args.adjusted), ('SuSiE-RSS', 'ERV_TRUE', args.rss),
               ('SuSiE-RSS', 'ERV_FALSE', args.rss_historical)]
    summaries, details, pooled = [], [], []
    for pipeline, mode, path in sources:
        data = prepare_pips(path, pve, causal)
        require(set(data.phenotype) == {'pheno{}'.format(k) for k in range(1, 11)},
                'Expected ten named phenotype replicates.')
        for method, group in data.groupby('GWAS_method'):
            calls = group.loc[group.call]
            record = dict(pipeline=pipeline, fit_mode=mode, method=method,
                          n_evaluated=len(group), n_causal_evaluated=int(group.causal.sum()),
                          n_high_pip=len(calls), n_high_pip_causal=int(calls.causal.sum()),
                          n_high_pip_null=int(calls['null'].sum()),
                          n_high_pip_central_null=int(calls.central_null.sum()),
                          n_central_null_evaluated=int(group.central_null.sum()))
            match = metrics.loc[(metrics.pipeline == pipeline) & (metrics.fit_mode == mode) &
                                (metrics.method == method)]
            require(len(match) == 1, 'Expected one canonical pooled metric row.')
            for key, value in record.items():
                require(match.iloc[0][key] == value, 'Canonical pooled recount mismatch: ' + key)
            pooled.append(record)
        summary, detail = pve_group_summaries(data, pipeline, mode)
        summaries.append(summary)
        details.append(detail)
    summary = pd.concat(summaries, ignore_index=True)
    differences = {}
    if args.supervisor_groups:
        old = pd.read_csv(args.supervisor_groups, sep='\t')
        reproduced = summary.loc[(summary.finemap == 'SuSiE-adj') | (summary.fit_mode == 'ERV_FALSE')]
        compared = old.merge(reproduced, on=['method', 'finemap', 'group'],
                             suffixes=('_upload', '_recomputed'), validate='one_to_one')
        require(len(compared) == len(old) == 24, 'Unexpected supervisor group-table inventory.')
        for key in ['n', 'calls', 'pooled', 'power', 'lo', 'hi', 'min_n']:
            differences[key] = float(np.max(np.abs(compared[key + '_upload'] - compared[key + '_recomputed'])))
            require(differences[key] < 1e-12, 'Supervisor estimator reproduction failed: ' + key)
    args.output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output / 'pve_groups_all_runs.tsv', sep='\t', index=False, float_format='%.15g')
    pd.concat(details, ignore_index=True).to_csv(args.output / 'pve_groups_by_phenotype.tsv',
                                               sep='\t', index=False, float_format='%.15g')
    pd.DataFrame(pooled).to_csv(args.output / 'independently_recounted_pooled_metrics.tsv', sep='\t', index=False)
    evidence = []
    final = summary.loc[summary.fit_mode == 'ERV_TRUE']
    for label in ['<0.5×', 'Non-causal', '≥4×']:
        block = final.loc[final.group == label]
        evidence.append(dict(statistic='mean_phenotype_high_pip_rate', group=label,
                             minimum=float(block.power.min()), maximum=float(block.power.max()), units='fraction'))
    evidence.append(dict(statistic='supplied_recovered_pve', group='all_causal_pairs',
                         minimum=float(pve.recovered_pve.min()), maximum=float(pve.recovered_pve.max()), units='fraction'))
    pd.DataFrame(evidence).to_csv(args.output / 'claim_evidence.tsv', sep='\t', index=False, float_format='%.15g')
    inputs = [('recovered_effects', args.effects), ('fixed_causal_truth', args.causal),
              ('pooled_metrics', args.pooled_metrics)]
    inputs += [(pipeline + '_' + mode, path) for pipeline, mode, path in sources]
    if args.supervisor_groups:
        inputs.append(('supervisor_group_baseline', args.supervisor_groups))
    manifest = dict(status='complete', pooled_counts_match_saved_metrics=True,
                    supervisor_table_reproduced=bool(args.supervisor_groups),
                    supervisor_table_max_absolute_differences=differences,
                    grouping='Left-closed PVE / (h2/m) bins; equal-phenotype means; sample SD; rounded t9 critical value 2.262157',
                    n_phenotypes=10, expected_pve_fraction=0.5 / len(causal),
                    maximum_recovered_pve_percent=float(100 * pve.recovered_pve.max()),
                    effect_source='Author-supplied recovered effects; recovery generator not supplied',
                    detailed_outputs_written=False,
                    inputs=[fingerprint(role, path) for role, path in inputs], source_sha256=sha256(__file__))
    (args.output / 'finemapping_pve_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('Pooled counts verified; recovered-PVE aggregates written.')


if __name__ == '__main__':
    main()
