#!/usr/bin/env python3
"""Summarize author-supplied recovered PVE at canonical lead SNPs.

Inputs remain restricted. Outputs contain only aggregate summaries and fingerprints.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import causal_ids, effect_long, fingerprint, require, sha256


def summarize_leads(effect, loci, causal, n_phenotypes=10, h2=0.5):
    expected = h2 / len(causal)
    long = effect_long(effect, causal, n_phenotypes)
    long['ratio_to_expected'] = long.recovered_pve_fraction / expected
    # Average ranks, including ties, among every supplied causal SNP-phenotype pair.
    long['pooled_causal_percentile'] = long.ratio_to_expected.rank(pct=True, method='average') * 100
    flags = effect[['snp', 'in_perfectLD_group']].rename(columns={'snp': 'SNP'})
    long = long.merge(flags, on='SNP', validate='many_to_one')
    methods = ['SLEMM', 'GCTA_MLMA', 'BOLT-LMM', 'REGENIE_LOCO']
    loci = loci.loc[loci.configuration.isin(methods)].copy()
    require(not loci.duplicated(['configuration', 'phenotype', 'window_id']).any(), 'Repeated canonical loci.')
    # Use IDs to establish truth rather than relying on CSV boolean parsing.
    is_causal = loci.lead_snp.isin(causal)
    leads = loci.merge(long, left_on=['lead_snp', 'phenotype'], right_on=['SNP', 'phenotype'],
                       how='left', validate='many_to_one')
    require(leads.loc[is_causal.to_numpy(), 'recovered_pve_fraction'].notna().all(),
            'Causal lead missing from recovered-effect table.')
    rows = []
    for method, group in leads.groupby('configuration', sort=False):
        selected = group.loc[group.lead_snp.isin(causal)]
        rows.append(dict(method=method, n_loci=len(group), n_causal_leads=len(selected),
                         causal_lead_percent=100 * len(selected) / len(group),
                         median_recovered_pve_fraction=selected.recovered_pve_fraction.median(),
                         median_recovered_pve_percent=100 * selected.recovered_pve_fraction.median(),
                         median_ratio_to_expected=selected.ratio_to_expected.median(),
                         median_pooled_causal_percentile=selected.pooled_causal_percentile.median(),
                         percentile_q25=selected.pooled_causal_percentile.quantile(0.25),
                         percentile_q75=selected.pooled_causal_percentile.quantile(0.75),
                         n_causal_leads_perfectLD_flagged=int(selected.in_perfectLD_group.eq(True).sum())))
    result = pd.DataFrame(rows).set_index('method').reindex([m for m in methods if m in set(leads.configuration)])
    return result.reset_index(), long


def normalization_check(effect, reference, n_phenotypes=10):
    rows = []
    reference = reference.set_index('PHENOTYPE')
    for k in range(1, n_phenotypes + 1):
        alpha = effect['alpha{}'.format(k)].to_numpy(float)
        pve = effect['pve{}'.format(k)].to_numpy(float)
        require(np.isfinite(alpha).all(), 'Recovered alpha must be finite.')
        require(np.any(pve > 0), 'No positive PVE available for normalization check.')
        implied = alpha[pve > 0] ** 2 / pve[pve > 0]
        variance = float(reference.loc[k, 'VAR_Y'])
        require(np.isfinite(variance) and variance > 0, 'Invalid reference phenotype variance.')
        require(np.allclose(implied, variance, rtol=1e-10, atol=1e-10), 'Recovered normalization differs from reference.')
        rows.append(dict(phenotype=k, implied_phenotype_variance=float(np.median(implied)),
                         historical_31301_variance=variance,
                         maximum_absolute_difference=float(np.max(np.abs(implied - variance)))))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ['effects', 'loci', 'causal', 'reference-variances', 'output']:
        parser.add_argument('--' + option, type=Path, required=True)
    args = parser.parse_args()
    effect = pd.read_csv(args.effects, sep='\t', float_precision='round_trip')
    truth = causal_ids(args.causal)
    require(len(truth) == 10000, 'Expected the fixed 10,000 causal SNPs.')
    summary, long = summarize_leads(effect, pd.read_csv(args.loci, sep='\t'), truth)
    normalization = normalization_check(effect, pd.read_csv(args.reference_variances, sep='\t'))
    args.output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output / 'recovered_lead_pve_summary.tsv', sep='\t', index=False, float_format='%.15g')
    normalization.to_csv(args.output / 'pve_normalization_check.tsv', sep='\t', index=False, float_format='%.15g')
    effect_rows = []
    for phenotype, group in list(long.groupby('phenotype')) + [('pooled', long)]:
        row = dict(phenotype=phenotype, n_causal_pairs=len(group))
        for name in ['mean', 'median', 'min', 'max']:
            row[name + '_recovered_pve_fraction'] = float(getattr(group.recovered_pve_fraction, name)())
            row[name + '_ratio_to_expected'] = float(getattr(group.ratio_to_expected, name)())
        effect_rows.append(row)
    pd.DataFrame(effect_rows).to_csv(args.output / 'recovered_effect_summary.tsv',
                                    sep='\t', index=False, float_format='%.15g')
    manifest = dict(status='complete', n_causal_snps=len(truth), n_causal_snp_phenotype_pairs=len(long),
                    expected_pve_fraction=0.5 / len(truth), expected_pve_percent=100 * 0.5 / len(truth),
                    median_all_causal_ratio=float(long.ratio_to_expected.median()),
                    mean_all_causal_ratio=float(long.ratio_to_expected.mean()),
                    percentile_definition='Average ranks across all 100,000 supplied causal SNP-phenotype PVE ratios, rank/n*100',
                    lead_selection='Input canonical fixed 500 kb loci; minimum P, ties by position then SNP ID',
                    effect_source='Author-supplied recovered effects; recovery generator not supplied',
                    normalization='Supplied pve=alpha^2/phenotype variance of 31,301 simulation animals',
                    n_perfectLD_flagged_causal_snps=int(effect.in_perfectLD_group.sum()),
                    detailed_outputs_written=False,
                    inputs=[fingerprint(role, path) for role, path in [
                        ('recovered_effects', args.effects), ('canonical_loci', args.loci),
                        ('fixed_causal_truth', args.causal), ('reference_variances', args.reference_variances)]],
                    source_sha256=sha256(__file__))
    (args.output / 'recovered_lead_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
