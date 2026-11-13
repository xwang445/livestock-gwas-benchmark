"""Cross-check current submission summaries against independent analysis tables."""
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'derived_data/submission_2026-11-12'


class SubmissionResultChecks(unittest.TestCase):
    def test_table_s1_uses_pooled_canonical_counts(self):
        canonical = pd.read_csv(ROOT / 'derived_data/revision_2026-11-05/association/association_by_phenotype.tsv', sep='\t')
        actual = pd.read_csv(SOURCE / 'tables/Table_S1_locus_distribution.tsv', sep='\t').set_index('configuration')
        for row in actual.itertuples():
            group = canonical[canonical.configuration.eq(row.Index)]
            self.assertEqual(len(group), 10)
            self.assertEqual(int(group.loci.sum()), row.pooled_loci)
            for name in ['effect', 'null', 'central_null']:
                count = int(group[name + '_loci'].sum())
                self.assertEqual(count, getattr(row, 'pooled_' + name + '_loci'))
                self.assertAlmostEqual(100 * count / row.pooled_loci, getattr(row, 'pooled_' + name + '_percent'), places=10)
        self.assertEqual(actual.loc['SLEMM', 'pooled_effect_percent'], 93.75)
        self.assertAlmostEqual(actual.loc['GCTA_MLMA', 'pooled_effect_percent'], 100 * 14 / 15)

    def test_spillover_ci_includes_all_phenotypes(self):
        by = pd.read_csv(SOURCE / 'tables/Table_S1_by_phenotype.tsv', sep='\t')
        summary = pd.read_csv(SOURCE / 'tables/Table_S1_central_null_spillover.tsv', sep='\t').set_index('configuration')
        for method, group in by.groupby('configuration'):
            rates = 100 * group.central_null_loci / group.central_null_windows
            half = t.ppf(.975, 9) * rates.std(ddof=1) / np.sqrt(10)
            self.assertEqual(len(rates), 10)
            self.assertAlmostEqual(summary.loc[method, 'spillover_ci95_lower_percent'], rates.mean() - half, places=8)
            self.assertAlmostEqual(summary.loc[method, 'spillover_ci95_upper_percent'], rates.mean() + half, places=8)

    def test_recovered_pve_groups_match_final_phenotype_rates(self):
        per = pd.read_csv(ROOT / 'derived_data/revision_2026-11-11/finemapping/pve_groups_by_phenotype.tsv', sep='\t')
        per = per[per.fit_mode.eq('ERV_TRUE')]
        summary = pd.read_csv(SOURCE / 'figure08_pve_group_summary.tsv', sep='\t')
        for row in summary.itertuples():
            group = per[per.method.eq(row.method) & per.pipeline.eq(row.finemap) & per.group.eq(row.group)]
            self.assertEqual(len(group), 10)
            self.assertEqual(int(group.n.sum()), row.n)
            self.assertEqual(int(group.calls.sum()), row.calls)
            rates = group.calls / group.n
            half = 2.262157 * rates.std(ddof=1) / np.sqrt(10)
            self.assertAlmostEqual(row.power, rates.mean(), places=12)
            self.assertAlmostEqual(row.lo, rates.mean() - half, places=12)
            self.assertAlmostEqual(row.hi, rates.mean() + half, places=12)

    def test_figure8_metrics_use_current_rss(self):
        data = pd.read_csv(SOURCE / 'figure08_pooled_metrics.tsv', sep='\t')
        self.assertEqual(len(data), 4)
        self.assertEqual(set(data.fit_mode), {'ERV_TRUE'})
        rss = data[data.pipeline.eq('SuSiE-RSS')].set_index('method')
        self.assertEqual(rss.loc['BOLT-LMM', 'n_high_pip'], 2021)
        self.assertEqual(rss.loc['REGENIE-LOCO', 'n_high_pip'], 2072)
        self.assertTrue(np.allclose(data.random_precision_baseline, data.n_causal_evaluated / data.n_evaluated))


if __name__ == '__main__':
    unittest.main()
