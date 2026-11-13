"""Current public scope, empty-set semantics and source-path privacy regressions."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import current_auxiliary_summaries as auxiliary

spec = importlib.util.spec_from_file_location('public_summary_support',
    ROOT / 'scripts/revision_2026-11-05/association/derive_public_summary_support.py')
support = importlib.util.module_from_spec(spec)
spec.loader.exec_module(support)


class CurrentAuxiliarySummaries(unittest.TestCase):
    def test_cross_replicate_empty_sets_are_explicit_and_ids_do_not_escape(self):
        loci = pd.DataFrame(dict(configuration=['SYNTHETIC'] * 4,
                                 phenotype=[2, 2, 3, 3], window_id=[11, 12, 12, 13]))
        result = support.cross_replicate_stability(loci, ['SYNTHETIC'], [1, 2, 3])
        result = result.set_index(['phenotype_a', 'phenotype_b'])
        self.assertEqual(len(result), 9)
        self.assertTrue(np.isnan(result.loc[(1, 2), 'jaccard']))
        self.assertEqual(result.loc[(1, 2), 'jaccard_set_theoretic'], 0.)
        self.assertTrue(np.isnan(result.loc[(1, 1), 'jaccard_set_theoretic']))
        self.assertEqual(result.loc[(2, 3), 'n_intersection'], 1)
        self.assertEqual(result.loc[(2, 3), 'n_union'], 3)
        self.assertAlmostEqual(result.loc[(2, 3), 'jaccard'], 1 / 3)
        self.assertNotIn('window_id', result.columns)

    def test_source_namespaces_require_explicit_authorized_roots(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            benchmark, gcta = base / 'benchmark', base / 'gcta'
            self.assertEqual(support.safe_source_path(benchmark / 'SLEMM_qf90_window3000/run.txt', benchmark, gcta),
                             'benchmark/SLEMM_qf90_window3000/run.txt')
            self.assertEqual(support.safe_source_path(gcta / 'run.mlma', benchmark, gcta), 'gcta/run.mlma')
            with self.assertRaisesRegex(ValueError, 'outside'):
                support.safe_source_path(base / 'unrelated/run.txt', benchmark, gcta)

    def test_current_fine_mapping_and_theory_scopes(self):
        tables = auxiliary.build_auxiliary_tables(ROOT)
        summary = tables['finemapping_summary/high_pip_summary.tsv']
        rss = summary.loc[summary.pipeline.eq('SuSiE-RSS')].set_index('method')
        self.assertEqual(rss.loc['BOLT-LMM', 'n_high_pip'], 2021)
        self.assertEqual(rss.loc['REGENIE_LOCO', 'n_high_pip'], 2072)
        self.assertEqual(rss.loc['BOLT-LMM', 'n_causal'], 261)
        self.assertEqual(rss.loc['REGENIE_LOCO', 'n_causal'], 254)
        self.assertEqual(set(summary.fit_mode), {'ERV_TRUE'})
        self.assertEqual(len(tables['finemapping_summary/high_pip_by_pheno.tsv']), 40)
        theory = tables['equation13_summary/equation13_summary.tsv']
        self.assertEqual(set(theory.method), {'BOLT-LMM', 'REGENIE_LOCO'})
        self.assertEqual(set(theory.region), {'Genome-wide'})
        self.assertEqual(len(theory), 40)
        self.assertEqual(len(tables['equation13_summary/equation13_by_pheno.tsv']), 400)
        self.assertEqual(int(tables['ld_summary/ld_decay_10kb.tsv'].n_pairs.sum()), 78082232)


if __name__ == '__main__':
    unittest.main()
