"""Check current public views against the complete canonical method registry."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import unittest

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('current_public_exports', ROOT / 'scripts/publish_current_summaries.py')
EXPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORT)


class CurrentSummaryExports(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = EXPORT.read_table(ROOT / 'config/revision_methods.tsv')
        cls.methods = set(cls.registry.configuration.replace(EXPORT.ALIASES))
        cls.canonical = EXPORT.read_table(ROOT / 'derived_data/revision_2026-11-05/association/association_by_phenotype.tsv')
        cls.canonical = cls.canonical[cls.canonical.configuration.isin(cls.registry.configuration)].copy()
        cls.canonical['method'] = cls.canonical.configuration.replace(EXPORT.ALIASES)
        cls.canonical = cls.canonical.set_index(['method', 'phenotype']).sort_index()

    def test_all_root_method_views_cover_current_registry(self):
        for relative in ['method_summary/significant_snps_by_pheno.tsv',
                         'lambda_gc_summary/lambda_gc_by_pheno.tsv',
                         'spillover_summary/spillover_by_pheno.tsv', 'input_manifest.tsv']:
            with self.subTest(path=relative):
                frame = EXPORT.read_table(ROOT / 'derived_data' / relative)
                self.assertEqual(set(frame.method), self.methods)
                self.assertEqual(len(frame), 100)
                self.assertFalse(frame.duplicated(['method', 'phenotype']).any())
                self.assertTrue(frame.groupby('method').phenotype.nunique().eq(10).all())

    def test_significant_counts_and_lambda_match_current_fits(self):
        counts = EXPORT.read_table(ROOT / 'derived_data/method_summary/significant_snps_by_pheno.tsv').set_index(
            ['method', 'phenotype']).sort_index()
        self.assertTrue(np.array_equal(counts.n_significant_snps, self.canonical.significant_snps))
        self.assertTrue(np.array_equal(counts.n_significant_loci, self.canonical.loci))
        self.assertTrue(np.array_equal(counts.strict_significant_snps, self.canonical.strict_significant_snps))
        self.assertEqual(int(counts.loc['GCTA_MLMA'].n_significant_snps.sum()), 58)
        self.assertEqual(int(counts.loc['SLEMM'].n_significant_loci.sum()), 16)
        lam = EXPORT.read_table(ROOT / 'derived_data/lambda_gc_summary/lambda_gc_by_pheno.tsv').set_index(
            ['method', 'phenotype']).sort_index()
        self.assertTrue(np.allclose(lam.lambda_gc, self.canonical.lambda_gc_central_null_historical, rtol=1e-14))
        self.assertTrue(np.allclose(lam.lambda_gc_p_based, self.canonical.lambda_gc_central_null, rtol=1e-14))

    def test_power_grid_retains_precise_cutoffs_and_all_methods(self):
        power = EXPORT.read_table(ROOT / 'derived_data/power_summary/power_by_pheno.tsv')
        self.assertEqual(set(power.method), self.methods)
        self.assertTrue(power.groupby(['method', 'phenotype']).alpha.nunique().eq(34).all())
        self.assertTrue(np.array_equal(power.alpha, power.alpha_hex.map(float.fromhex)))
        self.assertTrue(np.allclose(power.power, power.n_significant_causal / power.n_tested_causal))
        for method, group in power.groupby('method'):
            self.assertEqual(len(group[group.alpha.eq(5e-8)]), 10, method)

    def test_spatial_geometry_and_current_hit_counts_agree(self):
        spatial = EXPORT.read_table(ROOT / 'derived_data/power_summary/spatial_enrichment_by_pheno.tsv')
        self.assertEqual(set(spatial.method), self.methods)
        self.assertEqual(len(spatial), 4000)
        self.assertTrue(spatial.groupby('distance_index').n_tested.nunique().eq(1).all())
        self.assertTrue(np.allclose(spatial.pooled_observed_expected,
                                    spatial.n_significant / (spatial.n_tested * 5e-8), rtol=1e-13))

    def test_crossreplicate_overlap_is_complete_and_symmetric(self):
        crossrep = EXPORT.read_table(ROOT / 'derived_data/overlap_summary/cross_replicate_stability.tsv')
        self.assertEqual(set(crossrep.method), self.methods)
        self.assertEqual(len(crossrep), 1000)
        keys = ['method', 'phenotype_a', 'phenotype_b']
        reverse = crossrep.rename(columns={'phenotype_a': 'phenotype_b', 'phenotype_b': 'phenotype_a'})
        pd.testing.assert_frame_equal(crossrep.set_index(keys).sort_index(), reverse.set_index(keys).sort_index())
        diagonal = crossrep[crossrep.phenotype_a.eq(crossrep.phenotype_b)].set_index(['method', 'phenotype_a'])
        self.assertTrue(np.array_equal(diagonal.sort_index().n_union, self.canonical.loci))

    def test_heritability_distinguishes_reported_and_rounded_estimates(self):
        h2 = EXPORT.read_table(ROOT / 'derived_data/method_summary/heritability_from_logs.tsv')
        self.assertEqual(set(h2.method), {'fastGWA', 'SLEMM', 'GCTA_MLMA'})
        self.assertTrue(h2.groupby('method').phenotype.nunique().eq(10).all())
        gcta = h2[h2.method.eq('GCTA_MLMA')]
        self.assertTrue(np.allclose(gcta.h2, gcta.Vg_printed / (gcta.Vg_printed + gcta.Ve_printed)))
        self.assertTrue(gcta.precision_note.str.contains('5 decimal').all())
        self.assertTrue(h2[h2.method.eq('SLEMM')].source_path.str.contains('SLEMM_qf90_window3000').all())

    def test_check_mode_is_deterministic_and_does_not_write(self):
        expected = EXPORT.rendered_outputs(ROOT)
        self.assertEqual(expected, EXPORT.rendered_outputs(ROOT))
        before = {name: (ROOT / 'derived_data' / name).read_bytes() for name in expected}
        subprocess.run([sys.executable, str(ROOT / 'scripts/publish_current_summaries.py'), '--repo', str(ROOT), '--check'],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(before, {name: (ROOT / 'derived_data' / name).read_bytes() for name in expected})


if __name__ == '__main__':
    unittest.main()
