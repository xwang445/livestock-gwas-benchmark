"""Independent checks of recovered-effect estimands and restricted-input helpers."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts/revision_2026-11-11'
sys.path.insert(0, str(SCRIPTS))


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


lead = load('calculate_recovered_lead_pve')
fm = load('check_finemapping_paragraph')
gcta = load('calculate_gcta_mlma_pve')
distance = load('summarize_null_distance')


class RecoveredEffectEstimators(unittest.TestCase):
    def test_percentiles_use_full_pooled_truth_and_average_tie_ranks(self):
        effect = pd.DataFrame(dict(snp=['A', 'B', 'C'], pve1=[1., 2., 4.], pve2=[1., 2., 8.],
                                   in_perfectLD_group=[False, True, False]))
        loci = pd.DataFrame(dict(configuration=['SLEMM'] * 3, phenotype=[1, 2, 1],
                                 window_id=[0, 1, 2], lead_snp=['A', 'B', 'NONCAUSAL']))
        summary, long = lead.summarize_leads(effect, loci, {'A', 'B', 'C'}, n_phenotypes=2, h2=6.)
        row = summary.iloc[0]
        self.assertEqual(row.n_loci, 3)
        self.assertEqual(row.n_causal_leads, 2)
        self.assertAlmostEqual(row.median_ratio_to_expected, 0.75)
        # The two causal leads rank 1.5/6 and 3.5/6 in the full six-pair truth.
        self.assertAlmostEqual(row.median_pooled_causal_percentile, 100 * 2.5 / 6)
        self.assertEqual(row.n_causal_leads_perfectLD_flagged, 1)
        self.assertEqual(len(long), 6)

    def test_recovered_normalization_is_checked_against_supplied_variance(self):
        effect = pd.DataFrame(dict(pve1=[1., 4.], alpha1=[2., 4.]))
        reference = pd.DataFrame(dict(PHENOTYPE=[1], VAR_Y=[4.]))
        self.assertAlmostEqual(lead.normalization_check(effect, reference, 1).iloc[0].implied_phenotype_variance, 4.)
        reference.VAR_Y = 3.
        with self.assertRaisesRegex(ValueError, 'normalization'):
            lead.normalization_check(effect, reference, 1)

    def test_phenotype_mean_is_distinct_from_pooled_rate(self):
        records = []
        for k in range(1, 11):
            # One call among one pair in phenotype 1; zero among nine elsewhere.
            records += [dict(GWAS_method='SYNTHETIC', phenotype='pheno{}'.format(k), group='<0.5×',
                             call=(k == 1), recovered_pve=1e-6)] * (1 if k == 1 else 9)
        summary, per = fm.pve_group_summaries(pd.DataFrame(records), 'SuSiE-adj', 'ERV_TRUE')
        row = summary.iloc[0]
        self.assertAlmostEqual(row.power, 0.1)
        self.assertAlmostEqual(row.pooled, 1 / 82)
        expected_half_width = 2.262157 * np.std([1.] + [0.] * 9, ddof=1) / np.sqrt(10)
        self.assertAlmostEqual(row.hi - row.power, expected_half_width)
        self.assertAlmostEqual(row.power - row.lo, expected_half_width)
        self.assertEqual(len(per), 10)

    def test_bed_variance_ignores_padding_and_mean_imputes_missing(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / 'SYNTHETIC'
            # Five samples, two bytes per marker. Padding uses dosage-2 code 00.
            codes = [[0, 2, 3, 0, 3], [1, 2, 0, 3, 1]]
            packed = bytearray(b'\x6c\x1b\x01')
            for marker in codes:
                packed.append(sum(code << (2 * i) for i, code in enumerate(marker[:4])))
                packed.append(marker[4])
            Path(str(prefix) + '.bed').write_bytes(packed)
            result = gcta.genotype_variance(prefix, 5, 2)
            expected = [np.var([2, 1, 0, 2, 0], ddof=1), np.var([1, 1, 2, 0, 1], ddof=1)]
            np.testing.assert_allclose(result.genotype_sample_variance, expected)
            self.assertEqual(result.n_genotypes_missing.tolist(), [0, 2])

    def test_distance_shares_use_full_call_and_pool_denominators(self):
        data = pd.DataFrame(dict(GWAS_method=['SYNTHETIC'] * 6,
                                 PIP=[0.9, 0.99, 0.95, 0.9, 0.1, 0.1],
                                 is_null_window=[1, 1, 1, 0, 1, 0],
                                 is_central_null=[0, 0, 1, 0, 1, 0],
                                 dist_to_boundary=[199999, 500000, 1000000, 1500000, 200000, 0]))
        row = distance.summarize(data, 'SuSiE-adj', 'ERV_TRUE').iloc[0]
        self.assertEqual(row.n_high_pip, 4)
        self.assertEqual(row.n_evaluated, 6)
        self.assertAlmostEqual(row.high_pip_null_fraction, 3 / 4)
        self.assertAlmostEqual(row.high_pip_null_ge_200kb_fraction, 2 / 4)
        self.assertAlmostEqual(row.high_pip_null_ge_500kb_fraction, 2 / 4)
        self.assertAlmostEqual(row.high_pip_null_ge_1000kb_fraction, 1 / 4)
        self.assertAlmostEqual(row.evaluated_null_ge_200kb_fraction, 3 / 6)
        self.assertAlmostEqual(row.high_pip_central_null_fraction, 1 / 4)
        self.assertAlmostEqual(row.evaluated_central_null_fraction, 2 / 6)

    def test_spatial_boundary_convention_is_explicit(self):
        ids = ['1_2000000_A', '1_3000000_A', '1_5000000_A', '1_6000000_A', '1_700000_A', '1_600000_A']
        pve = pd.DataFrame(dict(SNP=ids[:5], phenotype=['pheno1'] * 5,
                                recovered_pve=[0.01, 0.075, 0.15, 0.3, 0.5]))
        raw = pd.DataFrame(dict(SNP=ids, phenotype=['pheno1'] * 6,
                                GWAS_method=['SYNTHETIC'] * 6, PIP=[0.9] * 6))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'synthetic.tsv'
            raw.to_csv(path, sep='\t', index=False)
            result = fm.prepare_pips(path, pve, set(ids[:5])).set_index('SNP').loc[ids[:4]]
            self.assertEqual(result['null'].tolist(), [True, True, True, False])
            self.assertEqual(result.central_null.tolist(), [False, True, False, False])

    def test_public_manifests_have_no_source_paths(self):
        import json
        for path in (ROOT / 'derived_data/revision_2026-11-11').rglob('*manifest.json'):
            manifest = json.loads(path.read_text())
            self.assertFalse(manifest['detailed_outputs_written'])
            for source in manifest['inputs']:
                self.assertEqual(set(source), {'role', 'bytes', 'sha256'})
                self.assertEqual(len(source['sha256']), 64)


if __name__ == '__main__':
    unittest.main()
