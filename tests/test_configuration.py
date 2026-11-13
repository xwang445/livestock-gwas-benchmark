from pathlib import Path
import hashlib,json,sys,tempfile,unittest
import yaml
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from build_result_config import build
from check_config import public_paths,resolve_profile

class Configuration(unittest.TestCase):
    def test_results_match_executed_tables(self):
        actual=yaml.safe_load((ROOT/'config/verified_results.yaml').read_text())
        self.assertEqual(actual,build(ROOT))
        self.assertEqual(actual['historical_study_checks']['gwas_panels'],90)
        self.assertEqual(actual['study_checks']['planned_configurations'],10)
        self.assertEqual(actual['current_summary_exports']['completed_panels'],100)
        self.assertEqual(actual['current_summary_exports']['completed_configurations'],10)
        self.assertEqual(actual['significance']['rule'],'P <= 5e-8')
        self.assertEqual(len(actual['finemapping']['pooled']),4)
    def test_profile_availability_and_missing_input_detection(self):
        with tempfile.TemporaryDirectory() as folder:
            profile=Path(folder)/'profiles.yaml'
            profile.write_text(yaml.safe_dump({'repository_root':str(ROOT),'profiles':{
                'mac':{'benchmark_root':folder},'hpcc':{'benchmark_root':folder,'ssh_host':'synthetic-host'}}}))
            mac=resolve_profile(ROOT,'mac',profile)
            hpcc=resolve_profile(ROOT,'hpcc',profile)
            self.assertIn('equation13_cached_theory',mac['inputs'])
            self.assertNotIn('ld_pairs_10mb',mac['inputs'])
            self.assertIn('ld_pairs_10mb',hpcc['inputs'])
            self.assertNotIn('equation13_cached_theory',hpcc['inputs'])
            self.assertEqual(mac['input_check'],'NOT_RUN')
            with self.assertRaisesRegex(ValueError,'Missing configured inputs'):
                resolve_profile(ROOT,'mac',profile,check_inputs=True)
    def test_public_configuration_is_portable(self):
        config,methods=public_paths(ROOT)
        self.assertEqual(len(config['historical_result_paths']),27)
        self.assertTrue(all(p.startswith('derived_data/historical_2026-11-04/') for p in config['historical_result_paths'].values()))
        self.assertEqual(config['repository_paths']['derived_data'],'derived_data')
        self.assertEqual(len(methods),10)
        self.assertTrue(methods['SLEMM']['output_patterns']['hpcc'].startswith('SLEMM_qf90_window3000/'))
        for name in ['paths.yaml','verified_results.yaml','methods.yaml']:
            text=(ROOT/'config'/name).read_text()
            self.assertNotIn('/Users/',text)
            self.assertNotIn('/mnt/archive2/',text)
            self.assertNotIn('/path/to/',text)

    def test_historical_snapshot_content_is_preserved(self):
        manifest=json.loads((ROOT/'derived_data/historical_2026-11-04/snapshot_manifest.json').read_text())
        for row in manifest['records']:
            path=ROOT/row['archive_path']
            self.assertEqual(path.stat().st_size,row['bytes'])
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),row['sha256'])

if __name__=='__main__':unittest.main()
