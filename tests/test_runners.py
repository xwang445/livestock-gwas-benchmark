"""Mock external binaries to verify paths/arguments, without fitting any model."""
import json,os,subprocess,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class RunnerArguments(unittest.TestCase):
    def test_fastgwa_relative_output_not_nested(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder);mock=p/'gcta';log=p/'calls'
            mock.write_text('#!/usr/bin/env python3\nimport os,sys,json\nwith open(os.environ["MOCK_LOG"],"a") as f:f.write(json.dumps(sys.argv[1:])+"\\n")\n');mock.chmod(0o755)
            env=dict(os.environ,GENO_DIR='restricted/genotypes',PHENO_DIR='restricted/phenotypes',OUTPUT_DIR='new_results',GCTA_BIN=str(mock),MOCK_LOG=str(log))
            subprocess.run(['bash',str(ROOT/'scripts/04_gwas/run_fastGWA.sh')],cwd=p,env=env,check=True)
            calls=[json.loads(x) for x in log.read_text().splitlines()];self.assertEqual(len(calls),12)
            for call in calls:
                self.assertEqual(Path(call[call.index('--out')+1]).parent,(p/'new_results').resolve())
            refused=subprocess.run(['bash',str(ROOT/'scripts/04_gwas/run_fastGWA.sh')],cwd=p,env=env,capture_output=True)
            self.assertNotEqual(refused.returncode,0)
    def test_mph_historical_weight_seed_and_replicates(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder);mock=p/'mph';log=p/'calls'
            mock.write_text('#!/usr/bin/env python3\nimport os,sys,json\nwith open(os.environ["MOCK_LOG"],"a") as f:f.write(json.dumps(sys.argv[1:])+"\\n")\n');mock.chmod(0o755)
            for name in ['input.bed','input.bim','input.fam','causal.csv']:(p/name).touch()
            env=dict(os.environ,MPH_BIN=str(mock),MOCK_LOG=str(log));env.pop('NUM_PHENOTYPES',None);env.pop('HERITABILITY',None)
            subprocess.run(['bash',str(ROOT/'scripts/03_simulation/simulate_polygenic_phenotypes.sh'),'input','causal.csv','simulation'],cwd=p,env=env,check=True,capture_output=True)
            calls=[json.loads(x) for x in log.read_text().splitlines()];self.assertEqual(len(calls),2)
            call=calls[1];self.assertEqual(call[call.index('--num_phenotypes')+1],'1000');self.assertEqual(call[call.index('--seed')+1],'0')
            self.assertEqual((p/'simulation/grmlist.txt').read_text().split()[-1],'0.5')
if __name__=='__main__':unittest.main()
