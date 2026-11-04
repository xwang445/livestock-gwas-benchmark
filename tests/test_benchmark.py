import importlib.util
from pathlib import Path
import sys, tempfile, unittest
import numpy as np
import pandas as pd
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import benchmark_core as b
spec=importlib.util.spec_from_file_location('synthetic',ROOT/'example/synthetic/generate.py')
synthetic=importlib.util.module_from_spec(spec);spec.loader.exec_module(synthetic)

class ScientificDefinitions(unittest.TestCase):
    def test_output_cannot_modify_historical_tree(self):
        with self.assertRaises(ValueError):b.protect_source('/tmp/source','/tmp/source/outputs')
        with self.assertRaises(ValueError):b.protect_source('/tmp/source','/tmp/source')
        b.protect_source('/tmp/source','/tmp/new_outputs')
    def test_locus_threshold_and_boundaries(self):
        d=pd.DataFrame({'chr':[1,1,1,2],'bp':[499999,500000,1000000,500000],'p':[1e-9,1e-9,5e-8,1e-9]})
        self.assertEqual(b.locus_keys(d),{(1,0),(1,1),(2,1)})
        self.assertIn((1,2),b.locus_keys(d,inclusive=True))
    def test_central_null_closed_coordinates(self):
        d=pd.DataFrame({'chr':[1]*5,'bp':[3000000,3000001,3500000,5000000,5000001],'p':[1e-9]*5})
        null=pd.DataFrame({'chr':[1],'start':[2000000],'end':[6000000]})
        mask,keys=b.central_mask_and_keys(d,null)
        self.assertEqual(mask.tolist(),[False,True,True,True,False]);self.assertEqual(len(keys),2)
    def test_spatial_shared_boundary_counts_twice(self):
        d=pd.DataFrame({'chr':[1],'bp':[2050000],'p':[1e-9]})
        null=pd.DataFrame({'chr':[1],'start':[2000000],'end':[6000000]})
        x=b.spatial_enrichment(d,null)
        self.assertEqual(x.n_tested.sum(),2);self.assertEqual(x.n_significant.sum(),2)
    def test_lambda_historical_finite_tail_exclusion(self):
        p=np.array([0.,.5,.9]);d=pd.DataFrame({'method':['BLINK_NYC']*3,'p':p})
        observed,n,_=b.historical_lambda(d)
        self.assertEqual(n,2)
        self.assertAlmostEqual(observed,np.median(chi2.ppf(1-p[1:],1))/chi2.ppf(.5,1))
        self.assertAlmostEqual(b.lambda_gc([.5,.5,.5]),1.)
        with self.assertRaises(ValueError):b.lambda_gc([-.1])
    def test_overlap_empty_and_union(self):
        self.assertEqual(b.jaccard({1,2},{2,3}),1/3)
        self.assertTrue(np.isnan(b.jaccard(set(),{1})))
        self.assertTrue(np.isnan(b.jaccard(set(),set())))
    def test_equation13_prediction_only_stable_bins(self):
        self.assertEqual(b.equal_count_bins([2,1,1,3,4],2).tolist(),[1,1,1,2,2])
        with self.assertRaises(ValueError):b.equal_count_bins([np.nan]*20)
        ncp=np.array([0.,1.,5.]);self.assertTrue(np.array_equal(1+ncp,[1,2,6]))

class SyntheticSmoke(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.base=Path(cls.temp.name)
        cls.input=synthetic.generate(cls.base/'input')
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()
    def test_all_nine_harmonizers(self):
        for method in b.METHODS:
            d=b.harmonize(b.input_path(self.input,method,1),method,1,36)
            self.assertTrue(set(b.SCHEMA).issubset(d));self.assertEqual(len(d),36)
            self.assertEqual(set(d.chr),{1,2});self.assertTrue(d.p.between(0,1).all())
            if method.startswith(('BLINK_','FarmCPU_')):self.assertTrue(d.beta.isna().all())
    def test_validation_rejects_corruption(self):
        d=b.harmonize(b.input_path(self.input,'BOLT-LMM',1),'BOLT-LMM',1,36)
        for column,value in [('p',1.01),('bp',None),('chr',19),('method','TYPO'),('phenotype',11)]:
            corrupt=d.copy();corrupt.loc[0,column]=value
            with self.assertRaises(ValueError):b.validate_harmonized(corrupt)
        with self.assertRaises(ValueError):b.validate_harmonized(pd.concat([d,d.iloc[:1]]))
    def test_full_downstream_and_resume(self):
        out=self.base/'outputs';cache=self.base/'cache'
        b.summarize(self.input,out,cache,expected_snps=36,phenotypes=[1,2])
        first=(out/'method_summary/significant_snps_by_pheno.tsv').read_bytes()
        b.summarize(self.input,out,cache,expected_snps=36,phenotypes=[1,2])
        self.assertEqual(first,(out/'method_summary/significant_snps_by_pheno.tsv').read_bytes())
        d=pd.read_csv(out/'method_summary/method_summary.tsv',sep='\t')
        self.assertEqual(set(d.method),set(b.METHODS));self.assertTrue(d.n_phenotypes.eq(2).all())
        self.assertTrue((out/'overlap_summary/locus_overlap_by_pheno.tsv').is_file())

if __name__=='__main__':unittest.main()
