"""Boundary regressions for the scientific decision rule and source-scale plots."""
import importlib.util
from pathlib import Path
import sys,unittest
from unittest.mock import patch
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts/revision_2026-11-05'))
sys.path.insert(0,str(ROOT/'scripts'))
from association.significance import ALPHA,significant,strict_for_historical_audit
import benchmark_core
spec=importlib.util.spec_from_file_location('current_manhattan',ROOT/'scripts/revision_2026-11-05/plot_manhattan.py')
manhattan=importlib.util.module_from_spec(spec);spec.loader.exec_module(manhattan)

class InclusiveAssociationBoundary(unittest.TestCase):
    def test_exact_raw_boundary_without_tolerance(self):
        p=[0.,np.nextafter(ALPHA,0.),ALPHA,np.nextafter(ALPHA,np.inf),1.,np.nan,-1.]
        self.assertEqual(significant(p).tolist(),[True,True,True,False,False,False,False])
        self.assertFalse(strict_for_historical_audit([ALPHA])[0])

    def test_native_log_boundary_without_exponentiation_or_tolerance(self):
        boundary=-np.log10(ALPHA)
        scores=np.array([np.nextafter(boundary,-np.inf),boundary,np.nextafter(boundary,np.inf),np.nan,-1.])
        self.assertEqual(significant(scores,is_log10=True).tolist(),[False,True,True,False,False])
        self.assertFalse(strict_for_historical_audit([boundary],is_log10=True)[0])
        d=pd.DataFrame({'chr':[1]*3,'bp':[1,500000,1000000],
            'p':10.**(-scores[:3]),'source_log10p':scores[:3]})
        self.assertEqual(benchmark_core.locus_keys(d),{(1,1),(1,2)})

    def test_manhattan_reader_preserves_native_decisions_and_display(self):
        boundary=-np.log10(ALPHA);n=135358
        scores=np.zeros(n);scores[:3]=[boundary,np.nextafter(boundary,np.inf),np.nextafter(boundary,-np.inf)]
        raw=pd.DataFrame({'CHROM':np.ones(n,dtype=int),'GENPOS':np.arange(1,n+1),
            'ID':['SYNTHETIC_VARIANT_'+str(i) for i in range(n)],'LOG10P':scores})
        with patch.object(manhattan.pd,'read_csv',return_value=raw):d=manhattan.load('synthetic.regenie','REGENIE_LOCO')
        self.assertEqual(d.significant.iloc[:3].tolist(),[True,True,False])
        np.testing.assert_array_equal(d.minus_log10p.iloc[:3].to_numpy(),scores[:3])

    def test_strict_summary_cannot_run_as_current_default(self):
        with self.assertRaisesRegex(ValueError,'legacy audit'):
            benchmark_core.summarize('/nonexistent','/nonexistent-output','/nonexistent-cache')

if __name__=='__main__':unittest.main()
