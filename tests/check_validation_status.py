"""Check corrected claims while retaining the original manuscript discrepancies."""
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
d=pd.read_csv(ROOT/'audit/manuscript_validation.tsv',sep='\t')
h=pd.read_csv(ROOT/'audit/historical_manuscript_validation.tsv',sep='\t')
assert set(h[h.status=='FAIL'].manuscript_item)=={'Results REGENIE maximum SNP count','Figure4 maximum minus log10 P approximately 102'}
assert not d.status.eq('FAIL').any(),d[d.status=='FAIL'].to_string()
for item,expected in [('Results REGENIE maximum SNP count',20461),('Figure4 corrected maximum minus log10 P',107.517)]:
 row=d.set_index('manuscript_item').loc[item];assert row.status=='PASS' and float(row.expected)==expected
assert d[d.manuscript_item=='Equation12 covariance prediction generator'].status.iloc[0]=='REQUIRES_RECOMPUTATION'
assert d[d.manuscript_item=='Maximum SNP PVE 0.063%-0.093%'].status.iloc[0]=='AMBIGUOUS'
assert d[d.manuscript_item.str.startswith(('Sample order ','Sample set ','Saved MPH GRM ','Phenotype chain'))].status.eq('PASS').all()
print('Corrected claims pass; original discrepancies and deferred provenance remain recorded.')
