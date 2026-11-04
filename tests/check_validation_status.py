from pathlib import Path
import pandas as pd
p=Path(__file__).resolve().parents[1]/'audit/manuscript_validation.tsv'
d=pd.read_csv(p,sep='\t')
assert set(d[d.status=='FAIL'].manuscript_item)=={'Results REGENIE maximum SNP count','Figure4 maximum minus log10 P approximately 102'}
assert len(d[d.status=='PASS'])==95
assert d[d.manuscript_item=='Equation12 covariance prediction generator'].status.iloc[0]=='REQUIRES_RECOMPUTATION'
print('Known discrepancy and unresolved provenance remain explicitly reported.')
