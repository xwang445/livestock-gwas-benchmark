#!/usr/bin/env python3
"""Fingerprint public aggregate and figure artifacts; no restricted inputs read."""
import argparse,csv,hashlib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',default=Path(__file__).resolve().parents[1]);a=p.parse_args();root=Path(a.repo)
rows=[]
for folder in ['derived_data','figures/reconstructed','audit/reference_values','config','example/schemas']:
    for f in sorted((root/folder).rglob('*')):
        if f.is_file():
            if folder=='config' and f.name.startswith('local'):continue
            rows.append(dict(path=str(f.relative_to(root)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),classification='A_PUBLIC_CONFIG' if folder in ['config','example/schemas'] or f.suffix=='.md' else 'B_AGGREGATE_REVIEW_CANDIDATE'))
with (root/'audit/artifact_manifest.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['path','bytes','sha256','classification'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
print(f'Fingerprinted {len(rows)} review artifacts.')
