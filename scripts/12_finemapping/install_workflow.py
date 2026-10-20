#!/usr/bin/env python3
"""Install code into the existing benchmark layout; never overwrite changed files."""
import argparse
from pathlib import Path
import shutil

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('benchmark_root',type=Path,help='Directory containing BOLTLMM/, REGENIE/, simu_pheno/, etc.')
    args=parser.parse_args()
    source=Path(__file__).resolve().parent
    target=args.benchmark_root.expanduser().resolve()
    trees=['SuSiE_adj_finemapping','SuSiE_finemapping']
    files=[p for name in trees for p in (source/name).rglob('*')
           if p.is_file() and p.suffix in {'.py','.R','.sh','.json'} and '__pycache__' not in p.parts]
    for p in files:
        dest=target/p.relative_to(source)
        if dest.exists() and dest.read_bytes()!=p.read_bytes():
            parser.error('Refusing to overwrite different existing code: '+str(dest)+'. Use a separate benchmark workspace.')
    for p in files:
        dest=target/p.relative_to(source)
        dest.parent.mkdir(parents=True,exist_ok=True)
        if not dest.exists():shutil.copy2(str(p),str(dest))
    directories={
        'SuSiE_adj_finemapping':['qc','logs','tmp','status','genotype','grm','loci'],
        'SuSiE_adj_finemapping/loco':['qc','logs','status','grm','ld','results','evaluation','figures','susie_results'],
        'SuSiE_finemapping':['qc','logs','status','ld','results','figures','susie_results','publication']}
    for root,names in directories.items():
        for name in names:(target/root/name).mkdir(parents=True,exist_ok=True)
    print('Installed {} code/configuration files under {}'.format(len(files),target))

if __name__=='__main__':main()
