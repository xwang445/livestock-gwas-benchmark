#!/usr/bin/env python3
"""Copy ten requested SLEMM GWAS outputs; verify originals and copies bytewise."""
import argparse,csv,hashlib,shutil
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);rows=[]
    for ph in range(1,11):
        fn='spg.even.100k.pheno%d.allchr.gwas.txt'%ph;src=a.source/fn;dest=a.output/fn
        h=sha(src)
        if not dest.exists():shutil.copy2(str(src),str(dest))
        if sha(dest)!=h:raise ValueError('Handoff file conflicts with source: '+str(dest))
        with dest.open() as f:
            header=f.readline().strip();n=sum(1 for _ in f)
        if n!=135358:raise ValueError('Incomplete output: '+str(dest))
        rows.append(dict(phenotype=ph,source_path=str(src),handoff_path=str(dest),sha256=h,size_bytes=dest.stat().st_size,n_variants=n,
            configuration='SLEMM --num_qf_markers 90 --window_size 3000',status='copied_and_checksum_verified'))
    with (a.output/'PI_handoff_manifest.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (a.output/'README.txt').write_text('Ten complete updated SLEMM merged chromosome 1–18 GWAS outputs. Source files are unchanged. No individual-level genotypes or sample records are included. Each output has 135,358 unique markers; checksums and copy validation are in PI_handoff_manifest.tsv. These are observed variant association outputs, not causal generating effects. Keep restricted until suitability for public release is reviewed.\n')

if __name__=='__main__':main()
