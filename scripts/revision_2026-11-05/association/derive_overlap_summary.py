#!/usr/bin/env python3
"""Summarize the canonical per-phenotype table; no GWAS parsing or fitting."""
import argparse
from pathlib import Path
import pandas as pd
from derive_association_tables import save,overlap_summary

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--tables',type=Path,required=True);a=ap.parse_args()
    d=pd.read_csv(a.tables/'pairwise_locus_overlap.tsv',sep='\t',keep_default_na=False,na_values=['NA'])
    save(overlap_summary(d),a.tables,'pairwise_locus_overlap_summary.tsv')
