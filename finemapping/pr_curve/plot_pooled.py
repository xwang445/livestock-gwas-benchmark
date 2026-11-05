#!/usr/bin/env python3
"""Redraw final Figure 8 from distributed pooled aggregates and safe vertices."""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'derived_data/revision_2026-11-05/susie'
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True,help='A new figure directory')
p.add_argument('--display-vertices',type=Path,default=DATA/'Figure_8_display_vertices.tsv')
p.add_argument('--metrics',type=Path,default=DATA/'pooled_finemapping_metrics.tsv')
a=p.parse_args()
subprocess.run([sys.executable,str(ROOT/'scripts/revision_2026-11-05/susie/plot_figure8_pooled.py'),
    '--display-vertices',str(a.display_vertices),'--metrics',str(a.metrics),
    '--output',str(a.output)],check=True)
