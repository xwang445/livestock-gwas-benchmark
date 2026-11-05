#!/usr/bin/env python3
"""Run the canonical association and lead-SNP summary; no duplicate estimator."""
from pathlib import Path
import runpy

runpy.run_path(str(Path(__file__).resolve().parents[2] /
    'scripts/revision_2026-11-05/association/derive_association_tables.py'),run_name='__main__')
