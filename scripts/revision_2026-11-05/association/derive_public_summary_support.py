#!/usr/bin/env python3
"""Export aggregate support for current public summary views from private audits.

Detailed loci supply only cross-replicate intersection/union counts. No locus IDs,
SNP IDs, individual records or absolute source paths are written.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

METHODS = ['BOLT-LMM', 'REGENIE_LOCO', 'REGENIE_WO_LOCO', 'fastGWA',
           'FarmCPU_NYC', 'FarmCPU_Kansas', 'BLINK_NYC', 'BLINK_Kansas', 'SLEMM', 'GCTA_MLMA']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def cross_replicate_stability(loci, methods=METHODS, phenotypes=range(1, 11)):
    require(not loci.duplicated(['configuration', 'phenotype', 'window_id']).any(), 'Repeated detected locus.')
    rows = []
    for method in methods:
        sets = {p: set(loci.loc[loci.configuration.eq(method) & loci.phenotype.eq(p), 'window_id'])
                for p in phenotypes}
        for a in phenotypes:
            for b in phenotypes:
                first, second = sets[a], sets[b]
                n_intersection, n_union = len(first & second), len(first | second)
                rows.append(dict(configuration=method, phenotype_a=a, phenotype_b=b,
                                 n_intersection=n_intersection, n_union=n_union,
                                 jaccard=n_intersection / n_union if first and second else np.nan,
                                 jaccard_set_theoretic=n_intersection / n_union if n_union else np.nan))
    return pd.DataFrame(rows)


def safe_source_path(path, benchmark_root, gcta_root):
    path = Path(path).resolve()
    for role, root in [('benchmark', benchmark_root), ('gcta', gcta_root)]:
        try:
            relative = path.relative_to(Path(root).resolve())
        except ValueError:
            continue
        return role + '/' + relative.as_posix()
    raise ValueError('Source is outside the explicitly supplied benchmark/GCTA roots.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ['qc', 'detected-loci', 'input-manifest', 'benchmark-root', 'gcta-root', 'output']:
        parser.add_argument('--' + option, type=Path, required=True)
    args = parser.parse_args()
    qc = pd.read_csv(args.qc, sep='\t')
    qc = qc.loc[qc.configuration.isin(METHODS), ['configuration', 'phenotype', 'rows', 'missing_p', 'zero_p']].copy()
    manifest = pd.read_csv(args.input_manifest, sep='\t')
    manifest = manifest.loc[manifest.configuration.isin(METHODS)].copy()
    expected = {(method, phenotype) for method in METHODS for phenotype in range(1, 11)}
    for table in [qc, manifest]:
        require(not table.duplicated(['configuration', 'phenotype']).any(), 'Repeated configuration/phenotype.')
        require(set(zip(table.configuration, table.phenotype)) == expected, 'Incomplete current method coverage.')
        require(table.rows.eq(135358).all(), 'Unexpected marker panel size.')
    require(manifest.status.eq('derived_complete').all(), 'Unfinished source analysis.')
    require((qc.missing_p >= 0).all() and (qc.zero_p >= 0).all(), 'Negative source QC counts.')
    manifest['path'] = manifest.path.map(lambda path: safe_source_path(path, args.benchmark_root, args.gcta_root))
    manifest = manifest[['configuration', 'phenotype', 'path', 'sha256', 'rows']]
    loci = pd.read_csv(args.detected_loci, sep='\t', usecols=['configuration', 'phenotype', 'window_id'])
    loci = loci.loc[loci.configuration.isin(METHODS)]
    require(set(loci.configuration).issubset(METHODS), 'Unknown method in detected loci.')
    stability = cross_replicate_stability(loci)
    args.output.mkdir(parents=True, exist_ok=True)
    outputs = {'association_qc.tsv': qc.sort_values(['configuration', 'phenotype']),
               'current_input_manifest.tsv': manifest.sort_values(['configuration', 'phenotype']),
               'cross_replicate_stability.tsv': stability}
    for name, data in outputs.items():
        data.to_csv(args.output / name, sep='\t', index=False, na_rep='NA', float_format='%.17g')
    sources = [('association_qc', args.qc), ('canonical_detected_loci', args.detected_loci),
               ('source_input_manifest', args.input_manifest)]
    provenance = dict(status='complete', configurations=10, phenotypes=10,
                      cross_replicate_rows=len(stability), detailed_outputs_written=False,
                      jaccard='NA if either set is empty; set-theoretic alternative separately labeled',
                      source_paths='Explicit benchmark/ and gcta/ relative namespaces',
                      inputs=[dict(role=role, sha256=sha(path), bytes=path.stat().st_size) for role, path in sources],
                      source_sha256=sha(__file__))
    (args.output / 'public_summary_support_manifest.json').write_text(json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    main()
