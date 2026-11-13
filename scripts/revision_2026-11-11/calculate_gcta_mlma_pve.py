#!/usr/bin/env python3
"""Compute GWAS-estimated single-SNP variance fractions for completed GCTA MLMA.

These are estimates from MLMA coefficients, not saved generating-effect PVE.
Use the exact analysis animals and total observed phenotype sample variance.
All paths must be supplied explicitly. Only aggregate outputs are written.
"""
import argparse
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from common import fingerprint

import numpy as np
import pandas as pd


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def genotype_variance(prefix, n_samples, n_markers):
    """Exact A1 dosage sums from SNP-major BED, including missing-call counts."""
    path = Path(str(prefix) + '.bed')
    with path.open('rb') as handle:
        require(handle.read(3) == b'\x6c\x1b\x01', 'BED must be SNP-major PLINK binary.')
    stride = (n_samples + 3) // 4
    require(path.stat().st_size == 3 + stride * n_markers, 'Unexpected BED size.')
    raw = np.memmap(str(path), dtype=np.uint8, mode='r', offset=3,
                    shape=(n_markers, stride))
    codes = (np.arange(256, dtype=np.uint16)[:, None] >> (2 * np.arange(4))) & 3
    # PLINK codes 00,01,10,11 correspond to A1 dosages 2,missing,1,0.
    dosage = np.array([2, 0, 1, 0], dtype=np.int16)[codes]
    lut_sum = dosage.sum(axis=1)
    lut_square = (dosage * dosage).sum(axis=1)
    lut_missing = (codes == 1).sum(axis=1)
    sums = np.empty(n_markers, dtype=np.int64)
    squares = np.empty(n_markers, dtype=np.int64)
    missing = np.empty(n_markers, dtype=np.int64)
    for start in range(0, n_markers, 2048):
        end = min(start + 2048, n_markers)
        chunk = raw[start:end]
        sums[start:end] = lut_sum[chunk].sum(axis=1)
        squares[start:end] = lut_square[chunk].sum(axis=1)
        missing[start:end] = lut_missing[chunk].sum(axis=1)
        if n_samples % 4:
            padding = codes[chunk[:, -1], n_samples % 4:]
            pad_dosage = np.array([2, 0, 1, 0], dtype=np.int16)[padding]
            sums[start:end] -= pad_dosage.sum(axis=1)
            squares[start:end] -= (pad_dosage * pad_dosage).sum(axis=1)
            missing[start:end] -= (padding == 1).sum(axis=1)
    observed = n_samples - missing
    require((observed > 1).all(), 'A SNP has fewer than two observed genotypes.')
    centered_squares = squares - sums.astype(float) ** 2 / observed
    require((centered_squares >= -1e-8).all(), 'Negative dosage variance.')
    # Missing calls, if present, are imputed to the observed SNP dosage mean.
    variance = np.maximum(centered_squares, 0) / (n_samples - 1)
    return pd.DataFrame({
        'n_genotypes_observed': observed,
        'n_genotypes_missing': missing,
        'dosage_A1_sum': sums,
        'dosage_A1_sum_squares': squares,
        'freq_A1_observed': sums / (2.0 * observed),
        'genotype_sample_variance': variance,
    })


def summarize(values, phenotype, group):
    values = np.asarray(values, dtype=float)
    row = dict(phenotype=phenotype, group=group, n_snps=len(values))
    for statistic, function in [('mean', np.mean), ('median', np.median),
                                ('min', np.min), ('max', np.max)]:
        value = float(function(values)) if len(values) else None
        row[statistic + '_pve_fraction'] = value
        row[statistic + '_pve_percent'] = None if value is None else 100 * value
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ['genotype-prefix', 'reference-bim', 'causal', 'loci', 'phenotype-dir', 'gcta-dir', 'output']:
        parser.add_argument('--' + option, type=Path, required=True)
    args = parser.parse_args()
    prefix, gcta_dir, output = args.genotype_prefix, args.gcta_dir, args.output
    names = ['Chr', 'SNP', 'cm', 'bp', 'A1', 'A2']
    bim = pd.read_csv(str(prefix) + '.bim', sep=r'\s+', header=None, names=names)
    original_bim = pd.read_csv(args.reference_bim, sep=r'\s+', header=None, names=names)
    require(bim.equals(original_bim), 'Analysis subset and original GCTA BIM differ.')
    require(len(bim) == 135358 and bim.SNP.is_unique, 'Unexpected marker inventory.')
    fam = pd.read_csv(str(prefix) + '.fam', sep=r'\s+', header=None, dtype=str)
    analysis_ids = set(zip(fam[0], fam[1]))
    require(len(fam) == len(analysis_ids) == 20000, 'Expected exactly 20,000 unique samples.')
    causal_path = args.causal
    causal = set(pd.read_csv(causal_path).SNP)
    require(len(causal) == 10000 and causal.issubset(set(bim.SNP)), 'Unexpected causal SNP inventory.')
    loci_path = args.loci
    loci = pd.read_csv(loci_path, sep='\t')
    loci = loci.loc[loci.configuration == 'GCTA_MLMA'].copy()
    require(len(loci) == 15, 'Expected 15 canonical GCTA lead SNPs.')

    print('Computing observed dosage variances for the matched 20,000 animals...', flush=True)
    genotype = genotype_variance(prefix, len(fam), len(bim))
    genotype = pd.concat([bim[['Chr', 'SNP', 'bp', 'A1', 'A2']], genotype], axis=1)
    require((genotype.n_genotypes_missing == 0).all(),
            'Missing genotypes present: verify GCTA imputation convention before proceeding.')
    output.mkdir(parents=True, exist_ok=True)
    summary_rows, phenotype_rows = [], []
    pooled = defaultdict(list)
    inputs = [Path(str(prefix) + ext) for ext in ('.bed', '.bim', '.fam')]
    inputs += [causal_path, loci_path, args.reference_bim]
    frequency_difference = 0.0
    expected_pve = 0.5 / len(causal)
    for phenotype in range(1, 11):
        phenotype_path = args.phenotype_dir / 'pheno.{}.txt'.format(phenotype)
        pheno = pd.read_csv(phenotype_path, sep=r'\s+', dtype={'FID': str, 'IID': str})
        require(set(zip(pheno.FID, pheno.IID)) == analysis_ids and len(pheno) == 20000,
                'Phenotype sample mismatch: {}'.format(phenotype))
        y = pheno['pheno{}'.format(phenotype)].to_numpy(float)
        require(np.isfinite(y).all() and not (y == -9).any(), 'Invalid phenotypes.')
        y_variance = float(np.var(y, ddof=1))
        mlma_path = gcta_dir / 'spg.even.100k.pheno{}.mlma'.format(phenotype)
        log_path = gcta_dir / 'spg.even.100k.pheno{}.log'.format(phenotype)
        log = log_path.read_text()
        require('Analysis finished' in log and '20000 individuals are in common' in log,
                'Missing successful completion or wrong sample count: {}'.format(log_path))
        gwas = pd.read_csv(mlma_path, sep=r'\s+')
        for column in ('Chr', 'SNP', 'bp', 'A1', 'A2'):
            require(gwas[column].equals(bim[column]), 'GWAS/BIM mismatch: ' + column)
        require(np.isfinite(gwas[['Freq', 'b', 'se', 'p']].to_numpy()).all(), 'Nonfinite GWAS output.')
        require(gwas.p.between(0, 1).all() and (gwas.se > 0).all(), 'Invalid GWAS statistics.')
        difference = float(np.max(np.abs(gwas.Freq - genotype.freq_A1_observed)))
        frequency_difference = max(frequency_difference, difference)
        require(difference < 1e-6, 'GCTA allele frequencies mismatch the matched genotype data.')
        lead = loci.loc[loci.phenotype == phenotype].set_index('lead_snp')
        require(lead.index.is_unique, 'Repeated lead SNP within a phenotype.')
        gwas.insert(0, 'phenotype', phenotype)
        gwas['is_causal'] = gwas.SNP.isin(causal)
        gwas['is_significant'] = gwas.p <= 5e-8
        gwas['is_lead'] = gwas.SNP.isin(lead.index)
        gwas['genotype_sample_variance'] = genotype.genotype_sample_variance
        gwas['phenotype_sample_variance'] = y_variance
        gwas['gwas_estimated_pve_fraction'] = gwas.b ** 2 * genotype.genotype_sample_variance / y_variance
        gwas['gwas_estimated_pve_percent'] = 100 * gwas.gwas_estimated_pve_fraction
        gwas['hwe_approx_pve_fraction'] = 2 * gwas.Freq * (1 - gwas.Freq) * gwas.b ** 2 / y_variance
        gwas['ratio_to_design_expected_pve'] = gwas.gwas_estimated_pve_fraction / expected_pve
        require(np.isfinite(gwas.gwas_estimated_pve_fraction).all(), 'Nonfinite estimated PVE.')
        selected = gwas.loc[gwas.is_lead].copy()
        selected['window_id'] = selected.SNP.map(lead.window_id)
        selected['region_type'] = selected.SNP.map(lead.region_type)
        require(np.allclose(selected.p, selected.SNP.map(lead.lead_p), rtol=1e-10, atol=0),
                'Canonical lead P values mismatch GCTA output.')
        require(np.allclose(selected.b, selected.SNP.map(lead.lead_beta), rtol=1e-10, atol=0),
                'Canonical lead betas mismatch GCTA output.')
        groups = dict(all_tested=np.ones(len(gwas), dtype=bool),
                      causal_tested=gwas.is_causal,
                      significant=gwas.is_significant,
                      all_leads=gwas.is_lead,
                      causal_leads=gwas.is_lead & gwas.is_causal,
                      noncausal_leads=gwas.is_lead & ~gwas.is_causal)
        for group, mask in groups.items():
            values = gwas.loc[mask, 'gwas_estimated_pve_fraction'].to_numpy()
            summary_rows.append(summarize(values, phenotype, group))
            pooled[group].append(values)
        phenotype_rows.append(dict(phenotype=phenotype, n=len(y), phenotype_variance=y_variance,
                                   n_tested=len(gwas), n_significant=int(gwas.is_significant.sum()),
                                   n_leads=len(selected), n_causal_leads=int(selected.is_causal.sum())))
        inputs.extend([phenotype_path, mlma_path, log_path])
        print('Phenotype {} complete: {} leads, {} causal.'.format(
            phenotype, len(selected), int(selected.is_causal.sum())), flush=True)
    require(sum(row['n_leads'] for row in phenotype_rows) == 15 and
            sum(row['n_causal_leads'] for row in phenotype_rows) == 10,
            'Unexpected pooled lead inventory.')
    for group, arrays in pooled.items():
        summary_rows.append(summarize(np.concatenate(arrays), 'pooled', group))
    summaries = pd.DataFrame(summary_rows)
    summaries.to_csv(output / 'estimated_pve_summary.tsv', sep='\t', index=False, float_format='%.15g')
    pd.DataFrame(phenotype_rows).to_csv(output / 'phenotype_variances.tsv', sep='\t', index=False,
                                       float_format='%.15g')
    manifest = dict(
        status='complete', created_utc=datetime.now(timezone.utc).isoformat(),
        method='GCTA-MLMA', estimand='GWAS-estimated single-SNP component PVE; not generating-effect PVE',
        formula='MLMA_b^2 * Var(observed_A1_dosage, ddof=1) / Var(observed_phenotype, ddof=1)',
        n_samples=20000, n_variants=135358, n_phenotypes=10,
        n_leads=15, n_causal_leads=10, n_missing_genotypes=int(genotype.n_genotypes_missing.sum()),
        maximum_allele_frequency_difference=frequency_difference,
        design_h2=0.5, n_simulated_causal=10000, design_expected_pve_fraction=expected_pve,
        locus_definition='Fixed 500 kb windows; canonical min-P leads and deterministic tie rule',
        limitations=[
            'MLMA coefficients are estimated association effects; these are not saved simulation effects.',
            'Selection of significant lead SNPs can inflate their estimated effect magnitudes.',
            'Single-SNP variance components must not be summed across LD-correlated SNPs.',
            'Ratios to h2/m compare an estimated component with the design expectation; they do not recover realized generating-effect ratios.',
        ],
        detailed_outputs_written=False,
        inputs=[fingerprint('restricted_input_{:02d}'.format(i), path) for i, path in enumerate(inputs, 1)],
        source_sha256=sha256(__file__), numpy=np.__version__, pandas=pd.__version__,
    )
    (output / 'gcta_estimated_pve_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    report = [
        'GCTA-MLMA GWAS-estimated single-SNP component PVE',
        'Formula: b_MLMA^2 * sample Var(G) / sample Var(y), both on matched 20,000 animals.',
        'Observed genotype variance is used; no missing genotype calls were found.',
        'PVE fraction 0.01 corresponds to 1% of phenotypic variance.',
        '',
    ]
    for _, row in summaries.loc[(summaries.phenotype == 'pooled') &
                                summaries.group.isin(['all_leads', 'causal_leads', 'noncausal_leads'])].iterrows():
        report.append('{}: n={}; median={:.8g} ({:.6g}%); range={:.6g}%–{:.6g}%.'.format(
            row.group, int(row.n_snps), row.median_pve_fraction,
            row.median_pve_percent, row.min_pve_percent, row.max_pve_percent))
    report += ['', 'These estimates cannot replace generating-effect PVE ratios or causal-effect percentiles.',
               'Do not sum these components across correlated SNPs.']
    (output / 'REPORT.txt').write_text('\n'.join(report) + '\n')
    print('\n'.join(report), flush=True)


if __name__ == '__main__':
    main()
