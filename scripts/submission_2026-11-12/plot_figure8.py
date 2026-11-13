#!/usr/bin/env python3
"""Updated PR and recovered-PVE figure, using final ERV_TRUE fine-mapping."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

from figure_style import (ROOT, DATA_ROOT, PIPELINE_COLORS, PIPELINE_MARKERS,
                          PIPELINE_LINESTYLES, configure, save_figure)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DATA_ROOT)
    parser.add_argument('--output', type=Path, default=ROOT)
    args = parser.parse_args()
    output = args.output.resolve()
    for name in ['captions', 'audit']:
        (output / name).mkdir(parents=True, exist_ok=True)
    configure(8.5)
    curve_path = args.source / 'figure08_pr_display_vertices.tsv'
    metric_path = args.source / 'figure08_pooled_metrics.tsv'
    group_path = args.source / 'figure08_pve_group_summary.tsv'
    curves = pd.read_csv(curve_path, sep='\t')
    metrics = pd.read_csv(metric_path, sep='\t')
    metrics = metrics.loc[metrics.fit_mode == 'ERV_TRUE'].copy()
    groups = pd.read_csv(group_path, sep='\t')
    assert set(curves.fit_mode) == {'ERV_TRUE'} and set(groups.fit_mode) == {'ERV_TRUE'}
    assert set(metrics.pipeline) == {'SuSiE-adj', 'SuSiE-RSS'} and len(metrics) == 4
    fig, axes = plt.subplots(2, 2, figsize=(170 / 25.4, 155 / 25.4))
    fig.subplots_adjust(left=0.105, right=0.96, bottom=0.155, top=0.895, hspace=0.62, wspace=0.32)
    pipelines = ['SuSiE-adj', 'SuSiE-RSS']
    labels = ['<0.5×', '0.5–1×', '1–2×', '2–4×', '≥4×']
    for column, method in enumerate(['BOLT-LMM', 'REGENIE-LOCO']):
        title = 'BOLT-LMM' if column == 0 else 'REGENIE_LOCO'
        ax = axes[0, column]
        baseline = None
        annotations = []
        for pipeline in pipelines:
            curve = curves.loc[(curves.method == method) & (curves.pipeline == pipeline)]
            curve = curve.loc[np.isfinite(curve.threshold) & np.isfinite(curve.precision)]
            curve = curve.sort_values('recall', kind='mergesort')
            assert len(curve) > 100
            ax.plot(curve.recall, curve.precision, color=PIPELINE_COLORS[pipeline],
                    linestyle=PIPELINE_LINESTYLES[pipeline], linewidth=1.25)
            record = metrics.loc[(metrics.method == method) & (metrics.pipeline == pipeline)].iloc[0]
            baseline = float(record.random_precision_baseline)
            annotations.append('{}: {:.2f}%'.format(pipeline, 100 * record.average_precision))
        ax.axhline(baseline, color='#555555', linestyle=':', linewidth=0.9)
        ax.text(0.97, 0.96, 'AP\n' + '\n'.join(annotations), transform=ax.transAxes,
                ha='right', va='top', fontsize=8.0)
        ax.text(0.98, baseline + 0.0035, 'Random: {:.2f}%'.format(100 * baseline),
                ha='right', va='bottom', fontsize=8.0, color='#444444')
        ax.set(xlim=(0, 1), ylim=(0, 0.16), xlabel='Recall', ylabel='Precision' if column == 0 else '')
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_yticks([0, 0.04, 0.08, 0.12, 0.16])
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
        ax.set_title(title, fontsize=9, fontweight='bold', loc='left', pad=5)
        ax.grid(axis='y')

        ax = axes[1, column]
        counts = None
        for index, pipeline in enumerate(pipelines):
            block = groups.loc[(groups.method == method) & (groups.finemap == pipeline)].set_index('group')
            selected = block.loc[labels]
            if counts is not None:
                assert np.array_equal(counts, selected.n.to_numpy())
            counts = selected.n.to_numpy()
            dx = -0.10 if index == 0 else 0.10
            y = selected.power.to_numpy()
            ax.errorbar(np.arange(5) + dx, y,
                        yerr=np.vstack([y - selected.lo.to_numpy(), selected.hi.to_numpy() - y]),
                        color=PIPELINE_COLORS[pipeline], marker=PIPELINE_MARKERS[pipeline],
                        linestyle='none', markersize=4.1, markeredgewidth=0.5,
                        markeredgecolor='white', elinewidth=0.9, capsize=2.0, capthick=0.7)
            ax.axhline(float(block.loc['Non-causal', 'power']), color=PIPELINE_COLORS[pipeline],
                       linestyle=':', linewidth=0.9, alpha=0.8)
        ax.set_xticks(np.arange(5))
        ax.set_xticklabels(['{}\n({:,})'.format(label, int(n)) for label, n in zip(labels, counts)], fontsize=8.0)
        ax.set(xlim=(-0.5, 4.5), ylim=(0, 0.03),
               xlabel='Recovered PVE / expected PVE',
               ylabel='SNPs with PIP ≥ 0.9' if column == 0 else '')
        ax.set_yticks([0, 0.01, 0.02, 0.03])
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
        ax.grid(axis='y')
    handles = [Line2D([], [], color=PIPELINE_COLORS[name], linestyle=PIPELINE_LINESTYLES[name],
                      marker=PIPELINE_MARKERS[name], markersize=3.5, linewidth=1.1, label=name)
               for name in pipelines]
    handles += [Line2D([], [], color='#555555', linestyle=':', linewidth=0.9, label='Random ranking (a)')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.55, 0.995),
               ncol=3, frameon=False, columnspacing=1.6, handlelength=2.1)
    fig.text(0.015, 0.91, 'a', fontweight='bold', fontsize=11)
    fig.text(0.015, 0.475, 'b', fontweight='bold', fontsize=11)
    fig.text(0.55, 0.052, 'Parentheses: causal SNP–phenotype pairs', ha='center', fontsize=8.0)
    fig.text(0.55, 0.024, 'Dotted lines in b: non-causal SNPs', ha='center', fontsize=8.0)
    save_figure(fig, 'Figure_08', height_mm=155, out_root=output)
    plt.close(fig)
    caption = '''Figure 8. Fine-mapping ranking and recovery by recovered causal PVE.
(a) Precision–recall curves for BOLT-LMM and REGENIE_LOCO loci, ranking SNP–phenotype pairs pooled across ten phenotypes by PIP. SuSiE-adj and SuSiE-RSS both use estimate_residual_variance = TRUE. Each SNP contributes once per phenotype and pipeline, retaining its maximum PIP across loci. AP denotes pooled stepwise average precision; dotted gray lines mark random-ranking precision, equal to causal prevalence in each analyzed pool. Precision is displayed from 0% to 16%; all empirical curve values fall within this range. The artificial zero-call endpoint is omitted. (b) Fraction of causal SNPs reaching PIP ≥0.9, grouped by author-supplied recovered PVE relative to the design expectation h²/m = 0.5/10,000 = 0.005%. Points are equal-phenotype mean call rates, and error bars are 95% t-based confidence intervals across ten phenotypes. Colored dotted lines give the corresponding non-causal call rates for each pipeline. Parentheses show analyzed causal SNP–phenotype counts. Recovered PVE is distinct from PVE estimated using GWAS association coefficients. Standard RSS uses ordinary in-sample LD and N = 20,000; SuSiE-adj uses LOCO-adjusted LD and locus-specific effective sample sizes. Both use L = 10, 95% credible sets, minimum absolute correlation 0.5 and at most 1,000 iterations. All 6,732 BOLT-LMM and 6,569 REGENIE_LOCO fits converged in each pipeline. Exact pooled metrics and phenotype-specific PVE-group summaries accompany the figure.'''
    (output / 'captions/Figure_08_caption.txt').write_text(caption + '\n')
    manifest = dict(status='complete', final_RSS='ERV_TRUE', adjusted='ERV_TRUE',
                    row_a='PR curves; BOLT-LMM left, REGENIE_LOCO right',
                    row_b='Recovered-PVE groups; BOLT-LMM left, REGENIE_LOCO right',
                    sources=[dict(path=path.name, sha256=sha(path)) for path in
                             [curve_path, metric_path, group_path]],
                    source_sha256=sha(__file__))
    (output / 'audit/figure08_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('Figure 8 exported using final TRUE fits and recomputed recovered-PVE groups.')


if __name__ == '__main__':
    main()
